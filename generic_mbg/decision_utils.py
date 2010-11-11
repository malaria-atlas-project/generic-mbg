# Copyright (C) 2010 Anand Patil
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from prediction_utils import *
from inference_utils import close
import tables as tb

def kldiv(p, mu_pri, V_pri, likefn, norms):
    mu_like = p[0]
    V_like = p[1]
    if V_like < 0:
        return np.inf
    mu_post = (mu_like*V_pri + mu_pri*V_like)/(V_pri+V_like)
    V_post = V_pri*V_like/(V_pri+V_like)
    
    norms = norms*np.sqrt(V_post) + mu_post
    mean_loglike = likefn(norms)+pm.normal_like(norms, mu_pri, 1./V_pri)
    mean_logapprox_like = pm.normal_like(norms, mu_post, 1./V_post)
    
    return (mean_logapprox_like - mean_loglike)/len(norms)
    
def find_approx_params(mu_pri, V_pri, likefn, norms, match_moments, optimfn = None):
    """
    Returns the 'likelihood' mean and variance. 
    Matches posterior moments if match_moments=True,
    minimizes KL divergence with posterior if False.
    match_moments=False is better, but slower.
    """
    snorms = mu_pri+np.sqrt(V_pri)*norms
    l = np.array([likefn(s) for s in snorms])
    w = np.exp(l-l.max())
    m1 = np.sum(snorms*w)/w.sum()
    m2 = np.sum(snorms**2*w)/w.sum()
    
    mu_post_init = m1
    v_post_init = m2-m1**2
    
    mu_pri_ = np.mean(snorms)
    V_pri_ = np.var(snorms)
    
    v_like_init = v_post_init*V_pri/(V_pri_ - v_post_init)
    mu_like_init = (mu_post_init*(v_like_init+V_pri_)-mu_pri_*v_like_init)/V_pri_
    if v_like_init<0:
        raise RuntimeError, 'Crap.'
    
    if match_moments:
        p = [mu_like_init, v_like_init]
    else:
        f = close(kldiv, mu_pri=mu_pri, V_pri=V_pri, likefn=likefn, norms=norms)

        if optimfn is None:
            from scipy import optimize
            optimfn = optimize.fmin
        p = optimfn(f, [mu_like_init, v_like_init],disp=0)

    return p
    
def obs_corrections(mu_pri, C_pri, like_m, like_v, i):
    C_pri_scale = C_pri[i,:]/(C_pri[i,i]+like_v)
    C_corr = - np.outer(C_pri[i,:], C_pri_scale)
    mu_corr = (like_m-mu_pri[i])*C_pri_scale
    return np.ravel(mu_corr), C_corr

def calc_mean_like(like_m,like_v,mu_pri,v_pri):
    A = 1./v_pri + 1./like_v
    B = mu_pri/v_pri + like_m/like_v
    return -.5*(-B**2/A+like_m**2/like_v+mu_pri**2/v_pri)-.5*np.log(A*v_pri)

def calc_norm_const(norms, like_m, like_v, mu_pri, v_pri, likefn):
    s_pri = np.sqrt(v_pri)
    mean_like = pm.flib.logsum([likefn(n*s_pri+mu_pri) for n in norms])-np.log(len(norms))
    log_norm_const = mean_like-calc_mean_like(like_m,like_v,mu_pri,v_pri)
    return log_norm_const

def mTcImT(m,s):
    m_ = np.asarray(m.copy('F')).squeeze()
    pm.flib.dtrsm_wrap(s,m_,side='l',transa='n',uplo='l')
    return np.dot(m_,m_)

def logdet(s):
    return np.sum(np.log(np.diag(s)))*2

def calc_imp_weight(m_pri, c_pri, m_like, v_like, m_post, c_post):
    s_pri = np.linalg.cholesky(c_pri)
    s_post = np.linalg.cholesky(c_post)

    log_weight = -.5*(mTcImT(m_pri,s_pri)+np.sum(m_like*m_like/v_like)-mTcImT(m_post,s_post))+.5*(logdet(s_post)-logdet(s_pri))

    return log_weight

def mean_reduce_with_hdf(hf, n_reps):
    def mean_reduce_(sofar, next, name, ind, hf=hf, n_reps=n_reps):
        """A function to be used with hdf5_to_samps"""
        if hasattr(hf.root, name+'_mean'):
            hfa = getattr(hf.root, name+'_mean')
        else:
            hfa = hf.createCArray('/',name+'_mean',shape=(n_reps,)+next.shape,atom=tb.FloatAtom(),filters=tb.Filters(complevel=1,complib='zlib'))
        if sofar is None:
            hfa[ind] = next
        else:
            hfa[ind] = hfa[ind] + next
        return 1
    return mean_reduce_

def var_reduce_with_hdf(hf, n_reps):
    def var_reduce_(sofar, next, name, ind, hf=hf, n_reps=n_reps):
        """A function to be used with hdf5_to_samps"""
        if hasattr(hf.root, name+'_var'):
            hfa = getattr(hf.root, name+'_var')
        else:
            hfa = hf.createCArray('/',name+'_var',shape=(n_reps,)+next.shape,atom=tb.FloatAtom(),filters=tb.Filters(complevel=1,complib='zlib'))
        if sofar is None:
            hfa[ind] = next**2
        else:
            hfa[ind] = hfa[ind] + next**2
        return 1
    return var_reduce_

def histogram_reduce_with_hdf(bins, binfn, hf, n_reps):
    """Produces an accumulator to be used with hdf5_to_samps"""
    def hr(sofar, next, name, ind, hf=hf, n_reps=n_reps):
        if hasattr(hf.root, name+'_histogram'):
            hfa = getattr(hf.root, name+'_histogram')
        else:
            hfa = hf.createCArray('/',name+'_histogram',shape=(n_reps,)+next.shape+(len(bins),),atom=tb.FloatAtom(),filters=tb.Filters(complevel=1,complib='zlib'))
        if sofar is None:
            hfa[ind] = 0.
        # Call to Fortran function multiinc
        hist_ind = binfn(next)
        sofar = hfa[ind]
        multiinc(sofar,hist_ind)
        hfa[ind] = sofar
        return 1
    return hr

def find_joint_approx_params(mu_pri, C_pri, likefns, match_moments, approx_param_fn = None, tol=1.e-3, maxiter=100):
    if approx_param_fn is None:
        norms = np.random.normal(size=1000)
    else:
        raise NotImplementedError

    delta_m = np.ones_like(mu_pri)*np.inf
    delta_v = np.ones_like(mu_pri)*np.inf

    like_means = np.zeros_like(mu_pri)
    like_vars = np.ones_like(mu_pri)*np.inf
    
    norm_consts = np.ones_like(mu_pri)

    mu = mu_pri.copy()
    C = C_pri.copy()

    # Skip this to mock
    iter = 0
    while (np.any(np.abs(delta_m)>tol) or np.any(np.abs(delta_v/like_vars)>tol)) and iter < maxiter:
        iter += 1
        for i in xrange(len(mu_pri)):
            
            if not np.isinf(like_vars[i]):
                mu_corr, C_corr = obs_corrections(mu, C, like_means[i], -like_vars[i], i)            
                mu += mu_corr
                C += C_corr
            
            if approx_param_fn is None:
                new_like_mean, new_like_var = find_approx_params(mu[i], C[i,i], likefns[i], norms, match_moments)
            else:
                new_like_mean, new_like_var = approx_param_fn(mu[i], C[i,i], likefns[i])
    
            delta_m[i] = new_like_mean-like_means[i]
            delta_v[i] = new_like_var- like_vars[i]
            like_means[i] = new_like_mean
            like_vars[i] = new_like_var
            
            norm_consts[i] = calc_norm_const(norms, like_means[i], like_vars[i], mu[i], C[i,i], likefns[i])
                        
            mu_corr, C_corr = obs_corrections(mu, C, like_means[i], like_vars[i], i)
            mu += mu_corr
            C += C_corr
            
    if iter==maxiter:
        raise RuntimeError, 'EP algorithm failed to converge.'

    log_imp_weight = calc_imp_weight(mu_pri, C_pri, like_means, like_vars, mu, C)+np.sum(norm_consts)

    return like_means, like_vars, norm_consts, mu, C, log_imp_weight

def hdf5_to_survey_eval(M, x, nuggets, burn, thin, total, fns, postprocs, pred_covariate_dict, survey_x, survey_data, survey_covariate_dict, survey_likelihood, survey_plan, match_moments, finalize=None, continue_past_npd=False):
    
    hf=M.db._h5file
    gp_submods = list(set(filter(lambda c: isinstance(c,pm.gp.GPSubmodel), M.containers)))
    if len(gp_submods)>1:
        raise NotImplementedError, 'Not implemented for multiple fields yet.'
    
    f_labels = [gps.name for gps in gp_submods]

    products, postproc_args, extra_postproc_args = get_args(postprocs, fns, f_labels, M)
    sl_args,extra_sl_args = get_one_args(survey_likelihood, f_labels, M)
    extra_sl_args -= set(['i','data','survey_plan'])
        
    iter = np.arange(burn,all_chain_len(hf),thin)

    M_obs = {}
    C_obs = {}
    nugs = {}

    if len(iter)==0:
        raise ValueError, 'You asked for %i burnin iterations with thinnnig %i but the chains are only %i iterations long.'%(burn, thin, all_chain_len(hf))
    n_per = total/len(iter)+1
    actual_total = n_per * len(iter)
    time_count = -np.inf
    time_start = time.time()
    log_imp_weights = np.zeros((len(iter), len(survey_data)))
    
    for k in xrange(len(iter)):
        
        i = iter[k]
        # Restore the i'th cache fram
        all_chain_remember(M,i)
        
        # Add the covariate values on the prediction mesh to the appropriate caches.
        covariate_covariances = []
        for d in M.deterministics:
            if isinstance(d.value, pm.gp.Covariance):
                if isinstance(d.value.eval_fun,CovarianceWithCovariates):
                    d.value.eval_fun.add_values_to_cache(x,pred_covariate_dict)
                    d.value.eval_fun.add_values_to_cache(survey_x,survey_covariate_dict)
                    d.value.eval_fun.add_values_to_cache(np.vstack((d.value.eval_fun.meshes[0], survey_x[:,:2])), dict([(key, np.hstack((d.value.eval_fun.dicts[0][key], survey_covariate_dict[key]))) for key in survey_covariate_dict.iterkeys()])) 
        
        if time.time() - time_count > 10:
            print ((k*100)/len(iter)), '% complete',
            time_count = time.time()      
            if k > 0:      
                print 'expect results %s (in %s hours)'%(time.ctime((time_count-time_start)*len(iter)/float(k)+time_start),(time_count-time_start)*(len(iter)-float(k))/float(k)/3600)
            else:
                print

        base_M_preds = {}
        base_S_preds = {}
        norms = dict([(s, np.random.normal(size=n_per)) for s in gp_submods])
        for s in gp_submods:
            M_obs[s] = pm.utils.value(s.M_obs)
            C_obs[s] = pm.utils.value(s.C_obs)
            nugs[s] = pm.utils.value(nuggets[s])
            base_M_preds[s], base_V_pred = pm.gp.point_eval(M_obs[s], C_obs[s], x)
            base_S_preds[s] = np.sqrt(base_V_pred + nugs[s])
            
            apply_postprocs_and_reduce(M, n_per, base_M_preds, base_S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=-1, norms=norms[s])
        try:
            for l in xrange(len(survey_data)):   
                M_preds = {}
                S_preds = {}
                for s in gp_submods:

                    mu_pri = M_obs[s](survey_x)
                    C_pri = C_obs[s](survey_x, survey_x)+nugs[s]*np.eye(len(survey_x))
                    
                    # FIXME: This will only work for single fields currently.
                    closure_dict = {'data': survey_data[l], 'survey_plan': survey_plan}
                    
                    for extra_arg in extra_sl_args:
                        closure_dict[extra_arg] = pm.utils.value(getattr(M, extra_arg))
                    optim_fns = []
                    for ii in xrange(survey_data.shape[1]):
                        closure_dict_ = {'i':ii}
                        closure_dict_.update(closure_dict)
                        optim_fns.append(close(survey_likelihood,**closure_dict_))

                    like_means, like_vars, norm_consts, mu_post, C_post, log_imp_weight = find_joint_approx_params(mu_pri, C_pri, optim_fns, match_moments)
                    log_imp_weights[k,l]=log_imp_weight

                    M_obs_ = copy.copy(M_obs[s])
                    C_obs_ = copy.copy(C_obs[s])

                    pm.gp.observe(M_obs_, C_obs_, obs_mesh=survey_x, obs_vals = like_means, obs_V = like_vars)

                    M_preds[s], V_pred = pm.gp.point_eval(M_obs_, C_obs_, x)

                    if np.any(V_pred<0):
                        if continue_past_npd:
                            warnings.warn('Some elements of V_pred were negative. Assuming non-positive definiteness but not checking yet.')
                            actual_total -= n_per
                            log_imp_weights[k,:] = -np.inf
                            raise np.linalg.LinAlgError
                    S_preds[s] = np.sqrt(V_pred + nugs[s])    

                    # This is the time-consuming step. Makes sense I guess.
                    # TODO: Try it with only a 3dmap reduce.

                    apply_postprocs_and_reduce(M, n_per, M_preds, S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=l, norms=norms[s])
                
        
        except np.linalg.LinAlgError:
            continue
    
    for l in xrange(len(survey_data)):
        # Normalize importance weights
        log_imp_weights[:,l] -= pm.flib.logsum(log_imp_weights[:,l])        

    return actual_total, log_imp_weights