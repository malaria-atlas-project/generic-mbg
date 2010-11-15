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
import time
import warnings
from histogram_utils import impw, meanl, obsc, linit
import scipy
from scipy import stats

def kldiv(p, mu_pri, V_pri, likefn, norms):
    """
    Approximates the Kullback-Liebler divergence
    D( N(mu_pri,V_pri)*N(mu_like,V_like) || N(mu_pri,V_pri)*likefn )
    up to a constant of proportionality.
    """
    mu_like = p[0]
    V_like = p[1]
    if V_like < 0:
        return np.inf
    mu_post = (mu_like*V_pri + mu_pri*V_like)/(V_pri+V_like)
    V_post = V_pri*V_like/(V_pri+V_like)
    
    norms = norms*np.sqrt(V_post) + mu_post
    mean_loglike = likefn(norms).sum()+pm.normal_like(norms, mu_pri, 1./V_pri)
    mean_logapprox_like = pm.normal_like(norms, mu_post, 1./V_post)
    
    return (mean_logapprox_like - mean_loglike)/len(norms)
    
def kld2(like_mean, like_var, mu_pri, v_pri, likefn, norms):
    """
    Approximates the Kullback-Liebler divergence
    D( N(mu_pri,V_pri)*N(mu_like,V_like) || N(mu_pri,V_pri)*likefn ).
    """
    mu_post = (mu_pri*like_var + like_mean*v_pri)/(like_var+v_pri)
    v_post = v_pri*like_var/(v_pri+like_var)
    
    h = .5*np.log(2.*np.pi*np.exp(1)*v_post)
    
    pri_norms = norms*np.sqrt(v_pri)+mu_pri
    post_norms = norms*np.sqrt(v_post)+mu_post
    
    lr = likefn(pri_norms)
    lp = likefn(post_norms)
    pri_logdens = -(post_norms-mu_pri)**2/2./v_pri
    
    t1 = pm.flib.logsum(lr)-np.log(len(norms))
    t2 = np.mean(lp)
    t3 = np.log(1./np.sqrt(2.*np.pi*v_pri))+np.mean(pri_logdens)
    return h-t1+t2+t3

def plotcompare(mu_pri, V_pri, likefn, mu_like, V_like, xlo=-10, xhi=10):
    """
    Call from find_approx_params to check the approximation.
    """
    import pylab as pl
    x = np.linspace(xlo,xhi,101)
    pri = np.exp(-(mu_pri-x)**2/2./V_pri)
    like = np.exp(likefn(x))
    post = like*pri
    apost = np.exp(-(mu_like-x)**2/2./V_like)*pri
    for p_ in [pri , like , post , apost]:
        p_ /= p_.max()
    pl.clf()
    pl.plot(x,pri,'b-.')
    pl.plot(x,like,'r-.')
    pl.plot(x,post,'r-')
    pl.plot(x,apost,'g-')
    from IPython.Debugger import Pdb
    Pdb(color_scheme='Linux').set_trace() 
    
def find_approx_params(mu_pri, V_pri, likefn, norms, match_moments, optimfn = None, debug=False):
    """
    Returns the 'likelihood' mean and variance. 
    Matches posterior moments if match_moments=True,
    minimizes KL divergence with posterior if False.
    match_moments=False is better, but slower.
    """
    # TODO: Fortran.
    snorms = mu_pri+np.sqrt(V_pri)*norms
    l = likefn(snorms)

    mu_like_init, v_like_init = linit(l,snorms)

    if match_moments:
        p = [mu_like_init, v_like_init]
    else:
        f = close(kldiv, mu_pri=mu_pri, V_pri=V_pri, likefn=likefn, norms=norms)

        if optimfn is None:
            from scipy import optimize
            optimfn = optimize.fmin
        p = optimfn(f, [mu_like_init, v_like_init],disp=0)

    if debug:
        plotcompare (mu_pri, V_pri, likefn, p[0], p[1])

    return p
    
def calc_norm_const(norms, like_m, like_v, mu_pri, v_pri, likefn):
    """
    Finds the normalizing constant associated with the exponentiated quadratic
    approximation to the likelihood.
    """
    if v_pri <= 0:
        # This is evidence of extreme lack of fit.
        warnings.warn('Got a negative v_pri in calc_norm_const.')
        return -np.inf
    else:
        s_pri = np.sqrt(v_pri)
        mean_like = pm.flib.logsum(likefn(norms*s_pri+mu_pri))-np.log(len(norms))
        log_norm_const = mean_like-meanl(like_m,like_v,mu_pri,v_pri)
        return log_norm_const

def mean_reduce_with_hdf(hf, n_reps):
    """Produces an accumulator to be used with hdf5_to_samps"""    
    def mean_reduce_(sofar, next, name, ind, hf=hf, n_reps=n_reps):
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
    """Produces an accumulator to be used with hdf5_to_samps"""
    def var_reduce_(sofar, next, name, ind, hf=hf, n_reps=n_reps):
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

def find_joint_approx_params(mu_pri, C_pri, likefns, match_moments, approx_param_fn = None, tol=1.e-3, maxiter=100, debug=False):

    if approx_param_fn is not None:
        raise NotImplementedError
        
    norms = np.random.normal(size=1000)

    delta_m = np.ones_like(mu_pri)*np.inf
    delta_v = np.ones_like(mu_pri)*np.inf

    like_means = np.zeros_like(mu_pri)
    like_vars = np.ones_like(mu_pri)*np.inf
    
    norm_consts = np.ones_like(mu_pri)

    mu = mu_pri.copy('F')
    C = C_pri.copy('F')
    
    # Estimate the probability of this simulated dataset given this MCMC sample.
    # If it is radically improbable, the importance weight will be set to zero
    # and the algorithm for finding the exponentiated quadratic approximation
    # will be skipped.
    init_evidences = np.empty_like(mu_pri)
    init_like_means = np.empty_like(mu_pri)
    init_like_vars = np.empty_like(mu_pri)
    for i in xrange(len(mu_pri)):
        init_like_means[i], init_like_vars[i] = find_approx_params(mu[i], C[i,i], likefns[i], norms, match_moments, debug=debug)
        init_evidences[i] = pm.flib.logsum(likefns[i](mu[i]+C[i,i]*norms))-np.log(len(mu[i]))

    if init_evidences.min()>-20:
    
        for iter in xrange(maxiter):
            # Break the loop at convergence.
            if (np.all(np.abs(delta_m)<tol) and np.any(np.abs(delta_v/like_vars)<tol)):
                break
            # if iter % 200 == 0:
            #     print '200 iterations, dropping you into debugger next time.'
            #     debug = True
            for i in xrange(len(mu_pri)):
                
                # 'Unobserve' simulated datapoint i.
                if not np.isinf(like_vars[i]):
                    obsc(mu, C, like_means[i], -like_vars[i], i)

                if np.any(np.diag(C)<0):
                    warnings.warn('Negative element in diagonal of C. Assuming nonconvergence and returning early')
                    return init_like_means, init_like_vars, -np.inf
            
                # Find the exponentiated quadratic approximation for datapoint i.
                if approx_param_fn is None:
                    new_like_mean, new_like_var = find_approx_params(mu[i], C[i,i], likefns[i], norms, match_moments, debug=debug)
                else:
                    new_like_mean, new_like_var = approx_param_fn(mu[i], C[i,i], likefns[i])
                
                if np.isnan(new_like_var) or np.isnan(new_like_mean):
                    warnings.warn('Nan in like mean or var. Assuming nonconvergence and returning early.')
                    return init_like_means, init_like_vars, -np.inf
                
                delta_m[i] = new_like_mean-like_means[i]
                delta_v[i] = new_like_var- like_vars[i]
                like_means[i] = new_like_mean
                like_vars[i] = new_like_var
            
                # Compute the normalizing constant of the exponentiated quadratic approximation.
                norm_consts[i] = calc_norm_const(norms, like_means[i], like_vars[i], mu[i], C[i,i], likefns[i])
            
                # 'Re-observe' simulated datapoint i.
                if not np.isinf(like_vars[i]):
                    obsc(mu, C, like_means[i], like_vars[i], i)

        # After maximum number of iterations, check that the approximation is 'good enough.'
        if iter==maxiter:        
            klds = []
            for i in xrange(len(mu_pri)):
                obsc(mu,C,like_means[i],-like_vars[i],i)
                klds.append(kld2(like_means[i], like_vars[i], mu[i], C[i,i], likefns[i], norms))
                obsc(mu,C,like_means[i],like_vars[i],i)
            # A KL divergence of 0.1 is about OK.
            max_kld = np.max(klds)
            if max_kld > 0.2:
                warnings.warn('Maximum iterations used. Maximum KL divergence %f. Assuming nonconvergence and returning early.'%np.max(klds))
                return init_like_means, init_like_vars, -np.inf
                    
        log_imp_weight = impw(mu_pri,C_pri,like_means,like_vars,mu,C) + np.sum(norm_consts)    
    else:
        # This simulated observation is radically improbable given
        # this MCMC sample, so assign it an importance weight of zero.
        # It will have no role to play in predictions, and the algorithm
        # for finding the exponentiated quadratic approximation is likely
        # to fail, so don't bother.
        warnings.warn('Evidence very low (%f), returning early.'%init_evidences.min())
        return init_like_means, init_like_vars, -np.inf
    
    if np.isnan(log_imp_weight):
        raise RuntimeError
    
    return like_means, like_vars, log_imp_weight


def kloop_init(iter, k, M, x, pred_covariate_dict, survey_x, survey_covariate_dict, time_count, time_start):
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

     # Print a message about how much remains to be done.
     return time_msg(time_count, k, iter, time_start)


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
    like_means = np.zeros((len(iter),)+survey_data.shape)
    like_vars = np.zeros((len(iter),)+survey_data.shape)
    log_imp_weights = np.zeros((len(iter), len(survey_data)))
    
    print 'Looping through and getting the exponentiated quadratic approximations.'
    for k in xrange(len(iter)):
        
        time_count = kloop_init(iter, k, M, x, pred_covariate_dict, survey_x, survey_covariate_dict, time_count, time_start)

        base_M_preds = {}
        base_S_preds = {}
        # Draw some normal variates to use 
        # norms = dict([(s, scipy.stats.norm.ppf(np.linspace(1e-3,1-1e-3,1000))) for s in gp_submods])
        norms = dict([(s, np.random.normal(size=1000)) for s in gp_submods])

        # Accumulate for the 'current' maps.
        for s in gp_submods:
            M_obs[s] = pm.utils.value(s.M_obs)
            C_obs[s] = pm.utils.value(s.C_obs)
            nugs[s] = pm.utils.value(nuggets[s])
            base_M_preds[s], base_V_pred = pm.gp.point_eval(M_obs[s], C_obs[s], x)
            base_S_preds[s] = np.sqrt(base_V_pred + nugs[s])
            
        # The 'ind=-1' means 'given current data only'.
        apply_postprocs_and_reduce(M, n_per, base_M_preds, base_S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=-1, norms=norms[s])

        try:
            
            for l in xrange(len(survey_data)):

                M_preds = {}
                S_preds = {}

                for s in gp_submods:

                    # Assimilate this simulated dataset.
                    mu_pri = M_obs[s](survey_x)
                    C_pri = C_obs[s](survey_x, survey_x)+nugs[s]*np.eye(len(survey_x))
                    
                    closure_dict = {'data': survey_data[l], 'survey_plan': survey_plan}
                    
                    for extra_arg in extra_sl_args:
                        closure_dict[extra_arg] = pm.utils.value(getattr(M, extra_arg))
                    optim_fns = []
                    for ii in xrange(survey_data.shape[1]):
                        closure_dict_ = {'i':ii}
                        closure_dict_.update(closure_dict)
                        optim_fns.append(close(survey_likelihood,**closure_dict_))

                    # The importance weight expresses how probable simulated dataset l is
                    # given MCMC iteration k.
                    like_means[k,l], like_vars[k,l], log_imp_weights[k,l] = find_joint_approx_params(mu_pri, C_pri, optim_fns, match_moments)
        
        except np.linalg.LinAlgError:
            continue

    for l in xrange(len(survey_data)):
        # Normalize importance weights
        log_imp_weights[:,l] -= pm.flib.logsum(log_imp_weights[:,l])

    print 'Having resampled according to importance weights, producing maps conditional on execution of survey plan.'
    time_start = time.time()
    for k in xrange(len(iter)):        
        time_count = kloop_init(iter, k, M, x, pred_covariate_dict, survey_x, survey_covariate_dict, time_count, time_start)
        
        norms = dict([(s, np.random.normal(size=1000)) for s in gp_submods])

        # Accumulate for the 'current' maps.
        for s in gp_submods:
            M_obs[s] = pm.utils.value(s.M_obs)
            C_obs[s] = pm.utils.value(s.C_obs)
            nugs[s] = pm.utils.value(nuggets[s])

        try:
            for l in xrange(len(survey_data)):

                M_preds = {}
                S_preds = {}

                for s in gp_submods:

                    M_obs_ = copy.copy(M_obs[s])
                    C_obs_ = copy.copy(C_obs[s])
                    
                    # Accumulate for the 'conditional' maps.
                    pm.gp.observe(M_obs_, C_obs_, obs_mesh=survey_x, obs_vals = like_means[k,l], obs_V = like_vars[k,l])

                    M_preds[s], V_pred = pm.gp.point_eval(M_obs_, C_obs_, x)

                    if np.any(V_pred<0):
                        if continue_past_npd:
                            warnings.warn('Some elements of V_pred were negative. Assuming non-positive definiteness but not checking yet.')
                            actual_total -= n_per
                            log_imp_weights[k,:] = -np.inf
                            raise np.linalg.LinAlgError
                    S_preds[s] = np.sqrt(V_pred + nugs[s])
                    
                # The 'ind=l' means 'given simulated dataset l'.
                apply_postprocs_and_reduce(M, n_per, M_preds, S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=l, norms=norms[s])

        
        except np.linalg.LinAlgError:
            continue
    

    return actual_total, log_imp_weights
