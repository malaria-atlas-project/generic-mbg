from prediction_utils import *
from inference_utils import close

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
    
def find_approx_params(mu_pri, V_pri, likefn, norms):
    "Returns the 'likelihood' mean and variance minimizing K-L divergence with the likelihood."
    f = close(kldiv, mu_pri=mu_pri, V_pri=V_pri, likefn=likefn, norms=norms)

    from scipy import optimize
    p = optimize.fmin(f, [mu_pri, V_pri], disp=0)

    return p
    
def obs_corrections(mu_pri, C_pri, like_m, like_v, i):
    C_pri_scale = C_pri[i,:]/(C_pri[i,i]+like_v)
    C_corr = - np.outer(C_pri[i,:], C_pri_scale)
    mu_corr = (like_m-mu_pri[i])*C_pri_scale
    return mu_corr, C_corr

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

def find_joint_approx_params(mu_pri, C_pri, likefns, tol=1.e-3):
    norms = np.random.normal(size=1000)

    delta_m = np.ones_like(mu_pri)*np.inf
    delta_v = np.ones_like(mu_pri)*np.inf

    like_means = np.zeros_like(mu_pri)
    like_vars = np.ones_like(mu_pri)*np.inf

    mu = mu_pri.copy()
    C = C_pri.copy()

    mu_corrs = np.zeros(shape=(mu_pri.shape[0],)*2)
    C_corrs = np.zeros(shape=(mu_pri.shape[0],)*3)

    while np.any(np.abs(delta_m)>tol) or np.any(np.abs(delta_v/like_vars)>tol):
        for i in xrange(len(mu_pri)):
            mu -= mu_corrs[i]
            C -= C_corrs[i]

            new_like_mean, new_like_var = find_approx_params(mu[i], C[i,i], likefns[i], norms)

            delta_m[i] = new_like_mean-like_means[i]
            delta_v[i] = new_like_var- like_vars[i]
            like_means[i] = new_like_mean
            like_vars[i] = new_like_var

            mu_corrs[i], C_corrs[i] = obs_corrections(mu, C, like_means[i], like_vars[i], i)

            mu += mu_corrs[i]
            C += C_corrs[i]

    return like_means, like_vars, mu, C

def hdf5_to_survey_eval(M, x, nuggets, burn, thin, total, fns, postprocs, pred_covariate_dict, survey_x, survey_data, survey_covariate_dict, survey_likelihood, survey_plan, finalize=None, continue_past_npd=False):
    
    hf=M.db._h5file
    gp_submods = list(set(filter(lambda c: isinstance(c,pm.gp.GPSubmodel), M.containers)))
    f_labels = [gps.name for gps in gp_submods]

    products, postproc_args, extra_postproc_args = get_args(postprocs, fns, f_labels, M)
        
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
                print 'expect results '+time.ctime((time_count-time_start)*len(iter)/float(k)+time_start)
            else:
                print

        base_M_preds = {}
        base_S_preds = {}
        for s in gp_submods:
            M_obs[s] = pm.utils.value(s.M_obs)
            C_obs[s] = pm.utils.value(s.C_obs)
            nugs[s] = pm.utils.value(nuggets[s])
            base_M_preds[s], base_V_pred = pm.gp.point_eval(M_obs[s], C_obs[s], x)
            base_S_preds[s] = np.sqrt(base_V_pred + nugs[s])
            
            apply_postprocs_and_reduce(n_per, base_M_preds, base_S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=-1)
        
        try:
            for l in xrange(len(survey_data)):    
                M_preds = {}
                S_preds = {}
                for s in gp_submods:
                    mu_pri = M_obs[s](survey_x)
                    C_pri = C_obs[s](survey_x, survey_x)+nugs[s]*np.eye(len(survey_x))
                    
                    
                    # FIXME: This will only work for single fields currently.
                    like_means, like_vars, mu_post, C_post = find_joint_approx_params(mu_pri, C_pri, [close(survey_likelihood, **{'data': survey_data[l], 'survey_plan': survey_plan, 'i': ii}) for ii in xrange(survey_data.shape[1])])
                
                    M_obs_ = copy.copy(M_obs[s])
                    C_obs_ = copy.copy(C_obs[s])

                    pm.gp.observe(M_obs_, C_obs_, obs_mesh=survey_x, obs_vals = like_means, obs_V = like_vars)

                    M_preds[s], V_pred = pm.gp.point_eval(M_obs_, C_obs_, x)

                    if np.any(V_pred<0):
                        if continue_past_npd:
                            warnings.warn('Some elements of V_pred were negative. Assuming non-positive definiteness but not checking yet.')
                            actual_total -= n_per
                            raise np.linalg.LinAlgError
                    S_preds[s] = np.sqrt(V_pred + nugs[s])
                    
                    apply_postprocs_and_reduce(n_per, M_preds, S_preds, postprocs, fns, products, postproc_args, extra_postproc_args, joint=False, ind=l)

        except np.linalg.LinAlgError:
            continue

    return actual_total