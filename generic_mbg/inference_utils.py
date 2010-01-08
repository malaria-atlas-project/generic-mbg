# Copyright (C) 2009 Anand Patil
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

import pymc as pm
import numpy as np
import time
import tables as tb
#from st_cov_fun import my_st

def spatial_mean(x, m_const):
    return m_const*np.ones(x.shape[0])
    
def zero_mean(x):
    return np.zeros(x.shape[:-1])

def st_mean_comp(x, m_const, t_coef):
    lon = x[:,0]
    lat = x[:,1]
    t = x[:,2]
    return m_const + t_coef * t

def combine_spatial_inputs(lon,lat):
    # Convert latitude and longitude from degrees to radians.
    lon = lon*np.pi/180.
    lat = lat*np.pi/180.
    
    # Make lon, lat tuples.
    data_mesh = np.vstack((lon, lat)).T 
    return data_mesh
    
def combine_st_inputs(lon,lat,t):
    # Convert latitude and longitude from degrees to radians.
    lon = lon*np.pi/180.
    lat = lat*np.pi/180.

    # Make lon, lat, t triples.
    data_mesh = np.vstack((lon, lat, t)).T 
    return data_mesh

def chains(hf):
    return [gr for gr in hf.listNodes("/") if gr._v_name[:5]=='chain']

def all_chain_len(hf):
    return np.sum([len(chain.PyMCsamples) for chain in chains(hf)])
    
def all_chain_trace(hf, name):
    return np.concatenate([np.ravel(chain.PyMCsamples.col(name)) for chain in chains(hf)])
    
def all_chain_getitem(hf, name, i, vl=False):
    c = chains(hf)
    lens = [len(chain.PyMCsamples) for chain in c]
    if i >= np.sum(lens):
        raise IndexError, 'Index out of bounds'
    s = 0
    j = i

    for k in xrange(len(lens)):
        s += lens[k]
        if i<s:
            if vl:
                return getattr(c[k].group0, name)[j]
            else:
                return getattr(c[k].PyMCsamples.cols,name)[j]
        else:
            j -= lens[k]
    
def add_standard_metadata(M, x_labels, *others):
    """
    Adds the standard metadata to an hdf5 archive.
    """
    
    hf = M.db._h5file
    hf.createGroup('/','metadata')
        
    hf.createVLArray(hf.root.metadata, 'covariates', tb.ObjectAtom())
    hf.root.metadata.covariates.append(M.covariate_dicts)
    
    for label in set(x_labels.itervalues()):
        hf.createArray(hf.root.metadata, label, getattr(M, label))
        
    for label in set(others):
        vla=hf.createVLArray(hf.root.metadata, label, tb.ObjectAtom())
        vla.append(getattr(M,label))
        
        
    
def cd_and_C_eval(covariate_values, C, data_mesh, ui=slice(None,None,None), fac=1e6):
    """
    Returns a {name: value, prior variance} dictionary
    and an evaluated covariance with covariates incorporated.
    """
    covariate_dict = {}
    # Set prior variances of covariate coefficients. They're huge, and scaled.
    means = []
        
    for cname, cval in covariate_values.iteritems():
        cov_var = np.var(cval)
        cov_mean = np.abs(np.mean(cval))+np.sqrt(cov_var)
        means.append(cov_mean)
        covariate_dict[cname] = (cval, cov_var*fac)
        
    # Constant term
    covariate_dict['m'] = (np.ones(data_mesh.shape[0]), (np.sum(np.array(means)**2) + 1)*fac)
    logp_mesh = data_mesh[ui]
                    
    # The evaluation of the Covariance object, plus the nugget.
    @pm.deterministic(trace=False)
    def C_eval(C=C,ui=ui,cd=covariate_dict):
        out = C(logp_mesh, logp_mesh)
        for val,var in cd.itervalues():
            valu = val[ui]
            out += np.outer(valu,valu)*var
        return out
    
    return covariate_dict, C_eval
    
def trivial_means(lpm):
    """
    Returns a trivial mean function and an evaluating node.
    """
    # The mean of the field
    @pm.deterministic(trace=True)
    def M():
        return pm.gp.Mean(zero_mean)
    
    # The mean, evaluated  at the observation points, plus the covariates    
    @pm.deterministic(trace=False)
    def M_eval(M=M, lpm=lpm):
        return M(lpm)
    return M, M_eval

def basic_spatial_submodel(lon, lat, covariate_values):
    """
    A stem for building spatial models.
    """
    logp_mesh = combine_spatial_inputs(lon,lat)

    # =====================
    # = Create PyMC model =
    # =====================  

    inc = pm.CircVonMises('inc', 0,0)
    sqrt_ecc = pm.Uniform('sqrt_ecc',0,.95)
    ecc = pm.Lambda('ecc', lambda s=sqrt_ecc: s**2)
    amp = pm.Exponential('amp',.1,value=1.)
    scale_shift = pm.Exponential('scale_shift',1./.08,value=1./.08)
    scale = pm.Lambda('scale',lambda ss=scale_shift:ss+.01)
    diff_degree = pm.Uniform('diff_degree',.01,3)
    
    M, M_eval = trivial_means(logp_mesh)

    @pm.deterministic(trace=True)
    def C(amp=amp,scale=scale,inc=inc,ecc=ecc,diff_degree=diff_degree):
        return pm.gp.FullRankCovariance(pm.gp.cov_funs.matern.aniso_geo_rad, amp=amp, scale=scale, inc=inc, ecc=ecc, diff_degree=diff_degree)
    
    covariate_dict, C_eval = cd_and_C_eval(covariate_values, C, logp_mesh)
    
    return locals()


def basic_st_submodel(lon, lat, t, covariate_values, cpus):
    """
    A stem for building spatiotemporal models.
    """
        
    logp_mesh = combine_st_inputs(lon,lat,t)

    inc = pm.CircVonMises('inc', 0,0)
    sqrt_ecc = pm.Uniform('sqrt_ecc',0,.95)
    ecc = pm.Lambda('ecc', lambda s=sqrt_ecc: s**2)
    amp = pm.Exponential('amp',.1,value=1.)
    scale = pm.Exponential('scale',1.,value=1.)
    scale_t = pm.Exponential('scale_t',.1,value=.1)
    t_lim_corr = pm.Uniform('t_lim_corr',0,.95)

    @pm.stochastic(__class__ = pm.CircularStochastic, lo=0, hi=1)
    def sin_frac(value=.1):
        return 0.

    M, M_eval = trivial_means(logp_mesh)
        
    # A constraint on the space-time covariance parameters that ensures temporal correlations are 
    # always between -1 and 1.
    @pm.potential
    def st_constraint(sd=.5, sf=sin_frac, tlc=t_lim_corr):    
        if -sd >= 1./(-sf*(1-tlc)+tlc):
            return -np.Inf
        else:
            return 0.

    # A Deterministic valued as a Covariance object. Uses covariance my_st, defined above. 
    @pm.deterministic
    def C(amp=amp,scale=scale,inc=inc,ecc=ecc,scale_t=scale_t, t_lim_corr=t_lim_corr, sin_frac=sin_frac):
        return pm.gp.FullRankCovariance(my_st, amp=amp, scale=scale, inc=inc, ecc=ecc,st=scale_t, sd=.5,
                                        tlc=t_lim_corr, sf = sin_frac, n_threads=cpus)

    covariate_dict, C_eval = cd_and_C_eval(covariate_values, C, logp_mesh)
        
    return locals()

def sample_covariates(covariate_dict, C_eval, d):
    """
    Samples covariates back in when they have been marginalized away.
        - covariate_dict : {name : value-on-input, prior-variance}
        - M_eval : array. Probably zeros, unless you did something fancy in the mean.
        - C_eval : covariance of d | covariates, m
        - d : current deviation from mean of covariates' immediate child.
    """
    # Extract keys to list to preserve order.
    n = covariate_dict.keys()
    
    cvv = [covariate_dict[k] for k in n]
    x = np.asarray([v[0] for v in cvv])
    prior_var = np.diag([v[1] for v in cvv])
    
    prior_offdiag = np.dot(prior_var,x).T
    prior_S = np.linalg.cholesky(np.asarray(C_eval) + np.dot(prior_offdiag, x))
    pm.gp.trisolve(prior_S, prior_offdiag, uplo='L', transa='N', inplace=True)
    post_C = prior_var - np.dot(prior_offdiag.T, prior_offdiag)
    post_mean = np.dot(prior_offdiag.T, pm.gp.trisolve(prior_S, d, uplo='L', transa='N'))
    
    new_val = pm.rmv_normal_cov(post_mean, post_C).squeeze()

    return dict(zip(n, new_val))

def get_d_C_eval(hf, f_label, nugget_label, i, mesh):
    """Utility fn"""
    if type(f_label) == type('str'):
        d = all_chain_getitem(hf, f_label, i, vl=False)
    else:
        d = f_label

    C = all_chain_getitem(hf, 'C', i, vl=True)
    if nugget_label is not None:
        nug = all_chain_getitem(hf, nugget_label, i, vl=False)

    C_eval = C(mesh, mesh) + nug*np.eye(np.sum(mesh.shape[:-1]))
    return d, C_eval

def covariate_trace(hf, f_label, nugget_label=None, burn=0, thin=1):
    """
    Produces a covariate trace from an existing hdf5 chain.
        - chain : hdf5 group
        - meta : hdf5 group
        - f_label : string or array
        - nugget_label : string
    """
    meta = hf.root.metadata
    
    covariate_dict = meta.covariates[0]

    out = dict.fromkeys(covariate_dict)
    for k in out.keys():
        out[k] = []
        
    if nugget_label is None:
        mesh = meta.logp_mesh[:]
    else:
        mesh = meta.data_mesh[:]

    n = all_chain_len(hf)
    time_count = -np.inf
    time_start = time.time()
        
    for i in xrange(burn,n,thin):

        if time.time() - time_count > 10:
            print ((i*100)/n), '% complete',
            time_count = time.time()     
            if i > 0:       
                print 'expect results '+time.ctime((time_count-time_start)*n/float(i)+time_start)        
            else:
                print

        d, C_eval = get_d_C_eval(hf, f_label, nugget_label, i, mesh)
        cur_vals = sample_covariates(covariate_dict, C_eval, d)

        for k, v in cur_vals.iteritems():
            out[k].append(v)

    return dict([(k,np.array(v)) for k,v in out.iteritems()])
