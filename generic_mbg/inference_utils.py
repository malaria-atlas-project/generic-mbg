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
from st_cov_fun import my_st
from histogram_utils import iinvlogit, iamul, iasq, icsum, subset_eq, iasadd
from pylab import csv2rec

def maybe_convert(ra, field, dtype):
    """
    Tries to cast given field of given record array to given dtype. 
    Raises helpful error on failure.
    """
    arr = ra[field]
    try:
        return arr.astype(dtype)
    except:
        for i in xrange(len(arr)):
            try:
                np.array(arr[i],dtype=dtype)
            except:
                raise ValueError, 'Input column %s, element %i (starting from zero) is %s,\n which cannot be cast to %s'%(field,i,arr[i],dtype)


def invlogit(x):
    """A shape-preserving, in-place, threaded inverse logit function."""
    if np.prod(np.shape(x))<10000:
        return pm.flib.invlogit(x)
    if not x.flags['F_CONTIGUOUS']:
        raise ValueError, 'x is not Fortran-contiguous'
    cmin, cmax = pm.thread_partition_array(x)        
    pm.map_noreturn(iinvlogit, [(x,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return x

def fast_inplace_mul(a,s):
    """Multiplies a by s in-place and returns a."""
    a = np.atleast_2d(a)
    s = np.atleast_2d(s)
    cmin, cmax = pm.thread_partition_array(a)
    pm.map_noreturn(iamul, [(a,s,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return a
    
def fast_inplace_scalar_add(a,s):
    """Adds s to a in-place and returns a. s should be a scalar."""
    a = np.atleast_2d(a)
    cmin, cmax = pm.thread_partition_array(a)
    pm.map_noreturn(iasadd, [(a,s,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return a

def fast_inplace_square(a):
    """Squares a in-place and returns it."""
    cmin, cmax = pm.thread_partition_array(a)
    pm.map_noreturn(iasq, [(a,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return a
    
def crossmul_and_sum(c,x,d,y):
    """Returns C + \sum_i d_i*outer(x[i,:],y[i,:])"""
    cmin, cmax = pm.thread_partition_array(y)
    pm.map_noreturn(icsum, [(c,x,d,y,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return c

def create_model(mod,db,input=None):
    """
    Takes either:
    - A module and a database object, or:
    - A module, a filename and a record array
    and returns an MCMC object, with step methods assigned.
    """
    
    if isinstance(db,pm.database.hdf5.Database):
        if input is not None:
            raise ValueError, 'Input provided with preexisting db.'
        input = db._h5file.root.input_csv[:]
        prev_db = db
    else:
        if input is None:
            raise ValueError, 'No input or preexisting db provided.'
        prev_db=None
        hfname = db
        
    lon = maybe_convert(input, 'lon', 'float')
    lat = maybe_convert(input, 'lat', 'float')
    mod_inputs = (lon,lat)
    if 't' in input.dtype.names:
        t = maybe_convert(input, 't', 'float')
        x = combine_st_inputs(lon,lat,t)
        mod_inputs = mod_inputs + (t,)
    else:
        x = combine_spatial_inputs(lon,lat)
    
    non_cov_columns = {}
    if hasattr(mod, 'non_cov_columns'):
        non_cov_coltypes = mod.non_cov_columns
    else:
        non_cov_coltypes = {}
    non_cov_colnames = non_cov_coltypes.keys()

    covariate_dict = {}
    for n in input.dtype.names:
        if n not in ['lon','lat','t']:
            if n in non_cov_colnames:
                non_cov_columns[n] = maybe_convert(input, n, non_cov_coltypes[n])
            else:
                covariate_dict[n]=maybe_convert(input, n, 'float')

    mod_inputs = mod_inputs + (covariate_dict,)

    # Create MCMC object, add metadata, and assign appropriate step method.

    if prev_db is None:
        M = pm.MCMC(mod.make_model(*mod_inputs,**non_cov_columns),db='hdf5',dbname=hfname,complevel=1,complib='zlib')
    else:
        M = pm.MCMC(mod.make_model(*mod_inputs,**non_cov_columns),db=prev_db)
        M.restore_sampler_state()
        
    # Pass MCMC object through the module's mcmc_init function.
    mod.mcmc_init(M)
    
    # Restore step method state if possible.
    M.assign_step_methods() 
    if prev_db is not None:
        M.restore_sm_state()

    return M

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

def add_standard_metadata(M, *labels):
    """
    Adds the standard metadata to an hdf5 archive.
    """
    
    hf = M.db._h5file
    hf.createGroup('/','metadata')
        
    for label in set(labels):
        vla=hf.createVLArray(hf.root.metadata, label, tb.ObjectAtom())
        vla.append(getattr(M,label))
        
class CachingCovariateEvaluator(object):
    """
    Evaluate this object on an array. If it has already been
    told what its value on the array is, it will return that value.
    Otherwise, it will throw an error.
    """
    def __init__(self, mesh, value, shift, scale):
        self.meshes = [mesh]
        self.values = [value]
        if np.any(np.isnan(value)):
            raise ValueError, 'NaN in covariate values'
        self.shift = shift
        self.scale = scale
    def add_value_to_cache(self, mesh, value):
        if np.any(np.isnan(value)):
            raise ValueError, 'NaN in covariate values'
        # elif len(set(value))==1:
        #     raise ValueError, 'Only one value in covariate. Not going to learn much here.'
        self.meshes.append(mesh)
        self.values.append((value-self.shift)/self.scale)
    def __call__(self, mesh):
        for i,m in enumerate(self.meshes):
            start,stop = subset_eq(m,mesh)
            if start>-1 and stop>-1:
                if stop-start != mesh.shape[0]:
                    raise ValueError
                return self.values[i][start:stop]
        raise RuntimeError, 'The given mesh is not in the cache.'
        
class CovarianceWithCovariates(object):
    """
    A callable that adds some covariates to the given covariance function C.
    
    The output is:
        ~C(x,y) = C(x,y) + fac(m + \sum_{n\in names} outer(~cov[n](x), ~cov[n](y)))

    m is the number of covariates.
        
    ~cov[n](x) is (cov[n](x)-mu[n])/sig[n], where mu[n] and sig[n] are the mean and
    standard deviation of the values passed into the init method.
    """
    
                        
    def __init__(self, cov_fun, mesh, cv, fac=1e6, ampsq_is_diag=False):
        self.cv = cv
        self.labels = self.cv.keys()
        self.m = len(cv)
        self.means = dict([(k,np.mean(v)) for k,v in cv.iteritems()])
        self.stds = dict([(k,np.std(v) or 1) for k,v in cv.iteritems()])
        self.evaluators = dict([(k,CachingCovariateEvaluator(mesh[:,:2], v, self.means[k], self.stds[k])) for k,v in cv.iteritems()])
        self.cov_fun = cov_fun
        self.fac = fac
        self.privar = np.ones(len(self.labels))*self.fac

    def diag_call(self,x, *args, **kwds):
        # Evaluate with one argument:
        x_evals = self.eval_covariates(x)        
        if hasattr(self.cov_fun, 'diag_call'):
            Cbase = self.cov_fun.diag_call(x,*args,**kwds)
        else:
            Cbase = self.cov_fun(x,y=None,*args,**kwds)
        if len(self.labels)>0:
            C = Cbase + np.sum(self.privar * x_evals.T**2, axis=1) + self.fac*self.m
        else:
            C = Cbase + self.fac*self.m
        return C
        
    def eval_covariates(self, x):
        out = np.asarray([self.evaluators[k](x[:,:2]) for k in self.labels], order='F')
        return out
        # return np.asarray([np.ones(len(x)) for k in self.labels], order='F')
        
    def add_values_to_cache(self, mesh, new_cv):
        for k,v in self.evaluators.iteritems():
            v.add_value_to_cache(mesh[:,:2], new_cv[k])

    def __call__(self, x, y, *args, **kwds):
        
        # Evaluate with both arguments:
        Cbase = self.cov_fun(x,y,*args,**kwds)
        
        if len(self.evaluators) > 0:
            x_evals = self.eval_covariates(x)
            if x is y:
                y_evals = x_evals
            else:
                y_evals = self.eval_covariates(y)
            C = crossmul_and_sum(Cbase, x_evals, self.privar, y_evals)
            C = fast_inplace_scalar_add(C, self.fac*self.m)
        else:
            C = Cbase
        return C

def sample_covariates(covariate_dict, C_eval, d):
    """
    Samples covariates back in when they have been marginalized away.
        - covariate_dict : {name : value-on-input, prior-variance}
        - M_eval : array. Probably zeros, unless you did something fancy in the mean.
        - C_eval : covariance of d | covariates, m
        - d : current deviation from mean of covariates' immediate child.
    """
    raise NotImplementedError, 'This has not been ported to PyMC 2.1 yet.'
    
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
