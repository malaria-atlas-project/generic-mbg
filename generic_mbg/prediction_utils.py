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
import tables as tb
import numpy as np
from map_utils import asc_to_ndarray, get_header, exportAscii
from scipy import ndimage, mgrid
from histogram_utils import *
from inference_utils import all_chain_len, all_chain_getitem
import time
import os

__all__ = ['grid_convert','mean_reduce','var_reduce','invlogit','hdf5_to_samps','vec_to_asc','asc_to_locs',
            'display_asc','display_datapoints','histogram_reduce','histogram_finalize','maybe_convert','sample_reduce',
            'sample_finalize','asc_to_vals']
            
memmax = 2.5e8

def thread_partition_array(x):
    "Partition work arrays for multithreaded addition and multiplication"
    n_threads = int(os.environ['OMP_NUM_THREADS'])
    if len(x.shape)>1:
        maxind = x.shape[1]
    else:
        maxind = x.shape[0]
    bounds = np.array(np.linspace(0, maxind, n_threads+1),dtype='int')
    cmin = bounds[:-1]
    cmax = bounds[1:]
    return cmin,cmax

def invlogit(x):
    """A shape-preserving, in-place, threaded inverse logit function."""
    if np.prod(np.shape(x))<10000:
        return pm.flib.invlogit(x)
    if not x.flags['F_CONTIGUOUS']:
        raise ValueError, 'x is not Fortran-contiguous'
    cmin, cmax = thread_partition_array(x)        
    pm.map_noreturn(iinvlogit, [(x,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return x

def fast_inplace_square(a):
    cmin, cmax = thread_partition_array(a)
    pm.map_noreturn(iasq, [(a,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return a
    
def square_and_sum(a,s):
    cmin, cmax = thread_partition_array(a)
    pm.map_noreturn(asqs, [(a,s,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return a
    
def crossmul_and_sum(c,x,d,y):
    cmin, cmax = thread_partition_array(y)
    pm.map_noreturn(icsum, [(c,x,d,y,cmin[i],cmax[i]) for i in xrange(len(cmax))])
    return c

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

def validate_format_str(st):
    for i in [0,2]:
        if not st[i] in ['x','y']:
            raise ValueError, 'Directions must be x or y'
    for j in [1,3]:
        if not st[j] in ['-', '+']:
            raise ValueError, 'Orders must be + or -'
            
    if st[0]==st[2]:
        raise ValueError, 'Directions must be different'
    
    
def grid_convert(g, frm, to, validate=False):
    """Converts a grid to a new layout.
      - g : 2d array
      - frm : format string
      - to : format string
      
      Example format strings:
        - x+y+ (the way Anand does it) means that 
            - g[i+1,j] is west of g[i,j]
            - g[i,j+1] is north of g[i,j]
        - y-x+ (map view) means that 
            - g[i+1,j] is south of g[i,j]
            - g[i,j+1] is west of g[i,j]"""
    
    # Validate format strings
    if validate:
        for st in [frm, to]:
            validate_format_str(st)
        
    # Transpose if necessary
    if not frm[0]==to[0]:
        g = g.T
                
    first_dir = to[1]
    if not first_dir == frm[frm.find(to[0])+1]:
        g=g[::-1,:]
        
    sec_dir = to[3]
    if not sec_dir == frm[frm.find(to[2])+1]:
        g=g[:,::-1]
        
    # print first_dir, sec_dir
    return g

def buffer(arr, n=5):
    """Creates an n-pixel buffer in all directions."""
    out = arr.copy()
    if n==0:
        return out
    for i in xrange(1,n+1):
        out[i:,:] += arr[:-i,:]
        out[:-i,:] += arr[i:,:]
        out[:,i:] += arr[:,:-i]
        out[:,:-i] += arr[:,i:]
    return out

def asc_to_locs(fname, path='', thin=1, bufsize=1):
    """Converts an ascii grid to a list of locations where prediction is desired."""
    lon,lat,data = asc_to_ndarray(fname,path)
    data = grid_convert(data,'y-x+','x+y+')
    unmasked = buffer(True-data[::thin,::thin].mask, n=bufsize)
    
    # unmasked = None
    lat,lon = [dir[unmasked] for dir in np.meshgrid(lat[::thin],lon[::thin])]
    # lon,lat = [dir.ravel() for dir in np.meshgrid(lon[::thin],lat[::thin])]
    return np.vstack((lon,lat)).T*np.pi/180., unmasked

def asc_to_vals(fname, path='', thin=1, unmasked=None):
    """
    Converts an ascii grid to a list of values where prediction is desired.
    If the unmasked argument is provided, the mask of the ascii grid is
    checked against that mask.
    """
    lon,lat,data = asc_to_ndarray(fname, path)
    data = grid_convert(data,'y-x+','x+y+')    
    if unmasked is not None:
        input_unmasked = True-data[::thin,::thin].mask
        if not (unmasked == input_unmasked).all():
            where_mismatch = np.where(input_unmasked != unmasked)
            import pylab as pl
            pl.clf()
            pl.plot(lon[where_mismatch[0]],lat[where_mismatch[1]],'k.',markersize=2)
            pl.savefig('mismatch.pdf')
            msg = '%s: mask does not match input mask at the following pixels (in decimal degrees):\n'%fname
            for i,j in zip(*where_mismatch):
                msg += "\t%f, %f\n"%(lon[i],lat[j])
            msg += 'Image of mismatched points saved as mismatch.pdf'
            raise ValueError, msg
    
    return data.data[unmasked]
    
def display_asc(fname, path='', radians=True, *args, **kwargs):
    """Displays an ascii file as a pylab image."""
    from pylab import imshow
    lon,lat,data = asc_to_ndarray(fname,path)
    if radians:
        lon *= np.pi/180.
        lat *= np.pi/180.
    imshow(grid_convert(data,'y-x+','y+x+'),extent=[lon.min(),lon.max(),lat.min(),lat.max()],*args,**kwargs)
    
def display_datapoints(h5file, path='', cmap=None, *args, **kwargs):
    """Adds as hdf5 archive's logp-mesh to an image."""
    hf = tb.openFile(path+h5file)
    lpm = hf.root.metadata.logp_mesh[:]
    hf.close()
    from pylab import plot
    if cmap is None:
        plot(lpm[:,0],lpm[:,1],*args,**kwargs)
        
def mean_reduce(sofar, next):
    """A function to be used with hdf5_to_samps"""
    if sofar is None:
        return next
    else:
        return sofar + next
        
def var_reduce(sofar, next):
    """A function to be used with hdf5_to_samps"""
    if sofar is None:
        return next**2
    else:
        return sofar + next**2
        
def moments_finalize(prod, n):
    """Finalizes accumulated moments in human-interpreted surfaces."""
    mean = prod[mean_reduce] / n
    var = prod[var_reduce] / n - mean**2
    std = np.sqrt(var)
    std_to_mean = std/mean
    out = {'mean': mean, 'var': var, 'std': std, 'std-to-mean':std_to_mean}
    return out
        
def sample_reduce(sofar, next):
    """A function to be used with hdf5_to_samps. Keeps all samples with no data loss."""
    if sofar is None:
        return [next]
    else:
        sofar.append(next)
        return sofar
        
def sample_finalize(prod, n):
    """Converts accumulated samples to an array for output."""
    return np.array(prod[sample_reduce])

def histogram_reduce(bins, binfn):
    """Produces an accumulator to be used with hdf5_to_samps"""
    def hr(sofar, next):
        if sofar is None:
            sofar = np.zeros(next.shape+(len(bins),), dtype=int, order='F')
        # Call to Fortran function multiinc
        ind = binfn(next)
        multiinc(sofar,ind)
        
        return sofar
    return hr
        
def histogram_finalize(bins, q, hr):
    """Converts accumulated histogram raster to desired quantiles"""
    def fin(products, n, bins=bins, q=q, hr=hr):
        out = {}
        hist = products[hr]
        # Call to Fortran function qextract
        quantile_surfs = qextract(hist,n,q,bins)
        for i in xrange(len(q)):
            out['quantile-%s'%q[i]] = quantile_surfs[i]

        return out
    return fin



def hdf5_to_samps(hf, x, burn, thin, total, fns, f_label, f_has_nugget, x_label, pred_cv_dict=None, nugget_label=None, postproc=None, finalize=None, diag_safe=False, **non_cov_columns):
    """
    Parameters:
        hf : PyTables file
            The trace file from which predictions should be made.
        x : array
            The lon, lat locations at which predictions are desired.
        burn : int
            Burnin iterations to discard.
        thin : int
            Number of iterations between ones that get used in the predictions.
        total : int
            Total number of iterations to use in thinning.
        fns : list of functions
            Each function should take four arguments: sofar, next, cols and i.
            Sofar may be None.
            The functions will be applied according to the reduce pattern.
        f_label : string
            The name of the hdf5 node containing f
        f_has_nugget : boolean
            Whether f is nuggeted.
        x_label : string
            The name of the hdf5 node containing the input mesh associated with f
            in the metadata.
        pred_cv_dict : dictionary
            {name : value on x}
        nugget_label : string (optional)
            The name of the hdf5 node giving the nugget variance
        postproc : function (optional)
            This function is applied to the realization before it is passed to
            the fns.
        finalize : function (optional)
            This function is applied to the products before returning. It should
            take a second argument which is the actual number of realizations
            produced.
    """
    
    # Add constant mean
    if pred_cv_dict is None:
        pred_cv_dict = {}
    pred_cv_dict['m'] = np.ones(x.shape[0])
    
    products = dict(zip(fns, [None]*len(fns)))
    iter = np.arange(burn,all_chain_len(hf),thin)
    if len(iter)==0:
        raise ValueError, 'You asked for %i burnin iterations with thinnnig %i but the chains are only %i iterations long.'%(burn, thin, all_chain_len(hf))
    n_per = total/len(iter)+1
    actual_total = n_per * len(iter)
    
    x_obs = getattr(hf.root.metadata,x_label)[:]
    
    # Avoid memory errors
    # max_chunksize = 1.e8 / x_obs.shape[0]
    # n_chunks = int(x.shape[0]/max_chunksize+1)
    # splits = np.array(np.linspace(0,x.shape[0],n_chunks+1)[1:-1],dtype='int')
    # x_chunks = np.split(x,splits)
    # i_chunks = np.split(np.arange(x.shape[0]), splits)
    
    # If postproc is not None, close on non-covariate columns.
    if postproc is not None:
        if len(non_cov_columns) > 0:
            postproc = postproc(**non_cov_columns)
    
    time_count = -np.inf
    time_start = time.time()

    for k in xrange(len(iter)):
        
        i = iter[k]
        
        if time.time() - time_count > 10:
            print ((k*100)/len(iter)), '% complete',
            time_count = time.time()      
            if k > 0:      
                print 'expect results '+time.ctime((time_count-time_start)*len(iter)/float(k)+time_start)
            else:
                print
        
        M_pred, S_pred = predictive_mean_and_std(hf, i, f_label, x_label, x, f_has_nugget, pred_cv_dict, nugget_label, diag_safe)
        cmin, cmax = thread_partition_array(M_pred)
        
        # Postprocess if necessary: logit, etc.
        norms = np.random.normal(size=n_per)
        
        for j in xrange(n_per):
            # surf = M_pred
            # surf = M_pred + S_pred * norms[j]

            surf = M_pred.copy('F')
            pm.map_noreturn(iaaxpy, [(norms[j], S_pred, surf, cmin[l], cmax[l]) for l in xrange(len(cmax))])
            
            if postproc is not None:
                surf = postproc(surf)
                        
            # Reduction step
            for f in fns:
                products[f] = f(products[f], surf)

    
    if finalize is not None:
        return finalize(products, actual_total)
    else:          
        return products
        
def normalize_for_mapcoords(arr, max):
    "Used to create inputs to ndimage.map_coordinates."
    arr /= arr.max()
    arr *= max
    
def vec_to_asc(vec, fname, out_fname, unmasked, path=''):
    """
    Converts a vector of outputs on a thin, unmasked, ravelled subset of an
    ascii grid to an ascii file matching the original grid.
    """
    
    if np.any(np.isnan(vec)):
        raise ValueError, 'NaN in vec'  
    
    header, headlines = get_header(fname,path)
    lon,lat,data = asc_to_ndarray(fname,path)
    data = grid_convert(data,'y-x+','x+y+')
    data_thin = np.zeros(unmasked.shape)
    data_thin[unmasked] = vec
    
    mapgrid = np.array(mgrid[0:data.shape[0],0:data.shape[1]], dtype=float)
    normalize_for_mapcoords(mapgrid[0], data_thin.shape[0]-1)
    normalize_for_mapcoords(mapgrid[1], data_thin.shape[1]-1)
    
    if data_thin.shape == data.shape:
        out = np.ma.masked_array(data_thin, mask=data.mask)
    else:
        out = np.ma.masked_array(ndimage.map_coordinates(data_thin, mapgrid), mask=data.mask)
        
    if np.any(np.isnan(out)):
        raise ValueError, 'NaN in output'
    
    # import pylab as pl
    # pl.close('all')
    # pl.figure()
    # pl.imshow(out, interpolation='nearest', vmin=0.)
    # pl.colorbar()
    # pl.title('Resampled')
    # 
    # pl.figure()
    # pl.imshow(np.ma.masked_array(data_thin, mask=True-unmasked), interpolation='nearest', vmin=0)
    # pl.colorbar()
    # pl.title('Original')
    
    out_conv = grid_convert(out,'x+y+','y-x+')
    
    header['NODATA_value'] = -9999
    exportAscii(out_conv.data,out_fname,header,True-out_conv.mask)
    
    return out
    

def predictive_mean_and_std(hf, i, f_label, x_label, x, f_has_nugget=False, pred_cv_dict=None, nugget_label=None, diag_safe=False):
    """
    Computes marginal (pointwise) predictive mean and variance for f(x).
    Expects input from an hdf5 datafile.
        - hf : hdf5 file.
        - i : integer.
        - f_label : string or array.
        - x : numpy array
        - pred_cv_dict : {name : value-on-predmesh}
        - nugget_label : string
    """

    meta = hf.root.metadata

    if pred_cv_dict is None:
        raise ValueError, 'No pred_cv_dict provided. You always have the constant term. Tell Anand.'
            
    n = pred_cv_dict.keys()

    covariate_dict = meta.covariates[0]
    prior_covariate_variance = np.diag([covariate_dict[key][1] for key in n])

    pred_covariate_values = np.empty((len(pred_cv_dict), x.shape[0]), order='F')
    input_covariate_values = np.empty((len(pred_cv_dict), len(covariate_dict[key][0])), order='F')
    for j in xrange(len(n)):
        pred_covariate_values[j,:] = pred_cv_dict[n[j]]
        input_covariate_values[j,:] = covariate_dict[n[j]][0]
    
    # How many times must a man condition a multivariate normal
    M = all_chain_getitem(hf, 'M', i, True)
    C = all_chain_getitem(hf, 'C', i, True)

    logp_mesh = np.asarray(getattr(meta,x_label)[:], order='F')    
    M_input = M(logp_mesh)
    M_pred = M(x)
    
    V_out = np.empty(M_pred.shape)
    M_out = np.empty(M_pred.shape)
    
    try:
        f = all_chain_getitem(hf, f_label, i)
    except:
        f = getattr(meta,f_label)[:]
        
    C_input = C(logp_mesh, logp_mesh)

    if pred_cv_dict is not None:
        C_input += np.dot(np.dot(input_covariate_values.T, prior_covariate_variance), input_covariate_values)
    if nugget_label is not None:
        nug = all_chain_getitem(hf, 'V', i)
        if f_has_nugget:
            C_input += nug*np.eye(np.sum(logp_mesh.shape[:-1]))
    else:
        nug = 0.
    
    try:
        S_input = np.linalg.cholesky(C_input)
        piv = None
    except np.linalg.LinAlgError:
        U, rank, piv = pm.gp.incomplete_chol.ichol_full(c=C_input, reltol=1.e-10)
        print 'Warning, full conditional covariance was not positive definite at iteration %i. Rank is %i of %i.'%(i, rank, C_input.shape[0])        
        if rank<=0:
            raise ValueError, "Matrix does not appear to be positive semidefinite. Tell Anand."
        else:
            S_input = np.asarray(U[:rank,:rank].T, order='F')
            logp_mesh = logp_mesh[piv[:rank]]
            f = f[piv[:rank]]
            M_input = M_input[piv[:rank]]
            input_covariate_values = np.asarray(input_covariate_values[:,piv[:rank]], order='F')
                
    max_chunksize = memmax / 8 / logp_mesh.shape[0]
    n_chunks = int(x.shape[0]/max_chunksize+1)
    splits = np.array(np.linspace(0,x.shape[0],n_chunks+1),dtype='int')
    x_chunks = np.split(x,splits[1:-1])
    i_chunks = [slice(splits[i],splits[i+1],None) for i in xrange(n_chunks)]


    for k in xrange(n_chunks):

        i_chunk = i_chunks[k]
        x_chunk = x_chunks[k]
        
        pcv = pred_covariate_values[:,i_chunk]
        C_cross = C(logp_mesh, x_chunk) 
        if diag_safe:
            if np.any(C_cross>C.params['amp']**2):
                raise ValueError, 'Off-diagonal elements of C are too large. Tell Anand.'
        
        for mat in [input_covariate_values, pcv, C_cross, S_input]:
            if not mat.flags['F_CONTIGUOUS']:
                raise ValueError, 'Matrix is not Fortran-contiguous. Tell Anand.'
        
        if diag_safe:
            V_pred = C.params['amp']**2 + nug
        else:
            V_pred = C(x_chunk) + nug
        
        
        if pred_cv_dict is not None:
            C_cross = crossmul_and_sum(C_cross, input_covariate_values, np.diag(prior_covariate_variance), pcv)
            V_pred_adj = V_pred + np.sum(np.dot(np.sqrt(prior_covariate_variance), pcv)**2, axis=0)
                        
        SC_cross = pm.gp.trisolve(S_input,C_cross,uplo='L',inplace=True)

        for mat in S_input, C_cross, SC_cross:
            if not mat.flags['F_CONTIGUOUS']:
                raise ValueError, 'Matrix is not Fortran-contiguous. Tell Anand.'
        
        scc = np.empty(x_chunk.shape[0])
        square_and_sum(SC_cross, scc)
        
        V_out[i_chunk] = V_pred_adj - scc
        M_out[i_chunk] = M_pred[i_chunk] + np.asarray(np.dot(SC_cross.T,pm.gp.trisolve(S_input, (f-M_input), uplo='L'))).squeeze()


    if np.any(np.isnan(np.sqrt(V_out))) or np.any(np.isnan(M_out)):
        raise ValueError, 'Some predictive samples were NaN. Keep all your input files and tell Anand.'

    return M_out, np.sqrt(V_out)