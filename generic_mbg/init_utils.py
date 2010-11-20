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

import imp, os, sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
matplotlib.interactive(False)
import pymc as pm
from pylab import csv2rec,rec2csv
import decimal
from map_utils import import_raster, export_raster

delayed_import_str ="""
import matplotlib
from map_utils import *
import tables as tb
import numpy as np
import pymc as pm
import pylab as pl"""

def get_nuggets_and_obs(mod, mod_name, M):
    try:
        nugget_labels = getattr(mod, 'nugget_labels')
        nuggets = dict([(getattr(M,k),getattr(M,v)) for k,v in mod.nugget_labels.iteritems()])
        obs_labels = getattr(mod, 'obs_labels')
    except:
        cls,inst,tb = sys.exc_info()
        new_inst = cls('Could not import nugget_labels or obs_labels from %s. Tell Anand. Original error message:\n\n\t%s'%(mod_name,inst.message))
        raise cls,new_inst,tb
    return nuggets, obs_labels

def parse_quantiles(o):
    if len(o.quantile_list) == 0:
        q = []
    else:
        q = map(decimal.Decimal, o.quantile_list.split(' '))
    return q

def get_reduce_finalize(mod):
    if hasattr(mod, 'extra_reduce_fns'):
        extra_reduce_fns = mod.extra_reduce_fns
        extra_finalize = mod.extra_finalize
    else:
        extra_reduce_fns = []
        extra_finalize = None
    return extra_reduce_fns, extra_finalize
        
def init_output_dir(o, suffix):
    hf_path, hf_basename  = os.path.split(o.hf_name)
    base_outname = os.path.splitext(hf_basename)[0]
    output_dir = os.path.join(hf_path, base_outname+'-%s'%suffix)
    try:
        os.mkdir(output_dir)
    except OSError:
        pass
    os.chdir(output_dir) 
    
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

def get_circle(n):
    cover = np.arange(-n,n+1)
    coversq = np.dstack(np.meshgrid(cover,cover))
    covercirc = np.where(np.sum(coversq**2, axis=-1)<=n**2)
    out =  coversq[covercirc]
    return out

def buffer(arr, n=5):
    """Creates an n-pixel buffer in all directions."""
    arr = np.asarray(arr,order='F',dtype='bool')
    if n > 0:
        circ_ind = get_circle(n)
        out = bufster(arr, circ_ind)
    else:
        out = arr.copy('F')
    return out.astype('bool')

def raster_to_locs(name, path='.', thin=1, bufsize=1):
    """Converts a raster grid to a list of locations where prediction is desired."""
    lon,lat,data,type = import_raster(name,path)
    data = grid_convert(data,'y-x+','x+y+')
    unmasked = buffer(True-data[::thin,::thin].mask, n=bufsize)
    
    # unmasked = None
    lat,lon = [dir[unmasked] for dir in np.meshgrid(lat[::thin],lon[::thin])]
    if np.any(np.abs(lon)>180.) or np.any(np.abs(lat)>90.):
        raise ValueError, 'Ascii file %s has incorrect header. Lower left corner is (%f,%f); upper right corner is (%f,%f).'%(fname,lat.min(),lon.min(),lat.max(),lon.max())
    # lon,lat = [dir.ravel() for dir in np.meshgrid(lon[::thin],lat[::thin])]
    return np.vstack((lon,lat)).T*np.pi/180., unmasked, type

def raster_to_vals(name, path='.', thin=1, unmasked=None):
    """
    Converts a raster grid to a list of values where prediction is desired.
    If the unmasked argument is provided, the mask of the ascii grid is
    checked against that mask.
    """
    lon,lat,data,type = import_raster(name, path)
    data = grid_convert(data,'y-x+','x+y+')    
    if unmasked is not None:
        input_unmasked = True-data[::thin,::thin].mask
        if not np.all(unmasked == input_unmasked):
            if unmasked.shape == input_unmasked.shape:
                where_mismatch = np.where(input_unmasked != unmasked)
                import pylab as pl
                pl.clf()
                pl.plot(lon[where_mismatch[0]],lat[where_mismatch[1]],'k.',markersize=2)
                pl.savefig('mismatch.pdf')
                msg = '%s: covariate raster\'s mask does not match mask at the following pixels (in decimal degrees):\n'%name
                for i,j in zip(*where_mismatch):
                    msg += "\t%f, %f\n"%(lon[i],lat[j])
                msg += 'Image of mismatched points saved as mismatch.pdf'
            else:
                msg = '%s: covariate raster is not same shape as mask.'%name
            raise ValueError, msg
    
    return data.data[::thin,::thin][unmasked]

def get_mask_t(o, hf):
    x, unmasked, output_type = raster_to_locs(o.mask_name, thin=o.raster_thin, bufsize=o.bufsize, path=o.raster_path)
    if o.year is None:
        if 't' in hf.root.input_csv.colnames:
            raise ValueError, 'No year provided, but the model appears to be temporal.'
    else:
        x = np.vstack((x.T,o.year*np.ones(x.shape[0]))).T
    return unmasked, x


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

def get_x_from_recarray(input):
    lon = maybe_convert(input, 'lon', 'float')
    lat = maybe_convert(input, 'lat', 'float')
    if 't' in input.dtype.names:
        t = maybe_convert(input, 't', 'float')
        x = combine_st_inputs(lon,lat,t)
    else:
        x = combine_spatial_inputs(lon,lat)
    return x

def get_covariate_dict(M, o, unmasked):
    all_covariate_keys = M.covariate_keys
    covariate_dict = {}
    for k in all_covariate_keys:
        if k != 'm':
            try:
                covariate_dict[k] = raster_to_vals(k, path=o.raster_path, thin=o.raster_thin, unmasked=unmasked)
            except IOError:
                raise IOError, 'Covariate raster %s not found in path %s.'%(k+'.asc',o.raster_path)
    return covariate_dict

def reload_module(module):
    mod_path, mod_name = os.path.split(module)
    mod_basename, mod_ext = os.path.splitext(mod_name)
    mod_search_path = [mod_path, os.getcwd()] + sys.path
    mod = imp.load_module(mod_basename, *imp.find_module(mod_basename, mod_search_path))
    return mod, mod_name

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
        hfname = db._h5file.filename
    else:
        if input is None:
            raise ValueError, 'No input or preexisting db provided.'
        prev_db=None
        hfname = db
    
    rec2csv(input,'%s-input-data.csv'%hfname)

    x = get_x_from_recarray(input)
    x[:,:2]*=180./np.pi
    mod_inputs = tuple(x.T)
    
    non_cov_columns = {}
    if hasattr(mod, 'non_cov_columns'):
        non_cov_coltypes = mod.non_cov_columns
    else:
        non_cov_coltypes = {}
    non_cov_colnames = non_cov_coltypes.keys()

    covariate_keys = []
    for n in input.dtype.names:
        if n not in ['lon','lat','t']:
            if n in non_cov_colnames:
                non_cov_columns[n] = maybe_convert(input, n, non_cov_coltypes[n])
            else:
                covariate_keys.append(n)

    mod_inputs = mod_inputs + ('%s-input-data.csv'%hfname,covariate_keys,)

    # Create MCMC object, add metadata, and assign appropriate step method.

    if prev_db is None:
        M = pm.MCMC(mod.make_model(*mod_inputs,**non_cov_columns),db='hdf5',dbname=hfname,dbcomplevel=1,dbcomplib='zlib')
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


def reload_model(mod, hf):
    mod, mod_name = reload_module(mod)
    db = pm.database.hdf5.load(hf)
    return create_model(mod,db), db._h5file, mod, mod_name

def check_args(args, n):
    if len(args) != n:
        raise ValueError, 'You must supply exactly %i positional arguments. You supplied %i.'%(n,len(args))