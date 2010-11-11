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

delayed_import_str ="""import matplotlib
matplotlib.use('PDF')
matplotlib.interactive(False)

from map_utils import *
import tables as tb
import numpy as np
import pymc as pm
import pylab as pl"""

def get_nuggets_and_obs(mod, mod_name, M):
    try:
        nugget_labels = getattr(mod, nugget_labels)
        nuggets = dict([(getattr(M,k),getattr(M,v)) for k,v in mod.nugget_labels.iteritems()])
        obs_labels = getattr(mod, obs_labels)
    except:
        cls,inst,tb = sys.exc_info()
        new_inst = cls('Could not import %s from %s. Tell Anand. Original error message:\n\n\t%s'%(n,mod_name,inst.message))
        raise cls,new_inst,tb
    return nugget_labels, obs_labels

def parse_quantiles(o):
    if len(o.quantile_list) == 0:
        q = []
    else:
        q = map(float, o.quantile_list.split(' '))
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

def get_mask_t(o, hf):
    x, unmasked, output_type = raster_to_locs(o.mask_name, thin=o.raster_thin, bufsize=o.bufsize, path=o.raster_path)
    if o.year is None:
        if 't' in hf.root.input_csv.colnames:
            raise ValueError, 'No year provided, but the model appears to be temporal.'
    else:
        x = np.vstack((x.T,o.year*np.ones(x.shape[0]))).T
    return unmasked, x
    
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
    if hasattr(input, 't'):
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

def reload_model(mod, hf):
    mod, mod_name = reload_module(mod)
    return create_model(mod,pm.database.hdf5.load(hf)), hf, mod, mod_name

def check_args(args, n):
    if len(args) != n:
        raise ValueError, 'You must supply exactly %i positional arguments. You supplied %i.'%(n,len(args))