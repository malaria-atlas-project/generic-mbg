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
import pymc as pm
from pylab import csv2rec,rec2csv
import decimal

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