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
    
def get_reduce_finalize(mod):
    if hasattr(mod, 'extra_reduce_fns'):
        extra_reduce_fns = mod.extra_reduce_fns
        extra_finalize = mod.extra_finalize
    else:
        extra_reduce_fns = []
        extra_finalize = None

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

def parse_and_check_args(p, n):
    (o, args) = p.parse_args()
    if len(args) != n:
        raise ValueError, 'You must supply exactly %i positional arguments. You supplied %i.'%(n,len(args))
    return o, args