# Author: Anand Patil
# Copyright (c) 2009 Anand Patil
# License: GPL

from setuptools import setup
from numpy.distutils.misc_util import Configuration
import os
config = Configuration('generic_mbg',parent_package=None,top_path=None)

config.add_extension(name='histogram_utils',sources=['generic_mbg/histogram_utils.f'])

prefix='/usr/local/bin'

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="0.1",
            description="The Malaria Atlas Project's generic MBG code.",
            author="Peter Gething and Anand Patil", 
            author_email="map@map.ox.ac.uk",
            url="www.map.ox.ac.uk",
            packages=['generic_mbg'],
            license="Creative commons BY-NC-SA",
            **(config.todict()))
    for ex_fname in ['mbg-infer','mbg-map','mbg-validate','mbg-scalar-priors','mbg-realize-prior']:
        os.system('chmod ugo+x %s'%ex_fname)
        os.system('cp %s %s'%(ex_fname,prefix))

