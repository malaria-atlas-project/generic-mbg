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

from setuptools import setup
from numpy.distutils.misc_util import Configuration
import os
import subprocess
import sys

prefix='/usr/local/bin'
for v in sys.argv:
    if v.find('--executable-dir')>-1:
        sys.argv.remove(v)
        prefix=v.split('=')[1].strip()
        

config = Configuration('generic_mbg',parent_package=None,top_path=None)

config.add_extension(name='histogram_utils',sources=['generic_mbg/histogram_utils.f'])

def get_syscall_output(str):
    process = subprocess.Popen(str, stdout=subprocess.PIPE, shell=True)
    os.waitpid(process.pid, 0)
    return process.stdout.read().strip()

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
    for ex_fname in ['mbg-map',
                    'mbg-3dmap',
                    'mbg-validate',
                    'mbg-scalar-priors',
                    'mbg-realize-prior',
                    'mbg-covariate-traces',
                    'mbg-decluster',
                    'get_declustered_sample.R',
                    'mbg-infer',
                    'mbg-describe-tracefile',
                    'mbg-init-specializing-module',
                    'mbg-init-user-account',
                    'mbg-evaluate-survey',
                    'mbg-areal-predict',
                    'mbg-text-traces',
                    'mbg-geweke-diagnostic']:

        commit = get_syscall_output('git show --pretty=format:"%H" --quiet')
        pythonpath = get_syscall_output('which python')

        path_fname = os.path.join(prefix,ex_fname)
        os.system('rm %s'%path_fname)
        file(path_fname,'w').write('#!%s\ngeneric_commit = "%s"\n'%(pythonpath, commit))
        os.system('cat %s >> %s'%(ex_fname,path_fname))
        os.system('chmod ugo+x %s'%path_fname)

