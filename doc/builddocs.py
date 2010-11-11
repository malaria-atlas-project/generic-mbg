# Copyright (C) 2010 Anand Patil
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

from docutils.core import publish_file
import os
userstr = ''.join(publish_file(file('../README.rst'),writer_name='latex').\
            split("""%___________________________________________________________________________""")[1:]).\
            replace("""\end{document}""",'').\
            replace('section*','section').\
            replace('\label{','\label{sec:')

userstr2 = []            
insec = False
for line in userstr.splitlines():
    if line=='}':
        insec=False
    if not insec:
        userstr2.append(line)
    if line.find("""section{""")>-1:
        insec=True

userstr = '\n'.join(userstr2)

usertex=file('user.tex','w')
usertex.write(r"""\chapter{User's guide} 
\label{chap:user} 

This section is generated automatically from \texttt{README.rst}, the same file that produces the \href{github.com/malaria-atlas-project/generic-mbg}{GitHub documentation}. 
""")

userlines = []
for line in userstr.splitlines():
    if line.find('pdfbookmark')==-1 and line.find('hypertarget')==-1:
        userlines.append(line)
    else:
        userlines.append('')

usertex.write('\n'.join(userlines))

usertex.close()

for i in xrange(2):
    os.system('pdflatex manual.tex')