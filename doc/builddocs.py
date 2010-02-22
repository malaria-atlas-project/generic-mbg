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