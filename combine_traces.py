import numpy as np
import tables as tb
import sys

# handle system parameters
output_tracefile_filename = sys.argv[1]

input_tracefile_filenames = []
for ii in np.arange(1,10):
    try:
        input_tracefile_filenames = input_tracefile_filenames + [sys.argv[ii]]
    except:
        print('')


print("output_tracefile_filename")
print(output_tracefile_filename)

print("input_tracefile_filenames")
print(input_tracefile_filenames)


#
#
#
#input_tracefile_filenames = ['1_eaf_MCT1_FULL','1_eaf_MCT1_FULL_chain0']
#output_tracefile_filename = '1_eaf_MCT1_FULL_chain02'
#
#hf_out = tb.openFile(output_tracefile_filename,'w')
##hf_in = tb.openFile('1_eaf_MCT1_FULL')
#
#Ninputs = len(input_tracefile_filenames)
#
##for i,f in enumerate(input_tracefile_filenames):
#for i in np.arange(0,Ninputs):
#    
#    print(str(i)+' = working on tracefile:'+input_tracefile_filenames[i])
#    hf_in = tb.openFile(input_tracefile_filenames[i])
#    
#    # for first chain of first file, set up new output tacefile with appropriate componenets
#    if i==0:
#        hf_out.createGroup('/','metadata')
#        for key,node in hf_in.root.metadata._v_children.iteritems():
#            a=hf_out.createVLArray('/metadata',key,tb.ObjectAtom())
#            a.append(node[0])
#        hf_out.createTable('/','input_csv',hf_in.root.input_csv[:])
#        for key in ['generic_commit','input_filename','mod_commit','mod_name']:
#            hf_out.root.input_csv.setAttr(key, hf_in.root.input_csv.attrs[key])
#        
#        hf_out.createGroup('/','chain0')
#        t=hf_out.createTable('/chain0','PyMCsamples',hf_in.root.chain0.PyMCsamples[:])
#        a=hf_out.createVLArray('/chain0', '_state_', tb.ObjectAtom())
#        hf_out.createGroup('/chain0', 'group0')
#        
#        for key in hf_in.root.chain0.group0._v_children.iterkeys():
#            node_out = hf_out.createVLArray('/chain0/group0', key, tb.ObjectAtom())
#    
#    # build list of extant chains in this tracefile
#    chains=[]
#    for cc in np.arange(0,10):
#        try:
#            hf_in.root._g_checkHasChild('chain'+str(cc))
#            chains=chains+ [hf_in.root._f_getChild('chain'+str(cc))]
#        except:
#            print('')
#    Nchains = len(chains)
#    print('Found '+str(Nchains)+' chains in tracefile '+input_tracefile_filenames[i]+', merging to single chain0 in output')
#    
#    # append this chain to chain0 of output (only procede for first tracefile if >1 chain) 
#    
#    if(i==0): seq=np.arange(1,Nchains)
#    if(i>0): seq=np.arange(0,Nchains) 
#    if(len(seq)>0):
#        for cc in seq:
#            t.append(chains[cc].PyMCsamples[:])
#            a.append(chains[cc]._state_[0])
#            for key in hf_in.root.chain0.group0._v_children.iterkeys():
#                node_in = getattr(chains[cc].group0,key)
#                for val in node_in:
#                    node_out.append(val)
#
#hf_in.close()
#hf_out.close()
#
#
