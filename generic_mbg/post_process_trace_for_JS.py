import tables as tb
import numpy as np
import os
import time

##############################################################################################
# utiliy function
def check_if_chain (child_name):
    if(child_name.find('chain') != -1):
        return (True)
    return (False)
check_if_chain_v=np.vectorize(check_if_chain)
#############################################################################################
def collapse_trace_to_chain0(hf_path):

    # import hdf5 tracefile
    hf = tb.openFile(hf_path,mode="a")

    # how many chains present in trace?
    child_list=np.array(hf.root._v_children.keys())
    ischain = check_if_chain_v(child_list)
    nchains=sum(ischain)

    # if we don;t have more than one chain, simply return
    if(nchains==0):
        raise ValueError ('No chains found')

    if(nchains==1):
        print('No additional chains found, returning')
        return(0)

    # if we have multiple chains, define list of their names above chain0
    extra_chain_names=child_list[ischain]
    extra_chain_names=extra_chain_names[extra_chain_names!='chain0']

    # loop through additional chains
    for extra_chain_name in extra_chain_names:

        chain=hf.root._f_getChild(extra_chain_name)

        # first deal with PyMCsamples table
        PyMCsamples=chain.PyMCsamples

        # check matches those of Chain0
        if(PyMCsamples.colnames!=hf.root.chain0.PyMCsamples.colnames):
            raise ValueError ('colnames of '+extra_chain_name+' do not match chain0')

#        # loop through 'columns' of PyMCsamlpes table, appending to those of chain0
#        for colname in PyMCsamples.colnames:
#
#            # get this column of chain0 (check no hanging dimensions)
#            col_0=hf.root.chain0.PyMCsamples.col(colname)
#            col_0=col_0.squeeze()
#
#            # get this column of this additional chain (check no hanging dimensions)
#            col_new=PyMCsamples.col(colname)
#            col_new=col_new.squeeze()
#
#            # check columns have same dimensionality
#            if( len(col_0.shape) != len(col_new.shape) ):
#                raise ValueError ('column '+colname+' dimension mis-match: chain0 = '+str(col_0.shape)+'; '+extra_chain_name+' = '+str(col_new.shape))
#
#            # concatenate new col to chain0
#            if(len(col_0.shape)==1):
#                col_conc=np.hstack((col_0,col_new))
#
#            if(len(col_new.shape)==2):
#                col_conc=np.vstack((col_0,col_new))        
#
#            if((len(col_new.shape)!=1) & (len(col_new.shape)!=2)):
#                raise ValueError('Something fishy going on: column '+colname+' is not 1-d or 2-d')
#
#            # overwrite this col of chain0 with this concatenated version
#            col_0=col_conc

        # initialise row iterator
        trow=hf.root.chain0.PyMCsamples.row

        # loop through number of rows (realisations) in this additional Chain PyMCsamples table
        
        jj=0
        for ii in xrange(0,100):        #PyMCsamples.nrows
            jj=jj+1
            if (jj==10):
                print('copying PyMCsamples  on row '+str(ii)+' of '+str(PyMCsamples.nrows))
                jj=0

            # loop through columns and copy across this row of each
            for colname in PyMCsamples.colnames:
                trow[colname]=hf.root.chain1.PyMCsamples.col(colname)[ii]

            # append this row to target table
            trow.append()

        # flush table to implement changes 
        hf.root.chain0.PyMCsamples.flush()


        # next deal with group0 table
        group0=chain.group0

        # check nodes matches those of Chain0.group0
        if(group0._v_children.keys()!=hf.root.chain0.group0._v_children.keys()):
            raise ValueError ('children of '+extra_chain_name+' do not match chain0')

        # loop through nodes of this chain's group0 group, appending to those of chain0
        for nodename in group0._v_children.keys():

            # get this column of chain0
            node_0=hf.root.chain0.group0._f_getChild(nodename) 

            # get this column of this additional chain
            node_new=group0._f_getChild(nodename)
            
            # check that the ..input-data.csv file exists (necessary to invoke group0.C or group0.sp_sub_f) and if not run mbg-describe-tracefile
            try:
                temp=group0.C[0]
            except IOError:
                print "IOError when invoking group0.C: running mbg-describe-tracefile"
                cmd = 'mbg-describe-tracefile '+hf_path
                os.sytem(cmd)
                time.sleep(10)
                temp=group0.C[0]

            # check columns have same dimensionality
            if( len(node_0.shape) != len(node_new.shape) ):
                raise ValueError ('node '+nodename+' dimension mis-match: chain0 = '+str(node_0.shape)+'; '+extra_chain_name+' = '+str(node_new.shape))

            # append contents of this node on this chain to that on chain0
            jj=0
            for index in np.arange(0,100):  #
            jj=jj+1
            if (jj==10):
                print('copying group0  on row '+str(ii)+' of '+str(PyMCsamples.nrows))
                jj=0
                node_0.append(node_new[index])

        # remove chain from tracefile
        chain._f_remove(recursive=True)

#############################################################################################