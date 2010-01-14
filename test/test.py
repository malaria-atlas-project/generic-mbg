from generic_mbg import asc_to_locs, asc_to_vals

l1,u1 = asc_to_locs('africa.asc',thin=1,bufsize=0)
l5,u5 = asc_to_locs('africa.asc',thin=5,bufsize=5)

africa1 = asc_to_vals('africa.asc',unmasked=u1)
africa5 = asc_to_vals('africa.asc',thin=5,unmasked=u5)