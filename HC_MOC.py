'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from opr import opr_E

def HC(C0, k_mtx, g_mtx, ne, no, ns, Z):
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns

    ###########################################################

    D = numpy.zeros([no,no,ns,ns])
    for r in xrange(no):
        for s in xrange(no):
            for k in xrange(ns):
                for j in xrange(ns):
#                    print opr_E(ne,no,r,s,kb,jb,Z) 
                    D[r,s,:,k] += opr_E(ne,no,r,s,k,j,Z) * C0[:,j] 
                    D[r,s,k,:] += opr_E(ne,no,r,s,k,j,Z) * C0[j,:]
#            print 'D(',r,',',s,') matrix'
#            print D[r,s,:,:]

    G = numpy.zeros([no,no,ns,ns])
    for p in xrange(no):
        for q in xrange(no):
            for r in xrange(no):
                for s in xrange(no):
                    G[p,q,:,:] += 0.5 * g_mtx[p*no+q,r*no+s] * D[r,s,:,:]
#            print 'G(',p,',',q,') matrix'
#            print G[p,q,:,:]

    sig2 = numpy.zeros([ns,ns])
    for ia in range(ns):
        for ib in range(ns):
            for p in range(no):
                for q in range(no):
                    for ka in range(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ia,ka,Z) * G[p,q,ka,ib]
                    for kb in range(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ib,kb,Z) * G[p,q,ia,kb]

    sig1 = numpy.zeros([ns,ns])
    for ia in range(ns):
        for ib in range(ns):
            for p in range(no):
                for q in range(no):
                    sig1[ia,ib] += k_mtx[p,q] * D[p,q,ia,ib]

    sig = numpy.zeros([ns,ns])
    sig = sig1 + sig2
    #-----------------------------------------------------------------------
    return sig
