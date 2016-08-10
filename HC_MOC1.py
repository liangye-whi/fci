'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from opr import formOccu
from opr_old import opr_E

def HC(C0, k_mtx, g_mtx, ne, no, ns, Z, occu):
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns
    #print occu
    ###########################################################
    sig1 = numpy.zeros([ns,ns])
    for i in xrange(ns):
        # diagonal elem
        res='0'*(no+2-len(bin(occu[0][i])))+bin(occu[0][i])[2:]
        cu = [k for k in xrange(no) if res[k] == '1']
        sumkpp = sum([k_mtx[p,p] for p in cu])
        sig1[i,:] += sumkpp * C0[i,:] # Ia
        sig1[:,i] += sumkpp * C0[:,i] # Ib
        
        # nondiagonal elem
        ucu = [k for k in xrange(no) if res[k] == '0']
        for p in cu:
            for q in ucu:
                assert p != q ##debug##
                joc = occu[0][i] ^ ((1 << (no-q-1)) + (1 << (no-p-1))) # XOR
                j = occu[1][str(joc)] # Epq
                sig1[i,:] += k_mtx[p,q] * C0[j,:] # Ia
                sig1[:,i] += k_mtx[p,q] * C0[:,j] # Ib
    
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

    sig = numpy.zeros([ns,ns])
    sig = sig1 + sig2
    print 'sig2'
    print sig2
    #-----------------------------------------------------------------------
    return sig
