'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from pyscf import gto, scf, ao2mo
from opr_E import constructZ, opr_E, math_C

def FCI(mol):
    m = scf.RHF(mol)
    m.kernel()

    print 'mo_coeff'
    print(m.mo_coeff)
    print 'mo_occ'
    print(m.mo_occ)
    print 'mo_energy'
    print(m.mo_energy)

    ne = mol.nelectron/2 #electron per string
    no = len(m.mo_energy)
    ns = math_C(no,ne)
    Z = constructZ(ne,no)
    #print 'Z'
    #print Z
    print 'ne =', ne, 'no =', no, 'ns =', ns

    ###########################################################

    h_ao_mtx = scf.hf.get_hcore(mol)
    h_mtx = numpy.dot(numpy.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
    print 'h matrix'
    print h_mtx

    g_mtx = ao2mo.kernel(mol, m.mo_coeff, compact=False)
    print 'g matrix'
    print g_mtx[0,0]

    k_mtx = numpy.zeros([no,no])
    for i in xrange(no):
        for j in xrange(no):
            for r in xrange(no):
                k_mtx[i,j] -= g_mtx[i*no+r,r*no+j]
            k_mtx[i,j] *= 0.5
    k_mtx = k_mtx + h_mtx

    #---------------------------------------------------------

    #===criterion of convergence===
    crt = 1e-9
    #==============================

    C0 = numpy.ones([ns,ns])
    #C0x = C0 ** 2
    #A = sum(sum(C0x))
    #C0 =  C0 / (A**0.5)
    C0 =  C0 / ns
    print 'C0'
    print C0

    E0 = 0
    iter_num = 0;

    ######################################################################

    while 1:
        iter_num += 1
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
    #    print 'sig'
    #    print sig

        #####################################################################################

        # Davidson step
        E1 = E0
        E0 = sum(sum(C0 * sig))
        cdav = - (1 - E0)**(-1) * (sig - E0 * C0)
        res = sum(sum(cdav**2))**.5
        print 'Iter ',iter_num, '=============='
        print 'E =', E0
        print 'dE =', E0 - E1
        print '|res| =', res
        if res < crt:
            print 'Iteration Converged.'
            break

        print 'Cnew'
        C0 = C0 + cdav
        # normalization
        A = sum(sum(C0**2))
    #    print 'normalization factor = ', A
        C0 =  C0 / (A**0.5)
        print C0
    #    C1x = C1 ** 2
    #    A = sum(sum(C1x))
    #    print 'normalization factor =', A
        
    #End of while 1
    #-----------------------------------------------------------------------
    return E0
