'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from pyscf import gto, scf, ao2mo
from opr import constructZ, opr_E, math_C
from HC_NR import HC

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
    for i in xrange(no):
        for j in xrange(i+1,no):
            print i,j, k_mtx[i,j],  k_mtx[j,i] , k_mtx[i,j]-  k_mtx[j,i]
    #---------------------------------------------------------

    C0 = numpy.ones([ns,ns])
    #C0x = C0 ** 2
    #A = sum(sum(C0x))
    #C0 =  C0 / (A**0.5)
    C0 =  C0 / ns
    print 'C0'
    print C0

    E0 = 0
    #===criterion of convergence===
    crt = 1e-9
    #==============================
    iter_num = 0
    iter_limit = 500
    ######################################################################

    while iter_num < iter_limit:
        # Davidson step
        iter_num += 1
        ##  ***  ##
        sig = HC(C0, k_mtx, g_mtx, ne, no, ns,Z)
    #    print 'sig'
    #    print sig

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
    else:
        print 'Iteration Max (',iter_max,'). Unconverged.'

    #End of while 1
    #-----------------------------------------------------------------------
    return E0