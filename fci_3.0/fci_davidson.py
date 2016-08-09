'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from pyscf import gto, scf, ao2mo
from opr import constructZ, formOccu, math_C
from HC_MOC2 import HC

import pdb

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
    print g_mtx[0,0] # g[p,q,r,s] => g[p*no+q,r*no+s]


    k_mtx = numpy.zeros([no,no])
    for i in xrange(no):
        for j in xrange(no):
            for r in xrange(no):
                k_mtx[i,j] -= g_mtx[i*no+r,r*no+j]
            k_mtx[i,j] *= 0.5
    k_mtx = k_mtx + h_mtx

    #---------------------------------------------------------

    occu = formOccu(ne,no,ns,Z)

    #---------------------------------------------------------
    
    C0 = numpy.zeros([ns,ns])
    C0[0,0] = 1.
    #C0 =  C0 / ns
    print 'C0'
    print C0
    sig = HC(C0, k_mtx, g_mtx, ne, no, ns, Z, occu)
    E0 = sum(sum(C0 * sig))
    print 'E0'
    print E0
    
    Hd = numpy.zeros([ns,ns])
    for i in xrange(ns):
        print i
        for j in xrange(ns):
            Ca = numpy.zeros([ns,ns])
            Ca[i,j] = 1.0
            siga = HC(Ca, k_mtx, g_mtx, ne, no, ns, Z, occu)
            Ea = siga[i,j]
            Hd[i,j] = Ea 
        Hd[0,0] +=  abs(Hd[0,0]) * 1e-2
    #print Hd
    #pdb.set_trace()
    #===criterion of convergence===
    crt = 1e-9
    #==============================
    iter_num = 0
    iter_limit = 10000
    ######################################################################
    print 'Now start iteration...'
    while iter_num < iter_limit:
        # Davidson step
        iter_num += 1

        cdav = - (Hd - E0)**(-1) * (sig - E0 * C0)
        res = sum(sum(cdav**2))**.5
        C0 = C0 + cdav
        # normalization
        A = sum(sum(C0**2))
    #    print 'normalization factor = ', A
        C0 =  C0 / (A**0.5)
        sig = HC(C0, k_mtx, g_mtx, ne, no, ns, Z, occu)
        E1 = sum(sum(C0 * sig))
        dE = E1 - E0
        E0 = E1
        print 'Iter ',iter_num, '=============='
        print 'E =', E0
        print 'dE =', dE
        print '|res| =', res
        print 'Cnew'
        print C0
        if res < crt:
            print 'Iteration Converged.'
            break
    else:
        print 'Iteration Max (',iter_limit,'). Unconverged.'
    #End of while 1
    #-----------------------------------------------------------------------
    return E0
