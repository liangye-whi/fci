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

def FCI(mol):
    m = scf.RHF(mol)
    m.kernel()

    #print 'mo_coeff'
    #print(m.mo_coeff)
    #print 'mo_occ'
    #print(m.mo_occ)
    #print 'mo_energy'
    #print(m.mo_energy)

    ne = mol.nelectron/2 #electron per string
    no = len(m.mo_energy)
    ns = math_C(no,ne)
    Z = constructZ(ne,no)
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns

    ###########################################################

    h_ao_mtx = scf.hf.get_hcore(mol)
    h_mtx = numpy.dot(numpy.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
    #print 'h matrix'
    #print h_mtx

    g_mtx = ao2mo.kernel(mol, m.mo_coeff, compact=False)
    #print 'g matrix'
    #print g_mtx[0,0] # g[p,q,r,s] => g[p*no+q,r*no+s]


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
    ######################################################################
    print 'Now start iteration...'
    H = numpy.zeros([ns*ns,ns*ns])
    for i in xrange(ns):
        for j in xrange(ns):
            C0 = numpy.zeros([ns,ns])
            C0[i,j] = 1
            ##  ***  ##
            sig = HC(C0, k_mtx, g_mtx, ne, no, ns, Z, occu)
            sig = sig.reshape(-1)
            H[i*ns+j] = sig
    #print H
    #numpy.savetxt('test.out', H)
    e, c = numpy.linalg.eig(H)
    #-----------------------------------------------------------------------
    return min(e)

