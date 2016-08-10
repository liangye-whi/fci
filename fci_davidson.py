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
import time

import pdb

def FCI(mol):
    start_time = time.time()
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
    HD = numpy.zeros(ns*ns)
    for i in xrange(ns*ns):
        Ca = numpy.matrix(numpy.zeros(ns*ns)).T
        Ca[i] = 1.0
        siga = HC(Ca, k_mtx, g_mtx, ne, no, ns, Z, occu)
        HD[i] = siga[i,0]
        #Hd[0,0] +=  1e-3
    #---------------------------------------------------------
    L = 1           #number of initial guess basis vectors
    K = 0           #ground state
    B = numpy.matrix(numpy.eye(ns*ns,L))
    I = numpy.matrix(numpy.eye(ns*ns))
    AB = numpy.matrix(HC(B, k_mtx, g_mtx, ne, no, ns, Z, occu))
    print AB.shape
    print AB
    A = B.T * AB
    A = (A + A.T) / 2.
    A[0,0] += A[0,0]*0.01
    
    dt = numpy.matrix(I)
    for i in xrange(L):
        dt = dt - B[:,i] * B[:,i].T

    e, c = numpy.linalg.eigh(A)
    idx = e.argsort()[K]
    ek = e[idx]
    ck = c[:,idx]
    print 'E =',ek
    print ck
    
    #===criterion of convergence===
    crt = 1e-9
    #==============================
    iter_num = 0
    iter_limit = 10000
    ######################################################################
    print 'Now start iteration...'
    start_davidson = time.time()
    for M in xrange(L, iter_limit):
        q = (AB - ek * B) * ck
        #D.
        epsi = numpy.matrix(numpy.diag((ek - HD) ** (-1))) * q
        #E.
        d = dt * epsi
        res = numpy.linalg.norm(d)
        print 'res =', res
        if res < crt: 
            print 'Iteration Converged.'
            break

        b = d / res
        #F.
        B = numpy.c_[B,b]
        #G.
        AB = numpy.c_[AB, HC(b, k_mtx, g_mtx, ne, no, ns, Z, occu)]
        #H.
        a = B.T * AB[:,M]
        A = numpy.c_[A,a[:M,0]]
        A = numpy.r_[A,a.T]
        print A.shape
        dt = dt - b * b.T

        e, c = numpy.linalg.eigh(A)
        idx = e.argsort()[K]
        ek = e[idx]
        ck = c[:,idx]
        print 'E =',ek

    else:
        print 'Iteration Max (',iter_limit,'). Unconverged.'
    end_davidson = time.time()
    print 'Davidson time =', end_davidson - start_davidson, 'seconds.'
    #End of while 1
    #-----------------------------------------------------------------------
    return ek
