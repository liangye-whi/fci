'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy as np
from pyscf import gto, scf, ao2mo
from opr import construct_string_data, construct_matrix, math_C
from HC import HC
import time

from pdb import set_trace as stop
def FCI(mol):
    start_time = time.time()
    m = scf.RHF(mol)
    m.kernel()

#    print 'mo_coeff'
#    print(m.mo_coeff)
#    print 'mo_occ'
#    print(m.mo_occ)
#    print 'mo_energy'
#    print(m.mo_energy)

    ne = mol.nelectron/2 #electron per string
    no = len(m.mo_energy)
    ns = math_C(no,ne)
    N = (ne,no,ns)
    print 'ne =', ne, 'no =', no, 'ns =', ns

    string_data = construct_string_data(N)

    ###########################################################

    h_ao_mtx = scf.hf.get_hcore(mol)
    h_mtx = np.dot(np.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
    del h_ao_mtx
    print 'h matrix done.' 

    g_mtx = ao2mo.kernel(mol, m.mo_coeff, compact=False)
    g_mtx = g_mtx.reshape(no,no,no,no)
    print 'g matrix done.'

    k_mtx, K_MOC, G_MOC, E_pq, G_pq, HD = construct_matrix(h_mtx, g_mtx, N, string_data)
    del h_mtx
    #---------------------------------------------------------
    K = 0           #ground state
    B = np.zeros(ns*ns)
    B[0] = 1.0
    B = B.reshape(-1,1)
    AB = HC(B, k_mtx, g_mtx, K_MOC, G_MOC, E_pq, G_pq, N, string_data)
    A = np.dot(B.T, AB)
    HD[0] += abs(HD[0]*0.001)
    
    e, c = np.linalg.eigh(A)
    idx = e.argsort()[K]
    ek = e[idx]
    ck = c[:,idx]
    print 'E =',ek
    
    #===criterion of convergence===
    crt = 1e-8
    #==============================
    iter_num = 0
    iter_limit = 100
    ######################################################################
    print 'Now start iteration...'
    start_davidson = time.time()
    for M in xrange(1, iter_limit):
        b = np.dot((AB - ek * B), ck)
        #D.
        b = (ek - HD) ** (-1) * b.reshape(-1)
        b.shape = (-1,1)
        #E.
        b = b - np.einsum('i,ji->j', np.dot(b.T, B).reshape(-1), B).reshape(-1,1)
        # orthogonalization
        res = np.linalg.norm(b)
        print 'res =', res
        if res < crt: 
            print 'Iteration Converged.'
            break
        b /= res
        #F.
        B = np.c_[B,b]
        #G.
        AB = np.c_[AB, HC(b, k_mtx, g_mtx, K_MOC, G_MOC, E_pq, G_pq, N, string_data)]
        #H.
        a = np.dot(B.T, AB[:,M]).reshape(-1,1)
        A = np.c_[A,a[:M,0]]
        A = np.r_[A,a.T]
        print A.shape
        e, c = np.linalg.eigh(A)
        idx = e.argsort()[K]
        ek = e[idx]
        ck = c[:,idx]
        print('E = %.12f' % (ek))

    else:
        print 'Iteration Max (',iter_limit,'). Unconverged.'
    end_davidson = time.time()
    print 'Davidson time =', end_davidson - start_davidson, 'seconds.'
    #-----------------------------------------------------------------------
    return ek
