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
    #E = np.eye(ns*ns)
    #H_mtx = np.zeros([ns*ns,1])
    #for i in xrange(ns):
    #    for j in xrange(ns):
    #        B = np.zeros(ns*ns)
    #        B[i*ns + j] = 1.0
    #        B = B.reshape(-1,1)

    #        AB = HC(B, k_mtx, g_mtx, K_MOC, G_MOC, E_pq, G_pq, N, string_data)
    #        A = np.dot(E.T, AB)
    #        H_mtx = np.c_[H_mtx,A]
    #H_mtx = H_mtx[:,1:]
    #np.set_printoptions(threshold='nan')
    #print H_mtx
    #print H_mtx.shape
    i, j = 0,0
    p, q = 20,0
    print string_data[0]
    print string_data[0][i]
    print string_data[0][j]
    print string_data[0][p]
    print string_data[0][q]
    print 'k', k_mtx[1,9]
    print 'g', .5*(g_mtx[3,1,9,2]+g_mtx[3,2,9,1])
    B = np.zeros(ns*ns)
    B[i*ns + j] = 1.0
    B = B.reshape(-1,1)
    E = np.zeros(ns*ns)
    E[p*ns + q] = 1.0
    E = E.reshape(-1,1)
    AB = HC(B, k_mtx, g_mtx, K_MOC, G_MOC, E_pq, G_pq, N, string_data)
    A = np.dot(E.T, AB)
    print A
#    e, c = np.linalg.eigh(H_mtx)
    #-----------------------------------------------------------------------
    return 0
