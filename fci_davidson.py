'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy as np
from pyscf import gto, scf, ao2mo
from opr import construct_string_data, math_C
from HC_MOC import HC, make_hdiag
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
    string_data = construct_string_data(N)
    occ,vir,aclist,Jlist,signlist = string_data
    print 'ne =', ne, 'no =', no, 'ns =', ns

    ###########################################################

    h_ao_mtx = scf.hf.get_hcore(mol)
    h_mtx = np.dot(np.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
    del h_ao_mtx
    print 'h matrix done.' 

    g_mtx = ao2mo.kernel(mol, m.mo_coeff, compact=False)
    g_mtx = g_mtx.reshape(no,no,no,no)
    print 'g matrix done.'


    k_mtx = np.einsum('irrj->ij',g_mtx)
    k_mtx *= -0.5
    k_mtx = k_mtx + h_mtx
    del h_mtx
    print 'k matrix done.'
    #=========================================================

    # diagonal of K
    K_MOC = np.diag([np.einsum('i->',k_mtx[i,i]) for i in occ])
    # off-diagonal of K
    for i,jlist in enumerate(Jlist): # j is a list
        K_MOC[i][jlist] = np.array([k_mtx[k] for k in aclist[i]]) \
                    * np.array(signlist[i])
    #------

    g_tmp = np.einsum('iijj->ij',g_mtx)
    G_MOC = np.diag([np.einsum('ij->',g_tmp[list(i)][:,list(i)]) \
            for i in occ])
    g_tmp = np.einsum('ijji->ij',g_mtx)
    G_MOC += np.diag([np.einsum('ij->', \
            g_tmp[list(occ[i])][:,list(vir[i])]) \
            for i in xrange(ns)])
    # off-diagonal of G
    for i,jlist in enumerate(Jlist):
        g_rrpq = np.einsum('rrpq->pq', \
                g_mtx[list(occ[i])][:,list(occ[i])])
        g_pqrr = np.einsum('pqrr->pq', \
                g_mtx[:,:,list(occ[i])][...,list(occ[i])])
        g_prrq = np.einsum('prrq->pq', \
                g_mtx[:,list(vir[i])][:,:,list(vir[i])])
        g_rqpr = np.einsum('rqpr->pq', \
                g_mtx[list(occ[i])][...,list(occ[i])])
        g_tmp = g_rrpq + g_pqrr + g_prrq - g_rqpr
        G_MOC[i][jlist] += np.array([g_tmp[k] for k in aclist[i]]) \
                    * np.array(signlist[i])
    del g_tmp
    for i,jl in enumerate(Jlist):
        for j_index,j in enumerate(jl):
            #----
            ijindex = Jlist[j].index(i)
            kl = Jlist[j][:]
            d = [aclist[i][j_index]+k for k in aclist[j]]
            ssign = list(signlist[i][j_index] \
                    * np.array(signlist[j]))
            for k, kac in enumerate(d):
                if len(set(kac)) < 4:
                    ssign[k] = 0
            g_tmp0 = [g_mtx[k] for k in d]
            G_MOC[i][kl] += np.array(g_tmp0) * np.array(ssign)
    
    #--------
    
    E_pq = np.zeros([no,no,ns,ns])
    for ia, p in enumerate(occ):
        E_pq[p,p,ia,ia] = 1
    for ia, jl in enumerate(Jlist):
        for j_index, pq in enumerate(aclist[ia]):
            p, q=pq
            E_pq[p,q,ia,jl[j_index]] = signlist[ia][j_index]
    
    #--------

    G_pq = np.zeros([no,no,ns,ns])
    for ib, rlist in enumerate(occ):
        g_tmppq = np.einsum('pqrr->pq', \
                g_mtx[:,:,rlist][...,rlist])
        G_pq[:,:,ib,ib] = g_tmppq
    del g_tmppq
    for ib, jl in enumerate(Jlist):
        for j_index, rs in enumerate(aclist[ib]):
            r, s = rs
            G_pq[:,:,ib,jl[j_index]] = g_mtx[:,:,r,s] \
                    * signlist[ib][j_index]

    #=========================================================
    #---------------------------------------------------------

    start_HD_time = time.time()
    HD = make_hdiag(k_mtx,g_mtx,N,string_data)
    print 'H0 matrix done.'
    end_HD_time = time.time()
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
    print 'HD time =', end_HD_time - start_HD_time, 'seconds.'
    print 'Davidson time =', end_davidson - start_davidson, 'seconds.'
    #-----------------------------------------------------------------------
    return ek
