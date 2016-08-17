'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy as np
from opr import sign
from pdb import set_trace as stop
#from opr import formOccu, sign

def HC(Cx, k_mtx, g_mtx, N, string_data):
    ne,no,ns = N
    occ,vir,aclist,Jlist,signlist = string_data
    # unpacking
    result = np.matrix(np.zeros(Cx.shape))

    for clm in xrange(Cx.shape[1]):
        C0 = np.array(Cx[:,clm]).reshape(ns,ns)
        ###########################################################
        
        # diagonal of K
        K = np.diag([np.einsum('i->',k_mtx[i,i]) for i in occ])
        # off-diagonal of K
        for i,jlist in enumerate(Jlist): # j is a list
            K[i][jlist] = np.array([k_mtx[k] for k in aclist[i]]) \
                        * np.array(signlist[i])
        sig = np.dot(K,C0) + np.dot(C0,K.T)
#        siga = np.dot(K,C0)
#        sigb = np.dot(C0,K.T) # = np.dot(K,C0.T).T
#        sig1 = siga + sigb

        ###########################################################
        # diagonal of G
        g_tmp = np.einsum('iijj->ij',g_mtx)
        G = np.diag([np.einsum('ij->',g_tmp[list(i)][:,list(i)]) \
                for i in occ])
        g_tmp = np.einsum('ijji->ij',g_mtx)
        G += np.diag([np.einsum('ij->', \
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
#            g_pqqq = np.einsum('pqqq->pq',g_mtx[list(occ[i])][:,list(vir[i])][:,:,list(vir[i])][...,list(vir[i])])
#            g_pqpp = np.einsum('pqpp->pq',g_mtx)
            g_tmp = g_rrpq + g_pqrr + g_prrq - g_rqpr
            G[i][jlist] += np.array([g_tmp[k] for k in aclist[i]]) \
                        * np.array(signlist[i])
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
                G[i][kl] += np.array(g_tmp0) * np.array(ssign)
        sig += 0.5 * (np.dot(G,C0) + np.dot(C0,G.T))    #########to be replaced
        #==================
        E_pq = np.zeros([no,no,ns,ns])
        for ia, p in enumerate(occ):
            E_pq[p,p,ia,ia] = 1
        for ia, jl in enumerate(Jlist):
            for j_index, pq in enumerate(aclist[ia]):
                p, q=pq
                E_pq[p,q,ia,jl[j_index]] = signlist[ia][j_index]
        D_pq = np.zeros([no,no,ns,ns])
        for p in xrange(no):
            for q in xrange(no):
                D_pq[p,q] = np.dot(E_pq[p,q], C0)
        G_pq = np.zeros([no,no,ns,ns])
        for ib, rlist in enumerate(occ):
            g_tmppq = np.einsum('pqrr->pq', \
                    g_mtx[:,:,rlist][...,rlist])
            G_pq[:,:,ib,ib] = g_tmppq
        for ib, jl in enumerate(Jlist):
            for j_index, rs in enumerate(aclist[ib]):
                r, s = rs
                G_pq[:,:,ib,jl[j_index]] = g_mtx[:,:,r,s] \
                        * signlist[ib][j_index]
        sig22pq = np.zeros([no,no,ns,ns])
        for p in xrange(no):
            for q in xrange(no):
                sig22pq[p,q] = np.dot(D_pq[p,q],G_pq[p,q].T)
        sig += np.einsum('pqij->ij',sig22pq)

        result[:,clm] = np.matrix(sig).reshape(ns*ns,1)
    #-----------------------------------------------------------------------
    return result 

def make_hdiag(k_mtx, g_mtx, N, string_data):
    ne, no, ns = N
    occ,vir,aclist,Jlist,signlist = string_data

    diagj = np.einsum('iijj->ij',g_mtx)
    diagk = np.einsum('ijji->ij',g_mtx)
    hdiag = []
    for ai,ao in enumerate(occ):
        for bi,bo in enumerate(occ):
            aocc = list(ao)
            avir = list(vir[ai])
            bocc = list(bo)
            bvir = list(vir[bi])
            e1 = k_mtx[aocc,aocc].sum() + k_mtx[bocc,bocc].sum()
            e2 = diagj[aocc][:,aocc].sum() + diagj[aocc][:,bocc].sum() \
               + diagj[bocc][:,aocc].sum() + diagj[bocc][:,bocc].sum() \
               + diagk[aocc][:,avir].sum() + diagk[bocc][:,bvir].sum()
            hdiag.append(e1 + e2*.5)
    return np.array(hdiag)
