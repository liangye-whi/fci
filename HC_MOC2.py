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
    occ,vir,spstr,sps2i,aclist,Jlist,signlist = string_data
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
        sig1 = np.dot(K,C0) + np.dot(C0,K.T)
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
                ssign = list(signlist[i][j_index] * np.array(signlist[j]))
                for k, kac in enumerate(d):
                    if len(set(kac)) < 4:
                        ssign[k] = 0
                g_tmp0 = [g_mtx[k] for k in d]
                G[i][kl] += np.array(g_tmp0) * np.array(ssign)
        sig21 = 0.5 * (np.dot(G,C0) + np.dot(C0,G.T))    #########to be replaced
        #==================
        sig22 = np.zeros([ns,ns])
        
        for ia in xrange(ns):
            for ib in xrange(ns):
                #ia=ja,ib=jb
                for p in occ[ia]:
                    for r in occ[ib]:
                        sig22[ia,ib] += g_mtx[p,p,r,r] * C0[ia,ib]
                #ia=ja,ib!=jb
                for p in occ[ia]:
                    for r in occ[ib]:
                        for s in vir[ib]:
                            joc = spstr[ib] ^ (1 << r|1 << s) # XOR
                            j = sps2i[joc] # Epq
                            sgn = sign(spstr[ib],joc)
                            sig22[ia,ib] += sgn * g_mtx[p,p,r,s] * C0[ia,j]
                #ia!=ja,ib=jb
                for p in occ[ia]:
                    for q in vir[ia]:
                        for r in occ[ib]:
                            joc = spstr[ia] ^ (1 << p|1 << q) # XOR
                            j = sps2i[joc] # Epq
                            sgn = sign(spstr[ia],joc)
                            sig22[ia,ib] += sgn * g_mtx[p,q,r,r] * C0[j,ib]
                for p in occ[ia]:
                    for q in vir[ia]:
                        joc1 = spstr[ia] ^ (1 << p|1 << q) # XOR
                        ja = sps2i[joc1] # Epq
                        sgna = sign(spstr[ia],joc1)
                        for r in occ[ib]:
                            for s in vir[ib]:
                                joc2 = spstr[ib] ^ (1 << r|1 << s)
                                jb = sps2i[joc2] # Epq
                                sgn = sign(spstr[ib],joc2) * sgna
                                sig22[ia,ib] += sgn * g_mtx[p,q,r,s] * C0[ja,jb]
                                 
        sig = np.zeros([ns,ns])
        sig = sig1 + sig21 + sig22
    result[:,clm] = np.matrix(sig).reshape(ns*ns,1)
    #-----------------------------------------------------------------------
    return result 
