import numpy as np
import itertools as it
import operator

def math_C(n,k):  
        return reduce(operator.mul, range(n - k + 1, n + 1)) / reduce(operator.mul, range(1, k +1))  

def sign(string0, string1):
    ss = string1 - string0
    def count_bit1(n):
        # see Hamming weight problem and K&R C program
        return bin(n).count('1')
    if ss > 0:
        # string1&ss gives the number of 1s between two strings
        return (-1) ** (count_bit1(string1&ss))
    elif ss == 0:
        return 1
    else:
        return (-1) ** (count_bit1(string0&(-ss)))

def construct_string_data(N):
    ne,no,ns = N
    #=====================
    # occ, vir
    occ = tuple(it.combinations(range(no),ne))
    vir = tuple(tuple(set(range(no)).difference(i)) for i in occ)
    # occ and vir for i_th string
    #=====================
    # spstr, sps2i
    spstr = []
    # i_th spin string
    sps2i = {}
    # spin string to its order i
    for i,iocc in enumerate(occ):
        sps_tmp = 0
        for j in iocc:
            sps_tmp |= 1 << j
        spstr.append(sps_tmp)
        sps2i[sps_tmp] = i
    #=====================
    # aclist
    # (annihilation, creation)
    aclist = [tuple(it.product(occ[i],vir[i])) for i in xrange(ns)]
    #=====================
    # Jlist
    Jlist = []
    for i,iac in enumerate(aclist):
        J0 = []
        str0 = spstr[i]
        for a,c in iac:
            str1 = str0 ^ (1<<a|1<<c)
            J0.append(sps2i[str1])
        Jlist.append(J0)

    for i, ac in enumerate(aclist):
        str0 = spstr[i]
        for p,q in ac:
            str1 = str0 ^ (1 << p | 1 << q)
    #=====================
    # signlist
    signlist = []
    for i,jlist in enumerate(Jlist):
        s_list = []
        for j in jlist:
            s_list.append(sign(spstr[i],spstr[j]))
        signlist.append(s_list)
    #=====================
    return occ, vir, aclist, Jlist, signlist

def construct_matrix(h_mtx, g_mtx, N, string_data):
    ne, no, ns = N
    occ,vir,aclist,Jlist,signlist = string_data
    #=========================================
    # k_mtx
    k_mtx = np.einsum('irrj->ij',g_mtx)
    k_mtx *= -0.5
    k_mtx = k_mtx + h_mtx
    del h_mtx
    print 'k matrix done.'
    #=========================================
    # K_MOC
    # diagonal of K
    K_MOC = np.diag([np.einsum('i->',k_mtx[i,i]) for i in occ])
    # off-diagonal of K
    for i,jlist in enumerate(Jlist): # j is a list
        K_MOC[i][jlist] = np.array([k_mtx[k] for k in aclist[i]]) \
                    * np.array(signlist[i])
    #=========================================
    # G_MOC
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
    
    #=========================================
    # E_pq
    E_pq = np.zeros([no,no,ns,ns])
    for ia, p in enumerate(occ):
        E_pq[p,p,ia,ia] = 1
    for ia, jl in enumerate(Jlist):
        for j_index, pq in enumerate(aclist[ia]):
            p, q=pq
            E_pq[p,q,ia,jl[j_index]] = signlist[ia][j_index]
    
    #=========================================
    # G_pq

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
    # HD
    diagj = np.einsum('iijj->ij',g_mtx)
    diagk = np.einsum('ijji->ij',g_mtx)
    hdiag = []
    for ai,ao in enumerate(occ):
        for bi,bo in enumerate(occ):
            aocc = list(ao)
            avir = list(vir[ai])
            bocc = list(bo)
            bvir = list(vir[bi])
            e1 = np.einsum('i->',k_mtx[aocc,aocc]) + np.einsum('i->',k_mtx[bocc,bocc])
            e2 = np.einsum('ij->',diagj[aocc][:,aocc]) + np.einsum('ij->',diagj[aocc][:,bocc]) \
               + np.einsum('ij->',diagj[bocc][:,aocc]) + np.einsum('ij->',diagj[bocc][:,bocc]) \
               + np.einsum('ij->',diagk[aocc][:,avir]) + np.einsum('ij->',diagk[bocc][:,bvir])
            hdiag.append(e1 + e2*.5)
    #=========================================

    return k_mtx, K_MOC, G_MOC, E_pq, G_pq, np.array(hdiag)
