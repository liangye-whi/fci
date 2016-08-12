'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy as np
from opr import sign
#from opr import formOccu, sign

def HC(Cx, k_mtx, g_mtx, N, string_data):
    ne,no,ns = N
    occ,vir,spstr,dict_sps2i,aclist,Jlist,signlist = string_data
    # unpacking
    result = np.matrix(np.zeros(Cx.shape))

    for clm in xrange(Cx.shape[1]):
        C0 = np.array(Cx[:,clm]).reshape(ns,ns)
        ###########################################################

        K = np.diag([np.einsum('i->',k_mtx[i,i]) for i in occ])
        for i,j in enumerate(Jlist): # j is a list
            K[i][j] = np.array([k_mtx[k] for k in aclist[i]]) * np.array(signlist[i])

        siga = np.dot(K,C0)
        sigb = np.dot(C0,K.T) # = np.dot(K,C0.T).T
        sig1 = siga + sigb

        ###########################################################
#        g_tmp = np.einsum('iijj->ij',g_mtx) + np.einsum('ijji'->'ij',g_mtx)


        #--------------
        sig21 = np.zeros([ns,ns])
        for i in xrange(ns):
            # leap 0, Ia=Ja
            sumgqqss = sum([g_mtx[q*no+q,s*no+s] for q in occ[i] for s in occ[i]])
            sumgsqqs = sum([g_mtx[s*no+q,q*no+s] for q in vir[i] for s in occ[i]])
            sig21[i,:] += (sumgqqss + sumgsqqs) * C0[i,:] # Ia
            sig21[:,i] += (sumgqqss + sumgsqqs) * C0[:,i] # Ib
            # leap 1,
        for i in xrange(ns):
            for p in occ[i]: 
                for q in vir[i]: #p->q
                    jstr = spstr[i] ^ (1 << p|1 << q) # XOR
                    j = dict_sps2i[jstr] # Epq
                    sgn = sign(spstr[i],jstr)
                    for r in occ[i]: #rrpq
                        sig21[i,:] += sgn * g_mtx[r*no+r,p*no+q] * C0[j,:] # Ia
                        sig21[:,i] += sgn * g_mtx[r*no+r,p*no+q] * C0[:,j] # Ib
                    for r in occ[i]: #pqrr
                        if r != p:
                            sig21[i,:] += sgn * g_mtx[p*no+q,r*no+r] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+q,r*no+r] * C0[:,j] # Ib
                        else:
                            sig21[i,:] += sgn * g_mtx[p*no+q,q*no+q] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+q,q*no+q] * C0[:,j] # Ib
                    for r in vir[i]: #prrq
                        if r != q:
                            sig21[i,:] += sgn * g_mtx[p*no+r,r*no+q] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+r,r*no+q] * C0[:,j] # Ib
                    for r in occ[i]: #rqpr
                        if r != p:
                            sig21[i,:] += - sgn * g_mtx[r*no+q,p*no+r] * C0[j,:] # Ia
                            sig21[:,i] += - sgn * g_mtx[r*no+q,p*no+r] * C0[:,j] # Ib
                        
            for p in occ[i]: #pq??
                for q in vir[i]:
                    jstr = spstr[i] ^ (1 << p|1 << q) # XOR
                    for r in occ[i]:
                        if r != p:
                            for s in vir[i]:
                                if s != q:
                                    jstr1 = jstr ^ (1 << r|1 << s) # XOR
                                    sgn = sign(spstr[i],jstr) * sign(jstr,jstr1)
                                    j = dict_sps2i[jstr1] # Epq##############
                                    sig21[i,:] += sgn * g_mtx[p*no+q,r*no+s] * C0[j,:] # Ia
                                    sig21[:,i] += sgn * g_mtx[p*no+q,r*no+s] * C0[:,j] # Ib
        sig21 *= 0.5
        #==================
        sig22 = np.zeros([ns,ns])
        
        for ia in xrange(ns):
            for ib in xrange(ns):
                #ia=ja,ib=jb
                for p in occ[ia]:
                    for r in occ[ib]:
                        sig22[ia,ib] += g_mtx[p*no+p,r*no+r] * C0[ia,ib]
                #ia=ja,ib!=jb
                for p in occ[ia]:
                    for r in occ[ib]:
                        for s in vir[ib]:
                            joc = spstr[ib] ^ (1 << r|1 << s) # XOR
                            j = dict_sps2i[joc] # Epq
                            sgn = sign(spstr[ib],joc)
                            sig22[ia,ib] += sgn * g_mtx[p*no+p,r*no+s] * C0[ia,j]
                #ia!=ja,ib=jb
                for p in occ[ia]:
                    for q in vir[ia]:
                        for r in occ[ib]:
                            joc = spstr[ia] ^ (1 << p|1 << q) # XOR
                            j = dict_sps2i[joc] # Epq
                            sgn = sign(spstr[ia],joc)
                            sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+r] * C0[j,ib]
                for p in occ[ia]:
                    for q in vir[ia]:
                        joc1 = spstr[ia] ^ (1 << p|1 << q) # XOR
                        ja = dict_sps2i[joc1] # Epq
                        sgna = sign(spstr[ia],joc1)
                        for r in occ[ib]:
                            for s in vir[ib]:
                                joc2 = spstr[ib] ^ (1 << r|1 << s)
                                jb = dict_sps2i[joc2] # Epq
                                sgn = sign(spstr[ib],joc2) * sgna
                                sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+s] * C0[ja,jb]
                                 
        sig = np.zeros([ns,ns])
        sig = sig1 + sig21 + sig22
    result[:,clm] = np.matrix(sig).reshape(ns*ns,-1)
    #-----------------------------------------------------------------------
    return result 
