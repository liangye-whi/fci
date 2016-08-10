'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

!!!! C0 is a (ns*ns, 1) matrix now.

'''
#-------------------------------------------

import numpy
from opr import formOccu, sign

def HC(Cx, k_mtx, g_mtx, ne, no, ns, Z, occu):
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns
    #print occu
    rst = numpy.matrix(numpy.zeros([ns*ns,Cx.shape[1]]))
    for clm in xrange(Cx.shape[1]):

        C0 = Cx[:,clm]
        #print Cx.shape,C0.shape
        C0 = C0.reshape([ns,ns])
    ###########################################################
        sig1 = numpy.matrix(numpy.zeros([ns,ns]))
        for i in xrange(ns):
            # diagonal elem
            sumkpp = sum(numpy.diag(k_mtx)[occu[2][i]])
            sig1[i,:] += sumkpp * C0[i,:] # Ia
            sig1[:,i] += sumkpp * C0[:,i] # Ib
            
            # nondiagonal elem
            for p in occu[2][i]:
                for q in occu[3][i]:
                    assert p != q ##debug##
                    joc = occu[0][i] ^ ((1 << (no-q-1)) + (1 << (no-p-1))) # XOR
                    j = occu[1][str(joc)] # Epq
                    sgn = sign(occu[2][i],p,q)
                    sig1[i,:] += sgn * k_mtx[p,q] * C0[j,:] # Ia
                    sig1[:,i] += sgn * k_mtx[p,q] * C0[:,j] # Ib
        
        ###########################################################

        sig21 = numpy.matrix(numpy.zeros([ns,ns]))
        for i in xrange(ns):
            # leap 0, Ia=Ja
            sumgqqss = sum([g_mtx[q*no+q,s*no+s] for q in occu[2][i] for s in occu[2][i]])
            sumgsqqs = sum([g_mtx[s*no+q,q*no+s] for q in occu[3][i] for s in occu[2][i]])
            sig21[i,:] += (sumgqqss + sumgsqqs) * C0[i,:] # Ia
            sig21[:,i] += (sumgqqss + sumgsqqs) * C0[:,i] # Ib
            # leap 1,
        for i in xrange(ns):
            ocu = occu[2][i]
            ucu = occu[3][i]
            for p in ocu: 
                for q in ucu: #p->q
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                    j = occu[1][str(joc)] # Epq
                    sgn = sign(ocu,p,q)
                    for r in ocu: #rrpq
                        sig21[i,:] += sgn * g_mtx[r*no+r,p*no+q] * C0[j,:] # Ia
                        sig21[:,i] += sgn * g_mtx[r*no+r,p*no+q] * C0[:,j] # Ib
                    for r in ocu: #pqrr
                        if r != p:
                            sig21[i,:] += sgn * g_mtx[p*no+q,r*no+r] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+q,r*no+r] * C0[:,j] # Ib
                        else:
                            sig21[i,:] += sgn * g_mtx[p*no+q,q*no+q] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+q,q*no+q] * C0[:,j] # Ib
                    for r in ucu: #prrq
                        if r != q:
                            sig21[i,:] += sgn * g_mtx[p*no+r,r*no+q] * C0[j,:] # Ia
                            sig21[:,i] += sgn * g_mtx[p*no+r,r*no+q] * C0[:,j] # Ib
                    for r in ocu: #rqpr
                        if r != p:
                            sig21[i,:] += - sgn * g_mtx[r*no+q,p*no+r] * C0[j,:] # Ia
                            sig21[:,i] += - sgn * g_mtx[r*no+q,p*no+r] * C0[:,j] # Ib
                        
            for p in ocu: #pq??
                for q in ucu:
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                    resj='0'*(no+2-len(bin(joc)))+bin(joc)[2:]
                    for r in ocu:
                        if r != p:
                            for s in ucu:
                                if s != q:
                                    sgn = sign(ocu,p,q) * sign(resj,r,s)
                                    joc1 = joc ^ ((1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                                    j = occu[1][str(joc1)] # Epq##############
                                    sig21[i,:] += sgn * g_mtx[p*no+q,r*no+s] * C0[j,:] # Ia
                                    sig21[:,i] += sgn * g_mtx[p*no+q,r*no+s] * C0[:,j] # Ib
        sig21 *= 0.5
        #==================
        sig22 = numpy.matrix(numpy.zeros([ns,ns]))
        
        for ia in xrange(ns):
            ocu1 = occu[2][ia] 
            ucu1 = occu[3][ia]
            for ib in xrange(ns):
                ocu2 = occu[2][ib]
                ucu2 = occu[3][ib]
                #ia=ja,ib=jb
                for p in ocu1:
                    for r in ocu2:
                        sig22[ia,ib] += g_mtx[p*no+p,r*no+r] * C0[ia,ib]
                #ia=ja,ib!=jb
                for p in ocu1:
                    for r in ocu2:
                        for s in ucu2:
                            joc = occu[0][ib] ^ ((1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                            j = occu[1][str(joc)] # Epq
                            sgn = sign(ocu2,r,s)
                            sig22[ia,ib] += sgn * g_mtx[p*no+p,r*no+s] * C0[ia,j]
                #ia!=ja,ib=jb
                for p in ocu1:
                    for q in ucu1:
                        for r in ocu2:
                            joc = occu[0][ia] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                            j = occu[1][str(joc)] # Epq
                            sgn = sign(ocu1,p,q)
                            sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+r] * C0[j,ib]
                for p in ocu1:
                    for q in ucu1:
                        joc1 = occu[0][ia] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                        ja = occu[1][str(joc1)] # Epq
                        sgna = sign(ocu1,p,q)
                        for r in ocu2:
                            for s in ucu2:
                                joc2 = occu[0][ib] ^ ((1 << (no-r-1)) + (1 << (no-s-1)))
                                jb = occu[1][str(joc2)] # Epq
                                sgn = sign(ocu2,r,s) * sgna
                                sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+s] * C0[ja,jb]
                                 
                sig = numpy.zeros([ns,ns])
                sig = sig1 + sig21 + sig22
            rst[:,clm] = sig.reshape(-1).T
    #-----------------------------------------------------------------------
    return rst
