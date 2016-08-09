'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from opr import formOccu, sign

def HC(Cx, k_mtx, g_mtx, ne, no, ns, Z, occu):
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns
    #print occu
    result = numpy.matrix(numpy.zeros([ns*ns,Cx.shape[1]]))
    for clm in xrange(Cx.shape[1]):
        C0 = numpy.array(Cx[:,clm]).reshape(ns,ns)
    ###########################################################
        sig1 = numpy.zeros([ns,ns])
        for i in xrange(ns):
            # diagonal elem
            res='0'*(no+2-len(bin(occu[0][i])))+bin(occu[0][i])[2:]
            cu = [k for k in xrange(no) if res[k] == '1']
            sumkpp = sum([k_mtx[p,p] for p in cu])
            sig1[i,:] += sumkpp * C0[i,:] # Ia
            sig1[:,i] += sumkpp * C0[:,i] # Ib
            
            # nondiagonal elem
            ucu = [k for k in xrange(no) if res[k] == '0']
            for p in cu:
                for q in ucu:
                    assert p != q ##debug##
                    joc = occu[0][i] ^ ((1 << (no-q-1)) + (1 << (no-p-1))) # XOR
                    j = occu[1][str(joc)] # Epq
                    sgn = sign(res,p,q)
                    sig1[i,:] += sgn * k_mtx[p,q] * C0[j,:] # Ia
                    sig1[:,i] += sgn * k_mtx[p,q] * C0[:,j] # Ib
        
        ###########################################################

        sig21 = numpy.zeros([ns,ns])
        for i in xrange(ns):
            res='0'*(no+2-len(bin(occu[0][i])))+bin(occu[0][i])[2:]
            cu = [k for k in xrange(no) if res[k] == '1']
            ucu = [k for k in xrange(no) if res[k] == '0']
            # leap 0, Ia=Ja
            sumgqqss = sum([g_mtx[q*no+q,s*no+s] for q in cu for s in cu])
            sumgsqqs = sum([g_mtx[s*no+q,q*no+s] for q in ucu for s in cu])
            sig21[i,:] += (sumgqqss + sumgsqqs) * C0[i,:] # Ia
            sig21[:,i] += (sumgqqss + sumgsqqs) * C0[:,i] # Ib
            # leap 1,
        for i in xrange(ns):
            res='0'*(no+2-len(bin(occu[0][i])))+bin(occu[0][i])[2:]
            cu = [k for k in xrange(no) if res[k] == '1']
            ucu = [k for k in xrange(no) if res[k] == '0']
            for p in cu: 
                for q in ucu: #p->q
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                    j = occu[1][str(joc)] # Epq
                    sgn = sign(res,p,q)
                    for r in cu: #rrpq
                        sig21[i,:] += sgn * g_mtx[r*no+r,p*no+q] * C0[j,:] # Ia
                        sig21[:,i] += sgn * g_mtx[r*no+r,p*no+q] * C0[:,j] # Ib
                    for r in cu: #pqrr
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
                    for r in cu: #rqpr
                        if r != p:
                            sig21[i,:] += - sgn * g_mtx[r*no+q,p*no+r] * C0[j,:] # Ia
                            sig21[:,i] += - sgn * g_mtx[r*no+q,p*no+r] * C0[:,j] # Ib
                        
            for p in cu: #pq??
                for q in ucu:
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                    resj='0'*(no+2-len(bin(joc)))+bin(joc)[2:]
                    for r in cu:
                        if r != p:
                            for s in ucu:
                                if s != q:
                                    sgn = sign(res,p,q) * sign(resj,r,s)
                                    joc1 = joc ^ ((1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                                    j = occu[1][str(joc1)] # Epq##############
                                    sig21[i,:] += sgn * g_mtx[p*no+q,r*no+s] * C0[j,:] # Ia
                                    sig21[:,i] += sgn * g_mtx[p*no+q,r*no+s] * C0[:,j] # Ib
        sig21 *= 0.5
        #==================
        sig22 = numpy.zeros([ns,ns])
        
        for ia in xrange(ns):
            res1='0'*(no+2-len(bin(occu[0][ia])))+bin(occu[0][ia])[2:]
            cu1 = [k for k in xrange(no) if res1[k] == '1']
            ucu1 = [k for k in xrange(no) if res1[k] == '0']
            for ib in xrange(ns):
                res2='0'*(no+2-len(bin(occu[0][ib])))+bin(occu[0][ib])[2:]
                cu2 = [k for k in xrange(no) if res2[k] == '1']
                ucu2 = [k for k in xrange(no) if res2[k] == '0']
                #ia=ja,ib=jb
                for p in cu1:
                    for r in cu2:
                        sig22[ia,ib] += g_mtx[p*no+p,r*no+r] * C0[ia,ib]
                #ia=ja,ib!=jb
                for p in cu1:
                    for r in cu2:
                        for s in ucu2:
                            joc = occu[0][ib] ^ ((1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                            j = occu[1][str(joc)] # Epq
                            sgn = sign(res2,r,s)
                            sig22[ia,ib] += sgn * g_mtx[p*no+p,r*no+s] * C0[ia,j]
                #ia!=ja,ib=jb
                for p in cu1:
                    for q in ucu1:
                        for r in cu2:
                            joc = occu[0][ia] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                            j = occu[1][str(joc)] # Epq
                            sgn = sign(res1,p,q)
                            sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+r] * C0[j,ib]
                for p in cu1:
                    for q in ucu1:
                        joc1 = occu[0][ia] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                        ja = occu[1][str(joc1)] # Epq
                        sgna = sign(res1,p,q)
                        for r in cu2:
                            for s in ucu2:
                                joc2 = occu[0][ib] ^ ((1 << (no-r-1)) + (1 << (no-s-1)))
                                jb = occu[1][str(joc2)] # Epq
                                resjb='0'*(no+2-len(bin(joc2)))+bin(joc2)[2:]
                                sgn = sign(res2,r,s) * sgna
                                sig22[ia,ib] += sgn * g_mtx[p*no+q,r*no+s] * C0[ja,jb]
                                 
        sig = numpy.zeros([ns,ns])
        sig = sig1 + sig21 + sig22
    result[:,clm] = numpy.matrix(sig).reshape(ns*ns,-1)
    #-----------------------------------------------------------------------
    return result 
