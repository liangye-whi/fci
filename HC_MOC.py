'''

This is a program calculating the ground state energy of 
a small system with Full CI method.

Published Edtion.

by Ye @ 26JUL2016.

'''
#-------------------------------------------

import numpy
from opr import formOccu

def HC(C0, k_mtx, g_mtx, ne, no, ns, Z, occu):
    #print 'Z'
    #print Z
    #print 'ne =', ne, 'no =', no, 'ns =', ns
    #print occu
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
                sig1[i,:] += k_mtx[p,q] * C0[j,:] # Ia
                sig1[:,i] += k_mtx[p,q] * C0[:,j] # Ib
    
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
        for p in cu:
            for r in cu:
                for s in ucu:
                    joc = occu[0][i] ^ ((1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                    j = occu[1][str(joc)]
                    sig21[i,:] += g_mtx[p*no+p,r*no+s] * C0[j,:]
                    sig21[:,i] += g_mtx[p*no+p,r*no+s] * C0[:,j]
        for p in cu:
            for r in cu:
                for q in ucu:
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1))) # XOR
                    j = occu[1][str(joc)]
                    sig21[i,:] += g_mtx[p*no+q,r*no+r] * C0[j,:]
                    sig21[:,i] += g_mtx[p*no+q,r*no+r] * C0[:,j]
        for p in cu:
            for q in ucu:
                for s in ucu:
                    joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-s-1))) # XOR
                    j = occu[1][str(joc)]
                    sig21[i,:] += g_mtx[p*no+q,q*no+s] * C0[j,:]
                    sig21[:,i] += g_mtx[p*no+q,q*no+s] * C0[:,j]
        for p in cu:
            for q in ucu:
                for r in cu:
                    joc = occu[0][i] ^ ((1 << (no-r-1)) + (1 << (no-q-1))) # XOR
                    j = occu[1][str(joc)]
                    sig21[i,:] += g_mtx[p*no+q,r*no+p] * C0[j,:]
                    sig21[:,j] += g_mtx[p*no+q,r*no+p] * C0[:,j]
        # leap 2,
        for p in cu:
            for q in ucu:
                for r in cu:
                    if p != r:
                        for s in ucu:
                            if q != s:
                                joc = occu[0][i] ^ ((1 << (no-p-1)) + (1 << (no-q-1)) + (1 << (no-r-1)) + (1 << (no-s-1))) # XOR
                                j = occu[1][str(joc)]
                                sig21[i,:] += g_mtx[p*no+q,r*no+s] * C0[j,:]
                                sig21[:,j] += g_mtx[p*no+q,r*no+s] * C0[:,j]
    sig21 *= 0.5

    sig22 = numpy.zeros([ns,ns])
    D = numpy.zeros([ns,ns,no,no])
    G = numpy.zeros([ns,ns,no,no])
    for ia in xrange(ns):
        res='0'*(no+2-len(bin(occu[0][ia])))+bin(occu[0][ia])[2:]
        cu = [k for k in xrange(no) if res[k] == '1']
        ucu = [k for k in xrange(no) if res[k] == '0']
        for jb in xrange(ns):
            for p in cu:
                D[ia,jb,p,p] += ne * C0[ia,jb]
                for q in ucu:
                    joc = occu[0][ia] ^ ((1 << (no-p-1)) + (1 << (no-q-1)))
                    ja = occu[1][str(joc)]
                    D[ia,jb,p,q] += C0[ja,jb]
    for ib in xrange(ns):
        res='0'*(no+2-len(bin(occu[0][ib])))+bin(occu[0][ib])[2:]
        cu = [k for k in xrange(no) if res[k] == '1']
        ucu = [k for k in xrange(no) if res[k] == '0']
        for r in cu:
            for p in xrange(no):
                for q in xrange(no):
                    G[ib,ib,p,q] += g_mtx[p*no+q,r*no+r]
            for s in ucu:
                joc = occu[0][ib] ^ ((1 << (no-r-1)) + (1 << (no-s-1)))
                jb = occu[1][str(joc)]
                for p in xrange(no):
                    for q in xrange(no):
                        G[ib,jb,p,q] += C0[ja,jb]
                        G[ia,jb,p,q] += g_mtx[p*no+q,r*no+s]
    sig = numpy.zeros([ns,ns])
    sig = sig1 + sig21 + sig22
    #-----------------------------------------------------------------------
    return sig
