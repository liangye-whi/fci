'''

This is a program calculating the ground state energy of 
a diatom system with Full CI method.

'''
#-------------------------------------------
import numpy
from opr_E import constructZ, getOrder, opr_E, math_C

def Hc(C0,ne,no,ns,k_mtx,g_mtx,Z):
    while 1:
        D = numpy.zeros([no,no,ns,ns])
        for r in range(no):
            for s in range(no):
                for ka in range(ns):
                    for kb in range(ns):
                        for jb in range(ns):
    #                        print opr_E(ne,no,r,s,kb,jb,Z) 
                            D[r,s,ka,kb] += opr_E(ne,no,r,s,kb,jb,Z) * C0[ka,jb] 
                        for ja in range(ns):
                            D[r,s,ka,kb] += opr_E(ne,no,r,s,ka,ja,Z) * C0[ja,kb]
    #            print 'D(',r,',',s,') matrix'
    #            print D[r,s,:,:]

        G = numpy.zeros([no,no,ns,ns])
        for p in range(no):
            for q in range(no):
                for r in range(no):
                    for s in range(no):
                        G[p,q,:,:] += 0.5 * g_mtx[p,q,r,s] * D[r,s,:,:]
    #            print 'G(',p,',',q,') matrix'
    #            print G[p,q,:,:]

        sig2 = numpy.zeros([ns,ns])
        for ia in range(ns):
            for ib in range(ns):
                for p in range(no):
                    for q in range(no):
                        for ka in range(ns):
                            for kb in range(ns):
                                sig2[ia,ib] += opr_E(ne,no,r,s,ia,ka,Z)*G[p,q,ka,ib]
                                sig2[ia,ib] += opr_E(ne,no,r,s,ib,kb,Z)*G[p,q,ia,kb]

        sig1 = numpy.zeros([ns,ns])
        for ia in range(ns):
            for ib in range(ns):
                for p in range(no):
                    for q in range(no):
                        sig1[ia,ib] += k_mtx[p,q]*D[p,q,ia,ib]

        sig = numpy.zeros([ns,ns])
        sig = sig1 + sig2
        return sig
#    print 'sig'
#    print sig
#    print sigx
#    print A
