'''

This is a program calculating the ground state energy of 
a diatom system with Full CI method.

'''
#-------------------------------------------
import numpy
from pyscf import gto, scf, ao2mo
from opr_E import constructZ, getOrder, opr_E, math_C

R=1.1
mol = gto.M(atom='H 0 0 0; H 0 0 '+str(R), basis='6-31g')

m = scf.RHF(mol)
m.kernel()

print 'mo_coeff'
print(m.mo_coeff)
print 'mo_occ'
print(m.mo_occ)
print 'mo_energy'
print(m.mo_energy)

ne = mol.nelectron/2 #electron per string
no = len(m.mo_energy)
ns = math_C(no,ne)
Z = constructZ(ne,no)
#print 'Z'
#print Z
print 'ne =', ne, 'no =', no, 'ns =', ns

###########################################################
h_ao_mtx = scf.hf.get_hcore(mol)
h_mtx = numpy.dot(numpy.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
print 'h matrix'
print h_mtx
#eri = ao2mo.kernel(mol, (m.mo_coeff,)*4, compact=False)
#eri = ao2mo.kernel(mol, m.mo_coeff, compact=False)
#print 'eri'
#print(eri.shape)
#print eri
g_mtx = ao2mo.kernel(mol, m.mo_coeff, compact=False)
k_mtx = numpy.zeros([no,no])
for i in xrange(no):
    for j in xrange(no):
        for r in xrange(no):
            k_mtx[i,j] -= g_mtx[i*no+r,r*no+j]
        k_mtx[i,j] *= 0.5
k_mtx = k_mtx + h_mtx
print 'g matrix'
print g_mtx[0,0]
################################################################
crt = 1e-9

C0 = numpy.ones([ns,ns])
#C0[1,0] , C0[0,1] , C0[1,1]= 0, 0, 0.4
C0x = C0 ** 2
A = sum(sum(C0x))
C0 =  C0 / (A**0.5)
print 'C0'
print C0

E0 = 0 
#H0 = numpy.zeros([ns,ns,ns,ns]) #H0[kakb,jajb]
#for i in xrange(ns):
#    for j in xrange(ns):
#        H0[i,j,i,j] = 1
######################################################################
while 1:
    D = numpy.zeros([no,no,ns,ns])
    for r in xrange(no):
        for s in xrange(no):
            for k in xrange(ns):
                for j in xrange(ns):
#                    print opr_E(ne,no,r,s,kb,jb,Z) 
                    D[r,s,:,k] += opr_E(ne,no,r,s,k,j,Z) * C0[:,j] 
                    D[r,s,k,:] += opr_E(ne,no,r,s,k,j,Z) * C0[j,:]
#            print 'D(',r,',',s,') matrix'
#            print D[r,s,:,:]

    G = numpy.zeros([no,no,ns,ns])
    for p in xrange(no):
        for q in xrange(no):
            for r in xrange(no):
                for s in xrange(no):
                    G[p,q,:,:] += 0.5 * g_mtx[p*no+q,r*no+s] * D[r,s,:,:]
#            print 'G(',p,',',q,') matrix'
#            print G[p,q,:,:]

    sig2 = numpy.zeros([ns,ns])
    for ia in xrange(ns):
        for ib in xrange(ns):
            for p in xrange(no):
                for q in xrange(no):
                    for ka in xrange(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ia,ka,Z)*G[p,q,ka,ib]
                    for kb in xrange(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ib,kb,Z)*G[p,q,ia,kb]

    sig1 = numpy.zeros([ns,ns])
    for ia in xrange(ns):
        for ib in xrange(ns):
            for p in xrange(no):
                for q in xrange(no):
                    sig1[ia,ib] += k_mtx[p,q]*D[p,q,ia,ib]

    sig = numpy.zeros([ns,ns])
    sig = sig1 + sig2
#    print 'sig'
#    print sig
#    print sigx
#    print A

#####################################################################################
    cdav = - (1 - E0)**(-1) * (sig - E0 * C0)
    C1 = C0 + cdav
    print 'Cnew  ===================='
#    print C1
    C1x = C1 ** 2
    A = sum(sum(C1x))
#    print 'normalization factor = ', A
    C1 =  C1 / (A**0.5)
    print C1
#    C1x = C1 ** 2
#    A = sum(sum(C1x))
#    print 'normalization factor =', A
    
    E1 = sum(sum(C1 * sig))
#    E = sum(sum((C1**2) * e_mtx))
    print ' E =', E1
    if E0 == 0:
        E0 = E1
        C0 = C1
        continue

    print 'dE =', E0-E1
    if abs(E1-E0) > crt:
        E0 = E1
        C0 = C1
    else:
        break
######################################################################################
a = 1.889726133
Ra = R*a #angstrom
nre = 1/Ra
print 'nre'
print nre
print 'groundstate energy'
print E1+nre
