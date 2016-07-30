'''

This is a program calculating the ground state energy of 
a diatom system with Full CI method.

'''
#-------------------------------------------
import numpy
from pyscf import gto, scf, ao2mo
from opr_E import constructZ, opr_E, math_C

R = 1.1
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
#print 'ns = ',ns
Z = constructZ(ne,no)
#print 'Z'
#print Z
print 'ne =', ne, 'no =', no
mocc = m.mo_coeff[:,m.mo_occ>0]
mvir = m.mo_coeff[:,m.mo_occ==0]
print 'mocc'
print mocc
print 'mvir'
print mvir

#print '################## h_mtx'
h_ao_mtx = scf.hf.get_hcore(mol)
h_mtx = numpy.dot(numpy.dot(m.mo_coeff.T,h_ao_mtx),m.mo_coeff)
print 'h matrix'
print h_mtx
#eri = ao2mo.kernel(mol, (m.mo_coeff,)*4, compact=False)
#eri = ao2mo.kernel(mol, m.mo_coeff, compact=False)
#print 'eri'
#print(eri.shape)
#print eri
g_mtx = numpy.zeros([no,]*4)
k_mtx = numpy.zeros([no,no])
for i in range(no):
    for j in range(no):
        for r in range(no):
            for s in range(no):
                g_mtx[i,r,s,j] = ao2mo.kernel(mol, (m.mo_coeff[:,i:i+1],m.mo_coeff[:,r:r+1],m.mo_coeff[:,s:s+1],m.mo_coeff[:,j:j+1]), compact=False)

            k_mtx[i,j] -= g_mtx[i,r,r,j]
        k_mtx[i,j] *= 0.5
k_mtx = k_mtx + h_mtx
print 'g matrix'
print g_mtx[0,0,0,0]
#print '################## k_mtx'
#print k_mtx

################################################################
crt = 1e-9

C0 = numpy.ones([ns,ns])
#C0x = C0 ** 2
#A = sum(sum(C0x))
#C0 = C0 / (A**0.5)
C0 = C0 / ns
print 'C0'
print C0

E0 = 0 
#H0 = numpy.zeros([ns,ns,ns,ns]) #H0[kakb,jajb]
#for i in range(ns):
#    for j in range(ns):
#        H0[i,j,i,j] = 1
######################################################################
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
            for ka in range(ns):
                for kb in range(ns):
                    for r in range(no):
                        for s in range(no):
                            G[p,q,ka,kb] += 0.5 * g_mtx[p,q,r,s] * D[r,s,ka,kb]
#            print 'G(',p,',',q,') matrix'
#            print G[p,q,:,:]

    sig2 = numpy.zeros([ns,ns])
    for ia in range(ns):
        for ib in range(ns):
            for p in range(no):
                for q in range(no):
                    for ka in range(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ia,ka,Z) * G[p,q,ka,ib]
                    for kb in range(ns):
                        sig2[ia,ib] += opr_E(ne,no,p,q,ib,kb,Z) * G[p,q,ia,kb]

    sig1 = numpy.zeros([ns,ns])
    for ia in range(ns):
        for ib in range(ns):
            for p in range(no):
                for q in range(no):
                    sig1[ia,ib] += k_mtx[p,q] * D[p,q,ia,ib]

    sig = numpy.zeros([ns,ns])
    sig = sig1 + sig2
#    print 'sig'
#    print sig

    #####################################################################################
    E0 = sum(sum(C0 * sig))
    cdav = - (1 - E0)**(-1) * (sig - E0 * C0)
    res = sum(sum(cdav**2))
    print 'new iter =============='
    print 'E =', E0
    print '|res| =', res
    if res < crt:
        print 'Iteration Converged.'
        break
    C1 = C0 + cdav
    print 'Cnew'
#    print C1
    C1x = C1 ** 2
    A = sum(sum(C1x))
#    print 'normalization factor = ', A
    C0 =  C1 / (A**0.5)
    print C0
#    C1x = C1 ** 2
#    A = sum(sum(C1x))
#    print 'normalization factor =', A
    
#    E1 = sum(sum(C1 * sig))
#    E = sum(sum((C1**2) * e_mtx))
#    print ' E =', E1
#    if E0 == 0:
#        E0 = E1
#        C0 = C1
#        continue

#    print 'dE =', E0-E1
#    if abs(E1-E0) > crt:
#        E0 = E1
#        C0 = C1
#    else:
#        break
######################################################################################
print 'groundstate energy'
print E0 + mol.energy_nuc()
