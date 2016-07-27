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
g_mtx = numpy.zeros([no,]*4)
k_mtx = numpy.zeros([no,no])
for i in xrange(no):
    for j in xrange(no):
        for r in xrange(no):
            for s in xrange(no):
                g_mtx[i,r,s,j] = ao2mo.kernel(mol, (m.mo_coeff[:,i:i+1],m.mo_coeff[:,r:r+1],m.mo_coeff[:,s:s+1],m.mo_coeff[:,j:j+1]), compact=False)

print 'g matrix'
print g_mtx[0,0,0,0]
print g_mtx
g_all = ao2mo.kernel(mol, m.mo_coeff, compact=False)
print g_all

################################################################

for i in xrange(no):
    for j in xrange(no):
        for r in xrange(no):
            for s in xrange(no):
                if g_mtx[i,j,r,s] == g_all[no*i+j,no*r+s]:
                    print 'TRUE'
                else:
                    print 'FALSE'

