'''

This is a program calculating the ground state energy of 
a diatom system with Full CI method.

At present it can only do primitive diagonalizaion of the
Hamiltonian.

'''
#-------------------------------------------
import h5py
import numpy
from pyscf import gto, scf, ao2mo

R=1.1
mol = gto.M(atom='H 0 0 0; H 0 0 '+str(R), basis='sto-3g')

m = scf.RHF(mol)
m.kernel()

print 'mo_coeff'
print(m.mo_coeff)
print 'mo_occ'
print(m.mo_occ)
print 'mo_energy'
print(m.mo_energy)

mocc = m.mo_coeff[:,m.mo_occ>0]
mvir = m.mo_coeff[:,m.mo_occ==0]
print 'mocc'
print mocc
print 'mvir'
print mvir

ao2mo.general(mol, (mocc,mocc,mocc,mocc), 'tmp.h5', compact=False)
feri = h5py.File('tmp.h5')
J11 = numpy.array(feri['eri_mo'])
print 'J11'
print(J11.shape)
print J11

ao2mo.general(mol, (mocc,mocc,mvir,mvir), 'tmp.h5', compact=False)
feri = h5py.File('tmp.h5')
J12 = numpy.array(feri['eri_mo'])
print 'J12'
print(J12.shape)
print J12

ao2mo.general(mol, (mvir,mvir,mvir,mvir), 'tmp.h5', compact=False)
feri = h5py.File('tmp.h5')
J22 = numpy.array(feri['eri_mo'])
print 'J22'
print(J22.shape)
print J22

ao2mo.general(mol, (mocc,mvir,mocc,mvir), 'tmp.h5', compact=False)
feri = h5py.File('tmp.h5')
K12 = numpy.array(feri['eri_mo'])
print 'K12'
print(K12.shape)
print K12
#----------------------------------------------------
import scipy
#import numpy 
from scipy import linalg

a = 1.889726133
Ra = R*a #angstrom
nre = 1/Ra
print 'nre'
print nre

print 'Hamiltonian'
a = numpy.array([[2*m.mo_energy[0]-J11[0,0],        K12[0,0]],
    [K12[0,0],      2*m.mo_energy[1]-4*J12[0,0]+J22[0,0]+2*K12[0,0]]])
print a

x, y = linalg.eig(a)
print 'eigenvalue'
print x
print 'eigenvector'
print y
print 'groundstate energy'
print x[0].real+nre
