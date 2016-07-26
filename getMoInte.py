from pyscf import gto, scf, ao2mo

#mol = gto.M(atom='H 0 0 0; H 0 0 1.1', basis='cc-pvdz')
mol = gto.M(atom='H 0 0 0; H 0 0 1.1', basis='sto-3g')

m = scf.RHF(mol)
m.kernel()
print 'mo_coeff'
print(m.mo_coeff)
print 'mo_occ'
print(m.mo_occ)
print 'mo_energy'
print(m.mo_energy)
#print(m.make_rdm1())
#j, k = scf.hf.get_jk(mol,m.make_rdm1())
#eri = ao2mo.incore.general(m._eri, (m.mo_coeff,)*4, compact=False)
#eri = ao2mo.kernel(mol,m.mo_coeff,compact=False)

import h5py
import numpy
mocc = m.mo_coeff[:,m.mo_occ>0]
mvir = m.mo_coeff[:,m.mo_occ==0]
print 'mocc'
print mocc
print 'mvir'
print mvir

ao2mo.general(mol, (mocc,mvir,mocc,mvir), 'tmp.h5', compact=False)
feri = h5py.File('tmp.h5')
ovvv = numpy.array(feri['eri_mo'])
print(ovvv.shape)
print ovvv
'''
print(len(eri))
print(eri.shape)
print(eri)
'''
