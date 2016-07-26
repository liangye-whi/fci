from pyscf import gto, scf, ao2mo

#mol = gto.M(atom='H 0 0 0; H 0 0 1.1', basis='cc-pvdz')
mol = gto.M(atom='H 0 0 0; H 0 0 1.1', basis='sto-3g')

m = scf.RHF(mol)
m.kernel()

print(m.mo_coeff)
print(m.mo_occ)
print(m.make_rdm1())
j, k = scf.hf.get_jk(mol,m.make_rdm1())
#eri = ao2mo.incore.general(m._eri, (m.mo_coeff,)*4, compact=False)
eri = ao2mo.kernel(mol,m.mo_coeff,compact=False)

print(len(eri))
print(eri.shape)
print(eri)
