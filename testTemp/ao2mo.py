#/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

from pyscf import gto, scf, ao2mo

'''
A simple example to call integral transformation for given orbitals
'''

mol = gto.Mole()
mol.build(
    atom = 'H 0 0 0; H 0 0 1.1',  # in Angstrom
    #basis = 'sto-3g',
    basis = '6-31g',
    symmetry = True,
)

myhf = scf.RHF(mol)
myhf.kernel()

orb = myhf.mo_coeff
eri_4fold = ao2mo.kernel(mol, orb, compact=False)
#eri_4fold = ao2mo.general(mol, orb)
print('MO integrals (ij|kl) with 4-fold symmetry i>=j, k>=l have shape %s' %
      str(eri_4fold.shape))
#print(eri_4fold)
