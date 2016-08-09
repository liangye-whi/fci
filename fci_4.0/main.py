from pyscf import gto, scf, fci
from fci_davidson import FCI

R = 1.1
mol = gto.Mole()
mol.build(
        atom = 'He 0 0 0; He 0 0 '+str(R), 
        basis = '6-31G',
        symmetry = True,
        )

print 'My FCI'
print('E(FCI) = %.12f' % (FCI(mol) + mol.energy_nuc()))

print 'Benchmark'
myhf = scf.RHF(mol)
myhf.kernel()
cisolver = fci.FCI(mol, myhf.mo_coeff)
print('E(FCI) = %.12f' % (cisolver.kernel()[0] + mol.energy_nuc()))


