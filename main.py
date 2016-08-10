from pyscf import gto, scf, fci
from fci_my import FCI

R = 1.1
mol = gto.M(atom='He 0 0 0; He 0 0 '+str(R), basis='6-31g')
nre = mol.energy_nuc()
myhf = scf.RHF(mol)
myhf.kernel()

print 'My FCI'
print('E(FCI) = %.12f' % (FCI(mol) + nre))

print 'Benchmark'
cisolver = fci.FCI(mol, myhf.mo_coeff)
print('E(FCI) = %.12f' % (cisolver.kernel()[0] + nre))


