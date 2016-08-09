from pyscf import gto, scf, fci
from fci_davidson import FCI

R = 1.1

mol = gto.Mole()
mol.verbose = 0
mol.output = None
mol.atom = [
    ['H', ( 1.,-1.    , 0.   )],
    ['H', ( 0.,-1.    ,-1.   )],
    ['H', ( 1.,-0.5   ,-1.   )],
    ['H', ( 0.,-0.    ,-1.   )],
    ['H', ( 1.,-0.5   , 0.   )],
    ['H', ( 0., 1.    , 1.   )],
]
mol.basis = 'sto-3g'
mol.build()

print 'My FCI'
print('E(FCI) = %.12f' % (FCI(mol) + mol.energy_nuc()))

print 'Benchmark'
myhf = scf.RHF(mol)
myhf.kernel()
cisolver = fci.FCI(mol, myhf.mo_coeff)
print('E(FCI) = %.12f' % (cisolver.kernel()[0] + mol.energy_nuc()))


