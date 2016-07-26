'''

This is a program calculating the ground state energy of 
a diatom system with Full CI method.

At present it can only do primitive diagonalizaion of the
Hamiltonian.

'''
import scipy
import numpy as np
from scipy import linalg

a = 1.889726133
R = 1.1*a #angstrom
nre = 1/R
print 'nre'
print nre

e1 = -0.45421869
e2 =  0.39669591
J11 = 0.60917168
J22 = 0.63747988
J12 = 0.60733543
K12 = 0.20322223

print 'Hamiltonian'
a = np.array([[2*e1-J11,K12],[K12,2*e2-4*J12+J22+2*K12]])
print a

x, y = linalg.eig(a)
print 'eigenvalue'
print x
print 'eigenvector'
print y
print 'groundstate e'
print x[0]+nre
print '--------'

d = e2-e1+J11/2.0+J22/2.0-2*J12+K12
ee = 2*(e1-J11)+J11+d-(d**2+K12**2)**0.5
print ee
