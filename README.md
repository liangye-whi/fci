Full-CI programming
=====================

version 4.1
-----------

#### New Features

#### Bugs

- Still cannot run properly. Modification see 4.1.

Instructions
-----------

- main.py
    Molecular system information input. Calls fci_davidson.py to run davidson iteration.

- fci_davidson.py
    Initialization of the electron integrals and run the davidson diagonalization. Calls HC_MOC.py to calculate the matrix multiplication of HC.

- HC_MOC.py
    Minimal Operation-Count method for the multiplication of HC.

- Opr.py
    Some basic tools and functions.

