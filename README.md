Full-CI programming
=====================

Version 1.1
-----------

#### New Features

- Simplified Davidson.

- N-resolution works. MOC doesn't.

#### Bugs


Instructions
-----------

- **main.py**

    Molecular system information input. Calls **fci_davidson.py** to run davidson iteration.

- **fci_davidson.py**

    Initialization of the electron integrals and run the davidson diagonalization. Calls **HC_MOC.py** to calculate the matrix multiplication of HC.

- **HC_NR.py**

    N-Resolution method for the multiplication of HC.

- **HC_MOC.py**

    Minimal Operation-Count method for the multiplication of HC.

- **Opr.py**

    Some basic tools and functions.

