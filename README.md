Full-CI programming
=====================

Version 5.1
-----------

#### New Features (In Progress)

- Redefine the string data structure.

- Accelerate the matrix multipication with matrix operations.

- Rewrite the diagonal elements HD calc

#### Bugs

- Try using preconditioner.

Instructions
-----------

- **main.py**

    Molecular system information input. Calls **fci_davidson.py** to run davidson iteration.

- **fci_davidson.py**

    Initialization of the electron integrals and run the davidson diagonalization. Calls **HC_MOC.py** to calculate the matrix multiplication of HC.

- **HC_MOC.py**

    Minimal Operation-Count method for the multiplication of HC.

- **Opr.py**

    Some basic tools and functions.

