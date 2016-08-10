Full-CI programming
=====================

Version 4.1
-----------

#### New Features

- Rewrite the Davidson algorithm into a FULL style according to his original paper.

- Result correct and run fast.

#### Bugs

- The perturbation should be added to initial guess of basis, but not the diagonal H0.

- The davidson iteration may still have a little problem, because the iteration steps are too many during calc He2 CC-PVDZ system (24 iter).

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

