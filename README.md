Full-CI programming
=====================

Version 5.0 alpha phase 2_
-----------

#### New Features (In Progress)

- Redefine the string logic (simplify).

    **Done.**

- Think about accelerate the matrix multipication.

- The slowest part of diagonal elements calc!!

    **sig1 & matrix K done.**

- Take care of the davidson steps.

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

