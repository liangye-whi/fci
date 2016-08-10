Full-CI programming
=====================

Version 0.1
-----------

#### New Features

- Using simplified Davidson steps and N-Resolution method for HC calc.

- Result correct only for H2 sto-3g and 6-31g.

#### Bugs

- The sign of exchange operators is incorrected. See version 2.0. 

- NR method is too slow.

- The initial guess is coarse.

Instructions
-----------

- **main.py**

    Molecular system information input. Calls **fci_davidson.py** to run davidson iteration.

- **fci_my.py**

    Initialization of the electron integrals and run the davidson diagonalization. 
