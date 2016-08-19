Full-CI programming
=====================

Version 5.4
-----------

#### New Features

- Redefine the string data structure.

- Accelerate the matrix multipication with matrix operations.

- Rewrite the diagonal elements HD calc.

- Detailed optimization. Faster Davidson.

- Now is able to calc Li2 6-31g level

- Separate the matrix construction to opr.py.

Instructions
-----------

- **main.py**

    Molecular system information input. Calls **davidson.py** to run davidson iteration.

- **davidson.py**

    Initialization of the electron integrals and run the davidson diagonalization. Calls **HC.py** to calculate the matrix multiplication of HC.

- **HC.py**

    Minimal Operation-Count method for the multiplication of HC.

- **Opr.py**

    Some basic tools and functions, constructing the string and matrix.

