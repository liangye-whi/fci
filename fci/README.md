
program of Full CI
==================

To set the molecule, modify the file `main.py` and python it to run.

the submodules of the program:

- fci_davidson
    main.py call this module to run a davidson iteration to solve the groundstate problem.

- HC_MOC
    inside the davidson iteration fci_davidson call this procedure to calculate the $\sigma = Hc$ multiplication. MOC refers to the algorithm called the Minimal Operation-Count method introduced by the purple book.

- opr
    this package includes the operation functions such as constructing the Shacitt map and the dictionary of the index of any specific spin orbital string.
