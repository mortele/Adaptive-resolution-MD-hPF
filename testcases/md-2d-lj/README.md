### Test case: 2D Lennard-Jones fluid
This test case is testing the MD part of this code directly against LAMMPS output, assuming the latter to be verified correct. 

Reproducing the plots in `/figures` is done by running the python scripts in `/analysis`. A compiled (serial) LAMMPS binary called `lammps` is assumed to be in the path.