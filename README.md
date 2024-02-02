# Quantum Dots
## Parameterization of quantum dots using machine learning

This project works on the parameterization of a system of quantum dots.

The project consists of four programs

occup.f90 is the module for the theoretical simulation of the model, which is used to write the database.

h5main4.f90 is the program that creates and writes the database using the hdf5 library and calling functions from the occup.f90 module

The decision to write this database in Fortran is due to the computational speed of the language.

dataset.ypynb has the same theoretical model implemented as occup.90. Although it was not used to create the database, this program is useful when viewing different cases and being able to view them.

double_dot.ipynb is the program that imports the database, creates the architecture of the Neural Network and proceeds with its training.
