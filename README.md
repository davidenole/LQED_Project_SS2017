# LQED_Project_SS2017
Serial C++ implementation of a transition study in Lattice QED

Containing three files:
  - lattice.hpp
  - lattice.cpp
  - latticeTest.cpp
  
The first two contain the definition and implementation of the core functions needed in the simulation

The latter contains a test MCMC run of 10000 sweeps, can be called by

./executable n_of_sides_per_direction beta opt

opt = (double) -> cold start with links frozen at value opt

opt = R -> hot start (random links)
