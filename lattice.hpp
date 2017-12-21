#include <string>

#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

// returns the total dimension of the lattice
//		line 17
int dTot(int*);

// returns the "1D" position on the lattice
//		line 23	
int myPos(int*, int, int, int, int);
//  overloading for array of int
//		line 29
int myPos(int*, int*);

// initialisation of the lattice to constant
//		line 37
void initLattice(int*, double*, double);

// initialisation of the lattice to random in [-pi,pi]
//		line 53
void initLattice(int*, double*);

// copy initialisation (given another lattice)
//		line 69
void initLattice(int*, double*, double*);

// prints the links to std output
//		line 80
void printL(int*, double*);

// prints the links on file
//		line 104
void printL(int*, double*, std::string);

// prints an observable to a file
//		line 130
void printO(int, double*, std::string);

// saves current configuration in the given position (already initialized)
//		line 141
void saveConf(int*, double*, double*);

// loads the configuration from a saved one
//		line 148
void loadConf(int*, double*, double*);

// evaluates the plaquette action in the wanted directions 
//		line 155
double lPAction(int*, double*, int, int, int, int, int, int);
// overloading for input with an array of int
//		line 175
double lPAction(int*, double*, int, int, int*);

// evaluates the QED action on the whole lattice
//		line 199
double lAction(int*, double*);

// evaluates the Wilson loop starting from a given point, with given dimensions and in given directions
//		line 224
double lWLoop(int*, double*, int,int,int,int, int,int, int,int);

// evaluates the action contribution given by one point, given a direction
//		line 249
double lPPAction(int*, double*, int, int,int,int,int);

#endif
