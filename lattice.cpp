#ifndef __LATTICE_CU__
#define __LATTICE_CU__

#include "lattice.hpp"

#include <vector>
#include <cstdlib> 
#include <cmath>
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
using namespace std;

// returns the total dimension of the lattice
//  given the array of dimensions
int dTot(int *d){
  return 4*d[0]*d[1]*d[2]*d[3];
}

// returns the "1D" position on the lattice
//  accepts a 4-tuple of int
int myPos(int *dim, int i, int j, int k, int l){
  return (l*dim[2]*dim[1]*dim[0]+k*dim[1]*dim[0]+j*dim[0]+i)*4;
}

// returns the "1D" position on lattice
//  overloading for vector of int
int myPos(int *dim, int *n){
  return (n[3]*dim[2]*dim[1]*dim[0]+n[2]*dim[1]*dim[0]+n[1]*dim[0]+n[0])*4;
}

// initialisation of the lattice to constant
//  d is the number of points in a direction
//  h must be already allocated
//	val is the constant value of initialisation
void initLattice(int *d, double *h, double val){
  for(int l=0; l<d[3]; ++l){
    for(int k=0; k<d[2]; ++k){
      for(int j=0; j<d[1]; ++j){
        for(int i=0; i<d[0]; ++i){
          for(int dir=0; dir<4; ++dir){   
            h[myPos(d,i,j,k,l)+dir] = val;
          }
        }
      }
    }
  }
};

// initialisation of the lattice to random in [-pi,pi]
// see above def.
void initLattice(int *d, double *h){
  for(int l=0; l<d[3]; ++l){
    for(int k=0; k<d[2]; ++k){
      for(int j=0; j<d[1]; ++j){
        for(int i=0; i<d[0]; ++i){
          for(int dir=0; dir<4; ++dir){
            h[myPos(d,i,j,k,l)+dir] = -M_PI + 2*M_PI*rand()/RAND_MAX;
          }
        }
      }
    }
  }
};

// copy initialisation (given another lattice)
// see above def.
void initLattice(int *d, double *to, double *from){
  int N = dTot(d);
  for(int i=0; i<N; ++i){
    to[i] = from[i];
  }
};


// function to print the configuation on std output
//	prints the links of the given lattice h
//	 with given dimensions stored in d
void printL(int *d, double *h){
  int myP;
  cout << "dim = " << d[0] << "x" << d[1] 
       <<  "x" << d[2] << "x" << d[3] << endl;
  for(int l=0; l<d[3]; ++l){
    for(int k=0; k<d[2]; ++k){
      for(int j=0; j<d[1]; ++j){
        for(int i=0; i<d[0]; ++i){
          myP = myPos(d,i,j,k,l);
          cout << "(" << h[myP+0] << "," << h[myP+1] << "," << h[myP+2] << "," << h[myP+3] << ")  ";
        }
        cout << endl << endl;
      }
      cout << endl << endl << endl;
    }
    cout << "---" << endl;
  }
};

// function to print the configuation on a file
//	prints the links of the given lattice h
//	 with given dimensions stored in d
//	in the file fName
// - analogous to preceding function
void printL(int *d, double *h, string fName){
  ofstream w;
  int myP;
  w.open(fName.c_str());
  w << "dim = " << d[0] << "x" << d[1] 
    <<  "x" << d[2] << "x" << d[3] << endl;
  for(int l=0; l<d[3]; ++l){
    for(int k=0; k<d[2]; ++k){
      for(int j=0; j<d[1]; ++j){
        for(int i=0; i<d[0]; ++i){
          myP = myPos(d,i,j,k,l);
          w << "(" << h[myP+0] << "," << h[myP+1] << "," << h[myP+2] << "," << h[myP+3] << ")  ";
        }
        w << endl;
      }
      w << endl << endl;
    }
    w << "---" << endl;
  }
  w.close();
};

// prints an observable
//	d is the number of measurements
//	ob stores the measured values
//	fName is the destination file
void printO(const int d, double* ob,  string fName){
  ofstream w;
  w.open(fName.c_str());
  for(int i=0; i<d; ++i)
    w << ob[i] << endl;
  w.close();
}

// save current configuration in the given position (already initialized)
//  conf is the configuration
//  sH is the saving spot
void saveConf(int *d, double *conf, double *sH){
  initLattice(d, sH, conf);  
};

// loads the configuration from a saved one
//  conf is the configuration
//  sH is the saving spot
void loadConf(int *d, double *conf, double *sH){
  initLattice(d, conf, sH);  
};

// evaluates the plaquette action in the wanted directions 
//  starting from the wanted point
//  going in the positive verse
double lPAction(int *d, double *lat, int dir1, int dir2, int i, int j, int k, int l){
  int n[4] = {i,j,k,l};  
  int t1 = n[dir1];
  double phi;
  if(dir2!=dir1){
    phi = lat[myPos(d,n)+dir1]-lat[myPos(d,n)+dir2];
    n[dir1] = (n[dir1]+1)%d[dir1];
    phi += lat[myPos(d,n)+dir2];
    n[dir1] = t1;
    n[dir2] = (n[dir2]+1)%d[dir2];
    phi -= lat[myPos(d,n)+dir1];
    return cos(phi);
  }
  else
    return 0;
}

// evaluates the plaquette action in the wanted directions 
//  overloading for int-4-array of position
//  see preceding function
double lPAction(int *d, double *lat, int dir1, int dir2, int *n){
  double phi;
  int t1 = n[dir1];
  int t2 = n[dir2];
  if(dir2>dir1)
    return lPAction(d, lat, dir2, dir1, n);
  if(dir2!=dir1){
    phi = lat[myPos(d,n)+dir1]-lat[myPos(d,n)+dir2];
    n[dir1] = (n[dir1]+1)%d[dir1];
    phi += lat[myPos(d,n)+dir2];
    n[dir1] = t1;
    n[dir2] = (n[dir2]+1)%d[dir2];
    phi -= lat[myPos(d,n)+dir1];
    n[dir2] = t2;
    return cos(phi);   
  }
  else
    return 0;
}

// evaluates the lattice action on the whole lattice
//	input: dimensions of the lattice - d
//				 lattice - lat
//  sums all the contributions from the single plaquettes in the lattice
double lAction(int *d, double *lat){
  double S = 0.0;
    // looping on all points
  for(int l=0; l<d[3]; ++l){
    for(int k=0; k<d[2]; ++k){
      for(int j=0; j<d[1]; ++j){
        for(int i=0; i<d[0]; ++i){
            // looping on all directions
            for(int d1=0; d1<4; ++d1){
              for(int d2=d1+1; d2<4; ++d2){
                S -= lPAction(d,lat,d1,d2,i,j,k,l);
              }
            }
          }
        }
      }
    }
  return S;
}

// evaluates the Wilson Loop 
// starting from point i,j,k,l
// having dimensions n x m
// moving in directions d1 x d2
// only in the bulk
double lWLoop(int *d, double *lat, int i, int j, int k, int l, int m, int n, int d1, int d2){
  int N[4] = {i,j,k,l};
  double phi = 0;
  for(int I=0; I<m; ++I){
    phi += lat[myPos(d,N)+d1];
    N[d1]++;    
  }
  for(int I=0; I<n; ++I){
    phi += lat[myPos(d,N)+d2];
    N[d2]++;    
  }  
  for(int I=0; I<m; ++I){
    N[d1]--; 
    phi -= lat[myPos(d,N)+d1];    
  }  
  for(int I=0; I<n; ++I){
    N[d2]--; 
    phi -= lat[myPos(d,N)+d2];    
  } 
  return cos(phi);
}

// evaluates the action contribution given by one point,
//	considering one direction
//	used to evaluate the change in the action once a link is changed
double lPPAction(int* d, double *lat, int dir1, int i, int j, int k, int l){
  int n[4] = {i,j,k,l};  
  double S = 0.0;
  for(int d2=0; d2<4; ++d2){
    // starting from the point
    if(dir1!=d2){
      S -= lPAction(d, lat, dir1, d2, n);
      // including the point, but starting somewhere else
      if(n[d2]!=0){
        n[d2]--;
        S -= lPAction(d, lat, dir1, d2, n);
        n[d2]++;
      }
      else{
        n[d2] = d[d2]-1;
        S -= lPAction(d, lat, dir1, d2, n);
        n[d2] = 0;
      }
    }
  }
  return S;
}

#endif
