#include "lattice.hpp"

#include <ctime>
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <sstream>
using namespace std;
 
int main(int argc, char* argv[]){

	// initialisation of the RNG
  srand(time(NULL)); 
  
	// dimensions of the lattice and number of plaquettes
	int dim[4] = {atof(argv[1]),atof(argv[1]),atof(argv[1]),atof(argv[1])};
	int latsize = dTot(dim);
	int nP = 6*pow(dim[0],4);
	
  // declaration and initialisation of the lattice
  double lat[latsize];
  
  if((string)argv[3]=="R")
  	initLattice(dim, lat); 
  else
  	initLattice(dim, lat, atoi(argv[3])); 
	
	// temp variables to check the acc/rej condition
	double lOld, SOld, dS;
	
	// number of sweeps and constants
	int Nsw = 10000;
  double delta = 0.2;  
	double b = atof(argv[2]);
	
	// observables
  double S[Nsw];
  double E[Nsw];
	
	// file names
	ostringstream strs1;
	ostringstream strs2; 
  strs1 << "stest_" << argv[3] << "_" << argv[1] << "_" << b*1000<<".dat";
  strs2 << "etest_" << argv[3] << "_" << argv[1] << "_" << b*1000<<".dat";
  string str1 = strs1.str();
  string str2 = strs2.str();  

	// main Metropolis loop 
  for(int sw=0; sw<Nsw; ++sw){
  	
  	// measuring the observables
  
    S[sw] = lAction(dim,lat);
    E[sw] = S[sw]/nP;
    
    // loop on the directions
    for(int dir=0; dir<4; ++dir){ 
    
    	// loop(s) on all points    
      for(int l=0; l<dim[3]; ++l){
		  for(int k=0; k<dim[2]; ++k){
    	for(int j=0; j<dim[1]; ++j){
      for(int i=0; i<dim[0]; ++i){
            
              // save old number              
              lOld = lat[myPos(dim,i,j,k,l)+dir];
              // old action value
              SOld = lPPAction(dim, lat, dir, i,j,k,l);
              
              // update              
              lat[myPos(dim,i,j,k,l)+dir] += -delta*M_PI + 2*delta*M_PI*rand()/RAND_MAX;
              
              // difference in the action              
              dS = lPPAction(dim, lat, dir, i,j,k,l) - SOld;
              
              //check condition for rejection -> rejection                
              if(dS>0 && rand()*1.0/RAND_MAX>exp(-b*dS))
                lat[myPos(dim,i,j,k,l)+dir] = lOld;  
    
                
      }
      }
      }
      } // end loop(s) on points
      
    } // end loop on directions
    
  } // end main sweep loop

	// print observables
  printO(Nsw, S, str1);  
  printO(Nsw, E, str2);
  	
	return 0;
	
}
