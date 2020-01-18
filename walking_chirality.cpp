//Sely85
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include "chirality_function.hpp"
 
using namespace std;

//Define a class with specie, x, y, z
class coordinates
{
public:
  coordinates(std::string sp, double xxx, double yyy, double zzz): specie(sp),x(xxx), y(yyy), z(zzz)
  {}
  coordinates(): specie(), x(0.0), y(0.0), z(0.0)
  {}
  std::string specie;
  double x;
  double y;
  double z;
};


int main(int argc, char *argv[])
{
  unsigned const stacksize = 500;

  std::cout << " " << std::endl;
  std::cout << "Welcome to the WalkingChirality code!"  << std::endl;
  std::cout << "   Usage: " << std::endl;
  std::cout << "          ./WalkingChirality <xyz> <cutoff> <Na>" << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;

  ifstream infile (argv[1]);

  if ( !infile.is_open() )
    {
      std::cout << "ERROR: Could not open file" << std::endl;
      return 0;
    }

  if ( argc < 4) 
    {
      std::cout << "ERROR: Missing input parameter" << std::endl;
      return 0;
    }


  string null;
  int totatom;

  infile >> totatom;
  infile >> null;

  std::cout << " Atoms.........." << totatom << std::endl;

  std::vector<coordinates> coord(totatom);
  for (int i=0; i<totatom; i++)
    {
      infile >> coord[i].specie >> coord[i].x >> coord[i].y >> coord[i].z;
    }    

      double cutoff = atof(argv[2]); //cutoff
      std::cout << " Cutoff........." << cutoff << std::endl; 
      int Na = atof(argv[3]); //cutoff
      std::cout << " Na............." << Na << std::endl; 
      std::cout << " " << std::endl;

      std::cout << "You set Na (sequence of connected length) to " << Na << " and your cutoff is " << cutoff << std::endl;
      std::cout << " " << std::endl;

      //Hence i,j,k,l are taken x-2, x-1, x, x+1, x+2. There are 5 residues accounting for the computation of Ga for the x^th residue
      int cc=1;
      //      int step = (Na-1)/2;
      double res_chi_ss[totatom-Na];
      double totchi = 0.0;

      ofstream fileout("chirality_index.txt", ios::out);
      fileout << "#[1]Number [2]chirality_index " << std::endl;

      ofstream dumpout("chirality_index.dump", ios::out);
      dumpout << "ITEM: TIMESTEP " << std::endl;
      dumpout << "0 " << std::endl;
      dumpout << "ITEM: NUMBER OF ATOMS " << std::endl;
      dumpout << totatom-Na << std::endl;
      dumpout << "ITEM: BOX BOUNDS " << std::endl;
      dumpout << "-100 100" << std::endl;
      dumpout << "-100 100" << std::endl;
      dumpout << "-100 100" << std::endl;
      dumpout << "ITEM: ATOMS id type x y z ga" << std::endl;

      for (int i=0; i<totatom-Na; i++)
	{
	  //	  std::cout << " i " << i  << " last " << totatom-Na << std::endl;
	  double res_chi = 0.0;

	  for (int a=i; a<i+Na-3; a++)
	    {
	      for (int b=a+1; b<i+Na-2; b++)
		{
		  for (int c=b+1; c<i+Na-1; c++)
		    {
		      for (int d=c+1; d<i+Na; d++)
			{
			  double r1[3], r2[3], r3[3], r4[3];

			  r1[0] = coord[a].x; 
			  r1[1] = coord[a].y; 
			  r1[2] = coord[a].z; 

			  r2[0] = coord[b].x; 
			  r2[1] = coord[b].y; 
			  r2[2] = coord[b].z; 

			  r3[0] = coord[c].x; 
			  r3[1] = coord[c].y; 
			  r3[2] = coord[c].z; 

			  r4[0] = coord[d].x; 
			  r4[1] = coord[d].y; 
			  r4[2] = coord[d].z; 

			  res_chi = res_chi + allperm_chirality(r1, r2, r3, r4, Na, cutoff);
		      
			  cc++;
			}
		    }
		}
	    }		      

	  totchi = totchi + res_chi;
	  fileout << i+1 << " " <<  res_chi << std::endl; 
	  dumpout << i+1 << " 1 " << coord[i].x  << " "<< coord[i].y  << " " << coord[i].z  << " " << res_chi << std::endl; 
	}

      totchi = totchi/(totatom-Na);

      std::cout << " Average chirality:       " << totchi << std::endl;
      std::cout << " " << std::endl;
} 


