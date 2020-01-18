//Sely85
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <vector>

double* distance(double crd1[], double crd2[])
{

  double* distance = new double[3];
  for (int i=0; i<3; i++)
    {
      distance[i]=0.0;
    }

  distance[0] = crd2[0] - crd1[0];
  distance[1] = crd2[1] - crd1[1];
  distance[2] = crd2[2] - crd1[2];
  //std::cout << " distance " << distance[0] << " " << distance[1] << " " << distance[2] << std::endl;

  return distance;
}

double norm(double crd[])
{

  double norm=0.0;

  norm = sqrt(crd[0]*crd[0] + crd[1]*crd[1] + crd[2]*crd[2]);

  return norm;
}

double* crossprod(double crd1[], double crd2[])
{
  double* cross = new double[3];
  for (int i=0; i<3; i++)
    {
      cross[i]=0.0;
    }

  cross[0] = crd1[1]*crd2[2]-crd1[2]*crd2[1];
  cross[1] = crd1[2]*crd2[0]-crd1[0]*crd2[2];
  cross[2] = crd1[0]*crd2[1]-crd1[1]*crd2[0];

  return cross;
}

double dotprod(double crd1[], double crd2[])
{
  double dot;

  dot = crd1[0]*crd2[0] +  crd1[1]*crd2[1] + crd1[2]*crd2[2];

  return dot;
}

double single_chirality(double crd1[], double crd2[], double crd3[], double crd4[], int Na, double cutoff)
{
  double*  rij = distance(crd1, crd2);
  double*  rjk = distance(crd2, crd3);
  double*  rkl = distance(crd3, crd4);
  double*  ril = distance(crd1, crd4);

  double nij = norm(rij);
  double njk = norm(rjk);
  double nkl = norm(rkl);
  double nil = norm(ril);

  double singlechi = 0.0;


  if (nij < cutoff && njk < cutoff && nkl < cutoff && nil < cutoff)
    {
      double* cross = crossprod(rij, rkl);
      double triple = dotprod(cross, ril);
      double rijjk = dotprod(rij, rjk);
      double rjkkl = dotprod(rjk, rkl);
      double num = triple*rijjk*rjkkl;

      double den = (nij*njk*nkl)*(nij*njk*nkl)*nil;

      if (abs(num) < 0.0000000000000000001 || abs(den) < 0.0000000000000000001)
	singlechi = 0.0;
      else
	singlechi = num/den;
    }
  
  return singlechi;
}

double allperm_chirality(double crd1[], double crd2[], double crd3[], double crd4[], int Na, double cutoff)
{

  double ri[3], rj[3], rk[3], rl[3];
  double allchi = 0.0;

  //All 24 possible permutations
  int pcnt=0;
  for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
	{
	  if (i != j)
	    {
	      for (int k=0; k<4; k++)
		{
		  if (i != j && j != k && i != k)
		    {
		      for (int l=0; l<4; l++)
			{
			  if (i != j && i != k && i != l && j != k && j != l && k != l)
			    {
			      if (i == 0)
				{
				  ri[0] = crd1[0];
				  ri[1] = crd1[1];
				  ri[2] = crd1[2];
				}
			      else if (i == 1)
				{
				  rj[0] = crd1[0];
				  rj[1] = crd1[1];
				  rj[2] = crd1[2];
				}
			      else if (i == 2)
				{
				  rk[0] = crd1[0];
				  rk[1] = crd1[1];
				  rk[2] = crd1[2];
				}
			      else if (i == 3)
				{
				  rl[0] = crd1[0];
				  rl[1] = crd1[1];
				  rl[2] = crd1[2];
				}

			      if (j == 0)
				{
				  ri[0] = crd2[0];
				  ri[1] = crd2[1];
				  ri[2] = crd2[2];
				}
			      else if (j == 1)
				{
				  rj[0] = crd2[0];
				  rj[1] = crd2[1];
				  rj[2] = crd2[2];
				}
			      else if (j == 2)
				{
				  rk[0] = crd2[0];
				  rk[1] = crd2[1];
				  rk[2] = crd2[2];
				}
			      else if (j == 3)
				{
				  rl[0] = crd2[0];
				  rl[1] = crd2[1];
				  rl[2] = crd2[2];
				}

			      if (k == 0)
				{
				  ri[0] = crd3[0];
				  ri[1] = crd3[1];
				  ri[2] = crd3[2];
				}
			      else if (k == 1)
				{
				  rj[0] = crd3[0];
				  rj[1] = crd3[1];
				  rj[2] = crd3[2];
				}
			      else if (k == 2)
				{
				  rk[0] = crd3[0];
				  rk[1] = crd3[1];
				  rk[2] = crd3[2];
				}
			      else if (k == 3)
				{
				  rl[0] = crd3[0];
				  rl[1] = crd3[1];
				  rl[2] = crd3[2];
				}


			      if (l == 0)
				{
				  ri[0] = crd4[0];
				  ri[1] = crd4[1];
				  ri[2] = crd4[2];
				}
			      else if (l == 1)
				{
				  rj[0] = crd4[0];
				  rj[1] = crd4[1];
				  rj[2] = crd4[2];
				}
			      else if (l == 2)
				{
				  rk[0] = crd4[0];
				  rk[1] = crd4[1];
				  rk[2] = crd4[2];
				}
			      else if (l == 3)
				{
				  rl[0] = crd4[0];
				  rl[1] = crd4[1];
				  rl[2] = crd4[2];
				}
			      
			      allchi = allchi + single_chirality(ri, rj, rk, rl, Na, cutoff);
			      pcnt++;

			    }
			}
		    }
		}
	    }
	}
    }


  return allchi;

}



double fiveres_chirality(double r1[], double r2[], double r3[], double r4[], double r5[], int Na, double cutoff)
{

  double allchi5res = 0.0;

  //we choose 4 elements from a set of 5 --> combination without repetitions
  // n!/(r!(n-r)!) --> 5!/(4!1!) = 5

  //COMB 1: 1 2 3 4
  allchi5res = allchi5res + allperm_chirality(r1, r2, r3, r4, Na, cutoff);

  //COMB 2: 1 2 3 5
  allchi5res = allchi5res + allperm_chirality(r1, r2, r3, r5, Na, cutoff);

  //COMB 3: 1 2 4 5
  allchi5res = allchi5res + allperm_chirality(r1, r2, r4, r5, Na, cutoff);

  //COMB 4: 1 3 4 5
  allchi5res = allchi5res + allperm_chirality(r1, r3, r4, r5, Na, cutoff);

  //COMB 5: 2 3 4 5
  allchi5res = allchi5res + allperm_chirality(r2, r3, r4, r5, Na, cutoff);

  return allchi5res;
}
