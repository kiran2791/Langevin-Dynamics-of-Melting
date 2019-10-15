/*
Name: LD_woutput.cpp
Description: Simulates Polymer crystallization using nchain number of chains (Bead-Spring chain in Langevin Dynamics)
with harmonic bond, harmonic bond angle, multi-harmonic dihedral angle interactions and excluded volume interactions. This file contains all classes and functions and requires only standard libraries.
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "vector1.h"
#include "random1.h"
#include "readinput.h"
#include "initialization.h"
#include "operators.h"
#include "GOP_LOP.h"
//#include "angles.h"
//#include "dihedrals.h"

using namespace std;

// ****************Declare Functions **************************************
void calculateacceleration(vector** position, vector** velocity, vector** acceleration, ofstream &, int chainindex);
void interchainLJ(vector** position,vector** acceleration, ofstream &);

int t;
// ********************** Main Function ***********************************
int main(int argc, char *argv[])
{
  stringstream convert(argv[1]);
  int myint;
  if (!(convert >> myint))
    {
      myint = 0;
    }
  
// ********************* Declare Variables ********************************
  srand(time(NULL)+myint);

  vector** position = new vector*[nchain];
  vector** velocity = new vector*[nchain];
  vector** acceleration = new vector*[nchain];

  double rg, rendtoend[nchain];
  double kineticenergy;
  double *gacc = new double;
  double *sysacc = new double;
  vector momentum, centerofmass, centerofmassvelocity;
  vector** initialposition = new vector*[nchain];
  for (int i=0; i<nchain; i++)
    {
      position[i] = new vector[N];
      velocity[i] = new vector[N];
      acceleration[i] = new vector[N];
      initialposition[i] = new vector[N];
    }
  
  ofstream vmdwrite;	vmdwrite.open("coordinates.xyz");
  //ofstream velwrite; velwrite.open("velocities.dat");
  ofstream energyout; energyout.open("energy.dat");
  ofstream radiiout; radiiout.open("radii.dat");
  ofstream dispout; dispout.open("displacements.dat");
  ofstream opwrite; opwrite.open("orderparameters.dat");

  //Set files output style
  vmdwrite.setf(ios_base::showpoint);
  //velwrite.setf(ios_base::showpoint);
  energyout.setf(ios_base::showpoint);
  radiiout.setf(ios_base::showpoint);
  dispout.setf(ios_base::showpoint);
  opwrite.setf(ios_base::showpoint);
  
  

  // ************************ Initialize Monomers *****************************
  initialization(position, velocity);
  //calculateacceleration(position, velocity, acceleration, energyout);
  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  initialposition[j][i] = vector(position[j][i].x, position[j][i].y, position[j][i].z);
	}
    }
  
  vmdwrite << int(N*nchain) << endl << "Atoms" << endl;
  for(int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  vmdwrite << "C\t" << initialposition[j][i].x << "\t" << initialposition[j][i].y << "\t" << initialposition[j][i].z << endl;
	}
    }

  cout<< "Time\n";
  energyout << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
  centerofmassvelocity = vector(0.0,0.0,0.0);

  for(int j=0; j<nchain; j++)
    {
      for (int i=0;i<N;i++)
	{
	  centerofmassvelocity = add(centerofmassvelocity,velocity[j][i]);
	}
    }
  centerofmassvelocity = divide(centerofmassvelocity,(double)N);
  
  kineticenergy = 0;
  momentum = vector(0,0,0);
  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[j][i],centerofmassvelocity)), 2);
	  momentum = add(momentum, multiply(subtract(velocity[j][i],centerofmassvelocity), mass));
	}
    }
  
  kineticenergy /= (double)N;
  t = 0;
  energyout << t*h << "\t" << (2.0/3.0) * kineticenergy <<"\t"<< kineticenergy << "\t" ;

  // ******************************* MAIN LOOP ********************************
  for(t=0; t<neq; t++)
    {
      if(t%tscr == 0) cout << t+1 << "\n";
      /*
      Velocity Verlet Integration
      */
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      position[j][i] = add(position[j][i], multiply(velocity[j][i], h));
	      position[j][i] = add(position[j][i], multiply(acceleration[j][i], 0.5*h*h));
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, energyout, j);
	}
      interchainLJ(position, acceleration, energyout);

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      // Write Positions to VMD File
      if(t % tvmd == 0)
	{
	  vmdwrite << int(N*nchain) << endl << "Atom type\tx\ty\tz\tvx\tvy\tvz" << endl;
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  vmdwrite << "C\t" << position[j][i].x << "\t" << position[j][i].y << "\t" << position[j][i].z << "\t" << velocity[j][i].x << "\t" <<velocity[j][i].y << "\t" << velocity[j][i].z << endl;
		}
	    }
	}
      
      // Write Energy, Momentum to Energyout File, first by defining center of mass velocity
      if (t % tener == 0)
	{
	  centerofmassvelocity = vector(0.0,0.0,0.0);
	  for (int j=0; j<nchain; j++)
	    {
	      for (int i=0;i<N;i++)
		{
		  centerofmassvelocity = add(centerofmassvelocity,velocity[j][i]);
		}
	    }
	  centerofmassvelocity = divide(centerofmassvelocity,(double)N);
	  kineticenergy = 0;
	  momentum = vector(0,0,0);
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[j][i],centerofmassvelocity)), 2);
		  momentum = add(momentum, multiply(subtract(velocity[j][i],centerofmassvelocity), mass));
		}
	    }
	  kineticenergy /= (double)N;

	  energyout << (t+1)*h << "\t" << (2.0/3.0) * kineticenergy <<"\t"<< kineticenergy << "\t" ;
	}
      
      // Calculate Center of Mass Position
      centerofmass = vector(0,0,0);
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      centerofmass = add(centerofmass, position[j][i]);
	    }
	  centerofmass = divide(centerofmass, (double)N);
	}
      
      // Write Chain Size to radii
      for (int i=0; i<nchain; i++)
	{
	  rg = 0;
	  rendtoend[i] = magnitude( subtract(position[i][N-1], position[i][0]) );
	}

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      rg += pow( magnitude(subtract(position[j][i], centerofmass)), 2 );
	    }
	}
      rg /= (double)(N+1);
      rg = sqrt(rg);
      radiiout << (t+1)*h << "\t" << rg << "\t" << rendtoend << endl;
    }

  // ************************END OF EQUILIBRATION LOOP *************************

  Tstar = Tstar2;
  kT = Tstar2*epsilon;

  // ************************QUENCH LOOP****************************************
  
  for(t=neq; t<(neq+nquench); t++)
    {
      if(t%tscr == 0) cout << t+1 << "\n";
      /*
      Velocity Verlet Integration
      */
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      position[j][i] = add(position[j][i], multiply(velocity[j][i], h));
	      position[j][i] = add(position[j][i], multiply(acceleration[j][i], 0.5*h*h));
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, energyout, j);
	}
      interchainLJ(position, acceleration, energyout);

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      //Calculation of the local and global order parameters
      if(t % top == 0)
      	{
      	  OP(position,gacc,sysacc);
      	  opwrite << setprecision(9) << (t+1)*h << "\t" << setw(6) << *gacc << "\t" << setw(6) << *sysacc << "\n";
      	}
      
      // Write Positions to VMD File
      if(t % tvmd == 0)
	{
	  vmdwrite << int(N*nchain) << endl << "Atom type\tx\ty\tz\tvx\tvy\tvz" << endl;
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  vmdwrite << "C\t" << position[j][i].x << "\t" << position[j][i].y << "\t" << position[j][i].z << "\t" << velocity[j][i].x << "\t" <<velocity[j][i].y << "\t" << velocity[j][i].z << endl;
		}
	    }
	}
      
      // Write Energy, Momentum to Energyout File, first by defining center of mass velocity
      if (t % tener == 0)
	{
	  centerofmassvelocity = vector(0.0,0.0,0.0);
	  for (int j=0; j<nchain; j++)
	    {
	      for (int i=0;i<N;i++)
		{
		  centerofmassvelocity = add(centerofmassvelocity,velocity[j][i]);
		}
	    }
	  centerofmassvelocity = divide(centerofmassvelocity,(double)N);
	  kineticenergy = 0;
	  momentum = vector(0,0,0);
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[j][i],centerofmassvelocity)), 2);
		  momentum = add(momentum, multiply(subtract(velocity[j][i],centerofmassvelocity), mass));
		}
	    }
	  kineticenergy /= (double)N;

	  energyout << (t+1)*h << "\t" << (2.0/3.0) * kineticenergy <<"\t"<< kineticenergy << "\t" ;
	}
      
      // Calculate Center of Mass Position
      centerofmass = vector(0,0,0);
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      centerofmass = add(centerofmass, position[j][i]);
	    }
	  centerofmass = divide(centerofmass, (double)N);
	}
      
      // Write Chain Size to radii
      for (int i=0; i<nchain; i++)
	{
	  rg = 0;
	  rendtoend[i] = magnitude( subtract(position[i][N-1], position[i][0]) );
	}

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      rg += pow( magnitude(subtract(position[j][i], centerofmass)), 2 );
	    }
	}
      rg /= (double)(N+1);
      rg = sqrt(rg);
      radiiout << (t+1)*h << "\t" << rg << "\t" << rendtoend << endl;
    }

  // **************************END OF QUENCH LOOP *********************************

  // **************************HEATING LOOP **************************************

  for(t=(nquench + neq); t<(neq + nquench + nheat); t++)
    {
      Tstar = Tstar + (h*hrate);
      kT = Tstar * epsilon;
      if(t%tscr == 0) cout << t+1 << "\n";
      
      /*
      Velocity Verlet Integration
      */
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      position[j][i] = add(position[j][i], multiply(velocity[j][i], h));
	      position[j][i] = add(position[j][i], multiply(acceleration[j][i], 0.5*h*h));
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, energyout, j);
	}
      interchainLJ(position, acceleration, energyout);

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
	    }
	}

      //Calculation of the local and global order parameters
      if(t % top == 0)
      	{
      	  OP(position,gacc,sysacc);
      	  opwrite << setprecision(9) << (t+1)*h << "\t" << setw(6) << *gacc << "\t" << setw(6) << *sysacc << "\n";
      	}
      
      // Write Positions to VMD File
      if(t % tvmd == 0)
	{
	  vmdwrite << int(N*nchain) << endl << "Atom type\tx\ty\tz\tvx\tvy\tvz" << endl;
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  vmdwrite << "C\t" << position[j][i].x << "\t" << position[j][i].y << "\t" << position[j][i].z << "\t" << velocity[j][i].x << "\t" <<velocity[j][i].y << "\t" << velocity[j][i].z << endl;
		}
	    }
	}
      
      // Write Energy, Momentum to Energyout File, first by defining center of mass velocity
      if (t % tener == 0)
	{
	  centerofmassvelocity = vector(0.0,0.0,0.0);
	  for (int j=0; j<nchain; j++)
	    {
	      for (int i=0;i<N;i++)
		{
		  centerofmassvelocity = add(centerofmassvelocity,velocity[j][i]);
		}
	    }
	  centerofmassvelocity = divide(centerofmassvelocity,(double)N);
	  kineticenergy = 0;
	  momentum = vector(0,0,0);
	  for (int j=0; j<nchain; j++)
	    {
	      for(int i=0; i<N; i++)
		{
		  kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[j][i],centerofmassvelocity)), 2);
		  momentum = add(momentum, multiply(subtract(velocity[j][i],centerofmassvelocity), mass));
		}
	    }
	  kineticenergy /= (double)N;

	  energyout << (t+1)*h << "\t" << (2.0/3.0) * kineticenergy <<"\t"<< kineticenergy << "\t" ;
	}
      
      // Calculate Center of Mass Position
      centerofmass = vector(0,0,0);
      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      centerofmass = add(centerofmass, position[j][i]);
	    }
	  centerofmass = divide(centerofmass, (double)N);
	}
      
      // Write Chain Size to radii
      for (int i=0; i<nchain; i++)
	{
	  rg = 0;
	  rendtoend[i] = magnitude( subtract(position[i][N-1], position[i][0]) );
	}

      for (int j=0; j<nchain; j++)
	{
	  for(int i=0; i<N; i++)
	    {
	      rg += pow( magnitude(subtract(position[j][i], centerofmass)), 2 );
	    }
	}
      rg /= (double)(N+1);
      rg = sqrt(rg);
      radiiout << (t+1)*h << "\t" << rg << "\t" << rendtoend << endl;
    }
  
  
  // ********************Deallocate Variables**********************************
  vmdwrite.close();
  //velwrite.close();
  energyout.close();
  radiiout.close();
  dispout.close();
  opwrite.close();

  for (int i=0;i<nchain;i++)
    {
      delete [] position[i];
      delete [] velocity[i];
      delete [] acceleration[i];
      delete [] initialposition[i];
    }
  
      

  delete [] position;
  delete [] velocity;
  delete [] acceleration;
  delete [] initialposition;
  delete gacc;
  delete sysacc;

  return 0;
}

// *************************Define Functions**************************

void calculateacceleration(vector** position, vector** velocity, vector** acceleration, ofstream &energyout, int chainindex)
{
  double bacc,ljacc,ranacc,friacc,ebond,elj,eangle,edihedral,etotal;
  ebond = 0.0;
  elj = 0.0;
  eangle = 0.0;
  edihedral = 0.0;
  double velocity_sigma = sqrt(kT/mass);
  // calculate random force
  //ranacc = 0.0;
  for(int i=0; i<N; i++)
    {
      double force_sigma = sqrt(2*kT*zeta/h);
      double ax,ay,az;
      ax = force_sigma*nr()/mass;
      ay = force_sigma*nr()/mass;
      az = force_sigma*nr()/mass;
      acceleration[chainindex][i].x = ax;
      acceleration[chainindex][i].y = ay;
      acceleration[chainindex][i].z = az;
      //ranacc += sqrt(pow(ax,2) + pow(ay,2) + pow(az,2));
    }
  //ranacc /= N;

  // calculate drag force
  //friacc = 0.0;
  for(int i=0; i<N; i++)
    {
      acceleration[chainindex][i] = add(acceleration[chainindex][i], multiply(velocity[chainindex][i], -zeta/mass));
      friacc += (-magnitude(velocity[chainindex][i])*zeta/mass);
    }
  //friacc /= N;

  // calculate bond force
  vector dr1, dr2, f, f1, f2;
  //bacc = 0.0;
  for(int i=0; i<N; i++)
    {
      if(i==0)
	{
	  dr1 = subtract(position[chainindex][i], position[chainindex][i+1]);
	  if (magnitude(dr1)<0.5 || magnitude(dr1) > 3.5)
	    {
	      for (int j=0; j<N ; j++)
		{
		  velocity[chainindex][i] = vector( velocity_sigma*nr(), velocity_sigma*nr(), velocity_sigma*nr() );
		}
	    }
	  f = multiply(unit(dr1), -2.0*k*(magnitude(dr1)-b));
	  ebond += (k*pow((magnitude(dr1)-b),2));
	}
      else if(i == N-1)
	{
	  dr2 = subtract(position[chainindex][i], position[chainindex][i-1]);
	  if (magnitude(dr2)<0.5 || magnitude(dr2) > 3.5)
	    {
	      for (int j=0; j<N ; j++)
		{
		  velocity[chainindex][i] = vector( velocity_sigma*nr(), velocity_sigma*nr(), velocity_sigma*nr() );
		}
	    }
	  f = multiply(unit(dr2), -2.0*k*(magnitude(dr2)-b));
	}
      else
	{
	  dr1 = subtract(position[chainindex][i], position[chainindex][i+1]);
	  dr2 = subtract(position[chainindex][i], position[chainindex][i-1]);
	  if (magnitude(dr1)<0.5 || magnitude(dr1) > 3.5 || magnitude(dr1)< 0.5 || magnitude(dr2) > 3.5)
	    {
	      for (int j=0; j<N ; j++)
		{
		  velocity[chainindex][i] = vector( velocity_sigma*nr(), velocity_sigma*nr(), velocity_sigma*nr() );
		}
	    }
	  f1 = multiply(unit(dr1), -2.0*k*(magnitude(dr1)-b));
	  f2 = multiply(unit(dr2), -2.0*k*(magnitude(dr2)-b));
	  f = add(f1, f2);
	  ebond += (k*pow((magnitude(dr1)-b),2));
	}
      acceleration[chainindex][i] = add(acceleration[chainindex][i], divide(f, mass));
      //bacc += (magnitude(f));
    }
  //bacc /= N;


  /*
  Calculate the forces arising out of the Lennard Jones 12-6 potential
  U = epsilon * ((sigma/r)^12 - 2(sigma/r)^6)
  Force is therefore, (-12 epsilon/r) * ((sigma/r)^12 - (sigma/r)^6)
  */
  
  vector dr,F_LJ;
  double drsq,FLJ,r6inv;
  //ljacc = 0.0;
  for (int i=0;i<N-1;i++)
    {
      for(int j=(i+1);j<N;j++)
	{
	  dr = subtract(position[chainindex][i],position[chainindex][j]);
	  drsq = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
	  if (drsq < (rcutoff*rcutoff))
	    {
	      r6inv = 1.0/(drsq*drsq*drsq);
	      if (abs(i-j)<5)
		{
		  FLJ = (12.0*epsilon/drsq)*((pow(sigma1,12)*(r6inv*r6inv)) - (pow(sigma1,6)*r6inv));
		  elj += (epsilon * ((pow(sigma1,12)*(r6inv*r6inv)) - (pow(sigma1,6)*r6inv)));
		}
	      else
		{
		  FLJ = (12.0*epsilon/drsq)*((pow(sigma2,12)*(r6inv*r6inv)) - (pow(sigma2,6)*r6inv));
		  elj += (epsilon * ((pow(sigma2,12)*(r6inv*r6inv)) - (pow(sigma1,6)*r6inv)));
		}
	      acceleration[chainindex][i].x += (FLJ*dr.x/mass);
	      acceleration[chainindex][i].y += (FLJ*dr.y/mass);
	      acceleration[chainindex][i].z += (FLJ*dr.z/mass);
	      //ljacc += (FLJ*magnitude(dr)/mass);

	      acceleration[chainindex][j].x -= (FLJ*dr.x/mass);
	      acceleration[chainindex][j].y -= (FLJ*dr.y/mass);
	      acceleration[chainindex][j].z -= (FLJ*dr.z/mass);
	    }
	}
    }
  //ljacc /= N;
  //cout << ranacc << "\t" << friacc << "\t" << bacc << "\t" << ljacc << "\t" ;

  /*
  Calculate the forces due to the bond angle potential
  */
  vector d21,d32;
  vector v1,v2,v3,v4;
  double cs,fa,d_b2,d_b1,d_b1b2,facc,cd;
  //vector acc0[N-2],acc1[N-2],acc2[N-2];
  // double mainacc;
  // mainacc = 0.0;

  // for (int i=0;i<N-2;i++)
  //   {
  //     acc0[i] = vector(0.0,0.0,0.0);
  //     acc1[i] = vector(0.0,0.0,0.0);
  //     acc2[i] = vector(0.0,0.0,0.0);
  //   }
  
  for (int i=0;i<N-2;i++)
    {
      d21 = subtract(position[chainindex][i+1],position[chainindex][i]);
      d32 = subtract(position[chainindex][i+2],position[chainindex][i+1]);
      d_b2 = dot(d32,d32);
      d_b1 = dot(d21,d21);
      d_b1b2 = dot(d32,d21);
      cs = d_b1b2/sqrt(d_b2*d_b1);
      fa = -2.0*ktheta*(cs-cs0);
      eangle += ktheta * pow((cs-cs0),2);
      cd = sqrt(d_b1*d_b2);
      
      //Force on atom 1 of the triplet
      v1 = multiply(d21,(fa * cs/d_b1));
      v2 = multiply(d32,(-1.0 * fa/cd)); 
      f1 = add(v1,v2);
      acceleration[chainindex][i] = add(acceleration[chainindex][i],f1);

      //Force on atom 3 of the triplet
      v3 = multiply(d32,(-1.0 * cs * fa/d_b2));
      v4 = multiply(d21,(fa/cd));
      f2 = add(v3,v4);
      acceleration[chainindex][i+2] = add(acceleration[chainindex][i+2],f2);

      //Force on atom 2 of the triplet
      acceleration[chainindex][i+1] = add(acceleration[chainindex][i+1],multiply(add(f1,f2),-1.0));

      // acc0[i] = add(acc0[i],f1);
      // acc2[i] = add(acc2[i],f2);
      // acc1[i] = add(acc1[i],multiply(add(f1,f2),-1.0));
    }

  // for (int i=0;i< N-4 ; i++)
  //   {
  //     mainacc += magnitude(add(add(acc2[i],acc1[i+1]),acc0[i+2]));
  //   }
  // mainacc += magnitude(acc0[0]);
  // mainacc += magnitude(add(acc0[1],acc1[0]));
  // mainacc += magnitude(add(acc1[N-3],acc2[N-4]));
  // mainacc += magnitude(acc2[N-3]);
  // mainacc /= N;
  // cout << "\t" << mainacc << "\t";

  /*
  Calculate the forces arising out of the dihedral angle potential
  */
  vector r12,r32,r43,f4,f3,c1232,c3243; 
  double csphi,fd,m_c1232,m_c3243,m_r12,m_r32,m_r43,co1,co2;
  double A1,A2,A3,A4,A5,A6,A7,A8,alpha,beta;
  //  vector dacc0[N-3],dacc1[N-3],dacc2[N-3],dacc3[N-3],accsum;

  // for (int i=0;i<N-3;i++)
  //   {
  //     dacc0[i] = vector(0.0,0.0,0.0);
  //     dacc1[i] = vector(0.0,0.0,0.0);
  //     dacc2[i] = vector(0.0,0.0,0.0);
  //     dacc3[i] = vector(0.0,0.0,0.0);
  //   }
  
  for (int i=0; i<N-3 ; i++)
    {
      //Defining the various coefficients that will be required
      r12 = subtract(position[chainindex][i],position[chainindex][i+1]);
      r32 = subtract(position[chainindex][i+2],position[chainindex][i+1]);
      r43 = subtract(position[chainindex][i+3],position[chainindex][i+2]);
      c1232 = cross(r12,r32);
      c3243 = cross(r32,r43);
      //c2143 = cross(r21,r43);
      m_c1232 = magnitude(c1232);
      m_c3243 = magnitude(c3243);
      m_r12 = magnitude(r12);
      m_r32 = magnitude(r32);
      m_r43 = magnitude(r43);
      csphi = dot(c1232,c3243)/(m_c1232 * m_c3243);
      fd = k1 + (4.0 * k2 * csphi) + (12.0 * k3 * csphi * csphi) - (3.0 * k3);
      edihedral += (k1*(1.0 - csphi) + k2*(1.0 - (2.0*pow(csphi,2) - 1.0)) + k3*(1.0 - (4.0*pow(csphi,3) - 3.0*csphi)));

      co1 = fd * dot(c1232,c3243)/(m_c3243 * pow((pow(m_r12*m_r32,2) - pow(dot(r12,r32),2)),1.5));
      A1 = co1 * dot(r12,r32);
      A2 = co1 * -1.0*pow(m_r32,2);
      A3 = fd * dot(r32,r43)/(m_c1232 * m_c3243);
      A4 = fd * -1.0*pow(m_r32,2)/(m_c1232 * m_c3243);
      co2 = fd * dot(c1232,c3243)/(m_c1232 * pow((pow(m_r32*m_r43,2) - pow(dot(r32,r43),2)),1.5));
      A5 = dot(r32,r43) * co2;
      A6 = -1.0 * pow(m_r32,2) * co2;
      A7 = fd * dot(r12,r32)/(m_c1232 * m_c3243);
      A8 = fd * -1.0*pow(m_r32,2)/(m_c1232 * m_c3243);

      //Force on atom 1
      f1 = add(multiply(r32,A1),multiply(r12,A2));
      f1 = add(f1,multiply(r32,A3));
      f1 = add(f1,multiply(r43,A4));
      f1 = divide(f1,mass);

      //Force on atom 4
      f4 = add(multiply(r32,A5),multiply(r43,A6));
      f4 = add(f4,multiply(r32,A7));
      f4 = add(f4,multiply(r12,A8));

      //Force on atom 2
      alpha = dot(cross(r43,f4),f1)/dot(cross(r32,f4),f1);
      v1 = subtract(r32,r12);
      beta = -1.0 * dot(cross(v1,f1),f4)/(dot(cross(r32,f1),f4));
      f2 = add(multiply(f4,alpha),multiply(f1,beta));

      //Force on atom 3
      f3 = add(multiply(f1,-1.0),multiply(f2,-1.0));
      f3 = add(f3,multiply(f4,-1.0));      

      //Adding the forces onto each of the atoms
      acceleration[chainindex][i] = add(acceleration[chainindex][i],f1);
      acceleration[chainindex][i+1] = add(acceleration[chainindex][i+1],f2);
      acceleration[chainindex][i+2] = add(acceleration[chainindex][i+2],f3);
      acceleration[chainindex][i+3] = add(acceleration[chainindex][i+3],f4);
      
      
      //Adding to the accumulators
      // dacc0[i] = add(dacc0[i],f1);
      // dacc1[i] = add(dacc1[i],f2);
      // dacc2[i] = add(dacc2[i],f3);
      // dacc3[i] = add(dacc3[i],f4);
    }
  
  // mainacc = 0.0;
  
  // for (int i=0;i< N-4 ; i++)
  //   {
  //     accsum = vector(0.0,0.0,0.0);
  //     accsum = add(dacc0[i+3],dacc1[i+2]);
  //     accsum = add(accsum,dacc2[i+1]);
  //     accsum = add(accsum,dacc3[i]);
  //     mainacc += magnitude(accsum);
  //   }
  
  // mainacc += magnitude(dacc0[0]);
  // mainacc += magnitude(add(dacc0[1],dacc1[0]));
  // mainacc += magnitude(add(add(dacc0[2],dacc1[1]),dacc2[0]));

  // mainacc += magnitude(add(add(dacc1[N-4],dacc2[N-5]),dacc3[N-6]));
  // mainacc += magnitude(add(dacc2[N-4],dacc3[N-5]));
  // mainacc += magnitude(dacc3[N-4]);
  // mainacc /= N;
  // cout << "\t" << mainacc << "\n" ;

  etotal = ebond + elj + eangle + edihedral;
  if(t % tener == 0)
    {
      energyout << setw(6)<<ebond << "\t" <<  setw(6)<<elj << "\t" << setw(6)<<eangle << "\t" << setw(6)<<edihedral << "\t" << setw(6)<<etotal << "\t";
    }
}

void interchainLJ(vector** position, vector** acceleration, ofstream &energyout)
{
  /*
  Calculates the Lennard Jones Interactions between chains.
  Therefore, the inter-particle distance between monomers of different chains is calculated here.
  */
  vector dr,F_LJ;
  double drsq,FLJ,r6inv;
  double elj = 0.0;
  double rcutsq;
  rcutsq = rcutoff*rcutoff;
  //ljacc = 0.0;
  for (int k=0; k<nchain-1; k++)
    {
      for (int l=(k+1); l<nchain; l++)
	{
	  for (int i=0;i<N;i++)
	    {
	      for(int j=0;j<N;j++)
		{
		  dr = subtract(position[k][i],position[l][j]);
		  drsq = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
		  if (drsq < (rcutsq))
		    {
		      r6inv = 1.0/(drsq*drsq*drsq);
		      FLJ = (12.0*epsilon/drsq)*((pow(sigma2,12)*(r6inv*r6inv)) - (pow(sigma2,6)*r6inv));
		      elj += (epsilon * ((pow(sigma2,12)*(r6inv*r6inv)) - (pow(sigma2,6)*r6inv)));
		      acceleration[k][i].x += (FLJ*dr.x/mass);
		      acceleration[k][i].y += (FLJ*dr.y/mass);
		      acceleration[k][i].z += (FLJ*dr.z/mass);
		      //ljacc += (FLJ*magnitude(dr)/mass);

		      acceleration[l][j].x -= (FLJ*dr.x/mass);
		      acceleration[l][j].y -= (FLJ*dr.y/mass);
		      acceleration[l][j].z -= (FLJ*dr.z/mass);
		    }
		}
	    }
	}
    }
  if(t % tener == 0)
    {
      energyout << setw(6)<< elj << "\n";
    }
}
