/*
Name: LangevinDynamics_Chain.cpp
Description: Simulates Rouse Dynamics (Bead-Spring chain in Langevin Dynamics
with harmonic bond interactions and excluded volume interactions. This
file contains all classes and functions and requires only standard libraries.
*/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <mpi.h>
#include "vector1.h"
#include "random1.h"
#include "readinput.h"
#include "initialization.h"
#include "operators.h"
#include "GOP_LOP.h"
#include "mpiutilities.h"
//#include "angles.h"
//#include "dihedrals.h"

using namespace std;

// ****************Declare Functions **************************************
void calculateacceleration(vector** position, vector** velocity, vector** acceleration, int chainindex, double &tote, double kT, double &ULJ, double &Ubond, double &Uangle, double &Udihedral);
void writevmdoutput(vector** position, vector** velocity, ofstream &);
void writeenenergyoutput(vector** position, vector** velocity, ofstream &, vector momentum, vector centerofmass, vector centerofmassvelocity, double kineticenergy, double ebond, double elj, double eangle, double dihedral, double etotal);
void writeradiioutput(vector** position, vector** velocity, ofstream &, vector centerofmass, double rg, double rendtoend);
void interchainLJ(vector** position, vector** acceleration, double &tote);

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

  int initmpi,numberofprocesses,processid;
  int node1, node2, node3, node4;
  
  vector** position = new vector*[nchain];
  vector** velocity = new vector*[nchain];
  vector** acceleration = new vector*[nchain];
  //vector *trueposition = new vector[N];

  double rg, rendtoend;
  double kineticenergy, totalenergy;
  double tosendx,tosendy,tosendz,vtosendx,vtosendy,vtosendz;
  double *gacc = new double;
  double *sysacc = new double;
  vector momentum, centerofmass, centerofmassvelocity;
  stringstream filename;
  totalenergy = 0.0;
  ofstream vmdwrite1; vmdwrite1.open("coordinates_1.xyz");
  ofstream vmdwrite2; vmdwrite2.open("coordinates_2.xyz");
  ofstream vmdwrite3; vmdwrite3.open("coordinates_3.xyz");
  ofstream vmdwrite4; vmdwrite4.open("coordinates_4.xyz");
  ofstream vmdwrite5; vmdwrite5.open("coordinates_5.xyz");
  ofstream vmdwrite0; vmdwrite0.open("coordinates_0.xyz");
  ofstream energyout1; energyout1.open("energy_1.dat");
  ofstream energyout2; energyout2.open("energy_2.dat");
  ofstream energyout3; energyout3.open("energy_3.dat");
  ofstream energyout4; energyout4.open("energy_4.dat");
  ofstream energyout5; energyout5.open("energy_5.dat");
  ofstream energyout0; energyout0.open("energy_0.dat");
  ofstream radiiout1; radiiout1.open("radii_1.dat");
  ofstream radiiout2; radiiout2.open("radii_2.dat");
  ofstream radiiout3; radiiout3.open("radii_3.dat");
  ofstream radiiout4; radiiout4.open("radii_4.dat");
  ofstream radiiout5; radiiout5.open("radii_4.dat");
  ofstream opwrite1; opwrite1.open("orderparameters_1.dat");
  ofstream opwrite2; opwrite2.open("orderparameters_2.dat");
  ofstream opwrite3; opwrite3.open("orderparameters_3.dat");
  ofstream opwrite4; opwrite4.open("orderparameters_4.dat");
  ofstream opwrite5; opwrite5.open("orderparameters_5.dat");
  ifstream readcrystals; readcrystals.open("crystals.xyz"); // File stream for reading a file consisting purely of crystal structures
  string junk;
  double randomnumber1;
  vector** positiontoswap = new vector*[nchain];
  vector** velocitytoswap = new vector*[nchain];
  double totalenergy2, deltae12, deltabeta12, temptotalenergy;
  char a;
  int valueofprocesszero, valueofneighborofzero;
  int *processids = new int[nreplicas+1];
  double acceptanceprobability;
  double ULJ, Ubond, Uangle, Udihedral;
  //cout << kTreplica[3] << "\n";
  for (int j=0; j<(nreplicas+1); j++)
    {
      processids[j] = j;
      //cout << processids[j] << "\t";
    }
  for (int i=0; i<nchain; i++)
    {
      position[i] = new vector[N];
      velocity[i] = new vector[N];
      acceleration[i] = new vector[N];
      positiontoswap[i] = new vector[N];
      velocitytoswap[i] = new vector[N];
    }
  // cout << "\n";
  // random_shuffle(&processids[0], &processids[nreplicas+1]);
  // for (int j=0; j<(nreplicas+1); j++)
  //   {
  //     cout << processids[j] << "\t";
  //   }
  // cout << "\n";

  //Set files output style
  vmdwrite1.setf(ios_base::showpoint);
  vmdwrite2.setf(ios_base::showpoint);
  vmdwrite3.setf(ios_base::showpoint);
  vmdwrite4.setf(ios_base::showpoint);
  vmdwrite5.setf(ios_base::showpoint);
  vmdwrite0.setf(ios_base::showpoint);
  //velwrite.setf(ios_base::showpoint);
  energyout1.setf(ios_base::showpoint);
  energyout2.setf(ios_base::showpoint);
  energyout3.setf(ios_base::showpoint);
  energyout4.setf(ios_base::showpoint);
  energyout5.setf(ios_base::showpoint);
  energyout0.setf(ios_base::showpoint);
  radiiout1.setf(ios_base::showpoint);
  radiiout2.setf(ios_base::showpoint);
  radiiout3.setf(ios_base::showpoint);
  radiiout4.setf(ios_base::showpoint);
  radiiout5.setf(ios_base::showpoint);
  opwrite1.setf(ios_base::showpoint);
  opwrite2.setf(ios_base::showpoint);
  opwrite3.setf(ios_base::showpoint);
  opwrite4.setf(ios_base::showpoint);
  opwrite5.setf(ios_base::showpoint);
  
  //************************INITIALIZE MPI*************************************
  initmpi = MPI_Init(&argc, &argv);

  //***************************************************************************
  
  //*************************FIND NUMBER OF PROCESSES**************************
  initmpi = MPI_Comm_size(MPI_COMM_WORLD,&numberofprocesses);

  //***************************************************************************
  
  /* *************************************************************
     Initialize monomers from each xyz file, using process 0. Then
     send the position an velocity vectors to each slave node.
     *************************************************************
   */
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 0)
    {
      for (int j=0;j<nreplicas;j++)
	{
	  filename.str(string());
	  filename << "last_" << j+1 << ".xyz";
	  //cout << filename.str() << "\n";
	  initialization(position, velocity, filename);
	  for (int i=0; i<nchain; i++)
	    {
	      for (int k = 0; k<N; k++)
		{
		  //tosendx = position[k].x;
		  //tosendy = position[k].y;
		  //tosendz = position[k].z;
		  //MPI_Send(&tosendx,1,MPI_DOUBLE,j+1,1,MPI_COMM_WORLD);
		  //MPI_Send(&tosendy,1,MPI_DOUBLE,j+1,2,MPI_COMM_WORLD);
		  //MPI_Send(&tosendx,1,MPI_DOUBLE,j+1,3,MPI_COMM_WORLD);
		  sendvector(position,tosendx,tosendy,tosendz,j+1,i,k);
		  sendvector(velocity,vtosendx,vtosendy,vtosendz,j+1,i,k);
		}
	    }
	}
    }
  if (processid == 1)
    {
      for (int i=0; i<nchain; i++)
	{
	  
	  for (int k = 0; k<N; k++)
	    {
	      //MPI_Recv(&tosendx,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      //MPI_Recv(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      //MPI_Recv(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      //position[k].x = tosendx;
	      //position[k].y = tosendy;
	      //position[k].z = tosendz;
	      receivevector(position, tosendx, tosendy, tosendz, 0, i,k);
	  
	      //MPI_Recv(&vtosendx,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      //MPI_Recv(&vtosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      //MPI_Recv(&vtosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      //velocity[k].x = vtosendx;
	      //velocity[k].y = vtosendy;
	      //velocity[k].z = vtosendz;
	      receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, i, k);
	    }
	  //cout << "Received " << velocity[k].x << "\t" << velocity[k].y << "\t" << velocity[k].z << "in process " << processid << "\n";
	}
    }

  if (processid == 2)
    {
      for (int i=0; i<nchain; i++)
	{
	  for (int k = 0; k<N; k++)
	    {
	      receivevector(position, tosendx, tosendy, tosendz, 0, i, k);
	      receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, i,k);
	      //cout << "Received " << position[k].x << "\t" << position[k].y << "\t" << position[k].z << "in process " << processid << "\n";
	      //cout << "Received " << velocity[k].x << "\t" << velocity[k].y << "\t" << velocity[k].z << "in process " << processid << "\n";
	    }
	}
    }

  if (processid == 3)
    {
      for (int i=0; i<nchain; i++)
	{
	  for (int k = 0; k<N; k++)
	    {
	      receivevector(position, tosendx, tosendy, tosendz, 0, i, k);
	      receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, i,k);
	      //cout << "Received " << position[k].x << "\t" << position[k].y << "\t" << position[k].z << "in process " << processid << "\n";
	      //cout << "Received " << velocity[k].x << "\t" << velocity[k].y << "\t" << velocity[k].z << "in process " << processid << "\n";
	    }
	}
    }

  if (processid == 4)
    {
      for (int i=0; i<nchain; i++)
	{
	  for (int k = 0; k<N; k++)
	    {
	      receivevector(position, tosendx, tosendy, tosendz, 0, i, k);
	      receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, i,k);
	      //cout << "Received " << position[k].x << "\t" << position[k].y << "\t" << position[k].z << "in process " << processid << "\n";
	      //cout << "Received " << velocity[k].x << "\t" << velocity[k].y << "\t" << velocity[k].z << "in process " << processid << "\n";
	    }
	}
    }
  
  if (processid == 5)
    {
      for (int i=0; i<nchain; i++)
	{
	  for (int k = 0; k<N; k++)
	    {
	      receivevector(position, tosendx, tosendy, tosendz, 0, i, k);
	      receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, i,k);
	      //cout << "Received " << position[k].x << "\t" << position[k].y << "\t" << position[k].z << "in process " << processid << "\n";
	      //cout << "Received " << velocity[k].x << "\t" << velocity[k].y << "\t" << velocity[k].z << "in process " << processid << "\n";
	    }
	}
    }

  // *****************CALCULATES ACCELERATION FOR THE FIRST TIME*****************
  totalenergy = 0.0;
  ULJ = 0.0; Ubond = 0.0; Uangle = 0.0; Udihedral = 0.0;
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 1)
    {
      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	  totalenergy += temptotalenergy;
	}
      interchainLJ(position, acceleration, temptotalenergy);
      totalenergy += temptotalenergy;
    }
  if (processid == 2)
    {
      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	  totalenergy += temptotalenergy;
	}
      interchainLJ(position, acceleration, temptotalenergy);
      totalenergy += temptotalenergy;
    }
  if (processid == 3)
    {
      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	  totalenergy += temptotalenergy;
	}
      interchainLJ(position, acceleration, temptotalenergy);
      totalenergy += temptotalenergy;
    }
  if (processid == 4)
    {
      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	  totalenergy += temptotalenergy;
	}
      interchainLJ(position, acceleration, temptotalenergy);
      totalenergy += temptotalenergy;
    }
  if (processid == 5)
    {
      for (int j=0; j<nchain; j++)
	{
	  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	  totalenergy += temptotalenergy;
	}
      interchainLJ(position, acceleration, temptotalenergy);
      totalenergy += temptotalenergy;
    }

  // *****************WRITE FIRST POSITION AND VELOCITY TO VMD**********************
  // From 1st slave node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 1) writevmdoutput(position, velocity, vmdwrite1);
  // From 2nd slave node
  if (processid == 2) writevmdoutput(position, velocity, vmdwrite2);
  //From 3rd slave node
  if (processid == 3) writevmdoutput(position, velocity, vmdwrite3);
  // From 4th slave node
  if (processid == 4) writevmdoutput(position, velocity, vmdwrite4);
  if (processid == 5) writevmdoutput(position, velocity, vmdwrite5);

  // ***********************************************************************
  // *********************CALCULATING KINETIC ENERGY AND MOMENTUM **************
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 1 || processid == 2 || processid == 3 || processid == 4 || processid == 5)
    {
      centerofmassvelocity = vector(0.0,0.0,0.0);
      for (int i=0; i<nchain; i++)
	{
	  for (int k=0; k<N; k++)
	    {
	      centerofmassvelocity = add(centerofmassvelocity,velocity[i][k]);
	    }
	}
      centerofmassvelocity = divide(centerofmassvelocity,(double)(N*nchain));

      kineticenergy = 0.0;
      momentum = vector(0.0,0.0,0.0);
      for (int i=0; i<nchain; i++)
	{
	  for (int k=0; k<N; k++)
	    {
	      kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[i][k],centerofmassvelocity)),2);
	      momentum = add(momentum, multiply(subtract(velocity[i][k],centerofmassvelocity), mass));
	    }
	}
      kineticenergy /= (double)(N*nchain);	  
    }

  t = 0;

  // ***********SEND KINETIC ENERGY AND MOMENTUM TO HEAD NODE FOR WRITING ENERGY OUTPUT************
  // From 1st node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 1)
    {
      tosendx = momentum.x;
      tosendy = momentum.y;
      tosendz = momentum.z;
      MPI_Send(&kineticenergy,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(&tosendx,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      MPI_Send(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    }
  if (processid == 0)
    {
      MPI_Recv(&kineticenergy,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendx,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendy,1,MPI_DOUBLE,1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendz,1,MPI_DOUBLE,1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      momentum = vector(tosendx,tosendy,tosendz);
      energyout1 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
      energyout1 << t*h << "\t" << (2.0*kineticenergy)/3.0 << "\t" << kineticenergy << "\n";
    }

  // From 2nd node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 2)
    {
      tosendx = momentum.x;
      tosendy = momentum.y;
      tosendz = momentum.z;
      MPI_Send(&kineticenergy,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(&tosendx,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      MPI_Send(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    }
  if (processid == 0)
    {
      MPI_Recv(&kineticenergy,1,MPI_DOUBLE,2,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendx,1,MPI_DOUBLE,2,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendy,1,MPI_DOUBLE,2,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendz,1,MPI_DOUBLE,2,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      momentum = vector(tosendx,tosendy,tosendz);
      energyout2 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
      energyout2 << t*h << "\t" << (2.0*kineticenergy)/3.0 << "\t" << kineticenergy << "\n";
    }

  // From 3rd node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 3)
    {
      tosendx = momentum.x;
      tosendy = momentum.y;
      tosendz = momentum.z;
      MPI_Send(&kineticenergy,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(&tosendx,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      MPI_Send(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    }
  if (processid == 0)
    {
      MPI_Recv(&kineticenergy,1,MPI_DOUBLE,3,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendx,1,MPI_DOUBLE,3,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendy,1,MPI_DOUBLE,3,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendz,1,MPI_DOUBLE,3,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      momentum = vector(tosendx,tosendy,tosendz);
      energyout3 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
      energyout3 << t*h << "\t" << (2.0*kineticenergy)/3.0 << "\t" << kineticenergy << "\n";
    }

  // From 4th node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 4)
    {
      tosendx = momentum.x;
      tosendy = momentum.y;
      tosendz = momentum.z;
      MPI_Send(&kineticenergy,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(&tosendx,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      MPI_Send(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    }
  if (processid == 0)
    {
      MPI_Recv(&kineticenergy,1,MPI_DOUBLE,4,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendx,1,MPI_DOUBLE,4,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendy,1,MPI_DOUBLE,4,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendz,1,MPI_DOUBLE,4,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      momentum = vector(tosendx,tosendy,tosendz);
      energyout4 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
      energyout4 << t*h << "\t" << (2.0*kineticenergy)/3.0 << "\t" << kineticenergy << "\n";
    }

  // From 4th node
  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
  if (processid == 5)
    {
      tosendx = momentum.x;
      tosendy = momentum.y;
      tosendz = momentum.z;
      MPI_Send(&kineticenergy,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
      MPI_Send(&tosendx,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      MPI_Send(&tosendy,1,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
      MPI_Send(&tosendz,1,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    }
  if (processid == 0)
    {
      MPI_Recv(&kineticenergy,1,MPI_DOUBLE,5,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendx,1,MPI_DOUBLE,5,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendy,1,MPI_DOUBLE,5,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&tosendz,1,MPI_DOUBLE,5,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      momentum = vector(tosendx,tosendy,tosendz);
      energyout5 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
      energyout5 << t*h << "\t" << (2.0*kineticenergy)/3.0 << "\t" << kineticenergy << "\n";
      energyout0 << "Time\tTemperature\tKE\tEbond\tELJ\tEangle\tEDihedral\tEtotal\n";
    }


  /* 
    ***************************************************************************************************
    ****************************MAIN DYNAMICS LOOP*****************************************************
    NEEDS TO INCLUDE: 1. Dynamics on each node except the head node
                      2. Swap structures randomly whenever a particular multiple of a timestep is accessed.
		      3. Communicate the structure to the head node and write the output to file
    *****************************************************************************************************
   */
  for (t=0; t<neq; t++)
    {
      initmpi = MPI_Comm_rank(MPI_COMM_WORLD, &processid);
      if (processid == 0)
	{
	  // Reads from crystals.xyz on the head node, or outputs timestep to file
	  if (t%tscr == 0) cout << t+1 << "\n";
	}

      // Velocity Verlet Integration on all other nodes
      else
	{
	  if (t%tswap != 0)
	    {
	      for (int j=0; j<nchain; j++)
		{
		  for (int i=0; i<N; i++)
		    {
		      position[j][i] = add(position[j][i], multiply(velocity[j][i], h));
		      position[j][i] = add(position[j][i], multiply(acceleration[j][i], 0.5*h*h));
		      velocity[j][i] = add(velocity[j][i], multiply(acceleration[j][i], 0.5*h));
		    }
		}
	    }
	}

      //****************************SWAP STEP***************************************************
      // The swap step swaps structures at random between processes. The head node has
      // crystal structures in its memory, while the other nodes have their dynamics structures.
      // Structures and velocities need to be exchanged for these five nodes, on the basis of the
      // Metropolis algorithm.
      // The structure at node i will be exchanged with another in node j by
      // comparing 1 with exp(((1/kT_i) - (1/kT_j))(U_i - U_j)), and assigning that as the acceptance probability
      // Basically, if the minimum turns out to be 1, then the structures will be exchanged.
      // If the minimum turns out to be exp(((1/kT_i) - (1/kT_j))(U_i - U_j)), then a univariate random number
      // is drawn, and compared with exp(((1/kT_i) - (1/kT_j))(U_i - U_j)). The structures are exchanged if
      // the random number is less than exp(((1/kT_i) - (1/kT_j))(U_i - U_j)), and not exhanged if the random number is greater
      // than that.

      // To choose i and j as well, a univariate random number needs to be drawn. A range is made for each choice of
      // process. The random number is called twice and the two processes are chosen.

      if (t%tswap == 0)
       	{
	  MPI_Barrier(MPI_COMM_WORLD); // Makes all the nodes stop at t = (multiple of) tswap
	  //cout << "Barrier okay\n";
       	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
       	  // Task 1. Randomly arrange the processids using the random_shuffle method provided by algorithm
       	  if (processid == 0)
       	    {
	      totalenergy = 0.0;
	      ULJ = 0.0; Ubond = 0.0; Uangle = 0.0; Udihedral = 0.0;
       	      random_shuffle(&processids[0], &processids[nreplicas+1]);

       	      // Read crystal structure from file and calculate acceleration
	      getline(readcrystals,junk);
	      //cout << junk << "\n";
       	      getline(readcrystals,junk);
	      //cout << junk << "\n";
	      for (int j=0; j<nchain; j++)
		{
		  for (int k=0; k<N; k++)
		    {
		      readcrystals >> a >> position[j][k].x >> position[j][k].y >> position[j][k].z >> velocity[j][k].x >> velocity[j][k].y >> velocity[j][k].z;
		    }
		}
	      getline(readcrystals,junk);
	      for (int j=0; j<nchain; j++)
		{
		  calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
		  totalenergy += temptotalenergy;
		}
	      interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	      totalenergy += temptotalenergy;
	      writeenenergyoutput(position, velocity, energyout0, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	      //cout << "crystal file read correctly\n";
	      //cout << position[0].x << "\t" << position[N-1].x << "\n";
	      
	      // for (int k=0;k<(nreplicas+1);k++)
       	      // 	{
       	      // 	  cout << processids[k] << "\t";
       	      // 	}
       	      // cout << "\n";
       	    }
       	  for (int k=0;k<(nreplicas+1);k++)
       	    {
       	      MPI_Bcast(&processids[k],1,MPI_INT,0,MPI_COMM_WORLD);
       	    }

	  if (processid == 3)
	    {
	      cout << "Process " <<processid << "says okay\n";
	      for (int k=0;k<(nreplicas+1);k++)
       	   	{
       	   	  cout << processids[k] << "\t";
       	   	}
	      cout << "\n";
	    }

	  
       	  // Get the energy summation from the calculateacceleration function, and send both energies to
       	  // node 0 for comparison.
       	  // The first task is to find the node 0, and isolate it from the list of processids. Then, we
       	  // can swap across pairs separately, and then exchange the crystal structure separately
       	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD, &processid);
       	  if (processid == 0)
       	    {
       	      for (int k=0; k<(nreplicas+1); k++)
       		{
       		  if (processids[k] == 0)
       		    {
       		      valueofprocesszero = k;
       		      break;
       		    }
       		}
       	      if (valueofprocesszero%2 == 0) valueofneighborofzero = valueofprocesszero + 1;
       	      else valueofneighborofzero = valueofprocesszero - 1;
       	    }
       	  MPI_Bcast(&valueofprocesszero,1,MPI_INT,0,MPI_COMM_WORLD);
       	  MPI_Bcast(&valueofneighborofzero,1,MPI_INT,0,MPI_COMM_WORLD);
	  
       	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD, &processid);
	  if (processid == 3)
	    {
	      cout << "Process 0's index after shuffling is\t" << valueofprocesszero << " while its neighbor's index is\t"<< valueofneighborofzero << "\n";
	    }
       	  for (int k=0; k<(nreplicas+1); k=k+2)
       	    {
       	      //if (processid != processids[valueofprocesszero] && processid != processids[valueofneighborofzero])
	      //{
	      if (processid == processids[k+1])
		{
		  //cout << processids[k+1] << " sending" << totalenergy << " to process " << processids[k] << "\n";
		  MPI_Send(&totalenergy, 1 , MPI_DOUBLE, processids[k], 0, MPI_COMM_WORLD);
		}
	      if (processid == processids[k])
		{
		  MPI_Recv(&totalenergy2, 1, MPI_DOUBLE, processids[k+1], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  //cout << processids[k] << " received "<< totalenergy2 << " from " << processids[k+1] << "\n";
		  deltae12 = totalenergy - totalenergy2;
		  //cout << "delta e12 =" << deltae12 << "\n";
		  deltabeta12 = (1.0/kTreplica[processids[k]]) - (1.0/kTreplica[processids[k+1]]);
		  //cout << "deltabeta12 = " << deltabeta12 << "\n";
		  acceptanceprobability = exp(deltae12 * deltabeta12);
		  if (acceptanceprobability > 100.0) acceptanceprobability = 100;
		  //cout << acceptanceprobability << " is the acceptance probability on " << processid << " to be sent to " << processids[k+1] << "\n";
		  MPI_Send(&acceptanceprobability, 1, MPI_DOUBLE, processids[k+1], 1, MPI_COMM_WORLD);
		  //cout << "Sent " << acceptanceprobability << " to process " << processids[k+1] << "\n";
		  if (acceptanceprobability > 1 ) // The acceptance probability is now 1 and swapping can occur
		    {
		      // Send the coordinates and velocities to processids[k+1]
		      //cout << "Sending " << velocity[0].x << " to " << processids[k+1] << "from process " << processid << "\n";
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0; l<N; l++)
			    {
			      sendvector(position, tosendx, tosendy, tosendz, processids[k+1], p, l);
			      sendvector(velocity, vtosendx, vtosendy, vtosendz, processids[k+1], p, l);
			    }
			}
			  
			  
		      // Receive the incoming coordinates into positiontoswap, and velocity into velocitytoswap
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0;l<N;l++)
			    {
			      receivevector(positiontoswap, tosendx, tosendy, tosendz, processids[k+1], p, l);
			      receivevector(velocitytoswap, vtosendx, vtosendy, vtosendz, processids[k+1], p, l);
			    }
			}
		      //cout << "Received " << positiontoswap[0].x << "on process " << processid << " from process " << processids[k+1] << "\n";

		      // Assign the positiontoswap to the position vector, and velocitytoswap to velocity vector
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0; l<N; l++)
			    {
			      position[p][l] = vector(positiontoswap[p][l].x, positiontoswap[p][l].y, positiontoswap[p][l].z);
			      velocity[p][l] = vector(velocitytoswap[p][l].x, velocitytoswap[p][l].y, velocitytoswap[p][l].z);
			    }
			}
		      //cout << positiontoswap[0].x << " copied to " << position[0].x << " on process " << processid << "\n";
		      // Re-scale the velocities by the temperature factor
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0; l<N; l++) velocity[p][l] = multiply(velocity[p][l],sqrt(kTreplica[processids[k]]/kTreplica[processids[k+1]]));
			}
		      //cout << velocitytoswap[0].x << " re-scaled to " << velocity[0].x << " on process " << processid << "\n";
		    }
		  else
		    {
		      randomnumber1 = ur();
		      MPI_Send(&randomnumber1, 1, MPI_DOUBLE, processids[k+1], 2, MPI_COMM_WORLD);
		      cout << "Random number " << randomnumber1 << " sent from process " << processid << " to process " << processids[k+1] << "\n";
		      if (randomnumber1 <= acceptanceprobability)
			{
			  // Send the coordinates and velocities to processids[k+1]
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0; l<N; l++)
				{
				  sendvector(position, tosendx, tosendy, tosendz, processids[k+1], p, l);
				  sendvector(velocity, vtosendx, vtosendy, vtosendz, processids[k+1], p, l);
				}
			    }
 
			  // Receive the incoming coordinates into positiontoswap, and velocity to velocitytoswap
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0;l<N;l++)
				{
				  receivevector(positiontoswap, tosendx, tosendy, tosendz, processids[k+1], p, l);
				  receivevector(velocitytoswap, vtosendx, vtosendy, vtosendz, processids[k+1], p, l);
				}
			    }
			  // Assign the positiontoswap to the position vector, and the velocityswap to the velocity vector
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0; l<N; l++)
				{
				  position[p][l] = vector(positiontoswap[p][l].x, positiontoswap[p][l].y, positiontoswap[p][l].z);
				  velocity[p][l] = vector(velocitytoswap[p][l].x, velocitytoswap[p][l].y, velocitytoswap[p][l].z);
				}
			    }
				  
			      
			  // Re-scale the velocities by the temperature factor
			  for (int p=0;p<nchain; p++)
			    {
			      for (int l=0; l<N; l++) velocity[p][l] = multiply(velocity[p][l],sqrt(kTreplica[processids[k]]/kTreplica[processids[k+1]]));
			    }
			}
		    }
		}
	      if (processid == processids[k+1])
		{
		  //cout << "Preparing to accept acceptance probability on process " << processids[k+1] << " from process " << processids[k] << "\n";
		  MPI_Recv(&acceptanceprobability, 1, MPI_DOUBLE, processids[k], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  //cout << "Received " << acceptanceprobability << " from process " << processids[k] << "\n";
		  if (acceptanceprobability > 1)
		    {
		      // Receive the incoming coordinates from processids[k] and assign them to positions to swap
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0;l<N;l++)
			    {
			      receivevector(positiontoswap, tosendx, tosendy, tosendz, processids[k], p, l);
			      receivevector(velocitytoswap, vtosendx, vtosendy, vtosendz, processids[k], p, l);
			    }
			}
		      //cout << "Received " << velocitytoswap[0].x << " from process " << processids[k] << " on process " << processid << "\n" ;

		      // Send the positions and velocities back to processids[k]
		      for (int p=0;p<nchain; p++)
			{
			  for (int l=0;l<N;l++)
			    {
			      sendvector(position, tosendx, tosendy, tosendz, processids[k], p, l);
			      sendvector(velocity, vtosendx, vtosendy, vtosendz, processids[k], p, l);
			    }
			}
		      //cout << "Sending " << position[0].x << " to process " << processids[k] << " from process " << processid << "\n";

		      // Assign the positiontoswap to the position vector, and the velocitytoswap to the velocity vector
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0; l<N; l++)
			    {
			      position[p][l] = vector(positiontoswap[p][l].x, positiontoswap[p][l].y, positiontoswap[p][l].z);
			      velocity[p][l] = vector(velocitytoswap[p][l].x, velocitytoswap[p][l].y, velocitytoswap[p][l].z);
			    }
			}
		      //cout << positiontoswap[0].x << " copied to " << position[0].x << " on process " << processid << "\n";
		      // Re-scale the velocities by a factor of square root of ratio of temperatures
		      for (int p=0; p<nchain; p++)
			{
			  for (int l=0; l<N; l++) velocity[p][l] = multiply(velocity[p][l],sqrt(kTreplica[processids[k+1]]/kTreplica[processids[k]]));
			}
		      //cout << velocitytoswap[0].x << " re-scaled to " << velocity[0].x << " on process " << processid << "\n";
		    }
		  else
		    {
		      // Receive the randomnumber, which is to be measured against the acceptance probability
		      MPI_Recv(&randomnumber1, 1, MPI_DOUBLE, processids[k], 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      cout << "Random number " << randomnumber1 << " received from process " << processids[k] << " on process " << processid << "\n";

		      // Measure the randomnumber against the probability that was received before these if-else cases
		      if (randomnumber1 <= acceptanceprobability)
			{
			  // If the randomnumber is <= acceptance probability then receive the coordinates into positiontoswap
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0;l<N;l++)
				{
				  receivevector(positiontoswap, tosendx, tosendy, tosendz, processids[k], p, l);
				  receivevector(velocitytoswap, vtosendx, vtosendy, vtosendz, processids[k], p, l);
				}
			    }

			  // Send the positions and velocities back to processids[k]
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0;l<N;l++)
				{
				  sendvector(position, tosendx, tosendy, tosendz, processids[k], p, l);
				  sendvector(velocity, vtosendx, vtosendy, vtosendz, processids[k], p, l);
				}
			    }

			  // Assign the positiontoswap to position vectors
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0; l<N; l++)
				{
				  position[p][l] = vector(positiontoswap[p][l].x, positiontoswap[p][l].y, positiontoswap[p][l].z);
				  velocity[p][l] = vector(velocitytoswap[p][l].x, velocitytoswap[p][l].y, velocitytoswap[p][l].z);
				}
			    }

			  // Re-scale the velocities by a factor of square root of ratio of temperatures
			  for (int p=0; p<nchain; p++)
			    {
			      for (int l=0; l<N; l++) velocity[p][l] = multiply(velocity[p][l],sqrt(kTreplica[processids[k+1]]/kTreplica[processids[k]]));
			    }
			}
		    }
		}       		
       	    }
       	}
      
      

      // Calculation of acceleration for the next step. Has to be done individually
      totalenergy = 0.0;
      ULJ = 0.0; Ubond = 0.0; Uangle = 0.0; Udihedral = 0.0;
      initmpi = MPI_Comm_rank(MPI_COMM_WORLD, &processid);
      if (processid == 1)
	{
	  for (int j=0;j<nchain;j++)
	    {
	      calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	      totalenergy += temptotalenergy;
	    }
	  interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	  totalenergy += temptotalenergy;
	}
      if (processid == 2)
	{
	  for (int j=0;j<nchain;j++)
	    {
	      calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	      totalenergy += temptotalenergy;
	    }
	  interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	  totalenergy += temptotalenergy;
	}
      if (processid == 3)
	{
	  for (int j=0;j<nchain;j++)
	    {
	      calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	      totalenergy += temptotalenergy;
	    }
	  interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	  totalenergy += temptotalenergy;
	}
      if (processid == 4)
	{
	  for (int j=0;j<nchain;j++)
	    {
	      calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	      totalenergy += temptotalenergy;
	    }
	  interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	  totalenergy += temptotalenergy;
	}
      if (processid == 5)
	{
	  for (int j=0;j<nchain;j++)
	    {
	      calculateacceleration(position, velocity, acceleration, j, temptotalenergy, kTreplica[processid], ULJ, Ubond, Uangle, Udihedral);
	      totalenergy += temptotalenergy;
	    }
	  interchainLJ(position, acceleration, temptotalenergy); ULJ += temptotalenergy;
	  totalenergy += temptotalenergy;
	}


      // Outputs need to be written to the respective files
      // 1. VMD Outputs
      if(t % tvmd == 0)
	{
	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
	  if (processid == 0) writevmdoutput(position, velocity, vmdwrite0);
	  if (processid == 1) writevmdoutput(position, velocity, vmdwrite1);
	  if (processid == 2) writevmdoutput(position, velocity, vmdwrite2);
	  if (processid == 3) writevmdoutput(position, velocity, vmdwrite3);
	  if (processid == 4) writevmdoutput(position, velocity, vmdwrite4);
	  if (processid == 5) writevmdoutput(position, velocity, vmdwrite5);
	}

      // 2.Write Energy, Momentum to Energyout File, first by defining center of mass velocity. Write energy output does all this
      if (t % tener == 0)
	{
	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
	  if (processid == 1) writeenenergyoutput(position, velocity, energyout1, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	  if (processid == 2) writeenenergyoutput(position, velocity, energyout2, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	  if (processid == 3) writeenenergyoutput(position, velocity, energyout3, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	  if (processid == 4) writeenenergyoutput(position, velocity, energyout4, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	  if (processid == 5) writeenenergyoutput(position, velocity, energyout5, momentum, centerofmass, centerofmassvelocity, kineticenergy, Ubond, ULJ, Uangle, Udihedral, totalenergy);
	}

      // 3. Write radii output to file
      if (t % trad == 0)
	{
	  initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid);
	  if (processid == 1) writeradiioutput(position, velocity, radiiout1, centerofmass, rg, rendtoend);
	  if (processid == 2) writeradiioutput(position, velocity, radiiout2, centerofmass, rg, rendtoend);
	  if (processid == 3) writeradiioutput(position, velocity, radiiout3, centerofmass, rg, rendtoend);
	  if (processid == 4) writeradiioutput(position, velocity, radiiout4, centerofmass, rg, rendtoend);
	  if (processid == 5) writeradiioutput(position, velocity, radiiout5, centerofmass, rg, rendtoend);
	}
    }

  // ************************END OF EQUILIBRATION LOOP *************************

  // ************************FINALIZE MPI**************************************
  MPI_Finalize();
  // **************************************************************************
  
  // ********************Deallocate Variables**********************************
  //vmdwrite.close();
  //velwrite.close();
  //energyout.close();
  //radiiout.close();
  //dispout.close();
  //opwrite.close();
  for (int i=0;i<nchain;i++)
    {
      delete [] position[i];
      delete [] velocity[i];
      delete [] acceleration[i];
      delete [] positiontoswap[i];
      delete [] velocitytoswap[i];
    }
  delete [] positiontoswap;
  delete [] velocitytoswap;
  readcrystals.close();

  vmdwrite1.close();
  vmdwrite2.close();
  vmdwrite3.close();
  vmdwrite4.close();
  vmdwrite5.close();
  energyout1.close();
  energyout2.close();
  energyout3.close();
  energyout4.close();
  energyout5.close();
  radiiout1.close();
  radiiout2.close();
  radiiout3.close();
  radiiout4.close();
  radiiout5.close();
  opwrite1.close();
  opwrite2.close();
  opwrite3.close();
  opwrite4.close();
  opwrite5.close();    
  delete [] position;
  delete [] velocity;
  delete [] acceleration;
  delete [] processids;
  //delete [] initialposition;
  //delete [] trueposition;
  delete gacc;
  delete sysacc;

  return 0;
}

// *************************Define Functions**************************

void calculateacceleration(vector** position, vector** velocity, vector** acceleration, int chainindex, double &tote, double kT, double &ULJ, double &Ubond, double &Uangle, double &Udihedral)
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
      //friacc += (-magnitude(velocity[i])*zeta/mass);
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
  Ubond += ebond;
  //bacc /= N;
  
  // Calculate the forces arising out of the Lennard Jones 12-6 potential
  // U = epsilon * ((sigma/r)^12 - 2(sigma/r)^6)
  // Force is therefore, (-12 epsilon/r) * ((sigma/r)^12 - (sigma/r)^6)
  
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
  ULJ += elj;
  //ljacc /= N;
  //cout << ranacc << "\t" << friacc << "\t" << bacc << "\t" << ljacc << "\t" ;

  // Calculate the forces due to the bond angle potential
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
  Uangle += eangle;

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

  //Calculate the forces arising out of the dihedral angle potential
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
  Udihedral += edihedral;
  
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
  tote = etotal;
}

void writevmdoutput(vector** position, vector** velocity, ofstream &vmdwrite)
{
  vmdwrite << N*nchain << endl << "Atom type\tx\ty\tz\tvx\tvy\tvz" << endl;
  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  vmdwrite << "C\t" << position[j][i].x << "\t" << position[j][i].y << "\t" << position[j][i].z << "\t" << velocity[j][i].x << "\t" <<velocity[j][i].y << "\t" << velocity[j][i].z << endl;
	}
    }
}
  
void writeenenergyoutput(vector** position, vector** velocity, ofstream &energyout, vector momentum, vector centerofmass, vector centerofmassvelocity, double kineticenergy, double ebond, double elj, double eangle, double edihedral, double etotal)
{
  centerofmassvelocity = vector(0.0,0.0,0.0);
  for (int j=0; j<nchain; j++)
    {
      for (int i=0;i<N;i++)
	{
	  centerofmassvelocity = add(centerofmassvelocity,velocity[j][i]);
	}
    }
  centerofmassvelocity = divide(centerofmassvelocity,(double)(N*nchain));
	  
  kineticenergy = 0.0;
  momentum = vector(0,0,0);
  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  kineticenergy += 0.5*mass*pow(magnitude(subtract(velocity[j][i],centerofmassvelocity)), 2);
	  momentum = add(momentum, multiply(subtract(velocity[j][i],centerofmassvelocity), mass));
	}
    }
  kineticenergy /= (double)(N*nchain);

  energyout << (t+1)*h << "\t" << (2.0/3.0) * kineticenergy <<"\t"<< kineticenergy << "\t" ;
  energyout << setw(6)<<ebond << "\t" <<  setw(6)<<elj << "\t" << setw(6)<<eangle << "\t" << setw(6)<<edihedral << "\t" << setw(6)<<etotal << "\n";
}

void writeradiioutput(vector** position, vector** velocity, ofstream &radiiout, vector centerofmass, double rg, double rendtoend)
{
  // Calculate Center of Mass Position
  centerofmass = vector(0,0,0);
  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  centerofmass = add(centerofmass, position[j][i]);
	}
    }
  centerofmass = divide(centerofmass, (double)(N*nchain));

  // Write Chain Size to radii
  rg = 0;
  rendtoend = magnitude( subtract(position[0][N-1], position[nchain-1][0]) );

  for (int j=0; j<nchain; j++)
    {
      for(int i=0; i<N; i++)
	{
	  rg += pow( magnitude(subtract(position[j][i], centerofmass)), 2 );
	}
    }
  rg /= (double)((N*nchain)+1);
  rg = sqrt(rg);
  radiiout << (t+1)*h << "\t" << rg << "\t" << rendtoend << endl;
}

void interchainLJ(vector** position, vector** acceleration, double &tote)
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
  tote = elj;
}
