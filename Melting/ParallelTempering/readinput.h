#ifndef _READINPUT_H_
#define _READINPUT_H_

extern int nchain = 21;                 // Number of chains
extern int N = 200;			// chain length
extern double mass = 1.0;		// bead mass
extern double epsilon = 1.0;	// LJ well depth
extern double sigma1 = 1.0;	// LJ separation distance
extern double sigma2 = 2.941558442; //LJ separation for pairs farther than 1-5
extern double k = 2435.125;		 // harmonic spring constant
extern double ktheta = 535.7142857; //Cosine squared angle potential constant
extern double cs0 = 0.3255681545; //cos(theta0)
extern double b = 1.0;		// harmonic equilibrium bond distance
extern double Tstar = 12.0; // T* required for the equilibration
extern double Tstar2 = 9.0; // T* required for the crystallization
extern double hrate = 0.001; //Heating rate required
extern double kTreplica[6] = {9.0, 11.6, 11.8, 12.0, 12.2, 12.4}; // Array of the temperatures of the replicas in the system
extern int nreplicas = 5;                // Number of replicas of the system
extern double zeta = 1.0;		// drag coefficient
extern double k1 = 26.96428571;  //Coefficient of first term in dihedral potential
extern double k2 = -5.0;  //Coefficient of second term in dihedral potential
extern double k3 = 23.03571429; //Coefficient of third term in dihedral potential
  

extern double h = 0.001;	// step size
extern int neq = 5000000;	// no. time iterations for equilibration at high temperature
extern int nquench = 0; // No. Time iterations for quench
extern int nheat = 0; // No. Time iterations required for heating
extern double density = 1e-4;	// particle number density

extern double rcutoff = 6.5; // Cut-off distance for application of Lennard jones potential
extern double L = pow((N/density), (1.0/3.0)); // Box size

extern int tswap = 50000; // Frequency of swapping structures
extern int tvmd = 1000; //Output frequency for VMD
extern int tdisp = 10000; //Output frequency for display file
extern int tener = 10000000; //Output frequency for energy file
extern int tscr = 1000; //Output frequency for screen printing
extern int top = 1000000; //Output frequency for order parameter analysis in the heating phase
extern int trad = 1000000; //Output frequency for radius of gyration calculation
#endif
