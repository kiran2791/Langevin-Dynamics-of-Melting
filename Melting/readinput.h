#ifndef _READINPUT_H_
#define _READINPUT_H_

extern int nchain = 21;                  // Number of chains
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
extern double Tstar2 = 12.0; // T* required for the quench
extern double hrate = 0.0; //Heating rate required
extern double kT = Tstar*epsilon;		// thermal energy
extern double zeta = 1.0;		// drag coefficient
extern double k1 = 26.96428571;  //Coefficient of first term in dihedral potential
extern double k2 = -5.0;  //Coefficient of second term in dihedral potential
extern double k3 = 23.03571429; //Coefficient of third term in dihedral potential
  

extern double h = 0.001;	// step size
extern int neq = 0;	// no. time iterations for equilibration at high temperature
extern int nquench = 5000000; // No. Time iterations for quench
extern int nheat = 0; // No. Time iterations required for heating
extern double density = 1e-4;	// particle number density

extern double rcutoff = 6.5; // Cut-off distance for application of Lennard jones potential
extern double L = pow((N/density), (1.0/3.0)); // Box size

extern int tvmd = 1000; //Output frequency for VMD
extern int tdisp = 1000; //Output frequency for display file
extern int tener = 1000; //Output frequency for energy file
extern int tscr = 1000; //Output frequency for screen printing
extern int top = 100; //Output frequency for order parameter analysis in the heating phase
#endif
