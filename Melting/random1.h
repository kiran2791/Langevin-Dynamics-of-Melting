#ifndef _RANDOM1_H_
#define _RANDOM1_H_

// returns uniform distribution variate between 0 and 1
double ur()
{
  return (double)rand()/(double)RAND_MAX;
}

// returns standard Gaussian distribution variate (Joseph Leva)
double nr()
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if(iset==0)
    {
      do
	{
	  v1=2.0*ur()-1.0;
	  v2=2.0*ur()-1.0;
	  rsq=v1*v1+v2*v2;
	} while(rsq >=1.0 || rsq==0.0);

      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return(v2*fac);
    }
  else
    {
      iset=0;
      return gset;
    }
}

/*double nr(){
    double u, v, x, y, q;
    bool run = true;
    while(run == true){
        u = ur();
        v = 1.7156*(ur()-0.5);
        x = u - 0.449871;
        y = abs(v) + 0.386595;
        q = x*x + y*(0.19600*y-0.25472*x);
        if(q<0.27597){
            run = false;
        } else if(q<0.27846){
            if(v*v<-4.0*log(u)*u*u) run = false;
        }
    }
    return v/u;
}*/

#endif
