#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_

void initialization(vector *position, vector *velocity)
{
  double velocity_sigma = sqrt(kT/mass);

  // initialize monomer positions
  /*double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double x,y,z,r1,r2;
  srand(time(NULL));
  position[0] = vector(x0,y0,z0);
  for (int i=1;i<N;i++)
    {
      r1 = M_PI_2*(rand()/double(RAND_MAX));
      r2 = M_PI*(rand()/double(RAND_MAX));
      x = x0 + (b*cos(r1)*sin(r2));
      y = y0 + (b*sin(r1)*sin(r2));
      z = z0 + (b*cos(r2));

      x0 = x;
      y0 = y;
      z0 = z;
      position[i] = vector(x0,y0,z0);
    }
  */
  
  //double inc = b*0.01;
  /*for(int i=0; i<N; i++)
    {
      position[i] = vector(inc*(i+1), (L+0.01)/2.0, (L-0.01)/2.0);
    }
  */

  // Initialize monomer positions. Read from specific file
  std::ifstream fin;
  fin.open("last.xyz");
  std::string junk;
  getline(fin,junk);
  getline(fin,junk);
  char j;
  for (int k=0;k<N;k++)
    {
      fin>>j>>position[k].x>>position[k].y>>position[k].z>>velocity[k].x>>velocity[k].y>>velocity[k].z;
      //std::cout<<position[k].x<<"\t"<<position[k].y<<"\t"<<position[k].z<<"\n";
    }
  
  // initialize monomer velocities
  /* for(int i=0; i<N; i++) */
  /*   { */
  /*     velocity[i] = vector( velocity_sigma*nr(), velocity_sigma*nr(), velocity_sigma*nr() ); */
  /*   } */
  /* fin.close(); */
}

#endif
