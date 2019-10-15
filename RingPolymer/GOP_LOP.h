#ifndef _GOP_LOP_H_
#define _GOP_LOP_H_

void OP(vector *pos,double *gacc,double *sysacc)
{
  vector gl[N-2];
  vector dr;
  //std::ofstream fout;
  //fout.open ("lcount.dat");

  double cs2phi,lacc;

  double xmax,xmin,ymax,ymin,zmax,zmin;
  long int gcount,lcount,glcount;

  *gacc = 0.0;
  *sysacc = 0.0;
  gcount = 0;
  glcount = 0;
  
  //**************GLOBAL ORDER PARAMETER CALCULATION***********************
  //Calculation of orientation vectors
  for (int i=1;i<=N-2;i++)
    {
      gl[i-1] = subtract(pos[i+1],pos[i-1]);
    }

  //Calculation of the global order parameter for this time step
  for(int i=0;i<N-3;i++)
    {
      for(int j=i+1 ; j<N-2 ; j++)
	{
	  cs2phi = dot(gl[i],gl[j])*dot(gl[i],gl[j]);
	  cs2phi /= ((gl[i].x*gl[i].x + gl[i].y*gl[i].y + gl[i].z*gl[i].z)*(gl[j].x*gl[j].x + gl[j].y*gl[j].y + gl[j].z*gl[j].z));
	  *gacc += (3.0 * cs2phi - 1.0);
	  gcount ++;
	}
    }
  *gacc /= (2.0 * gcount);

  //*****************END OF GLOBAL ORDER PARAMETER CALCULATION*********************

  //**************LOCAL ORDER PARAMETER CALCULATION********************************
  double drsq;
  int aindex[50];
  int indexcount;
  for (int i=0 ; i<N ; i++)
    {
      lacc = 0.0;
      lcount = 0;
      indexcount = 0;
      for(int j=1 ; j<N-1; j++)
	{
	  if (j!=i)
	    {
	      dr = subtract(pos[i],pos[j]);
	      drsq = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
	      if (drsq <= sigma2*sigma2)
		{
		  aindex[indexcount] = j;
		  //std::cout << aindex[indexcount] << " ";
		  indexcount ++ ;
		}
	    }
	}
      
      for (int j=0; j<indexcount-1 ; j++)
	{
	  //std::cout<<"\n";
	  for(int k=j+1; k<indexcount ; k++)
	    {
	      cs2phi = dot(gl[aindex[j]],gl[aindex[k]])*dot(gl[aindex[j]],gl[aindex[k]]);
	      cs2phi /= ((gl[aindex[j]].x*gl[aindex[j]].x + gl[aindex[j]].y*gl[aindex[j]].y + gl[aindex[j]].z*gl[aindex[j]].z)
			 *(gl[aindex[k]].x*gl[aindex[k]].x + gl[aindex[k]].y*gl[aindex[k]].y + gl[aindex[k]].z*gl[aindex[k]].z));
	      if (isnan(cs2phi))
		{
		  cs2phi = 0.0;
		}
	      lacc += ((3.0 * cs2phi - 1.0)/2.0);
	      lcount ++;
	      //std::cout << cs2phi << " ";
	    }
	}

      //fout.seekp(std::ios_base::beg);
      lacc /= lcount;
      //fout << lcount << "\n";
      *sysacc += lacc;
      glcount ++;
    }
  //fout.close();
  *sysacc /= glcount;
}

#endif
