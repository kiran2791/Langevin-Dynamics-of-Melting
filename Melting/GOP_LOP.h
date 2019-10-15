#ifndef _GOP_LOP_H_
#define _GOP_LOP_H_

void OP(vector** pos,double *gacc,double *sysacc)
{
  vector gl[nchain][N-2];
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
  for (int j=0; j<nchain; j++)
    {
      for (int i=1;i<=N-2;i++)
	{
	  gl[j][i-1] = subtract(pos[j][i+1],pos[j][i-1]);
	}
    }

  //Calculation of the global order parameter for this time step.
  // Only averaged across all chain atoms, but calculated only from each single chain.
  for (int k=0; k<nchain; k++)
    {
      for(int i=0;i<N-3;i++)
	{
	  for(int j=i+1 ; j<N-2 ; j++)
	    {
	      cs2phi = dot(gl[k][i],gl[k][j])*dot(gl[k][i],gl[k][j]);
	      cs2phi /= ((gl[k][i].x*gl[k][i].x + gl[k][i].y*gl[k][i].y + gl[k][i].z*gl[k][i].z)*(gl[k][j].x*gl[k][j].x + gl[k][j].y*gl[k][j].y + gl[k][j].z*gl[k][j].z));
	      *gacc += (3.0 * cs2phi - 1.0);
	      gcount ++;
	    }
	}
    }
  *gacc /= (2.0 * gcount);

  //*****************END OF GLOBAL ORDER PARAMETER CALCULATION*********************

  //**************LOCAL ORDER PARAMETER CALCULATION********************************
  double drsq;
  int aindex[50];
  int indexcount;
  for (int l=0; l<nchain; l++)
    {
      for (int i=0 ; i<N ; i++)
	{
	  lacc = 0.0;
	  lcount = 0;
	  indexcount = 0;
	  for(int j=1 ; j<N-1; j++)
	    {
	      if (j!=i)
		{
		  dr = subtract(pos[l][i],pos[l][j]);
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
		  cs2phi = dot(gl[l][aindex[j]],gl[l][aindex[k]])*dot(gl[l][aindex[j]],gl[l][aindex[k]]);
		  cs2phi /= ((gl[l][aindex[j]].x*gl[l][aindex[j]].x + gl[l][aindex[j]].y*gl[l][aindex[j]].y + gl[l][aindex[j]].z*gl[l][aindex[j]].z)
			     *(gl[l][aindex[k]].x*gl[l][aindex[k]].x + gl[l][aindex[k]].y*gl[l][aindex[k]].y + gl[l][aindex[k]].z*gl[l][aindex[k]].z));
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
	  //fout << lcount << "\n";
	  *sysacc += (lacc/lcount);
	  glcount ++;
	}
      lacc /= lcount;
      //fout.close();
    }
  *sysacc /= glcount;
}

#endif
