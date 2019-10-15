#ifndef _MPIUTILITIES_H_
#define _MPIUTILITIES_H_

// MPI function to send a vector. Copies the contents of the vector
// Then executes a send command for all three vector components
void sendvector(vector** a, double tosendx, double tosendy, double tosendz, int recipient, int chainid, int index)
{
  tosendx = a[chainid][index].x;
  tosendy = a[chainid][index].y;
  tosendz = a[chainid][index].z;
  MPI_Send(&tosendx,1,MPI_DOUBLE,recipient,1,MPI_COMM_WORLD);
  MPI_Send(&tosendy,1,MPI_DOUBLE,recipient,2,MPI_COMM_WORLD);
  MPI_Send(&tosendz,1,MPI_DOUBLE,recipient,3,MPI_COMM_WORLD);
  //std::cout << "Sending" << tosendx << "\t" << tosendy << "\t" << tosendz << "\tto process" << recipient << "\n";
}

// MPI function to receive a vector. Receives the contents into
// double variables. Then, it assigns the received variables
// to the vector
void receivevector(vector** a, double tosendx, double tosendy, double tosendz, int sender, int chainid, int index)
{
   MPI_Recv(&tosendx,1,MPI_DOUBLE,sender,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Recv(&tosendy,1,MPI_DOUBLE,sender,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   MPI_Recv(&tosendz,1,MPI_DOUBLE,sender,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   a[chainid][index].x = tosendx;
   a[chainid][index].y = tosendy;
   a[chainid][index].z = tosendz;
}


// MPI function to swap an array of pointers between two vectors. The function requires
// the process IDs of the two nodes in question, then the vector pointers to each
// vector being sent. The function then uses node 0 as a buffer to collect
// the vector from the sender. Then the sender receives the vector from the
// recipient node. The node 0 then sends the sender's vector to the recipient.

/* void swapvectors(vector *position, vector *velocity,double tosendx, double tosendy, double tosendz, int sender, int recipient) */
/* { */
/*   //Identify the sender node in question. Then send the position and velocity */
/*   // to node 0 */
/*   int initmpi,numberofprocesses,processid; */
/*   initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid); */
/*   if (processid == sender) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  sendvector(position, tosendx, tosendy, tosendz, 0, k); */
/* 	  sendvector(velocity, vtosendx, vtosendy, vtosendz, 0, k); */
/* 	} */
/*     } */
/*   if (processid == 0) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  receivevector(positiontoswap,tosendx,tosendy,tosendz,sender,k); */
/* 	  receivevector(velocitytoswap, vtosendx, vtosendy, vtosendz, sender, k); */
/* 	} */
/*     } */

/*   // The sender has now sent the positions and the velocities to the head node */
/*   // The recipient can now send the positions and velocities to the sender node */
/*   MPI_Barrier(MPI_COMM_WORLD); */
/*   initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid); */
/*   if (processid == recipient) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  sendvector(position, tosendx, tosendy, tosendz, sender, k); */
/* 	  sendvector(velocity, vtosendx, vtosendy, vtosendz, sender, k); */
/* 	} */
/*     } */
/*   if (processid == sender) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  receivevector(position,tosendx,tosendy,tosendz,recipient,k); */
/* 	  receivevector(velocity, vtosendx, vtosendy, vtosendz, recipient, k); */
/* 	} */
/*     } */

/*   // The recipient has now sent the positions and the velocities to the sender */
/*   // The head node can now send positions and vectors to the recipient node */
/*   MPI_Barrier(MPI_COMM_WORLD); */
/*   initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid); */
/*   if (processid == 0) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  sendvector(positiontoswap, tosendx, tosendy, tosendz, recipient, k); */
/* 	  sendvector(velocitytoswap, vtosendx, vtosendy, vtosendz, recipient, k); */
/* 	} */
/*     } */
/*   if (processid == recipient) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  receivevector(position,tosendx,tosendy,tosendz,0,k); */
/* 	  receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, k); */
/* 	} */
/*     } */
/* } */

/* /\* ************************************************************************ */
/*     This  function swaps monomers in node1 with the crystal structures maintained */
/* in node 0. The node 0 need not receive a structure in return. */
/*     ***********************************************************************   *\/ */
/* void swapvectorwithcrystal(vector *a, double tosendx, double tosendy, double tosendz, int recipient, int index) */
/* { */
/*   // Identify the node1 in question and get the head node to send it the crystal structure */
/*   int initmpi,numberofprocesses,processid; */
/*   initmpi = MPI_Comm_rank(MPI_COMM_WORLD,&processid); */
/*   if (processid == 0) */
/*     { */
/*       for (int k=0;k<N;k++) */
/* 	{ */
/* 	  sendvector(position, tosendx, tosendy, tosendz, recipient, k); */
/* 	  sendvector(velocity, vtosendx, vtosendy, vtosendz, recipient, k); */
/* 	} */
/*     } */
/*   if (processid == recipient) */
/*     { */
/*       for (int k=0; k<N; k++) */
/* 	{ */
/* 	  receivevector(position,tosendx,tosendy,tosendz,0,k); */
/* 	  receivevector(velocity, vtosendx, vtosendy, vtosendz, 0, k); */
/* 	} */
/*     } */
/* } */
      


#endif
