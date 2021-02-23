#include <iostream>
#include <cmath>
#include <fstream>
#include <mpi.h>

using namespace std;

int main () {
  // Declare integration parameters
  int M = 602; // Number spacial points
  double Md = M;
  double dx = 2/(Md-2), dt = 0.2/(Md-2);
  double dtdx = dt*dt/(dx*dx);  // Used in calculations below
  double T = 1;


  // Initialise MPI interface
  int rank, size;
  MPI_Init(NULL, NULL);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
  int sizeSq = pow(size, 0.5);
  int J = (M-2)/sizeSq + 2;  // Edge length for each process
  double rankd = rank, Jd = J;
  // Initialize solution
  double un[J][J] = {};
  double unm1[J][J] = {};
  double utemp[J][J] = {};
  double uswap[J] = {};
  double uswapOut0[J] = {}, uswapOut1[J] = {};
  ofstream fileOut, initialOut, timeOut;  // Used to save files
  double uOut [J][M], uOutJ[J][J];
  if (rank == 0){
    initialOut.open("initial-cond.txt");
    
  }

  
  // State start and end points in x and y direction away from boundary points
  int startx = 0, endx = J;
  if (rank <  sizeSq){startx = 1;}  
  if (rank >= size-sizeSq){endx = J-1;}

  int starty = 0, endy = J;
  if (rank%sizeSq == 0){starty = 1;}  
  if ((rank+1)%sizeSq == 0){endy = J-1;}

  // Start time
  MPI_Barrier(comm);
  double t_start = MPI_Wtime();
  double time_write = 0;  // Store time spent writing to files
  int posx = (rank - rank%sizeSq)*(J-2)/sizeSq, posy = (rank%sizeSq)*(J-2);
  for (int i=startx; i<endx; ++i){
    for (int j=starty; j<endy; ++j){
      un[i][j] = exp(-40*(pow((posx+i)*dx - 1.4, 2) + pow((posy + j)*dx - 1,2)));
      unm1[i][j] = un[i][j];
    }
  }
  double tWrite0 = MPI_Wtime();
  if (rank > 0){MPI_Send(un, J*J, MPI_DOUBLE, 0, 0, comm);}
  if (rank == 0){
      for (int i=0; i<J; ++i){
	for (int j=0; j<J; ++j){
	  uOut[i][j] = un[i][j];
	}
      }
  
    for (int r = 1; r<size; ++r){
      int startpointx = 2, startpointy = 2;
      if (r < sizeSq){startpointx = 0;}
      if (r%sizeSq > 0.5){startpointy = 0;}
      int pos = (r%sizeSq)*J;
      MPI_Recv(uOutJ, J*J, MPI_DOUBLE, r, 0, comm, MPI_STATUS_IGNORE);
      for (int i=0; i<J; ++i){
	for (int j=pos; j < pos + J; ++j){
	  uOut[i][j] = uOutJ[i][j - pos];
	}
      }
      if ((r+1)%sizeSq == 0){
	for (int i=startpointx; i < J; ++i){
	  for (int j=startpointy; j < M; ++j){
	    initialOut << uOut[i][j] << " ";
	  }
	  initialOut << "\n";
	}
      }
    }
    initialOut.close();
  }
  double tWrite1 = MPI_Wtime();
  time_write += tWrite1 - tWrite0;
  
  // Iterate solution for 2 <= n <= N
  for (double t = 2*dt; t < T; t += dt){
    if (rank == 0){cout << "\n t = "<< t;}
    // perform update
    for (int i=1; i < J-1; ++i){
      for (int j=1; j < J-1; ++j){
	utemp[i][j] = dtdx*(un[i+1][j] + un[i-1][j] + un[i][j+1]
			    + un[i][j-1] - 4*un[i][j]) + 2*un[i][j] - unm1[i][j];
	  }
    }
    // Update parameters
    for (int i=1; i < J-1; ++i){
      for (int j=1; j < J-1; ++j){
	unm1[i][j] = un[i][j];
	un[i][j] = utemp[i][j];
      }
    }

    // halo-swapping
    // Send
    if (rank >= sizeSq){
      MPI_Send(un[1], J, MPI_DOUBLE, rank - sizeSq, 0, comm);
    }
    if (rank < size - sizeSq){
      MPI_Send(un[J-2], J, MPI_DOUBLE, rank + sizeSq, 1, comm);
    }
    if (rank%sizeSq > 0.5){
      for (int j=0; j < J; ++j){uswapOut0[j] = un[j][1];}
      MPI_Send(uswapOut0, J, MPI_DOUBLE, rank - 1, 2, comm);
    }
    if ((rank + 1)%sizeSq > 0.5){
      for (int j=0; j < J; ++j){uswapOut1[j] = un[j][J-2];}
      MPI_Send(uswapOut1, J, MPI_DOUBLE, rank + 1, 3, comm);
    }
  
    
    // Receive
    if (rank >= sizeSq){
      MPI_Recv(uswap, J, MPI_DOUBLE, rank - sizeSq, 1, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < J; ++m){un[0][m] = uswap[m];}
    }
    if (rank < size - sizeSq){
      MPI_Recv(uswap, J, MPI_DOUBLE, rank + sizeSq, 0, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < J; ++m){un[J-1][m] = uswap[m];}
    }
    if ((rank+1)%sizeSq > 0.5){
      MPI_Recv(uswap, J, MPI_DOUBLE, rank + 1, 2, comm, MPI_STATUS_IGNORE);
      for (int j=0; j < J; ++j){un[j][J-1] = uswap[j];}
    }
    if (rank%sizeSq > 0.5){
      MPI_Recv(uswap, J, MPI_DOUBLE, rank - 1, 3, comm, MPI_STATUS_IGNORE);
      for (int j=0; j < J; ++j){un[j][0] = uswap[j];}
    }
    

    // Save data if t = 1/3, 2/3, 1
    
    if ((pow(1/1. - t,2) <dt*dt/100)|(pow(1/3. - t,2) <dt*dt/100)|(pow(2/3. - t,2) <dt*dt/100)){
      // Time process of writing
      tWrite0 = MPI_Wtime();
      int tag = t*3;
      if (pow(1/3. - t,2) <dt*dt/100){fileOut.open("data1.txt");}
      if (pow(2/3. - t,2) <dt*dt/100){fileOut.open("data2.txt");}
      if (pow(3/3. - t,2) <dt*dt/100){fileOut.open("data3.txt");}
      if (rank>0){MPI_Send(un, J*J,  MPI_DOUBLE, 0, tag, comm);}
      if (rank == 0){
	for (int i=0; i < J; ++i){
	  for (int j=0; j < J; ++j){
	    uOut[i][j] = un[i][j];
	  }
	}
      for (int r = 1; r<size; ++r){
	int startpointx = 2, startpointy = 2;
	if (r < sizeSq){startpointx = 0;}
	if (r%sizeSq > 0.5){startpointy = 0;}
	int pos = (r%sizeSq)*J;
	MPI_Recv(uOutJ, J*J, MPI_DOUBLE, r, tag, comm, MPI_STATUS_IGNORE);
	for (int i=0; i<J; ++i){
	  for (int j=pos; j < pos + J; ++j){
	    uOut[i][j] = uOutJ[i][j - pos];
	  }
	}
	if ((r+1)%sizeSq == 0){
	  for (int i=startpointx; i < J; ++i){
	    for (int j=startpointy; j < M; ++j){
	      fileOut << uOut[i][j] << " ";
	    }
	    fileOut << "\n";
	  }
	}
      }
      fileOut.close();
      tWrite1 = MPI_Wtime();
      time_write += tWrite1 - tWrite0;
      }
    }
  }
  MPI_Barrier(comm);
  double t_end = MPI_Wtime();
  if (rank == 0){
    timeOut.open("runtimes.txt", ios_base::app);
    timeOut << size << " " << t_end - t_start << "\n";
    timeOut.close();
    timeOut.open("writetimes.txt", ios_base::app);
    timeOut << size << " " << time_write << "\n";
    timeOut.close();
  }
  MPI_Finalize();
  return 0;
}
