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
  int J = (M-2)/size + 2;  // Number of rows for each process
  double rankd = rank, Jd = J;
  // Initialize solution
  double un[M][J] = {};
  double unm1[M][J] = {};
  double utemp[M][J] = {};
  double uswap[M] = {}, usenddown[M] = {}, usendup[M] = {};
  ofstream fileOut, initialOut, timeOut;  // Used to save files
  double uOut [M][J];
  if (rank == 0){
    initialOut.open("initial-cond.txt");   
  }

  
  // State start and end points in x direction away from boundary points
  int start = 0, end = J;
  if (rank == 0){start = 1;}  
  if (rank == size-1){end = J-1;}

  // Start time
  MPI_Barrier(comm);
  double t_start = MPI_Wtime();
  double time_write = 0;  // Store time spent writing to files
  
  for (int i=1; i<M-1; ++i){
    for (int j=start; j<end; ++j){
      un[i][j] = exp(-40*(pow(i*dx - 1.4, 2) + pow((rankd*(Jd-2)+j)*dx - 1,2)));
      unm1[i][j] = un[i][j];
    }
  }
  double tWrite0 = MPI_Wtime();
  if (rank > 0){MPI_Send(un, J*M, MPI_DOUBLE, 0, 0, comm);}
  if (rank == 0){
      for (int j=0; j<J; ++j){
	for (int i=0; i<M; ++i){
	  initialOut << un[i][j] << " ";  // Transpose output for simplicity
	}
	initialOut << "\n";
      }
    for (int r = 1; r<size; ++r){
      MPI_Recv(uOut, J*M, MPI_DOUBLE, r, 0, comm, MPI_STATUS_IGNORE);
      for (int j=2; j<J; ++j){
	for (int i=0;i<M; ++i){
	  initialOut << uOut[i][j] << " ";
	}
	initialOut << "\n";
      }
    }
    initialOut.close();
  }
  double tWrite1 = MPI_Wtime();
  time_write += tWrite1 - tWrite0;
  
  // Iterate solution for 2 <= n <= N
  for (double t = 2*dt; t < T; t += dt){
    // perform update
    for (int i=1; i < M-1; ++i){
      for (int j=1; j < J-1; ++j){
	utemp[i][j] = dtdx*(un[i+1][j] + un[i-1][j] + un[i][j+1]
			    + un[i][j-1] - 4*un[i][j]) + 2*un[i][j] - unm1[i][j];
	  }
    }
    // Update parameters
    for (int i=1; i < M-1; ++i){
      for (int j=1; j < J-1; ++j){
	unm1[i][j] = un[i][j];
	un[i][j] = utemp[i][j];
      }
    }

    // halo-swapping

    // Send
    if (rank > 0){
      for (int j=0; j<M; ++j){
	usenddown[j] = un[j][1];
      }
      MPI_Send(usenddown, M, MPI_DOUBLE, rank - 1, 0, comm);
    }
    if (rank < size - 1){
      for (int j=0; j<M; ++j){
	usendup[j] = un[j][J-2];
      }
      MPI_Send(usendup, M, MPI_DOUBLE, rank + 1, 1, comm);
    }

    // Receive
    if (rank > 0){
      MPI_Recv(uswap, M, MPI_DOUBLE, rank - 1, 1, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < M; ++m){un[m][0] = uswap[m];}
    }
    if (rank < size - 1){
      MPI_Recv(uswap, M, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < M; ++m){un[m][J-1] = uswap[m];}
    }

    // Save data if t = 1/3, 2/3, 1
    
    if ((pow(1/1. - t,2) <dt*dt/100)|(pow(1/3. - t,2) <dt*dt/100)|(pow(2/3. - t,2) <dt*dt/100)){
      // Time process of writing
      tWrite0 = MPI_Wtime();
      int tag = t*3;
      if (pow(1/3. - t,2) <dt*dt/100){fileOut.open("data1.txt");}
      if (pow(2/3. - t,2) <dt*dt/100){fileOut.open("data2.txt");}
      if (pow(3/3. - t,2) <dt*dt/100){fileOut.open("data3.txt");}
      if (rank>0){MPI_Send(un, J*M,  MPI_DOUBLE, 0, tag, comm);}
      if (rank == 0){
	for (int j=0; j < J; ++j){
	  for (int i=0; i < M; ++i){
	    fileOut << un[i][j] << " ";
	  }
	  fileOut << "\n";
	}
	for (int r=1; r<size; ++r){
	  // Receive from process r
	  MPI_Recv(uOut, J*M, MPI_DOUBLE, r, tag, comm, MPI_STATUS_IGNORE);
	  // output soln
	  for (int j=2; j < J; ++j){
	    for (int i=0; i<M; ++i){
	      fileOut << uOut[i][j] << " ";
	    }
	    fileOut << "\n";
	  }
	}
	fileOut.close();
      }
      tWrite1 = MPI_Wtime();
      time_write += tWrite1 - tWrite0;
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
