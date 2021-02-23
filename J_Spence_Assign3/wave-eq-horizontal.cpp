#include <iostream>
#include <cmath>
#include <fstream>
#include <mpi.h>

using namespace std;

int main () {
  // Declare integration parameters
  int M = 2306; // Number spacial points
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
  double** un = new double* [J];
  double** unm1 = new double* [J];
  double** utemp = new double* [J];
  double uswap[M] = {};
  ofstream fileOut, initialOut, timeOut;  // Used to save files
  double** uOut  = new double* [J];
  for (int j=0; j < J; ++j){
    un[j] = new double [M];
    unm1[j] = new double [M];
    utemp[j] = new double [M];
    uOut[j] = new double [M];
  }
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
  
  for (int i=start; i<end; ++i){
    for (int j=1; j<M-1; ++j){
      un[i][j] = exp(-40*(pow((rankd*(Jd-2)+i)*dx - 1.4, 2) + pow(j*dx - 1,2)));
      unm1[i][j] = un[i][j];
    }
  }
  double tWrite0 = MPI_Wtime();
  
  if (rank > 0){
    for (int j = 0; j<J; ++j){MPI_Send(&un[j][0], M, MPI_DOUBLE, 0, j, comm);}
  }
  if (rank == 0){
      for (int i=0; i<J; ++i){
	for (int j=0; j<M; ++j){
	  initialOut << un[i][j] << " ";
	}
	initialOut << "\n";
      }
    for (int r = 1; r<size; ++r){
      for (int j=0; j < J; ++j){
	MPI_Recv(&uOut[j][0], M, MPI_DOUBLE, r, j, comm, MPI_STATUS_IGNORE);
      }
      for (int i=2; i<J; ++i){
	for (int j=0; j<M; ++j){
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
  int count = 0;
  for (double t = 2*dt; t < T; t += dt){
    ++count;
    if (rank == 0 && count%100 == 0){cout << "\n t = " << t;}
    // perform update
    for (int i=1; i < J-1; ++i){
      for (int j=1; j < M-1; ++j){
	utemp[i][j] = dtdx*(un[i+1][j] + un[i-1][j] + un[i][j+1]
			    + un[i][j-1] - 4*un[i][j]) + 2*un[i][j] - unm1[i][j];
	  }
    }
    // Update parameters
    for (int i=1; i < J-1; ++i){
      for (int j=1; j < M-1; ++j){
	unm1[i][j] = un[i][j];
	un[i][j] = utemp[i][j];
      }
    }

    // halo-swapping

    // Send
    if (rank > 0){
      MPI_Send(&un[1][0], M, MPI_DOUBLE, rank - 1, 0, comm);
    }
    if (rank < size - 1){
      MPI_Send(&un[J-2][0], M, MPI_DOUBLE, rank + 1, 1, comm);
    }

    // Receive
    if (rank > 0){
      MPI_Recv(&uswap, M, MPI_DOUBLE, rank - 1, 1, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < M; ++m){un[0][m] = uswap[m];}
    }
    if (rank < size - 1){
      MPI_Recv(&uswap, M, MPI_DOUBLE, rank + 1, 0, comm, MPI_STATUS_IGNORE);
      for (int m=0; m < M; ++m){un[J-1][m] = uswap[m];}
    }

    // Save data if t = 1/3, 2/3, 1
    
    if ((pow(1/1. - t,2) <dt*dt/100)|(pow(1/3. - t,2) <dt*dt/100)|(pow(2/3. - t,2) <dt*dt/100)){
      // Time process of writing
      tWrite0 = MPI_Wtime();
      int tag = t*3;
      if (pow(1/3. - t,2) <dt*dt/100){fileOut.open("data1.txt");}
      if (pow(2/3. - t,2) <dt*dt/100){fileOut.open("data2.txt");}
      if (pow(3/3. - t,2) <dt*dt/100){fileOut.open("data3.txt");}
      if (rank > 0){
	for (int j = 0; j<J; ++j){MPI_Send(&un[j][0], M, MPI_DOUBLE, 0, J*tag + j, comm);}
      }
      if (rank == 0){
	for (int i=0; i < J; ++i){
	  for (int j=0; j < M; ++j){
	    fileOut << un[i][j] << " ";
	  }
	  fileOut << "\n";
	}
	for (int r=1; r<size; ++r){
	  // Receive from process r
	  for (int j = 0; j<J; ++j) {
	    MPI_Recv(&uOut[j][0], M, MPI_DOUBLE, r, J*tag + j, comm, MPI_STATUS_IGNORE);
	  }
	  // output soln
	  for (int i=2; i < J; ++i){
	    for (int j=0; j<M; ++j){
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
    timeOut.open("runtimes-hor.txt", ios_base::app);
    timeOut << size << " " << t_end - t_start << "\n";
    timeOut.close();
    timeOut.open("writetimes.txt", ios_base::app);
    timeOut << size << " " << time_write << "\n";
    timeOut.close();
  }
  for (int j=0; j < J; ++j){
      delete[] un[j];
      delete[] unm1[j];
      delete[] utemp[j];
      delete[] uOut[j];
    }
  delete[] un;
  delete[] unm1;
  delete[] utemp;
  delete[] uOut;
  MPI_Finalize();
  return 0;
}
