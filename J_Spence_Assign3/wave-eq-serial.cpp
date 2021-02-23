#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main () {
  // Declare integration parameters
  int M = 2306; // Number spacial points
  double Md = M;
  double dx = 2/(Md-2), dt = 0.2/(Md-2);
  double dtdx = dt*dt/(dx*dx);  // Used in calculations below
  double T = 1;
  ofstream fileOut, initialOut;  // Used to save files
  fileOut.open("data.txt");
  initialOut.open("initial-cond.txt");

  // Initialize solution
  double un[M][M] = {};
  double unm1[M][M] = {};
  double utemp[M][M] = {};
  for (int i=1; i<M-1; ++i){
    for (int j=1; j<M-1; ++j){
      un[i][j] = exp(-40*(pow(i*dx - 1.4, 2) + pow(j*dx - 1,2)));
      unm1[i][j] = un[i][j];
    }
  }
  for (int i=0; i<M; ++i){
    for (int j=0; j<M; ++j){
      initialOut << un[i][j] << " ";
    }
    initialOut << "\n";
  }
  initialOut.close();
  // Iterate solution for 2 <= n <= N
  for (double t = 2*dt; t < T; t += dt){
    // perform update
    for (int i=1; i < M-1; ++i){
      for (int j=1; j < M-1; ++j){
	utemp[i][j] = dtdx*(un[i+1][j] + un[i-1][j] + un[i][j+1]
			    + un[i][j-1] - 4*un[i][j]) + 2*un[i][j] - unm1[i][j];
	  }
    }
    // Update parameters
    for (int i=1; i < M-1; ++i){
      for (int j=1; j < M-1; ++j){
	unm1[i][j] = un[i][j];
	un[i][j] = utemp[i][j];
      }
    }

    // Save data if t = 1/3, 2/3, 1
    if ((pow(1/3. - t,2) <dt*dt/100) | (pow(2./3.- t, 2) < dt*dt/100) | (pow(1. - t,2) < dt*dt/100)){
      cout << "TEST";
      for (int i=0; i < M; ++i){
        for (int j=0; j < M; ++j){
          fileOut << un[i][j] << " ";
        }
        fileOut << "\n";
      }
    }
  }
  fileOut.close();
  return 0;
}
