#include "mpi.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 

using namespace std;

int Nx = 20;			
int Ny = 10;			// grid
int Nt = 10;
double hx = 1.0 / Nx;
double hy = 0.5 / Ny;	// step sizes
double ht = 10 / Nt;

FILE* logFile;			// log file reference

double** allocateArray(int rows, int cols) {
	/* manually allocate a two-dimensional array that can be sent via MPI_Send */
	
	double* data = (double*) malloc(rows * cols * sizeof(double));
	double** array = (double**) malloc(rows * sizeof(double));
	
	for (int i = 0; i < rows; i++) {
		array[i] = &(data[cols * i]);
	}
	
	return array;
}

void assembleResult(double** T, int _nx, int _ny, int rank, int size) {
	/* send this rank's result of the current cycle to rank 0 for the final result assembly */
	
	MPI_Status status[size];
	
	MPI_Send (&(T[0][0]), _nx * _ny, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
}

void exchangeBorders(double* leftBorder, double* rightBorder, double* leftGhost, double* rightGhost, int rank, int size) {
	/* exchange borders with neighbors */
	
	MPI_Status status;
	
	if (rank % 2 == 0) {
		/*
		* even ranks send their borders first, odd ranks send their borders after receiving ghosts
		* ranks that have only one neighbor skip the extra send/receive
		*/
		if (rank < size - 1) {
			MPI_Send(rightBorder, Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
			MPI_Recv(rightGhost, Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
		}
		if (rank > 1) {
			MPI_Send(leftBorder, Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
			MPI_Recv(leftGhost, Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
		}
	} else {
		if (rank > 1) {
			MPI_Recv(leftGhost, Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(leftBorder, Ny, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
		}
		if (rank < size - 1) {
			MPI_Recv(rightGhost, Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
			MPI_Send(rightBorder, Ny, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
		}
	}
}

void overrideArray(double** dest, double** source, int _nx, int _ny) {
	/* override dest array with data from source array */
	
	for (int i = 0; i < _nx; i++) {
		for (int j = 0; j < _ny; j++) {
			dest[i][j] = source[i][j];
		}
	}
}

void calculate(double** T, int& rank, int& size) {
	
	/* initialization section	*/
	
	int _nx = Nx / (size - 1);
	int _ny = Ny;							// segment of the grid assigned to the current process
	
	int _shift = _nx * (rank - 1);
	/*
	int _left = _nx * (rank - 1);
	int _right = _nx * rank;
	*/
	double** T0 = allocateArray(_nx, _ny);
	double** T1 = allocateArray(_nx, _ny);
	
	double lambda[_nx][_ny];
	
	for (int i = 0; i < _nx; i++) {
		for (int j = 0; j < _ny; j++) {
			T0[i][j] = 300;
			T1[i][j] = 300;
			lambda[i][j] = (i * hx >= 0.25 && i * hx <= 0.65) && (j * hy >= 0.1 && j * hy <=0.25) ? 10e-2 : 10e-4;
		}
	}
	assembleResult(T1, _nx, _ny, rank, size);
	
	double leftBorder[_ny];					// grid borders
	double rightBorder[_ny];
	
	double leftGhost[_ny];					// neighboring grid borders
	double rightGhost[_ny];
	
	/* calculation loop */
	
	double alpha;
	double beta;
	double prevAlpha;
	double prevBeta;
	
	for (int t = 1; t < Nt; t++) {
		/* time loop */
		
		for (int j = 0; j < _ny; j++) {
			leftBorder[j] = T0[0][j];
			rightBorder[j] = T0[_nx - 1][j];
		}
		exchangeBorders(leftBorder, rightBorder, leftGhost, rightGhost, rank, size);
		
		for (int i = 0; i < _nx; i++) {
			for (int j = 0; j < _ny; j++) {
				/* grid loop */
				
				if (j == 0) {
					T1[i][j] = T[i + _shift][j];
					continue;
				}
				if (j == Ny - 1) {
					T1[i][j] = T[i + _shift][j];
					continue;
				}
				
				double tau = t * ht;
				double lambdaPlus;
				double lambdaMinus;
				double F;
				
				if (i == 0) {
					if (rank == 1) {
						T1[i][j] = T[i][j];
						continue;
					}
					
					lambdaPlus = lambda[i][j];
					lambdaMinus = lambda[i][j];
					F = T0[i][j] / tau + ( lambdaPlus * (T0[i + 1][j] - T0[i][j]) 
										- lambdaMinus * (T0[i][j] - leftGhost[j]) ) / (2 * pow(hx, 2));
				} else if (i == _nx - 1) {
					if (rank == size - 1) {
						T1[i][j] = T[Nx - 1][j];
						continue;
					}
					
					lambdaPlus = lambda[i][j];
					lambdaMinus = lambda[i][j];
					F = T0[i][j] / tau + ( lambdaPlus * (rightGhost[j] - T0[i][j])
											- lambdaMinus * (T0[i][j] - T0[i - 1][j]) ) / (2 * pow(hx, 2));
				} else {
					lambdaPlus = (lambda[i + 1][j] + lambda[i][j]) / 2;
					lambdaMinus = (lambda[i - 1][j] + lambda[i][j]) / 2;
					F = T0[i][j] / tau + ( lambdaPlus * (T0[i + 1][j] - T0[i][j]) 
											- lambdaMinus * (T0[i][j] - T0[i - 1][j]) ) / (2 * pow(hx, 2));
				}
				
				double A = - ( lambdaMinus / (2 * pow(hy, 2)) );
				double B = - ( lambdaPlus / (2 * pow(hy, 2)) );
				double C = 1 / tau - A - B;
				
				alpha = -B / (C + A * prevAlpha);
				beta = (F - A * prevBeta) / (C + A * prevAlpha);
				
				T1[i][j] = alpha * T0[i][j] + beta;
			}
		}
		
		assembleResult(T1, _nx, _ny, rank, size);
		
		overrideArray(T0, T1, _nx, _ny);

		prevAlpha = alpha;
		prevBeta = beta;
	}
}

void writeToLog(double** T, int cycle) {
	/* write the result of the current file to log */
	
	fprintf(logFile, "----- cycle %d -----\n", cycle);
	fprintf(logFile, "t = %f\n", cycle * ht);
	
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			fprintf(logFile, " %9.3f ", T[i][j]);
		}
		fprintf(logFile, "\n");
	}
	
	fprintf(logFile, "\n\n");
}

int main(int argc, char* argv[]) {
	
	logFile = fopen("log.txt", "w+");

	MPI_Init(&argc, &argv);

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status status;
	
	if (rank == 0) {
        double** T = new double*[Nx];
	    for (int i = 0; i < Nx; i++) {
            T[i] = new double[Nx];
        }
		double** data = allocateArray(Nx / (size - 1), Ny);

		for (int t = 0; t < Nt; t++) {
			for (int R = 1; R < size; R++) {
				MPI_Recv(&(data[0][0]), (Nx / (size - 1)) * Ny, MPI_DOUBLE, R, 1, MPI_COMM_WORLD, &status);
				
				int shift = (Nx / (size - 1)) * (R - 1);
				
				for (int i = 0; i < Nx / (size - 1); i++) {
					for (int j = 0; j < Ny; j++) {
						T[i + shift][j] = data[i][j];
					}
				}
			}
			
			writeToLog(T, t);
		}
	} else {
        double** T = new double*[Nx];
	    for (int i = 0; i < Nx; i++) {
            T[i] = new double[Nx];
        }

		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				
				if (i == 0) {
					T[i][j] = 600;
					continue;
				}
				
				if (i == Nx - 1) {
					T[i][j] = 1200;
					continue;
				}
				
				if (j == 0) {
					T[i][j] = 600 * (1 + i * hx);
					continue;
				}
				
				if (j == Ny - 1) {
					T[i][j] = 600 * (1 + pow(i * hx, 3));
					continue;
				}
				
				T[i][j] = 0;
			}
		}
		
		calculate(T, rank, size);
	}
	
	MPI_Finalize();
	fclose(logFile);
	return 0;
}
