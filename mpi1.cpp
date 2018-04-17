#include "mpi.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <string>

using namespace std;

int Nx = 20;			
int Ny = 10;			// grid
int Nt = 25;
double hx = 1.0 / (double) Nx;
double hy = 0.5 / (double) Ny;	// step sizes
double ht = 10.0 / (double) Nt;

FILE* logFile;			// log file reference

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

void writeToSeparateFile(double** T, int cycle) {
	/* write the result of the current file to log */
	
	FILE* file;
	
	char buff[20];
	snprintf(buff, sizeof(buff), "cycle%d.txt", cycle);
	string name = buff;
	
	file = fopen(name.c_str(), "w+");
	
	for (int j = 0; j < Ny; j++) {
		for (int i = 0; i < Nx; i++) {
			fprintf(file, " %9.3f ", T[i][j]);
		}
		fprintf(file, "\n");
	}
	
	fprintf(file, "\n\n");
	
	fclose(file);
}

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

void calculate(int& rank, int& size) {
	
	//см. ПУНКТ 2
	
	int _nx = Nx / (size - 1);
	int _ny = Ny;							// segment of the grid assigned to the current process
	
	int _shift = _nx * (rank - 1);
	
	double** T0 = allocateArray(_nx, _ny);
	double** T1 = allocateArray(_nx, _ny);
	double lambda = 10e-4;
	
	for (int i = 0; i < _nx; i++) {
		for (int j = 0; j < _ny; j++) {
			
			if (i == 0 && rank == 1) {
				T0[i][j] = (double) 600.0;
				T1[i][j] = (double) 600.0;
				continue;
			}
			if (i == _nx - 1 && rank == size - 1) {
				T0[i][j] = (double) 1200.0;
				T1[i][j] = (double) 1200.0;
				continue;
			}
			
			if (j == 0) {
				T0[i][j] = (double)(600.0 * (1 + (i + _shift + 1) * hx));
				T1[i][j] = (double)(600.0 * (1 + (i + _shift + 1) * hx));
				continue;
			}
			if (j == Ny - 1) {
				T0[i][j] = (double)(600.0 * (1 + pow((i + _shift + 1) * hx, 3)));
				T1[i][j] = (double)(600.0 * (1 + pow((i + _shift + 1) * hx, 3)));
				continue;
			}
			T0[i][j] = (double) 300.0;
			T1[i][j] = (double) 300.0;
		}
	}
	
	assembleResult(T1, _nx, _ny, rank, size);

	double leftBorder[_ny];					// grid borders
	double rightBorder[_ny];
	
	double leftGhost[_ny];					// neighboring grid borders
	double rightGhost[_ny];
	//см. ПУНКТ 3
	for (int t = 1; t < Nt; t++) {
		for (int j = 0; j < _ny; j++) {
			leftBorder[j] = T0[0][j];
			rightBorder[j] = T0[_nx - 1][j];
		}
		
		exchangeBorders(leftBorder, rightBorder, leftGhost, rightGhost, rank, size);
		
		for (int i = 0; i < _nx; i++) {
			for (int j = 0; j < _ny; j++) {
				if (i == 0 && rank == 1) {
					T1[i][j] = T0[i][j];
					continue;
				}
				if (i == _nx - 1 && rank == size - 1) {
					T1[i][j] = T0[i][j];
					continue;
				}
				if (j == 0) {
					T1[i][j] = T0[i][j];
					continue;
				}
				if (j == _ny - 1) {
					T1[i][j] = T0[i][j];
					continue;
				}
				
				double leftT0;
				double rightT0;
				
				if (i == 0) {
					leftT0 = leftGhost[j];
				} else {
					leftT0 = T0[i - 1][j];
				}
				
				if (i == _nx - 1) {
					rightT0 = rightGhost[j];
				} else {
					rightT0 = T0[i + 1][j];
				}
				
				T1[i][j] = T0[i][j] + lambda * ((double)ht / (double)pow(hx, 2.0)) * (leftT0 + T0[i][j - 1] + rightT0 + T0[i][j + 1] - 4.0 * T0[i][j]);
			}
		}
		
		assembleResult(T1, _nx, _ny, rank, size);
		
		overrideArray(T0, T1, _nx, _ny);
	}
}

int main(int argc, char* argv[]) {
	
	logFile = fopen("log.txt", "w+");

	MPI_Init(&argc, &argv);

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status status;
	
	// см. ПУНКТ 1
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
			
			if (t % 5 == 0) writeToSeparateFile(T, t);
		}
	} else {
		calculate(rank, size);
	}
	
	MPI_Finalize();
	fclose(logFile);
	return 0;
}
