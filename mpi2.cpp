#include "mpi.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <string>

#define ROOT 0
#define SIZE 4
using namespace std;

double** allocateArray(int rows, int cols) {
	/* manually allocate a two-dimensional array that can be sent via MPI_Send */
	
	double* data = (double*) malloc(rows * cols * sizeof(double));
	double** array = (double**) malloc(rows * sizeof(double));
	
	for (int i = 0; i < rows; i++) {
		array[i] = &(data[cols * i]);
	}
	
	return array;
}

int main(int argc, char* argv[]) {
	
	double** A;
	double b[SIZE] = { 2.0, 3.0, 4.0, 5.0 };
	double c[SIZE];
	
	MPI_Init(&argc, &argv);
	
	int mpiSize;
	int mpiRank;
	MPI_Comm mpiComm;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	
	if (SIZE % mpiSize > 0) {
		printf("Please specify the number of processes as a divisor of %d!\n", SIZE);
		MPI_Finalize();
        exit(EXIT_FAILURE);
	}
	
	if (mpiRank == ROOT) {
		A = allocateArray(SIZE, SIZE);
		for (int i = 0; i < SIZE; i++) {
			for (int j = 0; j < SIZE; j++) {
				A[i][j] = (i == j) ? 1.0 : 0.0;
			}
		}
	}
	
	int rowsInSubmatrix = SIZE / mpiSize;
	
	double resultBuffer[rowsInSubmatrix];
	double** submatrix = allocateArray(rowsInSubmatrix, SIZE);
	
	MPI_Scatter(&A[0][0], rowsInSubmatrix * SIZE, MPI_FLOAT, &submatrix[0][0], rowsInSubmatrix * SIZE, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
	
	if (mpiRank == 0) {
		for (int i = 0; i < rowsInSubmatrix; i++) {
			for (int j = 0; j < SIZE; j++) {
				printf("%f  ", submatrix[i][j]);
			}
			printf("\n");
		}
	}
	
	//printf("i'm rank %d and here's a number %7.1f\n", mpiRank, submatrix[0][0]);
	
	for (int i = 0; i < rowsInSubmatrix; i++) {
		double temp = 0.0;
		for (int j = 0; j < SIZE; j++) {
			
			
			//temp += b[j] + submatrix[i][j];
		}
		//resultBuffer[i] = temp;
	}
	
	/*MPI_Gather(resultBuffer, rowsInSubmatrix, MPI_FLOAT, c, rowsInSubmatrix, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
	
	if (mpiRank == ROOT) {
		for (int i = 0; i < SIZE; i++) {
			printf("%f\n", c[i]);
		}
	}
	*/
	MPI_Finalize();
	return 0;
}
