#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

int N = 50;
int R = 10;

int test(int x) {
	// Transform to a number beween 0 and 2 * pi.
	double xd = ((double)x * (2 * M_PI) / N);

	// Do a fixed amount of computation.
	int n = 25000000;
	for (int i = 0; i < n; i++) {
		xd = sin(xd);
		xd = asin(xd);
	}
	xd = sin(xd);

	// Test the result.
	return xd < -0.5;
}

int test_imbalanced(int x) {
	// Transform to a number beween 0 and 2 * pi.
	double xd = ((double)x * (2 * M_PI) / N);

	// Do computation, the amount depends on the value of x.
	int n = 100000 * x;

	for (int i = 0; i < n; i++) {
		xd = sin(xd);
		xd = asin(xd);
	}
	xd = sin(xd);

	// Test the result.
	return xd < -0.5;
}

int *allocate_mem(int N) {
	int *A = calloc(N, sizeof(int));
	if (!A)
		exit(0);
	return A;
}

void fill_random(int *A, int N) {
	srand(time(0));

	for (int i = 0; i < N; i++) {
		A[i] = (rand() % N);
	}
}

void fill_ascending(int *A, int N) {
	for (int i = 0; i < N; i++) {
		A[i] = i;
	}
}

int main(int argc, char *argv[]) {
	MPI_Init(NULL, NULL);
  	int world_size;
  	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	MPI_Status status;
	MPI_Request request;
	
	int* A; 
	A = allocate_mem(N);
	int results = 0, local_result = 0, flag = 1, i = 0, a_i;
	
	time_t start = time(NULL);
	
	printf("1");
	if (world_rank == 0) {
      		fill_random(A, N);
		for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
			MPI_Send(&A[i], 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
			i++;
		}
	} else {
		while (flag){
			printf("start process");
			MPI_Recv(&a_i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("Received a = %d",a_i);
			if (a_i != -1){
				local_result = test(a_i);
				MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			} else {
				flag = 0;
			}
		}
	}
	MPI_Finalize();
	return 0;
}
