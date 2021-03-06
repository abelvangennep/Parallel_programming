#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

int N = 500;
int R = 100;

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
	
	if (world_size != 1) {
    		fprintf(stderr, "World size must be one for %s\n", argv[0]);
    		MPI_Abort(MPI_COMM_WORLD, 1);
  	}
	
	int* A; 
	A = allocate_mem(N);
	int results = 0;
	
	fill_random(A, N);
	
	time_t start = time(NULL);
	
	for (int m = 0; m < N; m++) {
		if (test_imbalanced(A[m])){
			results++;
		}
// 		results += test(A[m]);
		
		if (results >= R){
			printf("%.2f\n", (double)(time(NULL) - start));
			printf("process is finished at itteration %d\n",m);
			MPI_Finalize();
    			return 0;
		}
	}
	printf("%.2f\n", (double)(time(NULL) - start));
	return 0;
}
