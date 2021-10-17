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
	  int results = 0, local_result = 0, flag = 0, i = 0;
	
	  time_t start = time(NULL);
	
	
	  if (world_rank == 0) {
      		fill_random(A, N);
		for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
			MPI_Send(&A[i], 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
			i++;
		 }
       		 int message_received = 1;
		
		  while (results < R && i < N) {
			  MPI_Status status;
			  MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &message_received, &status);
			
			  if (message_received) {    
				  printf("message reveiced");
				  MPI_Recv(&local_result, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
				  MPI_Send(&A[i], 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
				  i++;
				  results += local_result;
			 }
		}
	}
	MPI_Finalize();
	return 0;
}
