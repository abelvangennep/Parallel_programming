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
	
	MPI_Status status;
	MPI_Request request;
	
	
	int results = 0, local_result = 0, i = 0, message_received = 0, a_i;
	int* A; 
	A = allocate_mem(N);
	
	if (world_rank == 0) {
		fill_random(A, N);
		
		for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
        		MPI_Send(A, N, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(A, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	
	if (world_rank == 0) {
		time_t start = time(NULL);
		for (int l = 0; l < 5; l++) {
			results += test(A[l]);
		}
		int m = world_size * 5;
			
		for (m; m < N; m += 5) {
			MPI_Status status;

			MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

			MPI_Recv(&local_result, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
			results += local_result;
			if ( results >= R) {
				for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
					int flag = -1;
					MPI_Send(&flag, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				printf("After barrier%.2f\n", (double)(time(NULL) - start));
				printf("at itteration:%d\n",m);
				MPI_Finalize();
				return 0;
			}
			MPI_Send(&m, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);			
		}
	} else {
		for (int l; l < 5; l++) {
			local_result += test(A[l+world_rank*5]);
		}
		MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		local_result = 0;
		for (int k; k < N; k++) {
			MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			
			if (a_i != -1){
				local_result += test(A[i]);
				MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				local_result = 0;
			} else {
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
				return 0;
			}
		}
	}
	MPI_Finalize();
	return 0;
}
