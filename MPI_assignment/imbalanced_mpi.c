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
	
	srand(time(NULL) * world_rank);
	
	MPI_Status status;
	MPI_Request request;
	
	int* A; 
	A = allocate_mem(N);
	int results = 0;
	int local_result = 0;
	int flag = 0;
	int step_size = world_size - 1;
	int domain_start = world_rank - 1;
	
	time_t start = time(NULL);
	
	if (world_rank == 0) {
		fill_random(A, N);
		for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
        		MPI_Send(A, N, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(A, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	for (int m = 0; m < N; m += step_size) {
		if (world_rank == 0) {
			if (results >= R) {
				flag = 1;
				for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
        				MPI_Isend(&flag, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, &request);
				}
	
				MPI_Barrier(MPI_COMM_WORLD);
				printf("After barrier%.2f\n", (double)(time(NULL) - start));
				printf("process %d is finished\n",world_rank);
				printf("results %d\n", results);
				MPI_Finalize();
    				return 0;
			}
			for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {
				MPI_Recv(&local_result, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				results += local_result;
				printf("%d\n",results);
			}
    		} else {
			local_result = test(A[domain_start + m]);
			MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			local_result = 0;
			
			MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
			if (flag){
				MPI_Barrier(MPI_COMM_WORLD);
				printf("process %d is finished, at itteration %d\n",world_rank, m);
				MPI_Finalize();
    				return 0;
			}

		}
	}
	for (int k = 1; k <= N % step_size; k++) {
		if (world_rank == 0) {
			if (results >= R) {
				flag = 1;
				for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
        				MPI_Isend(&flag, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, &request);
				}
	
				MPI_Barrier(MPI_COMM_WORLD);
				printf("After barrier%.2f\n", (double)(time(NULL) - start));
				printf("process %d is finished\n",world_rank);
				printf("results %d\n", results);
				MPI_Finalize();
    				return 0;
			}
			for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {
				MPI_Recv(&local_result, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				results += local_result;
				printf("%d\n",results);
			}
    		} else if (world_rank == k){
			local_result = test(A[domain_start + m]);
			MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			local_result = 0;
			
			MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
			if (flag){
				MPI_Barrier(MPI_COMM_WORLD);
				printf("process %d is finished, at itteration %d\n",world_rank, m);
				MPI_Finalize();
    				return 0;
			}

		}
	}
	printf("Did not finish");
	MPI_Finalize();
	return 0;
}