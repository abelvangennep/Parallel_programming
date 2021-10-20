#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

int N = 2500;
int R = 500;

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

void decompose_domain(int N, int world_rank, int world_size, int* subdomain_start, int* subdomain_size) {
  	*subdomain_size = N / (world_size);
  	*subdomain_start = *subdomain_size * (world_rank);
	
  	if (world_rank == world_size - 1) {
    		// Give remainder to last process
    		*subdomain_size += N % world_size;
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
	int results = 0, local_result = 0, flag = 0;
	int subdomain_start, subdomain_size, maximum_sends_recvs;
	
	decompose_domain(N, world_rank, world_size, &subdomain_start, &subdomain_size);
	
	maximum_sends_recvs = N / (world_size) + N % world_size;
	
	if (world_rank == 0) {
		fill_ascending(A, N);
		
		for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {	
        		MPI_Send(A, N, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
		}
	} else {
		MPI_Recv(A, N, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	time_t start = time(NULL);
	for (int m = 0; m < maximum_sends_recvs; m++) {
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
			results += test(A[subdomain_start + m]);
			
			for (int partner_rank = 1; partner_rank < world_size; partner_rank++) {
				MPI_Recv(&local_result, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				results += local_result;
			}
    		} else {
			local_result = test(A[subdomain_start + m]);
			
			MPI_Send(&local_result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			
			local_result = 0;
			
			MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
			
			if (flag) {
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
