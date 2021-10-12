EXECS=test_mpi
MPICC?=mpicc

all: ${EXECS}

test_mpi:test_mpi.c
	${MPICC} -Wall -O3 -o test_mpi test_mpi.c -lm

clean:
	rm -f *.o ${EXECS} *~ *core
