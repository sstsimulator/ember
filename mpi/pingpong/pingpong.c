
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);

	int me = -1;
	int world = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &world);

	if( 0 == me ) {
		printf("MPI PingPong Pattern\n");
		printf("Info:\n");
		printf("- Total Ranks:     %8d\n", world);
	}

	MPI_Finalize();
}
