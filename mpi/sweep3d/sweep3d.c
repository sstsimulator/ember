
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include <errno.h>

void get_position(const int rank, const int pex, const int pey, int* myX, int* myY) {
	*myX = rank % pex;
        *myY = rank / pex;
}

void compute(long sleep) {
	struct timespec sleepTS;
        sleepTS.tv_sec = 0;
        sleepTS.tv_nsec = sleep;

        struct timespec remainTS;

	if( nanosleep( &sleepTS, &remainTS ) == EINTR ) {
        	while( nanosleep( &remainTS, &remainTS ) == EINTR );
        }
}

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int me = -1;
	int world = -1;

	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &world);

	int pex = -1;
	int pey = -1;
	int nx = 50;
	int ny = 50;
	int nz = 100;
	int kba = 10;
	int repeats = 1;

	int vars = 1;
	long sleep = 1000;

	for( int i = 0; i < argc; ++i ) {
		if( strcmp("-pex", argv[i]) == 0 ) {
			pex = atoi(argv[i+1]);
			i++;
		} else if( strcmp("-pey", argv[i]) == 0) {
			pey = atoi(argv[i+1]);
			i++;
		} else if( strcmp("-iterations", argv[i]) == 0) {
			repeats = atoi(argv[i+1]);
			i++;
		} else if( strcmp("-nx", argv[i]) == 0) {
			nx = atoi(argv[i+1]);
			i++;
		} else if( strcmp("-ny", argv[i]) == 0) {
			ny = atoi(argv[i+1]);
			i++;
		} else if( strcmp("-nz", argv[i]) == 0) {
			nz = atoi(argv[i+1]);
			i++;
		} else if (strcmp("-sleep", argv[i]) == 0) {
			sleep = atol(argv[i+1]);
			i++;
		} else if (strcmp("-vars", argv[i]) == 0) {
			vars = atoi(argv[i+1]);
			i++;
		} else if (strcmp("-kba", argv[i+1]) == 0) {
			kba = atoi(argv[i+1]);
			i++;
		}
	}

	if( kba == 0 ) {
		if( me == 0 ) {
			fprintf(stderr, "K-Blocking Factor must not be zero. Please specify -kba <value > 0>\n");
		}

		exit(-1);
	}

	if( nz % kba != 0 ) {
		if( me == 0 ) {
			fprintf(stderr, "KBA must evenly divide NZ, KBA=%d, NZ=%d, remainder=%d (must be zero)\n",
				kba, nz, (nz % kba));
		}

		exit(-1);
	}

	if( (pex * pey) != world ) {
		if( 0 == me ) {
			fprintf(stderr, "Error: processor decomposition (%d x %d) != number of ranks (%d)\n",
				pex, pey, world);
		}

		exit(-1);
	}

	if( me == 0 ) {
		printf("# Sweep3D Communication Pattern\n");
		printf("# Info:\n");
		printf("# Px:              %8d\n", pex);
		printf("# Py:              %8d\n", pey);
		printf("# Nx x Ny x Nz:    %8d x %8d x %8d\n", nx, ny, nz);
		printf("# KBA:             %8d\n", kba);
		printf("# Variables:       %8d\n", vars);
		printf("# Iterations:      %8d\n", repeats);
	}

	int myX = -1;
	int myY = -1;

	get_position(me, pex, pey, &myX, &myY);

	const int xUp   = (myX != (pex - 1)) ? me + 1 : -1;
	const int xDown = (myX != 0) ? me - 1 : -1;

	const int yUp   = (myY != (pey - 1)) ? me + pex : -1;
	const int yDown = (myY != 0) ? me - pex : -1;
	
	MPI_Status status;

	double* xRecvBuffer = (double*) malloc( sizeof(double) * nx * kba * vars );
	double* xSendBuffer = (double*) malloc( sizeof(double) * nx * kba * vars );

	double* yRecvBuffer = (double*) malloc( sizeof(double) * ny * kba * vars );
	double* ySendBuffer = (double*) malloc( sizeof(double) * ny * kba * vars );

	for( int i = 0; i < nx; ++i ) {
		xRecvBuffer[i] = 0;
		xSendBuffer[i] = i;
	}

	for( int i = 0; i < ny; ++i ) {
		yRecvBuffer[i] = 0;
		ySendBuffer[i] = i;
	}

	struct timeval start;
	struct timeval end;

	gettimeofday( &start, NULL );

	// We repeat this sequence twice because there are really 8 vertices in the 3D data domain and
	// we sweep from each of them, processing the top four first and then the bottom four vertices
	// next. 
	for( int i = 0; i < (repeats * 2); ++i ) {
		// Recreate communication pattern of sweep from (0,0) towards (Px,Py)
		for( int k = 0; k < nz; k += kba ) {
			if( xDown > -1 ) {
				MPI_Recv(xRecvBuffer, (nx * kba * vars), MPI_DOUBLE, xDown, 1000, MPI_COMM_WORLD, &status);
			}

			if( yDown > -1 ) {
				MPI_Recv(yRecvBuffer, (ny * kba * vars), MPI_DOUBLE, yDown, 1000, MPI_COMM_WORLD, &status);
			}

			compute(sleep);

			if( xUp > -1 ) {
				MPI_Send(xSendBuffer, (nx * kba * vars), MPI_DOUBLE, xUp, 1000, MPI_COMM_WORLD);
			}

			if( yUp > -1 ) {
				MPI_Send(ySendBuffer, (nx * kba * vars), MPI_DOUBLE, yUp, 1000, MPI_COMM_WORLD);
			}
		}

		// Recreate communication pattern of sweep from (Px,0) towards (0,Py)
		for( int k = 0; k < nz; k += kba ) {
			if( xUp > -1 ) {
				MPI_Recv(xRecvBuffer, (nx * kba * vars), MPI_DOUBLE, xUp, 2000, MPI_COMM_WORLD, &status);
			}

			if( yDown > -1 ) {
				MPI_Recv(yRecvBuffer, (ny * kba * vars), MPI_DOUBLE, yDown, 2000, MPI_COMM_WORLD, &status);
			}

			compute(sleep);

			if( xDown > -1 ) {
				MPI_Send(xSendBuffer, (nx * kba * vars), MPI_DOUBLE, xDown, 2000, MPI_COMM_WORLD);
			}

			if( yUp > -1 ) {
				MPI_Send(ySendBuffer, (nx * kba * vars), MPI_DOUBLE, yUp, 2000, MPI_COMM_WORLD);
			}
		}

		// Recreate communication pattern of sweep from (Px,Py) towards (0,0)
		for( int k = 0; k < nz; k += kba ) {
			if( xUp > -1 ) {
				MPI_Recv(xRecvBuffer, (nx * kba * vars), MPI_DOUBLE, xUp, 3000, MPI_COMM_WORLD, &status);
			}

			if( yUp > -1 ) {
				MPI_Recv(yRecvBuffer, (ny * kba * vars), MPI_DOUBLE, yUp, 3000, MPI_COMM_WORLD, &status);
			}

			compute(sleep);

			if( xDown > -1 ) {
				MPI_Send(xSendBuffer, (nx * kba * vars), MPI_DOUBLE, xDown, 3000, MPI_COMM_WORLD);
			}

			if( yDown > -1 ) {
				MPI_Send(ySendBuffer, (nx * kba * vars), MPI_DOUBLE, yDown, 3000, MPI_COMM_WORLD);
			}
		}

		// Recreate communication pattern of sweep from (0,Py) towards (Px,0)
		for( int k = 0; k < nz; k += kba ) {
			if( xDown > -1 ) {
				MPI_Recv(xRecvBuffer, (nx * kba * vars), MPI_DOUBLE, xDown, 4000, MPI_COMM_WORLD, &status);
			}

			if( yUp > -1 ) {
				MPI_Recv(yRecvBuffer, (ny * kba * vars), MPI_DOUBLE, yUp, 4000, MPI_COMM_WORLD, &status);
			}

			compute(sleep);

			if( xUp > -1 ) {
				MPI_Send(xSendBuffer, (nx * kba * vars), MPI_DOUBLE, xUp, 4000, MPI_COMM_WORLD);
			}

			if( yDown > -1 ) {
				MPI_Send(ySendBuffer, (nx * kba * vars), MPI_DOUBLE, yDown, 4000, MPI_COMM_WORLD);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday( &end, NULL );

	const double timeTaken = ( ((double) end.tv_sec) + ((double) end.tv_usec) * 1.0e-6 ) -
				( ((double) start.tv_sec) + ((double) start.tv_usec) * 1.0e-6 );
	const double bytesXchng = ((double) repeats) * (
		((double)( xUp > -1 ? sizeof(double) * nx * kba * vars * 2 : 0 )) +
		((double)( xDown > -1 ? sizeof(double) * nx * kba * vars * 2 : 0 )) +
		((double)( yUp > -1 ? sizeof(double) * ny * kba * vars * 2 : 0 )) +
		((double)( yDown > -1 ? sizeof(double) * ny * kba * vars * 2 : 0 )) );

	
	if( (myX == (pex/2)) && (myY == (pey/2)) ) {
		printf("# Results from rank: %d\n", me );
		printf("# %20s %20s %20s\n", "Time", "KBytesXchng/Rank-Max", "MB/S/Rank");
        	printf("  %20.6f %20.4f %20.4f\n",
                        timeTaken, bytesXchng / 1024.0, (bytesXchng / 1024.0) / timeTaken );
	}
	MPI_Finalize();

}
