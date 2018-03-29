// Copyright 2009-2018 Sandia Corporation. Under the terms
// of Contract DE-NA0003525 with Sandia Corporation, the U.S.
// Government retains certain rights in this software.
//
// Copyright (c) 2009-2018, Sandia Corporation
// All rights reserved.
//
// Portions are copyright of other developers:
// See the file CONTRIBUTORS.TXT in the top level directory
// the distribution for more information.
//
// This file is part of the SST software package. For license
// information, see the LICENSE file in the top level directory of the
// distribution.

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

	int me = -1;
	int world = -1;

	MPI_Comm_size(MPI_COMM_WORLD, &world);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);

	int repeats = 1;
	int msgsize = 1024;

	for(int i = 1; i < argc; ++i) {
		if( strcmp(argv[i], "-iterations") == 0 ) {
			repeats = atoi(argv[i+1]);
			i++;
		} else if( strcmp(argv[i], "-msgsize") == 0) {
			msgsize = atoi(argv[i+1]);
			i++;
		} else {
			if( 0 == me ) {
				fprintf(stderr, "Unknown option: %s\n", argv[i]);
				exit(-1);
			}
		}
	}

	if( 0 == me ) {
		printf("# Incast Communication Benchmark\n");
		printf("# Info:\n");
		printf("# Incast Ranks:     %8d\n", (world-1));
		printf("# Iterations:       %8d\n", repeats);
		printf("# Message Size:     %8d\n", msgsize);
	}

	double* recvBuffer = NULL;
	double* sendBuffer = NULL;

	if( me == (world-1)) {
		recvBuffer = (double*) malloc( sizeof(double) * msgsize * (world - 1));

		for(int i = 0; i < (msgsize * (world-1)); ++i) {
			recvBuffer[i] = 0.0;
		}
	}

	if( me != (world-1)) {
		sendBuffer = (double*) malloc( sizeof(double) * msgsize );

		for(int i = 0; i < msgsize; ++i) {
			sendBuffer[i] = (double) i;
		}
	}

	MPI_Request* requests = (MPI_Request*) malloc( sizeof(MPI_Request) * (world - 1));
	MPI_Status* status    = (MPI_Status*) malloc( sizeof(MPI_Status) * (world - 1));

	struct timeval start;
	struct timeval end;

	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday( &start, NULL );

	for(int i = 0; i < repeats; ++i) {
		if( me == (world - 1) ) {
			for( int r = 0; r < world - 1; ++r) {
				MPI_Irecv(&recvBuffer[r * msgsize], msgsize, MPI_DOUBLE, r, 1000,
					MPI_COMM_WORLD, &requests[r]);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if( me != (world - 1) ) {
			MPI_Send(sendBuffer, msgsize, MPI_DOUBLE, (world - 1), 1000, MPI_COMM_WORLD);
		} else {
			MPI_Waitall( (world-1), requests, status );
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	gettimeofday( &end, NULL );

	if( recvBuffer != NULL )
		free(recvBuffer);

	if( sendBuffer != NULL )
		free(sendBuffer);

	if( (world - 1) == me ) {
		const double timeTaken = (((double) end.tv_sec) + (double) end.tv_usec * 1.0e-6) -
			(((double) start.tv_sec) + (double) start.tv_usec * 1.0e-6);
		const double msgsRecv  = ((double) (repeats * (world-1)));
		const double dataRecv  = (((double) repeats) * ((double) msgsize) * ((double) (world-1))) /
			(1024.0 * 1024.0);

		printf("# Statistics:\n");
		printf("# %20s %20s %20s\n", "Time Taken", "Msgs/s", "MB/s");
		printf("  %20.6f %20f %20.6f\n", timeTaken, msgsRecv/timeTaken, dataRecv / timeTaken);
		
	}

	MPI_Finalize();

}
