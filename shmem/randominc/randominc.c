
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>

#include <shmem.h>

#ifdef USE_SNL_RANDOM
#include "xorrand.h"
#endif

#define GUPS_DEFAULT_DATA_ITEMS 32 * 1024 * 1024
#define GUPS_DEFAULT_UPDATES    4096
#define GUPS_DEFAULT_ITERATIONS 128

int main(int argc, char* argv[]) {

	// Start up the SHMEM world...
	start_pes(0);

	const int me = shmem_my_pe();
	const int pecount = shmem_n_pes();

	int data_size  = GUPS_DEFAULT_DATA_ITEMS;
	int updates    = GUPS_DEFAULT_UPDATES;
	int iterations = GUPS_DEFAULT_ITERATIONS;

	if( argc > 1 ) {
		data_size = atoi(argv[2]);
	}

	if( argc > 2 ) {
		updates = atoi(argv[3]);
	}

	if( 0 == me ) {
		printf("Allocating shared memory regions...\n");
		printf("Size per PE is: %d values (size of value = %d bytes)\n", data_size,
			(int)(sizeof(long long int)));
		printf("Iterations: %d\n", iterations);
	}

	long long int* data = (long long int*) shmalloc( sizeof(long long int) * data_size );

	if( NULL == data ) {
		fprintf(stderr, "Unable to allocate data array on PE: %d\n", me);
		exit(-1);
	}

	for( int i = 0; i < data_size; ++i ) {
		data[i] = 0;
	}

	//////////////////////////////////////////////////////////////////

	shmem_barrier_all();

	struct timeval start, end;
	gettimeofday( &start, NULL );

	if( 0 == me ) {
		printf("Performing %d updates per PE across the domain...\n", updates);
	}

	#ifdef USE_SNL_RANDOM
		snl_random randNum;
		seed_random(&randNum, (unsigned int) (start.tv_usec));
	#else
		srand( (unsigned int) (start.tv_usec) );
	#endif

	for( int j = 0; j < iterations; ++j) {
		for( int i = 0; i < updates; ++i ) {
			#ifdef USE_SNL_RANDOM
				int pe   = (int) (next_rand_u32(&randNum) % pecount);
			#else
				int pe   = rand() % pecount;
			#endif

			while( pe == me ) {
				#ifdef USE_SNL_RANDOM
					pe   = (int) (next_rand_u32(&randNum) % pecount);
				#else
					pe   = rand() % pecount;
				#endif
			}

			#ifdef USE_SNL_RANDOM
				const int addr = (int) (next_rand_u32(&randNum) % data_size);
			#else
				const int addr = rand() % data_size;
			#endif

			shmem_longlong_inc( &data[addr], pe );
		}
	}

	// Barrier to make sure all of the updates were completed
	shmem_barrier_all();

	gettimeofday( &end, NULL );

	if( 0 == me ) {
		double startSeconds = start.tv_sec + ( (double) start.tv_usec ) * 1.0e-6;
		double endSeconds = end.tv_sec + ( (double) end.tv_usec ) * 1.0e-6;
		double updateTotal = (double) iterations * (double) updates * (double) pecount;

		printf("GUpdates    = %25.6f\n", (updateTotal / (1000.0 * 1000.0 * 1000.0) ));
		printf("Time        = %25.6f\n", (endSeconds - startSeconds));
		printf("GUP/s       = %25.6f\n", ((updateTotal / (1000.0 * 1000.0 * 1000.0 )) /
			(endSeconds - startSeconds)));
	}

	//////////////////////////////////////////////////////////////////

	long long int mytotal = 0;

	for( int i = 0; i < data_size; ++i ) {
		mytotal += data[i];
	}

	for( int i = 0; i < pecount; ++i ) {
		shmem_barrier_all();

		if( me == i ) {
			printf("PE: %10d total is: %lld\n", i, mytotal);
		}

		shmem_barrier_all();
	}

	//////////////////////////////////////////////////////////////////

	shmem_barrier_all();

	if( 0 == me ) {
		printf("Deallocating shared memory regions.\n");
	}

	shfree( data );

	shmem_barrier_all();
	return 0;
}
