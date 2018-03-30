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

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <shmem.h>

#ifdef USE_SNL_RANDOM
#include "xorrand.h"
#endif

#define GUPS_DEFAULT_DATA_ITEMS 32 * 1024 * 1024
#define GUPS_DEFAULT_UPDATES 4096
#define GUPS_DEFAULT_ITERATIONS 32
#define HOTSPOT_PROB_MULTIPIER 4

#ifdef USE_SNL_RANDOM
int generate_next_pe(const int pecount, const int me, snl_random *randNum) {
#else
int generate_next_pe(const int pecount, const int me) {
#endif
  int pe = -1;

  // If I'm the last PE, then I need to make sure
  // I'm targetting the other PEs first.
  if (me == pecount - 1) {
#ifdef USE_SNL_RANDOM
    pe = (next_rand_u32(randNum) % (pecount - 1));
#else
    pe = rand() % (pecount - 1);
#endif
  } else {
    int pecountHS = pecount + HOTSPOT_PROB_MULTIPIER;
#ifdef USE_SNL_RANDOM
    pe = (next_rand_u32(randNum) % pecountHS);
#else
    pe = rand() % pecountHS;
#endif

    // If we generate a PE higher than we have
    // clamp ourselves to the highest PE
    if (pe >= pecount) {
      pe = pecount - 1;
    }
  }

  return pe;
}

int main(int argc, char *argv[]) {

  // Start up the SHMEM world...
  start_pes(0);

  const int me = shmem_my_pe();
  const int pecount = shmem_n_pes();

  int data_size = GUPS_DEFAULT_DATA_ITEMS;
  int updates = GUPS_DEFAULT_UPDATES;
  int iterations = GUPS_DEFAULT_ITERATIONS;

  if (argc > 1) {
    data_size = atoi(argv[2]);
  }

  if (argc > 2) {
    updates = atoi(argv[3]);
  }

  if (0 == me) {
    printf("Allocating shared memory regions...\n");
    printf("= Benchmark Summary: "
           "==============================================\n");
    printf("Data Region Size:      %10d items\n", data_size);
    printf("Update Size:           %10d bytes per item\n",
           (int)(sizeof(long long int)));
    printf("Updates:               %10d\n", updates);
    printf("Iterations:            %10d\n", iterations);
    printf("Hotspot Multiplier:    %10d\n", HOTSPOT_PROB_MULTIPIER);
    printf("==================================================================="
           "\n");
    fflush(stdout);
  }

  long long int *data =
      (long long int *)shmalloc(sizeof(long long int) * data_size);

  if (NULL == data) {
    fprintf(stderr, "Unable to allocate data array on PE: %d\n", me);
    exit(-1);
  }

  for (int i = 0; i < data_size; ++i) {
    data[i] = 0;
  }

  //////////////////////////////////////////////////////////////////

  shmem_barrier_all();

  struct timeval start, end;
  gettimeofday(&start, NULL);

  if (0 == me) {
    printf("Performing %d updates per PE across the domain...\n", updates);
  }

#ifdef USE_SNL_RANDOM
  snl_random randNum;
  seed_random(&randNum, (unsigned int)(start.tv_usec));
#else
  srand((unsigned int)(start.tv_usec));
#endif

  for (int j = 0; j < iterations; ++j) {
    for (int i = 0; i < updates; ++i) {
      int pe = me;

      while (pe == me) {
#ifdef USE_SNL_RANDOM
        pe = generate_next_pe(pecount, me, &randNum);
#else
        pe = generate_next_pe(pecount, me);
#endif
      }

#ifdef USE_SNL_RANDOM
      const int addr = (int)(next_rand_u32(&randNum) % data_size);
#else
      const int addr = rand() % data_size;
#endif

      fflush(stdout);
      shmem_longlong_inc(&data[addr], pe);
    }
  }

  // Barrier to make sure all of the updates were completed
  shmem_barrier_all();

  gettimeofday(&end, NULL);

  if (0 == me) {
    double startSeconds = start.tv_sec + ((double)start.tv_usec) * 1.0e-6;
    double endSeconds = end.tv_sec + ((double)end.tv_usec) * 1.0e-6;
    double updateTotal = (double)iterations * (double)updates * (double)pecount;

    printf("GUpdates    = %25.6f\n",
           (updateTotal / (1000.0 * 1000.0 * 1000.0)));
    printf("Time        = %25.6f\n", (endSeconds - startSeconds));
    printf("GUP/s       = %25.6f\n",
           ((updateTotal / (1000.0 * 1000.0 * 1000.0)) /
            (endSeconds - startSeconds)));
  }

  //////////////////////////////////////////////////////////////////

  long long int mytotal = 0;

  for (int i = 0; i < data_size; ++i) {
    mytotal += data[i];
  }

  for (int i = 0; i < pecount; ++i) {
    shmem_barrier_all();

    if (me == i) {
      printf("PE: %10d total is: %lld\n", i, mytotal);
    }

    shmem_barrier_all();
  }

  //////////////////////////////////////////////////////////////////

  shmem_barrier_all();

  if (0 == me) {
    printf("Deallocating shared memory regions.\n");
  }

  shfree(data);

  shmem_barrier_all();
  return 0;
}
