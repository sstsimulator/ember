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

#include <errno.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

void get_position(const int rank, const int pex, const int pey, const int pez,
                  int* myX, int* myY, int* myZ) {
  const int plane = rank % (pex * pey);
  *myY = plane / pex;
  *myX = (plane % pex) != 0 ? (plane % pex) : 0;
  *myZ = rank / (pex * pey);
}

int convert_position_to_rank(const int pX, const int pY, const int pZ,
                             const int myX, const int myY, const int myZ) {
  // Check if we are out of bounds on the grid
  if ((myX < 0) || (myY < 0) || (myZ < 0) || (myX >= pX) || (myY >= pY) ||
      (myZ >= pZ)) {
    return -1;
  } else {
    return (myZ * (pX * pY)) + (myY * pX) + myX;
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int me = -1;
  int world = -1;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &world);

  int pex = world;
  int pey = 1;
  int pez = 1;

  int nx = 10;
  int ny = 10;
  int nz = 10;

  int repeats = 100;
  int vars = 1;

  long sleep = 1000;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-nx") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -nx without a value.\n");
        }

        exit(-1);
      }

      nx = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-ny") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -ny without a value.\n");
        }

        exit(-1);
      }

      ny = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-nz") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -nz without a value.\n");
        }

        exit(-1);
      }

      nz = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-pex") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -pex without a value.\n");
        }

        exit(-1);
      }

      pex = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-pey") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -pey without a value.\n");
        }

        exit(-1);
      }

      pey = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-pez") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -pez without a value.\n");
        }

        exit(-1);
      }

      pez = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-iterations") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -iterations without a value.\n");
        }

        exit(-1);
      }

      repeats = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-vars") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -vars without a value.\n");
        }

        exit(-1);
      }

      vars = atoi(argv[i + 1]);
      ++i;
    } else if (strcmp(argv[i], "-sleep") == 0) {
      if (i == argc) {
        if (me == 0) {
          fprintf(stderr, "Error: specified -sleep without a value.\n");
        }

        exit(-1);
      }

      sleep = atol(argv[i + 1]);
      ++i;
    } else {
      if (0 == me) {
        fprintf(stderr, "Unknown option: %s\n", argv[i]);
      }

      exit(-1);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if ((pex * pey * pez) != world) {
    if (0 == me) {
      fprintf(stderr, "Error: rank grid does not equal number of ranks.\n");
      fprintf(stderr, "%7d x %7d x %7d != %7d\n", pex, pey, pez, world);
    }

    exit(-1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (me == 0) {
    printf("# MPI Nearest Neighbor Communication\n");
    printf("# Info:\n");
    printf("# Processor Grid:         %7d x %7d x %7d\n", pex, pey, pez);
    printf("# Data Grid (per rank):   %7d x %7d x %7d\n", nx, ny, nz);
    printf("# Iterations:             %7d\n", repeats);
    printf("# Variables:              %7d\n", vars);
    printf("# Sleep:                  %7ld\n", sleep);
  }

  int posX, posY, posZ;
  get_position(me, pex, pey, pez, &posX, &posY, &posZ);

  const int xFaceUp =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY, posZ);
  const int xFaceDown =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY, posZ);
  const int yFaceUp =
      convert_position_to_rank(pex, pey, pez, posX, posY + 1, posZ);
  const int yFaceDown =
      convert_position_to_rank(pex, pey, pez, posX, posY - 1, posZ);
  const int zFaceUp =
      convert_position_to_rank(pex, pey, pez, posX, posY, posZ + 1);
  const int zFaceDown =
      convert_position_to_rank(pex, pey, pez, posX, posY, posZ - 1);

  const int vertexA =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY - 1, posZ - 1);
  const int vertexB =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY - 1, posZ + 1);
  const int vertexC =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY + 1, posZ - 1);
  const int vertexD =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY + 1, posZ + 1);
  const int vertexE =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY - 1, posZ - 1);
  const int vertexF =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY - 1, posZ + 1);
  const int vertexG =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY + 1, posZ - 1);
  const int vertexH =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY + 1, posZ + 1);

  const int edgeA =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY - 1, posZ);
  const int edgeB =
      convert_position_to_rank(pex, pey, pez, posX, posY - 1, posZ - 1);
  const int edgeC =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY - 1, posZ);
  const int edgeD =
      convert_position_to_rank(pex, pey, pez, posX, posY - 1, posZ + 1);
  const int edgeE =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY, posZ + 1);
  const int edgeF =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY, posZ + 1);
  const int edgeG =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY, posZ - 1);
  const int edgeH =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY, posZ - 1);
  const int edgeI =
      convert_position_to_rank(pex, pey, pez, posX - 1, posY + 1, posZ);
  const int edgeJ =
      convert_position_to_rank(pex, pey, pez, posX, posY + 1, posZ + 1);
  const int edgeK =
      convert_position_to_rank(pex, pey, pez, posX + 1, posY + 1, posZ);
  const int edgeL =
      convert_position_to_rank(pex, pey, pez, posX, posY + 1, posZ - 1);

  double send_vertexA = 1.0;
  double send_vertexB = 1.0;
  double send_vertexC = 1.0;
  double send_vertexD = 1.0;
  double send_vertexE = 1.0;
  double send_vertexF = 1.0;
  double send_vertexG = 1.0;
  double send_vertexH = 1.0;

  double recv_vertexA = 1.0;
  double recv_vertexB = 1.0;
  double recv_vertexC = 1.0;
  double recv_vertexD = 1.0;
  double recv_vertexE = 1.0;
  double recv_vertexF = 1.0;
  double recv_vertexG = 1.0;
  double recv_vertexH = 1.0;

  int requestcount = 0;
  MPI_Status* status;
  status = (MPI_Status*)malloc(sizeof(MPI_Status) * 52);

  MPI_Request* requests;
  requests = (MPI_Request*)malloc(sizeof(MPI_Request) * 52);

  double* edgeASendBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeBSendBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeCSendBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeDSendBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeESendBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeFSendBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeGSendBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeHSendBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeISendBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeJSendBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeKSendBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeLSendBuffer = (double*)malloc(sizeof(double) * nx * vars);

  double* edgeARecvBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeBRecvBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeCRecvBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeDRecvBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeERecvBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeFRecvBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeGRecvBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeHRecvBuffer = (double*)malloc(sizeof(double) * ny * vars);
  double* edgeIRecvBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeJRecvBuffer = (double*)malloc(sizeof(double) * nx * vars);
  double* edgeKRecvBuffer = (double*)malloc(sizeof(double) * nz * vars);
  double* edgeLRecvBuffer = (double*)malloc(sizeof(double) * nx * vars);

  for (int i = 0; i < nz; ++i) {
    edgeASendBuffer[i] = (double)i;
    edgeARecvBuffer[i] = 0.0;
    edgeCSendBuffer[i] = (double)i;
    edgeCRecvBuffer[i] = 0.0;
    edgeISendBuffer[i] = (double)i;
    edgeIRecvBuffer[i] = 0.0;
    edgeKSendBuffer[i] = (double)i;
    edgeKRecvBuffer[i] = 0.0;
  }

  for (int i = 0; i < ny; ++i) {
    edgeESendBuffer[i] = (double)i;
    edgeERecvBuffer[i] = 0.0;
    edgeFSendBuffer[i] = (double)i;
    edgeFRecvBuffer[i] = 0.0;
    edgeGSendBuffer[i] = (double)i;
    edgeGRecvBuffer[i] = 0.0;
    edgeHSendBuffer[i] = (double)i;
    edgeHRecvBuffer[i] = 0.0;
  }

  for (int i = 0; i < nx; ++i) {
    edgeBSendBuffer[i] = (double)i;
    edgeBRecvBuffer[i] = 0.0;
    edgeDSendBuffer[i] = (double)i;
    edgeDRecvBuffer[i] = 0.0;
    edgeJSendBuffer[i] = (double)i;
    edgeJRecvBuffer[i] = 0.0;
    edgeLSendBuffer[i] = (double)i;
    edgeLRecvBuffer[i] = 0.0;
  }

  double* xFaceUpSendBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);
  double* xFaceUpRecvBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);

  double* xFaceDownSendBuffer =
      (double*)malloc(sizeof(double) * ny * nz * vars);
  double* xFaceDownRecvBuffer =
      (double*)malloc(sizeof(double) * ny * nz * vars);

  for (int i = 0; i < ny * nz * vars; i++) {
    xFaceUpSendBuffer[i] = i;
    xFaceUpRecvBuffer[i] = i;
    xFaceDownSendBuffer[i] = i;
    xFaceDownRecvBuffer[i] = i;
  }

  double* yFaceUpSendBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);
  double* yFaceUpRecvBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);

  double* yFaceDownSendBuffer =
      (double*)malloc(sizeof(double) * nx * nz * vars);
  double* yFaceDownRecvBuffer =
      (double*)malloc(sizeof(double) * nx * nz * vars);

  for (int i = 0; i < nx * nz * vars; i++) {
    yFaceUpSendBuffer[i] = i;
    yFaceUpRecvBuffer[i] = i;
    yFaceDownSendBuffer[i] = i;
    yFaceDownRecvBuffer[i] = i;
  }

  double* zFaceUpSendBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);
  double* zFaceUpRecvBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);

  double* zFaceDownSendBuffer =
      (double*)malloc(sizeof(double) * nx * ny * vars);
  double* zFaceDownRecvBuffer =
      (double*)malloc(sizeof(double) * nx * ny * vars);

  for (int i = 0; i < nx * ny * vars; i++) {
    zFaceUpSendBuffer[i] = i;
    zFaceUpRecvBuffer[i] = i;
    zFaceDownSendBuffer[i] = i;
    zFaceDownRecvBuffer[i] = i;
  }

  struct timeval start;
  struct timeval end;

  struct timespec sleepTS;
  sleepTS.tv_sec = 0;
  sleepTS.tv_nsec = sleep;

  struct timespec remainTS;

  gettimeofday(&start, NULL);

  for (int i = 0; i < repeats; ++i) {
    requestcount = 0;

    if (nanosleep(&sleepTS, &remainTS) == EINTR) {
      while (nanosleep(&remainTS, &remainTS) == EINTR)
        ;
    }

    if (xFaceUp > -1) {
      MPI_Irecv(xFaceUpRecvBuffer, ny * nz * vars, MPI_DOUBLE, xFaceUp, 1000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(xFaceUpSendBuffer, ny * nz * vars, MPI_DOUBLE, xFaceUp, 1000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (xFaceDown > -1) {
      MPI_Irecv(xFaceDownRecvBuffer, ny * nz * vars, MPI_DOUBLE, xFaceDown,
                1000, MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(xFaceDownSendBuffer, ny * nz * vars, MPI_DOUBLE, xFaceDown,
                1000, MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (yFaceUp > -1) {
      MPI_Irecv(yFaceUpRecvBuffer, nx * nz * vars, MPI_DOUBLE, yFaceUp, 2000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(yFaceUpSendBuffer, nx * nz * vars, MPI_DOUBLE, yFaceUp, 2000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (yFaceDown > -1) {
      MPI_Irecv(yFaceDownRecvBuffer, nx * nz * vars, MPI_DOUBLE, yFaceDown,
                2000, MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(yFaceDownSendBuffer, nx * nz * vars, MPI_DOUBLE, yFaceDown,
                2000, MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (zFaceUp > -1) {
      MPI_Irecv(zFaceUpRecvBuffer, nx * ny * vars, MPI_DOUBLE, zFaceUp, 4000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(zFaceUpSendBuffer, nx * ny * vars, MPI_DOUBLE, zFaceUp, 4000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (zFaceDown > -1) {
      MPI_Irecv(zFaceDownRecvBuffer, nx * ny * vars, MPI_DOUBLE, zFaceDown,
                4000, MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(zFaceDownSendBuffer, nx * ny * vars, MPI_DOUBLE, zFaceDown,
                4000, MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeA > -1) {
      MPI_Irecv(edgeARecvBuffer, nz * vars, MPI_DOUBLE, edgeA, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeASendBuffer, nz * vars, MPI_DOUBLE, edgeA, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeB > -1) {
      MPI_Irecv(edgeBRecvBuffer, nx * vars, MPI_DOUBLE, edgeB, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeBSendBuffer, nx * vars, MPI_DOUBLE, edgeB, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeC > -1) {
      MPI_Irecv(edgeCRecvBuffer, nz * vars, MPI_DOUBLE, edgeC, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeCSendBuffer, nz * vars, MPI_DOUBLE, edgeC, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeD > -1) {
      MPI_Irecv(edgeDRecvBuffer, nx * vars, MPI_DOUBLE, edgeD, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeDSendBuffer, nx * vars, MPI_DOUBLE, edgeD, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeE > -1) {
      MPI_Irecv(edgeERecvBuffer, ny * vars, MPI_DOUBLE, edgeE, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeESendBuffer, ny * vars, MPI_DOUBLE, edgeE, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeF > -1) {
      MPI_Irecv(edgeFRecvBuffer, ny * vars, MPI_DOUBLE, edgeF, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeFSendBuffer, ny * vars, MPI_DOUBLE, edgeF, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeG > -1) {
      MPI_Irecv(edgeARecvBuffer, ny * vars, MPI_DOUBLE, edgeG, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeASendBuffer, ny * vars, MPI_DOUBLE, edgeG, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeH > -1) {
      MPI_Irecv(edgeARecvBuffer, ny * vars, MPI_DOUBLE, edgeH, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeASendBuffer, ny * vars, MPI_DOUBLE, edgeH, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeI > -1) {
      MPI_Irecv(edgeIRecvBuffer, nz * vars, MPI_DOUBLE, edgeI, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeISendBuffer, nz * vars, MPI_DOUBLE, edgeI, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeJ > -1) {
      MPI_Irecv(edgeJRecvBuffer, nx * vars, MPI_DOUBLE, edgeJ, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeJSendBuffer, nx * vars, MPI_DOUBLE, edgeJ, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeK > -1) {
      MPI_Irecv(edgeKRecvBuffer, nz * vars, MPI_DOUBLE, edgeK, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeKSendBuffer, nz * vars, MPI_DOUBLE, edgeK, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    if (edgeL > -1) {
      MPI_Irecv(edgeLRecvBuffer, nx * vars, MPI_DOUBLE, edgeL, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
      MPI_Isend(edgeLSendBuffer, nx * vars, MPI_DOUBLE, edgeL, 8000,
                MPI_COMM_WORLD, &requests[requestcount++]);
    }

    MPI_Waitall(requestcount, requests, status);
    requestcount = 0;
  }

  gettimeofday(&end, NULL);

  MPI_Barrier(MPI_COMM_WORLD);

  free(xFaceUpRecvBuffer);
  free(xFaceDownRecvBuffer);
  free(yFaceUpRecvBuffer);
  free(yFaceDownRecvBuffer);
  free(zFaceUpRecvBuffer);
  free(zFaceDownRecvBuffer);

  if (convert_position_to_rank(pex, pey, pez, pex / 2, pey / 2, pez / 2) ==
      me) {
    printf("# Results from rank: %d\n", me);

    const double timeTaken =
        (((double)end.tv_sec) + ((double)end.tv_usec) * 1.0e-6) -
        (((double)start.tv_sec) + ((double)start.tv_usec) * 1.0e-6);
    const double bytesXchng =
        ((double)(xFaceUp > -1 ? sizeof(double) * ny * nz * 2 * vars : 0)) +
        ((double)(xFaceDown > -1 ? sizeof(double) * ny * nz * 2 * vars : 0)) +
        ((double)(yFaceUp > -1 ? sizeof(double) * nx * nz * 2 * vars : 0)) +
        ((double)(yFaceDown > -1 ? sizeof(double) * nx * nz * 2 * vars : 0)) +
        ((double)(zFaceUp > -1 ? sizeof(double) * nx * ny * 2 * vars : 0)) +
        ((double)(zFaceDown > -1 ? sizeof(double) * nx * ny * 2 * vars : 0));

    printf("# %20s %20s %20s\n", "Time", "KBytesXchng/Rank-Max", "MB/S/Rank");
    printf("  %20.6f %20.4f %20.4f\n", timeTaken, bytesXchng / 1024.0,
           (bytesXchng / 1024.0) / timeTaken);
  }

  MPI_Finalize();
}
