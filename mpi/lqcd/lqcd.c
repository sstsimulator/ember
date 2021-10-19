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


void initialize( int *nx, int *ny, int *nz, int *nt, int *n_ranks, int num_nodes, int my_rank )
{

    MPI_Comm_rank( &my_rank, MPI_COMM_WORLD );
    MPI_Comm_size( &nsize, MPI_COMM_WORLD);

    int sites_on_node = nx*ny*nz*nt/num_nodes;
    int even_sites_on_node = sites_on_node/2;
    int odd_sites_on_node = sites_on_node/2;

    // 3X3 matrix of complex values -> 18 (8 Byte) doubles
    // This represents the long links or fat links aka gluons aka force carriers
    sizeof_su3matrix = 18*8;

    // 3 component complex vector -> 6 (8 Byte) doubles
    // represents the matter fields (quarks)
    suzeof_su3vector = 6*8;

    // With data compression, a Xeon Phi (KNL) could get 300 GF/s
    // see "Staggered Dslash Performance on Intel Xeon Phi Architecture"
    uint64_t pe_flops = (uint64_t) 20ULL* 1000*1000*1000; // Subject to change

    const uint64_t total_sites = (uint64_t) (nx*ny*nz*nt);

    // The MILC code for CG gave this calculation sites_on_node*(offset*15+1205)
    // Offset is also referred to as "number of masses" with referenced values being 9 or 11.
    // Some algorithms compress data before writing it to memory.
    // The equation below accounts only for "useful" flops.
    // Turn this into an optional argument
    //const uint64_t flops_per_iter = sites_on_node*(11*15 + 1205)/2;
    const uint64_t flops_per_iter = 0;

    // From the portion of code we are looking at there is 157(sites_on_node)/2 for the resid calc
    // (12*num_offsets+12 +12 +1) * V/2
    const uint64_t flops_resid = sites_on_node*157/2;

    // Per SG: multiplication of 3X3 matrix by 3-vector:
    // 36 multiplies and 30 adds, i.e. 66 flops.
    // 3-vector add is 6 adds, so one call to this function is 4*72 flops.
    // Each call is, therefore, 288 flops and there is loop over sites,
    // but that loop is only over one parity.
    // But, there are two calls to mult_su3_mat_vec_sum_4dir, once with fatback4 and once with longback4.
    const uint64_t flops_mmvs4d = sites_on_node*288;

    // This is overly simplified, but divides the flops per iteration evenly across the compute segments
    // You have 4 approximately equal compute segments after the even/odd, pos/neg gathers,
    // then the resid calculation which we assign two segments to.
    // Getting computation segments to represent all architectures optimizations is not possible.
    // Arithmetic intensity can vary widely machine to machine.
    // This motif represents a starting point, given our profiled systems.

    // Converts FLOP/s into nano seconds of compute
    compute_nseconds_resid = ( flops_resid / ( pe_flops / 1000000000.0 ) );
    compute_nseconds_mmvs4d = ( flops_mmvs4d / ( pe_flops / 1000000000.0 ) );

    if(flops_per_iter){
        const uint64_t compute_nseconds = ( flops_per_iter / ( pe_flops / 1000000000.0 ) );
        nsCompute  = (uint64_t) params.find("arg.computetime", (uint64_t) compute_nseconds);
        output("nsCompute set:%" PRIu64 " ns\n", nsCompute);
    }
    else{
        nsCompute = 0;
        if ( my_rank == 0){
            output("LQCD compute time resid: %" PRIu64 " mmvs4d: %" PRIu64 " ns\n", compute_nseconds_resid, compute_nseconds_mmvs4d);
        }
    }

    iterations = (uint32_t) params.find("arg.iterations", 1);

    coll_start = 0;
    coll_time = 0;
    comp_start = 0;
    comp_time = 0;
    gather_pos_start = 0;
    gather_pos_time = 0;
    gather_neg_start = 0;
    gather_neg_time = 0;

    n_ranks[XDOWN] = -1;
    n_ranks[XUP]   = -1;
    n_ranks[YDOWN] = -1;
    n_ranks[YUP]   = -1;
    n_ranks[ZDOWN] = -1;
    n_ranks[ZUP]   = -1;
    n_ranks[TDOWN] = -1;
    n_ranks[TUP]   = -1;

    // configure();
}

//returns the rank that holds the site x,y,z,t
int get_node_index( int x, int y, int z, int t, int num_nodes)
{
    int i = -1;
    i = x+ nx * ( y + ny * (z+nz ) * t );
    return i % num_nodes;
}

//copies coordinates with an offset in the specified dimension
void EmberLQCDGenerator::neighbor_coords(const int coords[], int n_coords[], const int dim, const int offset){
    n_coords[XUP] = coords[XUP];
    n_coords[YUP] = coords[YUP];
    n_coords[ZUP] = coords[ZUP];
    n_coords[TUP] = coords[TUP];
    n_coords[dim] = n_coords[dim]+offset;
}


/*------------------------------------------------------------------*/
/* Convert rank to machine coordinates */
void EmberLQCDGenerator::lex_coords(int coords[], const int dim, const int size[], const uint32_t rank){
    int d;
    uint32_t r = rank;

    for(d = 0; d < dim; d++){
        coords[d] = r % size[d];
        r /= size[d];
    }
}


/*------------------------------------------------------------------*/
/* Parity of the coordinate */
int EmberLQCDGenerator::coord_parity(int r[]){
    return (r[0] + r[1] + r[2] + r[3]) % 2;
}

/*------------------------------------------------------------------*/
/* Convert machine coordinate to linear rank (inverse of
   lex_coords) */
//added two checks to return -1 if we are trying to access a neighbor that
//is beyond the scope of the lattice.

int EmberLQCDGenerator::lex_rank(const int coords[], int dim, int size[])
{
    int d;
    int x = coords[XUP];
    int y = coords[YUP];
    int z = coords[ZUP];
    int t = coords[TUP];
    size_t rank = coords[dim-1];
    if (x == -1 || y == -1 || z == -1 || t == -1){
        return(-1);
    }
    if (x >= nsquares[XUP] || y >= nsquares[YUP] || z >= nsquares[ZUP] || t >= nsquares[TUP]){
        return(-1);
    }

    for(d = dim-2; d >= 0; d--){
        rank = rank * size[d] + coords[d];
    }
    return rank;
}

void EmberLQCDGenerator::setup_hyper_prime(){
    int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
    int dir, i, j, k;
    int len_primes = sizeof(prime)/sizeof(int);
    /* Figure out dimensions of rectangle */
    squaresize[XUP] = nx; squaresize[YUP] = ny;
    squaresize[ZUP] = nz; squaresize[TUP] = nt;
    nsquares[XUP] = nsquares[YUP] = nsquares[ZUP] = nsquares[TUP] = 1;

    i = 1;      /* current number of hypercubes */
    while(i<num_nodes){
        /* figure out which prime to divide by starting with largest */
        k = len_primes-1;
        while( (num_nodes/i)%prime[k] != 0 && k>0 ) --k;
        /* figure out which direction to divide */

        /* find largest dimension of h-cubes divisible by prime[k] */
        for(j=0,dir=XUP;dir<=TUP;dir++)
            if( squaresize[dir]>j && squaresize[dir]%prime[k]==0 )
                j=squaresize[dir];

        /* if one direction with largest dimension has already been
           divided, divide it again.  Otherwise divide first direction
           with largest dimension. */
        for(dir=XUP;dir<=TUP;dir++)
            if( squaresize[dir]==j && nsquares[dir]>1 )break;
        if( dir > TUP)for(dir=XUP;dir<=TUP;dir++)
                if( squaresize[dir]==j )break;
        /* This can fail if I run out of prime factors in the dimensions */
        if(dir > TUP){
            if(rank() == 0)
                printf("LAYOUT: Can't lay out this lattice, not enough factors of %d\n"
                        ,prime[k]);
        }

        /* do the surgery */
        i*=prime[k]; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
    }
    if(0 == rank()) {
        output("LQCD nsquares size: x%" PRIu32 "y%" PRIu32 "z%" PRIu32 "t%" PRIu32"\n", nsquares[XUP], nsquares[YUP], nsquares[ZUP], nsquares[TUP]);
        output("LQCD square size: x%" PRIu32 "y%" PRIu32 "z%" PRIu32 "t%" PRIu32"\n", squaresize[XUP], squaresize[YUP], squaresize[ZUP], squaresize[TUP]);
    }
}
//returns the node number on which a site lives.
int EmberLQCDGenerator::node_number_from_lex(int x, int y, int z, int t) {
    int i;
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return(i);
}

//returns the surface volume to be transfered in a gather (1st neighbor)
int EmberLQCDGenerator::get_transfer_size(int dimension){
    int sitestrans = 1;
    int dim = dimension;
    // squaresize[neg] should == squaresize[pos]
    if (dimension >= TDOWN){
        dim = OPP_DIR(dim);
    }
    // for positive directions
    for (int d = XUP; d <= TUP; d++){
        if (d != dim){
            sitestrans *= squaresize[d];
        }
    }
    return sitestrans;
}

//returns the node number on which a site lives (given lattice coords).
/*------------------------------------------------------------------*/
int EmberLQCDGenerator::node_number(int x, int y, int z, int t) {
    int i;
    x /= squaresize[XUP]; y /= squaresize[YUP];
    z /= squaresize[ZUP]; t /= squaresize[TUP];
    i = x + nsquares[XUP]*( y + nsquares[YUP]*( z + nsquares[ZUP]*( t )));
    return( i );
}

void EmberLQCDGenerator::configure()
{

    //determine the problem size given to each node
    //code from MILC setup_hyper_prime()
    setup_hyper_prime();

    //the coordinates of this rank in pe_x, pe_y, pe_z, pe_t
    lex_coords(machine_coordinates, 4, nsquares, rank());

    /* Number of sites on node */
    sites_on_node = squaresize[XUP]*squaresize[YUP]*squaresize[ZUP]*squaresize[TUP];
    even_sites_on_node = sites_on_node/2;
    odd_sites_on_node = sites_on_node/2;
    if(0 == rank()) {
        output("LQCD problem size: x%" PRIu32 "y%" PRIu32 "z%" PRIu32 "t%" PRIu32"\n", nx, ny, nz, nt);
        output("LQCD gather num_elements (1st and Naik): %" PRIu64  " %" PRIu64 "\n", (uint64_t) nx*ny*nz*9, (uint64_t) nx*ny*nz*27);
        if (nsCompute != 0) output("LQCD compute time per segment: %" PRIu64 " ns\n", nsCompute);
        else output("LQCD compute time resid: %" PRIu64 " mmvs4d: %" PRIu64 " ns\n", compute_nseconds_resid, compute_nseconds_mmvs4d);
        output("LQCD iterations: %" PRIu32 "\n", iterations);
        output("LQCD sites per node: %" PRIu32 "\n", sites_on_node);
    }
    //determine rank of neighbors in each direction x,y,z,t
    int tmp_coords[4];
    //positive direction
    for (int d = XUP; d <= TUP; d++){
        neighbor_coords(machine_coordinates, tmp_coords, d, 1);
        n_ranks[d] = lex_rank(tmp_coords, 4, nsquares);
    }
    //negative direction
    for (int d = XUP; d <= TUP; d++){
        neighbor_coords(machine_coordinates, tmp_coords, d, -1);
        n_ranks[OPP_DIR(d)] = lex_rank(tmp_coords, 4, nsquares);
    }

    verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", World=%" PRId32 ", X=%" PRId32 ", Y=%" PRId32 ", Z=%" PRId32 ", T=%" PRId32 ", Px=%" PRId32 ", Py=%" PRId32 ", Pz=%" PRId32 ", Pt=%" PRId32 "\n",
            rank(), size(),
            machine_coordinates[XUP], machine_coordinates[YUP], machine_coordinates[ZUP], machine_coordinates[TUP],
            nsquares[XUP],nsquares[YUP],nsquares[ZUP],nsquares[TUP]);
    verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", X+: %" PRId32 ", X-: %" PRId32 "\n", rank(), n_ranks[XUP], n_ranks[XDOWN]);
    verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", Y+: %" PRId32 ", Y-: %" PRId32 "\n", rank(), n_ranks[YUP], n_ranks[YDOWN]);
    verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", Z+: %" PRId32 ", Z-: %" PRId32 "\n", rank(), n_ranks[ZUP], n_ranks[ZDOWN]);
    verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", T+: %" PRId32 ", T-: %" PRId32 "\n", rank(), n_ranks[TUP], n_ranks[TDOWN]);
}



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
   int num_nodes = -1;
   int nx = 64;
   int ny = 64;
   int nz = 64;
   int nt = 64;
   int repeats = 100;
   int vars = 1;
   long sleep = 1000;
   int compt;
   uint64_t peflops;
   int n_ranks[8];
   MPI_Request *posRequests;
   MPI_Request *negRequests;
   MPI_Status  *status_list;
   MPI_Status  status;
   double * UpRecvBuffer[8];
   double * UpSendBuffer[8];
   double * DownRecvBuffer[8];
   double * DownSendBuffer[8];

   MPI_Comm_rank(MPI_COMM_WORLD, &me);
   MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);



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
      } else if (strcmp(argv[i], "-nt") == 0) {
          if (i == argc) {
              if (me == 0) {
                  fprintf(stderr, "Error: specified -nz without a value.\n");
              }

              exit(-1);
          }

          nt = atoi(argv[i + 1]);
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
      } else if (strcmp(argv[i], "-peflops") == 0) {
         if (i == argc) {
            if (me == 0) {
               fprintf(stderr, "Error: specified -vars without a value.\n");
            }

            exit(-1);
         }

         peflops = atoi(argv[i + 1]);
         ++i;
      } else if (strcmp(argv[i], "-compt") == 0) {
         if (i == argc) {
            if (me == 0) {
               fprintf(stderr, "Error: specified -sleep without a value.\n");
            }

            exit(-1);
         }

         compt = atol(argv[i + 1]);
         ++i;
      } else {
         if (0 == me) {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
         }

         exit(-1);
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   initialize( &nx, &ny, &nz, &nt, n_ranks, &me, &num_nodes, &posRequests, &negRequests );
   // Even/Odd preconditioning, organizes sites into two sets (if enabled odd_even == 2)
   // If enabled we do two loops of this function with one-half transfer size
   // int even_odd = 1;
   int even_odd = 2;
   // Enqueue gather from positive and negative directions (8 total)  su3_vector

   int transsz = 500;

   for (int parity = 1; parity <= even_odd; parity++) {

       MPI_Irecv(xUpRecvBuffer, ny * nz * vars, MPI_DOUBLE, xUp, 1000,
                 MPI_COMM_WORLD, &requests[requestcount++]);
       MPI_Isend(xUpSendBuffer, ny * nz * vars, MPI_DOUBLE, xUp, 1000,
                 MPI_COMM_WORLD, &requests[requestcount++]);

       for (int d = XUP; d <= TUP; d++){
           if (n_ranks[d] != -1){
               // pos_requests.push_back(req);
               int transsz = get_transfer_size(d)/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "receiving pos elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Irecv( xUpRecvBuffer, n_ranks[d], sizeof_su3vector * transsz, 0, MPI_Comm_WORLD &posRequests[requestcount++]);  // GroupWorld
           }
       }
       for (int d = XUP; d <= TUP; d++){
           if ( n_ranks[d] != -1 ){
               transsz = get_transfer_size(d)/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "sending pos elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Send( xUpSendBuffer, n_ranks[d] ,sizeof_su3vector * transsz, 0, MPI_COMM_WORLD);
           }
       }

       // And start the 3-step gather (4 total)
       // In the Greg Bauer MILC performance model paper they list trans_sz as
       // Fetch Naik term or three link term (points 1, 2 and 3 planes into the next domain)
       // multiplied by 6 elements multiplied by surface volume.
       // We don't need to gather the first site twice (hence get_transfer_size(d)*2)
       for (int d = XUP; d <= TUP; d++){
           if (n_ranks[d] != -1){
               int transsz = get_transfer_size(d)*2/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "receiving pos (Naik) elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Irecv( xUpRecvBuffer, n_ranks[d], sizeof_su3vector * transsz, 0, MPI_Comm_WORLD &posRequests[requestcount++]);  // GroupWorld
           }
       }
       for (int d = XUP; d <= TUP; d++){
           if ( n_ranks[d] != -1 ){
               int transsz = get_transfer_size(d)*2/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "sending pos (Naik) elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Send( xUpSendBuffer, n_ranks[d] ,sizeof_su3vector * transsz, 0, MPI_COMM_WORLD);
           }
       }


       // Gather from negative directions (4 total)
       for (int d = TDOWN; d <= XDOWN; d++){
           if (n_ranks[d] != -1){
               //MessageRequest*  req  = new MessageRequest();
               // neg_requests.push_back(req);
               int transsz = get_transfer_size(d)/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "receiving neg elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Irecv( xDownRecvBuffer, n_ranks[d], sizeof_su3vector * transsz, 0, MPI_Comm_WORLD &negRequests[neg_requestcount++]);  // GroupWorld
           }
       }
       for (int d = TDOWN; d <= XDOWN; d++){
           if ( n_ranks[d] != -1 ){
               int transsz = get_transfer_size(d)/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "sending neg elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Send( xDownSendBuffer, n_ranks[d] ,sizeof_su3vector * transsz, 0, MPI_COMM_WORLD);
           }
       }

       // Time-wise 3-step gather (4 total)
       for (int d = TDOWN; d <= XDOWN; d++){
           if (n_ranks[d] != -1){
               // MessageRequest*  req  = new MessageRequest();
               // neg_requests.push_back(req);
               transsz = get_transfer_size(d)*2/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "receiving neg (Naik) elements: %" PRIu32 "\n", rank(), transsz);
               enQ_irecv( evQ, n_ranks[d],
                          sizeof_su3vector * transsz, 0, GroupWorld, req);
               MPI_Irecv( xDownRecvBuffer, n_ranks[d], sizeof_su3vector * transsz, 0, MPI_Comm_WORLD &negRequests[neg_requestcount++]);  // GroupWorld
           }
       }
       for (int d = TDOWN; d <= XDOWN; d++){
           if ( n_ranks[d] != -1 ){
               transsz = get_transfer_size(d)*2/even_odd;
               if (rank() == 0 || n_ranks[d] == 0)
                   verbose(CALL_INFO, 2, 0, "Rank: %" PRIu32 "sending neg (Naik) elements: %" PRIu32 "\n", rank(), transsz);
               MPI_Send( xDownSendBuffer, n_ranks[d] ,sizeof_su3vector * transsz, 0, MPI_COMM_WORLD);
           }
       }

       // Wait on positive gathers (8 total)
       verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", +Messages to wait on: %zu\n", rank(), pos_requests.size());
       MPI_Waitall(pos_requestcount, pos_requests, status);

       // Compute Matrix Vector Multiply Sum in 4 directions
       // C <- A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]
       // This is done for both positive "fat" and "long" links

#if 0
       comp_start = Simulation::getSimulation()->getCurrentSimCycle();
       //output("Rank %" PRIu32 ", Compute start: %" PRIu64 "\n", rank(), comp_start);
       if (nsCompute) enQ_compute(evQ, nsCompute);
       else enQ_compute(evQ, compute_nseconds_mmvs4d);
       comp_time += Simulation::getSimulation()->getCurrentSimCycle() - comp_start;
       //output("Rank %" PRIu32 ", Compute end: %" PRIu64 "\n", rank(), Simulation::getSimulation()->getCurrentSimCycle());
#endif
       //Wait on negative gathers (8 total)
       verbose(CALL_INFO, 1, 0, "Rank: %" PRIu32 ", -Messages to wait on: %zu\n", rank(), neg_requests.size());
       MPI_Waitall(neg_requestcount, neg_requests, status);
       // Compute Matrix Vector Multiply Sum in 4 directions
       // C <- A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]
       // This is done for both negative "fat" and "long" links

#if 0
       if (nsCompute) enQ_compute(evQ, nsCompute);
       else enQ_compute(evQ, compute_nseconds_mmvs4d);
#endif
       //output("Iteration on rank %" PRId32 " completed generation, %d events in queue\n",
       //    rank(), (int)evQ.size());
   } // END even-odd preconditioning loop

   // Allreduce MPI_DOUBLE
   // coll_start = Simulation::getSimulation()->getCurrentSimCycle();
   double colldata = 1.0;
   double colldata_out;
   MPI_Allreduce( &colldata, &colldata_out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   // coll_time += Simulation::getSimulation()->getCurrentSimCycle() - coll_start;

   // Compute ResidSq
   // Our profiling shows this is ~ 1/4th the FLOPS of the
   // compute segments above, but is more memory intensive (and takes longer).
   // Depending on whether the problem is small enough to utilize memory
   // effectively, the time spent here can vary widely.

   if (nsCompute) {
       // enQ_compute(evQ, nsCompute);
       // sleep();
   } else {
       //enQ_compute(evQ, compute_nseconds_resid);
       // sleep();
   }

    // Allreduce MPI_DOUBLE
    // coll_start = Simulation::getSimulation()->getCurrentSimCycle();
    //output("Rank %" PRIu32 ", Collective start: %" PRIu64 "\n", rank(), coll_start);
    MPI_Allreduce( &colldata, &colldata_out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    //coll_time += Simulation::getSimulation()->getCurrentSimCycle() - coll_start;
    //output("Rank %" PRIu32 ", Collective end: %" PRIu64 "\n", rank(), Simulation::getSimulation()->getCurrentSimCycle());


#if 0
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

   int xUp = convert_position_to_rank(pex, pey, pez, posX + 1, posY, posZ);
   int xDown = convert_position_to_rank(pex, pey, pez, posX - 1, posY, posZ);
   int yUp = convert_position_to_rank(pex, pey, pez, posX, posY + 1, posZ);
   int yDown = convert_position_to_rank(pex, pey, pez, posX, posY - 1, posZ);
   int zUp = convert_position_to_rank(pex, pey, pez, posX, posY, posZ + 1);
   int zDown = convert_position_to_rank(pex, pey, pez, posX, posY, posZ - 1);

   int requestcount = 0;
   MPI_Status* status;
   status = (MPI_Status*)malloc(sizeof(MPI_Status) * 4);

   MPI_Request* requests;
   requests = (MPI_Request*)malloc(sizeof(MPI_Request) * 4);

   double* xUpSendBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);
   double* xUpRecvBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);

   double* xDownSendBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);
   double* xDownRecvBuffer = (double*)malloc(sizeof(double) * ny * nz * vars);

   for (int i = 0; i < ny * nz * vars; i++) {
      xUpSendBuffer[i] = i;
      xUpRecvBuffer[i] = i;
      xDownSendBuffer[i] = i;
      xDownRecvBuffer[i] = i;
   }

   double* yUpSendBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);
   double* yUpRecvBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);

   double* yDownSendBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);
   double* yDownRecvBuffer = (double*)malloc(sizeof(double) * nx * nz * vars);

   for (int i = 0; i < nx * nz * vars; i++) {
      yUpSendBuffer[i] = i;
      yUpRecvBuffer[i] = i;
      yDownSendBuffer[i] = i;
      yDownRecvBuffer[i] = i;
   }

   double* zUpSendBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);
   double* zUpRecvBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);

   double* zDownSendBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);
   double* zDownRecvBuffer = (double*)malloc(sizeof(double) * nx * ny * vars);


   for (int i = 0; i < nx * ny * vars; i++) {
      zUpSendBuffer[i] = i;
      zUpRecvBuffer[i] = i;
      zDownSendBuffer[i] = i;
      zDownRecvBuffer[i] = i;
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

      if (xUp > -1) {
         MPI_Irecv(xUpRecvBuffer, ny * nz * vars, MPI_DOUBLE, xUp, 1000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(xUpSendBuffer, ny * nz * vars, MPI_DOUBLE, xUp, 1000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      if (xDown > -1) {
         MPI_Irecv(xDownRecvBuffer, ny * nz * vars, MPI_DOUBLE, xDown, 1000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(xDownSendBuffer, ny * nz * vars, MPI_DOUBLE, xDown, 1000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      MPI_Waitall(requestcount, requests, status);
      requestcount = 0;

      if (yUp > -1) {
         MPI_Irecv(yUpRecvBuffer, nx * nz * vars, MPI_DOUBLE, yUp, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(yUpSendBuffer, nx * nz * vars, MPI_DOUBLE, yUp, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      if (yDown > -1) {
         MPI_Irecv(yDownRecvBuffer, nx * nz * vars, MPI_DOUBLE, yDown, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(yDownSendBuffer, nx * nz * vars, MPI_DOUBLE, yDown, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      MPI_Waitall(requestcount, requests, status);
      requestcount = 0;


      if (yUp > -1) {
         MPI_Irecv(yUpRecvBuffer, nx * nz * vars, MPI_DOUBLE, yUp, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(yUpSendBuffer, nx * nz * vars, MPI_DOUBLE, yUp, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      if (yDown > -1) {
         MPI_Irecv(yDownRecvBuffer, nx * nz * vars, MPI_DOUBLE, yDown, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(yDownSendBuffer, nx * nz * vars, MPI_DOUBLE, yDown, 2000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      MPI_Waitall(requestcount, requests, status);
      requestcount = 0;

      if (zUp > -1) {
         MPI_Irecv(zUpRecvBuffer, nx * ny * vars, MPI_DOUBLE, zUp, 4000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(zUpSendBuffer, nx * ny * vars, MPI_DOUBLE, zUp, 4000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      if (zDown > -1) {
         MPI_Irecv(zDownRecvBuffer, nx * ny * vars, MPI_DOUBLE, zDown, 4000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
         MPI_Isend(zDownSendBuffer, nx * ny * vars, MPI_DOUBLE, zDown, 4000,
                   MPI_COMM_WORLD, &requests[requestcount++]);
      }

      MPI_Waitall(requestcount, requests, status);
      requestcount = 0;
   }

   gettimeofday(&end, NULL);

   MPI_Barrier(MPI_COMM_WORLD);

   free(xUpRecvBuffer);
   free(xDownRecvBuffer);
   free(yUpRecvBuffer);
   free(yDownRecvBuffer);
   free(zUpRecvBuffer);
   free(zDownRecvBuffer);
   if (convert_position_to_rank(pex, pey, pez, pex / 2, pey / 2, pez / 2) ==
       me) {
      printf("# Results from rank: %d\n", me);

      const double timeTaken =
            (((double)end.tv_sec) + ((double)end.tv_usec) * 1.0e-6) -
            (((double)start.tv_sec) + ((double)start.tv_usec) * 1.0e-6);
      const double bytesXchng =
            ((double)(xUp > -1 ? sizeof(double) * ny * nz * 2 * vars : 0)) +
            ((double)(xDown > -1 ? sizeof(double) * ny * nz * 2 * vars : 0)) +
            ((double)(yUp > -1 ? sizeof(double) * nx * nz * 2 * vars : 0)) +
            ((double)(yDown > -1 ? sizeof(double) * nx * nz * 2 * vars : 0)) +
            ((double)(zUp > -1 ? sizeof(double) * nx * ny * 2 * vars : 0)) +
            ((double)(zDown > -1 ? sizeof(double) * nx * ny * 2 * vars : 0));

      printf("# %20s %20s %20s\n", "Time", "KBytesXchng/Rank-Max", "MB/S/Rank");
      printf("  %20.6f %20.4f %20.4f\n", timeTaken, bytesXchng / 1024.0,
             (bytesXchng / 1024.0) / timeTaken);
   }
#endif
   MPI_Finalize();
}
