//
// Created by Teranishi, Keita on 9/13/21.
//

#ifndef LQCD_LQCD_H
#define LQCD_LQCD_H

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#define EVEN 0x02
#define ODD 0x01
#define EVENANDODD 0x03

#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7

#define NODIR -1  /* not a direction */

#define OPP_DIR(dir)    (7-(dir))       /* Opposite direction */
#define NDIRS 8                         /* number of directions */

/* defines for 3rd nearest neighbor (NAIK) stuff */
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15


#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir) ((dir)+8)
#define FORALL3UPDIR(dir) for(dir=X3UP; dir<=T3UP; dir++)

typedef struct __my_stats {
    int nx;
    int ny;
    int nz;
    int nt;
    int sizeof_su3matrix;
    int sizeof_su3vector;
    int nscompute;
    uint64_t pe_flops;
    uint64_t total_sites;
    uint64_t flops_per_iter;
    uint64_t flops_resid ;
    uint64_t flops_mmvs4d;
    double compute_nseconds_resid;
    double compute_nseconds_mmvs4d;

} lqParam;

#endif //LQCD_LQCD_H
