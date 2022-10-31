#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#include<sunlinsol/sunlinsol_spbcgs.h>
#include<sunlinsol/sunlinsol_spfgmr.h>
#include<sundials/sundials_nvector.h>
#include <cvodes/cvodes.h>                 /* prototypes for CVODE fcts., consts.  */
#include <cvodes/cvodes_spils.h>
#include <cvodes/cvodes_bandpre.h>

#ifndef USENVECOMP
#include <nvector/nvector_serial.h>        /* access to serial N_Vector            */
#else
#include <omp.h>
#include <nvector/nvector_openmp.h>
#endif

#include <sundials/sundials_types.h>       /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>        /* defs. of SUNRabs, SUNRexp, etc.      */
#ifndef NDEBUG
#define CVODECHKERR( flag ){\
    if (flag < 0) \
    {\
    printf("\nSUNDIALS_ERROR: function failed in file %s line %d with flag = %d\n\n",\
    __FILE__,__LINE__, flag);\
    return flag;\
    }\
    }
#else
#define CVODECHKERR(flag){while(0){}};
#endif

typedef void (*matvecfun)( int n_state, double t, double *x, double *y, void* user_data );

typedef void (*outputfun)( int i_tspan, double t, double *p, void* output_mem);

typedef void (*stoppingfun)( double t, double *p, void *stopping_data, int *stop);

typedef struct mydata {
    int n;
    matvecfun mvfun;
    void* mvdat;
    void* dmvdat;
};



/* Function to evaluate the RHS of the CME */
int CMERHSFun( double t, N_Vector y, N_Vector ydot, void *user_data );

/* Function to evaluate the action of the CME Jacobian on a vector */
int CMEJacVecFun( N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp );

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
assuming that the initial time is 0.
*/
int FspCVodeFromTimeZero( int n_tspan, double *tspan, int n_state, double *p0, matvecfun matvec, void *matvecdat,
                          outputfun out_fun,
                          void *out_mem, stoppingfun stop_fun, void *stop_mem );

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
with a non-zero initial time.
*/
int FspCVode( int n_tspan, double t_init, double *tspan, int n_state, double *p0, matvecfun matvec, void *matvecdat, outputfun out_fun,
              void *out_mem, stoppingfun stop_fun, void *stop_mem );
