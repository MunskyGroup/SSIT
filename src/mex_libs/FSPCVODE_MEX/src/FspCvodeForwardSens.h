#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USENVECOMP
    #include <omp.h>
#endif

#include<sunlinsol/sunlinsol_spbcgs.h>
#include<sunlinsol/sunlinsol_spfgmr.h>
#include<sundials/sundials_nvector.h>
#include <cvodes/cvodes.h>                 /* prototypes for CVODE fcts., consts.  */
#include <cvodes/cvodes_spils.h>


#ifdef USENEVECOMP
#include <nvector/nvector_openmp.h>
#else
#include <nvector/nvector_serial.h>        /* access to serial N_Vector            */
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

typedef void (*dmatdpvecfun)( int n_state, int i_par, double t, double *x, double *y, void* user_data );

typedef void (*outputfun)( int i_tspan, double t, double *p, double **dp, void* output_mem);

typedef void (*stoppingfun)( double t, double *p, void *stopping_data, int *stop);

typedef struct mydata {
    int n;
    matvecfun mvfun;
    dmatdpvecfun dmvfun;
    void* mvdat;
    void* dmvdat;
};

/* Function to evaluate the RHS of the CME */
int CMERHSFun( double t, N_Vector y, N_Vector ydot, void *user_data );

/* Function to evaluate the action of the CME Jacobian on a vector */
int CMEJacVecFun( N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp );

/* Function to evaluate the RHS of the sensitivity CME, one parameter at a time */
int SensCMERHS1Fun( int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *my_data,
                    N_Vector tmp1, N_Vector tmp2 );

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
where the functions to evaluate the Jacobian action and the partial derivatives wrt sensitivity parameters are given from Matlab.
*/
int FspCVodeForwardSens( int n_tspan, double t_init, double *tspan, int n_state, int n_par, double *p0, double *dp0,
                         matvecfun matvec, void *matvecdat, dmatdpvecfun dmatvec, void *dmatvecdat, outputfun out_fun,
                         void *out_mem, stoppingfun stop_fun, void *stop_mem );

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
assuming initial time is zero.
*/
int FspCVodeForwardSensFromTimeZero( int n_tspan, double *tspan, int n_state, int n_par, double *p0, double *dp0,
                         matvecfun matvec, void *matvecdat, dmatdpvecfun dmatvec, void *dmatvecdat, outputfun out_fun,
                         void *out_mem, stoppingfun stop_fun, void *stop_mem );
