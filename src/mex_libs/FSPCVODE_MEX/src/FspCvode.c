#include "FspCvode.h"

/* Function to evaluate the RHS of the CME */
int CMERHSFun( double t, N_Vector y, N_Vector ydot, void *user_data ) {
    struct mydata *f_data = ( struct mydata * ) user_data;
    f_data->mvfun( f_data->n, t, N_VGetArrayPointer( y ), N_VGetArrayPointer( ydot ), f_data->mvdat );
    return 0;
}

/* Function to evaluate the action of the CME Jacobian on a vector */
int CMEJacVecFun( N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp ) {
    struct mydata *f_data = ( struct mydata * ) user_data;
    realtype *vptr = N_VGetArrayPointer( v );
    realtype *Jvptr = N_VGetArrayPointer( Jv );
    f_data->mvfun( f_data->n, t, vptr, Jvptr, f_data->mvdat );
    return 0;
}

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
assuming that initial time is 0.
*/
int FspCVodeFromTimeZero( int n_tspan, double *tspan, int n_state, double *p0, matvecfun matvec, void *matvecdat,
                          outputfun out_fun,
                          void *out_mem, stoppingfun stop_fun, void *stop_mem ) {
    return FspCVode(n_tspan, 0.0, tspan, n_state, p0, matvec, matvecdat, out_fun, out_mem, stop_fun, stop_mem);
}

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
with a non-zero initial time.
*/
int FspCVode( int n_tspan, double t_init, double *tspan, int n_state, double *p0, matvecfun matvec, void *matvecdat, outputfun out_fun,
              void *out_mem, stoppingfun stop_fun, void *stop_mem ){
    /* CVODEs error code */
    int cvode_stat;
    /* CVODE error tolerances */
    double Atol = 1.0e-8; // Absolute tolerance for the solution
    double Rtol = 1.0e-4; // Relative tolerance for the solution
    /* CVODES object */
    void *cvode_mem;

    /* Sundials linear solver */
    SUNLinearSolver linear_solver;

    /* Initialize initial N_Vector solution */
    N_Vector p_nv;

    #ifndef USENVECOMP
    p_nv = N_VNew_Serial( n_state );
    #else
    int nthreads = omp_get_max_threads();
    p_nv = N_VNew_OpenMP(n_state, nthreads);
    #endif
    realtype *p_nv_dat = N_VGetArrayPointer( p_nv );
    for ( int i = 0; i < n_state; ++i ) {
        p_nv_dat[ i ] = p0[ i ];
    }

    /* Initialize data for user-supplied functions */
    struct mydata my_data;
    my_data.n = n_state;
    my_data.mvfun = matvec;
    my_data.mvdat = matvecdat;

    /* Initialize the main integrator for states */
    cvode_mem = CVodeCreate( CV_BDF );
    cvode_stat = CVodeInit( cvode_mem, &CMERHSFun, t_init, p_nv );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSetUserData( cvode_mem, ( void * ) &my_data );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSStolerances( cvode_mem, Rtol, Atol );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSetMaxNumSteps( cvode_mem, 10000 );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSetMaxConvFails( cvode_mem, 1000 );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSetMaxNonlinIters( cvode_mem, 1000 );
    CVODECHKERR( cvode_stat );

    /* Create the linear solver without preconditioning */
    linear_solver = SUNSPFGMR( p_nv, PREC_NONE, 50 );
    cvode_stat = CVSpilsSetLinearSolver( cvode_mem, linear_solver );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVSpilsSetJacTimes( cvode_mem, NULL, &CMEJacVecFun );
    CVODECHKERR( cvode_stat );

    /* Advance to the end, output vectors at times specified in tspan */
    double tret;
    double tnow = t_init;
    double tfinal = tspan[n_tspan-1];
    int istop = 0;
    int iout = 0;
    while ( tnow < tfinal ) {
        cvode_stat = CVode( cvode_mem, tfinal, p_nv, &tret, CV_ONE_STEP );
        CVODECHKERR( cvode_stat );
        if (tret > tfinal){
            tret = tfinal;
            cvode_stat = CVodeGetDky( cvode_mem, tfinal, 0, p_nv);
            CVODECHKERR(cvode_stat);
        }

        /* Check if we have to terminate early */
        if (stop_fun != NULL){
            stop_fun( tret, p_nv_dat, stop_mem, &istop );
            if ( istop == 1 ) {
                break;
            }
        }

        while (iout < n_tspan && tspan[iout] <= tret){
            cvode_stat = CVodeGetDky( cvode_mem, tspan[iout], 0, p_nv);
            CVODECHKERR( cvode_stat );
            out_fun( iout, tspan[ iout ], p_nv_dat, out_mem );
            iout++;
        }
        cvode_stat = CVodeGetDky( cvode_mem, tret, 0, p_nv);
        CVODECHKERR( cvode_stat );
        tnow = tret;
    }
    N_VDestroy( p_nv );
    SUNLinSolFree( linear_solver );
    CVodeSensFree( cvode_mem );
    CVodeFree( &cvode_mem );
    return 0;
}