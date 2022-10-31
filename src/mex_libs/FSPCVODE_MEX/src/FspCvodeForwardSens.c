#include "FspCvodeForwardSens.h"

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

/* Function to evaluate the RHS of the sensitivity CME, one parameter at a time */
int SensCMERHS1Fun( int Ns, double t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *my_data,
                    N_Vector tmp1, N_Vector tmp2 ) {
    struct mydata *f_data = ( struct mydata * ) my_data;
    f_data->mvfun( f_data->n, t, N_VGetArrayPointer( yS ), N_VGetArrayPointer( tmp1 ), f_data->mvdat );
    f_data->dmvfun( f_data->n, iS, t, N_VGetArrayPointer( y ), N_VGetArrayPointer( tmp2 ), f_data->dmvdat );
    N_VLinearSum( 1.0, tmp1, 1.0, tmp2, ySdot );
    return 0;
}

/*
This function is intended for solution and forward sensitivity analysis of ODEs systems of the form
p'(t) = A(t, theta)p(t)
p(0) = p0(theta)
where the functions to evaluate the Jacobian action and the partial derivatives wrt sensitivity parameters are given from Matlab.
*/
int FspCVodeForwardSens( int n_tspan, double t_init, double *tspan, int n_state, int n_par, double *p0, double *dp0,
                         matvecfun matvec, void *matvecdat, dmatdpvecfun dmatvec, void *dmatvecdat, outputfun out_fun,
                         void *out_mem, stoppingfun stop_fun, void *stop_mem ) {
    printf( "Calling solver\n" );
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
    #ifdef USENVECOMP
    int nthreads = omp_get_max_threads();
    p_nv = N_VNew_OpenMP(n_state, nthreads);
    #else
    p_nv = N_VNew_Serial( n_state );
    #endif


    realtype *p_nv_dat = N_VGetArrayPointer( p_nv );
    for ( int i = 0; i < n_state; ++i ) {
        p_nv_dat[ i ] = p0[ i ];
    }

    /* Initialize initial N_Vector sensitivities*/

//    N_Vector *dp_nv_array = N_VCloneVectorArray_Pthreads( n_par, p_nv );
#ifdef USENVECOMP
    N_Vector *dp_nv_array = N_VCloneVectorArray_OpenMP(n_par, p_nv);
#else
    N_Vector *dp_nv_array = N_VCloneVectorArray_Serial( n_par, p_nv );
#endif

    realtype **dp_nv_dat;
    dp_nv_dat = malloc( n_par * sizeof( realtype * ));
    for ( int i = 0; i < n_par; ++i ) {
        dp_nv_dat[ i ] = N_VGetArrayPointer( dp_nv_array[ i ] );
        for ( int j = 0; j < n_state; ++j ) {
            dp_nv_dat[ i ][ j ] = dp0[ i * n_state + j ];
        }
    }

    /* Initialize data for user-supplied functions */
    struct mydata my_data;
    my_data.n = n_state;
    my_data.mvfun = matvec;
    my_data.dmvfun = dmatvec;
    my_data.mvdat = matvecdat;
    my_data.dmvdat = dmatvecdat;

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
    cvode_stat = CVodeSetMaxConvFails( cvode_mem, 10000 );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSetMaxNonlinIters( cvode_mem, 10000 );
    CVODECHKERR( cvode_stat );

    /* Create the linear solver without preconditioning */
    linear_solver = SUNSPFGMR( p_nv, PREC_NONE, 50 );
    cvode_stat = CVSpilsSetLinearSolver( cvode_mem, linear_solver );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVSpilsSetJacTimes( cvode_mem, NULL, &CMEJacVecFun );
    CVODECHKERR( cvode_stat );

    /* Initialize cvodes forward sensitivity */
    cvode_stat = CVodeSensInit1( cvode_mem, n_par, CV_SIMULTANEOUS, &SensCMERHS1Fun, dp_nv_array );
    CVODECHKERR( cvode_stat );
    cvode_stat = CVodeSensEEtolerances( cvode_mem );
    CVODECHKERR( cvode_stat );

    /* Advance to the end, output vectors at times specified in tspan */
    double tret;
    double tnow = t_init;
    double tfinal = tspan[n_tspan - 1];
    int istop = 0;
    int iout = 0;
    while ( tnow < tfinal ) {
        cvode_stat = CVode( cvode_mem, tfinal, p_nv, &tret, CV_ONE_STEP );
        CVODECHKERR( cvode_stat );
        cvode_stat = CVodeGetSens( cvode_mem, &tret, dp_nv_array );
        CVODECHKERR( cvode_stat );
        if (tret > tfinal){
            tret = tfinal;
            cvode_stat = CVodeGetDky( cvode_mem, tfinal, 0, p_nv);
            CVODECHKERR(cvode_stat);
            cvode_stat = CVodeGetSensDky(cvode_mem, tret, 0, dp_nv_array);
            CVODECHKERR( cvode_stat );
        }

        /* Check if we have to terminate early */
        if ( stop_fun != NULL) {
            stop_fun( tret, p_nv_dat, stop_mem, &istop );
            if ( istop ) {
                break;
            }
        }

        while (tspan[iout] <= tret && iout < n_tspan) {
            cvode_stat = CVodeGetDky( cvode_mem, tspan[iout], 0, p_nv);
            CVODECHKERR(cvode_stat);
            cvode_stat = CVodeGetSensDky(cvode_mem, tspan[iout], 0, dp_nv_array);
            CVODECHKERR( cvode_stat );
            out_fun( iout, tspan[ iout ], p_nv_dat, dp_nv_dat, out_mem );
            iout ++;
        }

        tnow = tret;
        cvode_stat = CVodeGetDky( cvode_mem, tret, 0, p_nv);
        CVODECHKERR(cvode_stat);
        cvode_stat = CVodeGetSensDky(cvode_mem, tret, 0, dp_nv_array);
        CVODECHKERR( cvode_stat );
    }

    free( dp_nv_dat );
    N_VDestroy( p_nv );
    N_VDestroyVectorArray( dp_nv_array, n_par );
    SUNLinSolFree( linear_solver );
    CVodeSensFree( cvode_mem );
    CVodeFree( &cvode_mem );
    return 0;
}

int FspCVodeForwardSensFromTimeZero( int n_tspan, double *tspan, int n_state, int n_par, double *p0, double *dp0,
                                     matvecfun matvec, void *matvecdat, dmatdpvecfun dmatvec, void *dmatvecdat,
                                     outputfun out_fun, void *out_mem, stoppingfun stop_fun, void *stop_mem ) {
    return FspCVodeForwardSens(n_tspan, 0.0, tspan, n_state, n_par, p0, dp0, matvec,matvecdat, dmatvec, dmatvecdat, out_fun, out_mem, stop_fun, stop_mem);
}

