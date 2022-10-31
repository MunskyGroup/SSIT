#include <math.h>
#include <mex.h>
#include "FspCvodeForwardSens.h"

#define MXSetScalarLayout(p){\
    mxSetM(p, 1); mxSetN(p, 1);\
}

#define MXSetVecLength(p, len){\
    mxSetM(p, (mwSize) len);\
    mxSetN(p, 1);\
}

#define MXSetEmptyLayout(p){\
    mxSetM(p,0); mxSetN(p,0);\
    mxSetData(p, NULL);\
}


typedef struct f_out_ml_dat {
    mxArray *output_cell;
    mxArray *function_handle;
    mxArray *finput[3];
    mxArray **sens_vec_ptr;
    int n_par;
    int n_state;
};

typedef struct f_stop_ml_dat {
    mxArray *function_handle;
    mxArray *t;
    mxArray *p;
    mxArray *stop_cond;
    mxArray *stop_status;
    int n_state;

    mxArray* plhs[2];
    mxArray* prhs[4];

    double cput;
};

typedef struct mv_ml_dat {
    mxArray *function_handle;
    mxArray *input_args[2];
    int n_state;
};

typedef struct dmv_ml_dat {
    mxArray *function_handle_cell;
    int n_par;
    mxArray *input_args[2];
    int n_state;
};

void CreateFOutData(struct f_out_ml_dat *fdat, const mxArray *out_cell, const mxArray *fhandle, int n_par, int n_state){
    fdat->n_par = n_par;
    fdat->n_state = n_state;
    fdat->output_cell = out_cell;
    fdat->function_handle = fhandle;
    fdat->finput[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    fdat->finput[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    fdat->finput[2] = mxCreateCellMatrix(n_par, 1);
    // Initialize the pointers to sensitivity vectors
    fdat->sens_vec_ptr = mxMalloc(n_par*sizeof(mxArray*));
    for (int i = 0; i < n_par; ++i){
        fdat->sens_vec_ptr[i] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mxSetCell(fdat->finput[2], i, fdat->sens_vec_ptr[i]);
    }
}

void DestroyFOutData(struct f_out_ml_dat *fdat){
    mxDestroyArray(fdat->finput[0]);
    mxDestroyArray(fdat->finput[1]);
    mxDestroyArray(fdat->finput[2]);
    for (int i = 0; i < fdat->n_par; ++i){
        mxDestroyArray(fdat->sens_vec_ptr[i]);
    }
}

void CreateFStopData(struct f_stop_ml_dat *f_stop_dat, const mxArray *stop_cond, const mxArray *fhandle, int n_state){
    f_stop_dat->stop_cond = stop_cond;
    f_stop_dat->function_handle = fhandle;
    f_stop_dat->t = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    f_stop_dat->p = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    f_stop_dat->n_state = n_state;
    f_stop_dat->stop_status = NULL;

    f_stop_dat->prhs[0] = f_stop_dat->function_handle;
    f_stop_dat->prhs[1] = f_stop_dat->t;
    f_stop_dat->prhs[2] = f_stop_dat->p;
    f_stop_dat->prhs[3] = f_stop_dat->stop_cond;
    f_stop_dat->cput = 0.0;
}

void DestroyFStopData(struct f_stop_ml_dat *f_stop_dat){
    mxDestroyArray(f_stop_dat->t);
    mxDestroyArray(f_stop_dat->p);
}

void CreateMVData(struct mv_ml_dat *mv_dat, const mxArray *fhandle, int n_state){
    mv_dat->function_handle = fhandle;
    mv_dat->n_state = n_state;
    mv_dat->input_args[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    mv_dat->input_args[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
}

void DestroyMVData(struct mv_ml_dat *mv_dat){
    mxDestroyArray(mv_dat->input_args[0]);
    mxDestroyArray(mv_dat->input_args[1]);
}

void CreateDMVData(struct dmv_ml_dat *dmv_dat, const mxArray *fhandle_cell, int n_par, int n_state){
    dmv_dat->function_handle_cell = fhandle_cell;
    dmv_dat->n_par = n_par;
    dmv_dat->n_state = n_state;
    dmv_dat->input_args[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);   
    dmv_dat->input_args[1] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
}

void DestroyDMVData(struct dmv_ml_dat *dmv_dat){
    mxDestroyArray(dmv_dat->input_args[0]);
    mxDestroyArray(dmv_dat->input_args[1]);    
}

void mex_out_fun(int i_tspan, double t, double *p, double **dp, void* output_mem){
    struct f_out_ml_dat *f_dat = (struct f_out_ml_dat *) output_mem;
    // Put in input arrays into the arguments of the output function
    mxSetDoubles(f_dat->finput[0], &t); MXSetScalarLayout(f_dat->finput[0]);
    mxSetDoubles(f_dat->finput[1], p); MXSetVecLength(f_dat->finput[1], f_dat->n_state);
    // Set pointers to dp
    for (int i=0; i < f_dat->n_par; ++i){
        mxSetDoubles(f_dat->sens_vec_ptr[i], dp[i]);
        MXSetVecLength(f_dat->sens_vec_ptr[i], f_dat->n_state);
    }
    // Call the function handle
    int ierr;
    mxArray *lhs[1];
    mxArray *rhs[4];
    rhs[0] = f_dat->function_handle;
    rhs[1] = f_dat->finput[0];
    rhs[2] = f_dat->finput[1];
    rhs[3] = f_dat->finput[2];
    ierr = mexCallMATLAB(1, lhs, 4, rhs, "feval");
    if (ierr!=0){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Error when calling user-suplied output function.");
        exit(-1);
    }
    // Enter the result to the output cell    
    mxSetCell(f_dat->output_cell, i_tspan, lhs[0]);
    
    // Return dp to the calling process (i.e., remove references from sens_vec_ptr)
    MXSetEmptyLayout(f_dat->finput[0]);
    MXSetEmptyLayout(f_dat->finput[1]);
    for (int i=0; i < f_dat->n_par; ++i){        
        MXSetEmptyLayout(f_dat->sens_vec_ptr[i]);
    }    
}

void mex_stop_fun( double t, double *p, void *stopping_data, int *stop){
    struct f_stop_ml_dat *fstop_dat = (struct f_stop_ml_dat *) stopping_data;

    // Attach matlab function to the input arrays
    mxSetDoubles(fstop_dat->t, &t); MXSetScalarLayout(fstop_dat->t);
    mxSetDoubles(fstop_dat->p, p); MXSetVecLength(fstop_dat->p, fstop_dat->n_state);

    // Evaluate the stopping criteria
    int ierr;
    ierr = mexCallMATLAB(2, fstop_dat->plhs, 4, fstop_dat->prhs, "feval");
    if (ierr!=0){
        mexErrMsgIdAndTxt("MATLAB:FspCvodeMex",
                          "Error when calling user-suplied stop function.");
        exit(-1);
    }
    // Update the stop status and decision
    if (fstop_dat->stop_status != NULL){
        mxDestroyArray(fstop_dat->stop_status);
        fstop_dat->stop_status = NULL;
    }
    fstop_dat->stop_status = mxDuplicateArray(fstop_dat->plhs[1]);
    double *istop_db;
    if (!mxIsDouble(fstop_dat->plhs[0])){
        mexErrMsgIdAndTxt("MATLAB:FspCvodeMex",
                          "First ouput of the stopping evaluation function must be a real value of class double.");
        exit(-1);
    }
    istop_db = mxGetDoubles(fstop_dat->plhs[0]);
    *stop = (int) *istop_db;

    mxDestroyArray(fstop_dat->plhs[0]);
    mxDestroyArray(fstop_dat->plhs[1]);

    // Detach matlab function from the input arrays
    MXSetEmptyLayout(fstop_dat->t);
    MXSetEmptyLayout(fstop_dat->p);
}

void mex_matvec(int n_state, double t, double *x, double *y, void *user_data){    
    struct mv_ml_dat* mv_dat = (struct mv_ml_dat*) user_data;
    mxArray *lhs[1];
    mxArray *rhs[3];

    // Point mv_dat internal data to the input arrays
    mxSetDoubles(mv_dat->input_args[0], &t);
    mxSetDoubles(mv_dat->input_args[1], x);
    mxSetM(mv_dat->input_args[1], (mwSize) mv_dat->n_state); mxSetN(mv_dat->input_args[1], 1);

    // Form the rhs for feval
    rhs[0] = mv_dat->function_handle;
    rhs[1] = mv_dat->input_args[0];
    rhs[2] = mv_dat->input_args[1];    
    // Call the matrix action handle from MATLAB
    int ierr;
    ierr = mexCallMATLAB(1, lhs, 3, rhs, "feval");
    if (ierr!=0){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                          "Error when calling user-suplied matvec function.");
        exit(-1);
    }
    // Copy the result of matvec to the output array    
    mxDouble* out_ptr = mxGetDoubles(lhs[0]);
    for ( int i = 0; i < n_state; ++i ) {        
        y[i] = out_ptr[i];
    }
    mxDestroyArray(lhs[0]);    
    // Un-reference the mv_dat internal data from the input arrays
    mxSetDoubles(mv_dat->input_args[0], NULL);
    mxSetDoubles(mv_dat->input_args[1], NULL);
    mxSetM(mv_dat->input_args[1], 0); mxSetN(mv_dat->input_args[1], 0);
}

void mex_dmatvec(int n_state, int i_par, double t, double *x, double *y, void *user_data){        
    struct dmv_ml_dat* dmv_dat = (struct dmv_ml_dat*) user_data;
    mxArray *lhs[1];
    mxArray *rhs[3];

    // Point mv_dat internal data to the input arrays
    mxSetDoubles(dmv_dat->input_args[0], &t);
    mxSetDoubles(dmv_dat->input_args[1], x);
    mxSetM(dmv_dat->input_args[1], (mwSize) dmv_dat->n_state); mxSetN(dmv_dat->input_args[1], 1);

    // Form the rhs for feval
    rhs[0] = mxGetCell(dmv_dat->function_handle_cell, i_par);
    rhs[1] = dmv_dat->input_args[0];
    rhs[2] = dmv_dat->input_args[1];
    // Call the matrix action handle from MATLAB
    int ierr;
    ierr = mexCallMATLAB(1, lhs, 3, rhs, "feval");
    if (ierr!=0){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                          "Error when calling user-suplied matvec function.");
        exit(-1);
    }
    // Copy the result of matvec to the output array
    mxDouble* out_ptr = mxGetDoubles(lhs[0]);
    for ( int i = 0; i < n_state; ++i ) {
        y[i] = out_ptr[i];
    }
    mxDestroyArray(lhs[0]);
    // Un-reference the mv_dat internal data from the input arrays
    mxSetDoubles(dmv_dat->input_args[0], NULL);
    mxSetDoubles(dmv_dat->input_args[1], NULL);
    mxSetM(dmv_dat->input_args[1], 0); mxSetN(dmv_dat->input_args[1], 0);
}


/* Mex Function to evaluate a generic function of the CME solution and sensitivities
 * This function should be called in Matlab with the following syntax
 *       [y_out, stop_status] = FSPSensCVodeSolve(t_start, t_out, Av, dAv, p0, dp0, f_out, f_stop, stop_cond)
 * where
 *
 * t_start : initial time.
 *
 * t_out : vector of output times.
 *
 * Av : function handle to evaluate action of the (possibly time-dependent) CME matrix, callable in the form y = Av(t, x).
 *
 * dAv : cell of function handles to evaluate action of the partial derivative of the CME matrix, callable in the form y = dAv{i}(t, x)
 *where i is the index of the parameter with respect to which the partial derivative is taken.
 *
 * p0 : column vector of the initial probability vector.
 *
 * dp0 : column vector made by stacking together the initial sensitivity vectors.
 *
 * f_out : function handle to evaluate the outputs, callable in the form y = f_out(p, dp) where p is the probability vector and dp 
 * contains all sensitivity vectors at a time point.
 *
 * * f_stop : function handle to evaluate whether to stop the ODE integration before reaching final time (for example, when the FSP error exceeds tolerance).
 * This function should be callable in Matlab in the form [istop, stop_stat] = f_stop(t, p, stop_cond) where t is the current time for the ODE solver and p is the current solution. Here, istop = 0 if p(t) passes the stopping condition, and 1 otherwise. The second output contains further information (for example, sum of sink states in the FSP).
 *
 * stop_cond : an array/cell/struct compatible with f_stop, containing the information needed by f_stop to make decision.
 *
 * y_out : cell array with the same number of entries as t_out, used to store the results of f_out evaluations.
 *
 * stop_status : the status of the solver
 */
void 
mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray	*prhs[])
{
    /* Check that input arguments are correct */
    if (nlhs != 2){
        mexErrMsgIdAndTxt( "MATLAB:FSPSensCVodeSolve",
                "Not enough outputs.");
    }
    if (nrhs != 9){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Not enough inputs.");
    }
    if (!mxIsClass(prhs[0], "double")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Input t_start must be of class double.");
    }
    if (!mxIsClass(prhs[1], "double")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                          "Input tspan must be of class double.");
    }
    if (!mxIsClass(prhs[2], "function_handle")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Third argument must be a function handle.");
    }
    if (!mxIsClass(prhs[3], "cell")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Fourth argument must be a cell of function handles.");
    }
    if (!mxIsClass(prhs[4], "double")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Fifth argument must be a column vector of double type.");
    }
    if (!mxIsClass(prhs[5], "double")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Sixth argument must be a column vector of double type.");
    }
    if (!mxIsClass(prhs[6], "function_handle")){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                "Seventh argument must be a function handle.");
    }
    size_t n_tspan = mxGetNumberOfElements(prhs[1]);    
    size_t n_state = mxGetNumberOfElements(prhs[4]);
    size_t n_par = mxGetNumberOfElements(prhs[5])/n_state;
    plhs[0] = mxCreateCellMatrix(n_tspan, 1);
//     
//     printf("Evaluating for %d time points with %d sensitivity parameters and %d states.\n",
//             n_tspan, n_par, n_state);
//     
    struct f_out_ml_dat fout_dat;
    struct f_stop_ml_dat fstop_dat;
    struct mv_ml_dat mv_dat;
    struct dmv_ml_dat dmv_dat;
    
    CreateFOutData(&fout_dat, plhs[0], prhs[6], n_par, n_state);
    CreateMVData(&mv_dat, prhs[2], n_state);
    CreateDMVData(&dmv_dat, prhs[3], n_par, n_state);
    CreateFStopData(&fstop_dat, prhs[8], prhs[7], n_state);
    
    double* p0_ptr;
    double* dp0_ptr;
    double* t_start;
    double* tspan;
    p0_ptr = (double*)mxGetPr(prhs[4]);
    dp0_ptr = (double*)mxGetPr(prhs[5]);
    t_start = mxGetDoubles(prhs[0]);
    tspan = (double*)mxGetPr(prhs[1]);
    
    int ierr = FspCVodeForwardSens( n_tspan, *t_start, tspan, n_state, n_par, p0_ptr, dp0_ptr, &mex_matvec, ( void * ) &mv_dat,
                                    &mex_dmatvec, ( void * ) &dmv_dat, &mex_out_fun, ( void * ) &fout_dat,
                                    &mex_stop_fun, ( void * ) &fstop_dat );
    if (ierr != 0){
        mexErrMsgIdAndTxt("MATLAB:FSPSensCVodeSolve",
                          "CVode gives error status %d \n.", ierr);
    }
    plhs[1] = mxDuplicateArray(fstop_dat.stop_status);

    DestroyDMVData(&dmv_dat);
    DestroyMVData(&mv_dat);
    DestroyFStopData(&fstop_dat);
    DestroyFOutData(&fout_dat);
}