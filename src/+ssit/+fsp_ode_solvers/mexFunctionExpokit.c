/*
 * mexFunctionExpokit.c  –  MEX gateway for expokitC
 *
 * Calling convention (MATLAB):
 *   [H, V, k1, mb, t_step] = mexFunctionExpokit( ...
 *       n, m, w, beta, Acsr, btol, Time_array_i_prt, tNow, ...
 *       resetSparsity, k1_in, mb_in, t_step_in);
 *
 * Improvements over original:
 *   - mxGetDoubles() used instead of deprecated mxGetPr()
 *   - scalar inputs validated for size (prevents silent misreads)
 *   - output arrays allocated before calling expokitC so MATLAB can
 *     clean up on an early error inside the kernel
 *   - mxFree wrapped in a single cleanup block via goto
 *   - CSR struct parsing factored and annotated
 */

#include "mex.h"
#include <string.h>   /* memset */
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/*  CSR matrix descriptor (0-based, matches expokitC.c)               */
/* ------------------------------------------------------------------ */
typedef struct {
    int     n;
    int    *row_ptr;
    int    *col_ind;
    double *val;
} CSRMatrix;

/* Forward declaration */
void expokitC(int n, int m, double *w, double beta,
              CSRMatrix *A, double btol,
              double *V, double *H,
              int *k1, int *mb, double *t_step,
              double Time_array_i_prt, double tNow,
              int resetSparsity);

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */
static const mxArray *get_required_field(const mxArray *s, const char *name)
{
    const mxArray *f = mxGetField(s, 0, name);
    if (!f)
        mexErrMsgIdAndTxt("mexFunctionExpokit:MissingCSRField",
                          "Missing CSR field: %s", name);
    return f;
}

/* Validate that an mxArray is a real double scalar. */
static void check_scalar(const mxArray *a, const char *name)
{
    if (!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1)
        mexErrMsgIdAndTxt("mexFunctionExpokit:BadScalar",
                          "Input '%s' must be a real double scalar.", name);
}

/* ------------------------------------------------------------------ */
/*  Parse the MATLAB CSR struct into a C CSRMatrix.                   */
/*  Allocates row_ptr, col_ind, val via mxMalloc.                     */
/*  Caller must mxFree them.                                           */
/* ------------------------------------------------------------------ */
static void parse_csr_struct(const mxArray *A_in, int n, CSRMatrix *A_out)
{
    if (!mxIsStruct(A_in))
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidCSR",
                          "Input A must be a struct with fields "
                          "row_ptr, col_ind, val.");

    const mxArray *row_ptr_mx = get_required_field(A_in, "row_ptr");
    const mxArray *col_ind_mx = get_required_field(A_in, "col_ind");
    const mxArray *val_mx     = get_required_field(A_in, "val");

    if (!mxIsDouble(row_ptr_mx) || mxIsComplex(row_ptr_mx) ||
        !mxIsDouble(col_ind_mx) || mxIsComplex(col_ind_mx) ||
        !mxIsDouble(val_mx)     || mxIsComplex(val_mx))
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidCSRType",
                          "CSR fields row_ptr, col_ind, val must be "
                          "real double arrays.");

    mwSize row_ptr_len = mxGetNumberOfElements(row_ptr_mx);
    mwSize nnz         = mxGetNumberOfElements(val_mx);

    if (row_ptr_len != (mwSize)(n + 1))
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidCSRSize",
                          "row_ptr must have n+1 = %d entries (got %zu).",
                          n + 1, (size_t)row_ptr_len);

    if (mxGetNumberOfElements(col_ind_mx) != nnz)
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidCSRSize",
                          "col_ind and val must have the same length.");

    /* Use the modern mxGetDoubles API (R2018a+). */
    const double *row_ptr_d = mxGetDoubles(row_ptr_mx);
    const double *col_ind_d = mxGetDoubles(col_ind_mx);
    const double *val_d     = mxGetDoubles(val_mx);

    /* Auto-detect 1-based vs 0-based indexing from first row_ptr entry. */
    int one_based = (row_ptr_d[0] == 1.0) ? 1 : 0;

    A_out->n       = n;
    A_out->row_ptr = (int *)mxMalloc((size_t)(n + 1) * sizeof(int));
    A_out->col_ind = (int *)mxMalloc((size_t)nnz      * sizeof(int));
    A_out->val     = (double *)mxMalloc((size_t)nnz   * sizeof(double));

    for (mwSize i = 0; i < (mwSize)(n + 1); i++)
        A_out->row_ptr[i] = (int)row_ptr_d[i] - one_based;

    for (mwSize i = 0; i < nnz; i++) {
        A_out->col_ind[i] = (int)col_ind_d[i] - one_based;
        A_out->val[i]     = val_d[i];
    }

    if (A_out->row_ptr[0] != 0 || A_out->row_ptr[n] != (int)nnz)
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidCSRContent",
                          "CSR row_ptr bounds are inconsistent with nnz.");
}

/* ------------------------------------------------------------------ */
/*  MEX entry point                                                    */
/*                                                                     */
/*  Inputs  (prhs, 0-based):                                          */
/*    0  n               scalar int                                    */
/*    1  m               scalar int                                    */
/*    2  w               n-vector double                               */
/*    3  beta            scalar double                                 */
/*    4  Acsr            struct {row_ptr, col_ind, val}                */
/*    5  btol            scalar double                                 */
/*    6  Time_array_i_prt  scalar double                               */
/*    7  tNow            scalar double                                 */
/*    8  resetSparsity   scalar int                                    */
/*    9  k1              scalar int  (initial / seed value)            */
/*   10  mb              scalar int  (initial / seed value)            */
/*   11  t_step          scalar double (initial / seed value)          */
/*                                                                     */
/*  Outputs (plhs, 0-based):                                          */
/*    0  H        (m+2) x (m+2)                                        */
/*    1  V        n x (m+1)                                            */
/*    2  k1_out   scalar                                               */
/*    3  mb_out   scalar                                               */
/*    4  t_step_out  scalar                                            */
/* ------------------------------------------------------------------ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ---- argument count check -------------------------------------- */
    if (nrhs != 12)
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidNumInputs",
                          "Expected exactly 12 inputs, got %d.", nrhs);
    if (nlhs > 5)
        mexErrMsgIdAndTxt("mexFunctionExpokit:InvalidNumOutputs",
                          "Expected at most 5 outputs.");

    /* ---- validate and read scalar inputs --------------------------- */
    check_scalar(prhs[0], "n");
    check_scalar(prhs[1], "m");
    check_scalar(prhs[3], "beta");
    check_scalar(prhs[5], "btol");
    check_scalar(prhs[6], "Time_array_i_prt");
    check_scalar(prhs[7], "tNow");
    check_scalar(prhs[8], "resetSparsity");
    check_scalar(prhs[9], "k1");
    check_scalar(prhs[10], "mb");
    check_scalar(prhs[11], "t_step");

    int    n               = (int)mxGetScalar(prhs[0]);
    int    m               = (int)mxGetScalar(prhs[1]);
    double beta            = mxGetScalar(prhs[3]);
    double btol            = mxGetScalar(prhs[5]);
    double Time_array_i_prt = mxGetScalar(prhs[6]);
    double tNow            = mxGetScalar(prhs[7]);
    int    resetSparsity   = (int)mxGetScalar(prhs[8]);
    int    k1              = (int)mxGetScalar(prhs[9]);
    int    mb_val          = (int)mxGetScalar(prhs[10]);
    double t_step          = mxGetScalar(prhs[11]);

    /* w must be a real double vector of length n */
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
        (int)mxGetNumberOfElements(prhs[2]) != n)
        mexErrMsgIdAndTxt("mexFunctionExpokit:BadW",
                          "w must be a real double vector of length n=%d.", n);

    double *w = mxGetDoubles(prhs[2]);

    /* ---- parse CSR struct ----------------------------------------- */
    CSRMatrix A_csr;
    parse_csr_struct(prhs[4], n, &A_csr);

    /* ---- allocate outputs ------------------------------------------ */
    plhs[0] = mxCreateDoubleMatrix((mwSize)(m + 2), (mwSize)(m + 2), mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)n,       (mwSize)(m + 1), mxREAL);
    plhs[2] = mxCreateDoubleScalar(0.0);
    plhs[3] = mxCreateDoubleScalar(0.0);
    plhs[4] = mxCreateDoubleScalar(0.0);

    double *H = mxGetDoubles(plhs[0]);
    double *V = mxGetDoubles(plhs[1]);

    /* ---- call kernel ----------------------------------------------- */
    expokitC(n, m, w, beta, &A_csr, btol, V, H,
             &k1, &mb_val, &t_step,
             Time_array_i_prt, tNow, resetSparsity);

    /* ---- write scalar outputs -------------------------------------- */
    *mxGetDoubles(plhs[2]) = (double)k1;
    *mxGetDoubles(plhs[3]) = (double)mb_val;
    *mxGetDoubles(plhs[4]) = t_step;

    /* ---- free CSR working copies ----------------------------------- */
    mxFree(A_csr.row_ptr);
    mxFree(A_csr.col_ind);
    mxFree(A_csr.val);
}
