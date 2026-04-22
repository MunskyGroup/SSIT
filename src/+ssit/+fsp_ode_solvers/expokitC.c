/*
 * expokitC.c  –  Arnoldi-based Krylov subspace approximation (Expokit style)
 *
 * Improvements over original:
 *   - restrict-qualified pointers throughout to allow compiler auto-vectorisation
 *   - csr_matvec uses a local accumulator (avoids repeated memory writes)
 *   - manual zeroing replaced with memset
 *   - scalar division replaced with a single reciprocal multiply
 *   - V(:,j+1) = p/s  done with cblas_dscal after memcpy (one BLAS call)
 *   - cross-platform BLAS header (Apple Accelerate / OpenBLAS / MKL)
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>   /* memset, memcpy */

#ifdef __APPLE__
#  include <Accelerate/Accelerate.h>
#elif defined(MKL_BLAS)
#  include <mkl_cblas.h>
#else
#  include <cblas.h>   /* OpenBLAS or any other CBLAS implementation */
#endif

/* ------------------------------------------------------------------ */
/*  CSR matrix (0-based indices)                                       */
/* ------------------------------------------------------------------ */
typedef struct {
    int     n;
    int    *row_ptr;
    int    *col_ind;
    double *val;
} CSRMatrix;

/* ------------------------------------------------------------------ */
/*  Sparse matrix-vector product  y = A * x                           */
/*  All pointers are restrict: no aliasing between A's arrays and x/y */
/* ------------------------------------------------------------------ */
static void csr_matvec(const CSRMatrix * restrict A,
                       const double    * restrict x,
                       double          * restrict y)
{
    const int     n       = A->n;
    const int    *row_ptr = A->row_ptr;
    const int    *col_ind = A->col_ind;
    const double *val     = A->val;

    for (int i = 0; i < n; i++) {
        double acc = 0.0;                       /* local accumulator  */
        const int row_end = row_ptr[i + 1];
        for (int k = row_ptr[i]; k < row_end; k++) {
            acc += val[k] * x[col_ind[k]];
        }
        y[i] = acc;
    }
}

/* ------------------------------------------------------------------ */
/*  Main Arnoldi / Expokit kernel                                      */
/* ------------------------------------------------------------------ */
void expokitC(
    int      n,
    int      m,
    double  *w,                 /* in/out: initial vector (length n)   */
    double   beta,              /* norm of w                           */
    CSRMatrix *A,               /* sparse matrix in CSR format         */
    double   btol,              /* breakdown tolerance                 */
    double  *V,                 /* out: n x (m+1), column-major        */
    double  *H,                 /* out: (m+2) x (m+2), column-major    */
    int     *k1,
    int     *mb,
    double  *t_step,
    double   Time_array_i_prt,
    double   tNow,
    int      resetSparsity)
{
    /* ---- zero output arrays ---------------------------------------- */
    memset(V, 0, (size_t)n * (size_t)(m + 1) * sizeof(double));
    memset(H, 0, (size_t)(m + 2) * (size_t)(m + 2) * sizeof(double));

    /* ---- optional sparsity enforcement on w ------------------------ */
    if (resetSparsity) {
        for (int i = 0; i < n; i++) {
            if (w[i] < 1e-10) w[i] = 0.0;
        }
    }

    /* ---- V(:,0) = w / beta ----------------------------------------- */
    const double inv_beta = 1.0 / beta;
    for (int i = 0; i < n; i++)
        V[i] = w[i] * inv_beta;

    /* ---- working vector -------------------------------------------- */
    double *p = (double *)malloc((size_t)n * sizeof(double));
    if (!p) return;   /* allocation failure; caller sees unchanged outputs */

    const int ldH = m + 2;   /* leading dimension of H                */

    /* ---- Arnoldi loop ---------------------------------------------- */
    for (int j = 0; j < m; j++) {

        /* p = A * V(:,j) */
        csr_matvec(A, &V[j * n], p);

        /* Set values of p whose absolute value is less than 1e-10 to zero */
        // for (int i = 0; i < n; i++) {
        //     if (fabs(p[i]) < 1e-10) p[i] = 0.0;
        // }

        /* Modified Gram-Schmidt orthogonalisation */
        for (int i = 0; i <= j; i++) {
            double hij = cblas_ddot(n, p, 1, &V[i * n], 1);
            H[i + j * ldH] = hij;
            cblas_daxpy(n, -hij, &V[i * n], 1, p, 1);
        }

        double s = cblas_dnrm2(n, p, 1);

        /* Invariant subspace detected – early exit */
        if (__builtin_expect(s < btol, 0)) {
            *k1     = 0;
            *mb     = j + 1;
            *t_step = Time_array_i_prt - tNow;
            free(p);
            return;
        }

        H[(j + 1) + j * ldH] = s;

        /* V(:,j+1) = p / s
         * Copy p into the next column of V, then scale in-place.
         * One memcpy + one cblas_dscal is typically faster than a
         * scalar loop for large n.                                   */
        memcpy(&V[(j + 1) * n], p, (size_t)n * sizeof(double));
        cblas_dscal(n, 1.0 / s, &V[(j + 1) * n], 1);
    }

    free(p);
}
