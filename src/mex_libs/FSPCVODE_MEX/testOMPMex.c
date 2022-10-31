#include <math.h>
#include <mex.h>
#include <omp.h>
#include <time.h>

void 
mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray	*prhs[])
{
    int nthreads;
    nthreads = omp_get_max_threads();
    mexPrintf("nthreads = %d \n", nthreads);
}