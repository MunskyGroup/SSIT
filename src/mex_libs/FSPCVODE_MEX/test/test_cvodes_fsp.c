#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<FspCvodeForwardSens.c>

const int n_max = 10;
const double alpha = 5.0;

void matvec( int n_state, double t, double* x, double* y, void* data){
    for ( int i = 0; i < n_state; ++i ) {
        y[i] = -alpha*x[i];
    }
}

void dmatvec(int n_state, int i_par, double t, double* x, double* y, void* data){
    for ( int i = 0; i < n_state; ++i ) {
        y[i] = -1.0*x[i];
    }
}

void out_fun(int i_tspan, double t, double* p, double** dp, void* out_mem){
    for ( int i = 0; i < n_max; ++i ) {
        printf("p[%d] = %.2e, true value %.2e, dp[%d] = %.2e true value %.2e \n", i, p[i], exp(-alpha*2.0), i, dp[0][i], -2.0*exp(-alpha*2.0));
    }
    printf("---------\n");
}

int main(){
    printf("Test.\n");
    double *p0 = malloc(n_max*sizeof(double));
    double *dp0 = malloc(n_max*sizeof(double));
    double *tspan = malloc(sizeof(double));

    for ( int i = 0; i < n_max; ++i ) {
        p0[i] = 1.0;
        dp0[i] = 0.0;
    }

    tspan[0] = 2.0;

    FspCVodeForwardSens( 1, 0.0, tspan, n_max, 1, p0, dp0, &matvec, NULL, &dmatvec, NULL, &out_fun, NULL, NULL, NULL);

    free(tspan);
    free(p0);
    free(dp0);
    return 0;
}