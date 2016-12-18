#include "mex.h"  /* Always include this */
#include "lapack.h" /* To compute the Frobenius norm */
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[],      /* Output variables */
          int nrhs, const mxArray *prhs[]) /* Input variables  */
{
    #define A_IN prhs[0]
     
    mwSignedIndex M, N;
    
    int i, j;
    double k;
    
    double *A;
    
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    A = mxGetPr(prhs[0]);
    
    #pragma omp parallel for shared(A, i) private(j) reduction(+: k)
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            k += A[i + M*j];
            /* mexPrintf("%f\n", A[i + M*j]); */
        }
    }
    
    mexPrintf("Now computing the Frobenius norm\n");
    double frobnorm = 0;
    const char *norm = "F";
    frobnorm = dlange(norm, &M, &N, A, &M, NULL);
   
	/*create space for output*/
 	plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);

	/*Get pointer to output array*/
	double* y = mxGetPr(plhs[0]);
    double* y2 = mxGetPr(plhs[1]);

	/*Return sum to output array*/
	y[0] = k;
    y2[0] = frobnorm;
    
    return;
}
