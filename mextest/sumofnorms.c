#include "mex.h"  /* Always include this */

#if defined(MKL_ILP64)
    #include "mkl_blas.h"
    #include "mkl_lapack.h"
    #define ptrdiff_t MKL_INT
    #define indextype long long int
#else 
    #include "blas.h"
    #include "lapack.h"
    #define indextype mwSignedIndex
#endif

#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[],      /* Output variables */
          int nrhs, const mxArray *prhs[]) /* Input variables  */
{
    /* max_error(X, T, E, Uc, Ur) */
    
    #define X_IN prhs[0]
    /*#define T_IN prhs[1]
    #define E_IN prhs[2]
    #define Uc_IN prhs[3]
    #define Ur_IN prhs[4]*/
    
    /* Loop index and accumulator */
    int i;
    double k = 0;
    /* double* X, T, E, Uc, Ur; */
    
    /* Indicates wich norm we want to compute */
    const char *norm = "I";
    
    /* Sizes of the arrays should be consistent with that */
    const mwSize *D = mxGetDimensions(X_IN);
    
    const indextype M = D[0], N = D[1];
    
    /* Size of a page */
    const int STRIDE = M*N;
    
    /* Work space for infinity norm */
    double *work = (double *) mxMalloc(2*M * sizeof(double));
    
    /* Array and subarray that starts at index i*M*N */
    double *X, *Xi;
    X = mxGetPr(X_IN);
    /*T = mxGetPr(T_IN);
    E = mxGetPr(E_IN);
    Uc = mxGetPr(Uc_IN);
    Ur = mxGetPr(Ur_IN);*/
    
    #pragma omp parallel for shared(X, D) private(i, Xi) reduction(+: k)
    for (i = 0; i < D[2]; ++i) {
        /* [n0 + D[0]*(n1 + D[1]*n2)] */
        Xi = X + i*STRIDE;
        /* frobnorm = dlange(norm, &M, &N, Xi, &M, NULL); */
        k += dlange(norm, &M, &N, Xi, &M, work);
    }
   
	/*create space for output*/
 	 plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

	/*Get pointer to output array*/
	double* y = mxGetPr(plhs[0]);

	/*Return sum to output array*/
	y[0]=k;
    
    mxFree(work);
    
    return;
}
