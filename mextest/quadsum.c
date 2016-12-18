/* quadsum.c - Implements the quadruple for loop necessary for computing the
 * expectation of the norm in the bayesian model. Expect speedup if correc-
 * tly multithreaded.
 * 
 * WARNING: THIS CODE IS UNSAFE AND DOES NO CHECKING ON THE TYPE OR DIMENS-
 * ONS OF THE INPUTS.
 */

#include "mex.h"  /* Always include this */
#include <omp.h>

/* Function prototype: float bench2(float *A, float*B, float*T) */
void mexFunction(int nlhs, mxArray *plhs[],
          int nrhs, const mxArray *prhs[])
{
    #define A_IN prhs[0]
    #define B_IN prhs[1]
    #define T_IN prhs[2]
    #define S_IN prhs[3]
    
    /* Dimensions of T. All the dimensions are expected to match! */
    const mwSize *DT = mxGetDimensions(T_IN);

    mwSignedIndex r = DT[0];
    
    /* Retrieve the variables */
    double *A, *B, *T, *S;
    A = mxGetPr(A_IN);
    B = mxGetPr(B_IN);
    T = mxGetPr(T_IN);
    S = mxGetPr(S_IN);
    
    /* Accumulator */
    double sum = 0.0;
    
    /* Actual quadruple loop */
    #pragma omp parallel for shared(A, B, T) reduction(+: sum)
    for (mwSignedIndex i=0; i<r; ++i) {
        for (mwSignedIndex l=0; l<r; ++l) {
            for (mwSignedIndex k=0; k<r; ++k) {
                for (mwSignedIndex s=0; s<r; ++s) {
                    sum += A[i + r*s]*B[l + r*k]*
                            (
                                T[s + r*l]*T[i + r*k] +
                                S[s + r*l + r*r*( i + r*k )]
                            );
                }
            }
        }
    }
    
    /*create space for output*/
 	plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

	/*Get pointer to output array*/
	double* y = mxGetPr(plhs[0]);

	/*Return sum to output array*/
	y[0] = sum;
}