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

// #include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[],      /* Output variables */
          int nrhs, const mxArray *prhs[]) /* Input variables  */
{
    /* max_error(X, T, E, Uc, Ur)
     Compute the frobenius norm of Xn - UcTnUr' - En for each n
     Divide by ||Xn||_F
     Return the max
     */
    
    #define X_IN prhs[0]
    #define T_IN prhs[1]
    #define E_IN prhs[2]
    #define Uc_IN prhs[3]
    #define Ur_IN prhs[4]
    
    /* Duplicate input X since the content will be overwritten */
    mxArray *X_WORK = mxDuplicateArray(X_IN);
    
    /* Dimensions for X and E */
    const mwSize *DXE = mxGetDimensions(X_IN);
    /* Dimensions for T */
    const mwSize *DT = mxGetDimensions(T_IN);
    
    /* Size of a page for X, E and T */
    const indextype m = DXE[0], n = DXE[1], r = DT[0], p = DT[1];
    
    const indextype STRIDE_XE = m*n;
    const indextype STRIDE_T = r*p;
    
    /**************************************************/
    /* ACTUAL WORK */
    
    /* Indicates wich norm we want to compute */
    const char *norm = "F";
    /* Indicates if the matrix should be transposed */
    const char *trp = "T";
    const char *no_trp = "N";
    
    /* Work space for infinity norm */
    /* double *work = (double *) mxMalloc(2*M * sizeof(double)); */
    
    /* Arrays and subarrays that start at index i*M*N */
    double *X, *T, *E, *Uc, *Ur, *Xi, *Ti, *Ei;
    /* Constants for the matrix multiplication routine */
    const double one = 1.0, zero = 0.0, m_one = -1.0;
    const indextype index_one = 1;
   
    X = mxGetPr(X_WORK);
    T = mxGetPr(T_IN);
    E = mxGetPr(E_IN);
    Uc = mxGetPr(Uc_IN);
    Ur = mxGetPr(Ur_IN);
    
    double norm_X, norm_Rec, relative_error, max_error;
    mxArray *temp1 = mxCreateDoubleMatrix(m, p, mxREAL);
    mxArray *temp2 = mxCreateDoubleMatrix(m, n, mxREAL);
    double *C1 = mxGetPr(temp1);
    double *C2 = mxGetPr(temp2);
        
    /* Loop. First compute A = Uc*Tn, then B = A*Ur', then take norm of difference
     * and divide by the norm of X */ 
    max_error = -1.0;
//     #pragma omp parallel for shared(X, T, E, Uc, Ur) private(Xi, Ti, Ei) /* reduction(+: k) */
    for (int i = 0; i < DXE[2]; ++i) {
        Xi = X + i*STRIDE_XE;
        Ti = T + i*STRIDE_T;
        Ei = E + i*STRIDE_XE;
        
        norm_X = dlange(norm, &m, &n, Xi, &m, NULL);
        
        /* Uc*T then (Uc*T)*Ur' */
        dgemm (no_trp, no_trp, &m, &p, &r, &one, Uc, &m, Ti, &r, &zero, C1, &m);
        dgemm (no_trp, trp, &m, &n, &p, &one, C1, &m, Ur, &n, &zero, C2, &m);
        /* Xi <- -1*C2 + Xi */
        daxpy (&STRIDE_XE, &m_one, C2, &index_one, Xi, &index_one);
        /* Xi <- -1*Ei + Xi */
        daxpy (&STRIDE_XE, &m_one, Ei, &index_one, Xi, &index_one);
        
        /* Reconstruction error */
        norm_Rec = dlange(norm, &m, &n, Xi, &m, NULL);
        relative_error = norm_Rec / norm_X;
        
        if (relative_error > max_error)
            max_error = relative_error;
    }
    
    /**************************************************/    
   
    /* Free the memory */
    mxDestroyArray(temp1);
    mxDestroyArray(temp2);
    mxDestroyArray(X_WORK);
    
	/*create space for output*/
 	plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

	/*Get pointer to output array*/
	double* y = mxGetPr(plhs[0]);

	/*Return sum to output array*/
	y[0] = max_error;
    /*plhs[1] = X_WORK;
    plhs[2] = temp1;
    plhs[3] = temp2;*/
    
    return;
}
