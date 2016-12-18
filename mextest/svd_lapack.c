#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex m, n, lwork, info=0;
    double *A, *U, *S, *VT, *work;
    double workopt = 0;
    mxArray *in;

    /* verify input/output arguments */
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input matrix must be real double matrix.");
    }

    /* duplicate input matrix (since its contents will be overwritten) */
    in = mxDuplicateArray(prhs[0]);

    /* dimensions of input matrix */
    m = mxGetM(in);
    n = mxGetN(in);

    /* create output matrices */
    plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((m<n)?m:n, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);

    /* get pointers to data */
    A = mxGetPr(in);
    U = mxGetPr(plhs[0]);
    S = mxGetPr(plhs[1]);
    VT = mxGetPr(plhs[2]);

    /* query and allocate the optimal workspace size */
    lwork = -1;
    dgesvd("A", "A", &m, &n, A, &m, S, U, &m, VT, &n, &workopt, &lwork, &info);
    lwork = (mwSignedIndex) workopt;
    work = (double *) mxMalloc(lwork * sizeof(double));

    /* perform SVD decomposition */
    dgesvd("A", "A", &m, &n, A, &m, S, U, &m, VT, &n, work, &lwork, &info);

    /* cleanup */
    mxFree(work);
    mxDestroyArray(in);

    /* check if call was successful */
    if (info < 0) {
        mexErrMsgTxt("Illegal values in arguments.");
    } else if (info > 0) {
        mexErrMsgTxt("Failed to converge.");
    }
}