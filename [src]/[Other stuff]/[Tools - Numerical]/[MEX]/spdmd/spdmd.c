/* File:  spdmd.c                                                                   */
/* Compile:  mex spdmd.c                                                            */
/* Syntax  C = spdmd(D1,M,D2)                                                       */
/* Does  C = D1 * M * D2                                                            */
/* where M  = double real sparse NxN matrix                                         */
/*       D1 = double real N element full vector representing diagonal NxN matrix    */
/*       D2 = double real N element full vector representing diagonal NxN matrix    */
/*       C  = double real sparse NxN matrix                                         */
/* Programmer:  James Tursa                                                         */
/* Date: 2/24/2021                                                                  */

/* Includes ----------------------------------------------------------- */

#include "mex.h"
#include "string.h"

/* Gateway ------------------------------------------------------------ */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSize m, n, j, nrow;
    double *Mpr, *D1pr, *D2pr, *Cpr;
    mwIndex *Mir, *Mjc, *Cir, *Cjc;

/* Argument checks */
    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many outputs");
    }
    if( nrhs != 3 ) {
        mexErrMsgTxt("Need exactly three inputs");
    }
    if (!mxIsDouble(prhs[1]) || !mxIsSparse(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgTxt("2nd argument must be real double sparse matrix");
    }
    if( !mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) ||
	mxGetNumberOfDimensions(prhs[0]) != 2 || (mxGetM(prhs[0]) != 1 && mxGetN(prhs[0]) != 1)) {
        mexErrMsgTxt("1st argument must be real double full vector");
    }
    if (!mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]) ||
	mxGetNumberOfDimensions(prhs[2]) != 2 || (mxGetM(prhs[2]) != 1 && mxGetN(prhs[2]) != 1)) {
	mexErrMsgTxt("3rd argument must be real double full vector");
    }
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (m != n || mxGetNumberOfElements(prhs[0]) != n || mxGetNumberOfElements(prhs[2]) != n) {
        mexErrMsgTxt("Matrix dimensions must agree.");
    }

/* Sparse info */
    Mir = mxGetIr(prhs[1]);
    Mjc = mxGetJc(prhs[1]);

/* Create output */
    plhs[0] = mxCreateSparse( m, n, Mjc[n], mxREAL);
    
/* Get data pointers */
    Mpr = (double *) mxGetData(prhs[1]);
    D1pr = (double *) mxGetData(prhs[0]);
    D2pr = (double *) mxGetData(prhs[2]);
    Cpr = (double *) mxGetData(plhs[0]);
    Cir = mxGetIr(plhs[0]);
    Cjc = mxGetJc(plhs[0]);

/* Fill in sparse indexing */
    memcpy(Cjc, Mjc, (n+1) * sizeof(mwIndex));
    memcpy(Cir, Mir, Cjc[n] * sizeof(mwIndex));

/* Calculate result */
    for( j=0; j<n; j++ ) {
        nrow = Mjc[j+1] - Mjc[j];  /* Number of row elements for this column */
        while( nrow-- ) {
            *Cpr++ = *Mpr++ * (D2pr[j] * D1pr[*Cir++]);  
	}
    }
}