/*****************************************************************************/
/* sqrttimes.c                                                               */
/* Syntax: sqrttimes(A,B)                                                    */
/* Does the following:                                                       */
/* A = sqrt(A)       in-place first, then                                    */
/* B = diag(A) * B   in-place second.                                        */
/* Where A = a double vector of size 1 x M or M x 1                          */
/*       B = a double matrix of size M x N                                   */
/*                                                                           */
/* The idea is to call this function as follows:                             */
/*                                                                           */
/*   >> sqrttimes(A,B);                                                      */
/*   >> C = B' * B;                                                          */
/*                                                                           */
/* The resulting C matrix will be the same as if you had done this:          */
/*                                                                           */
/*   >> C = B' * diag(A) * B                                                 */
/*                                                                           */
/* But using sqrttimes will result in practically 0 extra temporary memory   */
/* usage so it is very memory efficient. This can be important if the A and  */
/* B matrices are very large. Note that the original A and B matrices have   */
/* been changed in-place, which is necessary to avoid extra memory usage.    */
/*                                                                           */
/* Programmer: James Tursa                                                   */
/* Date: 2011-Apr-28                                                         */
/*****************************************************************************/
#include <math.h>
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Apr, *Bpr, *Apr0;
  mwSize i, j, k, l, m, n;
/* Check inputs */
  if( nrhs != 2 ) {
    mexErrMsgTxt("Need exactly two inputs.");
  }
  if( nlhs != 0 ) {
    mexErrMsgTxt("Too many outputs.");
  }
  if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || mxIsSparse(prhs[0]) || mxIsSparse(prhs[1]) ) {
    mexErrMsgTxt("Inputs must be double non-sparse.");
  }
  if( mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2 ) {
    mexErrMsgTxt("Inputs must be 2D.");
  }
  m = mxGetM(prhs[0]);
  k = mxGetN(prhs[0]);
  l = mxGetM(prhs[1]);
  n = mxGetN(prhs[1]);
  if( m != 1 && k != 1 ) {
    mexErrMsgTxt("1st Input must be a vector.");
  }
  m = m * k;
  if( m != l ) {
    mexErrMsgTxt("1st vector input must have same number of elements as 1st dim of 2nd input.");
  }
/* Get pointers to the data */
  Apr = Apr0 = mxGetPr(prhs[0]);
  Bpr = mxGetPr(prhs[1]);
/* Do A = sqrt(A) in-place */
  for( i=0; i<m; i++ ) {
    *Apr = sqrt(*Apr);
    Apr++;
  }
/* Do B = diag(A)*B in-place */
  for( j=0; j<n; j++ ) {
    Apr = Apr0;
    for( i=0; i<m; i++ ) {
      *Bpr++ *= *Apr++;
    }
  }
}