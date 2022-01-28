LURP: GAUSSIAN ELIMINATION WITH ROOK PIVOTING
USING EFFICIENT FORTRAN AND MATLAB MEX CODE 
  
The LURP package has efficient Fortran code to compute
LU factorizations of dense, unsymmetric matrices using
Gaussian elimination with rook pivoting. Also the package 
includes a MATLAB mex interface so that the routines
can be called directly from MATLAB as well as MATLAB
code to install and test the package.
  
Gaussian elimination with rook pivoting produces an
LU factorization of a matrix A: PAQ = LU = LDV, where
P and Q are permutation matrices, L is unit lower
triangular or trapezoidal matrix, U is an upper triangular
or trapezoidal matrix, D is a diagonal matrix and V is an
unit upper triangular or trapezoidal matrix. The trapezoidal 
matrices appear when A is not square.  In practice rook
pivoting produces a rank revealing factorization which can
be used to construct fundamental subspaces of a matrix,
to solve systems of equations involving rank deficient 
matrices, for basis repair in optimization algorithms and
for other uses. (See Nicholas Higham's text "Accuracy and
Stability of Numerical Algorithms", 2nd edition, pp. 
159-160, 170, 188, 219-220 and its references for additional
information.)
  
The Basic Linear Algebra Subroutines (BLAS) are sets 
of linear algebra primitives that are highly optimized
for specific computer architectures.  In developing 
efficient linear algebra software it is critical to
incorporate the BLAS and, in particular the level three
BLAS (BLAS-3), as much as possible. Code using BLAS-3
routines is more efficient than code using BLAS-1 and 
BLAS-2 routines.  The Fortran routines supplied with the
LURP package use BLAS-3 routines whenever possible. 
(See Higham's text pp. 578-579 for more information 
about the BLAS.)
  
Standard rook pivoting results in a factorization with
the largest magnitude elements in L and in V no larger 
than one.  The code also includes threshold rook 
pivoting which produces factorizations with the largest
magnitude elements in L and in V no larger than a user
selected pivoting tolerance TOL. 
  
Here is comparison of speed of the LURP routine with
alternatives for factoring a 4000 by 4000 random matrix
on 3 computers (A, B and C):
 

rook pivoting, [L,U]= lurp(A) ------ A: 8.8s, B: 3.3s, C: 3.7s
threshold rook, [L,U] = lurp(A,2) -- A: 7.8s, B: 2.9s, C: 3.1s
partial pivoting, [L,U] = lu(A) ---- A: 5.8s, B: 2.5s, C: 2.0s
QR factorization, [Q,R,P] = qr(A) -- A: 69s,  B: 62s,  C: 24s
singular value, [U,D,V] = svd(A) --- A: 163s, B: 126s, C: 38s
unblocked rook pivoting code ------- A: 109s, B: 79s,  C: 28s
      
Computer A: 2 cores, Intel SU4100 processor, MATLAB 7.14
Computer B: 4 cores, Intel Xeon E5404 processor, MATLAB 7.14
Computer C: 8 cores, Intel Xeon E5620 processor, MATLAB 7.12
  
The results show that MATLAB's partial pivoting code is 30%
to 85% faster than rook pivoting for these matrices and 
computers.  Note that partial pivoting does not reliably
reveal rank.  The results also show that rook pivoting 
is faster, by more than a factor of 6, than other MATLAB
dense matrix routines that produce rank revealing 
factorizations in practice.  The unblocked algorithm 
uses BLAS-1 and BLAS-2 but not BLAS-3 routines.  It is 
much slower than the BLAS-3 rook pivoting code.

To install the package:

   download and uncompress the zip file containing the
       package
   start MATLAB and move from inside MATLAB to the folder
       containing the uncompressed files
   set the proper compiler using MATLAB's mex utility: type
       mex -setup 
       and following the instructions
   type
       lurp_install
       
Lurp_install requires either gfortran for Linux computers
or Intel Fortran on Windows computers. See lurp_install 
("help lurp_install") for more discussion of the 
installation.

The files in this package include:
  lurp_install.m -- installs the LURP package
  lurp.m -- "help lurp" describes the use of LURP
  lurp.F -- the Mex gateway routine for LURP package
  lurpf.F -- a Mex "computational routine" called by lurp.F
  dgerpf.F -- Fortran code implementing rook pivoting
  dgerp2.F, dgerp3.F, dlaswq.F - called by dgerpf.F
  ilaenv.F --  Fortran version of LAPACK code that returns
     the block size. Included since the Fortran version
     returns a good block size whereas the version of
     ILAENV that comes with MATLAB's LAPACK library may not
     return a good block size.
  lurp_test.m -- quick test of the compiled LURP
  lurp_test_SJ.m -- test of LURP using matrices from 
     the San Jose State University Sparse Matrix Database
  lurp_test_UF.m -- test of LURP using matrices from the
     University of Florida Sparse Matrix Collection
  readme.txt -- describes the installation and use of the
     package
     
Other sources of rook pivoting software include:

BCSLIB-EXT (Access Analytics, PO Box 981, Redmond, WA 98073)
    Rook pivoting for symmetric matrices.  The software is 
    proprietary.

HSL_MA74 (HSL Software Library, STFC Rutherford Appleton 
    Laboratory Harwell, Oxford, DIDCOT OX11 0PE 
    United Kingdom ) Threshold rook pivoting for dense,
    unsymmetric matrices.  The software is proprietary.
    
LUSOL (Systems Optimization Laboratory, MS & E Department,
    Stanford University, Stanford, CA 94305-4121) Threshold
    rook pivoting for sparse matrices. The software is 
    governed by the Common Public License (CPL).
    

The LURP package was created by L. Foster, Department of
Mathematics, San Jose State University.  7-31-2012, copr.
