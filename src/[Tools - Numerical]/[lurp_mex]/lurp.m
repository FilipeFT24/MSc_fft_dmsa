%LURP     LU factorization using rook pivoting, using efficient code
%
%    [L,U] = LURP(A) stores a "psychologically lower triangular matrix" 
%    (i.e. a product of a permutation matrix and a unit lower triangular
%    matrix) in L and a "psychologically upper triangular matrix" (i.e. a
%    product of upper triangular matrix and a permutation matrix) in U,so
%    that A = L*U. If A is rectangular, L or U will be trapezoidal.
% 
%    [L,U,P] = LURP(A) returns a unit lower triangular matrix L, a
%    "psychologically upper triangular matrix"  matrix U, and a
%    permutation matrix P so that P*A = L*U.
%
%    [L,U,P,Q] = LURP(A) returns a unit lower triangular matrix L, an
%    upper triangular matrix U, a pemutation matrix P, and a permutation
%    matrix Q so that P*A*Q = L*U.
% 
%    [L,U,p] = LURP(A,'vector') or [L,U,p,q] = LURP(A,'vector') returns
%    the permutation information as vectors instead of matrices.  That is,
%    for [L,U,p,q] = LURP(A,'vector') for example, p and q are row vectors
%    so that  A(p,q) = L*U.  Similarly, [L,U,P] = LURP(A,'matrix') or
%    [L,U,P,Q] = LURP(A,'matrix') returns permutation matrices.  This is 
%    the default behavior.
%
%    A call to LURP(A,TOL), LURP(A,'vector',TOL) or LURP(A,'matrix',TOL)
%    where TOL is a real number calculates an LU factorization of A with
%    threshold rook pivoting.  Standard rook pivoting has TOL = 1.0 and
%    is the default.
%    if TOL >= 1 do threshold pivoting on rows and columns.
%       The maximum magnitude entry of L will be <= TOL and the maximum
%       magnitude entry of D^(-1)*U will be <= TOL where D is the
%       diagonal of U (when U is triangular). If D has an entry equal to
%       zero, the corresponding row of U is all zero for finite TOL.
%       TOL = Inf forces Gaussian elimination with no pivoting. A large
%       or infinite value of TOL can produce an inaccurate factorization.
%    if TOL < -1 do threshold pivoting on rows only.
%       The maximum magnitude entry of L will be <= 1 and the maximum
%       magnitude entry of D^(-1)*U will be <= |TOL| where D is the
%       diagonal of U (when U is triangular). If D has an entry equal to
%       zero, the corresponding row of U is all zero for finite TOL.
%       TOL = -Inf forces Gaussian elimation with partial pivoting.
%    if -1 <= TOL <= 1 use TOL = 1 (standard rook pivoting)
%
%    Y = LURP(A) returns the output from the Fortran routine DGERPF. 
%
%    Example:
%    A = randn(200,200);
%    [L,U]= lurp(A);
%    norm(A - L*U,1) / norm(A,1)
%    b = randn(200,1);
%    c = L \ b;
%    x = U \ c;
%    norm(b - A*x) / norm(b)

%    This routine calls the Fortran subroutine dgerpf.F.  Dgerpf.F is
%       written so that most of the computations are done using BLAS-3
%       routines to reduce, where possible, memory access and
%       increase efficiency.
%    associated files:
%      lurp_install.m -- a function that installs the lurp package
%      lurp.m -- "help lurp" describes the use of the lurp code
%      lurp.F -- the gateway routine for lurp package
%      lurpf.F -- a "computational routine" called by lurp.F
%      dgerpf.F -- Fortran code implementing rook pivoting
%      dgerp2.F, dgerp3.F, dlaswq.F - called by dgerpf.F
%      ilaenv.F --  Fortran version of LAPACK code that returns block size.
%         Included since the Fortran version returns a good block size
%         whereas the version of ilaenv that comes with MATLAB's LAPACK
%         library may not return a good block size.
%      lurp_test.m -- quick test of the compiled lurp
%      lurp_test_SJ.m -- test of lurp using matrices from the San Jose
%         State University Sparse Matrix Database
%      lurp_test_UF.m -- test of lurp using matrices from the University
%         of Florida Sparse Matrix Collection
%      readme.txt -- describes the installation and use of the package
%  created or modified by L. Foster, 7-11-2012, copr.

