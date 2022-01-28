c
c lurpf.F - Computational routine for the lurp code which calls 
c          dgerpf.F, a BLAS-3 implementation of rook pivoting
c associated files:
c   lurp_install.m -- a function that installs the lurp package
c   lurp.m -- "help lurp" describes the use of the lurp code
c   lurp.F -- the gateway routine for lurp package
c   lurpf.F -- a "computational routine" called by lurp.F
c   dgerpf.F -- Fortran code implementing rook pivoting
c   dgerp2.F, dgerp3.F, dlaswq.F - called by dgerpf.F
c   ilaenv.F --  Fortran version of LAPACK code that returns block size.
c      Included since the Fortran version returns a good block size
c      whereas the version of ilaenv that comes with MATLAB's LAPACK
c      library may not return a good block size.
c   lurp_test.m -- quick test the compiled lurp
c   lurp_test_SJ.m -- test of lurp using matrices from the San Jose
c       State University Sparse Matrix Database
c   lurp_test_UF.m -- test of lurp using matrices from the University
c       of Florida Sparse Matrix Collection
c   readme.txt -- describes the installation and use of the package
c  created or modified by L. Foster, 7-9-2012, copr.
c
c****************************************************************
c
c
c****************************************************************
c
      subroutine lurpf( m,n,A,lda,L,ldl,U,ldu,P,ldp,ipiv,ipiv2,Q,
     &               ldq,jpiv,jpiv2,work,lwork,iwork,tol,nb,vec,nlhs)
c
c  Purpose: computational routine that is part of mex files
c           to call dgerpf.F, a BLAS-3 implementation of rook pivoting
c
c  Purpose: computational routine that is part of mex files to call
c     dgerpf.F, a BLAS-3 implementation of rook pivoting. The routine
c     (1) makes a copy of the input matrix A so that it is not changed
c     (2) calls dgerpf
c     (3) copies the output of dgerpf into lower and upper triangular
c         matrices
c     (4) interchanges rows of L and columns of U, if necessary
c     (5) changes the pivot information from LAPACK style to Matlab style
c       
c     .. Scalar Arguments ..
****      implicit none
      double precision   tol
      character*6        vec
      integer            m,n,lda,ldl,ldu,ldp,ldq,lwork,iwork,nb
      integer*4 nlhs
c     ..
c     .. Array Arguments ..
      double precision   A( lda, * ), L( ldl, * ), U(ldu,  * )
      double precision   P(ldp, * ), Q(ldq, * ), work(*)
      integer ipiv(*), ipiv2(*), jpiv(*), jpiv2(*)
c
c  Arguments
c  =========
c
c  m       (input) integer
c          The number of rows of the matrix A.  m >= 0.
c  n       (input) integer
c          The number of columns of the matrix A.  n >= 0.
c  A       (input) double precision array, dimension (lda,n)
c          on entry, the m-by-n matrix to be factored.
c  lda     (input) integer
c          The leading dimension of the array A.  lda >= max(1,m).
c  L       (output) double precision array, dimension (ldl,*)
c          if nlhs = 1 the output of dgerpf
c          if nlhs > 1 contains L in the LU factorization of A
c  ldl     (input) integer
c          The leading dimension of the array L.  ldl >= max(1,m)
c  U       (output) double precision array, dimension (ldu,*)
c          if nlhs = 1 not used
c          if nlhs > 1 contains U in the LU factorization of A
c  ldu     (input) integer
c          The leading dimension of the array U.
c          if nlhs = 1, ldu = 1
c          if nlhs >= 2, ldu = min(m,n)
c  P       (output) double precision array, dimension (ldp,m)
c          if nlhs < 3 P is not used
c          if vec .eq. 'vector' P is a 1 by m permutation vector
c          if vec .ne. 'vector' P is a m by m permutation matrix
c  ldp     (input) integer
c          The leading dimension of the array P (1 or m).
c  ipiv    (output) integer array, dimension (min(m,n))
c          the pivot indices; for 1 <= i <= min(m,n), row i of the
c          matrix was interchanged with row ipiv(i).
c  ipiv2   (output) integer array, dimension m
c          the pivot indices in MATLAB style. The pivoting of A
c          has moved row i to row ipiv2(i)
c  Q       (output) double precision array, dimension (ldq,n)
c          if nlhs < 4 Q is not used
c          if vec .eq. 'vector' Q is a 1 by n permutation vector
c          if vec .ne. 'vector' Q is a n by n permutation matrix
c  ldq     (input) integer
c          The leading dimension of the array Q (1 or n).
c  jpiv    (output) integer array, dimension (min(m,n))
c          the pivot indices; for 1 <= j <= min(m,n), column j of the
c          matrix was interchanged with column jpiv(j).
c  jpiv2   (output) integer array, dimension n
c          the pivot indices in MATLAB style. The pivoting of A
c          has moved column j to column jpiv2(j)
c  work    (workspace/output) double precision array, dimension (lwork)
c          used only for workspace here, not for output
c  lwork   (input) integer
c          lwork is set to max(5, m + n + nb*n) in lurp.F
c  iwork   (workspace) integer array
c          integer workspace, set to size m in lurp.F
c  tol     (input) double precision
c          The threshold used in rook pivoting.  See dgerpf.F or lurp.m
c          for a description
c  nb      (input) integer
c          the block size for the BLAS-3 code,  set to the size 
c          recommended for LAPACK's dgetrf in lurp.F
c  vec     (input) character string of length 6
c          if vec .eq. 'vector' P is a 1 by m permutation vector
c          if vec .ne. 'vector' P is a m by m permutation matrix
c  nlhs    (input) integer
c          the number of left hand side parameters in MATLAB call
c          to lupp (1 to 3)
c  =======      
c     
c     local variables
c
      integer i, j, itemp, jtemp, minmn, info
      double precision ratio
****      character*120 string1
c
c     make a copy of the input matrix A so that A is unchanged
      if ( ( m .ge. n ) .or. ( nlhs .le. 1 ) ) then 
c        L is the size of A so copy A into L
         do 50 j = 1,n
            do 40 i = 1,m
               L(i,j) = A(i,j)
   40       continue
   50    continue
      else
c        U is the size of A so copy A into U
         do 70 j = 1,n
            do 60 i = 1,m
               U(i,j) = A(i,j)
   60       continue
   70    continue
      end if
      
***** example of a write statement for debugging        
*       write (string1,'(i12,6g12.4)') 2, m, n, ldl, ldu, nb, tol
*       call mexprintf( string1 //CHAR(10) )
c
c     call rook pivoting routine dgerpf:
c     forces the default value of a parameter (ratio) in dgerpf:      
      ratio = -1.0
*ddd:    For no storage of block row in block col: lwork = max(5,m+n)
*ddd:    To switch to unblocked code while factoring: ratio=.01 
*ddd:    To force unblocked code:  nb = 1
      if ( ( m .ge. n ) .or. ( nlhs .le. 1 ) ) then
         call dgerpf( m, n, L, ldl, ipiv, jpiv, info, work, lwork,
     $      iwork, tol, nb, ratio)
      else
         call dgerpf( m, n, U, ldl, ipiv, jpiv, info, work, lwork,
     $      iwork, tol, nb, ratio)
      end if
      if ( info .lt. 0 ) then
c        illegal input parameter in dgerpf (should never happen)      
         return
      end if
c
c     modify LAPACK output to be suitable for MATLAB use:
      if ( nlhs .ge. 2 ) then
c         dgerpf incorporates both the lower and upper triangular
c            matrices into its output matrix.  The following code
c            puts the lower and upper triangular parts in 
c            separate matrices
          if ( m .ge. n ) then
             do 120 j = 1,n
                do 110 i = 1,j
                   U(i,j) = L(i,j)
                   L(i,j) = 0
  110           continue
                L(j,j) = 1.0
c               note: no need to zero out the lower triangle of U
c                  since Matlab's mxcreatedoublematrix fills the 
c                  array with zeros
  120        continue
          else
             do 140 j = 1,m
c               note: no need to zero out the upper triangle of L
                L(j,j) = 1.0
                do 130 i = j+1,m
                   L(i,j) = U(i,j)
                   U(i,j) = 0
  130           continue
  140        continue
          end if
       end if
c
       if ( nlhs .eq. 2 ) then
          minmn = min(m,n)
c         Use LAPACK's dlaswp to swap the rows of L.
c         L will become "psychologically" lower triangular
          call dlaswp(minmn,L,ldl,1,minmn,ipiv,-1)
       end if
       if ( ( nlhs. eq. 2) .or. (nlhs .eq. 3) ) then
          minmn = min(m,n)
c         Use dlaswq (similar to dlaswp) to swap the columns of U.
c         U will become "psychologically" upper triangular
          call dlaswq(minmn,U,ldu,1,minmn,jpiv,-1)
       end if
c       
       if ( (nlhs .eq. 3) .or. (nlhs .eq. 4) ) then
c          format the pivot information of P
c
c          change the form of the pivot vector from LAPACK 
c          style (ipiv) to  Matlab style (ipiv2)
           minmn = min(m,n)
           do 150 i = 1,m
              ipiv2(i) = i
  150      continue
           do 160 i = 1,minmn
               itemp = ipiv2(i)
               ipiv2( i ) = ipiv2( ipiv(i) )
               ipiv2( ipiv(i) ) = itemp
  160      continue
c
           if ( vec .eq. 'vector' ) then
c             format pivot informaton as a pivot vector:
              do 170 i = 1,m
                 P(1,i) = ipiv2(i)
  170         continue
           else
c             format pivot information as a pivot matrix:
c             note that mxcreatedoublematrix initializes P to zero
              do 180 i = 1,m
                 P(i,ipiv2(i) ) = 1.0
  180         continue
           end if
      end if
c      
       if ( nlhs .eq. 4 ) then
c          format the pivot information of Q
c
c          change the form of the pivot vector from LAPACK 
c          style (jpiv) to  Matlab style (jpiv2)
           minmn = min(m,n)
           do 250 j = 1,n
              jpiv2(j) = j
  250      continue
           do 260 j = 1,minmn
               jtemp = jpiv2(j)
               jpiv2( j ) = jpiv2( jpiv(j) )
               jpiv2( jpiv(j) ) = jtemp
  260      continue
c
           if ( vec .eq. 'vector' ) then
c             format pivot informaton as a pivot vector:
              do 270 j = 1,n
                 Q(1,j) = jpiv2(j)
  270         continue
           else
c             format pivot information as a pivot matrix:
c             note that mxcreatedoublematrix initializes P to zero
              do 280 j = 1,n
                 Q(jpiv2(j), j ) = 1.0
  280         continue
           end if
      end if
c
      return
      end
