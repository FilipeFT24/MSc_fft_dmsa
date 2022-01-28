#include <fintrf.h>
c
c lurp.F - Gateway function for the lurp code which calls 
c          dgerpf.f, a BLAS-3 implementation of rook pivoting
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
c This is an example of the FORTRAN code required for interfacing
c a .MEX file to MATLAB and to calling LAPACK and (indirectly) BLAS
c routines..
c
c This subroutine is the main gateway to MATLAB.  When a MEX function
c  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
c  MEX file.  
c
c Copyright 1984-2006 The MathWorks, Inc.
c $Revision: 1.9.2.3 $
c
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
c-----------------------------------------------------------------------
c     (pointer) Replace integer by integer*8 on 64-bit platforms
c
****      implicit none
      mwpointer plhs(*), prhs(*)
c
c-----------------------------------------------------------------------
c

      INTEGER*4 NLHS, NRHS
c
c-----------------------------------------------------------------------
c     (pointer) Replace integer by integer*8 on 64-bit platforms
c
      mwpointer mxcreatedoublematrix, mxgetpr
c
c-----------------------------------------------------------------------
c
      mwsize mxgetm, mxgetn
      integer mxGetString, mxIsChar
c
c-----------------------------------------------------------------------
c     (pointer) Replace integer by integer*8 on 64-bit platforms
c
      mwpointer A, L, U, P, Q, ipiv,ipiv2,jpiv,jpiv2,work,iwork,tol_pr
c
c-----------------------------------------------------------------------
c
      mwsize m, n, lda, ldl, ldu, ldp, ldq, nb, lwork, minmn
      integer ilaenv, status
      integer*4 mxIsInf
      double precision tol
      character*6 vec
      
c     useful for debugging      
      character*120 string1
      
#if defined MSWIND
c For Windows only!
c This resets the floating point exception to allow divide by zero,
c overflow and invalid numbers. 
c
   INTEGER(2) CONTROL
   CALL GETCONTROLFPQQ(CONTROL)
   CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
   CONTROL = CONTROL .OR. FPCW$INVALID
   CONTROL = CONTROL .OR. FPCW$OVERFLOW
   CALL SETCONTROLFPQQ(CONTROL)
#endif
c
c check for proper number of arguments
c
      if (nrhs .gt. 3) then
        call mexerrmsgtxt('lurp requires 1 to 3 input arguments')
      elseif (nlhs .gt. 4) then
        call mexerrmsgtxt(
     $           'lurp requires at most four output arguments')
      endif
      
c
c assign pointers to the various parameters
c
      A = mxgetpr(prhs(1))
      m = mxgetm(prhs(1))
      n = mxgetn(prhs(1))
      lda = max(m,1)
      minmn = min(m, n)
c     The second input parameters .
      if ( nrhs .eq. 1) then
c        use default values for tol and vec
         tol = 1.0
         vec = 'matrix'
      endif
      if (nrhs .eq. 2) then
         if( mxIsChar(prhs(2)) .ne. 1) then
c           if second parameter is not a character string
c           let second parameter be the rook pivoting tolerance
            tol_pr = mxgetpr(prhs(2))
            call mxCopyPtrToReal8(tol_pr,tol,1)
            vec = 'matrix'
         else
c           second parameter controls the permutation format
            tol = 1.0
            status = mxGetString(prhs(2),vec,6)
         endif
      endif
      if (nrhs .eq. 3) then
c        second parameter controls permutation format
c        and third parameter is the rook pivoting tolerance
         if( mxIsChar(prhs(2)) .ne. 1) then
            call mexErrMsgTxt(
     $         'With 3 input parameters, 2nd must be a string.')
         endif
         status = mxGetString(prhs(2),vec,6)
         tol_pr = mxgetpr(prhs(3))
         call mxCopyPtrToReal8(tol_pr,tol,1)
      endif
c
c     Choose the block size the same as block size LAPACK recommends
c        for partial pivoting (a reasonable rook pivoting choice)
      nb = ilaenv( 1, 'dgetrf', ' ', m, n, -1, -1 )
c     allocate space for dgerpf work spaces
      lwork = max(5, m + n + nb*n)
      work = mxgetpr( mxcreatedoublematrix(lwork,1,0 ) )
      iwork = mxgetpr( mxcreatedoublematrix(m,1,0 ) )

******* example of a write statement for debugging purposes
*       write (string1,'(i12,5g12.4)') 1, m ,n , nb, lwork, tol
*       call mexprintf( string1 //CHAR(10) )

c     allocate work space for vectors containing pivot information
      ipiv = mxgetpr( mxcreatedoublematrix(minmn,1,0))
      ipiv2 = mxgetpr( mxcreatedoublematrix(m,1,0))
      jpiv = mxgetpr( mxcreatedoublematrix(minmn,1,0))
      jpiv2 = mxgetpr( mxcreatedoublematrix(n,1,0))
c
c     allocate space for output
      if ( nlhs .le. 1 ) then
          plhs(1) = mxcreatedoublematrix(m,n,0)
          L = mxgetpr( plhs(1) )
          ldl = max(1,m)
          U = mxgetpr( mxcreatedoublematrix(1,1,0) )
          ldu = 1
      else
c        choose L and U to be the proper size
c            if m >= n:  L is m by n and U is n by n
c            if m < n:   L is m by m and U is m by n
         if ( m .ge. n) then
            plhs(1) = mxcreatedoublematrix(m,n,0)
            L = mxgetpr( plhs(1) )
            ldl = max(1,m)
            plhs(2) = mxcreatedoublematrix(n,n,0)
            U = mxgetpr( plhs(2) )
            ldu = max(1,n)
         else
            plhs(1) = mxcreatedoublematrix(m,m,0)
            L = mxgetpr( plhs(1) )
            ldl = max(1,m)
            plhs(2) = mxcreatedoublematrix(m,n,0)
            U = mxgetpr( plhs(2) )
            ldu = max(1,m)
         endif
      endif
c     allocate space for permutation output
      if ( nlhs .ge. 3 ) then
c        choose P the proper size
c            if nlhs < 3 P is not used (a dummy is allocated)
c            if nlhs = 3 and
c                if vec is 'vector': P is 1 by m
c                if vec is not 'vector': P is m by m
         if ( vec .eq. 'vector' ) then
            plhs(3) = mxcreatedoublematrix(1,m,0)
            P = mxgetpr( plhs(3) )
            ldp = 1
         else         
            plhs(3) = mxcreatedoublematrix(m,m,0)
            P = mxgetpr( plhs(3) )
            ldp = m
         endif
      else
c        P not used set it to a "dummy"
         P = mxgetpr(mxcreatedoublematrix(1,1,0) )
         ldp = 1
      endif
      
      if ( nlhs .eq. 4 ) then
c        choose Q the proper size
c            if nlhs < 4 Q is not used (a dummy is allocated)
c            if nlhs = 4 and
c                if vec is 'vector': Q is 1 by n
c                if vec is not 'vector': Q is n by n
         if ( vec .eq. 'vector' ) then
            plhs(4) = mxcreatedoublematrix(1,n,0)
            Q = mxgetpr( plhs(4) )
            ldq = 1
         else         
            plhs(4) = mxcreatedoublematrix(n,n,0)
            Q = mxgetpr( plhs(4) )
            ldq = n
         endif
      else
c        Q not used set it to a "dummy"
         Q = mxgetpr(mxcreatedoublematrix(1,1,0) )
         ldq = 1
      endif
      if ( mxIsInf(tol) .eq. 1 .and. tol .gt. 0) then
          call mexWarnMsgTxt(CHAR(10) //'Input tolerance for rook '//
     $       'pivoting is infinity. ' // CHAR(10) //
     $       'This forces Gaussian elimination with no pivoting. ' //
     $       CHAR(10) // 'The factorization may be incorrect. '//
     $       CHAR(10) )
      endif

c     call wrapper for rook pivoting routine dgerpf

      call lurpf(m,n,%val(A),lda,%val(L),ldl,%val(U),ldu,%val(P),
     $           ldp,%val(ipiv),%val(ipiv2),%val(Q),ldq,%val(jpiv), 
     $           %val(jpiv2),%val(work),lwork,%val(iwork),tol,
     $           nb,vec,nlhs)
c
      return
      end
