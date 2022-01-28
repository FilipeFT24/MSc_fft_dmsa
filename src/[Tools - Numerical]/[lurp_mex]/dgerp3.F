      SUBROUTINE DGERP3( M, N, A, LDA, IPIV, JPIV, INFO, AROW, ACOL,
     $                   WORK, IWORK, FWORK, TOL, NB, RMXFLP, RNRCSR,
     $                   RNCOMP, RNEFLP )
*
*     This code is part of a package implementing Gaussian Elimination
*     with rook pivoting.  Written by L. Foster, Department of
*     Mathematics, San Jose State University, San Jose, CA, 95192,
*     foster@math.sjsu.edu.  November, 2005.  copr.
*
*     .. Scalar Arguments ..
****      IMPLICIT           NONE
      LOGICAL            FWORK
      INTEGER            INFO, LDA, M, N, NB
      DOUBLE PRECISION   RMXFLP, RNCOMP, RNEFLP, RNRCSR, TOL
*
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), ACOL( * ), AROW( * ), WORK( N, * )
*     ..
*
*  Purpose
*  =======
*
*  DGERP3 computes an LU factorization of NB columns and
*  NB rows of a general m-by-n matrix A using threshold
*  rook pivoting.
*
*  The (full) factorization has the form
*     A = P * L * U * Q
*  where P, Q are permutation matrices, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  Since DGERP3 only factors NB columns and rows the routine calculates
*    o   NB columns of L
*    o   NB rows of U
*    o   the rook pivoting permutations need to determine L and U
*  The equation A = P * L * U * Q is valid only for the first NB
*  columns and NB rows of A.  The elements of A outside of these
*  rows and columns are unchanged except that the elements are
*  permuted as indicated by P and Q.
*
*  This is a level 3  BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M by N matrix to be factored.
*          On exit,
*          o   the elements below the diagonal in the first NB columns
*              of A contain the first NB columns of L, below the
*              diagonal
*          o   the elements on and to the right of the diagonal in the
*              first NB rows of A contain the first NB rows of U, on
*              and to the right of the diagonal
*          o   the elements of A in rows NB+1 to M and columns NB+1 to
*              N are the original elements of A.  However rows and
*              columns have been moved around as indicated by P and Q
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(NB,M,N))
*          The pivot indices for row interchanges; for
*          1 <= i <= min(NB,M,N), row i of the matrix has been
*          interchanged with row IPIV(i).
*
*  JPIV    (output) INTEGER array, dimension (min(NB,M,N))
*          The pivot indices for column interchanges, for
*          1 <= i <= min(NB,M,N), column i of the matrix has been
*          interchanged with column IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  AROW    (work space) DOUBLE PRECISION array, dimension (N)
*           space to store a row of A
*
*  ACOL    (work space) DOUBLE PRECISION array, dimension (M)
*           space to store a column of A
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N,NB)
*          Space to store NB rows of U as NB columns.  The use of
*          WORK increases cache memory hits. Not used if FWORK is
*          false.
*
*  IWORK   (workspace) INTEGER array, dimension (M).
*          Space to store row interchanges (in a different format
*          than used in IPIV).  Current row i of A was originally
*          stored in row IWORK(i).  Not used if FWORK is false.
*
*  FWORK   (input) LOGICAL
*          If FWORK is true then use workspaces WORK and IWORK to
*             enhance cache memory hits.
*          If FWORK is false workspaces WORK and IWORK are not used.
*
*  TOL     (input) DOUBLE PRECISION
*          if TOL >= 1 do threshold pivoting on rows and columns.
*             Pivot only if:
*                magnitude of pivot entry < (magnitude of the largest
*                entry in the corresponding row or column ) / TOL
*             The maximum magnitude entry of L will be <= TOL and
*             the maximum magnitude entry of D^(-1)*U will be <= TOL
*             where D is the diagonal of U. If D has an entry equal
*             to zero, the corresponding row of U is all zero (for 
*             finite TOL).
*          if TOL < -1 do threshold pivoting on rows only. Use
*             standard pivoting (i.e. TOL = 1 ) on columns and pivot
*             in a row only if:
*                magnitude of pivot entry < (magnitude of the largest
*                entry in the corresponding row ) / | TOL |
*             The maximum magnitude entry of L will be <= 1 and
*             the maximum magnitude entry of D^(-1)*U will be <= |TOL|
*             where D is the diagonal of U. If D has an entry equal
*             to zero, the corresponding row of U is all zero (for
*             finite TOL).
*          if -1 <= TOL <= 1 use TOL = 1 (standard rook pivoting)
*
*  NB      (input/output) INTEGER
*          On entry, the number of columns and rows to attempt to factor
*          On exit, if the block factorization is requiring too many
*             flops (when RNEFLP > RMXFLP) then the block size is
*             reduced. The return value is the actual number of columns
*             factored.  This is the same as the entry value of NB when
*             RNEFLP <= RMXFLP.
*
*  RMXFLP   (input) DOUBLE PRECISION
*            The maximum number of  extra flops required by the block
*            algorithm beyond those in the unblocked algorithm.  Stop
*            the factoriation when RNEFLP > RMXFLP.
*
*  RNRCSR  (input / output) DOUBLE PRECISION
*          On entry, the number of rows and columns searched prior to
*          this call.  On exit, the number of rows and columns searched
*          at the completion of this call to the subroutine.
*
*  RNCOMP  (input / output) DOUBLE PRECISION
*          On entry, the number of comparisons in the algorithm prior
*          to this call. On exit, the number of comparisons at the
*          completion of this call to the subroutine.
*
*  RNEFLP  (input / output) DOUBLE PRECISION
*          On entry, the number of extra flops required by the block
*          algorithm (beyond those required by the unblocked algorithm)
*           prior to this call.
*          On exit, the number of extra flops required by the block
*          algorithm (beyond those required by the unblocked algorithm)
*          at the completion of the subroutine call.
*
*  =====================================================================
*
*     .. Parameters ..
*
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FCOL, SKIPRS
      INTEGER            I, INEW, IOLD, IP, IPT, ISRCH, ISRCHT, ISW
      INTEGER            J, JMAX, JMAXB, JP, JPT, MINMN
      DOUBLE PRECISION   AT, ATNEW, ATNEWV, ATV, RTOLC, RTOLR, TEMP
      character*120 string1
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DLASWP, DRSCL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( NB.LT.0 ) THEN
         INFO = -14
      ELSE IF( RMXFLP.LT.ZERO ) THEN
         INFO = -15
      ELSE IF( RNCOMP.LT.ZERO ) THEN
         INFO = -16
      ELSE IF( RNRCSR.LT.ZERO ) THEN
         INFO = -17
      ELSE IF( RNEFLP.LT.ZERO ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERP3', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      IF( TOL.GT.ONE ) THEN
*         do threshold pivoting on rows and columns
         RTOLC = ONE / TOL
         RTOLR = ONE / TOL
      ELSE IF( TOL.LT.-ONE ) THEN
*         do threshold pivoting on rows only
         RTOLC = ONE
         RTOLR = ONE / ABS( TOL )
      ELSE
*        do standard rook pivoting 
         RTOLR = ONE
         RTOLC = ONE
      END IF
      
*     since Fortran does not guarantee short-circuit evaluation:
      SKIPRS = .FALSE.
      IF ( TOL .NE. ZERO ) THEN
         IF ( ( ONE / TOL ) .EQ. ZERO ) THEN
            SKIPRS = .TRUE.
         END IF
      END IF
*     SKIPRS is used to skip searching rows when 1,0 / TOL is zero
*
      IF( FWORK ) THEN
*        initialize IWORK
         DO 10 I = 1, M
            IWORK( I ) = I
   10    CONTINUE
      END IF
*
      MINMN = MIN( M, N )
      JMAX = MIN( MINMN, NB )
      JMAXB = JMAX
      DO 60 J = 1, JMAX
*
*        Find pivot and test for singularity.
*
*        Generate column A(J:M,J) of updated A:
*           ACOL = A(J:M,J)- L(J:M,1:J-1) * U(1:J-1, J)
*        and search the column for the largest magnitude entry
*
         IF( FWORK ) THEN
*
*           use columns of workspace WORK to store rows of U
*
*           copy column of A to ACOL
            CALL DCOPY( M, A( 1, J ), 1, ACOL, 1 )
*           switch entries in ACOL
            DO 20 IPT = 1, J - 1
               ISW = IPIV( IPT )
               IF( IPT.NE.ISW ) THEN
                  TEMP = ACOL( IPT )
                  ACOL( IPT ) = ACOL( ISW )
                  ACOL( ISW ) = TEMP
               END IF
   20       CONTINUE
            CALL DGEMV( 'N', M-J+1, J-1, -ONE, A( J, 1 ), LDA,
     $                  WORK( J, 1 ), N, ONE, ACOL( J ), 1 )
            ATV = ACOL( J )
            AT = ABS( ATV )
            IP = IDAMAX( M-J+1, ACOL( J ), 1 )
            ATNEWV = ACOL( J+IP-1 )
            ATNEW = ABS( ATNEWV )
            IP = J - 1 + IP
*
*           count comparisons and rows/columns searched
            RNCOMP = RNCOMP + M - J
            RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*           check threshold pivoting condition
            IF( AT.GE.RTOLC*ATNEW ) THEN
*              do not pivot
               IP = J
            ELSE
               AT = ATNEW
               ATV = ATNEWV
            END IF
*
         ELSE
*
*           use rows of A to store rows of U
*
*           copy column of A to ACOL
            CALL DCOPY( M-J+1, A( J, J ), 1, ACOL, 1 )
            CALL DGEMV( 'N', M-J+1, J-1, -ONE, A( J, 1 ), LDA,
     $                  A( 1, J ), 1, ONE, ACOL, 1 )
            ATV = ACOL( 1 )
            AT = ABS( ATV )
            IP = IDAMAX( M-J+1, ACOL, 1 )
            ATNEWV = ACOL( IP )
            ATNEW = ABS( ATNEWV )
            IP = J - 1 + IP
*
*           count comparisons and rows/columns searched
            RNCOMP = RNCOMP + M - J
            RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*           check threshold pivoting condition
            IF( AT.GE.RTOLC*ATNEW ) THEN
*              do not pivot
               IP = J
            ELSE
               AT = ATNEW
               ATV = ATNEWV
            END IF
*
         END IF
*
*        Generate row A(IP,J:N) of updated A:
*           AROW = A(IP,J:N) - L(IP,1:J-1) * U(1:J-1,J:N)
*
         IF( FWORK ) THEN
*
*           use columns of workspace WORK to store rows of U
*
            JPT = IWORK( IP )
            CALL DCOPY( N-J+1, A( JPT, J ), LDA, AROW, 1 )
            CALL DGEMV( 'N', N-J+1, J-1, -ONE, WORK( J, 1 ), N,
     $                  A( IP, 1 ), LDA, ONE, AROW, 1 )
         ELSE
*
*           use rows of A to store rows of U
*
            CALL DCOPY( N-J+1, A( IP, J ), LDA, AROW, 1 )
            CALL DGEMV( 'T', J-1, N-J+1, -ONE, A( 1, J ), LDA,
     $                  A( IP, 1 ), LDA, ONE, AROW, 1 )
         END IF
*
         IF ( SKIPRS ) THEN
*           do not search row (so only search column J)
*           when 1.0 / TOL is zero or TOL is infinity
            JP = 1
            AROW( JP ) = ATV
            ATNEWV = ATV
            ATNEW = AT
         ELSE
*
*           search row
            JP = IDAMAX( N-J+1, AROW, 1 )
*
*           count comparisons and rows/columns searched
            RNCOMP = RNCOMP + N - J + 1
            RNRCSR = RNRCSR + MIN( N-J, 1 )
*
            ATNEWV = AROW( JP )
            ATNEW = ABS( ATNEWV )
*
*           check threshold pivoting condition
            IF( AT.GE.RTOLR*ATNEW ) THEN
*              do not pivot (and stop search for pivot)
               JP = 1
               AROW( JP ) = ATV
*              insure that element in intersection of col and row
*              (which we be the pivot element) will be the same
*              as calculated in the row and column
*              this assignment fixes a potential computer arithmetic
*              problem when A is essentially singular
            ELSE
               AT = ATNEW
               ATV = ATNEWV
            END IF
         END IF
*
         JP = J - 1 + JP
*
*        if JP = J search done, continue only if JP .NE. J
         IF( JP.EQ.J ) THEN
*
            JPIV( J ) = J
*
         ELSE
*
            ISRCH = IP
            IOLD = J
            INEW = JP
            FCOL = .TRUE.
*
*           start of while INEW .NE. IOLD loop
*
   30       CONTINUE
            IF( INEW.EQ.IOLD )
     $         GO TO 50
*
            IOLD = ISRCH
            ISRCH = INEW
*
            IF( FCOL ) THEN
*
*              Generate column A(J:M,ISRCH) of updated A:
*                 ACOL = A(J:M,ISRCH)- L(J:M,1:J-1) * U(1:J-1, ISRCH)
*              and search the column for the largest magnitude entry
*
               IF( FWORK ) THEN
*
*                 use columns of workspace WORK to store rows of U
*
                  CALL DCOPY( M, A( 1, ISRCH ), 1, ACOL, 1 )
*                 switch entries in ACOL
                  DO 40 IPT = 1, J - 1
                     ISW = IPIV( IPT )
                     IF( IPT.NE.ISW ) THEN
                        TEMP = ACOL( IPT )
                        ACOL( IPT ) = ACOL( ISW )
                        ACOL( ISW ) = TEMP
                     END IF
   40             CONTINUE
                  CALL DGEMV( 'N', M-J+1, J-1, -ONE, A( J, 1 ), LDA,
     $                        WORK( ISRCH, 1 ), N, ONE, ACOL( J ), 1 )
*
                  INEW = IDAMAX( M-J+1, ACOL( J ), 1 )
                  ATNEWV = ACOL( J+INEW-1 )
                  ATNEW = ABS( ATNEWV )
                  INEW = INEW + J - 1
*
*                 count comparisons and rows/columns searched
                  RNCOMP = RNCOMP + M - J
                  RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*                 count the extra flops due to use of a block algorithm
                  RNEFLP = RNEFLP + 2*( M-J+1 )*( J-1 )
*
*                 check threshold pivoting condition
                  IF( AT.GE.RTOLC*ATNEW ) THEN
*                    stop search
                     INEW = IOLD
                     ACOL( IOLD ) = ATV
*                    insures that element in intersection of col and
*                    row (which we be the pivot element) will be the
*                    same as calculated in the row and column
                  ELSE
                     ATV = ATNEWV
                     AT = ATNEW
                  END IF
                  FCOL = .FALSE.
*
               ELSE
*
*                 use rows of A to store rows of U
*
                  CALL DCOPY( M-J+1, A( J, ISRCH ), 1, ACOL, 1 )
                  CALL DGEMV( 'N', M-J+1, J-1, -ONE, A( J, 1 ), LDA,
     $                        A( 1, ISRCH ), 1, ONE, ACOL, 1 )
*
                  INEW = IDAMAX( M-J+1, ACOL, 1 )
                  ATNEWV = ACOL( INEW )
                  ATNEW = ABS( ATNEWV )
                  INEW = INEW + J - 1
*
*                 count comparisons and rows/columns searched
                  RNCOMP = RNCOMP + M - J
                  RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*                 count the extra flops due to use of a block algorithm
                  RNEFLP = RNEFLP + 2*( M-J+1 )*( J-1 )
*
*                 check threshold pivoting condition
                  IF( AT.GE.RTOLC*ATNEW ) THEN
*                    stop search
                     INEW = IOLD
                     ACOL( IOLD-J+1 ) = ATV
*                    insure that element in intersection of col and
*                    row (which we be the pivot element) will be the
*                    same as calculated in the row and column
                  ELSE
                     ATV = ATNEWV
                     AT = ATNEW
                  END IF
                  FCOL = .FALSE.
*
               END IF
*
            ELSE
*
*              Generate row A(ISRCH,J:N) of updated A:
*                 AROW = A(ISRCH,J:N) - L(ISRCH,1:J-1) * U(1:J-1,J:N)
*              and search row ISRCH for the largest magnitude entry
*
               IF( FWORK ) THEN
*
*                 use columns of workspace WORK to store rows of U
*
                  ISRCHT = IWORK( ISRCH )
                  CALL DCOPY( N-J+1, A( ISRCHT, J ), LDA, AROW, 1 )
                  CALL DGEMV( 'N', N-J+1, J-1, -ONE, WORK( J, 1 ), N,
     $                        A( ISRCH, 1 ), LDA, ONE, AROW, 1 )
               ELSE
*
*                 use rows of A to store rows of U
*
                  CALL DCOPY( N-J+1, A( ISRCH, J ), LDA, AROW, 1 )
                  CALL DGEMV( 'T', J-1, N-J+1, -ONE, A( 1, J ), LDA,
     $                        A( ISRCH, 1 ), LDA, ONE, AROW, 1 )
*
               END IF
*
               INEW = IDAMAX( N-J+1, AROW, 1 )
               ATNEWV = AROW( INEW )
               ATNEW = ABS( ATNEWV )
               INEW = INEW + J - 1
*
*              count comparisons and rows/columns searched
               RNCOMP = RNCOMP + N - J
               RNRCSR = RNRCSR + MIN( N-J, 1 )
*
*              count the extra flops due to use of a block algorithm
               RNEFLP = RNEFLP + 2*( N-J+1 )*( J-1 )
*
*              check threshold pivoting condition
               IF( AT.GE.RTOLR*ATNEW ) THEN
*                 stop search
                  INEW = IOLD
                  AROW( IOLD-J+1 ) = ATV
*                 insure that element in intersection of col and row
*                 (which becomes the pivot element) will be the same
*                 as calculated in the row and column
               ELSE
                  ATV = ATNEWV
                  AT = ATNEW
               END IF
               FCOL = .TRUE.
            END IF
*
            GO TO 30
*
*           end of while INEW .NE. IOLD loop
*
   50       CONTINUE
*
            IF( FCOL ) THEN
               IP = ISRCH
               JP = IOLD
            ELSE
               IP = IOLD
               JP = ISRCH
            END IF
*           IP is a row and JP is a column index
*
*           save JPIV(j), switch columns in A, switch AROW
            JPIV( J ) = JP
            IF( FWORK ) THEN
               CALL DCOPY( M, A( 1, J ), 1, A( 1, JP ), 1 )
               CALL DSWAP( J-1, WORK( J, 1 ), N, WORK( JP, 1 ), N )
            ELSE
               CALL DSWAP( M, A( 1, J ), 1, A( 1, JP ), 1 )
            END IF
            AROW( JP-J+1 ) = AROW( 1 )
            AROW( 1 ) = ATV
*
*
         END IF
*
         IF( FWORK ) THEN
*
*           use columns of workspace WORK to store rows of U
*
*           Copy ACOL to column j
            ACOL( IP ) = ATV
            CALL DCOPY( M-J+1, ACOL( J ), 1, A( J, J ), 1 )
*
            IPIV( J ) = IP
*           update IWORK
            IF( IP.NE.J ) THEN
               IPT = IWORK( J )
               IWORK( J ) = IWORK( IP )
               IWORK( IP ) = IPT
            END IF
*
            IF( AT.NE.ZERO ) THEN
*
*              Apply the interchange to columns 1:j.
*              do  rest of columns later
*
               IF( IP.NE.J ) THEN
                  CALL DSWAP( J, A( J, 1 ), LDA, A( IP, 1 ), LDA )
               END IF
*
*              Copy AROW to row J
               CALL DCOPY( N-J+1, AROW, 1, WORK( J, J ), 1 )
*
*              Compute elements J+1:M of J-th column.
*
               IF( J.LT.M )
     $            CALL DRSCL( M-J, WORK( J, J ), A( J+1, J ), 1 )
*
            ELSE
*
*              Copy AROW to row J
               CALL DCOPY( N-J+1, AROW, 1, WORK( J, J ), 1 )
*
               IF( INFO.EQ.0 ) THEN
                  INFO = J
               END IF
            END IF
*
         ELSE
*
*           use rows of A to store rows of U
*
*           Copy ACOL to column J
            ACOL( IP-J+1 ) = ATV
            CALL DCOPY( M-J+1, ACOL, 1, A( J, J ), 1 )
            IPIV( J ) = IP
*
            IF( AT.NE.ZERO ) THEN
*
*              Apply the interchange to columns 1:N.
               IF( IP.NE.J ) THEN
                  CALL DSWAP( N, A( J, 1 ), LDA, A( IP, 1 ), LDA )
               END IF
*
*              Copy AROW to row J
               CALL DCOPY( N-J+1, AROW, 1, A( J, J ), LDA )
*
*              Compute elements J+1:M of J-th column.
*
               IF( J.LT.M )
     $            CALL DRSCL( M-J, A( J, J ), A( J+1, J ), 1 )
*
            ELSE
*
*              Copy AROW to row J
               CALL DCOPY( N-J+1, AROW, 1, A( J, J ), LDA )
*
               IF( INFO.EQ.0 ) THEN
                  INFO = J
               END IF
            END IF
*
         END IF
*
         IF( RNEFLP.GT.RMXFLP ) THEN
*
*           the block algorithm is taking too many operations
*           reduce the block size to the current number of columns
*           factored and exit the do loop early
*
            NB = J
            JMAXB = J
            GO TO 70
*
*           exit loop early
*
         END IF
*
   60 CONTINUE
*
   70 CONTINUE
*
*
      IF( FWORK ) THEN
*        apply row interchanges to columns JMAXB+1 to N
         CALL DLASWP( N-JMAXB, A( 1, JMAXB+1 ), LDA, 1, JMAXB, IPIV, 1 )
*
*        copy workspace WORK back into A
         JMAX = MIN( JMAX, NB )
         DO 80 J = 1, JMAX
            CALL DCOPY( J, WORK( J, 1 ), N, A( 1, J ), 1 )
   80    CONTINUE
         DO 90 J = JMAX + 1, N
            CALL DCOPY( NB, WORK( J, 1 ), N, A( 1, J ), 1 )
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of DGERP3
*
      END
