      SUBROUTINE DGERP2( M, N, A, LDA, IPIV, JPIV, INFO, TOL, RNRCSR,
     $                   RNCOMP )
*
*     This code is part of a package implementing Gaussian Elimination
*     with rook pivoting.  Written by L. Foster, Department of
*     Mathematics, San Jose State University, San Jose, CA, 95192,
*     foster@math.sjsu.edu.  November, 2005.  copr.
*
*     .. Scalar Arguments ..
****      IMPLICIT           NONE
      INTEGER            INFO, LDA, M, N
      DOUBLE PRECISION   RNCOMP, RNRCSR, TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*
*  Purpose
*  =======
*
*  DGERP2 computes an LU factorization of a general m-by-n matrix A
*  using partial rook pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U * Q
*  where P, Q are permutation matrices, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is a Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M by N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N).  Prior to
*          applying elimination to column i, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  JPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N).  Prior to
*          applying elimination to column i, column i of the
*          matrix was interchanged with column JPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
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
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*
*     .. Local Scalars ..
      LOGICAL            FCOL, SKIPRS
      INTEGER            INEW, IOLD, IP, ISRCH, J, JP
      DOUBLE PRECISION   RTOLC, RTOLR
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DRSCL, DSWAP, XERBLA
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
      ELSE IF( RNRCSR.LT.ZERO ) THEN
         INFO = -9
      ELSE IF( RNCOMP.LT.ZERO ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERP2', -INFO )
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
*
*     since Fortran does not guarantee short-circuit evaluation:
      SKIPRS = .FALSE.
      IF ( TOL .NE. ZERO ) THEN
         IF ( ( ONE / TOL ) .EQ. ZERO ) THEN
            SKIPRS = .TRUE.
         END IF
      END IF
*     SKIPRS is used to skip searching rows when 1,0 / TOL is zero
*
      DO 30 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
*        search column
         IP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
*
*        count comparisons and rows/columns searched
         RNCOMP = RNCOMP + ( M-J )
         RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*        check threshold condition
         IF( ABS( A( J, J ) ).GE.RTOLC*ABS( A( IP, J ) ) ) THEN
*           do not pivot
            IP = J
         END IF
*
         IF ( SKIPRS ) THEN
*           do not search row (so only search column J)
*           when 1.0 / TOL is zero or TOL is infinity
            JP = J
         ELSE
*
*           search row
            JP = J - 1 + IDAMAX( N-J+1, A( IP, J ), LDA )
*
*           count comparisons and rows/columns searched
            RNCOMP = RNCOMP + ( N-J ) + 1
            RNRCSR = RNRCSR + MIN( N-J, 1 )
*
*           check threshold condition
            IF( ABS( A( IP, J ) ).GE.RTOLR*ABS( A( IP, JP ) ) ) THEN
*              do not pivot
               JP = J
            END IF
         END IF
*
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
   10       CONTINUE
            IF( INEW.EQ.IOLD )
     $         GO TO 20
*
            IOLD = ISRCH
            ISRCH = INEW
*
            IF( FCOL ) THEN
*
*              search column isrch for the largest magnitude entry
*
               INEW = IDAMAX( M-J+1, A( J, ISRCH ), 1 ) + J - 1
*
*              count comparisons and rows/columns searched
               RNCOMP = RNCOMP + ( M-J ) + 1
               RNRCSR = RNRCSR + MIN( M-J, 1 )
*
*              check threshold pivoting condition
               IF( ABS( A( IOLD, ISRCH ) ).GE.RTOLC*
     $             ABS( A( INEW, ISRCH ) ) ) THEN
*                 do not pivot
                  INEW = IOLD
               END IF
               FCOL = .FALSE.
            ELSE
*
*              search row isrch for the largest magnitude entry
*
               INEW = IDAMAX( N-J, A( ISRCH, J+1 ), LDA ) + J
*
*              count comparisons and rows/columns searched
               RNCOMP = RNCOMP + N - J + 1
               RNRCSR = RNRCSR + MIN( N-J, 1 )
*
*              check threshold pivoting condition
               IF( ABS( A( ISRCH, IOLD ) ).GE.RTOLR*
     $             ABS( A( ISRCH, INEW ) ) ) THEN
*                 do not pivot
                  INEW = IOLD
               END IF
               FCOL = .TRUE.
            END IF
            GO TO 10
*
*           end of while INEW .NE. IOLD loop
*
   20       CONTINUE
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
*           switch columns and save JPIV(J)
*
            JPIV( J ) = JP
            CALL DSWAP( M, A( 1, J ), 1, A( 1, JP ), 1 )
*
         END IF
*
         IPIV( J ) = IP
         IF( A( IP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( IP.NE.J ) THEN
               CALL DSWAP( N, A( J, 1 ), LDA, A( IP, 1 ), LDA )
            END IF
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN
               CALL DRSCL( M-J, A( J, J ), A( J+1, J ), 1 )
*              DRSCL is safer than DSCAL
*
            END IF
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
*
   30 CONTINUE
*
      RETURN
*
*     End of DGERP2
*
      END
