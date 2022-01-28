      SUBROUTINE DGERPF( M, N, A, LDA, IPIV, JPIV, INFO, WORK, LWORK,
     $                   IWORK, TOL, NB, RATIO )
*
*     This code is part of a package implementing Gaussian Elimination
*     with rook pivoting.  Written by L. Foster, Department of
*     Mathematics, San Jose State University, San Jose, CA, 95192,
*     foster@math.sjsu.edu.  November, 2005. copr.
*
*     .. Scalar Arguments ..
****      IMPLICIT           NONE
      INTEGER            INFO, LDA, LWORK, M, N, NB
      DOUBLE PRECISION   RATIO, TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGERPF computes an LU factorization of a general M-by-N matrix A
*  using threshold rook pivoting.
*
*  The factorization has the form
*            A = P * L * U * Q
*  where P, Q are permutation matrices, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the Level 3 BLAS version of the algorithm.  Therefore the
*  run time of this code is only a modest amount larger than code using
*  LAPACK's partial pivoting routine DGETRF.
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
*          On entry, the M-by-N matrix to be factored.
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
*          The pivot indices;  for 1 <= i <= min(M,N).  Prior to
*          applying elimination to column i, column i of the
*          matrix was interchanged with column JPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if LWORK = -1, WORK(1) returns the optimal LWORK.
*          On exit, if LWORK is not -1:
*               WORK(1) is the number of row interchanges plus the
*                        column interchanges.
*               WORK(2) is the number of element interchanges.
*               WORK(3) is the number of rows and columns searched.
*               WORK(4) is the number of comparisons required to find
*                       the pivot elements.
*               WORK(5) is the number of extra floating point
*                       operations or flops required to carry out
*                       the block algorithm (beyond the approximately
*                       n*n*(m-n) + 2*n^3 / 3 for m >= n
*                       or m*m*(n-m) + 2*m^3 / 3 for n > m required
*                       by standard Gaussian elimination)
*               The above counts are approximate and are relatively
*               more accurate for large M and N.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= MAX(5,M+N).
*          For optimal performance, for sufficiently large matrices,
*          LWORK >= MAX(5, M + N + NB*N), where NB is the blocksize.
*          In this case a copy of a block row of U is stored column
*          wise which will reduce cache misses.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array.
*          IF LWORK >= MAX(5, M + N + NB*N) IWORK must have dimension
*               at least M.
*          If LWORK < MAX(5, M + N + NB*N) IWORK is not accessed.
*
*  TOL     (input) DOUBLE PRECISION
*          if TOL >= 1 do threshold pivoting on rows and columns.
*             Pivot only if:
*                magnitude of pivot entry < (magnitude of the largest
*                entry in the corresponding row or column ) / TOL
*             The maximum magnitude entry of L will be <= TOL and
*             the maximum magnitude entry of D^(-1)*U will be <= TOL
*             where D is the diagonal of U. If D has an entry equal
*             to zero, the corresponding row of U is all zero (for.
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
*  NB      (input/output)  If NB > 0, NB is the block size.  NB columns
*          and rows of A are factored prior to updating the lower right
*          corner of the current A matrix.
*          If NB <= 0 NB is set to the default block size (from
*          LAPACK's ILAENV).
*
*  RATIO   (input) If RATIO > 0 the blocked algorithm will switch to
*          block size 1 if the "excess" floating point operations due
*          to using a blocked algorithm becomes more than RATIO times
*          the flops required for the unblocked algorithm.
*          If RATIO <= 0 the default value 0.25 is a chosen. In
*          practice a switch to the unblocked algorithm should
*          almost never occur when using a ratio of 0.25.
*
*  =====================================================================
*
*   Algorithm overview:
*          
*          For larger matrices a blocked algorithm is used: DGERP3 is
*          called to factor a block column and block row of the matrix.
*          In DGERP3 rows and columns that need to be searched are
*          updated as needed.
*
*          After a block column and block row are factored, the lower
*          right hand corner of the matrix is updated using the BLAS-3
*          routine DGEMM.
* 
*          The above calculations are repeated on the updated lower
*          right hand corner of the matrix.
*
*          There are matrices where the blocked algorithm requires
*          O(NB * M * N * N) floating point operations (flops).  The code 
*          calculates the number of flops due to use of a blocked 
*          algorithm and if this becomes large the code swithces to
*          DGERP2, which implements unblocked rook pivoting.
*
*          In practice the fraction of floating point operations done
*          using BLAS-1 and BLAS-2 routines relative to the flops done
*          using BLAS-3 routines approaches zero as M and N become large.
*          However, DGERPF can be slower than a good implementation of
*          partial pivoting since DGERPF does more BLAS-2 calculations
*          than required in partial pivoting and DGERPF requires more
*          frequent access to elements of the matrix A.
*
*  =====================================================================
*
*   DGERPF requires the routines DGERP2 and DGERP3 described above and
*          an auxiliary routine DLASWQ that interchanges columns of a
*          matrix.
*   DGERPF also requires the following Basic Linear Algebra Routines: 
*          DCOPY, DGEMM, DGEMV, DGER, DSWAP, IDAMAX.
*   DGERPF requires the LAPACK routines DLASWP, DRSCL, ILAENV.  These
*          routines may call other LAPACK or BLAS routines.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FWORK, LQUERY
      INTEGER            I, INFO2, J, JB, JN
      DOUBLE PRECISION   RFLPUB, RM, RMXFLP, RN, RNCOMP, RNCSW, RNEFLP
      DOUBLE PRECISION   RNRCSR, RNRCSW, RNRSW, RNSWCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGERP2, DGERP3, DLASWP, DLASWQ, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( ( LWORK.LT.MAX( 5, M+N ) ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERPF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
*
*        workspace query
*
         IF( NB.LE.0 )
     $      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
         WORK( 1 ) = MAX(5, M + N + NB*N)
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     initialize counters for work required
*
      RNRCSR = ZERO
      RNCOMP = ZERO
      RNEFLP = ZERO
      RN = N
      RM = M
*     calculate flops required for unblocked algorithm
      IF( M.GE.N ) THEN
         RFLPUB = RN*( 2*RM-RN-1 ) / 2 + RN*( RN-1 )*( 3*RM-RN-1 ) / 3
      ELSE
         RFLPUB = RM*( RM-1 ) / 2 + RM*( RM-1 )*( 3*RN-RM-1 ) / 3
      END IF
      IF ( RATIO .GT. ZERO ) THEN
         RMXFLP = RATIO*RFLPUB
      ELSE
         RMXFLP = 0.25*RFLPUB
      END IF
*     RMXFLP is the maximum extra flops allowed in the block algorithm
*
*     if NB <= 0, determine the block size for this environment.
*
      IF( NB.LE.0 ) THEN
         NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      END IF
*
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N )) THEN
*
*        Use unblocked code.
*
         CALL DGERP2( M, N, A, LDA, IPIV, JPIV, INFO, TOL, RNRCSR,
     $                RNCOMP )
*
      ELSE
*
*        Use blocked code.
*
*        Determine if the workspace is large enough to store a block
*        row or U columnwise (in order to enhance cache hits)
*
         IF( LWORK.GE.M+N+NB*N ) THEN
            FWORK = .TRUE.
         ELSE
            FWORK = .FALSE.
         END IF
*
         DO 50 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor block column and block row and test for exact
*           singularity.
*
            CALL DGERP3( M-J+1, N-J+1, A( J, J ), LDA, IPIV( J ),
     $                   JPIV( J ), INFO2, WORK( 1 ), WORK( N+1 ),
     $                   WORK( M+N+1 ), IWORK, FWORK, TOL, JB, RMXFLP,
     $                   RNRCSR, RNCOMP, RNEFLP )
*
*           Adjust INFO and row pivot indices.
*
            IF( INFO.EQ.0 .AND. INFO2.GT.0 )
     $         INFO = INFO2 + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
*
*           Adjust the column pivot indices
*
            DO 20 I = J, MIN( N, J+JB-1 )
               JPIV( I ) = J - 1 + JPIV( I )
   20       CONTINUE
*
*           Apply column interchanges to row 1:J-1
*
            CALL DLASWQ( J-1, A, LDA, J, J+JB-1, JPIV, 1 )
*

*
            IF( J+JB.LE.N ) THEN
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
*
               END IF
            END IF
*
*           CHECK TO SEE IF BLOCKED ALGORITHM IS TAKING TOO MANY FLOPS
*
            IF( RNEFLP.GT.RMXFLP ) THEN
*
*              THE BLOCK ALGORITHM IS TAKING TOO MANY OPERATIONS,
*              REDUCE THE BLOCK SIZE TO ONE AND FINISH THE
*              FACTORIZATION
*
               JN = J + JB
               JB = MIN( M, N ) - JN + 1
*
               CALL DGERP2( M-JN+1, N-JN+1, A( JN, JN ), LDA,
     $                      IPIV( JN ), JPIV( JN ), INFO2, TOL, RNRCSR,
     $                      RNCOMP )
*
*              Adjust INFO and the row pivot indices.
*
               IF( INFO.EQ.0 .AND. INFO2.GT.0 )
     $            INFO = INFO2 + JN - 1
               DO 30 I = JN, MIN( M, JN+JB-1 )
                  IPIV( I ) = JN - 1 + IPIV( I )
   30          CONTINUE
*
*              Apply interchanges to columns 1:JN-1.
*
               CALL DLASWP( JN-1, A, LDA, JN, JN+JB-1, IPIV, 1 )
*
*              Adjust the column pivot indices
*
               DO 40 I = JN, MIN( N, JN+JB-1 )
                  JPIV( I ) = JN - 1 + JPIV( I )
   40          CONTINUE
*
*              Apply column interchanges to row 1:JN-1
*
               CALL DLASWQ( JN-1, A, LDA, JN, JN+JB-1, JPIV, 1 )
*
*              EXIT THE LOOP EARLY.
*
               GO TO 60
*
            END IF
*
   50    CONTINUE
*
   60    CONTINUE
*
      END IF
*
*     count the number of row and column interchanges
*
      RNRCSW = ZERO
      RNRSW = ZERO
      RNCSW = ZERO
      DO 70 I = 1, MIN( M, N )
         IF( JPIV( I ).NE.I ) THEN
            RNCSW = RNCSW + 1
         END IF
         IF( IPIV( I ).NE.I ) THEN
            RNRSW = RNRSW + 1
         END IF
   70 CONTINUE
      RNRCSW = RNRSW + RNCSW
*
*     calculate the number of elements switched
*
      RNSWCH = RN*RNRSW + RM*RNCSW
*
*     store the counts of the work required in WORK
*
      WORK( 1 ) = RNRCSW
      WORK( 2 ) = RNSWCH
      WORK( 3 ) = RNRCSR
      WORK( 4 ) = RNCOMP
      WORK( 5 ) = RNEFLP
*
      RETURN
*
*     End of DGERPF
*
      END
