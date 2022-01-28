      SUBROUTINE DLASWQ( M, A, LDA, K1, K2, JPIV, INCX )
*
*     Modification of:
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*     To do column interchanges rather than row interchanges.
*     L. Foster
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            JPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWQ performs a series of column interchanges on the matrix A.
*  One column interchange is initiated for each of columns K1 through
*  K2 of A.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,*)
*          On entry, the matrix  to which the row interchanges will
*          be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of JPIV for which a column interchange
*          will be done.
*
*  K2      (input) INTEGER
*          The last element of JPIV for which a column interchange 
*          will be done.
*
*  JPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of JPIV are accessed.
*          JPIV(K) = L implies columns K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of JPIV.  If INCX
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange column I with column JPIV(I) for each of columns
*     K1 to K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = JPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( M, A( 1, I ), 1, A( 1, IP ), 1 )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = JPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( M, A( 1, I ), 1, A( 1, IP ), 1 )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = JPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( M, A( 1, I ), 1, A( 1, IP ), 1 )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWQ
*
      END
