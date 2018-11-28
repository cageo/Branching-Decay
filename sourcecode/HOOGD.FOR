C
C     *********************
      SUBROUTINE HOOGD(T,F)
C     *********************
C     DRIVER SUBROUTINE FOR THE DE HOOG INVERSION
C     DOUBLE PRECISION VERSION
C
C     DECLARATION OF VARIABLES
C     ========================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      ALPHA  = 0.D0
      RELERR = 1.D-4
C
      BIGT   = 0.8D0*T
      ATERM  = ALPHA-(DLOG(RELERR)/(2.D0*BIGT))
      M      = 15
C
      CALL HQDD(BIGT,ATERM,M,T,F)
      RETURN
      END
C
C     *********************************
      SUBROUTINE HQDD(BIGT,ATERM,M,T,F)
C     *********************************
C     SUBROUTINE FOR NUMERICAL INVERSION OF LAPLACE TRANSFORMS
C     USING THE QUOTIENT-DIFFERENCE ALGORITHM OF DE HOOG ET AL. (1982)
C
C     IMPLEMENTED BY: C.J. NEVILLE
C                     JUNE 13, 1991
C
C     NOTES: 1. THIS IS A DOUBLE PRECISION VERSION
C            2. THIS VERSION IS DESIGNED TO INVERT ANALYTICAL LAPLACE
C               TRANSFORMED EXPRESSIONS
C            3. THIS VERSION HAS BEEN MODIFIED TO PREVENT UNDERFLOW
C
C     DECLARATION OF VARIABLES
C     ========================
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      DIMENSION D(0:40),WORK(0:40)
      DOUBLE PRECISION T,BIGT,ATERM,F,PI,FACTOR,RESULT
      DOUBLE PRECISION AR,AI,SMALL
	DOUBLE PRECISION DIMAG
C
      PI     = 4.D0*DATAN(1.D0)
      ZERO   = DCMPLX(0.D0,0.D0)
      ONE    = DCMPLX(1.D0,0.D0)
      TWO    = DCMPLX(2.D0,0.D0)
      SMALL  = 1.0D-15
      FACTOR = PI/BIGT
C
C     CHECK THAT M >= 1
C     =================
      IF(M.LT.1) THEN
         WRITE(6,100)
100      FORMAT(5X,'ERROR: NUMBER OF TERMS IN INVERSION MUST EXCEED 0')
         RETURN
      END IF
      M2=M*2
C
C     CALCULATE Z
C     ===========
      Z=DCMPLX(DCOS(T*FACTOR),DSIN(T*FACTOR))
C
C     CALCULATE THE PADE TABLE
C     ========================
      AOLD=FBAR(DCMPLX(ATERM,0.D0))/TWO
      A=FBAR(DCMPLX(ATERM,FACTOR))
C
C     INITIALIZE THE TABLE ENTRIES
C     ----------------------------
      D(0)=AOLD
      WORK(0)=ZERO
C
      AR=DABS(DBLE(AOLD))
      AI=DABS(DIMAG(AOLD))
      IF(AR.LT.SMALL.AND.AI.LT.SMALL) THEN
         AOLD=DCMPLX(SMALL,SMALL)
      END IF
C
      WORK(1)=A/AOLD
      D(1)=-WORK(1)
      AOLD=A
C
C     CALCULATE SUCCESSIVE DIAGONALS OF THE TABLE
C     -------------------------------------------
      DO 10 J=2,M2
C
C        INITIALIZE CALCULATION OF THE DIAGONAL
C        --------------------------------------
         OLD2=WORK(0)
         OLD1=WORK(1)
         A=FBAR(DCMPLX(ATERM,DBLE(J)*FACTOR))
C
C        CALCULATE NEXT TERM AND SUM OF POWER SERIES
C        -------------------------------------------
         WORK(0)=ZERO
C
         AR=DABS(DBLE(AOLD))
         AI=DABS(DIMAG(AOLD))
         IF(AR.LT.SMALL.AND.AI.LT.SMALL) THEN
            AOLD=DCMPLX(SMALL,SMALL)
         END IF
C
         WORK(1)=A/AOLD
         AOLD=A
C
C        CALCULATE DIAGONAL USING THE RHOMBUS RULES
C        ------------------------------------------
         DO 20 I=2,J
            OLD3=OLD2
            OLD2=OLD1
            OLD1=WORK(I)
C
C           QUOTIENT-DIFFERENCE ALGORITHM RULES
C           -----------------------------------
            IF((I/2)*2.EQ.I) THEN
C              I EVEN: DIFFERENCE FORM
C              -----------------------
               WORK(I)=OLD3+(WORK(I-1)-OLD2)
            ELSE
C              I ODD: QUOTIENT FORM
C              --------------------
               AR=DABS(DBLE(OLD2))
               AI=DABS(DIMAG(OLD2))
               IF(AR.LT.SMALL.AND.AI.LT.SMALL) THEN
                  OLD2=DCMPLX(SMALL,SMALL)
               END IF
C
               WORK(I)=OLD3*(WORK(I-1)/OLD2)
            END IF
20       CONTINUE
C
C        SAVE CONTINUED FRACTION COEFFICIENTS
C        ------------------------------------
         D(J)=-WORK(J)
10    CONTINUE
C
C     EVALUATE CONTINUED FRACTION
C     ===========================
C     INITIALIZE RECURRENCE RELATIONS
C     -------------------------------
      AOLD2=D(0)
      AOLD1=D(0)
      BOLD2=ONE
      BOLD1=ONE+(D(1)*Z)
C
C     USE RECURRENCE RELATIONS
C     ------------------------
      DO 30 J=2,M2
         A=AOLD1+D(J)*Z*AOLD2
         AOLD2=AOLD1
         AOLD1=A
         B=BOLD1+D(J)*Z*BOLD2
         BOLD2=BOLD1
         BOLD1=B
30    CONTINUE
C
C     RESULT OF QUOTIENT-DIFFERENCE ALGORITHM
C     =======================================
      AR=DABS(DBLE(A))
      AI=DABS(DIMAG(A))
      IF(AR.LT.SMALL.AND.AI.LT.SMALL) THEN
         A=DCMPLX(0.D0,0.D0)
      END IF
      AR=DABS(DBLE(B))
      AI=DABS(DIMAG(B))
      IF(AR.LT.SMALL.AND.AI.LT.SMALL) THEN
         RESULT=0.D0
      ELSE
         RESULT=DBLE(A/B)
      END IF
C
C     CALCULATE REQUIRED APPROXIMATE INVERSE
C     ======================================
      F=DEXP(ATERM*T)*RESULT/BIGT
      RETURN
      END
