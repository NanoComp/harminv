      SUBROUTINE CQZHES(NM,N,A,B,MATZ,Z)
!!
!! This and the two following subroutines  
!! were modified by V.Mandelshtam (Sep.1/1998)
!! to change the real arrays by complex
!!
      implicit none
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      COMPLEX*16 A(NM,N),B(NM,N),Z(NM,N)
      REAL*8 R,S,U2,XR,RHO
      REAL*8 DSQRT,CDABS,DABS
      LOGICAL MATZ
      COMPLEX*16 DCMPLX,U1,T,X,Y
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            Z(I,J) = 0d0
    2    CONTINUE
C
         Z(I,I) = (1d0,0)
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0d0
C
         DO 20 I = L, N
            S = S + DABS(dreal(B(I,L))) + DABS(dimag(B(I,L)))
   20    CONTINUE
C
         IF (S .EQ. 0d0) GO TO 100
         RHO = 0d0
C
         DO 25 I = L, N
            B(I,L) = B(I,L) / S
            RHO = RHO + dreal(B(I,L))**2 + dimag(B(I,L))**2
   25    CONTINUE
C
         R = DSQRT(RHO)
         XR = CDABS(B(L,L))
         IF (XR .EQ. 0d0) GO TO 27
         RHO = RHO + XR * R
         U1 = -B(L,L) / XR
         B(L,L) = (R / XR + 1.0) * B(L,L)
         GO TO 28
C
   27    B(L,L) = dcmplx(R,dimag(B(L,L)))
         U1 = (-1d0,0)
C
   28    DO 50 J = L1, N
            T = (0d0,0d0)
C
            DO 30 I = L, N
               T = T + dconjg(B(I,L)) * B(I,J)
   30       CONTINUE
C
            T = T / RHO
C
            DO 40 I = L, N
               B(I,J) = B(I,J) - T * B(I,L)
   40       CONTINUE
C
            B(L,J) = dconjg(U1) * B(L,J)
   50    CONTINUE
C
         DO 80 J = 1, N
            T = (0,0d0)
C
            DO 60 I = L, N
               T = T + dconjg(B(I,L)) * A(I,J)
   60       CONTINUE
C
            T = T / RHO
C
            DO 70 I = L, N
               A(I,J) = A(I,J) - T * B(I,L)
   70       CONTINUE
C
            A(L,J) = dconjg(U1) * A(L,J)
   80    CONTINUE
C
         B(L,L) = dcmplx(R * S,0d0)
C
         DO 90 I = L1, N
            B(I,L) = (0,0d0)
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (dimag(A(N,K)) .EQ. 0d0) GO TO 105
         R = CDABS(A(N,K))
         U1 = A(N,K) / R
         A(N,K) = dcmplx(R,0d0)
C
         DO 103 J = K1, N
            A(N,J) = dconjg(U1) * A(N,J)
  103    CONTINUE
C
         B(N,N) = dconjg(U1) * B(N,N)
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = DABS(dreal(A(L,K))) + DABS(dimag(A(L,K))) 
     &           + dreal(A(L1,K))
            IF (S .EQ. 0d0) GO TO 150
            U1 = A(L,K) / S
            U2 = dreal(A(L1,K)) / S
            R = DSQRT(dreal(U1)**2+dimag(U1)**2+U2*U2)
            U1 = U1 / R
            U2 = U2 / R
            A(L,K) = dcmplx(R * S,0d0)
            A(L1,K) = dcmplx(0d0,dimag(A(L1,K)))  !  ???????
C
            DO 110 J = K1, N
               X = A(L,J)
               Y = A(L1,J)
               A(L,J) = dconjg(U1) * X + U2 * Y
               A(L1,J) = U1 * Y  - U2 * X
  110       CONTINUE
C
            XR = dreal(B(L,L))
            B(L,L) = dconjg(U1) * XR
            B(L1,L) = dcmplx(-U2 * XR,dimag(B(L1,L)))
C
            DO 120 J = L1, N
               X = B(L,J)
               Y = B(L1,J)
               B(L,J) = dconjg(U1) * X + U2 * Y
               B(L1,J) = U1 * Y  - U2 * X
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = DABS(dreal(B(L1,L1))) + DABS(dimag(B(L1,L1))) 
     &           + DABS(dreal(B(L1,L)))
            IF (S .EQ. 0d0) GO TO 150
            U1 = B(L1,L1) / S
            U2 = dreal(B(L1,L)) / S
            R = DSQRT(dreal(U1)**2+dimag(U1)**2+U2**2)
            U1 = U1 / R
            U2 = U2 / R
            B(L1,L1) = dcmplx(R * S,0d0)
            B(L1,L) = dcmplx(0d0,dimag(B(L1,L)))   ! ???????
C
            DO 130 I = 1, L
               X = B(I,L1)
               Y = B(I,L)
               B(I,L1) = dconjg(U1) * X + U2 * Y
               B(I,L) = U1 * Y - U2 * X
  130       CONTINUE
C
            DO 140 I = 1, N
               X = A(I,L1)
               Y = A(I,L)
               A(I,L1) = dconjg(U1) * X + U2 * Y
               A(I,L) = U1 * Y - U2 * X
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               X = Z(I,L1)
               Y = Z(I,L)
               Z(I,L1) = dconjg(U1) * X + U2 * Y
               Z(I,L) = U1 * Y - U2 * X
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZVAL(NM,N,A,B,EPS1,ALF,BETA,MATZ,Z,IERR)
C
      implicit none
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      COMPLEX*16 A(NM,N),B(NM,N),ALF(N),Z(NM,N)
      REAL*8 BETA(N)
      REAL*8 R,S,A2,EP,U2,ANI,BNI,B11,XR,
     X       EPSA,EPSB,EPS1,ANORM,BNORM
      REAL*8 DSQRT,CDABS,DABS
      INTEGER MAX0
      LOGICAL MATZ
      COMPLEX*16 X,Y,Z3,U1,B33,B44,B3344,A1,A33,A34,A43,A44,SH
      COMPLEX*16 CDSQRT,DCMPLX
      REAL*8 DREAL,DIMAG
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0d0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0d0
      BNORM = 0d0
C
      DO 30 I = 1, N
         ANI = 0d0
         IF (I .NE. 1) ANI = DABS(dreal(A(I,I-1)))
         BNI = 0d0
C
         DO 20 J = I, N
            ANI = ANI + DABS(dreal(A(I,J))) + DABS(dimag(A(I,J)))
            BNI = BNI + DABS(dreal(B(I,J))) + DABS(dimag(B(I,J)))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. 0d0) ANORM = 1.0
      IF (BNORM .EQ. 0d0) BNORM = 1.0
      EP = EPS1
      IF (EP .GT. 0d0) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0
   40 EP = EP / 2.0
      IF (1.0 + EP .GT. 1.0) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (DABS(dreal(A(L,LM1))) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 A(L,LM1) = dcmplx(0d0,dimag(A(L,LM1)))   !   ??????
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = CDABS(B(L,L))
      IF (B11     .EQ. 0d0) GO TO 98
      U1 = B(L,L) / B11
C
      DO 97 J = L, ENORN
         A(L,J) = dconjg(U1) * A(L,J)
         B(L,J) = dconjg(U1) * B(L,J)
   97 CONTINUE
C
      B(L,L) = dcmplx(dreal(B(L,L)),0d0)           ! ??????
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALF(EN) = A(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      B(L,L) = dcmplx(0d0,dimag(B(L,L)))           ! ??????
      S = DABS(dreal(A(L,L))) + DABS(dimag(A(L,L))) 
     &     + DABS(dreal(A(L1,L)))
      U1 = A(L,L) / S
      U2 = dreal(A(L1,L)) / S
      R = DSQRT(dreal(U1)**2+dimag(U1)**2+U2**2)
      U1 = U1 / R
      U2 = U2 / R
      A(L,L) = dcmplx(R * S,0d0)
C
      DO 110 J = L1, ENORN
         X = A(L,J)
         Y = A(L1,J)
         A(L,J) = dconjg(U1) * X + U2 * Y
         A(L1,J) = U1 * Y - U2 * X
         X = B(L,J)
         Y = B(L1,J)
         B(L,J) = dconjg(U1) * X + U2 * Y
         B(L1,J) = U1 * Y - U2 * X
  110 CONTINUE
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = B(NA,NA)
      IF (CDABS(B33) .GE. EPSB) GO TO 122
      B33 = dcmplx(EPSB,0d0)
  122 B44 = B(EN,EN)
      IF (CDABS(B44) .GE. EPSB) GO TO 124
      B44 = dcmplx(EPSB,0d0)
  124 B3344 = B33 * B44
      A33 = A(NA,NA) * B44 
      A34 = A(NA,EN) * B33 - A(NA,NA) * B(NA,EN)
      A43 = dreal(A(EN,NA)) * B44
      A44 = A(EN,EN) * B33 - dreal(A(EN,NA)) * B(NA,EN)
      SH = A44
      X = A34 * A43
      IF (X .EQ. (0d0,0d0)) GO TO 140
      Y = (A33 - SH) * 0.5d0
      Z3 = CDSQRT(Y**2+X)
      U1 = Z3
      IF (dreal(Y)*dreal(U1)+dimag(Y)*dimag(U1).GE.0d0) GO TO 125
      U1 = -U1
  125 Z3 = (SH - X / (Y+U1)) / B3344
      SH = Z3
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = dcmplx(dreal(A(EN,NA)) + dreal(A(NA,ENM2)),0d0)
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = A(L,L) / B11 - SH
      A2 = dreal(A(L1,L)) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = A(K,KM1)
         A2 = dreal(A(K1,KM1))
  170    S = DABS(dreal(A1)) + DABS(dimag(A1)) + DABS(A2)
         U1 = A1 / S
         U2 = A2 / S
         R = DSQRT(dreal(U1)**2+dimag(U1)**2+U2**2)
         U1 = U1 / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            X = A(K,J)
            Y = A(K1,J)
            A(K,J) = dconjg(U1) * X + U2 * Y
            A(K1,J) = U1 * Y - U2 * X
            X = B(K,J)
            Y = B(K1,J)
            B(K,J) = dconjg(U1) * X + U2 * Y
            B(K1,J) = U1 * Y - U2 * X
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         A(K,KM1) = dcmplx(dreal(A(K,KM1)),0d0)        !   ??????
         A(K1,KM1) = (0d0,0d0)
C     ********** ZERO B(K+1,K) **********
  240    S = DABS(dreal(B(K1,K1))) + DABS(dimag(B(K1,K1))) 
     &        + DABS(dreal(B(K1,K)))
         U1 = B(K1,K1) / S
         U2 = dreal(B(K1,K)) / S
         R = DSQRT(dreal(U1)**2+dimag(U1)**2+U2**2)
         U1 = U1 / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = dreal(A(K2,K1))
         A(K2,K1) = dconjg(U1) * XR
         A(K2,K) = dcmplx(-U2 * XR,dimag(A(K2,K)))
C
  245    DO 250 I = LOR1, K1
            X = A(I,K1)
            Y = A(I,K)
            A(I,K1) = dconjg(U1) * X + U2 * Y
            A(I,K) = U1 * Y - U2 * X
            X = B(I,K1)
            Y = B(I,K)
            B(I,K1) = dconjg(U1) * X + U2 * Y
            B(I,K) = U1 * Y - U2 * X
  250    CONTINUE
C
         B(K1,K1) = dcmplx(dreal(B(K1,K1)),0d0)
         B(K1,K) = (0d0,0d0)
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            X = Z(I,K1)
            Y = Z(I,K)
            Z(I,K1) = dconjg(U1) * X + U2 * Y
            Z(I,K) = U1 * Y - U2 * X
  255    CONTINUE
C
  260 CONTINUE
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (dimag(A(EN,NA)) .EQ. 0d0) GO TO 70
      R = CDABS(A(EN,NA))
      U1 = A(EN,NA) / R
      A(EN,NA) = dcmplx(R,0d0)
C
      DO 270 J = EN, ENORN
         A(EN,J) = dconjg(U1) * A(EN,J)
         B(EN,J) = dconjg(U1) * B(EN,J)
  270 CONTINUE
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) B(N,1) = dcmplx(EPSB,dimag(B(N,1)))
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZVEC(NM,N,A,B,ALF,BETA,Z)
C
      implicit none
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      COMPLEX*16 A(NM,N),B(NM,N),ALF(N),Z(NM,N)
      REAL*8 BETA(N),BETM,EPSB
      REAL*8 CDABS
      COMPLEX*16 Z3,ALM,R,T
      COMPLEX*16 DCMPLX
      REAL*8 DREAL,DIMAG
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (N .LE. 1) GO TO 1001
      EPSB = dreal(B(N,1))
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALM = ALF(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = (0,0d0)
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * A(I,J) - ALM * B(I,J)
               IF (J .EQ. EN) GO TO 605
               T = T * B(J,EN)
  605          R = R + T
  610       CONTINUE
C
            T = ALM * BETA(I) - BETM * ALF(I)
            IF (T .EQ. (0,0d0)) T = EPSB
            Z3 = R / T
            B(I,EN) = Z3
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               Z(I,J) = Z(I,J) + Z(I,K) * B(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** DO NOT NORMALIZE THE EIGENVECTORS ***********
C
 1001 RETURN
C     ********** LAST CARD OF CQZVEC **********
      END
c
