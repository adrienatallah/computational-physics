      subroutine MatInv(A,AINV,N,NP,INDX)
C
C     takes N X N matrix A (physical dimension NP X NP) and 
C     computes A^-1 = AINV
C
C   calls subroutines
C       LUDCMP  -- LU decomposition
C          note: LUDCMP overwrites input 
C       LUBKSB  -- LU forward and back substitution
C

      implicit none
      integer n,np
      real A(NP,NP),AINV(NP,NP)
      integer INDX(Np)			! dummy
      integer i,j
      real d

      do i = 1,n
      	do j = 1,n
      		ainv(i,j) = 0.0
      	enddo
      	ainv(i,i) = 1.0
      enddo
      call ludcmp(a,n,np,indx,d)
      do j = 1,n
      	call lubksb(a,n,np,indx,ainv(1,j))
      enddo
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   subroutine LUbksb(A,n,np,indx,b)
C   applies forward- and backsubstitution to get A^-1 b
C   
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   subroutine ludcmp(A,n,np,indx,d)
C   reads in matrix A and does LU decomposition
C   from Numerical Recipes
C   
C   INPUT:
C       A(NP,NP) : real matrix; must be declared in calling routine
C      n    = dimension of A used 
C      np  = physical declared dimension of array A (np >= n)
C     indx(n) : an array declared by calling routine; used for decomposition
C     d = +/- 1; keeps track of exchanges of rows (needed to get sign of determinant) 
C       
C   OUTPUT:
C       L and U returned in matrix A; assumes L(i,i) = 1
C       note!  matrix A is destroyed in the process 
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=1000,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
