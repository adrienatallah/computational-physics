      program matrixexample
! C
! C   illustration inverting matrices
! C   using LU decomposition
! C
! C  CWJ SDSU March 2005
! C
! C  uses numerical receipes program LUdcmp and LUbksb
! C 

      implicit none

! C.............first dimension arrays..................

      integer size        ! declared dimension of arrays
      parameter(size = 10)

      real a(size,size)   ! the matrix in question
      real atemp(size,size) ! needed as disposable array
      real ainv(size,size)    ! = A^-1
      
      integer matdim   ! actual dimension used

      integer indx(size)  ! needed for LU routines 

!C..................dummy indices for copying matrices........

      integer i,j
!C........... get the matrix from a file

      call getmatrix(size,matdim,a)
      
      print*,' The original matrix is: '
      call printnicematrix(size,matdim,a)
      print*,' '

!C................ copy original file 

      do i =1,matdim
          do j = 1,matdim
              atemp(i,j) = a(i,j)
          enddo
      enddo
!C.................. invert..................

      call matinv(atemp,ainv,matdim,size,indx)
      
!C.................print out inverted.........

      print*,' inverted matrix = '
      call printnicematrix(size,matdim,ainv)
    
      print*,' ' 

!C.................test................

      print*,' As a test, compute A^-1 A  = '
      call multmat(size,matdim,ainv,a,atemp)
    
      call printnicematrix(size,matdim,atemp)

      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getmatrix(size,n,array)

      implicit none

!C...................... matrix info......................

      integer size
      integer n
      real array(size,size)

!C.....................file variables .......................

      logical flag
      integer ilast
      character filename*15
!C...................dummy indices.................

      integer i,j

!C...................get name of file................

       flag = .false.

      do while(.not.flag)
          write(6,*)' Enter name of .dat file where matrix is stored'
          write(6,*)'(Do not enter the .dat extension )'
          read(5,*)filename
          ilast = index(filename,' ')-1
          open(unit=1,file=filename(1:ilast)//'.dat',status='old', err=10)
          flag=.true.
10      continue
          if(.not.flag)print*,'That file does not exist '
      enddo
      read(1,*)n
      print*,' Dimension of matrix is ',n
   
      do i = 1,n
         read(1,*,err=101)(array(i,j),j=1,n)
      enddo

      return

  101 continue
         write(6,*)' not enough matrix elements '
          stop    

      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine printnicematrix(size,n,array)

      implicit none

      integer size,n
      real array(size,size)
      integer i,j

      do i =1,n
           write(6,100)(array(i,j),j=1,n)
       enddo
  100 format(10f8.3)
       return
       end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine multmat(size,n,a,b,c)
      implicit none
      
      integer size
      real a(size,size),b(size,size),c(size,size)
      integer n
      integer i,j,k
      real temp


      do i = 1,n
        do j = 1,n
           temp = 0.0
           do k = 1,n
              temp = temp+a(i,k)*b(k,j)
           enddo
           c(i,j)=temp
        enddo
      enddo
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine MatInv(A,AINV,N,NP,INDX)
! C
! C     takes N X N matrix A (physical dimension NP X NP) and 
! C     computes A^-1 = AINV
! C
! C   calls subroutines
! C       LUDCMP  -- LU decomposition
! C          note: LUDCMP overwrites input 
! C       LUBKSB  -- LU forward and back substitution
! C

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


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! C   subroutine LUbksb(A,n,np,indx,b)
! C   applies forward- and backsubstitution to get A^-1 b
! C   
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
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! C   subroutine ludcmp(A,n,np,indx,d)
! C   reads in matrix A and does LU decomposition
! C   from Numerical Recipes
! C   
! C   INPUT:
! C       A(NP,NP) : real matrix; must be declared in calling routine
! C      n    = dimension of A used 
! C      np  = physical declared dimension of array A (np >= n)
! C     indx(n) : an array declared by calling routine; used for decomposition
! C     d = +/- 1; keeps track of exchanges of rows (needed to get sign of determinant) 
! C       
! C   OUTPUT:
! C       L and U returned in matrix A; assumes L(i,i) = 1
! C       note!  matrix A is destroyed in the process 
! C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=1000,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) WRITE (*,*) 'Singular matrix.'
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
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
