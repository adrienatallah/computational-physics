Program project3
 implicit none
 integer::i, j, index(3) 
 double precision::data(15, 3), G(3, 3), Ginv(3,3), Gtemp(3, 3), Psi(3), R(3, 3)
    
 character :: continue = 'y' , view = 'n'

 do while (continue == 'y' .or. continue == 'Y')    

	call ReadData(data) !call subroutine to read in data from file
	
	call CalculateGmatrix(data, G) !call subroutine to create matrix G
		
		Gtemp = G

	call CalculatePsiMatrix(data, Psi) !call subroutine to create Matrix Psi
	
	call MatInv(Gtemp, Ginv, 3, 3, index)
	
	call multmat(3, 3, Ginv, G, R)
	
	 print *, ' '
	 print *,'G^-1 matrix:'
	 print *, ' '
 
	 do i = 1, 3     
	 	write(6, '(10f20.8)') (Ginv(i, j), j = 1, 3)
			!print G^-1 matrix to screen in 3x3 form with 2decimal places of accuracy 			
		enddo 
		
	print*, ' '
	
	print *, ' '
	print *,'R matrix:'
	print *, ' '
 
	 do i = 1, 3     
	 	write(6, '(10f20.2)') (R(i, j), j = 1, 3)
			!test to see if G^-1*G = Identity matrix		
		enddo 
		
	print*, ' '
	
!!!!!!!!!!!!!
!Ask user if they want to view data
!!!!!!!!!!!!!
	print *, ' '
100	print *, 'Would you like to view the data?'
	print *, 'type "y" for yes, "n" for no'
	print *, ' '
	read *, view	
	
	
	if (view /= 'n' .and. view /= 'N' .and. view /= 'y' .and. view /= 'Y') then
		print *, ' '
		print *, 'Invalid entry, please try again'
		print *, ' '
		goto 100
	endif
	
	if (view == 'y' .or. view == 'Y') then
		call PrintData(data) !call subroutine that prints data
	endif

!!!!!!!!!!!!!
!Ask user if they want to run the program again
!!!!!!!!!!!!!!	
	print *, ' '
10	print *, 'Would you like to run the program again?'
	print *, 'type "y" for yes, "n" for no'
	print *, ' '
	read *, continue	
	
	
	if (continue /= 'n' .and. continue /= 'N' .and. continue /= 'y' .and. continue /= 'Y') then
		print *, ' '
		print *, 'Invalid entry, please try again'
		print *, ' '
		goto 10
	endif
end do
 
END  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ReadData(data) Reads in data from file, 'dataset4.dat' wich contains climate 
! data for 15 years.  The first column contains the year, the second column
! is the temperature for that year and the 3rd column is the uncertainty for 
! the temperature measured that year.  This data is read into a 15x3 array, 
! 'data,' which the subroutine returns. 
!
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadData(data)
implicit none

 integer :: i, j 
 double precision, intent(out) :: data(15, 3)

 do i = 1, 15 !initialize 15x3 data matrix to 0
	do j = 1, 3
		data(i, j) = 0.0  		
	enddo
 enddo    

 OPEN(unit=1, FILE='dataset4.dat', status='old', form='formatted')  
	!open 'dataset4.dat,' set to unit = 1

 do i = 1, 15      
	read(1,*) (data(i ,j), j = 1, 3)		
		!read in and store the data on 'dataset4.dat' to 'data' array
 enddo 	 

close(1) !close data file	
		
return
end	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PrintData(data)
!
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintData(data)
implicit none

 double precision, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j 
 double precision :: year(15), temp(15), var(15) 

 do i = 1, 15
	 temp(i) = 0.0 !Initialize arrays to 0	
	 year(i) = 0.0
	 var(i) = 0.0		
 enddo    

!print data from file
print*, ' '	
print*, 'the data from the file is:'
print*, '1st column: year'
print*, '2nd column: temperature'
print*, '3rd column: Uncertainty'
print*, ' '

 do i = 1, 15      
	write(6, '(15f8.2)') (data(i, j), j = 1, 3)
		!print data to screen in 15 rows by year with 2decimal places of accuracy 			
 enddo 
 print*, ' '
 
 do i = 1, 15
	year(i) = data(i, 1)
	temp(i) = data(i, 2)
	var(i) = data(i, 3)
	print *, i, 'Year=', year(i)
	print *, i, 'Temperature=',temp(i)
	print *, i, 'Uncertainty=',var(i)
	print*, ' '
 enddo
	 
return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalculateGmatrix(data)
!
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculateGmatrix(data, Gmat)
implicit none

 double precision, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j, k 
 double precision :: gy, sum
 	!intermediate variables
 double precision, intent(out) :: Gmat(3, 3)
 	!output variable
 		
 sum = 0
 
 do i = 1, 3
 	do j = 1,3
	 	Gmat(i, j) = 0.0 !Initialize Gmat array to 0	 	
	 enddo
 enddo    

 do i = 1, 3
 	do j = 1,3
 		do k = 1, 15
 			sum = sum + gy(i, k, data)*gy(j, k, data) 			
 		enddo
 		
 		Gmat(i, j) = sum
 		sum = 0
	 enddo
 enddo  

 print *, ' '
 print *,'G matrix:'
 print *, ' '
 
 do i = 1, 3     
	write(6, '(10f20.2)') (Gmat(i, j), j = 1, 3)
		!print G matrix to screen in 3x3 form with 2decimal places of accuracy 			
 enddo 
 print*, ' '

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalculatePsiMatrix(data)
!
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculatePsiMatrix(data, PsiMat)
implicit none

 double precision, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j, k 
 double precision :: gy, sum
 	!Intermediate variables
 double precision, intent(out) :: PsiMat(3)
 	!output variable
 
 
 sum = 0
 
 do i = 1, 3 
	 PsiMat(i) = 0.0 !Initialize PsiMat array to 0	 	
 enddo    

 do i = 1, 3
 	do j = 1,15
 		sum = sum + data(j, 2)*gy(i, j, data) 			
 	enddo
 		
 	PsiMat(i) = sum
 	sum = 0
	
 enddo  

 print *, ' '
 print *,'Psi matrix:'
 print *, ' '
 
 do i = 1, 3
 	print*, PsiMat(i)
 enddo

 print*, ' '

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function gy(i, yj, data) returns g for a certain y value, using an index number
! called 'yj'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 double precision function gy(i, yj, data)
 implicit none
 
 double precision, intent(in) :: data(15, 3)
 integer, intent(in) :: i, yj
 
if (i == 1) then	
	gy = 1
	
else if ( i == 2) then
	gy = data(yj, 1) - 1900
	
else if ( i == 3) then
	gy = (data(yj, 1) - 1900)*(data(yj, 1)-1900)
	
end if

 return
 end 
 
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
      double precision A(NP,NP),AINV(NP,NP)
      integer INDX(Np)			! dummy
      integer i,j
      double precision d

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine 'multmat(size, n, a, b, c) multiplies two matrices
! a, b; returns the result.  From 'matrixexample.f' off website
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
      subroutine multmat(size,n,a,b,c)
      implicit none
      
      integer size
      double precision a(size,size),b(size,size),c(size,size)
      integer n
      integer i,j,k
      double precision temp


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
