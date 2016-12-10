!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 3
!
!This program takes in climate data from the file, 'dataset4.dat' and does a
! linear least squares fit using polynomials.  'dataset4.dat' contains 15 temperature 
! measurements, each corresponding to a different year.  The file is organized
! in 3 columns; the first being the year, the second being the temperature measured
! that year and the third column contains the uncertainty on that measurement.
! This program creates 3 main data files; the first, 'fort.7' is the temperature vs. the year
! calculated by a 2nd degree polynomial; the second, 'fort.8' is the same data but
! calculated by a 3rd degree polynomial and the third, 'fort.9' is calculated from
! a fourth degree polynomial.  Two other data files are also created; 'fort.3' is used
! to print the solutions for the coefficeints for each polynomial and 'fort.10' is data
! for an extrapolation to the year 2050.  Two plots are made from these data files; the
! first uses the gnuplot script, 'plot.p' and plots the data on 'dataset4.dat' with error
! bars along with the 2nd and 3rd degree polynomials using the data on 'fort.7' and
! 'fort.8.'  The second is the temperature extrapolated to the year 2050, plotted by
! gnuplot script, 'plot_extra.p' which uses the data on 'fort.10.'  
! 
! 
!	Intermediate variables: 
!		-'data': array that reads in data off 'dataset4.dat'
!		-'G2', 'G3', 'G4': arrays containg G matrix values for 2nd 3rd and 4th 
!					   order polynomials respectively
!		-'Gtemp2', 'Gtemp3', 'Gtemp4': used for temporary storage of G's
!		-'Ginv2', 'Ginv3', 'Ginv4': inverse of G matrix for each order polynomial
!		-'Psi2', 'Psi3', 'Psi4': arrays to store Psi matrix for each order polynomial
!		-'C2', 'C3', 'C4': arrays to store coefficient matrix for each order polynomail
!
! 	Subroutines used: 
!		-'ReadData(data)'
!			-Reads in data from 'dataset4.dat' and stores in 'data' array
!		-'CalculateGmatrix(data, G2, G3, G4)'
!			-takes in array 'data' and calculates G matricies
!		-'CalculatePsiMatrix(data, Psi2, Psi3, Psi4) '
!			-takes in array 'data' and calculates Psi matricies
!		-'MatInv(Gtemp2, Ginv2, 3, 3, index)'
!			-computes inverse of a square matrix, in this case 'Gtemp2'
!		-'CalculateCoeff(Ginv2, Ginv3, Ginv4, Psi2, Psi3, Psi4, C2, C3, C4)'
!			-takes in G inverse and Psi matrices for each order polynomial
!			 and calculates coefficient matrices for each order polynomial and
!			 also creates 'fort.3'
!		-'CalculateTemp(data, C2, C3, C4)'
!			-calcuates temperatures for each order polynomials and creates
!			 data files 'fort.7', 'fort.8', 'fort.9'
!		-'Extrapolate(C2)'
!			-uses 2nd order polynomial to extrapolate the temperature to 
!			 the year 2050 and creates 'fort.10'
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program project3
 implicit none
 integer :: i, j, index(3) 
 real :: data(15, 3), G2(3, 3), Gtemp2(3, 3), G3(4, 4), Gtemp3(4, 4), G4(5, 5), Gtemp4(5, 5)
 real :: Ginv2(3, 3), Ginv3(4, 4), Ginv4(5, 5), Psi2(3), Psi3(4), Psi4(5), C2(3), C3(4), C4(5)
    
 character :: continue = 'y' , view = 'n'

  do i = 1, 3
 	do j = 1,3
	 	Gtemp2(i, j) = 0.0 !Initialize Gtemp2 array to 0	 	
	 enddo
 enddo    
 
 do i = 1, 4
 	do j = 1,4
	 	Gtemp3(i, j) = 0.0 !Initialize Gtemp3 array to 0	 	
	 enddo
 enddo    
 
 do i = 1, 5
 	do j = 1,5
	 	Gtemp4(i, j) = 0.0 !Initialize Gtemp4 array to 0	 	
	 enddo
 enddo    
 
 do while (continue == 'y' .or. continue == 'Y')    

	call ReadData(data) !call subroutine to read in data from file
	
	call CalculateGmatrix(data, G2, G3, G4) !call subroutine to create G matricies
	
	Gtemp2 = G2
	Gtemp3 = G3
	Gtemp4 = G4
		!use temporary array Gtemp to send to Matinv since it will be destroyed

	call CalculatePsiMatrix(data, Psi2, Psi3, Psi4) !call subroutine to create matrix Psi	
	
	call MatInv(Gtemp2, Ginv2, 3, 3, index) !call subroutine to calculate G inverses
	call MatInv(Gtemp3, Ginv3, 4, 4, index)
	call MatInv(Gtemp4, Ginv4, 5, 5, index)
	
	call CalculateCoeff(Ginv2, Ginv3, Ginv4, Psi2, Psi3, Psi4, C2, C3, C4)
		!call subroutine to calculate coefficients for each order polynomial
			
	call CalculateTemp(data, C2, C3, C4)
		!call subroutine to create data files of year vs. temperature for each order polynomial
			
	call Extrapolate(C2)
		!call subroutine to extrapolate temperature to the year 2050
			
!!!!!!!!!!!!!
!Ask user if they want to view data
!!!!!!!!!!!!!
	print *, ' '
100	print *, 'Would you like to view the data on ''dataset4.dat''?'
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
!ReadData(data) reads in data from file, 'dataset4.dat' wich contains climate 
! data for 15 years.  The first column contains the year, the second column
! is the temperature for that year and the 3rd column is the uncertainty for 
! the temperature measured that year.  This data is read into a 15x3 array, 
! 'data,' which the subroutine returns. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadData(data)
implicit none

 integer :: i, j 
 real, intent(out) :: data(15, 3)

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
!PrintData(data) takes in data as read in by subroutine, 'ReadData(data)'
! and prints the data to screen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintData(data)
implicit none

 real, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j 
 real :: year(15), temp(15), var(15) 

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
!CalculateGmatrix(data, Gmat2, Gmat3, Gmat4) takes in array 'data' and returns 
! G matricies as arrays for each order polynomial.
!
! 	Functions used: 
!		-'gy(i, yj, data)': returns g for a certain y value, using an index number
! 					  called 'yj'
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculateGmatrix(data, Gmat2, Gmat3, Gmat4)
implicit none

 real, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j, k 
 real :: gy, sum
 	!intermediate variables
 real, intent(out) :: Gmat2(3, 3), Gmat3(4, 4), Gmat4(5, 5)
 	!output variable
 		
 sum = 0
 
 do i = 1, 3
 	do j = 1,3
	 	Gmat2(i, j) = 0.0 !Initialize Gmat2 array to 0	 	
	 enddo
 enddo    
 
  do i = 1, 4
 	do j = 1,4
	 	Gmat3(i, j) = 0.0 !Initialize Gmat3 array to 0	 	
	 enddo
 enddo    
 
  do i = 1, 5
 	do j = 1,5
	 	Gmat4(i, j) = 0.0 !Initialize Gmat4 array to 0	 	
	 enddo
 enddo    
 

 ! create array of G matrix for 2nd order ploynomial
 do i = 1, 3
 	do j = 1,3
 		do k = 1, 15
 			sum = sum + gy(i, k, data)*gy(j, k, data)*(1/(data(k,3)**2.0)) 			
 		enddo
 		
 		Gmat2(i, j) = sum 
 		sum = 0
	 enddo
 enddo  
 
 ! create array of G matrix for 3rd order ploynomial
 do i = 1, 4
 	do j = 1,4
 		do k = 1, 15
 			sum = sum + gy(i, k, data)*gy(j, k, data)*(1/(data(k,3)**2.0)) 			
 		enddo
 		
 		Gmat3(i, j) = sum
 		sum = 0
	 enddo
 enddo  
 
 ! create array of G matrix for 4th order ploynomial
 do i = 1, 5
 	do j = 1,5
 		do k = 1, 15
 			sum = sum + gy(i, k, data)*gy(j, k, data)*(1/(data(k,3)**2.0)) 			
 		enddo
 		
 		Gmat4(i, j) = sum
 		sum = 0
	 enddo
 enddo  

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalculatePsiMatrix(data, PsiMat2, PsiMat3, PsiMat4) takes in array 'data' and 
! returns Psi matricies as arrays for each order polynomial.
!
! 	Functions used: 
!		-'gy(i, yj, data)': returns g for a certain y value, using an index number
! 					  called 'yj'
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculatePsiMatrix(data, PsiMat2, PsiMat3, PsiMat4)
implicit none

 real, intent(in) :: data(15, 3)
 	!input variables
 integer :: i, j, k 
 real :: gy, sum
 	!Intermediate variables
 real, intent(out) :: PsiMat2(3), PsiMat3(4), PsiMat4(5)
 	!output variable 
 
 sum = 0
 
 do i = 1, 3 
	 PsiMat2(i) = 0.0 !Initialize PsiMat2 array to 0	 	
 enddo    
 
 do i = 1, 4 
	 PsiMat3(i) = 0.0 !Initialize PsiMat3 array to 0	 	
 enddo 
 
  do i = 1, 5 
	 PsiMat4(i) = 0.0 !Initialize PsiMat4 array to 0	 	
 enddo 

 do i = 1, 3
 	do j = 1,15
 		sum = sum + data(j, 2)*gy(i, j, data)*(1/(data(j,3)**2.0)) 	 
 			!data(j, 2) => Temperature; data(j, 3) => Variance
 	enddo
 		
 	PsiMat2(i) = sum
 	sum = 0	
 enddo  
 
  do i = 1, 4
 	do j = 1,15
 		sum = sum + data(j, 2)*gy(i, j, data)*(1/(data(j,3)**2.0)) 	 
 			!data(j, 2) => Temperature; data(j, 3) => Variance
 	enddo
 		
 	PsiMat3(i) = sum
 	sum = 0	
 enddo  
 
  do i = 1, 5
 	do j = 1,15
 		sum = sum + data(j, 2)*gy(i, j, data)*(1/(data(j,3)**2.0)) 	 
 			!data(j, 2) => Temperature; data(j, 3) => Variance
 	enddo
 		
 	PsiMat4(i) = sum
 	sum = 0	
 enddo  

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalculateCoeff(Ginv2, Ginv3, Ginv4, Psi2, Psi3, Psi4, C2, C3, C4) takes in G 
! inverse and Psi matrices for each order polynomial and returns coefficient
! matrices as arrays: C2, C3 and C4 for each order polynomial and creates 'fort.3'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculateCoeff(Ginv2, Ginv3, Ginv4, PsiMat2, PsiMat3, PsiMat4, C2, C3, C4)
implicit none

 real, intent(in) :: Ginv2(3, 3), Ginv3(4, 4), Ginv4(5, 5), PsiMat2(3), PsiMat3(4), PsiMat4(5)
 	!input variables
 integer :: i, j 
 real :: sum
 	!Intermediate variables
 real, intent(out) :: C2(3), C3(4), C4(5)
 	!output variable 
 
 sum = 0
 
 do i = 1, 3 
	 C2(i) = 0.0 !Initialize Coeffecient array for 2nd order to 0	 	
 enddo
 
 do i = 1, 4 
	 C3(i) = 0.0 !Initialize Coeffecient array for 3rd order to 0	 	
 enddo    
 
 do i = 1, 5 
	 C4(i) = 0.0 !Initialize Coeffecient array for 4th order to 0	 	
 enddo    

 do i = 1, 3
 	do j = 1, 3
 		sum = sum + Ginv2(i, j)*PsiMat2(j) 	 
 	enddo
 		
 	C2(i) = sum
 	sum = 0	
 enddo  
 sum = 0
 
 do i = 1, 4
 	do j = 1, 4
 		sum = sum + Ginv3(i, j)*PsiMat3(j) 	 
 	enddo
 		
 	C3(i) = sum
 	sum = 0	
 enddo  
 sum = 0
 
 do i = 1, 5
 	do j = 1, 5
 		sum = sum + Ginv4(i, j)*PsiMat4(j) 	 
 	enddo
 		
 	C4(i) = sum
 	sum = 0	
 enddo  
 
  do i = 1, 5 	 
 	 write(3,*)C2(i), C3(i), C4(i)
 enddo
 
  print*, ' '

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalculateTemp(data, C2, C3, C4) takes in 'data' and each C array to calcuate 
! temperatures for each order polynomials and create data files 'fort.7', 'fort.8', 'fort.9'
!
! 	Functions used: 
!		-'T(i, C2, C3, C4, y)': returns temperature for any order polynomial (order indicated
! 						  by integer i) as a function of the coefficient array and the year.	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalculateTemp(data, C2, C3, C4)
implicit none

 real, intent(in) :: data(15, 3), C2(3), C3(4), C4(5)
 	!input variables
 integer :: i, j 
 real :: sum, T, Temp2, Temp3, Temp4, year
 	!Intermediate variables
 
 sum = 0

 do i = 1, 15
 	year = data(i, 1) - 1900
	Temp2 = T(2, C2, C3, C4, year) 
	write(7,*)data(i, 1), Temp2 
		!write year versus best fit calcuated temperature for 2nd order to fort.7
			
	Temp3 = T(3, C2, C3, C4, year) 
	write(8,*)data(i, 1), Temp3
		!write year versus best fit calcuated temperature for 3rd order to fort.8
			
	Temp4 = T(4, C2, C3, C4, year) 	
	write(9,*)data(i, 1), Temp4
		!write year versus best fit calcuated temperature for 4th order to fort.9	
	Temp2 = 0
	Temp3 = 0
	Temp4 = 0
 enddo  

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'Extrapolate(C2)' uses 2nd order polynomial to extrapolate the temperature to 
! the year 2050 and creates 'fort.10'
!
! 	Functions used: 
!		-'T(i, C2, C3, C4, y)': returns temperature for any order polynomial (order indicated
! 						  by integer i) as a function of the coefficient array and the year.	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Extrapolate(C2)
implicit none

 real, intent(in) :: C2(3)
 	!input variables
 integer :: i, j 
 real :: sum, T, Temp, year, y, C3(4), C4(5)
 	!Intermediate variables
 
 sum = 0

 do i = 1, 75
 	year = 2*i ! y goes from 2 to 150 => years 1902 to 2050
 	
	Temp = T(2, C2, C3, C4, year) 
	y = year + 1900
	write(10,*)y, Temp 
		!write year versus best fit calcuated temperature for 2nd order to fort.7
	
 enddo  

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function gy(i, yj, data) returns g for a certain y value, using an index number
! called 'yj'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function gy(i, yj, data)
 implicit none
 
 real, intent(in) :: data(15, 3)
 integer, intent(in) :: i, yj
 
if (i == 1) then	
	gy = 1
	
else if ( i == 2) then
	gy = data(yj, 1) - 1900
	
else if ( i == 3) then
	gy = (data(yj, 1) - 1900)*(data(yj, 1)-1900)

else if ( i == 4) then
	gy = (data(yj, 1) - 1900)*(data(yj, 1)-1900)*(data(yj, 1)-1900)
	
else if ( i == 5) then
	gy = (data(yj, 1) - 1900)*(data(yj, 1)-1900)*(data(yj, 1)-1900)*(data(yj, 1)-1900)
	
end if

 return
 end 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function T(i, C2, C3, C4, y) returns temperature for any order polynomial (order indicated
! by integer i) as a function of the coefficient array and the year.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function T(i, C2, C3, C4, y) 
 implicit none
 
 integer, intent(in) :: i
 real, intent(in) :: C2(3), C3(4), C4(5), y
 
 if (i == 2) then	
	T = C2(1) + C2(2)*y + C2(3)*y**2
	
else if ( i == 3) then
	T = C3(1) + C3(2)*y + C3(3)*y**2 + C3(4)*y**3
	
else if ( i == 4) then
	T = C4(1) + C4(2)*y + C4(3)*y**2 + C4(4)*y**3 + C4(5)*y**4
	
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine 'multmat(size, n, a, b, c) multiplies two matrices
! a, b; returns the result.  From 'matrixexample.f' off website
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
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
