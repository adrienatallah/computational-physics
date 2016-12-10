!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 2
!
!This program calculates the flux of neutron emitting reactor measured by a detector.
! The user enters values for the spatial dimensions of the reactor, the distance from the 
! reactor to the detector as well as the number of lattice points in each direction.
! This data is sent to a subroutine,  'TripleInt(H, W, D, x0, y0, N, flux),' that sets up 
! and approximates a triple integral using Boole's Rule.  It also calls two other subroutines, 
! 'LatticeDependence(),' generates a data file, 'fort.8' with the number of lattice points in one
! column and the error in the other; 'InverseSquareDeviations(),' generates a data file, 'fort.9'
! containing the value of x0 in the first column, the flux from the reactor in the second 
! column, the flux from a point source in the third, the flux using Monte Carlo Integration in the fourth
! and the error of the Monte Carlo approximation in the fifth. these two data files will then be used to 
! create 2 plots using gnuplot.  Load 'plot-error.p' to plot 'fort.8' and 'plot-InverseSquare.p'
! to plot 'fort.9'
! 
!	Input variables: 
!		-H: Heigth of the reactor, inputted by user
!		-W: Width of the reactor, inputted by user
!		-D: Depth of the reactor, inputted by user
!		-x0: Distance from the front of the reactor to the detector, inputted by user
!		-y0: Distance from the side of the reactor to the detector, inputted by user
!		-N: Number of lattice points in each direction, inputted by user
!
!	Output variables:
!		-User_Values_Flux: the approximated value of the flux triple integral using 
!						 Boole's rule with the user inputted values
!
! 	Subroutines used: 
!		-'TripleInt(H, W, D, x0, y0, N, User_Values_Flux)'
!			-Approximates the triple integral using Boole's rule with the user
!			 inputted values
!		-'LatticeDependence()'
!			-generates fort.8 data file
!		-'InverseSquareDeviations() '
!			-generates fort.9 data file
!		-'MonteCarlo(H, W, D, x0, y0, Np, seed, flux, error) '
!			-computes flux using Monte Carlo integration method
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 program project2
 implicit none	
 integer :: N, Np, seed 	
 double precision :: H, W, D, x0, y0
	!Input variables; values inputted by user
		
 double precision :: User_Values_Flux, MonteCarlo_Flux, MonteCarlo_Error
 	!Output variable
 
 character :: continue = 'y'

 do while (continue == 'y' .or. continue == 'Y')

 	print *, ' '
	print *, 'Computational Project 2:'
	print *, 'Quadrature: Neutron Flux From A Reactor'
	print *, 'Written by: Adrien Atallah'
	print *, ' '

	print *, ' '
	print*, 'Enter the height (H) of the reactor in meters'
	print *, ' '
	read*, H
	
	print *, ' '
	print*, 'Enter the Width (W) of the reactor in meters'
	print *, ' '
	read*, W
	
	print *, ' '
	print*, 'Enter the Depth (D) of the reactor in meters'
	print *, ' '
	read*, D
	
	print *, ' '
	print*, 'Enter the distance (x0) from the reactor to detector in meters'
	print *, ' '
	read*, x0
	 
	print *, ' '
	print*, 'Enter the distance (y0) from the corner of the reactor to the detector'
	print *, ' '
	read*, y0
	
100	print *, ' '
	print*, 'Enter the number of lattice points (N) in each direction (Nx=Ny=Nz=N)'
	print*, 'Must be less than 1000 and a multiple of 4'
	read*, N
	
	if ( mod(N,4) /= 0 .or. N >= 1000) then
		print *, ' '
		print *, 'N value must be less than 1000 and a multiple of 4'
		print *, 'Please try again'
		print *, ' '
		goto 100
	endif
	
	print *, ' '
	print*, 'Enter the number of test points to use for Monte Carlo Integration:'
	print *, ' '
	read*, Np
	
	print *, ' '
	print*, 'Choose the Seed:'
	print *, ' '
	read*, seed

	call TripleInt(H, W, D, x0, y0, N, User_Values_Flux)
	
	call MonteCarlo(H, W, D, x0, y0, Np, seed, MonteCarlo_Flux, MonteCarlo_Error)
	
	print *, ' '
	print *, 'with the values you entered, the total flux (using Booles approximation) =  ', User_Values_Flux	
	print *, ' '
	print *, 'with the values you entered, the total flux (using Monte Carlo Integration) =  ', MonteCarlo_Flux	
	print *, 'The Error on the Monte Carlo Integration  =  ', MonteCarlo_Error	
	print *, ' '
	
	call LatticeDependence() 
		! generate data file with number of lattice points vs error
			
	call InverseSquareDeviations() 
		! generate data file with point source flux and detector flux as x0 varies
	
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

end program project2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TripleInt(H, W, D, x0, y0, N, flux) takes in the values for the spatial dimensions of 
! the reactor: H, W, D; the distance of the detector from the reactor: x0, y0;
! as well as the number of lattice points: N and approximates the value of the flux 
! triple integral using Boole's rule
!
!	Input variables: 
!		-H: Heigth of the reactor
!		-W: Width of the reactor
!		-D: Depth of the reactor
!		-x0: Distance from the front of the reactor to the detector
!		-y0: Distance from the side of the reactor to the detector
!		-N: Number of lattice points in each direction
!
!	Intermediate variables:
!		-ix, iy, iz: counters for the x, y, z loops respectively
!		-fxyz: array of integrand with varied z to be sent to Boole's rule subroutine
!		-gxy: array of integrand with varied y to be sent to Boole's rule subroutine
!		-fx: array of integrand with varied x to be sent to Boole's rule subroutine
!		-hx, hy, hz: step size in the x, y, z directions respectively
!		-x, y, z: cartesian spatial coordinates 
!		-sum: value returned by Boole's rule subroutine, represets the approximate
!			   integration of whatever integrand array was sent
!
!	Output variables:
!		-flux: the approximated value of the flux triple integral using Boole's rule
!
! 	Subroutines used: 
!		-'BRule(farray, h, N, total)'
!			-Calculates Boole's Rule approximation of farray integrand
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TripleInt(H, W, D, x0, y0, N, flux)
implicit none

 integer :: N		
 double precision :: H, W, D, x0, y0
	!Input variables
		
 double precision, dimension(0:1000) :: fxyz, gxy, fx !N must be less than 1000
 double precision :: hx, hy, hz, x, y, z, sum, pi
 integer :: ix, iy, iz
	!Intermediate variables
		
 double precision :: flux
	!Output variable
		
	pi = acos(-1.0)
	
	hx = D/N
	hy = W/N
	hz = H/N
	
	Do ix = 0, N
		x = ix*hx !vary x by stepsize i.e. x = x0 + ihx; x0 = 0
		
		Do iy = 0, N
			y = iy*hy
			
			Do iz = 0, N			
				z = iz*hz
			
				fxyz(iz) = ( 1.0/(4.0*pi) )*( 1.0/( (x+x0)**2 + (y-y0)**2 + (z)**2 ) ) 				
			
			EndDo
			
			Call BRule(fxyz, hz, N, sum)  ! Calculate approximate integral of z integrand using Boole's rule			
		
			gxy(iy) = sum !fill y integrand to send to Boole's rule
		
		EndDo
		
		Call BRule(gxy, hy, N, sum)  ! Calculate approximate integral of y integrand using Boole's rule
		
		fx(ix) = sum !fill y integrand to send to Boole's rule
			
	EndDo	
	
	Call BRule(fx, hx, N, flux) ! finally call subroutine to calculate approximate integral of x integrand and get flux value

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BRule(farray, h, N, total) takes in array of integrand to be integrated, stepsize, and 
! number lattice points, then applies Boole's Rule to approximate the integral
!
!	Input variables: 
!		-farray: array of integrand
!		-h: stepsize
!		-N: number of lattice points
!		
!	Intermediate variables:
!		-j: counter for loop
!		-I: array with an element for each iteration of Boole's rule sum
!		-k: constant used in Boole's rule
!
!	Output variables:
!		-total: result of Boole's rule applied to integrand array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BRule(farray, h, N, total)
implicit none
	
 integer :: N
 double precision :: farray(0:1000), h
 	!Input variables
 		
 double precision :: I(0:1000), k
 integer :: j
	!Intermediate variables
		
 double precision :: total
 	!Output variable
 	
total = 0
k = (2.0*h)/45.0

	 Do j = 0, N-4, 4	!Generate Boole's rule sum
	
		I(j) = k*( 7.0*farray(j) + 32.0*farray(j+1) + 12.0*farray(j+2) + 32.0*farray(j+3) + 7.0*farray(j+4) )		
		
		total = total + I(j)
	enddo	


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!LatticeDependence() uses fixed values for H, W, D, x0, y0 and varies the number
! of lattice points, N, to generate a data file, 'fort.8' with the difference between
! the exact value of the flux and the value with each number of lattice points in
! one column versus the number of lattice points in the other column.  The exact
! flux value is calculated using a large number of lattice points, N = 400.
!
!	Intermediate variables:
!		Nv: varied number of lattice points
!		exact: flux calculated using fixed values and N = 400
!		error: difference between exact and approximate at each N value
!		approx: approximated flux at each N value
!		x0, y0, H, W, D: fixed spatial dimensions and distance between
!					   reactor and detector (x0, y0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
 subroutine LatticeDependence()
 implicit none
 
 integer :: N, Nv !Nv will be varying number of lattice points
 double precision :: x0, y0, H, W, D, exact, error, Exact_Flux, approx
	! all intermediate variables 

 H = 10.0 ! fix the spatial dimensions of the reactor with these values
 W = 20.0
 D = 30.0
 
 x0 = 100.0 ! fix the the location of the detector with these x0 and y0 values
 y0 = 100.0
 
 ! to get an accurate approximation of the flux with these fixed values, I will use 400
 ! lattice points (N=400) and call it the exact value of the flux 
 N=400
 
 ! calculate the exact flux: 	 
 call TripleInt(H, W, D, x0, y0, N, Exact_Flux)
 exact = Exact_Flux
 
 ! now I will vary the number of lattice points (N) from 4 to 100 to generate a 
 ! data file, 'fort.8,' that will have the value of N in the first column and the difference
 ! between the exact flux value and the approximate in the second column.
 ! This file will be used to illustrate the error dependence on the lattice spacing
 ! with a plot of N vs. the difference in exact and approximate flux values
 
   Do Nv = 4, 100, 4 ! varying number of lattice points to illustrate the dependence 
  
  	 call TripleInt(H, W, D, x0, y0, Nv, approx)
  	 
  	 error = abs(approx - exact)
  	 
  	 write(8,*)Nv, error !write N vs error to fort.8
 
  EndDo
 
 return
 end
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! InverseSquareDeviations() uses fixed values for H, W, D, y0 and varies the distance
! from the front of the reactor to the detector, x0, to generate a data file, 'fort.9' 
! containing the value of x0 in the first column, the flux from the reactor in the second 
! column, the flux from a point source in the third, the flux using Monte Carlo Integration in the fourth
! and the error of the Monte Carlo approximation in the fifth.  This data will be used to 
! illustrate the devations from the inverse square law as the distance between
! the reactor and detector changes.
!
!	Intermediate variables:
!		N: number of lattice points
!		x0: varying distance between reactor and detector
!		reactor: flux from reactor calculated using fixed values 
!		ptSourceFlux: flux from a point source using same fixed values as reactor
!		H, W, D, y0: fixed spatial dimensions and distance between
!				     the side of the reactor and detector (y0)
! 
!	Subroutines used: 
!		-'MonteCarlo(H, W, D, x0, y0, Np, seed, flux, error) '
!			-computes flux using Monte Carlo integration method
!
!	Functions used:
!		'ptSource(H, W, D, x0, y0)'
!			-takes in spatial dimensions and distance between point source
!			 and detector; returns flux of point source
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine InverseSquareDeviations()	
 integer :: N, Np, seed
 double precision :: x0, y0, H, W, D, reactor, ptSourceFlux, Mflux, Merror
 double precision :: ptSource
 
 H = 10.0 ! fix the spatial dimensions of the reactor with these values
 W = 20.0
 D = 30.0
 
 N=40 ! fix the number of lattice points
 
 y0 = 10.0 ! fix the the location of the detector with respect to y0 
 
 Np = 1000 ! fix the number of test points for monte carlo integration
 seed = -1 ! choose seed value for random number generator 
 
 ! now I will vary x0, the distance from the reactor to the detector in order
 ! to illustrate the relationship between the flux from the reactor, the flux 
 ! from a point source and the flux using monte carlo integration as x0 gets larger
 
   Do j = 1, 100, 2 ! varying distance from reactor to detector, x0  	 
   	 x0 = j  	
  	 
   	 call TripleInt(H, W, D, x0, y0, N, reactor)  	 
  	 ptSourceFlux = ptSource(H, W, D, x0, y0)
  	 
  	 call MonteCarlo(H, W, D, x0, y0, Np, seed, Mflux, Merror)
  	 
  	 write(9,*)x0, reactor, ptSourceFlux, Mflux, Merror
  	 	!write x0 vs. flux of reactor, flux of pt source and monte carlo flux with error to 'fort.9'
 
  EndDo
 
 return
 end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function ptSource(H, W, D, x0, y0) returns the value of the flux from a point source 
! located at the center of the reactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 double precision function ptSource(H, W, D, x0, y0)
 implicit none
 double precision, intent(in) :: H, W, D, x0, y0
 double precision :: pi
 
 pi = acos(-1.0)
 
 ptSource = (W*H*D)/(4*pi) * (1/( (H/2.0)**2 + (x0 + (D/2.0))**2 + (y0 - (W/2.0))**2) )

 return
 end 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MonteCarlo(H, W, D, x0, y0, Np, seed, flux, error) 
!
!	Input variables: 
!		-H: Heigth of the reactor
!		-W: Width of the reactor
!		-D: Depth of the reactor
!		-x0: Distance from the front of the reactor to the detector
!		-y0: Distance from the side of the reactor to the detector
!		-Np: Number of test points
!		-seed: Integer for use in random number generation
!
!	Intermediate variables:
!		-i: counter for loop
!		-x, y, z: cartesian coordinates in flux integrand
!		-fsum: sum of evaluated integrand function in Monte Carlo sum
! 		-fsquared: store the square of integrand function
!		-fsquared_sum: sum of fsquared for Monte Carlo error
!		-var: variance of Monte Carlo integration
!		-rand: store random number generated by RAND3
!
!	Output variables:
!		-flux: flux using Monte Carlo integration
!		-error: standard deviation of Monte Carlo approximation
!
!
!	Functions used:
!		 'f(x0, y0, x, y, z)'
!			-evaluates integrand in flux integral at for x0, y0, x, y, z
!		'RAN3(IDUM)'
!			-takes in IDUM (seed) and returns a random number 
!			 of uniform distribution between 0.0 and 1.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine MonteCarlo(H, W, D, x0, y0, Np, seed, flux, error)	
 integer, intent(in) :: Np, seed
 double precision, intent(in) ::  x0, y0, H, W, D
 	! Input variables
 double precision, intent(out) :: flux, error
 	! Output variables	
 integer :: i
 double precision :: x, y, z, V, f, fsum, fsquared, fsquared_sum, var
 real :: RAN3, factor, rand1, rand2, rand3
 	! Intermediate variables 	
 
 V = H*W*D
 fsum = 0
 fsquared_sum = 0
 
 Do i = 1, Np
 	rand1 = RAN3(seed) !store random number generated by RAN3
 	rand2 = RAN3(seed)
 	rand3 = RAN3(seed)
 	
 	!if (rand > 1.0) then
 		!rand = rand - 0.1
 	!end if 
 	
 	x = rand1 * D ! Set x to random number between 0 and D
 	y = rand2 * W ! Set y to random number between 0 and W
 	z = rand3 * H ! Set z to random number between 0 and H
 
	fsum = fsum + f(x0, y0, x, y, z)
	
	fsquared = f(x0, y0, x, y, z)**2
	
	fsquared_sum = fsquared_sum + fsquared
 
 EndDo
 
 flux = (V/Np)*fsum
 var = (1.0/Np)*fsquared_sum - ( (1.0/Np) * fsum )**2
 
 error = V*((var/Np)**(1.0/2.0))
   	
 return
 end
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function f(x0, y0, x, y, z) returns the integrand in the flux integral evaluated at
! x0, y0, x, y, z for use in the Monte Carlo Integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 double precision function f(x0, y0, x, y, z)
 implicit none
 double precision, intent(in) :: x0, y0, x, y, z
 double precision :: pi
 
 pi = acos(-1.0)
 
 f = ( 1.0/(4.0*pi) )*( 1.0/( (x+x0)**2 + (y-y0)**2 + (z)**2 ) ) 

 return
 end 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  function ran3
!
!  generates uniformly distributed random numbers 
!  between 0.0 and 1.1
!  taken from Numerical Recipes (fortran version) 
!
!  common block need to store inext, etc 
!
! INPUT: IDUM  integer "seed" 
!   note: best if IDUM is negative first call, in order 
!         to initialize
! OUTPUT: ran3 
!
      real FUNCTION RAN3(IDUM)
!         IMPLICIT REAL*4(M)
!         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      common/randumb/inext,inextp,ma     ! needed on HP
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END
