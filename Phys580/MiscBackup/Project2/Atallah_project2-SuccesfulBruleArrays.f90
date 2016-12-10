!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 2
!
!This program uses the 3pt and 5pt formulas to numerically approximate the second derivative of a function 
!of x with respect to x.  It will generate two data files; fort.8 contains the values of the difference
!between the 3pt approximation and the exact values in one column and the stepsize in the other, fort.9 
!contains the values of the difference between the 5pt approximation and the exact values in one column 
!and the stepsize in the other.  The purpose of generating these two data files will be to make a plot
!of the difference between the exact and approximate values versus the step size for both the 3pt and 5pt
!formulas.  This plot is created using a gnuplot script, 'plot.p'
!
!Functions used:
!	'f(a)'  - used to define the function of which the derivatives will be taken, default = x*sin(x)
!		- input variables: 
!			-'a' used as x in f(x)
!		- returns the value of f(a); i.e. a*sin(a) by default
!
!	'f2p(b)' - used to define exact value for 2nd derivative of f(x), 
!		 - default = (x*sin(x))'' = 2*cos(x) - x*sin(x)
!		 - input variables: 
!			-'b' used as x in f''(x)
!		 - returns the value of f''(b); i.e. 2*cos(b) - b*sin(b) by default
!
!Subroutines used:
!	'df(x,dx,fp,f3,f5)' - calculates f''(x) for user inputted x and dx values with both 3pt and 5pt
!			      formulas and also generates fort.8 and fort.9
!			    - input variables:
!				  -'x' value for which the the derivatives of f(x) will be evaluated at;
!				   inputted by user
!				  -'dx' stepsize; inputted by user
!			    -output variables:
!				  -'fp' exact 2 derivatived as calculated by function 'f2p(b)' with user
!				   inputted value of x
!				  -'f3' numerically approximated 2nd derivative using 3pt formula with
!				    user inputted values for x and dx
!				  -'f5' numerically approximated 2nd derivative using 5pt formula with
!				    user inputted values for x and dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 program project2
 implicit none	
 integer :: N, ix, iy, iz		
 double precision :: H, W, D, total_flux
	!Input variables; values inputted by user
 double precision, dimension(0:1000) :: fxyz, gxy, fx, flux !N must be less than 1000
	!arrays passed back and forth from booles rule
 double precision :: hx, hy, hz, x0, y0, x, y, z, xsum, ysum, zsum, pi
	!intermediate variables - h: stepsize with respect to each coordinate
 character :: continue = 'y'
	
pi = acos(-1.0)


 do while (continue == 'y' .or. continue == 'Y')

	print *, 'Computational Project 2:'
	print *, 'Quadrature: Neutron Flux From A Reactor'
	print *, 'Written by: Adrien Atallah'
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
	print*, 'Must be a multiple of 4'
	read*, N
	
	if ( mod(N,4) /= 0 ) then
		print *, 'N value must be a multiple of 4'
		print *, 'Please try again'
		print *, ' '
		goto 100
	endif

	hx = D/N
	hy = W/N
	hz = H/N
	
	ix = 0
	iy = 0
	!Do ix = 0, N
		x = ix*hx !vary x by stepsize i.e. x = x0 + ihx; x0 = 0
		
		!Do iy = 0, N
			y = iy*hy
			
			Do iz = 0, N			
				z = iz*hz
			
			fxyz(iz) = ( 1.0/(4.0*pi) )*( 1.0/( (x+x0)**2 + (y-y0)**2 + (z)**2 ) ) 
			print *, 'fxyz(',iz,') =  ', fxyz(iz)	
			
			EndDo
			
		Call BRule(fxyz, hz, N, zsum)  ! Call subroutine to calculate approximate integral using booles rule
		print *, 'zsum =  ', zsum
		
		!fy( jy + 1 ) = gxy( jy+1 ) !fill y integrand to send to booles rule, only problem is gxy has N/4 elements when fy should have N		
		!EndDo
		
	!Call BRule(fy,hy,N,fxf)
		
	!fx( jy + 1 ) = fxf( jy+1 ) !fill x integrand to send to booles rule, only problem is fxf has N/4 elements when fx should have N
	!EndDo	
	
	!Call BRule(fx,hx,N,flux)

	!total_flux = Sum(flux)
	
	!print *, 'total flux =  ', total_flux

	
10	print *, 'Would you like to run the program again?'
	print *, 'type "y" for yes, "n" for no'
	print *, ' '
	read *, continue	
	
	if (continue /= 'n' .and. continue /= 'N' .and. continue /= 'y' .and. continue /= 'Y') then
		print *, 'Invalid entry, please try again'
		print *, ' '
		goto 10
	endif
end do


end program project2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BRule(farray, h, N, total)
implicit none
	
integer :: N, j	
	 !Input variables
 double precision :: farray(0:1000), k, total, h, I(0:1000)
	
total = 0
k = (2.0*h)/45.0

	 Do j = 0, N-4, 4	
	
		I(j) = k*( 7.0*farray(j) + 32.0*farray(j+1) + 12.0*farray(j+2) + 32.0*farray(j+3) + 7.0*farray(j+4) )		
		
		print *, 'I(', j, ') =  ', I(j)
		total = total + I(j)
	enddo	


return
end

