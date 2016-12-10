!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 1
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
double precision :: x0, xf 	
integer :: N		
	!Input variables; values inputted by user
double precision :: Integral
	!Output Variable

!double precision :: H, W, D, x0, y0, n 			
	!input variables


character :: continue = 'y'

do while (continue == 'y' .or. continue == 'Y')

	print *, 'Computational Project 2:'
	print *, 'Quadrature: Neutron Flux From A Reactor'
	print *, 'Written by: Adrien Atallah'
	print *, ' '


	print *, 'Test 1D Booles Rule Subroutine, take integral of f(x)=sinx'
	print *, 'Enter values for x0, xf, N'
	print *, '(N must be an integer multiple of 4)'
	print *, ' '
	read *, x0, xf, N

	
	
100	if ( mod(N,4) /= 0 ) then
		print *, 'N must be a multiple of 4!'
		print *, 'Please re-enter N'
		print *, ' '
		read *, N
		goto 100
	endif

	
	Call BRule1D(N,x0,xf,Integral)  ! Call subroutine to calculate and approximate 2nd derivatives of a function of x


	print *, 'Integral of f(x)=sinx : ', Integral

	
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

subroutine BRule1D(N,x0,xf,Integral)
implicit none

double precision :: x0, xf 	
integer :: N		
	!Input variables; values inputted by user
double precision :: Integral
	!Output Variable
double precision :: Int(N/4), xh, kz 			
	!Intermediate variables
double precision :: f
	!Declare function f = f(x)=sinx
integer i,j 
	!counters for loop and index of I(*)
j = 0
i = 0

xh = Abs(xf-x0)/N
kz = (2*xh)/45

Do While (i < (N-4))
	i = 4*j

	Int(j+1) = kz*( 7*f(i*xh) + 32*f((i+1)*xh) + 12*f((i+2)*xh) + 32*f((i+3)*xh) + 7*f((i+4)*xh) )
	j=j+1
end do

Integral=Sum(Int)





return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function f(x) used to define a generic function f(x), using a for x to avoid confusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function f(a)
implicit none
double precision, intent(in) :: a

f = sin(a) ! define generic f(x) here

return
end function f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function f2p(x) used to define exact value for 2nd derivative of a generic function f(x), using b for x to
! avoid confusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!double precision function f2p(b)
!implicit none
!double precision, intent(in) :: b
!
!f2p = 2*cos(b) - b*sin(b) ! define 2nd derivative of f(x) here
!
!return
!end function f2p

