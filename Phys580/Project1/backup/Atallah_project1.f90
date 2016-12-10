!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 1
!
!This program uses the 3pt and 5pt formulas to numerically approximate the second derivative of a function 
!of x with respect to x.  It will generate two data files; fort.8 contains the values of the difference
!between the 3pt approximation and the exact values in one column and the stepsize in the other, fort.9 
!contains the values of the difference between the 5pt approximation and the exact values in one column 
!and the stepsize in the other.  The purpose of generating these two data files will be to plot make a plot
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
!	'df(x,dx,fp,f3,f5)' - calculates f''(x) for user inputed x and dx values with both 3pt and 5pt
!			      formulas and also generates fort.8 and fort.9
!			    - input variables:
!				  -'x' value for which the the derivatives of f(x) will be evaluated at;
!				   inputed by user
!				  -'dx' stepsize; inputed by user
!			    -output variables:
!				  -'fp' exact 2 derivatived as calculated by function 'f2p(b)' with user
!				   inputed value of x
!				  -'f3' numerically approximated 2nd derivative using 3pt formula
!				  -'f5' numerically approximated 2nd derivative using 5pt formula
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program project1
implicit none
double precision :: x, dx 			
	!input variables
double precision :: f3, f5 			
	!numerically calculated 2nd derivatives of f(x); f3 is for 3pt formula, f5 is for 5pt
double precision :: fp			
	!variable to store exact formula for 2nd derivative of f(x) 
double precision :: fx, fxph, fxmh, fxp2h, fxm2h 		
	!fx = f(x), fxph=f(x+h), fxmh=f(x-h), fxp2h=f(x+2h), fxm2h=f(x-2h)
character :: continue = 'y'

do while (continue == 'y' .or. continue == 'Y')

	print *, 'This is Computational Project 1: Numerical Derivatives'
	print *, 'Written by: Adrien Atallah'
	print *, ' '


	print *, 'Compute the 2nd derivative of f(x)=xsinx'
	print *, 'Enter values of x and dx'
	print *, ' '
	read *, x, dx
	
	call df(x,dx,fp,f3,f5)  ! Call subroutine to calculate and approximate 2nd derivatives of a function of x

	
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

end program project1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine df(x,dx,fp,f3,f5)
implicit none

double precision :: x, dx 			
	!input variables
double precision :: f3, f5 			
	!numerically calculated 2nd derivatives of f(x); f3 is for 3pt formula, f5 is for 5pt
double precision :: fp			
	!variable to store exact formula for 2nd derivative of f(x) 
double precision :: fx, fxph, fxmh, fxp2h, fxm2h 		
	!fx = f(x), fxph=f(x+h), fxmh=f(x-h), fxp2h=f(x+2h), fxm2h=f(x-2h)
integer i

			!fx = x*sin(x)
fxph = f(x+dx)		!(x+dx)*sin(x+dx)
fxmh = f(x-dx)		!(x-dx)*sin(x-dx)
fxp2h = f(x+2*dx)	!(x+2*dx)*sin(x+2*dx)
fxm2h = f(x-2*dx)	!(x-2*dx)*sin(x-2*dx)

fp = 2*cos(x) - x*sin(x)
f3 = (-2*fx + fxph + fxmh)/(dx**2)
f5 = (-fxm2h + 16*fxmh - 30*fx + 16*fxph - fxp2h)/(12*dx**2)

print *, 'f"(',x,') = ', fp
print *, 'Three Point Formula = ', f3
print *, 'Five Point Formula = ', f5


! loop to generate data files in order to create plot with gnuplot
! after generating fort.8 and fort.9, load gnuplot script 'plot.p' to create plot as postscript file
do i = 1,7        

				!fx = x*sin(x)
	fxph = f(x+dx)		!(x+dx)*sin(x+dx)
	fxmh = f(x-dx)		!(x-dx)*sin(x-dx)
	fxp2h = f(x+2*dx)	!(x+2*dx)*sin(x+2*dx)
	fxm2h = f(x-2*dx)	!(x-2*dx)*sin(x-2*dx)

	fp = 2*cos(x) - x*sin(x)
	f3 = (-2*fx + fxph + fxmh)/(dx**2)
	f5 = (-fxm2h + 16*fxmh - 30*fx + 16*fxph - fxp2h)/(12*dx**2)

        write(8,*)dx,abs(f3-fp) !write difference for 3pt vs dx to fort.8

        write(9,*)dx,abs(f5-fp) !write difference for 5pt vs dx to fort.9

	dx=dx/3

enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function f(x) used to define a generic function f(x), using a for x to avoid confusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function f(a)

f = a*sin(a) ! define generic f(x) here

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function f2p(x) used to define exact value for 2nd derivative of a generic function f(x), using b for x to
! avoid confusion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function f2p(b)

f2p = 2*cos(b) - b*sin(b) ! define 2nd derivative of f(x) here

return
end

