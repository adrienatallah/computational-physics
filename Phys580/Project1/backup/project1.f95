program project1
implicit none
real :: x, dx 			!input variables
real :: f2, f5 			!numerically calculated 2nd derivatives of f(x); f2 is for 2pt formula, f5 is for 5pt
real :: fp			!variable to store exact formula for 2nd derivative of f(x) 
real :: fx, fxph, fxmh 		!fx = f(x), fxph=f(x+h), fxmh=f(x-h)
character :: continue = 'y'

do while (continue == 'y' .or. continue == 'Y')

	print *, 'Project 1: Numerical Derivatives'
	print *, 'Adrien Atallah'
	print *, ' '


	print *, 'Compute the 2nd derivative of f(x)=xsinx'
	print *, 'Enter values of x and dx'
	print *, ' '
	read *, x, dx
	
	fx = x*sin(x)
	fxph = (x+dx)*sin(x+dx)
	fxmh = (x-dx)*sin(x-dx)

	fp = 2*cos(x) - x*sin(x)
	f2 = (-2*fx + fxph + fxmh)/(dx**2)

	!print *, 'x=', x 
	!print *, 'dx=', dx
	!print *, 'fx=', fx
	!print *, 'fxph=', fxph
	!print *, 'fxmh=', fxmh
	print *, 'f"(',x,') = ', fp
	print *, 'Three Point Formula = ', f2


	
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
