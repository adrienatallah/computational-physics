program project1
implicit none
real :: x, dx 			!input variables
real :: f2, f5 			!numerically calculated 2nd derivatives of f(x); f2 is for 2pt formula, f5 is for 5pt
real :: fp			!variable to store exact formula for 2nd derivative of f(x) 
real :: fx, fxph, fxmh, fxp2h, fxm2h 		!fx = f(x), fxph=f(x+h), fxmh=f(x-h), fxp2h=f(x+2h), fxm2h=f(x-2h)
character :: continue = 'y'

do while (continue == 'y' .or. continue == 'Y')

	print *, 'This is Computational Project 1: Numerical Derivatives'
	print *, 'Written by: Adrien Atallah'
	print *, ' '


	print *, 'Compute the 2nd derivative of f(x)=xsinx'
	print *, 'Enter values of x and dx'
	print *, ' '
	read *, x, dx
	
	call xsinx(x,dx,fp,f2,f5)  
	
	print *, 'f"(',x,') = ', fp
	print *, 'Three Point Formula = ', f2
	print *, 'Five Point Formula = ', f5

	
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

subroutine xsinx(x,dx,fp,f2,f5)
implicit none

real :: x, dx 			!input variables
real :: f2, f5 			!numerically calculated 2nd derivatives of f(x); f2 is for 2pt formula, f5 is for 5pt
real :: fp			!variable to store exact formula for 2nd derivative of f(x) 
real :: fx, fxph, fxmh, fxp2h, fxm2h 		!fx = f(x), fxph=f(x+h), fxmh=f(x-h), fxp2h=f(x+2h), fxm2h=f(x-2h)

fx = x*sin(x)
fxph = (x+dx)*sin(x+dx)
fxmh = (x-dx)*sin(x-dx)
fxp2h = (x+2*dx)*sin(x+2*dx)
fxm2h = (x-2*dx)*sin(x-2*dx)

fp = 2*cos(x) - x*sin(x)
f2 = (-2*fx + fxph + fxmh)/(dx**2)
f5 = (-fxm2h + 16*fxmh - 30*fx + 16*fxph - fxp2h)/(12*dx**2)

return
end

