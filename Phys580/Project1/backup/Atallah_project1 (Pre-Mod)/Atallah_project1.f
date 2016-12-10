program project1
implicit none
real :: x, dx 			!input variables
real :: f3, f5 			!numerically calculated 2nd derivatives of f(x); f3 is for 3pt formula, f5 is for 5pt
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
	
	call xsinx(x,dx,fp,f3,f5)  ! Call subroutine to calculate 2nd derivatives of xsinx

	
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

subroutine xsinx(x,dx,fp,f3,f5)
implicit none

real :: x, dx 			!input variables
real :: f3, f5 			!numerically calculated 2nd derivatives of f(x); f3 is for 3pt formula, f5 is for 5pt
real :: fp			!variable to store exact formula for 2nd derivative of f(x) 
real :: fx, fxph, fxmh, fxp2h, fxm2h 		!fx = f(x), fxph=f(x+h), fxmh=f(x-h), fxp2h=f(x+2h), fxm2h=f(x-2h)
integer i

fx = x*sin(x)
fxph = (x+dx)*sin(x+dx)
fxmh = (x-dx)*sin(x-dx)
fxp2h = (x+2*dx)*sin(x+2*dx)
fxm2h = (x-2*dx)*sin(x-2*dx)

fp = 2*cos(x) - x*sin(x)
f3 = (-2*fx + fxph + fxmh)/(dx**2)
f5 = (-fxm2h + 16*fxmh - 30*fx + 16*fxph - fxp2h)/(12*dx**2)

print *, 'f"(',x,') = ', fp
print *, 'Three Point Formula = ', f3
print *, 'Five Point Formula = ', f5


! loop to generate data files in order to create plot with gnuplot
! after generating fort.8 and fort.9, load gnuplot script 'plot.p' to create plot as postscript file
do i = 1,7        

	fx = x*sin(x)
	fxph = (x+dx)*sin(x+dx)
	fxmh = (x-dx)*sin(x-dx)
	fxp2h = (x+2*dx)*sin(x+2*dx)
	fxm2h = (x-2*dx)*sin(x-2*dx)

	fp = 2*cos(x) - x*sin(x)
	f3 = (-2*fx + fxph + fxmh)/(dx**2)
	f5 = (-fxm2h + 16*fxmh - 30*fx + 16*fxph - fxp2h)/(12*dx**2)

        write(8,*)dx,abs(f3-fp) !write difference for 3pt vs dx to fort.8

        write(9,*)dx,abs(f5-fp) !write difference for 5pt vs dx to fort.9

	dx=dx/3

enddo

return
end

