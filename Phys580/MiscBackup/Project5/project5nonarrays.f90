!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 5
!
!This program 
!
! 	Input Variables:
!		-R: effective nuclear radius
!		-L: size of the box
!		-N: number of lattice points
! 
!	Intermediate variables: 
!		-dx: stepsize
!		-PEmat: Potential Energy matrix
!		-KEmat: Kinetic Energy matrix
!		-Hmat: Kinetec + Potential matrix
!		-KEeVal, PEeVal: KE and PE array of eigenvalues
!		-HeVal: array of eigenvalues of KE + PE matrix
!		-KEeVec, PEeVec: KE and PE matrix of eigenvectors
!		-HeVec: KE + PE eigenvectors matrix
!		!HeValWS, HeVecWS: eigenvalues/eigenvectors for Woods-Saxon Potential
!
! 	Subroutines used: 
!		-'KEmatrix(N, hbar, m, dx, KEmat)'	
!			-Creates KE matrix used for all scenarios
!		-'ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)'	
!			-Takes in KE matrix, calcuates eigenvectors/eigenvalues and creates
!			 text files to be plotted for the particle in a box scenario
!		-'HarmonicOscillator(N, KEmat, dx, L, hbar, m, Hmat, HeVal, HeVec)'
!			-Takes in KE matrix, Initializes PE matrix, calcuates eigenvectors/eigenvalues of
!			 their sum, then creates text files to be plotted for the harmonic oscillator scenario
!		-'WoodsSaxon(N, KEmat, dx, L, a, R, V0, HmatWS, HeValWS, HeVecWS)'
!			-Takes in KE matrix, Initializes Woods-Saxon PE matrix, calcuates eigenvectors/eigenvalues of
!			 their sum, then creates text files to be plotted for the Woods-Saxon potential
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program project5
 implicit none
 
 real :: M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2, dt, t
 integer :: N 
 	!input variables, given by user
 
 real :: fxt0(8), xyv(8), xyvstepped(8), k1(8), k2(8), k3(8), k4(8), E, Energy
 integer :: i, j
 	!Intermediate variables
    
 character :: continue = 'y' 
 
 
 do while (continue == 'y' .or. continue == 'Y')    


 	print *, ' '
	print *, 'Computational Project 5:'
	print *, 'Orbital Mechanics'
	print *, 'Written by: Adrien Atallah'
	print *, ' '

	print *, ' '
	print *, 'Enter the primary mass, M (kg):'
	print *, ' '
	read *, M
	
	print *, ' '
	print *, 'Enter the mass of the first orbiting body, m1 (kg):'
	print *, ' '
	read *, m1
	
	print *, ' '
	print *, 'Enter the mass of the second orbiting body, m2 (kg):'
	print *, ' '
	read *, m2
	
	print *, ' '
	print *, 'Enter the coordinates for the initial position of the first orbiting body, x1 and y1 (meters):'
	print *, ' '
	read*, x1, y1
	
	print *, ' '
	print *, 'Enter the coordinates for the initial position of the second orbiting body, x2 and y2 (meters):'
	print *, ' '
	read*, x2, y2
	
	print *, ' '
	print *, 'Enter the x and y components of the initial velocity of the first orbiting body (m/s):'
	print *, ' '
	read*, vx1, vy1
	
	print *, ' '
	print *, 'Enter the x and y components of the initial velocity of the second orbiting body (m/s):'
	print *, ' '
	read*, vx2, vy2
	
	print *, ' '
	print *, 'Enter the size of the time step, dt (seconds):'
	print *, ' '
	read *, dt
	
	print *, ' '
	print *, 'Enter the number of steps to take (integer):'
	print *, ' '
	read *, N

	E = 0
	t = 0
	
	open(unit = 1, file = 'locations.txt')
	write(1,*), '# Adrien Atallah'
	write(1,*), '# Project 5'
	write(1,*), '# Contains time vs. location for planet 1 and planet 2 plus the energy of the system'
	write(1,*), '# t, x1, y1, x2, y2, Energy'
	write(1,*), t, x1, y1, x2, y2, E
	
	!call fxt(M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2, fxt0)
	
	xyv(1) = x1
	xyv(2) = vx1
	xyv(3) = y1
	xyv(4) = vy1
	xyv(5) = x2
	xyv(6) = vx2
	xyv(7) = y2
	xyv(8) = vy2
	
	!call getks(dt, M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2, k1, k2, k3, k4)	
	!call RKstep(fxt0, k1, k2, k3, k4, xyv) !get first step
	
	do i = 1, N	
	
		!print *, 'xyv(1) =', xyv(1)
	
		call getks(dt, M, m1, m2, xyv, k1, k2, k3, k4)	
		
		E = Energy(M, m1, m2, xyv)
		
		call RKstep(xyv, k1, k2, k3, k4, xyvstepped)
		t = t + dt !after calling rkstep t increases by dt
		
		write (1 , *), t, xyvstepped(1), xyvstepped(3), xyvstepped(5), xyvstepped(7), E		
		
		do j = 1, 8
			xyv(j) = xyvstepped(j) !so previous step is fed back into RKstep
		enddo
		
	enddo
	
	close(1)
	
	print *, ' '
	print *, 'Planets location data stored in "locations.txt" '
	print *, ' '

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
!'ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)'	
!	-Takes in KE matrix, calcuates eigenvectors/eigenvalues and creates
!	  text files to be plotted for the particle in a box scenario. 
!
! 	Subroutines used: 
!		-'eig(KEmat, N, N, KEeVal, KEeVec, work)'	
!			-Returns NxN matrix of eigenvectors and N-dimensional array
!			 of eigenvalues of the KE matrix
!		-'eigsrt(KEeVal, KEeVec, N, N)'	
!			-Sorts eigenvalues in ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fxt(M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2, fxt0)
 implicit none
 
 real :: G !Declare gravitational constant as a parameter
 parameter (G = 6.673E-11) !units are in m^3*kg^-1*s^-2

 real, intent(in) :: M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2
 	!input variables
 
 integer :: i
 real :: r1, r2, r12
 	!intermediate variables
 		
 real, intent(out) :: fxt0(8)
 	!output variables 


r1 = sqrt( x1**2 + y1**2 ) !create r values for use in fxt
r2 = sqrt( x2**2 + y2**2 )
r12 = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )

 do i = 1, 8 !initialize fxt array to 0
	fxt0(i) = 0
enddo  

fxt0(1) = vx1!create array for the right hand side to use in RK4
fxt0(2) = ( (-G*M) / (r1**3.0) )*x1 - ( (G*m2) / (r12**3.0) )*(x1 - x2)
fxt0(3) = vy1
fxt0(4) = ( (-G*M) / (r1**3.0) )*y1 - ( (G*m2) / (r12**3.0) )*(y1 - y2)
fxt0(5) = vx2
fxt0(6) = ( (-G*M) / (r2**3.0) )*x2 - ( (G*m1) / (r12**3.0) )*(x2 - x1)
fxt0(7) = vy2
fxt0(8) = ( (-G*M) / (r2**3.0) )*y2 - ( (G*m1) / (r12**3.0) )*(y2 - y1)

  do i = 1, 8 
	! print*, 'fxt(', i, ')=', fxt0(i)
 enddo
 
  print *, ' '

return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)'	
!	-Takes in KE matrix, calcuates eigenvectors/eigenvalues and creates
!	  text files to be plotted for the particle in a box scenario. 
!
! 	Subroutines used: 
!		-'eig(KEmat, N, N, KEeVal, KEeVec, work)'	
!			-Returns NxN matrix of eigenvectors and N-dimensional array
!			 of eigenvalues of the KE matrix
!		-'eigsrt(KEeVal, KEeVec, N, N)'	
!			-Sorts eigenvalues in ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getks(dt, M, m1, m2, xyv, k1, k2, k3, k4)	
 implicit none 
 
 real, intent(in) :: M, m1, m2, xyv(8)
 real, intent(in) :: dt
 	!input variables
 
 integer :: i
 real :: x1temp, y1temp, x2temp, y2temp, vx1temp, vy1temp, vx2temp, vy2temp
 real :: fxt1(8), fxt2(8), fxt3(8), fxt4(8)
 	!intermediate variables
 		
 real, intent(out) :: k1(8), k2(8), k3(8), k4(8)
 	!output variables 
		
  do i = 1, 8 !initialize all arrays to 0
 		k1(i) = 0
 		k2(i) = 0
 		k3(i) = 0
 		k4(i) = 0
 		fxt1(i) = 0
 		fxt2(i) = 0
 		fxt3(i) = 0
 		fxt4(i) = 0
 enddo 
 
 x1temp = xyv(1)
 vx1temp = xyv(2)
 y1temp = xyv(3)
 vy1temp = xyv(4) 
 x2temp = xyv(5) 
 vx2temp = xyv(6)
 y2temp = xyv(7)
 vy2temp = xyv(8)
 
 ! print *, 'x1temp (for k1) =', x1temp
 ! print *, 'vx1temp (for k1) =', vx1temp
 ! print *, ' '
 
 	
 call fxt(M, m1, m2, x1temp, y1temp, x2temp, y2temp, vx1temp, vy1temp, vx2temp, vy2temp, fxt1) !get fxt1 array to use for getting k's
	
 do i = 1, 8 
	k1(i) = dt * fxt1(i)
	print *, 'k1(', i ,') =', k1(i)
 enddo 		
 

 x1temp = xyv(1) + 0.5*k1(1)
 	!k1(1) = vx1 = dx1/dt => k1(1) keeps track of how x changes 
 vx1temp = xyv(2) + 0.5*k1(2) 
 	!k1(2) = dvx1/dt => k1(2) keeps track of the change in vx
 y1temp = xyv(3) + 0.5*k1(3) !k1(3) = dy1/dt 
 vy1temp = xyv(4) + 0.5*k1(4) !k1(4) = dvy1/dt 
 x2temp = xyv(5) + 0.5*k1(5) !k1(5) = dx2/dt 
 vx2temp = xyv(6) + 0.5*k1(6) !k1(6) = dvx2/dt
 y2temp = xyv(7) + 0.5*k1(7) !k1(7) = dy2/dt 
 vy2temp = xyv(8) + 0.5*k1(8) !k1(8) = dvy2/dt
 	!step taken using k1, now call fxt to update array and get fxt2 for next step
 
 ! print *, 'vx1temp = xyv(2) + 0.5*k1(2)'
 ! print *, ' '		
  ! print *, 'xyv(2)=', xyv(2)
 ! print *, 'x1temp (for k2) =', x1temp
  ! print *, 'vx1temp (for k2) =', vx1temp
   ! print *, 'y1temp (for k2) =', y1temp
 		
 call fxt(M, m1, m2, x1temp, y1temp, x2temp, y2temp, vx1temp, vy1temp, vx2temp, vy2temp, fxt2)
 
  x1temp = 0
 vx1temp = 0
 y1temp = 0
 vy1temp = 0 
 x2temp = 0
 vx2temp = 0
 y2temp = 0
 vy2temp = 0
 
 do i = 1, 8 
	k2(i) = dt * fxt2(i)
	print *, 'k2(', i ,') =', k2(i)
 enddo 
 
 x1temp = xyv(1) + 0.5*k2(1)
 	!k1(1) = vx1 = dx1/dt => k1(1) keeps track of how x changes 
 vx1temp = xyv(2) + 0.5*k2(2) 
 	!k1(2) = dvx1/dt => k1(2) keeps track of the change in vx
 y1temp = xyv(3) + 0.5*k2(3) !k1(3) = dy1/dt 
 vy1temp = xyv(4) + 0.5*k2(4) !k1(4) = dvy1/dt 
 x2temp = xyv(5) + 0.5*k2(5) !k1(5) = dx2/dt 
 vx2temp = xyv(6) + 0.5*k2(6) !k1(6) = dvx2/dt
 y2temp = xyv(7) + 0.5*k2(7) !k1(7) = dy2/dt
 vy2temp = xyv(8) + 0.5*k2(8) !k1(8) = dvy2/dt
 	!step taken using k2, now call fxt to update array and get fxt3 for next step
 call fxt(M, m1, m2, x1temp, y1temp, x2temp, y2temp, vx1temp, vy1temp, vx2temp, vy2temp, fxt3)
 
  x1temp = 0
 vx1temp = 0
 y1temp = 0
 vy1temp = 0 
 x2temp = 0
 vx2temp = 0
 y2temp = 0
 vy2temp = 0
 
	
 do i = 1, 8 
	k3(i) = dt * fxt3(i)
	print *, 'k3(', i ,') =', k3(i)
 enddo 		
 
 x1temp = xyv(1) + k3(1)
 	!k1(1) = vx1 = dx1/dt => k1(1) keeps track of how x changes 
 vx1temp = xyv(2) + k3(2) 
 	!k1(2) = dvx1/dt => k1(2) keeps track of the change in vx
 y1temp = xyv(3) + k3(3) !k1(3) = dy1/dt 
 vy1temp = xyv(4) + k3(4) !k1(4) = dvy1/dt 
 x2temp = xyv(5) + k3(5) !k1(5) = dx2/dt 
 vx2temp = xyv(6) + k3(6) !k1(6) = dvx2/dt
 y2temp = xyv(7) + k3(7) !k1(7) = dy2/dt 
 vy2temp = xyv(8) + k3(8) !k1(8) = dvy2/dt
 	!step taken using k3, now call fxt to update array and get fxt4
 		
 call fxt(M, m1, m2, x1temp, y1temp, x2temp, y2temp, vx1temp, vy1temp, vx2temp, vy2temp, fxt4)
	
 do i = 1, 8 
	k4(i) = dt * fxt4(i)
	print *, 'k4(', i ,') =', k4(i)
 enddo 
	!now k1, k2, k3, k4 arrays have all been created 
 
return
end	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)'	
!	-Takes in KE matrix, calcuates eigenvectors/eigenvalues and creates
!	  text files to be plotted for the particle in a box scenario. 
!
! 	Subroutines used: 
!		-'eig(KEmat, N, N, KEeVal, KEeVec, work)'	
!			-Returns NxN matrix of eigenvectors and N-dimensional array
!			 of eigenvalues of the KE matrix
!		-'eigsrt(KEeVal, KEeVec, N, N)'	
!			-Sorts eigenvalues in ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RKstep(xyv_n, k1, k2, k3, k4, xyv_np1)
 implicit none

 real, intent(in) :: xyv_n(8), k1(8), k2(8), k3(8), k4(8)
 	!input variables
 
 integer :: i
 	!intermediate variables
 		
 real, intent(out) :: xyv_np1(8)
 	!output variables 

 do i = 1, 8 !initialize all arrays to 0
 	xyv_np1(i) = 0
 enddo    

 do i = 1, 8
 	xyv_np1(i) = xyv_n(i) + (1.0/6.0) * ( k1(i) + 2.0*k2(i) + 2.0*k3(i) + k4(i) ) !RK4 step
 	!print *, 'xyv_np1(', i, ') =', xyv_np1(i)
 enddo
 
return
end	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy(M, m1, m2, xyv) returns the total energy, which should be conserved
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function Energy(M, m1, m2, xyv)
 implicit none
 
 real, intent(in) :: M, m1, m2, xyv(8)
 	!input variables
 
 real :: G !Declare gravitational constant as a parameter
 parameter (G = 6.673E-11) !units are in m^3*kg^-1*s^-2
 real :: r1, r2, r12
 	!intermediate variables
 		
 !xyv(1) = x1
 !xyv(2) = vx1
 !xyv(3) = y1
 !xyv(4) = vy1
 !xyv(5) = x2
 !xyv(6) = vx2
 !xyv(7) = y2
 !xyv(8) = vy2
 	 
r1 = sqrt( xyv(1)**2 + xyv(3)**2 ) !create r values for use in fxt
r2 = sqrt( xyv(5)**2 + xyv(7)**2 )
r12 = sqrt( (xyv(1) - xyv(5))**2 + (xyv(3) - xyv(7))**2 )

Energy = (m1/2.0) * (xyv(2)**2 + xyv(4)**2) + (m2/2.0)*(xyv(6)**2 + xyv(8)**2)&
 & - ( (G*M*m1)/(r1) ) - ( (G*M*m2)/(r2) ) - ( (G*m1*m2)/(r12) )

 return
 end 


