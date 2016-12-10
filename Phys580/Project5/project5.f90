!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!Project 5
!
!This program uses Runge-Kutta integration (RK4) to plot the locations
! of various orbital configurations.  The output is saved to 'locations.txt' 
! which contains the time in the first column, the x1, y1 positions in the 
! 2nd and 3rd, the x2, y2 in the 4th and 5th and the total energy of the 
! system in the 6th column. Use 'proj5.script' in conjuction with various
! input files for various initial configurations.  'proj5basic1.input' contains the 
! initial values for the setup of the Earth orbiting the Sun with the veclocity
! that makes this orbit circular.  'proj5basic2.input' contains values for the 
! normal elliptical orbit of the Earth about the Sun.  'proj5advanced1.input'
! contains values for two masses orbiting a primary. 'proj5advanced2.input'
! contains values for the Moon orbiting the Earth. 'locations.txt' is to be plotted
! using a gnuplot script correlating to each input file.  'plotb1.p' plots 'locations.txt'
! when created using 'proj5basic1.input'; 'plotb2.p' plots 'locations.txt' as created
! from 'proj5basic2.input'; 'plota1.p' plots 'locations.txt' as created from 
! 'proj5advanced1.input'; 'plota2.p' plots 'locations.txt' as created from
! 'proj5advanced2.input'.  
!
! 	Input Variables:
!		-M: primary mass 
!		-m1, m2: mass of planets 1, 2
!		-x1, y1, x2, y2: initial positions of planets 1, 2
!		-vx1, vy1, vx2, vy2: initial velocities of planets 1, 2
!		-dt: time step
!		-N: number of steps
! 
!	Intermediate variables: 
!		-xyv(8), xyvstepped(8): arrays of positions/velocites
!		-k1(8), k2(8), k3(8), k4(8): k arrays for use in RK4
!		-E: energy value returned by function 'Energy' 
!		-t: time variable
!
! 	Subroutines used: 
!		-'getks(dt, M, m1, m2, xyv, k1, k2, k3, k4)	'	
!			-creates the k arrays, k1, k2, k3, k4, to be used in the RK4 routine
!		-'RKstep(xyv, k1, k2, k3, k4, xyvstepped)'	
!			-Takes the current array of position/velocities (n) and returns the array of 
!	  		 next values (n+1).  Uses the k arrays in RK4 to do this
!
!	Functions used:
!		-'Energy(M, m1, m2, xyv)'
!			-returns the value of the total energy, which should be conserved
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program project5
 implicit none
 
 real :: M, m1, m2, x1, y1, x2, y2, vx1, vy1, vx2, vy2, dt
 integer :: N 
 	!input variables, given by user
 
 real :: xyv(8), xyvstepped(8), k1(8), k2(8), k3(8), k4(8), E, Energy, t
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

	t = 0
	
	xyv(1) = x1
	xyv(2) = vx1
	xyv(3) = y1
	xyv(4) = vy1
	xyv(5) = x2
	xyv(6) = vx2
	xyv(7) = y2
	xyv(8) = vy2
	
	E = Energy(M, m1, m2, xyv)	
	
	open(unit = 1, file = 'locations.txt')
	write(1,*), '# Adrien Atallah'
	write(1,*), '# Project 5'
	write(1,*), '# Contains time vs. location for planet 1 and planet 2 plus the energy of the system'
	write(1,*), '# t, x1, y1, x2, y2, Energy'
	write(1,*), t, x1, y1, x2, y2, E
		
	do i = 1, N	
	
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
!'fxt(M, m1, m2, xyv, fxt0)'	
!	-Takes in the array of positions and velocites (xyv) and returns an array of
!	 the derivatives (fxt0) as given by the equations of motion.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fxt(M, m1, m2, xyv, fxt0)
 implicit none
 
 real :: G !Declare gravitational constant as a parameter
 parameter (G = 6.673E-11) !units are in m^3*kg^-1*s^-2

 real, intent(in) :: M, m1, m2, xyv(8)
 	!input variables
 
 integer :: i
 real :: r1, r2, r12
 	!intermediate variables
 		
 real, intent(out) :: fxt0(8)
 	!output variables 


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

 do i = 1, 8 !initialize fxt array to 0
	fxt0(i) = 0
enddo  

fxt0(1) = xyv(2)!create array for the right hand side to use in RK4
fxt0(2) = ( (-G*M) / (r1**3.0) )*xyv(1) - ( (G*m2) / (r12**3.0) )*(xyv(1) - xyv(5))
fxt0(3) = xyv(4)
fxt0(4) = ( (-G*M) / (r1**3.0) )*xyv(3) - ( (G*m2) / (r12**3.0) )*(xyv(3) - xyv(7))
fxt0(5) = xyv(6)
fxt0(6) = ( (-G*M) / (r2**3.0) )*xyv(5) - ( (G*m1) / (r12**3.0) )*(xyv(5) - xyv(1))
fxt0(7) = xyv(8)
fxt0(8) = ( (-G*M) / (r2**3.0) )*xyv(7) - ( (G*m1) / (r12**3.0) )*(xyv(7) - xyv(3))

return
end	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'getks(dt, M, m1, m2, xyv, k1, k2, k3, k4)'	
!	-creates the k arrays, k1, k2, k3, k4, to be used in the RK4 routine. 
!
! 	Subroutines used: 
!		-'fxt(M, m1, m2, xyv, fxt1)'	
!			-Returns the array of the derivatives of the positions and velocites
!			 as given by the equations of motion.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getks(dt, M, m1, m2, xyv, k1, k2, k3, k4)	
 implicit none 
 
 real, intent(in) :: M, m1, m2, xyv(8)
 real, intent(in) :: dt
 	!input variables
 
 integer :: i
 real :: xyv_n1(8), xyv_n2(8), xyv_n3(8), fxt1(8), fxt2(8), fxt3(8), fxt4(8)
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
 		xyv_n1(i) = 0
 		xyv_n2(i) = 0
 		xyv_n3(i) = 0
 enddo 
 
 do i = 1, 8
 	xyv_n1(i) = xyv(i)
 	xyv_n2(i) = xyv(i)
 	xyv_n3(i) = xyv(i)
 enddo
 
 call fxt(M, m1, m2, xyv, fxt1) !get fxt1 array to use for getting k's
	
 do i = 1, 8 
	k1(i) = dt * fxt1(i)
	
	xyv_n1(i) = xyv(i) + 0.5*k1(i)
 enddo 		
 

! xyv_n1(1) = xyv(1) + 0.5*k1(1) !k1(1) = vx1 = dx1/dt => k1(1) keeps track of how x changes  	
! xyv_n1(2) = xyv(2) + 0.5*k1(2) !k1(2) = dvx1/dt => k1(2) keeps track of the change in vx 	
! xyv_n1(3) = xyv(3) + 0.5*k1(3) !k1(3) = dy1/dt 
! xyv_n1(4) = xyv(4) + 0.5*k1(4) !k1(4) = dvy1/dt 
! xyv_n1(5) = xyv(5) + 0.5*k1(5) !k1(5) = dx2/dt 
! xyv_n1(6) = xyv(6) + 0.5*k1(6) !k1(6) = dvx2/dt
! xyv_n1(7) = xyv(7) + 0.5*k1(7) !k1(7) = dy2/dt 
! xyv_n1(8) = xyv(8) + 0.5*k1(8) !k1(8) = dvy2/dt
 	! !step taken using k1, now call fxt to update array and get fxt2 for next step

 call fxt(M, m1, m2, xyv_n1, fxt2)
 
 do i = 1, 8 
	k2(i) = dt * fxt2(i)
		
	xyv_n2(i) = xyv(i) + 0.5*k2(i)
 enddo 
 	!step taken using k2, now call fxt to update array and get fxt3 for next step
 		
 call fxt(M, m1, m2, xyv_n2, fxt3)
 	
 do i = 1, 8 
	k3(i) = dt * fxt3(i)
		
	xyv_n3(i) = xyv(i) + k3(i)
 enddo  
 	!step taken using k3, now call fxt to update array and get fxt4
 		
 call fxt(M, m1, m2, xyv_n3, fxt4)
	
 do i = 1, 8 
	k4(i) = dt * fxt4(i)
	
 enddo 
	!now k1, k2, k3, k4 arrays have all been created 
 
return
end	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'RKstep(xyv_n, k1, k2, k3, k4, xyv_np1)'	
!	-Takes the current array of position/velocities (n) and returns the array of 
!	  next values (n+1).  Uses the k arrays in RK4 to do this
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
 	 
r1 = sqrt( xyv(1)**2 + xyv(3)**2 ) 
r2 = sqrt( xyv(5)**2 + xyv(7)**2 )
r12 = sqrt( (xyv(1) - xyv(5))**2 + (xyv(3) - xyv(7))**2 )

Energy = (m1/2.0) * (xyv(2)**2 + xyv(4)**2) + (m2/2.0)*(xyv(6)**2 + xyv(8)**2)&
 & - ( (G*M*m1)/(r1) ) - ( (G*M*m2)/(r2) ) - ( (G*m1*m2)/(r12) )

 return
 end 


