      program diff_sine4
!
!  a 4th program to compare numerical differentiation
!  (somehwat) better documentation
!  compare asymmetric 2-pt formula against 
!  symmetric 3-pt formula
!  CWJ --SDSU - August 2006
!  

      implicit none

      real x    ! argument (in radians)
      real dx   ! delta x for numerical derivative 
      real dfdx  ! derivative of sine

      integer i  ! counter
      
100   continue
      print*,' enter x '
      read*,x

 
      dx = 1.0

      do i = 1,6
        dx = dx / sqrt(10.)

!........... 2 -point formula      
         dfdx = (sin(x+dx)-sin(x))/dx
         print*,dfdx,cos(x)
         write(8,*)dx,abs(dfdx-cos(x))
!............3 - point formula
         dfdx = (sin(x+dx) - sin(x-dx))/(2*dx)
         write(10,*)dx,abs(dfdx-cos(x))

         write(9,*)x,cos(x)
      enddo
      print*,' Approx written to fort.8 '
      print*,' Exact written to fort.9 '
      print*,' Improved approx to fort.10 '

      end

