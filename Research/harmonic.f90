!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Adrien Atallah
!
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
!		-'HarmonicOscillator(N, KEmat, dx, L, hbar, m, Hmat, HeVal, HeVec)'
!			-Takes in KE matrix, Initializes PE matrix, calcuates eigenvectors/eigenvalues of
!			 their sum, then creates text files to be plotted for the harmonic oscillator scenario
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program harmonic
 implicit none
 real :: c, hbar, V0, m, a !Constants
 	!Declare constant values as parameters:
 parameter (c = 1)
 parameter (hbar = 197.33) !c = 1 makes hbar = 197.33 MeV-fm = 197.33*10E-15 MeV-m
 parameter (V0 = 50) !units are in MeV
 parameter (m = 939) !units are in MeV
 parameter (a =0.2) !units are in femtometers or fermis, where 1fm = 10E-15m 
 
 real :: R, L
 integer :: N 
 	!input variables
 
 real :: x, xi, dx, PEmat(1000, 1000), KEmat(1000, 1000), KEeVal(1000), KEeVec(1000, 1000)
 real :: work(1000), Hmat(1000, 1000), HeVal(1000), HeVec(1000, 1000), HmatWS(1000, 1000)
 real :: HeValWS(1000), HeVecWS(1000, 1000)
 integer :: i
 	!Intermediate variables
    
 character :: continue = 'y' 
 
 
 do while (continue == 'y' .or. continue == 'Y')    


 	print *, ' '
	print *, ' '
	print *, 'Written by: Adrien Atallah'
	print *, ' '

	print *, ' '
	print *, 'Enter the size of the box, L'
	print *, ' '
	read *, L
	
	print *, ' '
	print *, 'Enter the number of lattice points (less than 1000)'
	print *, ' '
	read *, N
	
	! print *, ' '
	! print *, 'Enter the nuclear radius value (in femtometers or fermis)'
	! print *, ' '
	! read*, R
	
	dx = (2*L)/(N-1)
	
	call KEmatrix(N, hbar, m, dx, KEmat)
	
	!call ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)
	
	call HarmonicOscillator(N, KEmat, dx, L, hbar, m, Hmat, HeVal, HeVec)
	
	!call WoodsSaxon(N, KEmat, dx, L, a, R, V0, HmatWS, HeValWS, HeVecWS)	 

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
!KEmatrix(a, R, V0, x, dx, N, L, PEmat) initializes diagonal matrix for the Kinetic
! energy to be used in all three scenarios
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KEmatrix(N, hbar, m, dx, KEmat)
implicit none

 real, intent(in) :: hbar, m, dx
 integer, intent(in) :: N
 	!input variables
 		
 integer :: i, j 
 real :: k
 	!intermediate variables
 		
 real, intent(out) :: KEmat(N, N)
 	!output variables
 

 do i = 1, N !initialize NxN Potential Energy matrix to 0
 	do j = 1, N
		KEmat(i, j) = 0
	enddo
 enddo    

 k = (hbar)**2/( 2*m*dx**2)

 do i = 1, N !initialize NxN Kinetic Energy matrix 	
	KEmat(i, i) = 2*k	
	KEmat(i+1, i) = -k
	KEmat(i, i+1) = -k
 enddo  
 
		
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
! subroutine ParticleInABox(N, KEmat, dx, L, KEeVal, KEeVec, work)
! implicit none
! 
 ! real, intent(in) :: KEmat, dx, L
 ! integer, intent(in) :: N
 	! !input variables
 ! 
 ! real :: xi
 ! integer :: i, j 
 	! !intermediate variables
 		! 
 ! real, intent(out) :: KEeVal(N), KEeVec(N, N), work(N)
 	! !output variables 
! 
 ! do i = 1, N !initialize all arrays to 0
 	! do j = 1, N
 		! KEeVal(j) = 0
		! KEeVec(i, j) = 0
		! work(j) = 0
	! enddo
 ! enddo    
! 
 ! call eig(KEmat, N, N, KEeVal, KEeVec, work) !np = N
 ! call eigsrt(KEeVal, KEeVec, N, N) !np = N
  ! 
 ! open(unit = 1, file = 'GroundState.txt')
 ! open(unit = 2, file = '1stExState.txt')
 ! 
 ! do i = 1, N !create file with 1st and 2nd columns of eigenvector matrix 	
	! xi = -L + (i - 1)*dx
	! write(1,*) xi, KEeVec(i, 1)
	! write(2,*) xi, KEeVec(i, 2)
 ! enddo  	
 ! 
! return
! end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'HarmonicOscillator(N, KEmat, dx, L, hbar, m, Hmat, HeVal, HeVec)'
!	-Takes in KE matrix, Initializes PE matrix, calcuates eigenvectors/eigenvalues of
!	  their sum, then creates text files to be plotted for the harmonic oscillator scenario
!
! 	Subroutines used: 
!		-'eig(Hmat, N, N, HeVal, HeVec, work)'	
!			-Returns NxN matrix of eigenvectors and N-dimensional array
!			 of eigenvalues of the KE +PE matrix
!		-'eigsrt(HeVal, HeVec, N, N)'	
!			-Sorts eigenvalues in ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HarmonicOscillator(N, KEmat, dx, L, hbar, m, Hmat, HeVal, HeVec)
implicit none

 real, intent(in) :: KEmat(N, N), dx, L, hbar, m
 integer, intent(in) :: N
 	!input variables
 
 real :: xi, Vho, PEmat(N, N), work(N)
 integer :: i, j 
 	!intermediate variables
 		
 real, intent(out) :: Hmat(N, N), HeVal(N), HeVec(N, N)
 	!output variables
 

 do i = 1, N !initialize all arrays to 0
 	do j = 1, N
 		HeVal(j) = 0
		HeVec(i, j) = 0
		PEmat(i, j) = 0
		Hmat(i, j) = 0
		work(j) = 0
	enddo
 enddo    

  do i = 1, N !initialize PE for Harmonic oscillator and H = KE + PE matrices
  	xi = -L + (i - 1)*dx 
  	PEmat(i, i) = Vho(hbar, m, xi) !create diagonalized matrix for PE of Harmonic Oscillator
  	
 	do j = 1, N 		
		Hmat(i, j) = KEmat(i, j) + PEmat(i, j)
	enddo
	
 enddo    

 
 call eig(Hmat, N, N, HeVal, HeVec, work) !np = N
 call eigsrt(HeVal, HeVec, N, N) !np = N
  
 open(unit = 3, file = 'GroundStateHO.txt')
 open(unit = 4, file = '1stExStateHO.txt')
 open(unit = 11, file = '2ndExStateHO.txt')
 open(unit = 12, file = '3rdExStateHO.txt')
 open(unit = 13, file = '4thExStateHO.txt')
 open(unit = 14, file = '5thExStateHO.txt')
 
 do i = 1, N !create file with 1-6 columns of eigenvector matrix 	
	xi = -L + (i - 1)*dx
	write(3,*) xi, HeVec(i, 1) !ground state
	write(4,*) xi, HeVec(i, 2) !1st excited state
	write(11,*) xi, HeVec(i, 3) !2nd excited state
	write(12,*) xi, HeVec(i, 4) !3rd excited state
	write(13,*) xi, HeVec(i, 5) !4th excited state
	write(14,*) xi, HeVec(i, 6) !5th excited state
 enddo  	
 
return
end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!'WoodsSaxon(N, KEmat, dx, L, a, R, V0, HmatWS, HeValWS, HeVecWS)'
!	-Takes in KE matrix, Initializes Woods-Saxon PE matrix, calcuates eigenvectors/eigenvalues of
!	  their sum, then creates text files to be plotted for the Woods-Saxon potential
!
!	Subroutines used: 
!		-'eig(Hmat, N, N, HeVal, HeVec, work)'	
!			-Returns NxN matrix of eigenvectors and N-dimensional array
!			 of eigenvalues of the KE +PE matrix
!		-'eigsrt(HeVal, HeVec, N, N)'	
!			-Sorts eigenvalues in ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine WoodsSaxon(N, KEmat, dx, L, a, R, V0, Hmat, HeVal, HeVec)
! implicit none
! 
 ! real, intent(in) :: KEmat(N, N), dx, L, a, R, V0
 ! integer, intent(in) :: N
 	! !input variables
 ! 
 ! real :: xi, V, PEmat(N, N), work(N)
 ! integer :: i, j 
 	! !intermediate variables
 		! 
 ! real, intent(out) :: Hmat(N, N), HeVal(N), HeVec(N, N)
 	! !output variables
 ! 
! 
 ! do i = 1, N !initialize all arrays to 0
 	! do j = 1, N
 		! HeVal(j) = 0
		! HeVec(i, j) = 0
		! PEmat(i, j) = 0
		! Hmat(i, j) = 0
		! work(j) = 0
	! enddo
 ! enddo    
! 
  ! do i = 1, N !initialize PE for Harmonic oscillator and H = KE + PE matrices
  	! xi = -L + (i - 1)*dx 
  	! PEmat(i, i) = V(a, R, V0, xi) !create diagonalized matrix for Woods-Saxon PE
  	! 
 	! do j = 1, N 		
		! Hmat(i, j) = KEmat(i, j) + PEmat(i, j)
	! enddo
	! 
 ! enddo    
! 
 ! 
 ! call eig(Hmat, N, N, HeVal, HeVec, work) !np = N
 ! call eigsrt(HeVal, HeVec, N, N) !np = N
  ! 
 ! open(unit = 7, file = 'GroundStateWS.txt') !unit #5, 6 reserved
 ! open(unit = 8, file = '1stExStateWS.txt')
 ! open(unit = 9, file = '2ndExStateWS.txt') 
 ! open(unit = 10, file = '3rdExStateWS.txt')
 ! 
 ! do i = 1, N !create file with 1st, 2nd, 3rd and 4th columns of eigenvector matrix 	
	! xi = -L + (i - 1)*dx
	! write(7,*) xi, HeVec(i, 1) !write ground state vs x to 'GroundStateWS.txt'
	! write(8,*) xi, HeVec(i, 2) !write 1st excited state vs x to '1stExStateWS.txt'
	! write(9,*) xi, HeVec(i, 3) !write 2nd excited state vs x to '2ndExStateWS.txt'
	! write(10,*) xi, HeVec(i, 4) !write 3rd excited state vs x to '3rdExStateWS.txt'
 ! enddo  	
 ! 
! return
! end	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function Vho(a, R, V0, x) returns the potential for a harmonic oscillator as a function
!	of x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function Vho(hbar, m, x)
 implicit none
 
 real, intent(in) :: hbar, m, x
 	!input variables
 		
 real :: k
 	!intermediate variables
 		
 k = (hbar**2)/(m)
 
 Vho = 0.5*k*x**2

 return
 end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function V(a, R, V0, x) returns the woods-saxon potential for a certain x value, 
!	using constants a, R, V0 as defined in main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function V(a, R, V0, x)
 implicit none
 
 real, intent(in) :: a, R, V0, x
 	!Input variables
 
V = -V0/(1 + exp( (abs(x) - R)/a ) ) !Woods-Saxon Potential

 return
 end 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C   package EIGLIB   -- from numerical recipes
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine eig(A,n,np,e,vec,work)
! C
! C  INPUT: 
! C       A(i,j): symmetry n x n matrix (declared np x np)
! C      n = working dimension of A
! C     np = declared dimension of A
! C    work(i) = dummy array
! C
! C   OUTPUT:
! C     e(i) = vector of eigenvalues of A
! C     vec(i,j) = arrays of eigenvectors of A
! C
! C   SUBROUTINES CALLED:
! C     tred2: reduces A to tridiagonal
! C     tqli:   find eigenvalues, eigenvectors of tridiagonal
! C
      implicit none
      
      integer n,np
      real a(np,np)		! array to be diagonalized 
      real e(np)		! eigenvalues out 
      real vec(np,np)		! eigenvectors out
      real work(np)		!
      integer i,j
      
      do i =1,n
      	do j = 1,n
      	     vec(i,j)=a(i,j)
      	enddo
      enddo
      call tred2(vec,n,np,e,work)
      call tqli(e,work,n,np,vec)
      return
      end

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! C   subroutine TRED2
! C   
! C   reduces a matrix to tridiagonal; taken from numerical recipes
! C
! C  INPUT:
! C  A(i,j) : n x n symmetric matrix (this gets destroyed)
! C  n = working dimension of A
! C  np = declared dimension of A
! C  
! C  OUTPUT:
! C     d(i): diagonal elements
! C     e(i): off diagonal elements 
! C  NB: A is replaced by orthogonal transformation
! C

      SUBROUTINE tred2(a,n,np,d,e)
      INTEGER n,np
      REAL a(np,np),d(np),e(np)
      INTEGER i,j,k,l
      REAL f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.
        scale=0.
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(sqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.
            do 15 j=1,l
! C     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
! C     Omit following line if finding only eigenvalues.
      d(1)=0.
      e(1)=0.
      do 24 i=1,n
! C     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.)then
          do 22 j=1,l
            g=0.
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
!C     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
!C     Also delete lines from here ...
        a(i,i)=1.
        do 23 j=1,l
          a(i,j)=0.
          a(j,i)=0.
23      continue
!C     ... to here when finding only eigenvalues.
24    continue
      return
      END


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! C  subroutine TQLI
! C   
! C  solves tridiagonal matrix for eigenvalues through implicit QL shift
! C  (see numerical recipes)
! C
! C INPUT: 
! C   d(i): diagonal matrix elements
! C   e(i): off-diagonal matrix elements
! C   n = working dimension
! C   np = declared dimension
! C   z(i,j): output from TRED2; this will transform back to original basis
! C
! C  OUTPUT:
! C   d(i): will hold eigenvalues
! C   z(i,j) will hold eigenvectors in original basis  
! C 
! C FUNCTIONS CALLED: 
! C  pythag
! C

      SUBROUTINE tqli(d,e,n,np,z)
      INTEGER n,np
      REAL d(np),e(np),z(np,np)
!CU    USES pythag
      INTEGER i,iter,k,l,m
      REAL b,c,dd,f,g,p,r,s,pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.1000) print*,'too many iterations in tqli' !*** I changed from 30 to 300
          iter=iter+1
          g=(d(l+1)-d(l))/(2.*e(l))
          r=pythag(g,1.)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.
          c=1.
          p=0.
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.)then
              d(i+1)=d(i+1)-p
              e(m)=0.
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
!C     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
!C     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.
          goto 1
        endif
15    continue
      return
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION pythag(a,b)
      REAL a,b,pythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C
! C  subroutine eigsrt
! C  sorts eigenvalues into ASCENDING ORDER
! C
! C  INPUT:
! C  d(i): eigenvalues
! C   v(i,j): matrix of eigenvectors
! C   n = working dimension
! C  np = declared dimension   
! C
      SUBROUTINE eigsrt(d,v,n,np)
      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
      
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! c
! c	subroutine printnicematrix(size, n, array)
! c
! c	Prints a nice looking matrix
! c
! c	Input:
! c		integer size - maxsize of matrix
! c		integer n - actual size of matrix
! c		real array(size, size) - the matrix to print
! c

      subroutine printnicematrix(size,n,array)

      implicit none

      integer size,n
      real array(size,size)
      integer i,j

      do i =1,n
           write(6,100)(array(i,j),j=1,n)
       enddo
  100 format(10f8.3)
       return
       end



