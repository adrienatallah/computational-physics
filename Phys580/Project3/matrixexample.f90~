      program matrixexample
C
C   illustration inverting matrices
C   using LU decomposition
C
C  CWJ SDSU March 2005
C
C  uses numerical receipes program LUdcmp and LUbksb
C 

      implicit none

C.............first dimension arrays..................

      integer size        ! declared dimension of arrays
      parameter(size = 10)

      real a(size,size)   ! the matrix in question
      real atemp(size,size) ! needed as disposable array
      real ainv(size,size)    ! = A^-1
      
      integer matdim   ! actual dimension used

      integer indx(size)  ! needed for LU routines 

C..................dummy indices for copying matrices........

      integer i,j
C........... get the matrix from a file

      call getmatrix(size,matdim,a)
      
      print*,' The original matrix is: '
      call printnicematrix(size,matdim,a)
      print*,' '

C................ copy original file 

      do i =1,matdim
          do j = 1,matdim
              atemp(i,j) = a(i,j)
          enddo
      enddo
C.................. invert..................

      call matinv(atemp,ainv,matdim,size,indx)
      
C.................print out inverted.........

      print*,' inverted matrix = '
      call printnicematrix(size,matdim,ainv)
    
      print*,' ' 

C.................test................

      print*,' As a test, compute A^-1 A  = '
      call multmat(size,matdim,ainv,a,atemp)
    
      call printnicematrix(size,matdim,atemp)

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getmatrix(size,n,array)

      implicit none

C...................... matrix info......................

      integer size
      integer n
      real array(size,size)

C.....................file variables .......................

      logical flag
      integer ilast
      character filename*15
C...................dummy indices.................

      integer i,j

C...................get name of file................

       flag = .false.

      do while(.not.flag)
          write(6,*)' Enter name of .dat file where matrix is stored'
          write(6,*)'(Do not enter the .dat extension )'
          read(5,*)filename
          ilast = index(filename,' ')-1
          open(unit=1,file=filename(1:ilast)//'.dat',status='old',
     & err=10)
          flag=.true.
10      continue
          if(.not.flag)print*,'That file does not exist '
      enddo
      read(1,*)n
      print*,' Dimension of matrix is ',n
   
      do i = 1,n
         read(1,*,err=101)(array(i,j),j=1,n)
      enddo

      return

  101 continue
         write(6,*)' not enough matrix elements '
          stop    

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine multmat(size,n,a,b,c)
      implicit none
      
      integer size
      real a(size,size),b(size,size),c(size,size)
      integer n
      integer i,j,k
      real temp


      do i = 1,n
        do j = 1,n
           temp = 0.0
           do k = 1,n
              temp = temp+a(i,k)*b(k,j)
           enddo
           c(i,j)=temp
        enddo
      enddo
      return
      end
