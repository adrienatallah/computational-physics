  Program proj3
    implicit none
    INTEGER::i, j, k, l, m, n  
    
    REAL::datad(3, 15), year(15), temp(15), var(15) 
   ! temp = 0.0
    
    do i = 1, 3
    	do j = 1,15
    		datad(i,j) = 0.0   
		temp(j) = 0.0 		
		year(j) = 0.0
		var(j) = 0.0
    	enddo
    enddo    
    
    OPEN(unit=1, FILE='dataset4.dat', status='old', form='formatted')  

!print data from file
print*, ' '	
print*, 'the data from the file is:'
print*, '1st column: year'
print*, '2nd column: tempurature'
print*, '3rd column: variance'
print*, ' '

 do i = 1,3      
   	read(1,*) (datad(i,k), k=1,15)
   	!read(1,*) (year(k), k=1,15)
   	write(6,'(3f8.2)') (datad(i,k), k=1,15)
 enddo 
 print*, ' '
 
 do l = 1, 15
 	print *, datad(2, l)
 	datad(1, l) = year(l)
 	datad(2, l) = temp(l)
 	datad(3, l) = var(l)
 	!print *, year(l)
 enddo
 
 END  
