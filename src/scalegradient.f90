 !**********************************************************
 ! Scale gradient
 !**********************************************************
 
 SUBROUTINE scalegradient(grad,cost,m,perc,alpha)
 USE globalp
 INTEGER:: i,imax
 REAL*8:: perc,cost,gradmax,alpha
 REAL:: m(nop),grad(nop)
 
 gradmax=-1e20
 imax=-1
 
 Do i=1, nop
    IF (abs(grad(i)).GT.gradmax) THEN
       gradmax=abs(grad(i))
       imax=i
    ENDif
 ENDDO
 write(*,*)imax,m(imax)
 WRITE(*,*) 'GRADMAX PERC = ',gradmax,perc
 alpha=m(imax)*perc/gradmax
 WRITE(*,*) 'alpha = ',alpha
 
 END SUBROUTINE scalegradient
 
