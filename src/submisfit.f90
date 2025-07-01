!******************************************************
! misfit function for L2 norm
!******************************************************

SUBROUTINE submisfit(residuals,vel,velr,cost,cost1,cost2)  
 USE globalp
 INTEGER:: i,j,k
 REAL:: cost1,cost2,cost
 REAL*8:: residuals(nsrc,nrc),velr(0:nvp+1,0:nvt+1,0:nvr+1),vel(0:nvp+1,0:nvt+1,0:nvr+1)

 cost1=0.
 DO j=1, nsrc
 DO i=1, nrc
    cost1=cost1+residuals(j,i)**2
 ENDDO
 ENDDO
 
 cost2=0
 DO i=1, nvp
 DO j=1, nvt
 DO k=1, nvr
    cost2=cost2+((vel(i,j,k)-velr(i,j,k)))**2/cm(i,j,k)
 ENDDO
 ENDDO
 ENDDO
 
 write(*,*)'data variance',cost1,'model variance',cost2
 
 cost=0.5*cost1+0.5*epsilon*cost2
 
END SUBROUTINE submisfit
