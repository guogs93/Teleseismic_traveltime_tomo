!********************************************
! transformation between 3D and 1D data
!********************************************

SUBROUTINE subm2v(m,v)
 USE globalp
 INTEGER:: i,j,k,r
 REAL:: m(nop)
 REAL*8:: v(0:nvp+1,0:nvt+1,0:nvr+1)

 r=0
 
 DO i=0, nvp+1
 DO j=0, nvt+1
 DO k=0, nvr+1
    r=r+1
    v(i,j,k)=m(r)
 ENDDO
 ENDDO
 ENDDO
 
END SUBROUTINE subm2v
 
 
SUBROUTINE subv2m(v,m)
 USE globalp
 INTEGER:: i,j,k,r
 REAL:: m(nop)
 REAL*8:: v(0:nvp+1,0:nvt+1,0:nvr+1)
 
 r=0
 DO i=0, nvp+1
 DO j=0, nvt+1
 DO k=0, nvr+1
    r=r+1
    m(r)=v(i,j,k)
 ENDDO
 ENDDO
 ENDDO

END SUBROUTINE subv2m
