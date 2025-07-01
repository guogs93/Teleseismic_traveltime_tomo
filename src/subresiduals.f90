!*******************************************************
! calculate traveltime residuals at receiver locations
!*******************************************************

SUBROUTINE subresiduals(ishot,trcal,trobs,residuals)
 USE globalp
 INTEGER:: i,ishot
 REAL*8:: trcal(nrc),trobs(nrc),residuals(nrc)
 
 DO i=1, nrc
    IF(tstat(ishot,i).EQ.1)THEN
      residuals(i)=(trcal(i)-trobs(i))/cd(ishot,i)
    ELSE
      residuals(i)=0
    ENDIF
 ENDDO
END SUBROUTINE subresiduals
