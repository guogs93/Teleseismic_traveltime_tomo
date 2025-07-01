!********************************************
! add gaussian noise for the synthetic data
!********************************************
subroutine gaussnoise(tres,tresnoise)
use globalp 
REAL*8:: tres(nsrc,nrc), tresnoise(nsrc,nrc)
INTEGER:: is, ir
REAL*8, EXTERNAL:: gasdev

do is=1, nsrc
do ir=1, nrc
   tresnoise(is,ir)=tres(is,ir)+gasdev(rseed)*sdgn
enddo
enddo

end subroutine

REAL*8 FUNCTION gasdev(idum)
IMPLICIT NONE
INTEGER :: iset,i
INTEGER, PARAMETER :: imax=100000
INTEGER :: idum
REAL*8:: fac,rsq,v1,v2
REAL*8, SAVE :: gset
REAL*8, EXTERNAL :: ran1

iset=0
IF(iset.EQ.0)then
   DO i=1,imax
      v1=2.0*ran1(idum)-1.
      v2=2.0*ran1(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.LT.1.AND.rsq.NE.0.)EXIT
   ENDDO
   fac=sqrt(-2.0*LOG(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
ELSE
   gasdev=gset
   iset=0
ENDIF
END FUNCTION gasdev

REAL*8 FUNCTION ran1(idum)
IMPLICIT NONE
INTEGER :: idum
INTEGER, PARAMETER :: ia=16807,im=2147483647,iq=127773
INTEGER, PARAMETER :: ir=2836,ntab=32,ndiv=1+(im-1)/ntab
INTEGER :: j,k
INTEGER, SAVE :: iy
INTEGER, DIMENSION (:), SAVE ::iv(ntab)
REAL*8, PARAMETER :: eps=1.2e-7,rnmx=1.0-eps,am=1./im

iv=ntab*0
iy=0
IF(idum.LE.0.OR.iy.EQ.0)THEN
   DO j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      IF(idum.LT.0)idum=idum+im
      IF(j.LE.ntab)iv(j)=idum
   ENDDO
   iy=iv(1)
ENDIF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF(idum.LT.0)idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=min(am*iy,rnmx)

END FUNCTION ran1
