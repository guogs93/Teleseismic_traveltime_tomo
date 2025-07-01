!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module. 
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: checkstat
INTEGER, SAVE :: nnr,nnx,nnz,nrc,fom,ltfr
INTEGER, SAVE :: nvr,nvx,nvz,nrfr,nrfx,nrfz
INTEGER, SAVE :: nsnn,nsns,nsne,nsnw
REAL(KIND=i10), SAVE :: gor,gox,goz,dnr,dnx,dnz,snb,earth
REAL(KIND=i10), SAVE :: dvr,dvx,dvz
REAL(KIND=i10), SAVE :: goxd,gozd,dnxd,dnzd
REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: veln
REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: ttn 
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE, SAVE :: rcr,rcx,rcz
REAL(KIND=i10), PARAMETER :: pi=3.1415926535898
!
! nnr,nnx,nnz = Number of nodes of refined grid in r,x and z
! nvr,nvx,nvz = Number of B-spline vertices in r,x and z
! gor,gox,goz = Origin of grid (radius,theta,phi)
! dnr,dnx,dnz = Node separation of refined grid in r, x and z
! dvr,dvx,dvz = Node separation of B-spline grid in r, x and z
! nrfr,nrfx,nrfz = B-spline dicing level in r, x and z
! veln(i,j,k) = velocity values on a refined grid of nodes
! ttn(i,j,k) = traveltime field on the refined grid of nodes
! checkstat = check status of memory allocation
! fom = use first-order(0) or mixed-order(1) scheme
! snb = Maximum size of narrow band as fraction of nnx*nnz
! nrc = number of receivers
! rcr(i),rcx(i),rcz(i) = (r-earth,x,z) coordinates of receivers
! earth = radius of Earth (in km)
! goxd,gozd = gox,goz in radians
! dnxd,dnzd = dnx,dnz in radians
! dfr,dfx,dfz = B-spline dicing factor in r, theta and phi
! ltfr = Limit traveltime field to receiver array (0=no, 1=yes)
! nsnn,nsns = Latitude bounds (N-S) of surface receiver grid
! nsne,nsnw = Longitude bounds (E-W) of surface receiver grid
!
END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module contains all the subroutines used to calculate
! the first-arrival traveltime field through the grid.
! Subroutines are:
! (1) travel
! (2) setup
! (3) sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE traveltime
USE globalp
IMPLICIT NONE
INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: nsts
INTEGER imin(8),imax(8),istep(8),jmin(8),jmax(8),jstep(8),kmin(8),kmax(8),kstep(8)
REAL(KIND=i10), DIMENSION (:,:,:), ALLOCATABLE :: ttno 
! nsts(i,j,k) = node status (1=source,0= determined nodes) 


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using fast sweeping method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE travel(srcid)
IMPLICIT NONE
INTEGER :: i,j,k,sw,srcid
INTEGER :: isum,tnrn,nsw,ite
REAL(KIND=i10) :: rd1
REAL(KIND=i10) :: txmin,tzmin,trmin
REAL(KIND=i10) :: a,b,c,tbar,slown,ri,risti
REAL(KIND=i10) :: ww,vv,uu,hw,hv,hu
REAL(KIND=i10) :: cvg_err

! rd1 = substitution variable
! srcid = id number of source
! sw = switch (0 or 1)
! tnrn = total number of receiver nodes
!
! Allocate nsts and set all elements of array nsts equal to 
! (-1) since to begin with, all points are "big" points.
!
ALLOCATE(nsts(0:nnz+1,0:nnx+1,0:nnr+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: INTEGER nsts'
ENDIF
nsts=-1

ALLOCATE(ttno(nnz,nnx,nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE travel: real ttno'
ENDIF

!
! Read in traveltimes calculated a priori to base of grid,
! find point with minimum traveltime,
!
sw=0
ttn=1.0d10
DO i=1, nnz
   DO j=1, nnx
      READ(20)rd1
      IF(rd1.GT.0.0)THEN
         ttn(i,j,nnr)=rd1
         nsts(i,j,nnr)=0
         sw=1
      ENDIF
   ENDDO
ENDDO

IF(sw.EQ.0)THEN
   WRITE(6,*)'No initial traveltimes for source ',srcid
   WRITE(6,*)'Terminating program'
   STOP
ENDIF

CALL setup

ite=0
cvg_err=1.0d10
!write(6,*)'srcid=',srcid
!write(6,*)'cvg_err=',cvg_err
!pause

!DO WHILE((cvg_err.gt.1.0d-1).and.(ite.le.2000))
   !ttno=ttn
   IF(fom.EQ.0)THEN
     CALL sweep1
   ELSEIF(fom.EQ.1)THEN
     CALL sweep2
   ELSE
     CALL sweep3
   ENDIF
   !ite=ite+1
   !write(6,*)'ite=',ite
   !pause
   !cvg_err=0
   !cvg_err=sum(abs(ttno-ttn))/dble(nnr*nnx*nnz) 
   !write(6,*)'cvg_err=',cvg_err
   !pause
!ENDDO

DEALLOCATE(nsts, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE travel: nsts'
ENDIF

DEALLOCATE(ttno, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE travel: ttno'
ENDIF

END SUBROUTINE travel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is dertermine the direction of sweeping.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE setup
IMPLICIT NONE
imin(1)=1;imax(1)=nnr;istep(1)=1
jmin(1)=1;jmax(1)=nnx;jstep(1)=1
kmin(1)=1;kmax(1)=nnz;kstep(1)=1

imin(2)=1;imax(2)=nnr;istep(2)=1
jmin(2)=nnx;jmax(2)=1;jstep(2)=-1
kmin(2)=1;kmax(2)=nnz;kstep(2)=1

imin(3)=1;imax(3)=nnr;istep(3)=1
jmin(3)=1;jmax(3)=nnx;jstep(3)=1
kmin(3)=nnz;kmax(3)=1;kstep(3)=-1

imin(4)=1;imax(4)=nnr;istep(4)=1
jmin(4)=nnx;jmax(4)=1;jstep(4)=-1
kmin(4)=nnz;kmax(4)=1;kstep(4)=-1 

imin(5)=nnr;imax(5)=1;istep(5)=-1
jmin(5)=1;jmax(5)=nnx;jstep(5)=1
kmin(5)=1;kmax(5)=nnz;kstep(5)=1

imin(6)=nnr;imax(6)=1;istep(6)=-1
jmin(6)=1;jmax(6)=nnx;jstep(6)=1
kmin(6)=nnz;kmax(6)=1;kstep(6)=-1

imin(7)=nnr;imax(7)=1;istep(7)=-1
jmin(7)=nnx;jmax(7)=1;jstep(7)=-1
kmin(7)=nnz;kmax(7)=1;kstep(7)=-1

imin(8)=nnr;imax(8)=1;istep(8)=-1
jmin(8)=nnx;jmax(8)=1;jstep(8)=-1
kmin(8)=1;kmax(8)=nnz;kstep(8)=1
 

END SUBROUTINE setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 1-order fast sweeping method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sweep1
IMPLICIT NONE
INTEGER :: i,j,k,swp
REAL(KIND=i10) :: txmin,tzmin,trmin
REAL(KIND=i10) :: a,b,c,tbar,slown,ri,risti
REAL(KIND=i10) :: ww,vv,uu,hw,hv,hu

DO swp=1, 8
   DO i=imin(swp),imax(swp),istep(swp)
   DO j=jmin(swp),jmax(swp),jstep(swp)
   DO k=kmin(swp),kmax(swp),kstep(swp)  
      IF(nsts(k,j,i).NE.0)THEN       
        slown=1.0/veln(k,j,i)
        ri=gor-(i-1)*dnr+earth
        risti=ri*sin(gox+(j-1)*dnx)

        IF(i.GT.1.AND.i.LT.nnr)THEN
           trmin=MIN(ttn(k,j,i+1),ttn(k,j,i-1))
        ELSEIF(i.EQ.1)THEN
           trmin=ttn(k,j,i+1)           
        ELSE
           trmin=ttn(k,j,i-1)
        ENDIF      

        IF(j.GT.1.AND.j.LT.nnx)THEN
           txmin=MIN(ttn(k,j+1,i),ttn(k,j-1,i))      
        ELSEIF(j.EQ.1)THEN
           txmin=ttn(k,j+1,i)
        ELSE
           txmin=ttn(k,j-1,i)
        ENDIF

        IF(k.GT.1.AND.k.LT.nnz)THEN
           tzmin=MIN(ttn(k+1,j,i),ttn(k-1,j,i))      
        ELSEIF(k.EQ.1)THEN
           tzmin=ttn(k+1,j,i)
        ELSE
           tzmin=ttn(k-1,j,i)
        ENDIF

        ww=trmin
        vv=txmin
        uu=tzmin
        hw=dnr
        hv=ri*dnx
        hu=risti*dnz

        IF(ww.GT.vv)THEN
           call swap(ww,vv)
           call swap(hw,hv)
        ENDIF

        IF(ww.GT.uu)THEN
           call swap(ww,uu)
           call swap(hw,hu)
        ENDIF

        IF(vv.GT.uu)THEN
           call swap(vv,uu)
           call swap(hv,hu)
        ENDIF

        IF(.NOT.ww.LE.vv.AND.vv.LE.uu)THEN
           WRITE(6,*)"ERROR with sorting traveltimes"
        ENDIF

        tbar=ww+slown*hw
        IF(tbar.GT.vv)THEN
           a=hv**2+hw**2
           b=-2*(ww*hv**2+vv*hw**2)
           c=(hv*ww)**2+(hw*vv)**2-(slown*hv*hw)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF
        IF(tbar.GT.uu)THEN
           a=(hv*hu)**2+(hw*hu)**2+(hw*hv)**2
           b=-2*(ww*(hv*hu)**2+vv*(hw*hu)**2+uu*(hw*hv)**2)
           c=(ww*hv*hu)**2+(vv*hw*hu)**2+(uu*hw*hv)**2-(slown*hw*hv*hu)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF    

        ttn(k,j,i)=min(tbar,ttn(k,j,i))  
               
      ENDIF
  ENDDO
  ENDDO
  ENDDO
ENDDO
END SUBROUTINE sweep1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 2-order fast sweeping method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sweep2
IMPLICIT NONE
INTEGER :: i,j,k,swp
REAL(KIND=i10) :: txmin,tzmin,trmin,tr1,tr2,tx1,tx2,tz1,tz2
REAL(KIND=i10) :: a,b,c,tbar,slown,ri,risti
REAL(KIND=i10) :: ww,vv,uu,hw,hv,hu,hr,hx,hz

DO swp=1, 8
   DO i=imin(swp),imax(swp),istep(swp)
   DO j=jmin(swp),jmax(swp),jstep(swp)
   DO k=kmin(swp),kmax(swp),kstep(swp)  
      IF(nsts(k,j,i).NE.0)THEN       
        slown=1.0/veln(k,j,i)
        ri=gor-(i-1)*dnr+earth
        risti=ri*sin(gox+(j-1)*dnx)

        IF(i.GT.1.AND.i.LT.nnr)THEN
           trmin=MIN(ttn(k,j,i-1),ttn(k,j,i+1))
           hr=dnr
           IF(trmin.EQ.ttn(k,j,i-1))THEN
              IF(i-2.GE.1)THEN
                 IF(ttn(k,j,i-2).LT.ttn(k,j,i-1))THEN
                    trmin=ttn(k,j,i-1)+(ttn(k,j,i-1)-ttn(k,j,i-2))/3
                    hr=dnr*2/3
                 ENDIF
              ENDIF
           ELSE
              IF(i+2.LE.nnr)THEN
                 IF(ttn(k,j,i+2).LT.ttn(k,j,i+1))THEN
                    trmin=ttn(k,j,i+1)+(ttn(k,j,i+1)-ttn(k,j,i+2))/3
                    hr=dnr*2/3
                 ENDIF
              ENDIF
           ENDIF
        ELSEIF(i.EQ.1)THEN
           trmin=ttn(k,j,i+1)
           hr=dnr           
        ELSEIF(i.EQ.nnr)THEN
           trmin=ttn(k,j,i-1)
           hr=dnr
        ENDIF      

        IF(j.GT.1.AND.j.LT.nnx)THEN
           txmin=MIN(ttn(k,j-1,i),ttn(k,j+1,i))
           hx=dnx
           IF(txmin.EQ.ttn(k,j-1,i))THEN
              IF(j-2.GE.1)THEN
                 IF(ttn(k,j-2,i).LT.ttn(k,j-1,i))THEN
                    txmin=ttn(k,j-1,i)+(ttn(k,j-1,i)-ttn(k,j-2,i))/3
                    hx=dnx*2/3
                 ENDIF
              ENDIF
           ELSE
              IF(j+2.LE.nnx)THEN
                 IF(ttn(k,j+2,i).LT.ttn(k,j+1,i))then
                    txmin=ttn(k,j+1,i)+(ttn(k,j+1,i)-ttn(k,j+2,i))/3
                    hx=dnx*2/3
                 ENDIF
              ENDIF
           ENDIF    
        ELSEIF(j.EQ.1)THEN
           txmin=ttn(k,j+1,i)
           hx=dnx
        ELSEIF(j.EQ.nnx)THEN
           txmin=ttn(k,j-1,i)
           hx=dnx
        ENDIF

        IF(k.GT.1.AND.k.LT.nnz)THEN
           tzmin=MIN(ttn(k-1,j,i),ttn(k+1,j,i))
           hz=dnz
           IF(tzmin.EQ.ttn(k-1,j,i))THEN
              IF(k-2.GE.1)THEN
                IF(ttn(k-2,j,i).LT.ttn(k-1,j,i))THEN
                   tzmin=ttn(k-1,j,i)+(ttn(k-1,j,i)-ttn(k-2,j,i))/3
                   hz=dnz*2/3
                ENDIF
              ENDIF
           ELSE
              IF(k+2.LE.nnz)THEN
                IF(ttn(k+2,j,i).LT.ttn(k+1,j,i))THEN
                   tzmin=ttn(k+1,j,i)+(ttn(k+1,j,i)-ttn(k+2,j,i))/3
                   hz=dnz*2/3
                ENDIF
              ENDIF
           ENDIF 
        ELSEIF(k.EQ.1)THEN
           tzmin=ttn(k+1,j,i)
           hz=dnz
        ELSEIF(k.EQ.nnz)THEN
           tzmin=ttn(k-1,j,i)
           hz=dnz
        ENDIF

        ww=trmin
        vv=txmin
        uu=tzmin
        hw=hr
        hv=ri*hx
        hu=risti*hz

        IF(ww.GT.vv)THEN
           call swap(ww,vv)
           call swap(hw,hv)
        ENDIF

        IF(ww.GT.uu)THEN
           call swap(ww,uu)
           call swap(hw,hu)
        ENDIF

        IF(vv.GT.uu)THEN
           call swap(vv,uu)
           call swap(hv,hu)
        ENDIF

        IF(.NOT.ww.LE.vv.AND.vv.LE.uu)THEN
           WRITE(6,*)"ERROR with sorting traveltimes"
        ENDIF

        tbar=ww+slown*hw
        IF(tbar.GT.vv)THEN
           a=hv**2+hw**2
           b=-2*(ww*hv**2+vv*hw**2)
           c=(hv*ww)**2+(hw*vv)**2-(slown*hv*hw)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF
        IF(tbar.GT.uu)THEN
           a=(hv*hu)**2+(hw*hu)**2+(hw*hv)**2
           b=-2*(ww*(hv*hu)**2+vv*(hw*hu)**2+uu*(hw*hv)**2)
           c=(ww*hv*hu)**2+(vv*hw*hu)**2+(uu*hw*hv)**2-(slown*hw*hv*hu)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF    

        ttn(k,j,i)=min(tbar,ttn(k,j,i))  
               
      ENDIF
  ENDDO
  ENDDO
  ENDDO
ENDDO
END SUBROUTINE sweep2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 3-order fast sweeping method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sweep3
IMPLICIT NONE
INTEGER :: i,j,k,swp
REAL(KIND=i10) :: txmin,tzmin,trmin,tr1,tr2,tx1,tx2,tz1,tz2
REAL(KIND=i10) :: a,b,c,tbar,slown,ri,risti
REAL(KIND=i10) :: ww,vv,uu,hw,hv,hu,hr,hx,hz

DO swp=1, 8
   DO i=imin(swp),imax(swp),istep(swp)
   DO j=jmin(swp),jmax(swp),jstep(swp)
   DO k=kmin(swp),kmax(swp),kstep(swp)  
      IF(nsts(k,j,i).NE.0)THEN       
        slown=1.0/veln(k,j,i)
        ri=gor-(i-1)*dnr+earth
        risti=ri*sin(gox+(j-1)*dnx)

        IF(i.GT.1.AND.i.LT.nnr)THEN
           trmin=MIN(ttn(k,j,i-1),ttn(k,j,i+1))
           hr=dnr
           IF(trmin.EQ.ttn(k,j,i-1))THEN
              IF(i-2.GE.1)THEN
                 IF(ttn(k,j,i-2).LT.ttn(k,j,i-1))THEN
                    trmin=ttn(k,j,i-1)+(ttn(k,j,i-1)-ttn(k,j,i-2))/3
                    hr=dnr*2/3
                    IF(i-3.GE.1)THEN
                       IF(ttn(k,j,i-3).LT.ttn(k,j,i-2))THEN
                          trmin=ttn(k,j,i-1)+(7*ttn(k,j,i-1)-9*ttn(k,j,i-2)+2*ttn(k,j,i-3))/11
                          hr=dnr*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ELSE
              IF(i+2.LE.nnr)THEN
                 IF(ttn(k,j,i+2).LT.ttn(k,j,i+1))THEN
                    trmin=ttn(k,j,i+1)+(ttn(k,j,i+1)-ttn(k,j,i+2))/3
                    hr=dnr*2/3
                    IF(i+3.LE.nnr)THEN
                       IF(ttn(k,j,i+3).LT.ttn(k,j,i+2))THEN
                          trmin=ttn(k,j,i+1)+(7*ttn(k,j,i+1)-9*ttn(k,j,i+2)+2*ttn(k,j,i+3))/11
                          hr=dnr*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ELSEIF(i.EQ.1)THEN
           trmin=ttn(k,j,i+1)
           hr=dnr           
        ELSEIF(i.EQ.nnr)THEN
           trmin=ttn(k,j,i-1)
           hr=dnr
        ENDIF      

        IF(j.GT.1.AND.j.LT.nnx)THEN
           txmin=MIN(ttn(k,j-1,i),ttn(k,j+1,i))
           hx=dnx
           IF(txmin.EQ.ttn(k,j-1,i))THEN
              IF(j-2.GE.1)THEN
                 IF(ttn(k,j-2,i).LT.ttn(k,j-1,i))THEN
                    txmin=ttn(k,j-1,i)+(ttn(k,j-1,i)-ttn(k,j-2,i))/3
                    hx=dnx*2/3
                    IF(j-3.GE.1)THEN
                       IF(ttn(k,j-3,i).LT.ttn(k,j-2,i))THEN
                          txmin=ttn(k,j-1,i)+(7*ttn(k,j-1,i)-9*ttn(k,j-2,i)+2*ttn(k,j-3,i))/11
                          hx=dnx*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ELSE
              IF(j+2.LE.nnx)THEN
                 IF(ttn(k,j+2,i).LT.ttn(k,j+1,i))then
                    txmin=ttn(k,j+1,i)+(ttn(k,j+1,i)-ttn(k,j+2,i))/3
                    hx=dnx*2/3
                    IF(j+3.LE.nnx)THEN
                       IF(ttn(k,j+3,i).LT.ttn(k,j+2,i))THEN
                          txmin=ttn(k,j+1,i)+(7*ttn(k,j+1,i)-9*ttn(k,j+2,i)+2*ttn(k,j+3,i))/11
                          hx=dnx*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF    
        ELSEIF(j.EQ.1)THEN
           txmin=ttn(k,j+1,i)
           hx=dnx
        ELSEIF(j.EQ.nnx)THEN
           txmin=ttn(k,j-1,i)
           hx=dnx
        ENDIF

        IF(k.GT.1.AND.k.LT.nnz)THEN
           tzmin=MIN(ttn(k-1,j,i),ttn(k+1,j,i))
           hz=dnz
           IF(tzmin.EQ.ttn(k-1,j,i))THEN
              IF(k-2.GE.1)THEN
                IF(ttn(k-2,j,i).LT.ttn(k-1,j,i))THEN
                   tzmin=ttn(k-1,j,i)+(ttn(k-1,j,i)-ttn(k-2,j,i))/3
                   hz=dnz*2/3
                   IF(k-3.GE.1)THEN
                       IF(ttn(k-3,j,i).LT.ttn(k-2,j,i))THEN
                          tzmin=ttn(k-1,j,i)+(7*ttn(k-1,j,i)-9*ttn(k-2,j,i)+2*ttn(k-3,j,i))/11
                          hz=dnz*6/11
                       ENDIF               
                    ENDIF
                ENDIF
              ENDIF
           ELSE
              IF(k+2.LE.nnz)THEN
                IF(ttn(k+2,j,i).LT.ttn(k+1,j,i))THEN
                   tzmin=ttn(k+1,j,i)+(ttn(k+1,j,i)-ttn(k+2,j,i))/3
                   hz=dnz*2/3
                   IF(k+3.LE.nnz)THEN
                       IF(ttn(k+3,j,i).LT.ttn(k+2,j,i))THEN
                          tzmin=ttn(k+1,j,i)+(7*ttn(k+1,j,i)-9*ttn(k+2,j,i)+2*ttn(k+3,j,i))/11
                          hz=dnz*6/11
                       ENDIF               
                    ENDIF
                ENDIF
              ENDIF
           ENDIF 
        ELSEIF(k.EQ.1)THEN
           tzmin=ttn(k+1,j,i)
           hz=dnz
        ELSEIF(k.EQ.nnz)THEN
           tzmin=ttn(k-1,j,i)
           hz=dnz
        ENDIF

        ww=trmin
        vv=txmin
        uu=tzmin
        hw=hr
        hv=ri*hx
        hu=risti*hz

        IF(ww.GT.vv)THEN
           call swap(ww,vv)
           call swap(hw,hv)
        ENDIF

        IF(ww.GT.uu)THEN
           call swap(ww,uu)
           call swap(hw,hu)
        ENDIF

        IF(vv.GT.uu)THEN
           call swap(vv,uu)
           call swap(hv,hu)
        ENDIF

        IF(.NOT.ww.LE.vv.AND.vv.LE.uu)THEN
           WRITE(6,*)"ERROR with sorting traveltimes"
        ENDIF

        tbar=ww+slown*hw
        IF(tbar.GT.vv)THEN
           a=hv**2+hw**2
           b=-2*(ww*hv**2+vv*hw**2)
           c=(hv*ww)**2+(hw*vv)**2-(slown*hv*hw)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF
        IF(tbar.GT.uu)THEN
           a=(hv*hu)**2+(hw*hu)**2+(hw*hv)**2
           b=-2*(ww*(hv*hu)**2+vv*(hw*hu)**2+uu*(hw*hv)**2)
           c=(ww*hv*hu)**2+(vv*hw*hu)**2+(uu*hw*hv)**2-(slown*hw*hv*hu)**2
           IF(b**2-4*a*c.ge.0)THEN
              tbar=(-b+sqrt(b**2-4*a*c))/(2*a)
           ENDIF
        ENDIF    

        ttn(k,j,i)=min(tbar,ttn(k,j,i))  
               
      ENDIF
  ENDDO
  ENDDO
  ENDDO
ENDDO
END SUBROUTINE sweep3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is used to swap the values of two traveltimes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE swap(a,b)
REAL(KIND=i10) :: a,b,tem
tem=a
a=b
b=tem
END SUBROUTINE

END MODULE traveltime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM   
! CODE: FORTRAN 90
! This program is designed to implement the Fast Marching
! Method (FMM) for calculating first-arrival traveltimes
! through a 3-D continuous velocity medium. It is written in
! Fortran 90, although it is probably more accurately
! described as Fortran 77 with some of the Fortran 90
! extensions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM fmmin3d
USE globalp
USE traveltime
IMPLICIT NONE
CHARACTER (LEN=20) :: sources,itimes,receivers,gridv
CHARACTER (LEN=20) :: travelt,rtravel,wrays,cdum
CHARACTER (LEN=20) :: frechet,pht
INTEGER :: ii,i,j,k,l,nsrc,wttf,fsrt,wrgf,cfd,tnr,nra,idum
INTEGER :: awttf,awrgf
REAL(KIND=i10) :: cslat,cslong,abw,abe,abn,abs
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: scr,scx,scz
!
! sources = File containing source locations
! itimes = File containing initial traveltimes
! receivers = File containing receiver locations
! gridv = File containing grid of velocity vertices
! travelt = File name for storage of traveltime field
! wttf = Write traveltimes to file? (0=no,>0=source id)
! awttf = Receiver array ID for traveltimes
! fom = Use first-order(0) or mixed-order(1) scheme
! nsrc = number of sources
! scr,scx,scz = source location in r,x,z
! fsrt = find source-receiver traveltimes? (0=no,1=yes)
! rtravel = output file for source-receiver traveltimes
! cdum = dummy character variable
! wrgf = write ray geometries to file? (<0=all,0=no,>0=source id)
! awrgf = Receiver array ID for ray paths
! wrays = file containing raypath geometries
! cfd = calculate Frechet derivatives? (0=no,1=yes)
! frechet = output file containing matrix of frechet derivatives
! pht = phase type (e.g. P, S, pP etc.)
! tnr = total number of rays
! nra = number of receiver arrays
! cslat,cslong = Latidude and longitude width of receiver cushion
! abw,abe,abn,abs = Limits of receiver grid (W, E, N, S)
!
OPEN(UNIT=10,FILE='fm3dt.in',STATUS='old')
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)sources
READ(10,1)itimes
READ(10,1)receivers
READ(10,1)gridv
READ(10,*)nrfr,nrfx,nrfz
READ(10,*)earth
READ(10,*)fom
READ(10,*)snb
READ(10,*)ltfr
READ(10,*)cslat,cslong
READ(10,1)cdum
READ(10,1)cdum
READ(10,1)cdum
READ(10,*)fsrt
READ(10,1)rtravel
READ(10,*)cfd
READ(10,1)frechet
READ(10,*)wttf,awttf
READ(10,1)travelt
READ(10,*)wrgf,awrgf
READ(10,1)wrays
1   FORMAT(a20)
CLOSE(10)
!
! Call a subroutine which reads in the velocity grid
!
CALL gridder(gridv)
!
! Compute the total number of ray paths if
! necessary
!
IF(wrgf.NE.0)THEN
   OPEN(UNIT=10,FILE=sources,STATUS='old')
   OPEN(UNIT=20,FILE=receivers,STATUS='old')
   READ(10,*)nra
   READ(20,*)idum
   IF(idum.NE.nra)THEN
      WRITE(6,*)'ERROR!!!'
      WRITE(6,*)'Source and receiver files are'
      WRITE(6,*)'inconsistent!!!!'
      WRITE(6,*)'First line of each file should'
      WRITE(6,*)'be identical!!!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   tnr=0
   DO i=1,nra
      READ(10,*)nsrc
      READ(20,*)nrc
      IF(wrgf.GT.0)THEN
         IF(awrgf.EQ.i)tnr=nrc
      ELSE
         tnr=tnr+nsrc*nrc
      ENDIF
      DO j=1,nsrc
         READ(10,*)
         READ(10,*)
      ENDDO
      DO j=1,nrc
         READ(20,*)
      ENDDO
   ENDDO
   CLOSE(10)
   CLOSE(20)
ENDIF
!
! Read in the number of source types (i.e. the number
! of receiver arrays).
!
Open(UNIT=60,FILE=sources,STATUS='old')
READ(60,*)nra
!
! Open receiver file if required.
!
IF(fsrt.eq.1)THEN
   OPEN(UNIT=70,FILE=receivers,status='old')
   READ(70,*)idum
   IF(idum.NE.nra)THEN
      WRITE(6,*)'ERROR!!!'
      WRITE(6,*)'Source and receiver files are'
      WRITE(6,*)'inconsistent!!!!'
      WRITE(6,*)'First line of each file should'
      WRITE(6,*)'be identical!!!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
ENDIF
!
! Now work out, source by source, the first-arrival traveltime
! field plus traveltimes and ray paths from the base of the 3-D grid 
! if required. First, allocate memory to the
! traveltime field array
!
ALLOCATE(ttn(nnz,nnx,nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL ttn'
ENDIF
!
! Open file for source-receiver traveltime output if required.
!
IF(fsrt.eq.1)THEN
   OPEN(UNIT=10,FILE=rtravel,STATUS='unknown')
ENDIF
!
! Open file for ray path output if required
!
IF(wrgf.NE.0)THEN
   OPEN(UNIT=40,FILE=wrays,FORM='unformatted',STATUS='unknown')
   WRITE(40)tnr
ENDIF
!
! Open file for Frechet derivative output if required.
!
IF(cfd.EQ.1)THEN
   OPEN(UNIT=50,FILE=frechet,FORM='unformatted',STATUS='unknown')
ENDIF
!
! Open file containing traveltimes at base of model
!
OPEN(UNIT=20,FILE=itimes,FORM='unformatted',STATUS='old')
!
! Loop over each receiver array
!
DO ii=1,nra
   READ(60,*)nsrc
   ALLOCATE(scr(nsrc),scx(nsrc),scz(nsrc), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL scr,scx,scz'
   ENDIF
   DO i=1,nsrc
      READ(60,*)scr(i),scx(i),scz(i)
      READ(60,'(a8)')pht
!
!     Convert source coordinates in degrees to radians
!
      scx(i)=(90.0-scx(i))*pi/180.0
      scz(i)=scz(i)*pi/180.0
   ENDDO
   IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
      READ(70,*)nrc
      ALLOCATE(rcr(nrc),rcx(nrc),rcz(nrc), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM fmmin3d: REAL rcr,rcx,rcz'
      ENDIF
      DO i=1,nrc
         READ(70,*)rcr(i),rcx(i),rcz(i)
         IF(ltfr.EQ.1)THEN
            IF(i.EQ.1)THEN
               abe=rcz(i)
               abw=abe
               abn=rcx(i)
               abs=abn
            ELSE
               IF(rcz(i).LT.abw)THEN
                  abw=rcz(i)
               ELSE IF(rcz(i).GT.abe)THEN
                  abe=rcz(i)
               ENDIF
               IF(rcx(i).LT.abs)THEN
                  abs=rcx(i)
               ELSE IF(rcx(i).GT.abn)THEN
                  abn=rcx(i)
               ENDIF
            ENDIF
         ENDIF
!
!        Convert receiver coordinates in degrees to radians
!
         rcx(i)=(90.0-rcx(i))*pi/180.0
         rcz(i)=rcz(i)*pi/180.0
      ENDDO
      IF(ltfr.EQ.1)THEN
!
!        Add cushion values to bounds of receiver
!        array, and then locate surface grid which
!        is used to terminate the fast marching process.
!
         abw=(abw-cslong)*pi/180.0
         abe=(abe+cslong)*pi/180.0
         abs=abs-cslat
         abn=abn+cslat
         abs=(90.0-abs)*pi/180.0
         abn=(90.0-abn)*pi/180.0
         nsnn=INT((abn-gox)/dnx)
         IF(nsnn.LT.1)nsnn=1
         nsns=INT((abs-gox)/dnx)+2
         IF(nsns.GT.nnx)nsns=nnx
         nsne=INT((abe-goz)/dnz)+2
         IF(nsne.GT.nnz)nsne=nnz
         nsnw=INT((abw-goz)/dnz)
         IF(nsnw.LT.1)nsnw=1
      ENDIF
   ENDIF
   DO i=1,nsrc
!     
!     Call a subroutine that works out the first-arrival traveltime
!     field.
!
      CALL travel(i)
!
!     Find traveltimes from base of grid if required
!
      IF(fsrt.eq.1)THEN
         CALL srtimes
      ENDIF
!
!     Calculate raypath geometries and write to file if required.
!     Calculate Frechet derivatives with the same subroutine 
!     if required.
!
      IF((wrgf.EQ.i.AND.awrgf.EQ.ii).OR.wrgf.LT.0.OR.cfd.EQ.1)THEN
         !write(*,*)wrgf,i,awrgf,ii,cfd
         !pause
         CALL rpaths(wrgf,i,awrgf,ii,cfd)
      ENDIF
!
!     If required, write traveltime field to file
!
      IF(wttf.EQ.i.AND.awttf.EQ.ii)THEN
         OPEN(UNIT=30,FILE=travelt,FORM='unformatted',STATUS='unknown')
         WRITE(30)gor,goxd,gozd
         WRITE(30)nnr,nnx,nnz
         WRITE(30)dnr,dnxd,dnzd
         DO j=1,nnz
            DO k=1,nnx
               DO l=1,nnr
                  WRITE(30)ttn(j,k,l)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(30)
      ENDIF
   ENDDO
   DEALLOCATE(scr,scx,scz)
   IF(fsrt.EQ.1.or.ltfr.EQ.1)THEN
      DEALLOCATE(rcr,rcx,rcz)
   ENDIF
ENDDO
CLOSE(60)
!
! Close rtravel if required
!
IF(fsrt.eq.1)THEN
   CLOSE(70)
ENDIF
CLOSE(20)
IF(wrgf.NE.0)THEN
   CLOSE(40)
ENDIF
IF(cfd.EQ.1)THEN
   CLOSE(50)
ENDIF
DEALLOCATE (veln,ttn, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM fmmin3d: final deallocate'
ENDIF
WRITE(6,*)'Program fm3dt has finished successfully!'

STOP
END PROGRAM fmmin3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine is passed the name of the velocity
! grid file (grid) and reads in the velocity values.
! The gridded values are globally shared via
! a MODULE statement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gridder(gridv)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1
INTEGER :: conr,conx,conz,str,stx,stz
REAL(KIND=i10) :: u,sumi,sumj,sumk
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: velv
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi,wi
CHARACTER (LEN=20) :: gridv
!
! velv = B-spline velocity grid values
! u = Cubic spline independent variable
! ui,vi,wi = Cubic spline basis functions
! sumi,sumj,sumk = Summation variables for constructing spline
! conr,conx,conz = Counters for refining grid in r,theta,phi
! str,stx,stz = Refined grid location in r,theta,phi
! gridv = Input B-spline vertex file
!
! Read in B-spline grid
!
OPEN(UNIT=10,FILE=gridv,STATUS='old')
READ(10,*)nvr,nvx,nvz
READ(10,*)gor,goxd,gozd
READ(10,*)dvr,dvx,dvz
ALLOCATE(velv(0:nvr+1,0:nvx+1,0:nvz+1))
DO i=0,nvz+1
   DO j=0,nvx+1
      DO k=0,nvr+1
         READ(10,*)velv(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(10)
!
! Calculate total numer of refined nodes in r,theta,phi
! and the refined grid spacing.
!
nnr=(nvr-1)*nrfr+1
nnx=(nvx-1)*nrfx+1
nnz=(nvz-1)*nrfz+1
dnr=dvr/nrfr
dnxd=dvx/nrfx
dnzd=dvz/nrfz
ALLOCATE(veln(nnz,nnx,nnr))
!
! Calculate the values of the basis functions
!
ALLOCATE(ui(nvr+1,4))
DO i=1,nrfr+1
   u=nrfr
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(nvx+1,4))
DO i=1,nrfx+1
   u=nrfx
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
ALLOCATE(wi(nvz+1,4))
DO i=1,nrfz+1
   u=nrfz
   u=(i-1)/u
   wi(i,1)=(1.0-u)**3/6.0
   wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   wi(i,4)=u**3/6.0
ENDDO
!
! Calculate velocity values on refined grid
!
DO i=1,nvz-1
   conz=nrfz
   IF(i==nvz-1)conz=nrfz+1
   DO j=1,nvx-1
      conx=nrfx
      IF(j==nvx-1)conx=nrfx+1
      DO k=1,nvr-1
         conr=nrfr
         IF(k==nvr-1)conr=nrfr+1
         DO l=1,conz
            stz=nrfz*(i-1)+l
            DO m=1,conx
               stx=nrfx*(j-1)+m
               DO n=1,conr
                  str=nrfr*(k-1)+n
                  sumi=0.0
                  DO i1=1,4
                     sumj=0.0
                     DO j1=1,4
                        sumk=0.0
                        DO k1=1,4
                           sumk=sumk+ui(n,k1)*velv(k-2+k1,j-2+j1,i-2+i1)
                        ENDDO
                        sumj=sumj+vi(m,j1)*sumk
                     ENDDO
                     sumi=sumi+wi(l,i1)*sumj
                  ENDDO
                  veln(stz,stx,str)=sumi
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Finally, convert all parameters in degrees to radians
!
dnx=dnxd*pi/180.0
dnz=dnzd*pi/180.0
gox=(90.0-goxd)*pi/180.0
goz=gozd*pi/180.0
dvx=dvx*pi/180.0
dvz=dvz*pi/180.0
DEALLOCATE(velv)
DEALLOCATE(ui,vi,wi)
END SUBROUTINE gridder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates all receiver traveltimes for
! a given source and writes the results to file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE srtimes
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,irr,irx,irz,sw
REAL :: trr
REAL(KIND=i10) :: drr,drx,drz,produ
!
! irr,irx,irz = Coordinates of cell containing receiver
! trr = traveltime value at receiver
! produ = dummy multiplier
! drr,drx,drz = receiver distance from (i,j,k) grid node
!
! Determine source-receiver traveltimes one at a time.
!
DO i=1,nrc
! 
!  The first step is to locate the receiver in the grid.
!
   irr=INT((gor-rcr(i))/dnr)+1
   irx=INT((rcx(i)-gox)/dnx)+1
   irz=INT((rcz(i)-goz)/dnz)+1
   sw=0
   IF(irr.lt.1.or.irr.gt.nnr)sw=1
   IF(irx.lt.1.or.irx.gt.nnx)sw=1
   IF(irz.lt.1.or.irz.gt.nnz)sw=1
   IF(sw.eq.1)then
      WRITE(6,*)"Receiver lies outside model (rr,xr,zr)= ",irr,irx,irz
   ENDIF
   IF(irr.eq.nnr)irr=irr-1
   IF(irx.eq.nnx)irx=irx-1
   IF(irz.eq.nnz)irz=irz-1
!
!  Location of receiver successfully found within the grid. Now approximate
!  traveltime at receiver using trilinear interpolation from eight
!  surrounding grid points
!
   drx=(rcx(i)-gox)-(irx-1)*dnx
   drz=(rcz(i)-goz)-(irz-1)*dnz
   drr=(gor-rcr(i))-(irr-1)*dnr
   trr=0.0
   DO j=1,2
      DO k=1,2
         DO l=1,2
            produ=(1.0-ABS(((l-1)*dnz-drz)/dnz))*(1.0-ABS(((k-1)*dnx-drx)/dnx))
            produ=produ*(1.0-ABS(((j-1)*dnr-drr)/dnr))
            trr=trr+ttn(irz-1+l,irx-1+k,irr-1+j)*produ
         ENDDO
      ENDDO
   ENDDO
   WRITE(10,*)trr
ENDDO
END SUBROUTINE srtimes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine calculates ray path geometries for each
! receiver from the base of the model. It will also compute
! Frechet derivatives using these ray paths if required.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE rpaths(wrgf,csid,awrgf,caid,cfd)
USE globalp
IMPLICIT NONE
INTEGER,  PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
INTEGER :: i,j,k,l,m,n,ipr,ipx,ipz,nrp,sw
INTEGER :: wrgf,cfd,csid,ipro,ipxo,ipzo,caid,awrgf
INTEGER :: ivr,ivx,ivz,ivro,ivxo,ivzo,nhp
INTEGER :: ivrt,ivxt,ivzt,iprt,ipxt,ipzt,isum
INTEGER, DIMENSION (4) :: chp
INTEGER, PARAMETER :: maxrp=100000
REAL(KIND=i5), PARAMETER :: ftol=1.0e-6
REAL(KIND=i5) :: rayr,rayx,rayz
REAL(KIND=i10) :: dpl,crad,rd1,rd2,atio,ri,xi,zi,vel,velo
REAL(KIND=i10) :: u,v,w,rigz,rigx,rigr,dinc
REAL(KIND=i10) :: dtr,dtx,dtz,drr,drx,drz,produ
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: rgr,rgx,rgz
REAL(KIND=i5), DIMENSION (:,:,:), ALLOCATABLE :: fdm
REAL(KIND=i10), DIMENSION (4) :: vrat,ui,vi,wi,uio,vio,wio
!
! ipr,ipx,ipz = Coordinates of cell containing current point
! ipro,ipxo,ipzo = Coordinates of previous point
! rgr,rgx,rgz = (r,x,z) coordinates of ray geometry
! ivr,ivx,ivz = Coordinates of B-spline vertex containing current point
! ivro,ivxo,ivzo = Coordinates of previous point
! maxrp = maximum number of ray points
! nrp = number of points to describe ray
! dpl = incremental path length of ray
! crad = current radius of ray point
! atio = ratio for locating ray endpoint
! ri,xi,zi = edge of model coordinates
! dtr,dtx,dtz = components of gradT
! wrgf = Write out raypaths? (<0=all,0=no,>0=souce id)
! awrgf = Array ID if wrgf>0
! cfd = calculate Frechet derivatives? (0=no,1=yes)
! csid = current source id
! fdm = Frechet derivative matrix
! nhp = Number of ray segment-B-spline cell hit points
! vrat = length ratio of ray sub-segment
! chp = pointer to incremental change in r,x or z cell
! drr,drx,drz = distance from reference node of cell
! produ = variable for trilinear interpolation
! vel = velocity at current point
! velo = velocity at previous point
! u,v,w = local variables of r,x,z
! ui,vi,wi = B-spline basis functions at current point
! uio,vio,wio = ui,vi,wi for previous point
! ivrt,ivxt,ivzt = temporary ivr,ivx,ivz values
! rigr,rigx,rigz = end point of sub-segment of ray path
! iprt,ipxt,ipzt = temporary ipr,ipx,ipz values
! dinc = path length of ray sub-segment
! rayr,rayx,rayz = ray path coordinates in single precision
! caid = current receiver array id.
!
! Allocate memory to arrays for storing ray path geometry
!
ALLOCATE(rgr(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgr'
ENDIF
ALLOCATE(rgx(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgx'
ENDIF
ALLOCATE(rgz(nnx*nnr), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: SUBROUTINE rpaths: REAL rgz'
ENDIF
!
! Allocate memory to partial derivative array
!
IF(cfd.EQ.1)THEN
   ALLOCATE(fdm(0:nvz+1,0:nvx+1,0:nvr+1), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: SUBROUTINE gridder: REAL veln'
   ENDIF
ENDIF
!
! Set ray incremental path length equal to grid separation
! in depth.
!
  dpl=dnr
!
! Loop through all the receivers
!
DO i=1,nrc
   IF(cfd.EQ.1)THEN
      fdm=0.0
   ENDIF
!
!  The first step is to locate the receiver in the grid.
!
   ipr=INT((gor-rcr(i))/dnr)+1
   ipx=INT((rcx(i)-gox)/dnx)+1
   ipz=INT((rcz(i)-goz)/dnz)+1
   sw=0
   IF(ipr.lt.1.or.ipr.ge.nnr)sw=1
   IF(ipx.lt.1.or.ipx.ge.nnx)sw=1
   IF(ipz.lt.1.or.ipz.ge.nnz)sw=1
   IF(sw.eq.1)then
      WRITE(6,*)"Receiver lies outside model (rr,xr,zr)= ",ipr,ipx,ipz
   ENDIF
   IF(ipr.eq.nnr)ipr=ipr-1
   IF(ipx.eq.nnx)ipx=ipx-1
   IF(ipz.eq.nnz)ipz=ipz-1
!  
!  First point of the ray path is the receiver
!
   rgr(1)=rcr(i)
   rgx(1)=rcx(i)
   rgz(1)=rcz(i)
!
!  Now trace ray from receiver to "source"
!
   DO j=1,maxrp
!
!     Calculate traveltime gradient vector for current cell
!
      dtr=ttn(ipz,ipx,ipr+1)-ttn(ipz,ipx,ipr)
      dtr=dtr+ttn(ipz+1,ipx,ipr+1)-ttn(ipz+1,ipx,ipr)
      dtr=dtr+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx+1,ipr)
      dtr=dtr+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr)
      dtr=dtr/(4.0*dnr)
      dtx=ttn(ipz,ipx+1,ipr)-ttn(ipz,ipx,ipr)
      dtx=dtx+ttn(ipz+1,ipx+1,ipr)-ttn(ipz+1,ipx,ipr)
      dtx=dtx+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz+1,ipx,ipr+1)
      dtx=dtx+ttn(ipz,ipx+1,ipr+1)-ttn(ipz,ipx,ipr+1)
      dtx=dtx/(4.0*(earth+gor-(ipr-1)*dnr)*dnx)
      dtz=ttn(ipz+1,ipx,ipr)-ttn(ipz,ipx,ipr)
      dtz=dtz+ttn(ipz+1,ipx+1,ipr+1)-ttn(ipz,ipx+1,ipr+1)
      dtz=dtz+ttn(ipz+1,ipx+1,ipr)-ttn(ipz,ipx+1,ipr)
      dtz=dtz+ttn(ipz+1,ipx,ipr+1)-ttn(ipz,ipx,ipr+1)
      dtz=dtz/(4.0*(earth+gor-(ipr-1)*dnr)*SIN(rgx(j))*dnz)
!
!     Calculate the next ray path point
!
      crad=earth+rgr(j)
      rd1=SQRT(dtr**2+dtx**2+dtz**2)
      rgr(j+1)=rgr(j)+dpl*dtr/rd1
      rgx(j+1)=rgx(j)-dpl*dtx/(crad*rd1)
      rgz(j+1)=rgz(j)-dpl*dtz/(crad*SIN(rgx(j))*rd1)
!
!     Determine which cell the new ray endpoint
!     lies in.
!
      ipro=ipr
      ipxo=ipx
      ipzo=ipz
      ipr=INT((gor-rgr(j+1))/dnr)+1
      ipx=INT((rgx(j+1)-gox)/dnx)+1
      ipz=INT((rgz(j+1)-goz)/dnz)+1
      sw=0
      IF(ipr.LT.1.OR.ipr.GE.nnr)sw=1
      IF(ipx.LT.1.OR.ipx.GE.nnx)sw=1
      IF(ipz.LT.1.OR.ipz.GE.nnz)sw=1
!
!     If sw.ne.0 then we are done; find
!     the intersection point of the ray segment
!     with the boundary of the model.
!
      IF(sw.NE.0)THEN
         sw=0
         IF(ipr.LT.1.OR.ipr.GE.nnr)THEN
            IF(ipr.LT.1)THEN
               ri=gor
               ipr=1
            ELSE
               ri=gor-(nnr-1)*dnr
               ipr=nnr-1
            ENDIF
            atio=(ri-rgr(j))/(rgr(j+1)-rgr(j))
            sw=1
         ENDIF
         IF(ipx.LT.1.OR.ipx.GE.nnx)THEN
            IF(ipx.LT.1)THEN
               xi=gox
               ipx=1
            ELSE
               xi=(gox+(nnx-1)*dnx)
               ipx=nnx-1
            ENDIF
            rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
            IF(sw.eq.1)then
               IF(rd1.LT.atio)THEN
                  atio=rd1
                  sw=2
               ENDIF
            ELSE
               atio=rd1
               sw=2
            ENDIF
         ENDIF
         IF(ipz.LT.1.OR.ipz.GE.nnz)THEN
            IF(ipz.LT.1)THEN
               zi=goz
               ipz=1
            ELSE
               zi=goz+(nnz-1)*dnz
               ipz=nnz-1
            ENDIF
            rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
            IF(sw.NE.0)then
               IF(rd1.LT.atio)THEN
                  atio=rd1
                  sw=3
               ENDIF
            ELSE
               atio=rd1
               sw=3
            ENDIF
         ENDIF
         IF(sw.EQ.1)THEN
            ipx=ipxo
            ipz=ipzo
         ELSE IF(sw.EQ.2)THEN
            ipr=ipro
            ipz=ipzo
         ELSE IF(sw.EQ.3)THEN
            ipr=ipro
            ipx=ipxo
         ENDIF
         rgz(j+1)=rgz(j)+atio*(rgz(j+1)-rgz(j))
         rgx(j+1)=rgx(j)+atio*(rgx(j+1)-rgx(j))
         rgr(j+1)=rgr(j)+atio*(rgr(j+1)-rgr(j))
         nrp=j+1
         IF(cfd.EQ.1)THEN
            sw=1
         ELSE
            EXIT
         ENDIF
      ENDIF
!
!     Calculate the Frechet derivatives if required.
!
      IF(cfd.EQ.1)THEN
!
!        First determine which B-spline cell the refined cells
!        containing the ray path segment lies in. If they lie
!        in more than one, then we need to divide the problem
!        into separate parts (up to four).
!
         ivr=INT((ipr-1)/nrfr)+1
         ivx=INT((ipx-1)/nrfx)+1
         ivz=INT((ipz-1)/nrfz)+1
         ivro=INT((ipro-1)/nrfr)+1
         ivxo=INT((ipxo-1)/nrfx)+1
         ivzo=INT((ipzo-1)/nrfz)+1
!
!        Calculate up to three hit points between straight
!        ray segment and cell faces.
!
         nhp=0
         IF(ivr.NE.ivro)THEN
            nhp=nhp+1
            IF(ivr.GT.ivro)THEN
               ri=gor-(ivr-1)*dvr
            ELSE
               ri=gor-ivr*dvr
            ENDIF
            vrat(nhp)=(ri-rgr(j))/(rgr(j+1)-rgr(j))
            chp(nhp)=1
         ENDIF
         IF(ivx.NE.ivxo)THEN
            nhp=nhp+1
            IF(ivx.GT.ivxo)THEN
               xi=gox+(ivx-1)*dvx
            ELSE
               xi=gox+ivx*dvx
            ENDIF
            rd1=(xi-rgx(j))/(rgx(j+1)-rgx(j))
            IF(nhp.EQ.1)THEN
               vrat(nhp)=rd1
               chp(nhp)=2
            ELSE
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=2
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=2
               ENDIF
            ENDIF
         ENDIF
         IF(ivz.NE.ivzo)THEN
            nhp=nhp+1
            IF(ivz.GT.ivzo)THEN
               zi=goz+(ivz-1)*dvz 
            ELSE
               zi=goz+ivz*dvz
            ENDIF
            rd1=(zi-rgz(j))/(rgz(j+1)-rgz(j))
            IF(nhp.EQ.1)THEN
               vrat(nhp)=rd1
               chp(nhp)=3
            ELSE IF(nhp.EQ.2)THEN
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=3
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=3
               ENDIF
            ELSE
               IF(rd1.GE.vrat(nhp-1))THEN
                  vrat(nhp)=rd1
                  chp(nhp)=3
               ELSE IF(rd1.GE.vrat(nhp-2))THEN
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=rd1
                  chp(nhp-1)=3
               ELSE
                  vrat(nhp)=vrat(nhp-1)
                  chp(nhp)=chp(nhp-1)
                  vrat(nhp-1)=vrat(nhp-2)
                  chp(nhp-1)=chp(nhp-2)
                  vrat(nhp-2)=rd1
                  chp(nhp-2)=3
               ENDIF
            ENDIF
         ENDIF
         nhp=nhp+1
         vrat(nhp)=1.0
         chp(nhp)=0
!
!        Calculate the velocity, u,v and w values of the
!        first point
!
         drr=(gor-rgr(j))-(ipro-1)*dnr
         drx=(rgx(j)-gox)-(ipxo-1)*dnx
         drz=(rgz(j)-goz)-(ipzo-1)*dnz
         vel=0.0
         DO k=1,2
            DO l=1,2
               DO m=1,2
                  produ=(1.0-ABS(((m-1)*dnz-drz)/dnz))
                  produ=produ*(1.0-ABS(((l-1)*dnx-drx)/dnx))
                  produ=produ*(1.0-ABS(((k-1)*dnr-drr)/dnr))
         IF(ipzo-1+m.LE.nnz.AND.ipxo-1+l.LE.nnx.AND.ipro-1+k.LE.nnr)THEN
                  vel=vel+veln(ipzo-1+m,ipxo-1+l,ipro-1+k)*produ
         ENDIF
               ENDDO
            ENDDO
         ENDDO
         drr=(gor-rgr(j))-(ivro-1)*dvr
         drx=(rgx(j)-gox)-(ivxo-1)*dvx
         drz=(rgz(j)-goz)-(ivzo-1)*dvz
         u=drr/dvr
         v=drx/dvx
         w=drz/dvz
!
!        Calculate the 12 basis values at the point
!
         ui(1)=(1.0-u)**3/6.0
         ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
         ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
         ui(4)=u**3/6.0
         vi(1)=(1.0-v)**3/6.0
         vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
         vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
         vi(4)=v**3/6.0
         wi(1)=(1.0-w)**3/6.0
         wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
         wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
         wi(4)=w**3/6.0
         ivrt=ivro
         ivxt=ivxo
         ivzt=ivzo
!
!        Now loop through the one or more sub-segments of the
!        ray path segment and calculate partial derivatives
!
         DO k=1,nhp
            velo=vel
            uio=ui
            vio=vi
            wio=wi
            IF(k.GT.1)THEN
               IF(chp(k-1).EQ.1)THEN
                  ivrt=ivr
               ELSE IF(chp(k-1).EQ.2)THEN
                  ivxt=ivx
               ELSE IF(chp(k-1).EQ.3)THEN
                  ivzt=ivz
               ENDIF
            ENDIF
!
!           Calculate the velocity, u,v and w values of the
!           new point
!
            rigz=rgz(j)+vrat(k)*(rgz(j+1)-rgz(j))
            rigx=rgx(j)+vrat(k)*(rgx(j+1)-rgx(j))
            rigr=rgr(j)+vrat(k)*(rgr(j+1)-rgr(j))
            iprt=INT((gor-rigr)/dnr)+1
            ipxt=INT((rigx-gox)/dnx)+1
            ipzt=INT((rigz-goz)/dnz)+1
            drr=(gor-rigr)-(iprt-1)*dnr
            drx=(rigx-gox)-(ipxt-1)*dnx
            drz=(rigz-goz)-(ipzt-1)*dnz
            vel=0.0
            DO l=1,2
               DO m=1,2
                  DO n=1,2
                     produ=(1.0-ABS(((n-1)*dnz-drz)/dnz))
                     produ=produ*(1.0-ABS(((m-1)*dnx-drx)/dnx))
                     produ=produ*(1.0-ABS(((l-1)*dnr-drr)/dnr))
            IF(ipzt-1+n.LE.nnz.AND.ipxt-1+m.LE.nnx.AND.iprt-1+l.LE.nnr)THEN
                     vel=vel+veln(ipzt-1+n,ipxt-1+m,iprt-1+l)*produ
            ENDIF
                  ENDDO
               ENDDO
            ENDDO
            drr=(gor-rigr)-(ivrt-1)*dvr
            drx=(rigx-gox)-(ivxt-1)*dvx
            drz=(rigz-goz)-(ivzt-1)*dvz
            u=drr/dvr
            v=drx/dvx
            w=drz/dvz
!
!           Calculate the 12 basis values at the new point
!
            ui(1)=(1.0-u)**3/6.0
            ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
            ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
            ui(4)=u**3/6.0
            vi(1)=(1.0-v)**3/6.0
            vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
            vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
            vi(4)=v**3/6.0
            wi(1)=(1.0-w)**3/6.0
            wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
            wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
            wi(4)=w**3/6.0
!
!           Calculate the incremental path length
!     
            IF(k.EQ.1)THEN
               dinc=vrat(k)*dpl
            ELSE 
               dinc=(vrat(k)-vrat(k-1))*dpl
            ENDIF
!
!           Now compute the 64 contributions to the partial
!           derivatives.
!
            DO l=1,4
               DO m=1,4
                  DO n=1,4
                     rd1=ui(n)*vi(m)*wi(l)/vel**2
                     rd2=uio(n)*vio(m)*wio(l)/velo**2
                     rd1=-(rd1+rd2)*dinc/2.0
                     rd2=fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)
                     fdm(ivzt-2+l,ivxt-2+m,ivrt-2+n)=rd1+rd2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF(sw.EQ.1)EXIT
      ENDIF
   ENDDO
!
!  Write ray paths to output file
!
   IF((wrgf.EQ.csid.AND.awrgf.EQ.caid).OR.wrgf.LT.0)THEN
      WRITE(40)nrp
      DO j=1,nrp
         rayr=rgr(j)
         rayx=(pi/2-rgx(j))*180.0/pi
         rayz=rgz(j)*180.0/pi
         WRITE(40)rayr,rayx,rayz
      ENDDO
   ENDIF
!
!  Write partial derivatives to output file
!
   IF(cfd.EQ.1)THEN
!
!     Determine the number of non-zero elements.
!
      isum=0
      DO j=0,nvz+1
         DO k=0,nvx+1
            DO l=0,nvr+1
               IF(ABS(fdm(j,k,l)).GE.ftol)isum=isum+1
            ENDDO
         ENDDO
      ENDDO
      WRITE(50)isum
      isum=0
      DO j=0,nvz+1
         DO k=0,nvx+1
            DO l=0,nvr+1
               isum=isum+1
               IF(ABS(fdm(j,k,l)).GE.ftol)WRITE(50)isum,fdm(j,k,l)
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO
IF(cfd.EQ.1)THEN
   DEALLOCATE(fdm, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: fdm'
   ENDIF
ENDIF
DEALLOCATE(rgr,rgx,rgz, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE rpaths: rgr,rgx,rgz'
ENDIF
END SUBROUTINE rpaths

