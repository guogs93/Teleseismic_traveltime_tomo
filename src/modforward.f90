!******************************************
! Global paramater for froward simulation
!******************************************

MODULE globalfwd
 INTEGER, ALLOCATABLE:: nsts(:,:,:), lock(:,:,:)
 REAL*8, ALLOCATABLE:: velc(:,:,:), ttn(:,:,:)
 REAL*8, PARAMETER:: big=1.0d10
 INTEGER, DIMENSION(8):: imin,imax,istep,jmin,jmax,jstep,kmin,kmax,kstep
 
 !***************************************************
 ! nsts = status of calculate data 
 ! ttn = the temporary varible for traveltime field
 ! velc = the refined grids for forward model 
 ! imin, imax, istep = sweeping order 
 !***************************************************

END MODULE

!*************************************************************************
! Calculate traveltime field and data at receivers for difference sources
!*************************************************************************

SUBROUTINE subforward(velt,veln,ttc,dcal)
 USE globalp
 USE globalfwd
 INTEGER:: ishot,i,j,k
 REAL*8:: velt(0:nvp+1,0:nvt+1,0:nvr+1),dcal(1:nsrc,1:nrc),ttc(nsrc,nnp,nnt,nnr),veln(nnp,nnt,nnr)
 REAL*8, ALLOCATABLE:: dcalnoise(:,:)
 real:: dum
 !*******************************************************
 ! input 
 ! velt = current velocity model
 !
 ! output
 ! dcal = calculated traveltime data at receiver
 ! ttc = the traveltime field for different sources
 ! srcm = model traveltime residual mean for each source
 ! ntsrc = number of traveltimes for each source
 !*******************************************************
 
 ALLOCATE(ttn(nnp,nnt,nnr),velc(nnp,nnt,nnr),nsts(nnp,nnt,nnr),lock(nnp,nnt,nnr))

 !*********************************************
 ! Refine current velocity model velt -> velc
 !*********************************************
 
 CALL gridder(velt)
 
 !***************************************************
 ! read in crustal model 
 !***************************************************
 if(crt.eq.1)then
 open(unit=11,file=crtmodel,status='old')
 do k=1, ncrt
 do j=nnt, 1, -1
 do i=1, nnp
    read(11,*)dum,dum,dum,velc(i,j,k)
 enddo
 enddo
 enddo
 close(11)
 endif

 veln=velc
 
 !****************************************************
 ! Open file containing traveltimes at base of model
 !****************************************************
 
 OPEN(UNIT=20,FILE=itimes,FORM='unformatted',STATUS='old')
   
 DO ishot=1,  nsrc

    !************************************************************************
    ! Call a subroutine that works out the first-arrival traveltime field.
    !************************************************************************

    CALL travel(ishot)

    !**********************************************
    ! store traveltime field for different sources
    !**********************************************
   
    DO j=1, nnp
    DO k=1, nnt
    DO l=1, nnr
       ttc(ishot,j,k,l)=ttn(j,k,l)
    ENDDO
    ENDDO
    ENDDO
      
    !*******************************************************
    ! Find traveltimes at receivers from base of grid
    !*******************************************************
 
    CALL srtimes(ishot,dcal(ishot,:))
    
 ENDDO
 CLOSE(20)

 if(mode.eq.0)then
   open(unit=111,file='rtravel.out')
   do i=1, nsrc
   do j=1, nrc
      write(111,*) dcal(i,j) 
   enddo
   enddo
   close(111)
   write(*,*)"compute the reference traveltime finish ! rtravel.out !"
   stop
 endif 
 
 !*****************************************************
 ! Remove mean from model residuals if required.
 !*****************************************************

 IF(rmtr.EQ.1)THEN
    srcm=0.0
    ntsrc=0
    DO i=1, nsrc
    DO j=1, nrc
       IF(tstat(i,j).EQ.1)THEN
         srcm(i)=srcm(i)+(dcal(i,j)-dref(i,j))
         ntsrc(i)=ntsrc(i)+1
       ENDIF
    ENDDO
    ENDDO
    
    DO i=1, nsrc
    DO j=1, nrc
       IF(tstat(i,j).EQ.1)THEN
         dcal(i,j)=(dcal(i,j)-dref(i,j))-srcm(i)/REAL(ntsrc(i))
       ELSE
         dcal(i,j)=0
       ENDIF
    ENDDO
    ENDDO
   
 ENDIF
 
 if (mode.eq.1)then

    allocate(dcalnoise(nsrc,nrc))

    dcalnoise=dcal

    ! add the Gaussian noise to the data 
    if(agn.eq.1)then
      write(*,*) "add the Gaussian noise to the synthetic data", sdgn 
      call gaussnoise(dcal,dcalnoise)
    endif

    open(unit=111,file='otimes.ckb.dat')
    do i=1, nsrc
    do j=1, nrc
       write(111,'(I2,1X,F8.5,1X,F7.4)')tstat(i,j),dcalnoise(i,j),cd(i,j)
    enddo
    enddo
    close(111)

    deallocate(dcalnoise)

    stop

 endif
  
 DEALLOCATE (ttn,velc,nsts,lock)

END SUBROUTINE subforward

SUBROUTINE init_sweep
 USE globalfwd
 USE globalp

 imin(1)=1;imax(1)=nnr;istep(1)=1
 jmin(1)=1;jmax(1)=nnt;jstep(1)=1
 kmin(1)=1;kmax(1)=nnp;kstep(1)=1

 imin(2)=1;imax(2)=nnr;istep(2)=1
 jmin(2)=nnt;jmax(2)=1;jstep(2)=-1
 kmin(2)=1;kmax(2)=nnp;kstep(2)=1

 imin(3)=1;imax(3)=nnr;istep(3)=1
 jmin(3)=1;jmax(3)=nnt;jstep(3)=1
 kmin(3)=nnp;kmax(3)=1;kstep(3)=-1

 imin(4)=1;imax(4)=nnr;istep(4)=1
 jmin(4)=nnt;jmax(4)=1;jstep(4)=-1
 kmin(4)=nnp;kmax(4)=1;kstep(4)=-1 

 imin(5)=nnr;imax(5)=1;istep(5)=-1
 jmin(5)=1;jmax(5)=nnt;jstep(5)=1
 kmin(5)=1;kmax(5)=nnp;kstep(5)=1

 imin(6)=nnr;imax(6)=1;istep(6)=-1
 jmin(6)=1;jmax(6)=nnt;jstep(6)=1
 kmin(6)=nnp;kmax(6)=1;kstep(6)=-1

 imin(7)=nnr;imax(7)=1;istep(7)=-1
 jmin(7)=nnt;jmax(7)=1;jstep(7)=-1
 kmin(7)=nnp;kmax(7)=1;kstep(7)=-1

 imin(8)=nnr;imax(8)=1;istep(8)=-1
 jmin(8)=nnt;jmax(8)=1;jstep(8)=-1
 kmin(8)=1;kmax(8)=nnp;kstep(8)=1

END SUBROUTINE init_sweep 

!**********************************************************************************
! This subroutine transforms the inversed velocity model to refined velocity model
!**********************************************************************************

SUBROUTINE gridder(velt)
 USE globalp
 USE globalfwd
 REAL*8:: velt(0:nvp+1,0:nvt+1,0:nvr+1)
 INTEGER :: m,n,i1,j1,k1,i,j,k
 INTEGER :: conr,cont,conp,str,stt,stp
 REAL*8 :: u,sumi,sumj,sumk
 REAL*8, ALLOCATABLE :: ui(:,:),vi(:,:),wi(:,:)
 !
 ! u = Cubic spline independent variable
 ! ui,vi,wi = Cubic spline basis functions
 ! sumi,sumj,sumk = Summation variables for constructing spline
 ! conr,cont,conp = Counters for refining grid in r,theta,phi
 ! str,stt,stp = Refined grid location in r,theta,phi
 !
 ! Calculate the values of the basis functions
 !
 ! The basic function for p, z
 ALLOCATE(ui(nvp+1,4))
 DO i=1,nrfp+1
    u=nrfp
    u=(i-1)/u
    ui(i,1)=(1.0-u)**3/6.0
    ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
    ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
    ui(i,4)=u**3/6.0
 ENDDO
 
 ! The basic function for t, x 
 ALLOCATE(vi(nvt+1,4))
 DO i=1,nrft+1
    u=nrft
    u=(i-1)/u
    vi(i,1)=(1.0-u)**3/6.0
    vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
    vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
    vi(i,4)=u**3/6.0
 ENDDO
 
 ! The basic function for r, r 
 ALLOCATE(wi(nvr+1,4))
 DO i=1,nrfr+1
    u=nrfr
    u=(i-1)/u
    wi(i,1)=(1.0-u)**3/6.0
    wi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
    wi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
    wi(i,4)=u**3/6.0
 ENDDO

 !
 ! Calculate velocity values on refined grid
 !
 DO i=1,nvp-1
    conp=nrfp
    IF(i==nvp-1)conp=nrfp+1
    DO j=1,nvt-1
       cont=nrft
       IF(j==nvt-1)cont=nrft+1
       DO k=1,nvr-1
          conr=nrfr
          IF(k==nvr-1)conr=nrfr+1
          DO l=1,conp
             stp=nrfp*(i-1)+l
             DO m=1,cont
                stt=nrft*(j-1)+m
                DO n=1,conr
                   str=nrfr*(k-1)+n
                   sumi=0.0
                   DO i1=1,4
                      sumj=0.0
                      DO j1=1,4
                         sumk=0.0
                         DO k1=1,4
                            sumk=sumk+wi(n,k1)*velt(i-2+i1,j-2+j1,k-2+k1)
                         ENDDO
                         sumj=sumj+vi(m,j1)*sumk
                      ENDDO
                      sumi=sumi+ui(l,i1)*sumj
                   ENDDO
                   velc(stp,stt,str)=sumi
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 DEALLOCATE(ui,vi,wi)
 
END SUBROUTINE gridder

!*****************************************************************
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using fast sweeping method.
!*****************************************************************

SUBROUTINE travel(srcid)
 USE globalp
 USE globalfwd
 INTEGER :: sw,srcid,i,j
 INTEGER :: ite
 REAL*8 :: rd1, cvg_err
 REAL*8, ALLOCATABLE :: ttno(:,:,:) 
 !
 ! rd1 = substitution variable
 ! srcid = id number of source
 ! sw = switch (0 or 1)
 !
 ! Allocate nsts and set all elements of array nsts equal to 
 ! (-1) since to begin with, all points are "big" points.
 !
 nsts=-1
 lock=1
 ALLOCATE(ttno(nnp,nnt,nnr))
 !
 ! Read in traveltimes calculated a priori to base of grid, find point with minimum traveltime,
 !
 sw=0
 ttn=big
 DO i=1, nnp
    DO j=1, nnt
       READ(20)rd1
       IF(rd1.GT.0.0)THEN
          ttn(i,j,nnr)=rd1
          nsts(i,j,nnr)=0
          if((i-1).ge.1)lock(i-1,j,nnr)=0
          if((i+1).le.nnp)lock(i+1,j,nnr)=0
          if((j-1).ge.1)lock(i,j-1,nnr)=0
          if((j+1).le.nnt)lock(i,j+1,nnr)=0
          lock(i,j,nnr-1)=0
          sw=1
       ENDIF
    ENDDO
 ENDDO
 IF(sw.EQ.0)THEN
    WRITE(6,*)'No initial traveltimes for source ',srcid
    WRITE(6,*)'Terminating program'
    STOP
 ENDIF

 CALL init_sweep

 ite=1
 cvg_err=1.0d10
 DO WHILE((cvg_err.GT.1.0d-12).AND.(ite.LE.maxiter))
    ttno=ttn
    IF(fom.EQ.1)THEN
      CALL sweep1
    ELSEIF(fom.EQ.2)THEN
      CALL sweep2
    ELSEIF(fom.EQ.3)THEN
      CALL sweep3
    ENDIF
    ite=ite+1
    cvg_err=sum(abs(ttno-ttn))/dble(nnr*nnt*nnp) 
    !WRITE(6,*)'cvg_err=',cvg_err
 ENDDO

 DEALLOCATE(ttno)

END SUBROUTINE travel

!****************************************************************
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 1-order fast sweeping method.
!****************************************************************

SUBROUTINE sweep1
 USE globalfwd
 USE globalp
 INTEGER :: swp,i,j,k
 REAL*8 :: txmin,tzmin,trmin
 REAL*8 :: a,b,c,tbar,slown,ri,risti
 REAL*8 :: ww,vv,uu,hw,hv,hu

 DO swp=4, 8
    DO i=imin(swp),imax(swp),istep(swp)
    DO j=jmin(swp),jmax(swp),jstep(swp)
    DO k=kmin(swp),kmax(swp),kstep(swp)  
       IF(nsts(k,j,i).NE.0.and.lock(k,j,i).eq.0)THEN       
         slown=1.0/velc(k,j,i)
         ri=gor-(i-1)*dnr+earth
         risti=ri*sin(cgot+(j-1)*dnt)

         IF(i.GT.1.AND.i.LT.nnr)THEN
            trmin=MIN(ttn(k,j,i+1),ttn(k,j,i-1))
         ELSEIF(i.EQ.1)THEN
            trmin=ttn(k,j,i+1)           
         ELSE
            trmin=ttn(k,j,i-1)
         ENDIF      

         IF(j.GT.1.AND.j.LT.nnt)THEN
            txmin=MIN(ttn(k,j+1,i),ttn(k,j-1,i))      
         ELSEIF(j.EQ.1)THEN
            txmin=ttn(k,j+1,i)
         ELSE
            txmin=ttn(k,j-1,i)
         ENDIF

         IF(k.GT.1.AND.k.LT.nnp)THEN
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
         hv=ri*dnt
         hu=risti*dnp

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

         !ttn(k,j,i)=min(tbar,ttn(k,j,i))
         if(tbar.lt.ttn(k,j,i))then
            ttn(k,j,i)=tbar
            if(k+1.le.nnp)then
               if(ttn(k+1,j,i).gt.ttn(k,j,i))lock(k+1,j,i)=0
            endif
            if(k-1.ge.1)then
              if(ttn(k-1,j,i).gt.ttn(k,j,i))lock(k-1,j,i)=0
            endif
            if(j+1.le.nnt)then
              if(ttn(k,j+1,i).gt.ttn(k,j,i))lock(k,j+1,i)=0
            endif
            if(j-1.ge.1)then
              if(ttn(k,j-1,i).gt.ttn(k,j,i))lock(k,j-1,i)=0
            endif
            if(i+1.le.nnr)then
              if(ttn(k,j,i+1).gt.ttn(k,j,i))lock(k,j,i+1)=0
            endif
            if(i-1.ge.1)then
              if(ttn(k,j,i-1).gt.ttn(k,j,i))lock(k,j,i-1)=0
            endif
         endif
         lock(k,j,i)=1
      ENDIF
   ENDDO
   ENDDO
   ENDDO
 ENDDO
END SUBROUTINE sweep1

!*****************************************************************
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 2-order fast sweeping method.
!*****************************************************************

SUBROUTINE sweep2
 USE globalfwd
 USE globalp
 INTEGER :: swp,i,j,k
 REAL*8 :: txmin,tzmin,trmin
 REAL*8 :: a,b,c,tbar,slown,ri,risti
 REAL*8 :: ww,vv,uu,hw,hv,hu,hr,hx,hz

 DO swp=4, 8
    DO i=imin(swp),imax(swp),istep(swp)
    DO j=jmin(swp),jmax(swp),jstep(swp)
    DO k=kmin(swp),kmax(swp),kstep(swp)  
       IF(nsts(k,j,i).NE.0.and.lock(k,j,i).eq.0)THEN       
         slown=1.0/velc(k,j,i)
         ri=gor-(i-1)*dnr+earth
         risti=ri*sin(cgot+(j-1)*dnt) 

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

         IF(j.GT.1.AND.j.LT.nnt)THEN
            txmin=MIN(ttn(k,j-1,i),ttn(k,j+1,i))
            hx=dnt
            IF(txmin.EQ.ttn(k,j-1,i))THEN
               IF(j-2.GE.1)THEN
                  IF(ttn(k,j-2,i).LT.ttn(k,j-1,i))THEN
                    txmin=ttn(k,j-1,i)+(ttn(k,j-1,i)-ttn(k,j-2,i))/3
                    hx=dnt*2/3
                 ENDIF
              ENDIF
            ELSE
              IF(j+2.LE.nnt)THEN
                 IF(ttn(k,j+2,i).LT.ttn(k,j+1,i))then
                    txmin=ttn(k,j+1,i)+(ttn(k,j+1,i)-ttn(k,j+2,i))/3
                    hx=dnt*2/3
                 ENDIF
              ENDIF
           ENDIF    
         ELSEIF(j.EQ.1)THEN
           txmin=ttn(k,j+1,i)
           hx=dnt
         ELSEIF(j.EQ.nnt)THEN
           txmin=ttn(k,j-1,i)
           hx=dnt
         ENDIF

         IF(k.GT.1.AND.k.LT.nnp)THEN
           tzmin=MIN(ttn(k-1,j,i),ttn(k+1,j,i))
           hz=dnp
           IF(tzmin.EQ.ttn(k-1,j,i))THEN
              IF(k-2.GE.1)THEN
                IF(ttn(k-2,j,i).LT.ttn(k-1,j,i))THEN
                   tzmin=ttn(k-1,j,i)+(ttn(k-1,j,i)-ttn(k-2,j,i))/3
                   hz=dnp*2/3
                ENDIF
              ENDIF
           ELSE
              IF(k+2.LE.nnp)THEN
                IF(ttn(k+2,j,i).LT.ttn(k+1,j,i))THEN
                   tzmin=ttn(k+1,j,i)+(ttn(k+1,j,i)-ttn(k+2,j,i))/3
                   hz=dnp*2/3
                ENDIF
              ENDIF
           ENDIF 
        ELSEIF(k.EQ.1)THEN
           tzmin=ttn(k+1,j,i)
           hz=dnp
        ELSEIF(k.EQ.nnp)THEN
           tzmin=ttn(k-1,j,i)
           hz=dnp
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

          if(tbar.lt.ttn(k,j,i))then
            ttn(k,j,i)=tbar
            if(k+1.le.nnp)then
               if(ttn(k+1,j,i).gt.ttn(k,j,i))lock(k+1,j,i)=0
            endif
            if(k-1.ge.1)then
              if(ttn(k-1,j,i).gt.ttn(k,j,i))lock(k-1,j,i)=0
            endif
            if(j+1.le.nnt)then
              if(ttn(k,j+1,i).gt.ttn(k,j,i))lock(k,j+1,i)=0
            endif
            if(j-1.ge.1)then
              if(ttn(k,j-1,i).gt.ttn(k,j,i))lock(k,j-1,i)=0
            endif
            if(i+1.le.nnr)then
              if(ttn(k,j,i+1).gt.ttn(k,j,i))lock(k,j,i+1)=0
            endif
            if(i-1.ge.1)then
              if(ttn(k,j,i-1).gt.ttn(k,j,i))lock(k,j,i-1)=0
            endif
         endif
         lock(k,j,i)=1
          
       ENDIF
   ENDDO
   ENDDO
   ENDDO
 ENDDO
END SUBROUTINE sweep2

!******************************************************************
! This subroutine is passed the location of a source, and from
! this point the first-arrival traveltime field through the
! velocity grid is determined using 3-order fast sweeping method.
!******************************************************************

SUBROUTINE sweep3
 USE globalfwd
 USE globalp
 INTEGER :: swp,i,j,k
 REAL*8 :: txmin,tzmin,trmin
 REAL*8 :: a,b,c,tbar,slown,ri,risti
 REAL*8 :: ww,vv,uu,hw,hv,hu,hr,hx,hz

 DO swp=4, 8
    DO i=imin(swp),imax(swp),istep(swp)
    DO j=jmin(swp),jmax(swp),jstep(swp)
    DO k=kmin(swp),kmax(swp),kstep(swp)  
       IF(nsts(k,j,i).NE.0)THEN       
         slown=1.0/velc(k,j,i)
         ri=gor-(i-1)*dnr+earth
         risti=ri*sin(cgot+(j-1)*dnt)

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

        IF(j.GT.1.AND.j.LT.nnt)THEN
           txmin=MIN(ttn(k,j-1,i),ttn(k,j+1,i))
           hx=dnt
           IF(txmin.EQ.ttn(k,j-1,i))THEN
              IF(j-2.GE.1)THEN
                 IF(ttn(k,j-2,i).LT.ttn(k,j-1,i))THEN
                    txmin=ttn(k,j-1,i)+(ttn(k,j-1,i)-ttn(k,j-2,i))/3
                    hx=dnt*2/3
                    IF(j-3.GE.1)THEN
                       IF(ttn(k,j-3,i).LT.ttn(k,j-2,i))THEN
                          txmin=ttn(k,j-1,i)+(7*ttn(k,j-1,i)-9*ttn(k,j-2,i)+2*ttn(k,j-3,i))/11
                          hx=dnt*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ELSE
              IF(j+2.LE.nnt)THEN
                 IF(ttn(k,j+2,i).LT.ttn(k,j+1,i))then
                    txmin=ttn(k,j+1,i)+(ttn(k,j+1,i)-ttn(k,j+2,i))/3
                    hx=dnt*2/3
                    IF(j+3.LE.nnt)THEN
                       IF(ttn(k,j+3,i).LT.ttn(k,j+2,i))THEN
                          txmin=ttn(k,j+1,i)+(7*ttn(k,j+1,i)-9*ttn(k,j+2,i)+2*ttn(k,j+3,i))/11
                          hx=dnt*6/11
                       ENDIF               
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF    
        ELSEIF(j.EQ.1)THEN
           txmin=ttn(k,j+1,i)
           hx=dnt
        ELSEIF(j.EQ.nnt)THEN
           txmin=ttn(k,j-1,i)
           hx=dnt
        ENDIF

        IF(k.GT.1.AND.k.LT.nnp)THEN
           tzmin=MIN(ttn(k-1,j,i),ttn(k+1,j,i))
           hz=dnp
           IF(tzmin.EQ.ttn(k-1,j,i))THEN
              IF(k-2.GE.1)THEN
                IF(ttn(k-2,j,i).LT.ttn(k-1,j,i))THEN
                   tzmin=ttn(k-1,j,i)+(ttn(k-1,j,i)-ttn(k-2,j,i))/3
                   hz=dnp*2/3
                   IF(k-3.GE.1)THEN
                       IF(ttn(k-3,j,i).LT.ttn(k-2,j,i))THEN
                          tzmin=ttn(k-1,j,i)+(7*ttn(k-1,j,i)-9*ttn(k-2,j,i)+2*ttn(k-3,j,i))/11
                          hz=dnp*6/11
                       ENDIF               
                    ENDIF
                ENDIF
              ENDIF
           ELSE
              IF(k+2.LE.nnp)THEN
                IF(ttn(k+2,j,i).LT.ttn(k+1,j,i))THEN
                   tzmin=ttn(k+1,j,i)+(ttn(k+1,j,i)-ttn(k+2,j,i))/3
                   hz=dnp*2/3
                   IF(k+3.LE.nnp)THEN
                       IF(ttn(k+3,j,i).LT.ttn(k+2,j,i))THEN
                          tzmin=ttn(k+1,j,i)+(7*ttn(k+1,j,i)-9*ttn(k+2,j,i)+2*ttn(k+3,j,i))/11
                          hz=dnp*6/11
                       ENDIF               
                    ENDIF
                ENDIF
              ENDIF
           ENDIF 
        ELSEIF(k.EQ.1)THEN
           tzmin=ttn(k+1,j,i)
           hz=dnp
        ELSEIF(k.EQ.nnp)THEN
           tzmin=ttn(k-1,j,i)
           hz=dnp
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

!****************************************************************
! This subroutine is used to swap the values of two traveltimes 
!****************************************************************

SUBROUTINE swap(a,b)
 REAL*8 :: a,b,tem

 tem=a
 a=b
 b=tem

END SUBROUTINE

!******************************************************************
! This subroutine calculates all receiver traveltimes for
! a given source and writes the results to file.
!******************************************************************

SUBROUTINE srtimes(srcid,dmodc)
 USE globalfwd
 USE globalp
 REAL*8 :: dmodc(1:nrc)
 INTEGER :: irr,irt,irp,sw,srcid,i,j,k,l
 REAL :: trr
 REAL*8 :: drr,drt,drp,produ
 !
 ! irr,irt,irp = Coordinates of cell containing receiver
 ! trr = traveltime value at receiver
 ! produ = dummy multiplier
 ! drr,drt,drp = receiver distance from (i,j,k) grid node
 !
 ! Determine source-receiver traveltimes one at a time.
 !
 DO i=1, nrc
 ! 
 !  The first step is to locate the receiver in the grid.
 !
    irr=INT((gor-rcr(i))/dnr)+1
    irt=INT((rct(i)-cgot)/dnt)+1
    irp=INT((rcp(i)-gop)/dnp)+1
    sw=0

    !write(*,*)irr,irt,irp,i
    IF(irr.lt.1.or.irr.gt.nnr)sw=1
    IF(irt.lt.1.or.irt.gt.nnt)sw=1
    IF(irp.lt.1.or.irp.gt.nnp)sw=1
    IF(sw.eq.1)then
       WRITE(6,*)"Receiver lies outside model (rr,tr,pr)= ",irr,irt,irp
    ENDIF
    IF(irr.eq.nnr)irr=irr-1
    IF(irt.eq.nnt)irt=irt-1
    IF(irp.eq.nnp)irp=irp-1
 !
 !  Location of receiver successfully found within the grid. Now approximate
 !  traveltime at receiver using trilinear interpolation from eight surrounding grid points
 !
    drt=(rct(i)-cgot)-(irt-1)*dnt
    drp=(rcp(i)-gop)-(irp-1)*dnp
    drr=(gor-rcr(i))-(irr-1)*dnr
    trr=0.0
    DO j=1,2
       DO k=1,2
          DO l=1,2
             produ=(1.0-ABS(((l-1)*dnp-drp)/dnp))*(1.0-ABS(((k-1)*dnt-drt)/dnt))
             produ=produ*(1.0-ABS(((j-1)*dnr-drr)/dnr))
             trr=trr+ttn(irp-1+l,irt-1+k,irr-1+j)*produ
          ENDDO
       ENDDO
    ENDDO
    dmodc(i)=trr
 ENDDO

END SUBROUTINE srtimes
