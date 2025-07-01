!*******************************************
! Global parameter for adjoint calculation
!*******************************************

MODULE globaladj

 REAL*8, PARAMETER:: big=1.0d8
 INTEGER, ALLOCATABLE:: nsts(:,:,:)
 REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: ttn,adj
 REAL*8, EXTERNAL:: pos,neg
 INTEGER, DIMENSION(8) :: imin,imax,istep,jmin,jmax,jstep,kmin,kmax,kstep

 !****************************************************************
 ! ttn = the traveltime field for different sources
 ! adj = adjoint-state filed 
 ! nsts = status of node
 ! pos, neg = the function to calculate postive or negtive value
 ! imin, imax .... = sweeping order
 !**************************************************************** 

END MODULE

!***************************************************
! This subroutine computes adjoint-state 
!***************************************************

SUBROUTINE subadjoint(ishot,T,residuals,uu,adjn)
 USE globalp
 USE globaladj
 INTEGER:: ishot,ite,adjn,i,j,k
 REAL*8:: T(nnp,nnt,nnr),uu(nnp,nnt,nnr),residuals(nrc)
 REAL*8, ALLOCATABLE:: adjo(:,:,:)
 REAL*8:: cvg_err,eps,adjmax
 
 ! T = traveltime field for current source
 ! u = adjoint-state field for current source
 ! residuals = teleseismic traveltime residuals
 ! adjo = storage adjoint-state field

 ALLOCATE(adj(nnp,nnt,nnr),adjo(nnp,nnt,nnr))
 ALLOCATE(ttn(0:nnp+1,0:nnt+1,0:nnr+1),nsts(nnp,nnt,nnr))

 ttn=big

 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    ttn(i,j,k)=T(i,j,k)
 ENDDO
 ENDDO
 ENDDO

 !************************************
 ! Initialize sweeping order
 !************************************
 CALL init_sweep2

 !************************************
 ! Initialize ajoint at receivers
 !************************************

 CALL init_receiver(ishot,residuals,adjn)

 eps=1.0d-12
 ite=1   
 cvg_err=100000
      
 DO WHILE((cvg_err.GT.eps).AND.(ite.LE.maxiter))
    adjo=adj
    CALL adjoint_sweep
    ite=ite+1
    cvg_err=sum(abs(adjo(1:nnp,1:nnt,1:nnr)-adj(1:nnp,1:nnt,1:nnr)))/dble(nnr*nnt*nnp)    
 ENDDO

 uu=adj

 adjmax=maxval(abs(uu))
 
 IF(debug.AND.ishot.EQ.1)THEN
 open(unit=888,file='adjoint.dat')
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    WRITE(888,*)gopd+(i-1)*dnpd,gotd-(j-1)*dntd,gor-(k-1)*dnr,uu(i,j,k)/adjmax
 ENDDO
 ENDDO
 ENDDO
 close(888)
 ENDIF

 DEALLOCATE (ttn,nsts,adj,adjo)

END SUBROUTINE subadjoint

!*********************************************************
! This subroutine initializes the boundary conditions.
!*********************************************************

SUBROUTINE init_receiver(ishot,residuals,adjn)
 USE globaladj
 USE globalp
 REAL*8 :: residuals(nrc)
 INTEGER :: adjn,ishot,sw 
 INTEGER :: irr,irt,irp,l                                       
 REAL*8 :: par1,par2
 REAL*8 :: nr,dTr

 adj=0
 nsts=0

 DO l=1, nrc
    IF(tstat(ishot,l).EQ.1)THEN
      irr=INT((gor-rcr(l))/dnr)+1
      irt=INT((rct(l)-cgot)/dnt)+1
      irp=INT((rcp(l)-gop)/dnp)+1
      sw=0
      IF(irr.LT.1.OR.irr.GT.nnr)sw=1
      IF(irt.LT.1.OR.irt.GT.nnt)sw=1
      IF(irp.LT.1.OR.irp.GT.nnp)sw=1
      IF(sw.EQ.1)THEN
         WRITE(6,*)"Receiver lies outside model (rr,xr,zr)= ",irr,irt,irp,nnp,nnt,nnr
      ENDIF
      IF(irr.EQ.nnr)irr=irr-1
      IF(irt.EQ.nnt)irt=irt-1
      IF(irp.EQ.nnp)irp=irp-1

      nsts(irp,irt,irr)=1
   
      ri=gor-(i-1)*dnr+earth
      risti=ri*sin(cgot+(j-1)*dnt) 

      dTr=ttn(irp,irt,irr+1)-ttn(irp,irt,irr)
      dTr=dtr+ttn(irp-1,irt,irr+1)-ttn(irp-1,irt,irr)
      dTr=dtr+ttn(irp-1,irt-1,irr+1)-ttn(irp-1,irt-1,irr)
      dTr=dtr+ttn(irp,irt-1,irr+1)-ttn(irp,irt-1,irr)
      dTr=dtr/(4.0*dnr)

      nr=1

      par1=residuals(l)/cd(ishot,l)*(1-1./real(ntsrc(ishot)))
      
      if(ntsrc(ishot).eq.1)then
        write(*,*)'source has been recorded by one receiver',ishot
        stop
      endif
      
      par2=nr*dTr

      IF(abs(par1).GT.(1.0d-10).AND.abs(par2).GT.(1.0d-12))THEN
        IF(adjn.EQ.1)THEN
          adj(irp,irt,irr)=par1/par2
          adj(irp-1,irt,irr)=par1/par2
          adj(irp-1,irt-1,irr)=par1/par2
          adj(irp,irt-1,irr)=par1/par2
        ELSE
          adj(irp,irt,irr)=(1-1./real(ntsrc(ishot)))/cd(ishot,l)
          adj(irp-1,irt,irr)=(1-1./real(ntsrc(ishot)))/cd(ishot,l) 
          adj(irp-1,irt-1,irr)=(1-1./real(ntsrc(ishot)))/cd(ishot,l) 
          adj(irp,irt-1,irr)=(1-1./real(ntsrc(ishot)))/cd(ishot,l) 
        ENDIF
      ENDIF
      
    ENDIF
 ENDDO

END SUBROUTINE init_receiver

!*************************************************************
! This subroutine is dertermine the direction of sweeping.
!*************************************************************

SUBROUTINE init_sweep2
 USE globaladj
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

END SUBROUTINE init_sweep2

!**************************************************************
! This subroutine initializes the boundary conditions.
!**************************************************************

SUBROUTINE adjoint_sweep
 USE globaladj
 USE globalp
 INTEGER :: swp,i,j,k
 REAL*8 :: alphaw,alphae,betas,betan
 REAL*8 :: an,ap,bn,bp,cn,cp
 REAL*8 :: par1,par2,par3,par4,par5,ri,risti

 DO swp=1, 4
    DO i=imin(swp),imax(swp),istep(swp)
    DO j=jmin(swp),jmax(swp),jstep(swp)
    DO k=kmin(swp),kmax(swp),kstep(swp)  
       IF(nsts(k,j,i).EQ.0)THEN
         ri=gor-(i-1)*dnr+earth
         risti=ri*sin(cgot+(j-1)*dnt)    
 
         IF(max(ttn(k+1,j,i),ttn(k-1,j,i),ttn(k,j+1,i),ttn(k,j-1,i),ttn(k,j,i+1),ttn(k,j,i-1)).GE.big)THEN
           adj(k,j,i)=0
         ELSE
           an=-(ttn(k,j,i)-ttn(k,j,i-1))/(dnr)
           ap=-(ttn(k,j,i+1)-ttn(k,j,i))/(dnr)
           bn=-(ttn(k,j,i)-ttn(k,j-1,i))/(ri*dnt)
           bp=-(ttn(k,j+1,i)-ttn(k,j,i))/(ri*dnt)
           cn=-(ttn(k,j,i)-ttn(k-1,j,i))/(risti*dnp)
           cp=-(ttn(k+1,j,i)-ttn(k,j,i))/(risti*dnp)

           par1=(pos(an)*adj(k,j,i-1)-neg(ap)*adj(k,j,i+1))/(dnr)
           par2=(pos(bn)*adj(k,j-1,i)-neg(bp)*adj(k,j+1,i))/(ri*dnt)
           par3=(pos(cn)*adj(k-1,j,i)-neg(cp)*adj(k+1,j,i))/(risti*dnp)
           par4=(pos(ap)-neg(an))/(dnr)+(pos(bp)-neg(bn))/(ri*dnt)+(pos(cp)-neg(cn))/(risti*dnp)
          
           IF(abs(par4).GT.1.0d-12)THEN
             adj(k,j,i)=(par1+par2+par3)/par4
             !IF(abs(adj(k,j,i)).LT.1.0d-9)adj(k,j,i)=0
           ELSE
             adj(k,j,i)=0
           ENDIF

         ENDIF
       ENDIF      
    ENDDO
  ENDDO
  ENDDO
 ENDDO

END SUBROUTINE adjoint_sweep

!********************************************************
! Calculate the addition of an value and its absolute
!********************************************************

REAL*8 FUNCTION pos(par)
 REAL*8 :: par

 pos=(par+abs(par))/2.0

END FUNCTION pos

!*********************************************************
! Calculate the difference of an value and its absolute
!*********************************************************

REAL*8 FUNCTION neg(par)
 REAL*8 :: par

 neg=(par-abs(par))/2.0

END FUNCTION neg
