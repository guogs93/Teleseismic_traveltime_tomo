!*******************************************************************************************************
! PROGRAM TEAT3d
! Gaoshan Guo, CAS IGGCAS and CNRS GEOAZUR, Mar, 2021
! 3D teleseismic eikonal + adjoint tomograpahy (TEAT).
! A lock sweeping eikonal solver is implemented for the forward problem.
! Adjoint-state method is used to obtain gradient and preconditioner for least-square misfit function
!*******************************************************************************************************

!*********************************************************
! Global parameters for forward, adjoint and inversion.
!*********************************************************

MODULE globalp

 REAL*8:: dvr,dvtd,dvpd,dvt,dvp
 REAL*8:: dnr,dntd,dnpd,dnt,dnp
 REAL*8:: gor,gotd,gopd,cgot,gop
 REAL*8:: epsilon
 REAL*8,PARAMETER:: earth=6371.0
 INTEGER:: nrfr,nrft,nrfp,rmtr,fom,maxiter,preco
 INTEGER:: nvr,nvt,nvp,nnr,nnt,nnp,nsrc,nrc,nop,nwp,nwt,nwr
 INTEGER,ALLOCATABLE:: ntsrc(:)
 REAL*8,ALLOCATABLE:: srcm(:)
 REAL*8,ALLOCATABLE:: rcr(:),rct(:),rcp(:)
 INTEGER,ALLOCATABLE:: tstat(:,:)
 REAL*8,ALLOCATABLE:: dref(:,:),cd(:,:),cm(:,:,:)
 CHARACTER(LEN=20):: itimes
 REAL*8:: damping
 LOGICAL,PARAMETER:: debug=.false.
 REAL :: start, finish
 INTEGER:: crt, ncrt, mute, mode
 CHARACTER(LEN=160):: crtmodel
 INTEGER:: rseed, agn
 REAL*8:: sdgn, bound


 !***********************************************************
 ! nvr, nvt, nvp = Number of velocity nodes in radius, theta (latitude), phi (longitude)
 ! gor, gotd, gopd = Origin of velocity grid in radius, theta, phi (unit=degrees)
 ! dvr, dvtd, dvpd = Node separation of velocity grid in radius, theta, phi (unit=degrees)
 ! nrfr,nrfx, nrfz = B-spline dicing level in radius, theta, phi
 ! nnr, nnt, nnp = Number of nodes of refined grid in radius, theta and phi
 ! dnr, dntd, dnpd = Number of nodes of refined grid in radius, theta and phi
 ! cgot, gop = Origin of velocity grid in co-theta, phi (unit=radians)
 ! dvt, dvp = Node separation of velocity grid in radius, theta, phi (unit=radians)
 ! dnt, dnp = Node separation of refined velocity grid in theta and phi (unit=radians)
 ! epsilon = damping term
 ! fom = accuracy for FD stencil in eikonal solver
 ! nw = the length of Gaussian smoothing filter 
 ! tstat = status of observed time pick (0=none, 1=picked)
 ! rmtr = remove model traveltime residual (0=no,1=yes)
 ! cd = data covariance matrix
 ! cm = model covariance matrix
 ! maxiter = Maxminum iterations for eikonal and adjoint solver
 ! itimes = file containing initial traveltimes
 ! nrfr,nrft,nrfp = refined number in r, theta, phi
 !***********************************************************

END MODULE

!*****************************************************************
! The main program for 3D teleseismic eikonal adjoint tomography
!*****************************************************************

PROGRAM TEAT3D
 USE globalp
 INCLUDE 'optim_type.h'
 TYPE(optim_type):: optim
 CHARACTER(LEN=4):: FLAG
 CHARACTER(LEN=160):: gridi,gridr,otimes,rtimes,receivers,sources,gridm,namei
 CHARACTER(LEN=9):: pht
 INTEGER:: i,j,k,r,igrad,niter,optimalgo,maxinv
 REAL*8:: perc,alpha
 REAL*8,ALLOCATABLE:: scp(:),sct(:),scr(:)
 REAL*8,ALLOCATABLE:: vel0(:,:,:),vel(:,:,:),dobs(:,:),velr(:,:,:)
 REAL:: cost,cost1,cost2
 REAL,ALLOCATABLE:: vel1D(:),grad(:),pgrad(:),precon(:)
 INTEGER,PARAMETER:: lbytes=4,ibound=0
 REAL,PARAMETER:: pi=3.1415926535898,thresconv=1e-2,ngradlbfgs=5
 REAL*8,ALLOCATABLE:: velr1D(:)
  
 !*************************************************************
 ! gridi = file containing starting velocity field
 ! otimes = file containing observed traveltime residuals
 ! rtimes = file containing reference traveltimes
 ! rmtr = remove model traveltime residual (0=no,1=yes)
 ! sources = file containing source locations
 ! receivers = file containing receiver locations
 ! maxinv = maxminum iretations for inversion
 !*************************************************************
 ! parameter for seiscope optimization
 ! ibound = use bound resctriction or not
 ! thresconv = convergence threshold for inversion
 ! ngradlbfgs = number of gradient for L-BFGS inversion
 !*************************************************************

 !*************************************************************
 ! Read the input parameters
 !*************************************************************

 OPEN(UNIT=10,FILE="input.in",STATUS='old')
 READ(10,*)
 READ(10,*)
 READ(10,*)
 READ(10,*)mode
 READ(10,1)gridi
 READ(10,1)gridr
 READ(10,1)otimes
 READ(10,1)rtimes
 READ(10,*)itimes
 READ(10,*)rmtr
 READ(10,1)sources
 READ(10,1)receivers
 READ(10,*)maxiter
 READ(10,*)maxinv
 READ(10,*)epsilon
 READ(10,*)nrfr,nrft,nrfp
 READ(10,*)fom
 READ(10,*)nwp,nwt,nwr
 READ(10,*)optimalgo
 READ(10,*)perc,damping
 READ(10,*)preco
 READ(10,*)crt, crtmodel, ncrt
 READ(10,*)mute
 READ(10,*)agn, sdgn, rseed
 READ(10,*)bound
 CLOSE(10)
1 FORMAT(a25)
 CLOSE(10)
  
 !***********************************************
 ! Read reference = initial and reference grids
 !***********************************************
 OPEN(UNIT=11,FILE=gridi,status='old')
 READ(11,*)nvr,nvt,nvp
 READ(11,*)gor,gotd,gopd
 READ(11,*)dvr,dvtd,dvpd
 
 cgot=(90.0-gotd)*pi/180.0
 gop=gopd*pi/180.0

 dvt=dvtd*pi/180.0
 dvp=dvpd*pi/180.0

 !**************************************************************************************
 ! Calculate total numer of refined nodes in r, theta, phi and the refined grid spacing.
 !**************************************************************************************

 nnr=(nvr-1)*nrfr+1
 nnt=(nvt-1)*nrft+1
 nnp=(nvp-1)*nrfp+1
 dnr=dvr/nrfr
 dntd=dvtd/nrft
 dnpd=dvpd/nrfp
 
 write(*,*)'forward grid'
 write(*,*)nnr,nnt,nnp
 write(*,*)'refine factor'
 write(*,*)nrfr,nrft,nrfp
 write(*,*)'inverse grid'
 write(*,*)nvr,nvt,nvp

 !******************************************************
 ! convert space of refined grids in degrees to radians
 !******************************************************

 dnt=dntd*pi/180.0
 dnp=dnpd*pi/180.0

 !********************************************
 ! allocate memory to model arrays
 !********************************************
 
 ALLOCATE(vel0(0:nvp+1,0:nvt+1,0:nvr+1),velr(0:nvp+1,0:nvt+1,0:nvr+1),vel(0:nvp+1,0:nvt+1,0:nvr+1),cm(0:nvp+1,0:nvt+1,0:nvr+1))
 
 DO i=0, nvp+1
 DO j=0, nvt+1
 DO k=0, nvr+1
    READ(11,*)vel0(i,j,k),cm(i,j,k)
 ENDDO
 ENDDO
 ENDDO
 CLOSE(11)
 
 OPEN(UNIT=11,FILE=gridr,status='old')
 READ(11,*)
 READ(11,*)
 READ(11,*)
 DO i=0, nvp+1
 DO j=0, nvt+1
 DO k=0, nvr+1
    READ(11,*)velr(i,j,k),cm(i,j,k)
 ENDDO
 ENDDO
 ENDDO
 CLOSE(11)
 
 !***************************************
 ! number of points
 !***************************************
 
 nop=(nvr+2)*(nvp+2)*(nvt+2)
 WRITE(*,*) 'Number of points = ',nop
 
 ALLOCATE(vel1D(nop))
 ALLOCATE(grad(nop))
 ALLOCATE(pgrad(nop))
 ALLOCATE(precon(nop))
 ALLOCATE(velr1D(nop))
 
 CALL subv2m(velr,velr1D) 
 
 !*****************************************************
 ! read acquisition data
 !***************************************************
 
 OPEN(UNIT=60,FILE=sources,STATUS='old')
 OPEN(UNIT=70,FILE=receivers,STATUS='old')
 READ(60,*)
 READ(60,*)nsrc
 ALLOCATE(scr(nsrc),sct(nsrc),scp(nsrc))

 DO i=1, nsrc
    READ(60,*)sct(i),scp(i),scr(i)
    READ(60,'(a8)')pht
      
    !*********************************************************
    ! Convert source coordinates in degrees to radians
    !*********************************************************
      
    sct(i)=(90.0-sct(i))*pi/180.0
    scp(i)=scp(i)*pi/180.0
 ENDDO
 READ(70,*)
 READ(70,*)nrc
 ALLOCATE(rcr(nrc),rct(nrc),rcp(nrc))

 DO i=1, nrc
    READ(70,*)rcr(i),rct(i),rcp(i)
     
    !********************************************************
    ! Convert receiver coordinates in degrees to radians
    !********************************************************
     
    rct(i)=(90.0-rct(i))*pi/180.0
    rcp(i)=rcp(i)*pi/180.0
 ENDDO
 CLOSE(60)
 CLOSE(70)
 
 WRITE(*,*) 'Number of sources and receivers', nsrc, nrc

 !********************************************************************************
 ! Now read in the observed (dobs) and reference (dref) traveltime residuals
 !********************************************************************************
 ALLOCATE(dobs(nsrc,nrc),cd(nsrc,nrc),dref(nsrc,nrc),tstat(nsrc,nrc),srcm(nsrc),ntsrc(nsrc))

 OPEN(UNIT=10,FILE=otimes,STATUS='old')
 DO i=1, nsrc
 DO j=1, nrc
    READ(10,*)tstat(i,j),dobs(i,j),cd(i,j)
 ENDDO
 ENDDO
 CLOSE(10)
 
 if(mode.ne.0)then
   OPEN(UNIT=30,FILE=rtimes,STATUS='old')
   DO i=1, nsrc
   DO j=1, nrc
      READ(30,*)dref(i,j)
   ENDDO
   enddo
 endif

 CLOSE(10)
 CLOSE(30)

 open(unit=20,file='variance.dat')
 write(20,*)'data',' ','model'
 !****************************************************
 ! Calculate gradient for initial velocity
 !****************************************************
 
 CALL subgradient(dobs,velr,vel0,vel1D,grad,pgrad,precon,cost,cost1,cost2)
 
 write(20,*)cost1,cost2
 
 !*******************************************************
 ! Scale the gradient
 !*******************************************************
 
 CALL scalegradient(pgrad,cost,vel1D,perc,alpha)
 
 grad(:)=grad(:)*alpha
 pgrad(:)=pgrad(:)*alpha
 cost=cost*alpha
 
 WRITE(*,*)'scaled cost=',cost
 
 !*********************************************************
 ! Configuration for optimization
 !*********************************************************
 
 optim%niter_max = maxinv
 optim%conv = thresconv
 optim%print_flag = 1
 optim%bound = ibound
 optim%l = ngradlbfgs
 optim%debug = .false. 
 
 ! apply the bound contraint

 ALLOCATE(optim%ub(nop),optim%lb(nop))
 
 !*********************************************************
 ! Optimization loop
 !*********************************************************
 
 FLAG='INIT'       
 niter = 0
 igrad=0
 
 DO WHILE ( (FLAG .NE. 'CONV') .AND. (FLAG .NE. 'FAIL'))  

       IF (optimalgo==0) THEN
          WRITE(*,*) 'Steepest-descent optimization'
          
          CALL PSTD(nop,vel1D,cost,grad,pgrad,optim,FLAG)
          
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
          WRITE(*,*) 'FLAG = ',FLAG
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
       ELSEIF (optimalgo==1) THEN
          WRITE(*,*) 'LBFGS optimization'
          
          CALL PLBFGS(nop,vel1D,cost,grad,pgrad,optim,FLAG)
          
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
          WRITE(*,*) 'FLAG = ',FLAG
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
          
       ELSEIF (optimalgo==2) THEN
          WRITE(*,*) 'CG optimization'
          
          CALL PNLCG(nop,vel1D,cost,grad,pgrad,optim,FLAG)
          
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
          WRITE(*,*) 'FLAG = ',FLAG
          WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
       ENDIF
         
       IF (FLAG.EQ.'GRAD') THEN
          WRITE(*,*) 'COMPUTE GRAD'
          igrad = igrad + 1
 
          CALL subm2v(vel1D,vel)
          
          CALL subgradient(dobs,velr,vel,vel1D,grad,pgrad,precon,cost,cost1,cost2)
          
          grad(:)=grad(:)*alpha
          cost=cost*alpha
          pgrad(:)=grad(:)*alpha
          
          WRITE(*,*)'scaled cost=',cost

       !ELSE IF (FLAG.EQ.'PREC') THEN  
          
       !   WRITE(*,*) "PRECONDITIONER REQUIRED BY P-LBFGS 1"
       !   optim%q_plb(:)=optim%q_plb(:)/precon(:)
          
       ELSEIF (FLAG.EQ.'NSTE' .OR. FLAG.EQ.'CONV') THEN
          niter =  niter + 1 
          
          CALL subm2v(vel1D,vel)
          
          WRITE(*,*) "========================================"
          WRITE(*,*) 'WRITE FTDET MODEL',niter
          WRITE(*,*) "========================================"
       
          DO i=0, nvp+1
          DO j=0, nvt+1
          DO k=0, nvr+1
             if(vel(i,j,k).ge.(velr(i,j,k)*(1+bound)))then
               vel(i,j,k)=velr(i,j,k)*(1+bound)
             endif
             if(vel(i,j,k).le.(velr(i,j,k)*(1-bound)))then
               vel(i,j,k)=velr(i,j,k)*(1-bound)
             endif
          ENDDO
          ENDDO
          ENDDO
       
          !******************************************************************************
          ! 
          !******************************************************************************
       
          call subname(niter,namei)
       
          OPEN(UNIT=10,FILE='gridc_'//namei(1:LEN_TRIM(namei))//'.vtx',STATUS='unknown')
          WRITE(10,*)nvr,nvt,nvp
          WRITE(10,'(3f12.6)')gor,gotd,gopd
          WRITE(10,'(3f12.6)')dvr,dvtd,dvpd
          WRITE(10,'(1X)')
          DO i=0, nvp+1
          DO j=0, nvt+1
          DO k=0, nvr+1
             WRITE(10,'(2f12.8)')vel(i,j,k),cm(i,j,k)
          ENDDO
          WRITE(10,'(1X)')
          ENDDO
          WRITE(10,'(1X)')
          ENDDO
          CLOSE(10)
       
          write(20,*)cost1,cost2
       
       ENDIF

 ENDDO

 close(20)
 
 DEALLOCATE(vel0,vel,velr,cm)
 DEALLOCATE(vel1D,grad,pgrad,precon)
 DEALLOCATE(scr,sct,scp,rcr,rct,rcp)
 DEALLOCATE(dobs,tstat,dref,cd)
 DEALLOCATE(srcm,ntsrc)
 DEALLOCATE(optim%ub,optim%lb)

END PROGRAM TEAT3D         
    
