!**********************************************************************
! calculate the gradient matrix for least-square cost function 
!**********************************************************************

SUBROUTINE subgradient(dobs,velr,vel,vel1D,grad,pgrad,precon,cost,cost1,cost2)
 ! input                                          | output
 ! dobs: observed data at recerver locations      | grad: 1D gradient matrix
 ! velr: reference velocity model                  | cost: cost function
 ! vel: current velocity model                    
 ! vel1D: 1D current velocity model
 
 USE globalp
 INTEGER:: i,j,k,r,l,ishot
 REAL:: cost,cost1,cost2
 REAL:: vel1D(nop),grad(nop),pgrad(nop),precon(nop)
 REAL*8:: dobs(nsrc,nrc),damp,gradtmax
 REAL*8:: velr(0:nvp+1,0:nvt+1,0:nvr+1),vel(0:nvp+1,0:nvt+1,0:nvr+1)
 REAL*8,ALLOCATABLE:: T(:,:,:,:),adj(:,:,:),padj(:,:,:),veln(:,:,:)
 REAL*8,ALLOCATABLE:: dcal(:,:),residuals(:,:),gradt(:,:,:),pgradt(:,:,:),hesst(:,:,:),gradtt(:,:,:,:),hesstt(:,:,:,:)
 REAL*8,ALLOCATABLE:: gradtv(:,:,:),pgradtv(:,:,:),preconv(:,:,:)
 REAL*8,ALLOCATABLE:: gradsmo(:,:,:),pgradsmo(:,:,:)
 
 ! vel1D = current 1D velocity model
 ! T = Traveltime fields for difference sources at refined grid
 ! adj = adjoint fields for every source at velocity grid
 ! padj = preconditioned fields for every source at velocity grid
 ! gradt, hesst = gradient and hess for every source at velocity grid
 ! gradtt, hesstt = gradient and hess for different sources at velocity grid
 ! dcal = calculated data for different source and receiver couples
 ! residuals = traveltime residuals for different source and receive couples
 ! cost = cost function
 ! damp = damping factor for preconditioner 

 ALLOCATE(T(nsrc,nnp,nnt,nnr),adj(nnp,nnt,nnr),padj(nnp,nnt,nnr))
 ALLOCATE(gradtt(nsrc,nnp,nnt,nnr),hesstt(nsrc,nnp,nnt,nnr),gradt(nnp,nnt,nnr),hesst(nnp,nnt,nnr),pgradt(nnp,nnt,nnr))
 ALLOCATE(dcal(nsrc,nrc),residuals(nsrc,nrc))
 ALLOCATE(veln(nnp,nnt,nnr))
 ALLOCATE(gradtv(nvp,nvt,nvr),pgradtv(nvp,nvt,nvr),preconv(nvp,nvt,nvr))
 ALLOCATE(gradsmo(nnp,nnt,nnr),pgradsmo(nnp,nnt,nnr))
 
 !********************************************** 
 ! Transform the 3D velocity into 1D velocity  
 !**********************************************
 
 CALL subv2m(vel,vel1D)     

 !************************************************************
 ! Traveltime field modeling for all sources using TDEE solver  
 !************************************************************

 call cpu_time(start)
 
 CALL subforward(vel,veln,T,dcal)         

 call cpu_time(finish)
 
 write(*,*)'("Time = ",f6.3," seconds.")',finish-start
 
 if(debug)then
 
 open(unit=888,file='traveltime.dat')
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    WRITE(888,*)gopd+(i-1)*dnpd,gotd-(j-1)*dntd,gor-(k-1)*dnr,T(1,i,j,k)
 ENDDO
 ENDDO
 ENDDO
 close(888)
 
 endif
 
 gradtt=0
 hesstt=0
 
 !*************************************************************
 ! Calculate gradients for different sources
 !*************************************************************

 DO ishot=1,  nsrc

    !*********************************************************
    ! Calculate traveltime residuals at receiver locations  
    !*********************************************************
    CALL subresiduals(ishot,dcal(ishot,:),dobs(ishot,:),residuals(ishot,:))    

    !**********************************************************
    ! Calculate adjoint-state field     
    !**********************************************************
   
    CALL subadjoint(ishot,T(ishot,:,:,:),residuals(ishot,:),adj,1)  

    !***********************************
    ! Calculate preconditioned adjoint    
    !***********************************

    if(preco.eq.1)then
      CALL subadjoint(ishot,T(ishot,:,:,:),residuals(ishot,:),padj,2)
    endif
    !**********************************************************
    ! Calculate gradient for mistfit function  
    !**********************************************************

    gradtt(ishot,:,:,:)=adj(:,:,:)/veln(:,:,:)**3
    hesstt(ishot,:,:,:)=padj(:,:,:)/veln(:,:,:)**3

 ENDDO
 
 if(debug)then
 OPEN(UNIT=10,FILE='Initial_vp.dat')
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    WRITE(10,*)gopd+(i-1)*dnpd,gotd-(j-1)*dntd,gor-(k-1)*dnr,veln(i,j,k)
 ENDDO
 ENDDO
 ENDDO
 CLOSE(10)
 endif

 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    gradt(i,j,k)=sum(gradtt(1:nsrc,i,j,k))
    hesst(i,j,k)=sum(hesstt(1:nsrc,i,j,k))
 ENDDO
 ENDDO
 ENDDO
 
 call gauss_smooth(gradt,gradsmo)
 gradtmax=maxval(abs(gradsmo))

 !************************************
 ! fix the crustal structure
 !************************************
 if(crt.eq.1)then
 do k=1, ncrt
 do j=1, nnt
 do i=1, nnp
    gradt(i,j,k)=0
 enddo
 enddo
 enddo
 endif
 
 if(debug)then
 open(unit=888,file='gradient.dat')
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    WRITE(888,*)gopd+(i-1)*dnpd,gotd-(j-1)*dntd,gor-(k-1)*dnr,gradsmo(i,j,k)/gradtmax
 ENDDO
 ENDDO
 ENDDO
 close(888)
 write(*,*) 'write out gradient'
 endif
 
 damp=maxval(abs(hesst(1:nnp,1:nnt,1:nnr)))*damping
 
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    if(preco.eq.1)then
      hesst(i,j,k)=damp+hesst(i,j,k)
      pgradt(i,j,k)=gradt(i,j,k)/hesst(i,j,k)
    else
      hesst(i,j,k)=1
      pgradt(i,j,k)=gradt(i,j,k)
    endif
 ENDDO
 ENDDO
 ENDDO

 call gauss_smooth(pgradt,pgradsmo)
 pgradt=pgradsmo
 gradtmax=maxval(abs(pgradt))

 if(debug)then
 
 open(unit=888,file='pregradient.dat')
 DO i=1, nnp
 DO j=1, nnt
 DO k=1, nnr
    WRITE(888,*)gopd+(i-1)*dnpd,gotd-(j-1)*dntd,gor-(k-1)*dnr,pgradt(i,j,k)/gradtmax
 ENDDO
 ENDDO
 ENDDO
 close(888)
 write(*,*) 'write out pregradient'
 endif

 DO i=1, nvp
 DO j=1, nvt
 DO k=1, nvr
    gradtv(i,j,k)=gradt(nrfp*(i-1)+1,nrft*(j-1)+1,nrfr*(k-1)+1)+epsilon*(vel(i,j,k)-velr(i,j,k))/cm(i,j,k)
    pgradtv(i,j,k)=pgradt(nrfp*(i-1)+1,nrft*(j-1)+1,nrfr*(k-1)+1)+epsilon*(vel(i,j,k)-velr(i,j,k))/cm(i,j,k)
    preconv(i,j,k)=hesst(nrfp*(i-1)+1,nrft*(j-1)+1,nrfr*(k-1)+1)
 ENDDO
 ENDDO
 ENDDO
 
 r=0
 pgrad=0
 grad=0
 precon=0
 DO i=0, nvp+1
 DO j=0, nvt+1
 DO k=0, nvr+1
    r=r+1
    IF((i.EQ.0).OR.(j.EQ.0).OR.(k.EQ.0).OR.(i.EQ.(nvp+1)).OR.(j.EQ.(nvt+1)).OR.(k.EQ.(nvr+1)))THEN
      pgrad(r)=0
      grad(r)=0
      precon(r)=damp
    ELSE
      pgrad(r)=pgrad(r)+pgradtv(i,j,k)
      grad(r)=grad(r)+gradtv(i,j,k)
      precon(r)=precon(r)+preconv(i,j,k)
    ENDIF
 ENDDO
 ENDDO
 ENDDO
 
 ! normalize preconditioned gradient
 do r=1, nop
    pgrad(r)=pgrad(r) !*sum(abs(grad(1:nop)))/sum(abs(pgrad(1:nop)))
 enddo
 
 ! normalize preconditioner
 !do r=1, nop
    !precon(r)=1 !precon(r)/maxval(abs(precon(1:nop)))
 !enddo

 !******************************************
 ! calculate cost function for all sources
 !******************************************
 CALL submisfit(residuals,vel,velr,cost,cost1,cost2)    
 
 DEALLOCATE(adj,T,padj,gradt,gradtt,hesst,hesstt,pgradt)
 DEALLOCATE(dcal,residuals)
 DEALLOCATE(veln)
 DEALLOCATE(gradtv,pgradtv)
 DEALLOCATE(gradsmo,pgradsmo)

END SUBROUTINE subgradient


