!**********************************************************************
! smooth the gradient
!**********************************************************************
subroutine gauss_smooth(g,g1)
 use globalp
 real*8:: g(1:nnp,1:nnt,1:nnr),g1(1:nnp,1:nnt,1:nnr)
 real*8, allocatable:: w(:,:,:)
 real*8:: sigmap,sigma2p,sigmat,sigma2t,sigmar,sigma2r
 real*8:: mywet,mysum
 integer:: i,j,k,i1,j1,k1,i2,j2,k2,l1,l2,l3
 
 if(nwp.gt.0)then
 
 allocate(w(-nwp:nwp,-nwt:nwt,-nwr:nwr))

 sigmap=dble(nwp)/3.
 sigma2p=sigmap*sigmap

 sigmat=dble(nwt)/3.
 sigma2t=sigmat*sigmat
 
 sigmar=dble(nwr)/3.
 sigma2r=sigmap*sigmar
 
 do i=-nwp, nwp
 do j=-nwt, nwt
 do k=-nwr, nwr
    w(i,j,k)=exp(-0.5d0*(dble(i*i)/sigma2p+dble(j*j)/sigma2t+dble(k*k)/sigma2r))
 enddo
 enddo
 enddo

 w=w/sum(w)
 g1=g

 do i=1, nnp
 do j=1, nnt
 do k=1, nnr
    mywet=0.0
    mysum=0.0
    do i1=-nwp, nwp
    do j1=-nwt, nwt
    do k1=-nwr, nwr
       l1=i+i1
       l2=j+j1
       l3=k+k1
       if(l1.lt.1)   l1=1-l1
       if(l1.gt.nnp) l1=2*nnp-l1+1
       if(l2.lt.1)   l2=1-l2
       if(l2.gt.nnt) l2=2*nnt-l2+1
       if(l3.lt.1)   l3=1-l3
       if(l3.gt.nnr) l3=2*nnr-l3+1
       mywet=mywet+w(i1,j1,k1)
       mysum=mysum+w(i1,j1,k1)*g(l1,l2,l3)
    enddo
    enddo
    enddo
    g1(i,j,k)=mysum/mywet
 enddo
 enddo
 enddo
 
 deallocate(w)
 else
 g1=g
 endif

end subroutine
