!**************************************************************
! This subroutine dertermines the direction of sweeping.
!**************************************************************

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



