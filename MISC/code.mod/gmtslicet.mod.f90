!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to take depth slices through
! the FMM computed traveltime field and put them in a
! form suitable for contour plotting by GMT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i10) :: earthr
REAL(KIND=i10), PARAMETER :: pi=3.141592653589793
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: vela
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: veln
! vela = diced velocity values
! inta = diced interface values
! veln = velocity grid values
! intn = interface grid values
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! earthr = earth radius
END MODULE globalp

PROGRAM slice
USE globalp
IMPLICIT NONE
INTEGER :: i,j,k,l,m,n,i1,j1,k1,extrds,extrnss,extrews,extrgcs,idp,iew,ins
INTEGER :: checkstat,isw,conp,cont,conr,id1,id2
INTEGER :: nnr,nnt,nnp,pvgd,pvgns,pvgew,pvggc,pvgcs,prpd,prpns,prpew,nre,nr
INTEGER :: idm1,idm2,idm3,nnx,nnz,stp,stt,str
INTEGER :: ptfd,ptfns,ptfew,rnode,ngch,ngcd
INTEGER :: ddt,ddp,dnsr,dnst,dewr,dewp
INTEGER :: intid,indt,indp,ios,rnodeth,rnodeph,rnoder
REAL(KIND=i5) :: rdep,rlat,rlon
REAL(KIND=i10) :: gor,got,gop,rgsr,rgst,rgsp,sldep,slns,slew
REAL(KIND=i10) :: tt,ttt,ttb,sumi,sumj,sumk
REAL(KIND=i10) :: lft,rgt,btm,top
REAL(KIND=i10) :: rdm,rd1,rd2,rd3,u,v,w
REAL(KIND=i10), DIMENSION(:) :: wi(4)
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: ttn
REAL(KIND=i10), PARAMETER :: dtol=1.0e-6
REAL(KIND=i10), PARAMETER :: stol=1.0e-8
REAL(KIND=i10) :: deltas,deltaskm,azim
REAL(KIND=i10) :: deltanp,thgcp,phgcp,goti,gopi
REAL(KIND=i10) :: slgclat1,slgclon1,slgclat2,slgclon2
REAL(KIND=i10) :: slgclat1d,slgclon1d,slgclat2d,slgclon2d
CHARACTER (LEN=30) :: ifilet,ifilev,irfile,ifilevr,sep
CHARACTER (LEN=30) :: ofiledb,ofilensb,ofileewb
CHARACTER (LEN=30) :: ofiledv,ofilensv,ofileewv
CHARACTER (LEN=30) :: ofiledt,ofilenst,ofileewt
CHARACTER (LEN=30) :: ofiledr,ofilensr,ofileewr
CHARACTER (LEN=30) :: ofilergc,ofilesgc,ofilercgc,ofilegcb
CHARACTER (LEN=30) :: ofilegcv
!
! sldep = slice depth
! slns = NS slice
! slew = EW slice
! ifilet = input 3-D traveltime grid file
! ifilev = input velocity grid file
! ifilevr = input reference velocity grid file
! irfile = input ray path file
! ofiledt = output depth slice file for traveltimes
! ofilenst = output N-S slice file for traveltimes
! ofileewt = output E-W slice file for traveltimes
! ofiledb = bounds for output depth slice file
! ofilensb = bounds for output N-S slice file
! ofileewb = bounds for output E-W slice file
! ofilegcb = bounds for output great circle file
! ofiledv = output velocity file for depth slice
! ofilensv = output velocity file for N-S slice
! ofileewv = output velocity file for E-W slice
! ofiledr = output ray path file in depth
! ofilensr = output ray path file in N-S
! ofileewr = output ray path file in E-W
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! rgsr,rgst,rgsp = Refined node spacing in r,theta,phi
! idp = radial index of slice depth
! ins = NS index for latitude slice
! iew = EW index for longitude slice
! tt = traveltime at or near a node
! ttn = traveltime grid values
! ttt = traveltime of node above slice 
! ttb =  traveltime of node below slice
! lft,rgt,btm,top = plotting bounds
! pvgd = plot velocity depth slice (0=no,1=yes)
! pvgns = plot velocity N-S slice (0=no, 1=yes)
! pvgew = plot velocity E-W slice (0=no, 1=yes)
! ptfd = plot traveltime field depth slice? (0=no,1=yes)
! ptfns = plot traveltime field N-S slice? (0=no, 1=yes)
! ptfew = plot traveltime field E-W slice? (0=no,1=yes)
! veln = velocity grid values
! prpd = plot ray path in depth? (0=no,1=yes)
! prpns = plot ray path in N-S? (0=no,1=yes)
! prpew = plot ray path in E-W? (0=no,1=yes)
! nre = number of ray elements
! rdep,rlat,rlon = ray point depth, latitude, longitude
! nr = number of receivers
! sep = character marker for separating rays
! extrds = extract depth slice? (0=no,1=yes)
! extrnss = extract N-S slice? (0=no,1=yes)
! extrews = extract E-W slice? (0=no,1=yes)
! extrgcs = extract great circle slice? (0=no, 1=yes)
! ddt,ddp = dicing of velocity grid for depth slice
! dnsr,dnst = dicing of velocity grid for N-S slice
! dewr,dewp = dicing of velocity grid for E-W slice
! u,v,w = bspline independent parameters
! ui,vi = bspline basis functions
! vela = diced velocity values
! nnx,nnz = dimensions of vela
! conr,conp,cont = variables for edge of bspline grid
! str,stp,stt = counters for vela grid points
! rnode = reference node for slice
! slgclat1,slgclon1 = First point for great circle slice
! slgclat2,slgclon2 = Second point for great circle slice
! slgclat1,slgclon1 = First point for great circle slice (in degrees)
! slgclat2,slgclon2 = Second point for great circle slice (in degrees)
! ngch,ngcd = Number of sample points in horizontal and depth 
!             for great circle slice.
! stol = Slice tolerance for great circle slice
! azim = Azimuth of great circle slice
! deltanp = Angular distance to point along great circle
! thgcp,phgcp = Lat and long of point along great circle
! rnodeth,rnodeph,rnoder = Reference grid nodes for great circle slice
!
OPEN(UNIT=10,FILE='gmtslicet.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a26)')ifilev
READ(10,'(a26)')ifilevr
READ(10,'(a26)')ifilet
READ(10,'(a26)')irfile
READ(10,*) earthr
!
! Read in slice parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)extrds
READ(10,*)sldep
READ(10,'(a22)')ofiledb
READ(10,*)extrnss
READ(10,*)slns
READ(10,'(a22)')ofilensb
READ(10,*)extrews
READ(10,*)slew
READ(10,'(a22)')ofileewb
READ(10,*)extrgcs
READ(10,*)slgclat1d,slgclon1d
READ(10,*)slgclat2d,slgclon2d
READ(10,'(a27)')ofilegcb
!
! Calculate GMT bounds files. Start off by reading in
! velocity grid.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)nnr,nnt,nnp
READ(20,*)gor,got,gop
READ(20,*)rgsr,rgst,rgsp
ALLOCATE(veln(0:nnr+1,0:nnt+1,0:nnp+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
ENDIF
DO i=0,nnp+1
   DO j=0,nnt+1
      DO k=0,nnr+1
         READ(20,*)veln(k,j,i)
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
!
! Now read in velocity grid parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pvgd
READ(10,*)ddt,ddp
READ(10,'(a22)')ofiledv
READ(10,*)pvgns
READ(10,*)dnsr,dnst
READ(10,'(a22)')ofilensv
READ(10,*)pvgew
READ(10,*)dewr,dewp
READ(10,'(a22)')ofileewv
READ(10,*)ngch,ngcd
READ(10,'(a22)')ofilegcv
!
! Calculate GMT bounds file for depth slice if required.
!
IF(extrds.EQ.1)THEN
   lft=gop
   rgt=gop+(nnp-1)*rgsp
   btm=got-(nnt-1)*rgst
   top=got
   OPEN(UNIT=50,FILE=ofiledb,STATUS='unknown')
   WRITE(50,'(f16.11)')lft
   WRITE(50,'(f16.11)')rgt
   WRITE(50,'(f16.11)')btm
   WRITE(50,'(f16.11)')top
   id1=(nnp-1)*ddp+1
   id2=(nnt-1)*ddt+1
   WRITE(50,*)id1
   WRITE(50,*)id2
ENDIF
!
! Calculate GMT bounds file for N-S slice if required.
!
IF(extrnss.EQ.1)THEN
   rgt=got-(nnt-1)*rgst
   lft=got
   btm=gor-(nnr-1)*rgsr
   top=gor
   OPEN(UNIT=60,FILE=ofilensb,STATUS='unknown')
   WRITE(60,'(f16.11)')rgt
   WRITE(60,'(f16.11)')lft
   WRITE(60,'(f16.11)')btm
   WRITE(60,'(f16.11)')top
   id1=(nnt-1)*dnst+1
   id2=(nnr-1)*dnsr+1
   WRITE(60,*)id1
   WRITE(60,*)id2
ENDIF
!
! Calculate GMT bounds file for E-W slice if required.
!
IF(extrews.EQ.1)THEN
   lft=gop
   rgt=gop+(nnp-1)*rgsp
   btm=gor-(nnr-1)*rgsr
   top=gor
   OPEN(UNIT=70,FILE=ofileewb,STATUS='unknown')
   WRITE(70,'(f16.11)')lft
   WRITE(70,'(f16.11)')rgt
   WRITE(70,'(f16.11)')btm
   WRITE(70,'(f16.11)')top
   id1=(nnp-1)*dewp+1
   id2=(nnr-1)*dewr+1
   WRITE(70,*)id1
   WRITE(70,*)id2
ENDIF
!
! Calculate GMT bounds file for great circle slice if required.
!
IF(extrgcs.EQ.1)THEN
   slgclat1=slgclat1d*pi/180.0
   slgclat2=slgclat2d*pi/180.0
   slgclon1=slgclon1d*pi/180.0
   slgclon2=slgclon2d*pi/180.0
   IF(slgclon1d.EQ.slgclon2d)slgclon2=slgclon2+dtol
!
!  Compute angular distance between end points of slice
!
   deltas=SIN(slgclat2)*SIN(slgclat1)
   deltas=deltas+COS(slgclat2)*COS(slgclat1)*COS(slgclon2-slgclon1)
   deltas=ACOS(deltas)
   deltaskm=deltas*earthr
   lft=0.0
   rgt=deltaskm
   btm=gor-(nnr-1)*rgsr
   top=gor
   OPEN(UNIT=70,FILE=ofilegcb,STATUS='unknown')
   WRITE(70,'(f16.11)')lft
   WRITE(70,'(f16.11)')rgt
   WRITE(70,'(f16.11)')btm
   WRITE(70,'(f16.11)')top
   WRITE(70,'(f16.11)')ngch
   WRITE(70,'(f16.11)')ngcd
   CLOSE(70)
ENDIF
!
! Read in reference velocity grid file if required
!
IF(pvgd.EQ.1.OR.pvgns.EQ.1.OR.pvgew.EQ.1.OR.pvggc.EQ.1)THEN
   isw=0
   OPEN(UNIT=20,FILE=ifilevr,status='old')
   READ(20,*)idm1,idm2,idm3
   IF(idm1.NE.nnr.OR.idm2.NE.nnt.OR.idm3.NE.nnp)isw=1
   READ(20,*)rd1,rd2,rd3
   IF(rd1.NE.gor.OR.rd2.NE.got.OR.rd3.NE.gop)isw=1
   READ(20,*)rd1,rd2,rd3
   IF(rd1.NE.rgsr.OR.rd2.NE.rgst.OR.rd3.NE.rgsp)isw=1
   IF(isw.EQ.1)THEN
      WRITE(6,*)'ERROR! Actual velocity grid and reference'
      WRITE(6,*)'velocity grid have different dimensions or'
      WRITE(6,*)'different numbers of grid points!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   DO i=0,nnp+1
      DO j=0,nnt+1
         DO k=0,nnr+1
            READ(20,*)rd1
!
!           This gives the velocity anomaly.
!
            veln(k,j,i)=veln(k,j,i)-rd1
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Extract depth slice if required
!
IF(pvgd.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(sldep.GT.gor.OR.sldep.LT.gor-rgsr*(nnr-1))THEN
     WRITE(6,*)'Requested depth slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(nnt-1)*ddt+1
  nnz=(nnp-1)*ddp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(ddt+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,ddt+1
     u=ddt
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(ddp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,ddp+1
     u=ddp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given depth
!
  rnode=INT((gor-sldep)/rgsr)+1
  u=ABS(gor-sldep-(rnode-1)*rgsr)/rgsr
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,nnp-1
     conp=ddp
     IF(i==nnp-1)conp=ddp+1
     DO j=1,nnt-1
        cont=ddt
        IF(j==nnt-1)cont=ddt+1
        DO l=1,conp
           stp=ddp*(i-1)+l
           DO m=1,cont
              stt=ddt*(j-1)+m
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=wi(k1)*veln(rnode-2+k1,j-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+ui(m,j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(stt,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofiledv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
!
! Extract N-S slice if required
!
IF(pvgns.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slns.LT.gop.OR.slns.GT.gop+rgsp*(nnp-1))THEN
     WRITE(6,*)'Requested N-S slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(nnr-1)*dnsr+1
  nnz=(nnt-1)*dnst+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(dnsr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dnsr+1
     u=dnsr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dnst+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dnst+1
     u=dnst
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given longitude
!
  rnode=INT((slns-gop)/rgsp)+1
  u=ABS(slns-gop-(rnode-1)*rgsp)/rgsp
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO j=1,nnt-1
     cont=dnst
     IF(j==nnt-1)cont=dnst+1
     DO k=1,nnr-1
        conr=dnsr
        IF(k==nnr-1)conr=dnsr+1
        DO m=1,cont
           stt=dnst*(j-1)+m
           DO n=1,conr
              str=dnsr*(k-1)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*veln(k-2+k1,j-2+j1,rnode-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+vi(m,j1)*sumk
                 ENDDO
                 sumi=sumi+wi(i1)*sumj
              ENDDO
              vela(str,stt)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofilensv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
!
! Extract E-W slice if required
!
IF(pvgew.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slew.GT.got.OR.slew.LT.got-rgst*(nnt-1))THEN
     WRITE(6,*)'Requested E-W slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(nnr-1)*dewr+1
  nnz=(nnp-1)*dewp+1
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Compute the values of the basis functions
!
  ALLOCATE(ui(dewr+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
  ENDIF
  DO i=1,dewr+1
     u=dewr
     u=(i-1)/u
     ui(i,1)=(1.0-u)**3/6.0
     ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     ui(i,4)=u**3/6.0
  ENDDO
  ALLOCATE(vi(dewp+1,4), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
  ENDIF
  DO i=1,dewp+1
     u=dewp
     u=(i-1)/u
     vi(i,1)=(1.0-u)**3/6.0
     vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
     vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
     vi(i,4)=u**3/6.0
  ENDDO
!
! Work out the value of u for the given latitude
!
  rnode=INT((got-slew)/rgst)+1
  u=ABS(got-slew-(rnode-1)*rgst)/rgst
  wi(1)=(1.0-u)**3/6.0
  wi(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
  wi(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
  wi(4)=u**3/6.0
  DO i=1,nnp-1
     conp=dewp
     IF(i==nnp-1)conp=dewp+1
     DO k=1,nnr-1
        conr=dewr
        IF(k==nnr-1)conr=dewr+1
        DO l=1,conp
           stp=dewp*(i-1)+l
           DO n=1,conr
              str=dewr*(k-1)+n
              sumi=0.0
              DO i1=1,4
                 sumj=0.0
                 DO j1=1,4
                    sumk=0.0
                    DO k1=1,4
                       rdm=ui(n,k1)*veln(k-2+k1,rnode-2+j1,i-2+i1)
                       sumk=sumk+rdm
                    ENDDO
                    sumj=sumj+wi(j1)*sumk
                 ENDDO
                 sumi=sumi+vi(l,i1)*sumj
              ENDDO
              vela(str,stp)=sumi
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  OPEN(UNIT=30,FILE=ofileewv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,ui,vi, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,ui,vi'
  ENDIF
ENDIF
!
! Extract E-W slice if required
!
IF(pvgcs.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slgclat1d.GT.got+stol.OR.slgclat1d.LT.got-rgst*(nnt-1)-stol)THEN
     WRITE(6,*)'Requested latitude of point 1 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclat2d.GT.got+stol.OR.slgclat2d.LT.got-rgst*(nnt-1)-stol)THEN
     WRITE(6,*)'Requested latitude of point 2 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclon1d.LT.gop-stol.OR.slgclon1d.GT.gop+rgsp*(nnp-1)+stol)THEN
     WRITE(6,*)'Requested longitude of point 1 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclon2d.LT.gop-stol.OR.slgclon2d.GT.gop+rgsp*(nnp-1)+stol)THEN
     WRITE(6,*)'Requested longitude of point 2 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Compute azimuth of great circle cross-section
!
  azim=COS(slgclat1)*SIN(slgclat2)-COS(slgclat2)*SIN(slgclat1)*COS(slgclon2-slgclon1)
  azim=ACOS(azim/SIN(deltas))
!
! Allocate memory to velocity grid array
!
  nnx=ngcd
  nnz=ngch
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Now call a bspline subroutine to compute great circle slice
!
  DO j=1,nnz
     deltanp=deltas*(j-1)/REAL(nnz-1)
     thgcp=SIN(slgclat1)*COS(deltanp)+COS(slgclat1)*COS(azim)*SIN(deltanp)
     IF(thgcp.LT.-1.0)THEN
        thgcp=-1.0
     ELSE IF(thgcp.GT.1.0)THEN
        thgcp=1.0
     ENDIF
     thgcp=ASIN(thgcp)
     phgcp=COS(deltanp)-SIN(thgcp)*SIN(slgclat1)
     phgcp=phgcp/(COS(thgcp)*COS(slgclat1))
     IF(phgcp.LT.-1.0)THEN
        phgcp=-1.0
     ELSE IF(phgcp.GT.1.0)THEN
        phgcp=1.0
     ENDIF
     phgcp=ACOS(phgcp)
     IF(azim.GE.0.0.AND.azim.LE.pi)THEN
        phgcp=ABS(phgcp)
     ELSE
        phgcp=-ABS(phgcp)
     ENDIF
     IF(slgclon1.LT.slgclon2)THEN
        phgcp=slgclon1+phgcp
     ELSE
        phgcp=slgclon1-phgcp
     ENDIF
     thgcp=thgcp*180.0/pi
     phgcp=phgcp*180.0/pi
     !
     ! Make sure point lies within model
     !
     IF(thgcp.GT.got+stol.OR.thgcp.LT.got-rgst*(nnt-1)-stol)THEN
        WRITE(6,*)'Great circle slice latitude'
        WRITE(6,*)'lies outside model bounds!'
        WRITE(6,*)'TERMINATING PROGRAM!!!'
        STOP
     ELSE IF(phgcp.LT.gop-stol.OR.phgcp.GT.gop+rgsp*(nnp-1)+stol)THEN
        WRITE(6,*)'Great circle slice longitude'
        WRITE(6,*)'lies outside model bounds!'
        WRITE(6,*)'TERMINATING PROGRAM!!!'
        STOP
     ENDIF
     DO k=1,nnx
        rnodeth=INT((got-thgcp)/rgst)+1
        u=ABS(got-thgcp-(rnodeth-1)*rgst)/rgst
        IF(rnodeth.EQ.nnt)THEN
           rnodeth=rnodeth-1
           u=1.0
        ENDIF
        rnoder=INT((k-1)*(nnr-1)/REAL(nnx-1))+1
        v=(k-1)*(nnr-1)/REAL(nnx-1)-(rnoder-1)
        IF(rnoder.EQ.nnr)THEN
           rnoder=rnoder-1
           v=1.0
        ENDIF
        rnodeph=INT((phgcp-gop)/rgsp)+1
        w=ABS(phgcp-gop-(rnodeph-1)*rgsp)/rgsp
        IF(rnodeph.EQ.nnp)THEN
           rnodeph=rnodeph-1
           w=1.0
        ENDIF
        CALL vbspline(rnoder,rnodeth,rnodeph,u,w,v,k,j)
     ENDDO
  ENDDO
!
! Now write out great circle slice to file.
!
  OPEN(UNIT=30,FILE=ofilegcv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
        WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)

  DEALLOCATE(vela, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
ENDIF

DEALLOCATE(veln, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
ENDIF
!
! Read in input parameters for traveltime grid
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ptfd
READ(10,'(a22)')ofiledt
READ(10,*)ptfns
READ(10,'(a22)')ofilenst
READ(10,*)ptfew
READ(10,'(a22)')ofileewt
!
! Open traveltime field file if required
!
IF(ptfd.EQ.1.OR.ptfns.EQ.1.OR.ptfew.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifilet,FORM='unformatted',STATUS='old')
   READ(20)gor,got,gop
   READ(20)nnr,nnt,nnp
   READ(20)rgsr,rgst,rgsp
   ALLOCATE(ttn(nnp,nnt,nnr), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL ttn'
   ENDIF
   DO i=1,nnp
      DO j=1,nnt
         DO k=1,nnr
            READ(20)ttn(i,j,k)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Finish off GMT boundary files
!
WRITE(50,*)nnp
WRITE(50,*)nnt
WRITE(50,*)int(-sldep)

WRITE(60,*)nnt
WRITE(60,*)nnr

if(int(slns).ne.slns)then
WRITE(60,"(F4.1)")slns
else
write(60,*)int(slns)
endif

WRITE(70,*)nnp
WRITE(70,*)nnr
if(int(slew).ne.slew)then
WRITE(70,"(F4.1)")slew
else
write(70,*)int(slew)
endif
CLOSE(50)
CLOSE(60)
CLOSE(70)
!
! Extract depth slice if required.
!
IF(ptfd.EQ.1)THEN
   idp=INT((gor-sldep)/rgsr)+1
   IF(idp.LT.1.OR.idp.GT.nnr)THEN
      WRITE(6,*)'Depth requested lies outside model'
      STOP
   ENDIF
!
!  Write out slice to GMT xyz file
!
   OPEN(UNIT=30,FILE=ofiledt,STATUS='unknown')
   DO i=1,nnp
      DO j=nnt,1,-1
         ttt=ttn(i,j,idp)
         IF(idp.EQ.nnr)THEN
            ttb=ttt
         ELSE
            ttb=ttn(i,j,idp+1)
         ENDIF
!
!        Apply linear interpolation to get travltime
!
         tt=ttt+(ttb-ttt)*((gor-(idp-1)*rgsr)-sldep)/rgsr
         WRITE(30,*)tt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
!
! Extract N-S slice if required.
!
IF(ptfns.EQ.1)THEN
   ins=INT((slns-gop)/rgsp)+1
   IF(ins.LT.1.OR.ins.GT.nnp)THEN
      WRITE(6,*)'Longitude requested lies outside model'
      STOP
   ENDIF
   OPEN(UNIT=30,FILE=ofilenst,STATUS='unknown')
   DO i=1,nnt
      DO j=nnr,1,-1
         ttt=ttn(ins,i,j)
         IF(ins.EQ.nnp)THEN
            ttb=ttt
         ELSE
            ttb=ttn(ins+1,i,j)
         ENDIF
!
!        Apply linear interpolation to get velocity
!
         tt=ttt+(ttb-ttt)*(slns-(gop+(ins-1)*rgsp))/rgsp
         WRITE(30,*)tt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
!
! Extract E-W slice if required.
!
IF(ptfew.EQ.1)THEN
   iew=INT((got-slew)/rgst)+1
   IF(iew.LT.1.OR.iew.GT.nnp)THEN
      WRITE(6,*)'Latitude requested lies outside model'
      STOP
   ENDIF
   OPEN(UNIT=30,FILE=ofileewt,STATUS='unknown')
   DO i=1,nnp
      DO j=nnr,1,-1
         ttt=ttn(i,iew,j)
         IF(iew.EQ.nnt)THEN
            ttb=ttt
         ELSE
            ttb=ttn(i,iew+1,j)
         ENDIF
!
!        Apply linear interpolation to get velocity
!
         tt=ttt+(ttb-ttt)*((got-(iew-1)*rgst)-slew)/rgst
         WRITE(30,*)ttt
      ENDDO
   ENDDO
   CLOSE(30)
ENDIF
IF(ptfd.EQ.1.OR.ptfns.EQ.1.OR.ptfew.EQ.1)THEN
   DEALLOCATE(ttn, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL ttn'
   ENDIF
ENDIF
!
! Read in ray path parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)prpd
READ(10,'(a22)')ofiledr
READ(10,*)prpns
READ(10,'(a22)')ofilensr
READ(10,*)prpew
READ(10,'(a22)')ofileewr
CLOSE(10)
!
! Plot raypaths in depth.
!
IF(prpd.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofiledr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlon,rlat
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
!
! Plot raypaths in N-S.
!
IF(prpns.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofilensr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlat,rdep
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
!
! Plot raypaths in E-W.
!
IF(prpew.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofileewr,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      DO j=1,nre
         READ(20)rdep,rlat,rlon
         WRITE(30,'(2f10.4)')rlon,rdep
      ENDDO
      WRITE(30,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF
STOP
END PROGRAM slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine uses cubic B-spline interpolation to sample
! a point wthiin a 3-D velocity grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vbspline(rnode,gt,gp,u,v,w,stt,stp)
USE globalp
IMPLICIT NONE
INTEGER :: i1,j1,k1,stp,stt,gt,gp
INTEGER :: rnode,lid
REAL(KIND=i10), DIMENSION(4) :: ui,vi,wi
REAL(KIND=i10) :: u,v,w,sumi,sumj,sumk,rdm
!
! u,v,w = independent surface parameters
! ui,vi,wi = bspline basis functions
! stt,stp = global grid coordinate of diced node
! sumi,sumj,sumk = summation variables for spline
! lid= layer id.
! gt,gp = i,j coordinate of current cell
! rnode = Depth coordinate of current cell
!
! Compute the values of the basis functions
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
! Set parameters for slice type.
!
sumi=0.0
DO i1=1,4
   sumj=0.0
   DO j1=1,4
      sumk=0.0
      DO k1=1,4
         rdm=wi(k1)*veln(rnode-2+k1,gt-2+j1,gp-2+i1)
         sumk=sumk+rdm
      ENDDO
      sumj=sumj+ui(j1)*sumk
   ENDDO
   sumi=sumi+vi(i1)*sumj
ENDDO
vela(stt,stp)=sumi
END SUBROUTINE vbspline
