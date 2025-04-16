c
c
c                V. Magas & A. Reina- Barcelona - APRIL 2022


c     ****************************************************************************
c
c     PROGRAM: EFFECTIVE STRING ROPE MODEL EVENT-BY-EVENT SIMULATIONS
c
c     COMMENTS: This code calculate heavy ion collision event-by-event at high energies 
c     using the Effective String Rope Model (arXiv:hep-ph/0202085v2) and the 
c     Glauber Monte Carlo approach (Annu. Rev. Nucl. Part. Sci. 2007. 57:20543)
c     It prints a fort.100 file used as initial conditions for hydro calculations.
c 
c     ****************************************************************************

c     We are going to work in the center of rapidity frame --> rap01=-rap02      

c 
      PROGRAM ESRM
      IMPLICIT NONE
C     ----------------------------------------------------------------------
      CHARACTER(8)  :: DATE
      CHARACTER(10) :: time
      CHARACTER(5)  :: zone
      REAL start,finish
      INTEGER,DIMENSION(8) :: values
      REAL, PARAMETER :: pi=3.1415926535897932384626
      REAL, PARAMETER :: hbarc=197.326 ! MeV*fm
      INTEGER :: A1,A2,tnc,nteventos,vortin,for100,check,nn1,nn2,nx,ny
     1 ,nz,Nceldasx,Nceldasy,Nceldasz,ismax,ismin,Npuntosal,ixmin,iymin,
     2 iiig,lg,i,j,jj,k,l,iii,kkk,iy,neventos
     3 ,qqq,ll,kk,i7,N_lines
      REAL :: epsilon00,Mn,r1,r2,rho00,rn,b0,A,aws1,aws2,
     1 x1,y1,x2,y2,e00,gamma0,v0,rap00,b,xsize,ysize,zsize,dx,dy,dz,
     2 z1,z2,ymax,ymin,xmax,xmin,xsmax,xsmin,tbctp,tpn2,alpha_string
     3 ,Bag,Co2,alpha,M,prom,xplusC1,xminusC2,MMMT,MMAXT,MINT,
     4 trajectory111,trajectory222,zzz1,zzz2,Tfin,ksi1,ksi2,v12,gamma12
     5 ,w12,v21,gamma21,w21,v22,gamma22,w22,v11,gamma11,w11,v10,gamma10,
     6 w10,v20,gamma20,w20,q1,q2,v00,gamma00,tnpnaverage2,
     7 Nnptot,Nnp,Nnpaver,Nnp_str,Nnp_str_aver,ener_str,zmom_str,
     8 ener_str_aver,zmom_str_aver,tEtp,tPztp,tEpn2,tPzpn2,end_time,
     9 ener_str2,ener_str3,tPznpnaverage2,Enptot,nnhydro,eehydro,Enp
      REAL :: hvc, hv2, amp, amn, amnuc, ak, an0, bin, a2j,
     1 a3j, a4j, a5j, bfer, ph3, ent, ampi, ampi0, deg, bbag, alam, aq1,
     2 aq2, aq3, tcr, agk, degdelt, entdelt, amdelta, ddd, dcoef, alam4,
     3 aq4, aq5, zq1, aaq1, alim, anflav, atq1, atq2, pi2, hv3, qq9, op9
      INTEGER :: mqcdwt 
      REAL :: tq, eq1, pq, sq, anq1, amuq, op9anq, anq13
      REAL :: totclassvort,totrelvort,eednn
      REAL, ALLOCATABLE, DIMENSION(:) :: x,y,Cellcenterx,Cellcentery,
     1 mtotclassvort,mtotrelvort
      REAL, DIMENSION(1:4001) :: z     
      REAL, DIMENSION(1:197,1:1000) :: xab1,yab1
      REAL, DIMENSION(1:197,1:1000) :: xab2,yab2 
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: efin,thyfin,eadd,
     1 zmomadd,padd,wadd,ycorr,press,vort_fin,vort_rel_fin,nfin
      REAL, ALLOCATABLE, DIMENSION(:,:) :: Ne1,Ne2,Ne,l1,l2,
     1 rho01,rho02,chy,yf,dens,chyz,Ee,Pze
      REAL, DIMENSION(1:122,1:82) :: t0,e01,e02,Rpl,Rmi,X_lim,MINT_lim
     1 ,aa1,aa2,b1,b2,c1,c2,d1,d2,tauminus,tauplus,thy,sigma_s,rap02,
     2 xturn1,xturn2,efin0,rap01
      INTEGER, ALLOCATABLE,DIMENSION(:,:) :: fkmin,fkmax,tfkmin,tfkmax,
     1 kz1min,kz2min,kz1max,kz2max
      REAL, DIMENSION(1:4001,1:122,1:82) :: xpm0,tC1,zC1,tC2,
     1 zC2
      REAL, DIMENSION(1:122,1:495,1:82) :: ehydro,vhydro,pzhydro,
     1 presshydro,nhydro,u,v,pxhydro,pyhydro
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: epaper,vpaper,pzpaper,
     1 presspaper,npaper
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: efintot,thyfintot,
     1 presstot,zmomtot,nfintot,paddtot,eaddtot,zmomaddtot,waddtot,
     2 T00hydro,T0Zhydro,N0hydro,NZhydro
      REAL, DIMENSION(1:123,1:83) :: MAXT
C     ----------------------------------------------------------------------
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
      COMMON /nrd/ aws1,aws2,A1,A2,dx,dy,Nceldasx,Nceldasy,nx,ny,r1,r2,
     1 rn,x1,x2,y1,y2,z1,z2
      COMMON /inxplusC1/ aa2,b1,d1,rap01,tauminus,e01
      COMMON /inxminusC2/ aa1,b2,d2,rap02,tauplus,e02
      COMMON /intraj111/ thy,MAXT,xturn1,tC1,zC1,z,t0
      COMMON /intraj222/ xturn2,tC2,zC2
      COMMON /cond/ check,nteventos
      COMMON /uniq/ mqcdwt
      COMMON /eqqcd/ tq, eq1, pq, sq, anq1, amuq
      COMMON /fq3/ op9anq, anq13



c     ----------------------------------------------------------------------
c                  READ INITIAL PARAMETERS FROM INPUT FILE
c
      OPEN(UNIT=1,FILE='input/global_parameters_ESRM_code',STATUS='old'
     1 ,ACTION='read')
c
         READ(1,151) epsilon00
         READ(1,152) b0
         READ(1,153) A1
         READ(1,154) A2
         READ(1,156) r1
         READ(1,157) r2
         READ(1,158) aws1
         READ(1,159) aws2
         READ(1,161) Mn
         READ(1,171) rn
         READ(1,162) A
         READ(1,163) rho00
         READ(1,164) tnc
         READ(1,166) nteventos
         READ(1,167) vortin
         READ(1,168) for100
         READ(1,169) check
         READ(1,211) end_time
      
      CLOSE(UNIT=1)
      
      IF (A1 .NE. 197) THEN
      IF (A1 .EQ. 2) r1=2.095
      IF (A1 .EQ. 3) r1=1.976
      IF (A1 .EQ. 4) r1=1.696
      IF (A1 .EQ. 7) r1=2.390
      IF (A1 .EQ. 10) r1=2.450
      IF (A1 .EQ. 13) r1=2.440
      IF (A1 .EQ. 16) r1=2.718
      IF (A1 .GE. 20 .AND. A1 .LT. 100) r1=1.2*A1**(1./3.)
      IF (A1 .GT. 100) r1=1.1*A1**(1./3.)
      IF (A1 .EQ. 100) r1=(1.1*(A1+1)**(1./3.)+1.2*(A1-1)**(1./3.))/2.
      ENDIF
      IF (A2 .NE. 197) THEN
      IF (A2 .EQ. 2) r2=2.095
      IF (A2 .EQ. 3) r2=1.976
      IF (A2 .EQ. 4) r2=1.696
      IF (A2 .EQ. 7) r2=2.390
      IF (A2 .EQ. 10) r2=2.450
      IF (A2 .EQ. 13) r2=2.440
      IF (A2 .EQ. 16) r2=2.718
      IF (A2 .GE. 20 .AND. A2 .LT. 100) r2=1.2*A2**(1./3.)
      IF (A2 .GT. 100) r2=1.1*A2**(1./3.)
      IF (A2 .EQ. 100) r2=(1.1*(A2+1)**(1./3.)+1.2*(A2-1)**(1./3.))/2.
      ENDIF
      
      WRITE(*,*)
      WRITE(*,*) 'Initial energy per nucleon [GeV]: ',epsilon00
      WRITE(*,*) 'Beam parameter: ',b0
      WRITE(*,*) 'Number of nucleons target: ',A1
      WRITE(*,*) 'Number of nucleons projectile: ',A2
      WRITE(*,*) 'Target radius [fm]: ',r1
      WRITE(*,*) 'Projectile radius [fm]: ',r2
      WRITE(*,*) 'aws1 [fm]: ',aws1
      WRITE(*,*) 'aws2 [fm]: ',aws2
      WRITE(*,*) 'Nucleon mass [GeV]: ',Mn
      WRITE(*,*) 'Nucleon radius [fm]: ',rn
      WRITE(*,*) 'String tension parameter A [GeV]: ',A
      WRITE(*,*) 'Normal nuclear density [fm-3]: ',rho00
      WRITE(*,*) 'tnc parameter: ',tnc
      WRITE(*,*) 'Number of events: ',nteventos
      IF (vortin .EQ. 1) THEN
        WRITE(*,*) 'Vorticity implementation: YES'
      ELSE
        WRITE(*,*) 'Vorticity implementation: NO'
      ENDIF
      IF (for100 .EQ. 1) THEN
        WRITE(*,*) 'Print fort.100 file:  YES'
      ELSE
        WRITE(*,*) 'Print fort.100 file:  NO'
      ENDIF
      WRITE(*,*)
c       WRITE(*,*)'IF THESE INPUT PARAMETERS ARE CORRECT TYPE 1, OTHERWISE
c      1 TYPE 0'
c       READ(*,*) i
c
c       IF (i .EQ. 1) GOTO 455
c       IF (i .EQ. 0) STOP
      
  455 CONTINUE 
       
     
c     ----------------------------------------------------------------------
      CALL CPU_TIME(start)
c     ----------------------------------------------------------------------
c           DATE AND TIME IN WHICH THE SIMULATION HAS STARTED 
      CALL DATE_AND_TIME(DATE,time,zone,values)
      CALL DATE_AND_TIME(DATE=DATE,ZONE=zone)
      CALL DATE_AND_TIME(TIME=time)
      CALL DATE_AND_TIME(VALUES=values)
      PRINT '(a,2x,a,2x,a)', DATE, time, zone
      PRINT '(8i5)', values

      
c     ----------------------------------------------------------------------
c      
      OPEN(UNIT=2,FILE='output/output1',STATUS='REPLACE')
C
C     ----------------------------------------------------------------------
C
      WRITE(2,182) 
      WRITE(2,183) epsilon00
      WRITE(2,184) b0
      WRITE(2,186) A1
      WRITE(2,187) A2
      WRITE(2,188) r1
      WRITE(2,189) r2
      WRITE(2,192) aws1
      WRITE(2,193) aws2
      WRITE(2,194) Mn
      WRITE(2,196) rn
      WRITE(2,198) A
      WRITE(2,199) rho00
      WRITE(2,201) tnc
      WRITE(2,*)
      WRITE(2,202) 
      WRITE(2,203) nteventos
      WRITE(2,204) vortin
      WRITE(2,206) for100
      WRITE(2,207) check
      WRITE(2,208)
c
c     ----------------------------------------------------------------------
c                          INITIAL PARAMETERS
c
      e00=rho00*Mn ! Initial energy density(GeV*fm-3). It is the same 
C     for both nuclei: target and projectile
      gamma0=epsilon00/Mn ! gamma factor
      v0=SQRT(1.-1./gamma0**2) ! Initial velocity for each nucleon
      rap00=ATANH(v0) ! Initial rapidity for each nucleon (just an 
C     absolute value sing assumed to be negative automatically)
      b=b0*(r1+r2) ! Impact parameter
c
      WRITE(2,*)
      WRITE(2,173)
      WRITE(2,174) e00
      WRITE(2,176) gamma0
      WRITE(2,177) v0
      WRITE(2,178) rap00
      WRITE(2,179) b 
      WRITE(2,181)
c
      WRITE(*,*)
      WRITE(*,*) 'Impact parameter=',b, '[fm]'
c
c     ------------------------------------------------------------------
c                           CELL SIZE
c
      xsize=MAX(r1,r2)/FLOAT(tnc)
      ysize=xsize
      zsize=xsize/100.
      dx=xsize
      dy=ysize
      dz=zsize
c
      WRITE(2,*)
      WRITE(2,209)
      WRITE(2,*)
      WRITE(2,106)
      WRITE(2,107) dx,dy,dz
c
c     ----------------------------------------------------------------------
c                INITIAL POSITION FOR EACH NUCLEUS
c
C     The initial (x,y) coordinates of the center of each nucleus is read 
C     from the input file. Here only the z-coordinate of each nucleus
C     is calculated.
C9
      WRITE(2,*)
      WRITE(2,108)
c                                 TARGET
c     The right one stays in the origen of the xy plane and displaced
C     in the z-direction by
c
      x1=0.0
      y1=0.
      z1=2.*r1/gamma0
      nn1=NINT(z1/zsize+0.49)
      z1=zsize*nn1
c                                PROJECTIL
      x2=b
      y2=0.
      z2=2.*r2/gamma0
      nn2=NINT(z2/zsize+0.49)
      z2=-zsize*nn2
c
      WRITE(2,109) x1,y1,z1
      WRITE(2,111) x2,y2,z2
c
c     ----------------------------------------------------------------------
c                  TRANSVERSE PLANE DEFINITION

      ymax=2.*MAX(r1,r2)+1.1*rn
      ymin=-ymax
      xmax=MAX(x1+2.*r1,x2+2.*r2)+1.1*rn
      xmin=MIN(x1-2.*r1,x2-2.*r2)-1.1*rn
c
      WRITE(2,*)
      WRITE(2,112)
      WRITE(2,113) xmin,xmax,ymin,ymax 
c
c     ----------------------------------------------------------------------
c          NUMBER OF NODES USED IN EACH DIRECTION OF THE GRID
c
      nx=NINT((xmax-xmin)/xsize)+2
      ny=NINT((ymax-ymin)/ysize)+2
      nz=4001
c
      WRITE(2,*)
      WRITE(2,114)
      WRITE(2,116) nx,ny,nz
c
C     ----------------------------------------------------------------------
C                    POSITION OF EACH NODE IN THE GRID
c
      ALLOCATE(x(1:nx))
      ALLOCATE(y(1:ny))
c
      DO 5 i=1,nz
c
          z(i)=0.
          z(i)=-(FLOAT(nz-1)/2.)*zsize+FLOAT(i-1)*zsize
          write (222,*) z(i)
c
        IF (i .GT. ny) GO TO 10    
c
          y(i)=0.
          y(i)=ymin+FLOAT(i-1)*ysize 
c
   10   CONTINUE
c
        IF (i .GT. nx) GO TO 5
c
          x(i)=0.
          x(i)=xmin+FLOAT(i-1)*xsize
c
    5 CONTINUE
      
c 
c     ----------------------------------------------------------------------
c       POSITION OF THE CENTER OF EACH CELL IN THE TRANSVERSE PLANE
c
      Nceldasx=nx-1
      Nceldasy=ny-1  
      Nceldasz=nz-1
c
      WRITE(2,*)
      WRITE(2,117)
      WRITE(2,118) Nceldasx,Nceldasy,Nceldasz
c
      ALLOCATE(Cellcenterx(1:Nceldasx))
      ALLOCATE(Cellcentery(1:Nceldasy))
c
      DO 15 i=1,Nceldasx
      
          Cellcenterx(i)=0.
          Cellcenterx(i)=(x(i)+x(i+1))/2.
          
        IF (i .GT. Nceldasy) GO TO 15
      
          Cellcentery(i)=0.
          Cellcentery(i)=(y(i)+y(i+1))/2.
          
   15 CONTINUE
c
c     ---------------------------------------------------------------
c                       PARTICIPANT REGION
C     We define the limits of the traansverse plane in the x-direction
c
      xsmax=MIN(x1+2.*r1,x2+2.*r2)
      xsmin=MAX(x1-2.*r1,x2-2.*r2)
c     We calculate the nearest nodes to xsmax and xsmin
      DO 20 i=1,nx
      
        IF (ABS(x(i)-xsmax) .LT. xsize/2.) ismax=i ! MAXIMUM
        IF (ABS(x(i)-xsmin) .LT. xsize/2.) ismin=i ! MINIMUM
C  
   20 CONTINUE
c
      WRITE(2,*)
      WRITE(2,119)
      WRITE(2,*) 
      WRITE(2,121) ismin,ismax
      WRITE(2,122) xsmin,xsmax
c
c     ----------------------------------------------------------------------
c
      CLOSE(UNIT=2)
c
c     ----------------------------------------------------------------------
c
c     ######################################################################
c           ESRM GLAUBER-MONTE CARLO EVENT-BY-EVENT SIMULATIONS
c     ######################################################################
c
      ALLOCATE(efintot(1:nx,1:nz,1:ny))
      ALLOCATE(thyfintot(1:nx,1:nz,1:ny))
      ALLOCATE(nfintot(1:nx,1:nz,1:ny))
      ALLOCATE(zmomtot(1:nx,1:nz,1:ny))
      ALLOCATE(eaddtot(1:nx,1:nz,1:ny))
      ALLOCATE(paddtot(1:nx,1:nz,1:ny))
      ALLOCATE(zmomaddtot(1:nx,1:nz,1:ny))
      ALLOCATE(waddtot(1:nx,1:nz,1:ny))
      ALLOCATE(presstot(1:nx,1:nz,1:ny))
      ALLOCATE(tfkmin(1:Nceldasx,1:Nceldasy))
      ALLOCATE(tfkmax(1:Nceldasx,1:Nceldasy))
      ALLOCATE(vort_fin(1:nx,1:nz,1:ny))
      ALLOCATE(vort_rel_fin(1:nx,1:nz,1:ny))
c
c     ----------------------------------------------------------------------
c
      
      efintot=0.
      thyfintot=0.
      nfintot=0.
      zmomtot=0.
      presstot=0.
      eaddtot=0.
      paddtot=0.
      zmomaddtot=0.
      waddtot=0.   
      tfkmin=1
      tfkmax=1
c
      tnpnaverage2=0.
      tPznpnaverage2=0.
      Nnpaver=0.
      Nnp_str_aver=0.
      ener_str_aver=0.
      zmom_str_aver=0.
      neventos=0


C     ######################################################################
C                            START SIMULATIONS
c     ######################################################################
c
      
      DO 999 qqq=1,nteventos
c
c     ----------------------------------------------------------------------
c
        neventos=neventos+1
c
        WRITE(*,133) neventos   
c
c     ----------------------------------------------------------------------
c
  140   CONTINUE
  998   CONTINUE
        
c
c     ----------------------------------------------------------------------
c
        CALL nucl_rand_dist(xab1,yab1,xab2,yab2)
c
c     ----------------------------------------------------------------------
c
        IF (neventos .GT. 1) GO TO 54
c
          ALLOCATE(efin(1:nx,1:nz,1:ny))
          ALLOCATE(thyfin(1:nx,1:nz,1:ny))
          ALLOCATE(eadd(1:nx,1:nz,1:ny))
          ALLOCATE(zmomadd(1:nx,1:nz,1:ny))
          ALLOCATE(padd(1:nx,1:nz,1:ny))
          ALLOCATE(wadd(1:nx,1:nz,1:ny))
          ALLOCATE(ycorr(1:nx,1:nz,1:ny))
          ALLOCATE(press(1:nx,1:nz,1:ny))
          ALLOCATE(nfin(1:nx,1:nz,1:ny))
          ALLOCATE(Ne1(1:Nceldasx,1:Nceldasy))
          ALLOCATE(Ne2(1:Nceldasx,1:Nceldasy))
          ALLOCATE(Ne(1:Nceldasx,1:Nceldasy))
          ALLOCATE(Ee(1:Nceldasx,1:Nceldasy))
          ALLOCATE(Pze(1:Nceldasx,1:Nceldasy))
          ALLOCATE(l1(1:Nceldasx,1:Nceldasy))
          ALLOCATE(l2(1:Nceldasx,1:Nceldasy))
          ALLOCATE(kz1min(1:Nceldasx,1:Nceldasy))
          ALLOCATE(kz1max(1:Nceldasx,1:Nceldasy))
          ALLOCATE(kz2min(1:Nceldasx,1:Nceldasy))
          ALLOCATE(kz2max(1:Nceldasx,1:Nceldasy))
          ALLOCATE(fkmin(1:Nceldasx,1:Nceldasy))
          ALLOCATE(fkmax(1:Nceldasx,1:Nceldasy))
          ALLOCATE(rho01(1:Nceldasx,1:Nceldasy))
          ALLOCATE(rho02(1:Nceldasx,1:Nceldasy))
          ALLOCATE(chy(1:Nceldasx,1:Nceldasy))
          ALLOCATE(yf(1:Nceldasx,1:Nceldasy))
          ALLOCATE(dens(1:Nceldasx,1:Nceldasy))
          ALLOCATE(chyz(1:Nceldasx,1:Nceldasy))
          ALLOCATE(mtotclassvort(1:nteventos))
          ALLOCATE(mtotrelvort(1:nteventos))
c
   54   CONTINUE
    
c
c     ----------------------------------------------------------------------
c
        thyfin=0.
        nfin=0.
        efin=0.
        efin0=0.
        eadd=0.
        zmomadd=0.
        ycorr=0.
        padd=0.
        wadd=0.
                press=0
                Ne1=0.
                Ne2=0.
                l1=0.
                l2=0.
                Ne=0.
                Ee=0.
                Pze=0.
                sigma_s=0.    
                kz1min=2001
                kz1max=2001
                kz2min=2001-4*nn2
                kz2max=2001-4*nn2
                fkmin=0
                fkmax=0
                e01=0.
                e02=0.
                rap01=0.
                rho01=0.
                rap02=0.
                rho02=0.
                aa1=0.
                aa2=0.
                b1=0.
                b2=0.
                c1=0.
                c2=0.
                d1=0.
                d2=0.
                tauminus=0.
                tauplus=0.
                chy=0.
                thy=0.
                yf=0.
                dens=0.
                chyz=0.
                xturn1=0.
                xturn2=0.
                t0=0.
                xpm0=0.
                tC1=0.
                zC1=0.
                tC2=0.
                zC2=0.   

c
c     ----------------------------------------------------------------------
c      POSITION OF EACH FICTITIOUS PARTICLE IN THE TRANSVERSE PLANE           
c           
        Npuntosal=1000
c
        DO 85 k=1,MAX0(A1,A2) ! We sum over all nucleons within a nucleus 
          DO 90 l=1,Npuntosal ! We sum over all fictitious particles 
C     within each nucleon
c                            TARGET
c
            IF (k .GT. A1) GO TO 95
c
              i=FLOOR((-x(1)+xab1(k,l))/dx)+1 ! Position of each fictitious 
C     particle in the x-direction of the grid
              j=FLOOR((-y(1)+yab1(k,l))/dy)+1 ! Position of each fictitious 
C     particle in the y-direction of the grid
              Ne1(i,j)=Ne1(i,j)+1./FLOAT(Npuntosal) ! Baryon charge within each cell
C     corresponding to the target
c
   95       CONTINUE
c
c                          PROJECTILE
c
            IF (k .GT. A2) GO TO 90
c
              i=FLOOR((-x(1)+xab2(k,l))/dx)+1
              j=FLOOR((-y(1)+yab2(k,l))/dy)+1 
              Ne2(i,j)=Ne2(i,j)+1./FLOAT(Npuntosal) ! Baryon charge within each cell 
C     corresponding to the projectile
c
   90     CONTINUE
   85   CONTINUE
     
        tbctp = SUM(Ne1)+SUM(Ne2) ! Total baryon charge
        tEtp = epsilon00*tbctp  ! Total initial energy
        tPztp = -Mn*gamma0*v0*SUM(Ne1)+Mn*gamma0*v0*SUM(Ne2) ! Total initial momentum
        
        WRITE(*,*)
        WRITE(*,124) tbctp,tEtp,tPztp
        WRITE(*,*)
        
c
c     ---------------------------------------------------------------
c                  STRING TENSION AND LENGTH OF EACH STREAK
C
        DO 115 j=1,Nceldasy
          DO 120 i=ismin,ismax

            l1(i,j)=Ne1(i,j)/rho00/dx/dy/gamma0
            l2(i,j)=Ne2(i,j)/rho00/dx/dy/gamma0

            IF (l1(i,j) .GT. 0.05*zsize .AND. l2(i,j) .GT. 
     1 0.05*zsize) THEN      
C     
              Ne(i,j)=Ne1(i,j)+Ne2(i,j)
              Ee(i,j)=epsilon00*Ne(i,j)
              Pze(i,j)=-Mn*gamma0*v0*Ne1(i,j)+Mn*gamma0*v0*Ne2(i,j)
C
            ENDIF
C
  120     CONTINUE
  115   CONTINUE

         tpn2 = SUM(Ne) ! Number of participant nucleons
         tEpn2 = SUM(Ee) ! Total energy
         tPzpn2 = SUM(Pze) ! Total momentum
        
        WRITE(*,*)
        WRITE(*,144) tpn2,tEpn2,tPzpn2
        WRITE(*,*)
      
c      #####################################################################
c               GENERALIZED EFFECTIVE STRING ROPE MODEL
C      #####################################################################
C
c     ----------------------------------------------------------------------
c              INITIAL (COLLISIONS) COORDINATE FOR CELLS (z=t)
c
c     We calculate the cells, in z-direction, where are the extremes of 
c          each initial streak corresponding to each nucleus 
C
        DO 145 l=1,Nceldasy
          DO 150 iii=1,Nceldasx 
C
      KZ1:  IF (l1(iii,l) .GT. 0.) THEN
C
              k=2001-2*nn1   
C
  155         IF (z(k) .LT. z1-l1(iii,l)) THEN
C
                kz1min(iii,l)=k+1
                k=k+1           
                GO TO 155     
C
              ENDIF      
C
              k=kz1min(iii,l)
C
  160         IF (z(k) .LE. z1+l1(iii,l)) THEN
C
                kz1max(iii,l)=k
                k=k+1           
                GO TO 160           
C
              ENDIF  
C
            ELSE
C
              kz1max(iii,l)=kz1min(iii,l)-1    
C
            ENDIF KZ1
C
      KZ2:  IF (l2(iii,l) .GT. 0.) THEN
C
              k=2001-4*nn2    
C
  165         IF (z(k) .LT. z2-l2(iii,l)) THEN
C
                kz2min(iii,l)=k+1
                k=k+1
                GO TO 165     
C
              ENDIF       
C
              k=kz2min(iii,l)
C
  170         IF (z(k) .LE. z2+l2(iii,l)) THEN
C
                kz2max(iii,l)=k
                k=k+1              
                GO TO 170
C
              ENDIF    
C
            ELSE    
C
              kz2max(iii,l)=kz2min(iii,l)-1       
C
            ENDIF KZ2
C
  150     CONTINUE
  145   CONTINUE
C
c     ----------------------------------------------------------------------
c                           RECOIL MODEL
c
c     String production const depend on energy quadratically
C
        alpha_string=(epsilon00/Mn)**2
        Bag=0.330 ! Bag constant (GeV per fm**3)
        Co2=1./3. ! Speed of sound to square --> Take from EoS
        alpha=(1.-Co2)/(1.+Co2)
c
        DO 175 lg=1,Nceldasy
          DO 180 iiig=ismin,ismax  
c
      COND1:IF (l1(iiig,lg) .GT. 0.05*zsize .AND. 
     1 l2(iiig,lg) .GT. 0.05*zsize) THEN
c
c     Old definition of string tension:
c     sigma_s(iiig,lg)=alpha_string*A*rho00*sqrt(l1(iiig,lg)*l2(iiig,lg))

              sigma_s(iiig,lg)=(A*gamma0*(Ne1(iiig,lg)*Ne2(iiig,lg))
     1 **(1./4.))/SQRT(dx*dy) ! New definition of sigma
              t0(iiig,lg)=(l1(iiig,lg)+l2(iiig,lg))/2.
c
              DO 185 k=kz1min(iiig,lg),kz1max(iiig,lg)
c
                xpm0(k,iiig,lg)=2.*t0(iiig,lg)-ABS(z(k)) 
c
  185         CONTINUE
c
              DO 190 k=kz2min(iiig,lg),kz2max(iiig,lg)  
c
                xpm0(k,iiig,lg)=2.*t0(iiig,lg)-ABS(z(k))  
c
  190         CONTINUE
c
              e01(iiig,lg)=(rho00*Mn)/(1.+Co2)-sigma_s(iiig,lg)*rho00
     1 *EXP(-rap00)*l1(iiig,lg)/4./(epsilon00/Mn)**2/(1.+Co2)-
     2 (Bag+sigma_s(iiig,lg)**2/0.197/2.)*(l1(iiig,lg)+l2(iiig,lg))/
     3 (epsilon00/Mn)**2/(1.+Co2)/2./l1(iiig,lg)
              e02(iiig,lg)=(rho00*Mn)/(1.+Co2)-sigma_s(iiig,lg)*rho00
     1 *EXP(-rap00)*l2(iiig,lg)/4./(epsilon00/Mn)**2/(1.+Co2)-
     2 (Bag+sigma_s(iiig,lg)**2/0.197/2.)*(l1(iiig,lg)+l2(iiig,lg))/
     3 (epsilon00/Mn)**2/(1.+Co2)/2./l2(iiig,lg)
c
              IF (e01(iiig,lg) .LT. 0.) e01(iiig,lg)=0.
              IF (e02(iiig,lg) .LT. 0.) e02(iiig,lg)=0.
c
c     the right sign of rap01 taking into account manually
c
              rap01(iiig,lg)=rap00
              rho01(iiig,lg)=rho00
              rap02(iiig,lg)=rap00
              rho02(iiig,lg)=rho00 
c
c     According to my notes notations
c     Exact Physics
c
              c1(iiig,lg)=-sigma_s(iiig,lg)*rho00*EXP(rap00)/8./
     1 (1.+Co2)/(epsilon00/Mn)**2
      c2(iiig,lg)=c1(iiig,lg)
              aa1(iiig,lg)=c1(iiig,lg)-2.*sigma_s(iiig,lg)*rho00*
     1 EXP(rap00)*EXP(-2.*rap01(iiig,lg))+2.*sigma_s(iiig,lg)*rho00*
     2 EXP(rap00)
              aa2(iiig,lg)=c2(iiig,lg)-2.*sigma_s(iiig,lg)*rho00*
     1 EXP(rap00)*EXP(-2.*rap02(iiig,lg))+2.*sigma_s(iiig,lg)*rho00*
     2 EXP(rap00)
              d1(iiig,lg)=c1(iiig,lg)-2.*sigma_s(iiig,lg)*rho00*
     1 EXP(rap00)*EXP(-2.*rap01(iiig,lg))
              d2(iiig,lg)=c2(iiig,lg)-2.*sigma_s(iiig,lg)*rho00*
     1 EXP(rap00)*EXP(-2.*rap02(iiig,lg))
              b1(iiig,lg)=alpha*aa2(iiig,lg)+2.*sigma_s(iiig,lg)*
     1 rho00*EXP(rap00)
              b2(iiig,lg)=alpha*aa1(iiig,lg)+2.*sigma_s(iiig,lg)*
     1 rho00*EXP(rap00)
c
c     NEW FORMULAS !!!!!!!!!!!!!!!!!!!!!!!!!
c
              tauminus(iiig,lg)=e01(iiig,lg)*(1.+Co2)*
     1 EXP(2.*rap01(iiig,lg))/aa2(iiig,lg)
              tauplus(iiig,lg)=e02(iiig,lg)*(1.+Co2)*
     1 EXP(2.*rap02(iiig,lg))/aa1(iiig,lg)
c
c     Motion in common frame
c        FINAL RAPIDITY FOR AN HOMOGENEOUS DISTRIBUTION OF MATTER
c
      COND2:  IF (l2(iiig,lg) .NE. l1(iiig,lg)) THEN
c
                M=(l2(iiig,lg)+l1(iiig,lg))/(l2(iiig,lg)-l1(iiig,lg))
                prom=M**2*(1.+Co2)-2.*Co2
                chy(iiig,lg)=SQRT((prom+SQRT(prom**2+4.*Co2**2*
     1 (M**2-1.)))/2./(1.+Co2)/(M**2-1.)) ! Final rapidity using Landau's convection
c
               ELSEIF (l2(iiig,lg) .EQ. l1(iiig,lg)) THEN   
c
                chy(iiig,lg)=1. ! The streaks are at rest, yf=0 --> cosh(0)=1
c
               ENDIF COND2  
c
      COND3:  IF (l2(iiig,lg) .GT. l1(iiig,lg)) THEN ! Final streak moves to 
c     the right (yf>0) 
                thy(iiig,lg)=SQRT(chy(iiig,lg)**2-1.)/chy(iiig,lg) ! Final 
c     velocity (vf=tanh(yf)=thy)
c
              IF (thy(iiig,lg) .GE. 1.) THEN
c
                WRITE(*,*) 'thy>1',thy(iiig,lg)
                thy(iiig,lg)=0.99
              ENDIF
c
              IF (thy(iiig,lg) .LE. -1.) THEN
c
                WRITE(*,*) 'thy<-1',thy(iiig,lg)
                thy(iiig,lg)=-0.99
c
              ENDIF
c
              yf(iiig,lg)=ATANH(thy(iiig,lg)) ! Final rapidity
c
            ELSE ! Final streak moves to the left (yf<0)
c
              thy(iiig,lg)=-SQRT(chy(iiig,lg)**2-1.)/chy(iiig,lg)
c
              IF (thy(iiig,lg) .GE. 1.) THEN
c
                WRITE(*,*) 'thy>1',thy(iiig,lg)
                thy(iiig,lg)=0.99
c
              ENDIF
c
              IF (thy(iiig,lg) .LE. -1.) THEN
c
               WRITE(*,*) 'thy<-1',thy(iiig,lg)
               thy(iiig,lg)=-0.99
c
              ENDIF
c
              yf(iiig,lg)=ATANH(thy(iiig,lg))
c
            ENDIF COND3
c
            chyz(iiig,lg)=1./SQRT(1.-thy(iiig,lg)**2) ! Lorentz gamma factor (chyz=gamma=chy)
c        Position of the final streak formation, just before the final streaks 
c        starts to move like one object with uniform rapidity yf (see PRC64(014901))   
c
c         left extreme of each final streak
            xturn1(iiig,lg)=tauminus(iiig,lg)*(1.-((d1(iiig,lg)/
     1 b1(iiig,lg)+EXP(-2.*rap01(iiig,lg)))/(d1(iiig,lg)/b1(iiig,lg)+
     2 EXP(2.*yf(iiig,lg))))**(alpha*aa2(iiig,lg)/b1(iiig,lg)))
c          right extreme of each final streak
            xturn2(iiig,lg)=tauplus(iiig,lg)*(1.-((d2(iiig,lg)/
     1 b2(iiig,lg)+EXP(-2.*rap02(iiig,lg)))/(d2(iiig,lg)/b2(iiig,lg)+
     2 EXP(-2.*yf(iiig,lg))))**(alpha*aa1(iiig,lg)/b2(iiig,lg)))
c
c          TRAJECTORIES OF PARTONS (OR CELL ELEMENTS) FOR BOTH NUCLEI
c
            DO 195 k=kz1min(iiig,lg),kz1max(iiig,lg)
c
              tC1(k,iiig,lg)=0.5*(xplusC1(xturn1(iiig,lg),z(k))+
     1 xpm0(k,iiig,lg)+xturn1(iiig,lg))
c
              IF (ISNAN(tC1(k,iiig,lg)) .EQV. .TRUE.) tC1(k,iiig,lg)=0.
c
              zC1(k,iiig,lg)=0.5*(xplusC1(xturn1(iiig,lg),z(k))-
     1 xpm0(k,iiig,lg)-xturn1(iiig,lg)) ! Trajectories of partons for the target
c
              IF (ISNAN(zC1(k,iiig,lg)) .EQV. .TRUE.) zC1(k,iiig,lg)=0.
c
  195       CONTINUE
c
            DO 200 k=kz2min(iiig,lg),kz2max(iiig,lg)
c
              tC2(k,iiig,lg)=0.5*(xturn2(iiig,lg)+xpm0(k,iiig,lg)+
     1 xminusC2(xturn2(iiig,lg),z(k)))
c
              IF (ISNAN(tC2(k,iiig,lg)) .EQV. .TRUE.) tC2(k,iiig,lg)=0.
c
                zC2(k,iiig,lg)=0.5*(xturn2(iiig,lg)+xpm0(k,iiig,lg)-
     1 xminusC2(xturn2(iiig,lg),z(k))) ! Trajectories of partons for the projectile
c
              IF (ISNAN(zC2(k,iiig,lg)) .EQV. .TRUE.) zC2(k,iiig,lg)=0.
c
  200       CONTINUE
c
            ENDIF COND1
c
  180     CONTINUE
  175   CONTINUE
  
c
c     ######################################################################
c                         EXPANDING FINAL STREAKS
c     ######################################################################
c
        MMMT=0.
c
        DO 205 l=1,Nceldasy
          DO 210 iii=ismin,ismax   
c
      COND4:IF (l1(iii,l) .GT. 0.05*zsize .AND. l2(iii,l) .GT. 0.05*
     1 zsize) THEN
c
              MMAXT=0.
c
              DO 215 k=kz1min(iii,l),kz1max(iii,l) 
c
                IF (tC1(k,iii,l) .GT. MMAXT) MMAXT=tC1(k,iii,l) 
c
  215         CONTINUE
c
              DO 220 k=kz2min(iii,l),kz2max(iii,l)
c
                IF (tC2(k,iii,l) .GT. MMAXT) MMAXT=tC2(k,iii,l)  
c
  220         CONTINUE
c             
              MAXT(iii,l)=0.
              MAXT(iii,l)=MMAXT
c              
              IF (MMMT .LT. MMAXT) MMMT=MMAXT   
c
            ENDIF COND4    
c
  210     CONTINUE      
  205   CONTINUE
c
        WRITE(*,126) MMMT 
c
c     ----------------------------------------------------------------------
c
        MINT=1000.
        ixmin=0
        iymin=0
        Nnp_str=0.
c
C	  Volodya test Streak
c
        ener_str=0.
        ener_str2=0.
        ener_str3=0.
        zmom_str=0.
c
        DO 225 lg=1,Nceldasy
          DO 230 iiig=ismin,ismax
c
      COND5:IF (l1(iiig,lg) .GT. 0.05*zsize .AND. l2(iiig,lg) .GT. 
     1 0.05*zsize) THEN
c 
              zzz1=trajectory111(MAXT(iiig,lg),kz1min(iiig,lg))
              zzz2=trajectory222(MAXT(iiig,lg),kz2max(iiig,lg))
              fkmax(iiig,lg)=2001+NINT(MAX(zzz2,zzz1)/zsize)
              fkmin(iiig,lg)=2001+NINT(MIN(zzz2,zzz1)/zsize)
              Rpl(iiig,lg)=FLOAT((fkmax(iiig,lg)-2001))*zsize
              Rmi(iiig,lg)=FLOAT((fkmin(iiig,lg)-2001))*zsize
              X_lim(iiig,lg)=(Rpl(iiig,lg)*(1.-thy(iiig,lg)*
     1 SQRT(Co2))*(thy(iiig,lg)+SQRT(Co2))-Rmi(iiig,lg)*
     2 (1.+thy(iiig,lg)*SQRT(Co2))*(thy(iiig,lg)-SQRT(Co2)))/2./
     3 SQRT(Co2)/(1.-thy(iiig,lg)**2)
              MINT_lim(iiig,lg)=MAXT(iiig,lg)+(X_lim(iiig,lg)-
     1 Rpl(iiig,lg))*(1.-thy(iiig,lg)*SQRT(Co2))/(thy(iiig,lg)-
     2 SQRT(Co2))
              dens(iiig,lg)=(Ne1(iiig,lg)+Ne2(iiig,lg))*
     1 SQRT(1.-thy(iiig,lg)**2)/(Rpl(iiig,lg)-Rmi(iiig,lg))/xsize/ysize
c
              IF (MINT .GT. MINT_lim(iiig,lg)) THEN
c
                MINT=MINT_lim(iiig,lg)
                ixmin=iiig
                iymin=lg
c
              ENDIF
c
c            Final energy density, homogenious distribution
c
              efin0(iiig,lg)=((epsilon00/Mn)*rho00*epsilon00*
     1 (l2(iiig,lg)+l1(iiig,lg))-Bag*FLOAT((fkmax(iiig,lg)-
     2 fkmin(iiig,lg)))*zsize)/FLOAT((fkmax(iiig,lg)-fkmin(iiig,lg)))/
     3 ((1.+Co2)*chy(iiig,lg)**2-Co2)/zsize
c   
            ENDIF COND5 
c
  230     CONTINUE
  225   CONTINUE
c
c     ----------------------------------------------------------------------
c
      WRITE(*,131) MINT
      IF (MINT .LT. MMMT) WRITE(*,127)
c
c     ----------------------------------------------------------------------
c                       Give me Tfin >= MMMT
c                       Give me Tfin >= MMMT
c
        Tfin=end_time
C      
        IF (Tfin .LT. MMMT) WRITE(*,128)
C  
      COND6:IF (Tfin .lt. MINT) THEN
C      
              WRITE(*,129)
C      
        DO 235 lg=1,Nceldasy
          DO 240 iiig=ismin,ismax
C          
      COND7:IF (l1(iiig,lg) .GT. 0.05*zsize .AND.
     1 l2(iiig,lg) .GT. 0.05*zsize) THEN
C    
              zzz1=trajectory111(Tfin,kz1min(iiig,lg))
              zzz2=trajectory222(Tfin,kz2max(iiig,lg))
C                  
c              Final energy density, homogenious distribution
C
              fkmax(iiig,lg)=2001+NINT(MAX(zzz2,zzz1)/zsize)+1
              fkmin(iiig,lg)=2001+NINT(MIN(zzz2,zzz1)/zsize)-1
       IF (fkmax(iiig,lg) .GT. tfkmax(iiig,lg)) tfkmax(iiig,lg)=
     1 fkmax(iiig,lg)
       IF (fkmin(iiig,lg) .LT. tfkmin(iiig,lg)) tfkmin(iiig,lg)=
     1 fkmin(iiig,lg)
C                  
      COND8:  IF (Tfin .LE. MAXT(iiig,lg)) THEN
C                  
                DO 245 k=fkmin(iiig,lg),fkmax(iiig,lg)

c                       if (k .GT. 4001 .OR. k .LT.0) write (*,*)k
C                    
                   
                  efin(iiig,k,lg)=((epsilon00/Mn)*rho00*epsilon00*
     1 (l2(iiig,lg)+l1(iiig,lg))-Bag*FLOAT((fkmax(iiig,lg)-
     2 fkmin(iiig,lg)))*zsize)/FLOAT((fkmax(iiig,lg)-fkmin(iiig,lg)))/
     3 ((1.+Co2)*chy(iiig,lg)**2-Co2)/zsize
                  thyfin(iiig,k,lg)=thy(iiig,lg)
                  nfin(iiig,k,lg)=(Ne1(iiig,lg)+Ne2(iiig,lg))*
     1 SQRT(1.-thy(iiig,lg)**2)/FLOAT((fkmax(iiig,lg)-fkmin(iiig,lg)))
     2 /dz/xsize/ysize
                  press(iiig,k,lg)=Co2*efin(iiig,k,lg)-(1+Co2)*Bag
c
                  IF (press(iiig,k,lg).LT.0.)  press(iiig,k,lg)=0.
C
                  IF (for100 .NE. 1) GO TO 655
C
c                              USE TO PREPARE FORT.100
C
                    eadd(iiig,k,lg)=
     1 efin(iiig,k,lg)*((1.+Co2)/(1.-thyfin(iiig,k,lg)**2)-Co2)  ! T00
                    zmomadd(iiig,k,lg)=   
     1 efin(iiig,k,lg)*(1.+Co2)*thyfin(iiig,k,lg)/
     2 (1.-thyfin(iiig,k,lg)**2)                                   !T0Z
                    padd(iiig,k,lg)=   
     1 nfin(iiig,k,lg)/SQRT(1.-thyfin(iiig,k,lg)**2)              !N0
                    wadd(iiig,k,lg)=
     1 padd(iiig,k,lg)*thyfin(iiig,k,lg)                            !Nz
                    Nnp_str=Nnp_str+padd(iiig,k,lg)
C
C	                   Volodya test Streak
C
                    ener_str=ener_str+eadd(iiig,k,lg)
c      ener_str3=ener_str3+efin(iiig,k,lg)*chyz(iiig,lg)*dx*dy*dz
                    zmom_str=zmom_str+zmomadd(iiig,k,lg)
C
  655             CONTINUE     
C
  245           CONTINUE
C
              ELSE
c
                DO 250 k=fkmin(iiig,lg),fkmax(iiig,lg)
c      
                  ksi1=(z(k)-Rpl(iiig,lg))/(Tfin-MAXT(iiig,lg))
                  ksi2=(z(k)-Rmi(iiig,lg))/(Tfin-MAXT(iiig,lg)) 
c
      COND9:      IF (ksi1 .GT. 1.) THEN
c
                  ELSEIF (ksi2 .LT. -1.) THEN
c
                  ELSE  
c
      COND10:       IF (ksi1 .GT. (thy(iiig,lg)-
     1 SQRT(Co2))/(1.-thy(iiig,lg)*SQRT(Co2))) THEN
c
                      
                      efin(iiig,k,lg)=efin0(iiig,lg)*
     1 ((1.+thy(iiig,lg))*(1.-SQRT(Co2))*(1.-ksi1)/(1.-thy(iiig,lg))/
     2 (1.+SQRT(Co2))/(1.+ksi1))**((1.+Co2)/2./SQRT(Co2))
                      press(iiig,k,lg)=Co2*efin(iiig,k,lg)-(1+Co2)*Bag
C
                      IF (press(iiig,k,lg).LT.0.)  press(iiig,k,lg)=0.
C
                      thyfin(iiig,k,lg)=(ksi1+
     1 SQRT(Co2))/(1.+ksi1*SQRT(Co2))
C
                      IF (thyfin(iiig,k,lg) .GE. 1.) THEN
C
                        WRITE(*,*) 'thyfin>1',thyfin(iiig,k,lg)
                        thyfin(iiig,k,lg)=0.99
C
                      ENDIF
C
                      IF (thyfin(iiig,k,lg) .LE. -1.) THEN
C
                        WRITE(*,*) 'thyfin<-1',thyfin(iiig,k,lg)
                        thyfin(iiig,k,lg)=-0.99
C
                      ENDIF
C
                      ycorr(iiig,k,lg)=
     1 (1.-thyfin(iiig,k,lg))*(1.+thy(iiig,lg))/(1.-thy(iiig,lg))/
     2 (1.+thyfin(iiig,k,lg))
                      nfin(iiig,k,lg)=   
     1 dens(iiig,lg)*ycorr(iiig,k,lg)**(1./2./SQRT(Co2))    
C
                      IF (for100 .NE. 1) GO TO 255
C
c                     USE TO PREPARE FORT.100
C
                        eadd(iiig,k,lg)=
     1 efin(iiig,k,lg)*((1.+Co2)/(1.-thyfin(iiig,k,lg)**2)-Co2)  ! T00
                        zmomadd(iiig,k,lg)=   
     1 efin(iiig,k,lg)*(1.+Co2)*thyfin(iiig,k,lg)/
     2 (1.-thyfin(iiig,k,lg)**2)                                   !T0Z
                        padd(iiig,k,lg)=   
     1 nfin(iiig,k,lg)/SQRT(1.-thyfin(iiig,k,lg)**2)              !N0
                        wadd(iiig,k,lg)=
     1 padd(iiig,k,lg)*thyfin(iiig,k,lg)                            !Nz
                        Nnp_str=Nnp_str+padd(iiig,k,lg)
C
C	                      Volodya test Streak
C
                        ener_str=ener_str+eadd(iiig,k,lg)
c      ener_str3=ener_str3+efin(iiig,k,lg)*chyz(iiig,lg)*dx*dy*dz
                        zmom_str=zmom_str+zmomadd(iiig,k,lg)
C
  255                 CONTINUE  
C
                    ELSEIF (ksi2 .LT. (thy(iiig,lg)+
     1 SQRT(Co2))/(1.+thy(iiig,lg)*SQRT(Co2))) THEN
C
                      
                      efin(iiig,k,lg)=efin0(iiig,lg)*
     1 ((1.-thy(iiig,lg))*(1.-SQRT(Co2))*(1.+ksi2)/(1.+thy(iiig,lg))/(
     2 1.+SQRT(Co2))/(1.-ksi2))**((1.+Co2)/2./SQRT(Co2))
                      press(iiig,k,lg)=Co2*efin(iiig,k,lg)-(1+Co2)*Bag
C
                      IF (press(iiig,k,lg).LT.0.)  press(iiig,k,lg)=0.
C
                      thyfin(iiig,k,lg)=(ksi2-
     1 SQRT(Co2))/(1.-ksi2*SQRT(Co2))
C
                      IF (thyfin(iiig,k,lg) .GE. 1.) THEN
C
                        WRITE(*,*) 'thyfin>1',thyfin(iiig,k,lg)
                        thyfin(iiig,k,lg)=0.99
C
                      ENDIF
C
                      IF (thyfin(iiig,k,lg) .LE. -1.) THEN
C
                        WRITE(*,*) 'thyfin<-1',thyfin(iiig,k,lg)
                        thyfin(iiig,k,lg)=-0.99
C
                      ENDIF
C
                      ycorr(iiig,k,lg)=
     1 (1.+thyfin(iiig,k,lg))*(1.-thy(iiig,lg))/(1.+thy(iiig,lg))/
     2 (1.-thyfin(iiig,k,lg))
                      nfin(iiig,k,lg)=
     1 dens(iiig,lg)*ycorr(iiig,k,lg)**(1./2./SQRT(Co2))
C
                      IF (for100 .NE. 1) GO TO 260
C
c                    USE TO PREPARE FORT.100
C
                        eadd(iiig,k,lg)=
     1 efin(iiig,k,lg)*((1.+Co2)/(1.-thyfin(iiig,k,lg)**2)-Co2)  ! T00
                        zmomadd(iiig,k,lg)=
     1 efin(iiig,k,lg)*(1.+Co2)*thyfin(iiig,k,lg)/
     2 (1.-thyfin(iiig,k,lg)**2)
                        padd(iiig,k,lg)=
     1 nfin(iiig,k,lg)/SQRT(1.-thyfin(iiig,k,lg)**2) ! N0
                        wadd(iiig,k,lg)=
     1 padd(iiig,k,lg)*thyfin(iiig,k,lg)  
c
                        Nnp_str=Nnp_str+padd(iiig,k,lg)
C
C	                      Volodya test Streak
C
                        ener_str=ener_str+eadd(iiig,k,lg)
c      ener_str3=ener_str3+efin(iiig,k,lg)*chyz(iiig,lg)*dx*dy*dz
                        zmom_str=zmom_str+zmomadd(iiig,k,lg)
C
  260                 CONTINUE 
C
                    ELSE
C
                      
                      efin(iiig,k,lg)=efin0(iiig,lg) 
                      thyfin(iiig,k,lg)=thy(iiig,lg)
                      nfin(iiig,k,lg)=dens(iiig,lg)
                      press(iiig,k,lg)=Co2*efin(iiig,k,lg)-(1+Co2)*Bag
C
                      IF (press(iiig,k,lg).LT.0.)  press(iiig,k,lg)=0.
C
                      IF (for100 .NE. 1) GO TO 265  
C
c                              USE TO PREPARE FORT.100  
C   
                        eadd(iiig,k,lg)=
     1 efin0(iiig,lg)*((1.+Co2)*chyz(iiig,lg)**2-Co2)
                        zmomadd(iiig,k,lg)=
     1 efin0(iiig,lg)*(1.+Co2)*chyz(iiig,lg)**2*thy(iiig,lg)
                        wadd(iiig,k,lg)=
     1 dens(iiig,lg)*chyz(iiig,lg)*thy(iiig,lg)
                        padd(iiig,k,lg)=dens(iiig,lg)*chyz(iiig,lg)
                        Nnp_str=Nnp_str+padd(iiig,k,lg)
C
C	                      Volodya test Streak
C
                        ener_str=ener_str+eadd(iiig,k,lg)
c      ener_str3=ener_str3+efin(iiig,k,lg)*chyz(iiig,lg)*dx*dy*dz
                        zmom_str=zmom_str+zmomadd(iiig,k,lg)
C
  265                 CONTINUE   
C
                    ENDIF COND10
C
                  ENDIF COND9
C
  250           CONTINUE
C
              ENDIF COND8
C
            ENDIF COND7
C
  240     CONTINUE
  235   CONTINUE

  
        ELSE
C
          GO TO 140
          WRITE(*,132)
C
        ENDIF COND6
C
C     ######################################################################
C
        Nnp_str=Nnp_str*dx*dy*dz
        WRITE(*,*)'N total exp streaks ',Nnp_str
C
C	                          Volodya test Streak
C
        ener_str=ener_str*dx*dy*dz
        WRITE(*,*)'E total exp streaks ',ener_str
        write(*,*)'epn= ', ener_str/Nnp_str
        zmom_str=zmom_str*dx*dy*dz
        WRITE(*,*)'Pz total exp streaks ',zmom_str
        
      DO i=ismin,ismax
        DO j=1,Nceldasy
          DO k=tfkmin(i,j),tfkmax(i,j)

          if (i < 1 .or. k < 1) then
          write(*,*) i,k,tfkmin(i,j),tfkmax(i,j)
          stop 1254
          ENDIF
           
        IF (ISNAN(padd(i,k,j)) .EQV. .TRUE. .OR. ISNAN(wadd(i,k,j)) 
     1  .EQV. .TRUE. .OR. ISNAN(eadd(i,k,j)) .EQV. .TRUE. .OR.
     2  ISNAN(zmomadd(i,k,j)) .EQV. .TRUE. ) GO TO 998
     
        IF (ABS(padd(i,k,j)) .GT. HUGE(1.0D0) .OR. ABS(wadd(i,k,j)) .GT. 
     1   HUGE(1.0D0) .OR. ABS(eadd(i,k,j)) .GT. HUGE(1.0D0) .OR. 
     2   ABS(zmomadd(i,k,j)) .GT. HUGE(1.0D0)) THEN
     
         OPEN(UNIT=21,FILE='errors/Ne1.dat',STATUS='REPLACE')
C
                  WRITE(21,131) Ne1
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/Ne2.dat',STATUS='REPLACE')
C
                  WRITE(21,131) Ne2
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/l1.dat',STATUS='REPLACE')
C
                  WRITE(21,131) l1
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/l2.dat',STATUS='REPLACE')
C
                  WRITE(21,131) l2
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/sigma.dat',STATUS='REPLACE')
C
                  WRITE(21,131) sigma_s
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/tC1.dat',STATUS='REPLACE')
C
                  WRITE(21,131) tC1
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/tC2.dat',STATUS='REPLACE')
C
                  WRITE(21,131) tC2
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/efin.dat',STATUS='REPLACE')
C
                  WRITE(21,131) efin
C
         CLOSE(UNIT=21)
         
         OPEN(UNIT=21,FILE='errors/thyfin.dat',STATUS='REPLACE')
C
                  WRITE(21,131) thyfin
C
         CLOSE(UNIT=21)
         
         GO TO 998
         
        ENDIF     
        
          ENDDO
        ENDDO
      ENDDO
        
C
C     ######################################################################
C                               AVERAGE QUANTITIES
C     ######################################################################
C
        DO i=ismin,ismax
          DO j=1,Nceldasy
            DO k=tfkmin(i,j),tfkmax(i,j)
C
              paddtot(i,k,j)=(paddtot(i,k,j)*FLOAT(neventos-1)+
     1 padd(i,k,j))/FLOAT(neventos)
              waddtot(i,k,j)=(waddtot(i,k,j)*FLOAT(neventos-1)+
     1 wadd(i,k,j))/FLOAT(neventos)
              eaddtot(i,k,j)=(eaddtot(i,k,j)*FLOAT(neventos-1)+
     1 eadd(i,k,j))/FLOAT(neventos)
              zmomaddtot(i,k,j)=(zmomaddtot(i,k,j)*FLOAT(neventos-1)+
     1 zmomadd(i,k,j))/FLOAT(neventos)
C
C                 v,e,n de las leyes de conservación
C
              IF (paddtot(i,k,j) .EQ. 0. .AND. 
     1 eaddtot(i,k,j) .EQ. 0.) THEN
C
                thyfintot(i,k,j)=0.
                efintot(i,k,j)=0.
                presstot(i,k,j)=0.
                zmomtot(i,k,j)=0.
                nfintot(i,k,j)=0.
C
              ENDIF
C
              IF (paddtot(i,k,j) .EQ. 0. .AND. 
     1 eaddtot(i,k,j) .NE. 0.) THEN
C
                WRITE(*,*)'Averaging N=0,E>0'   
                thyfintot(i,k,j)= zmomaddtot(i,k,j)/eaddtot(i,k,j)! definición Landau ????????          
C      
                IF (thyfintot(i,k,j) .GE. 1.) THEN
C
                  WRITE(*,*) 'Velocidad>1',i,k,j,thyfintot(i,k,j)
                  thyfintot(i,k,j)=0.99
C
                ENDIF
C
                IF (thyfintot(i,k,j) .LE. -1.) THEN
C
                  WRITE(*,*) 'Velocidad<-1',i,k,j,thyfintot(i,k,j)
                  thyfintot(i,k,j)=-0.99
C
                ENDIF
C
                zmomtot(i,k,j)=zmomaddtot(i,k,j)   ! T0Z   ?????
                efintot(i,k,j)=eaddtot(i,k,j)/
     1 ((1+Co2)/(1.-thyfintot(i,k,j)**2)-Co2)
                presstot(i,k,j)=Co2*efintot(i,k,j)-(1+Co2)*Bag 
C
                IF (presstot(i,k,j).LT.0.)  presstot(i,k,j)=0.
C
                nfintot(i,k,j)=0.
C
              ENDIF
C
              IF (paddtot(i,k,j).NE. 0.) THEN
C
                thyfintot(i,k,j)= waddtot(i,k,j)/paddtot(i,k,j) ! definición Eckart 
C
                IF (thyfintot(i,k,j) .GE. 1.) THEN
C
                  WRITE(*,*) 'Velocidad>1',i,k,j,thyfintot(i,k,j)
                  thyfintot(i,k,j)=0.99
C
                ENDIF
C
                IF (thyfintot(i,k,j) .LE. -1.) THEN
C
                  WRITE(*,*) 'Velocidad<-1',i,k,j,thyfintot(i,k,j)
                  thyfintot(i,k,j)=-0.99
C
                ENDIF
C
                zmomtot(i,k,j)=zmomaddtot(i,k,j)   ! T0Z   ?????
                efintot(i,k,j)=eaddtot(i,k,j)/
     1 ((1+Co2)/(1.-thyfintot(i,k,j)**2)-Co2)
                presstot(i,k,j)=Co2*efintot(i,k,j)-(1+Co2)*Bag
C
                IF (presstot(i,k,j).LT.0.)  presstot(i,k,j)=0.
C
                nfintot(i,k,j)=paddtot(i,k,j)*SQRT(1.-
     1 thyfintot(i,k,j)**2)
C
              ENDIF
C
C               v,e,n de las leyes de conservación --END
C
C
            ENDDO
          ENDDO
        ENDDO   
C
      
        tnpnaverage2=(tnpnaverage2*FLOAT(neventos-1)+tpn2)/
     1 FLOAT(neventos)
        tPznpnaverage2=(tPznpnaverage2*FLOAT(neventos-1)+tPzpn2)/
     1 FLOAT(neventos)
        Nnp_str_aver=(Nnp_str_aver*FLOAT(neventos-1)+Nnp_str)/
     1 FLOAT(neventos)
C
C	                      Volodya test exp streaks
C
        ener_str_aver=(ener_str_aver*FLOAT(neventos-1)+ener_str)/
     1 FLOAT(neventos)
        zmom_str_aver=(zmom_str_aver*FLOAT(neventos-1)+zmom_str)/
     1 FLOAT(neventos)
       WRITE(*,*) Nnp_str_aver
       WRITE(*,*) ener_str_aver
       write(*,*)'epn (1438)= ', ener_str_aver/Nnp_str_aver

C
C     ######################################################################
C
C     ######################################################################
C                          SAVE A COPY OF THE RESULTS
C     ######################################################################
C
        IF (nteventos .LT. 5000) GO TO 330
C
        DO 335 j=1,100
C
          l=5000*j
C
          IF (neventos .EQ. l) THEN
          
C
            OPEN(UNIT=19,FILE='plots/inf.dat',STATUS='REPLACE')
C
                  WRITE(19,*) neventos
C
            CLOSE(UNIT=19) 
            
C
            OPEN(UNIT=20,FILE='plots/eaddtot.dat',STATUS='REPLACE')
C
                  WRITE(20,131) eaddtot
C
            CLOSE(UNIT=20) 
C
            OPEN(UNIT=21,FILE='plots/waddtot.dat',STATUS='REPLACE')
C
                  WRITE(21,131) waddtot
C
            CLOSE(UNIT=21)
            
            OPEN(UNIT=22,FILE='plots/zmomaddtot.dat',STATUS='REPLACE')
C
                  WRITE(22,131) zmomaddtot
C
            CLOSE(UNIT=22) 
C
            OPEN(UNIT=23,FILE='plots/paddtot.dat',STATUS='REPLACE')
C
                  WRITE(23,131) paddtot
C
            CLOSE(UNIT=23)
C
            
C
          ENDIF
C
  335   CONTINUE
C
  330   CONTINUE    
  
  
  
cc     ##################################################################

c      IF (vortin .NE. 1) GO TO 725
      
cc     ##################################################################
cc                         EVENT-BY-EVENT VORTICITY
cC     ##################################################################
cc
cc     Vfin is the velocity in z direction of each point of the discrete
cc     matter 
cc     Wy =1/2 (dxvz - dzvx)   ---->      vx = 0   ===>  Wy = 1/2 (dxvz)
cC     dx= partial derivative with repect to x
cc     Wy= vorticity in the reaction plane [xz]
cc     NOTATION 1 = previous point; 2 = next point; 0 = same point i.e: v12 =v-+
cc     gamma = 1/sqrt(1-v^2)
    
cC     ----------------------------------------------------------------------


c      vort_fin=0.
c      vort_rel_fin=0.

cc                             
c      DO iy=1,Nceldasy
c        DO iii=ismin,ismax
c            DO kkk=fkmin(iii,iy),fkmax(iii,iy)
cC                   
c               IF (efin(iii,kkk,iy) .NE. 0) THEN
cc                    Derivatives on x direction
c!
c!
c!      -----------------------------------------------------------------
cC 
c                       IF (efin(iii-1,kkk+1,iy) .EQ. 0.) THEN
cC                       
c                          v12=thyfin(iii,kkk+1,iy)
c                          gamma12=1./SQRT(1.-v12**2)
c                          w12 = 0.
cC                       
c                       ELSE
cC                       
c                          v12=thyfin(iii-1,kkk+1,iy)
c                          w12=1.
c                          gamma12= 1./SQRT(1.-v12**2)
cC                       
c                       ENDIF
ccC                       
c                       IF (efin(iii+1,kkk-1,iy) .EQ. 0.) THEN
cC                       
c                          v21=thyfin(iii,kkk-1,iy)
c                          gamma21=1./SQRT(1.-v21**2)
c                          w21=0.
cC                       
c                       ELSE
cC                       
c                          v21=thyfin(iii+1,kkk-1,iy)
c                          gamma21= 1./SQRT(1.-v21**2)
c                          w21=1.
cC                       
c                       ENDIF
ccC                       
c                       IF (efin(iii+1,kkk+1,iy) .EQ. 0.) THEN
cC                       
c                          v22 = thyfin(iii,kkk+1,iy)
c                          gamma22=1./SQRT(1.-v22**2)
c                          w22 = 0.
cC                       
c                       ELSE
cC                       
c                          v22 = thyfin(iii+1,kkk+1,iy)
c                          gamma22= 1./SQRT(1.-v22**2)
c                          w22 = 1.
cC                       
c                       ENDIF
cC                       
c                       IF (efin(iii-1,kkk-1,iy) .EQ. 0.) THEN
cC                       
c                           v11 = thyfin(iii,kkk-1,iy)
c                           gamma11=1./SQRT(1.-v11**2)
c                           w11 = 0.
cC                       
c                       ELSE
cC                       
c                           v11 = thyfin(iii-1,kkk-1,iy)
c                           gamma11= 1./SQRT(1.-v11**2)
c                           w11 = 1.
cC                       
c                       ENDIF
cC                       
cc                                  On z drection
cC 
cC 
c                       IF (efin(iii-1,kkk,iy) .EQ. 0.) THEN
cC                       
c                          v10 = thyfin(iii,kkk,iy)
c                          gamma10=1./SQRT(1.-v10**2)
c                          w10 = 0.
cC                       
c                       ELSE
cC                       
c                          v10 = thyfin(iii-1,kkk,iy)
c                          gamma10= 1./SQRT(1.-v10**2)
c                          w10=1.
cC                       
c                       ENDIF
cC                       
c                       IF (efin(iii+1,kkk,iy) .EQ. 0.) THEN
cC                       
c                          v20 = thyfin(iii,kkk,iy)
c                          gamma20=1./SQRT(1.-v20**2)
c                          w20 = 0.
cC                       
c                       ELSE
cC                       
c                          v20 = thyfin(iii+1,kkk,iy)
c                          gamma20= 1./SQRT(1.-v20**2)
c                          w20=1.
cC                       
c                       ENDIF
cC                       
cC                       
cc                              Boundary conditions
cC 
cC 
c                       IF (v22 .EQ. 0. .AND. w22 .EQ. 0.) THEN
cC                       
c                          v22 = thyfin(iii-1,kkk+1,iy)
c                          gamma22= 1./SQRT(1.-v22**2)
c                          w22=-1.
cC                       
c                       ENDIF
cC                       
c                       IF (v12 .EQ. 0. .AND. w12 .EQ. 0.) THEN
cC                       
c                          v12 = thyfin(iii+1,kkk+1,iy)
c                          gamma12= 1./SQRT(1.-v12**2)
c                          w12=-1.
cC                       
c                       ENDIF
cC                       
c                       IF (v21 .EQ. 0. .AND. w21 .EQ. 0.) THEN
cC                       
c                          v21 = thyfin(iii-1,kkk-1,iy)
c                          gamma21= 1./SQRT(1.-v21**2)
c                          w21=-1.
cC                       
c                       ENDIF
ccC                       
c                       IF (v11 .EQ. 0. .AND. w11 .EQ. 0.) THEN
cC                       
c                          v11 = thyfin(iii+1,kkk-1,iy)
c                          gamma11= 1./SQRT(1.-v11**2)
c                          w11=-1.
cC                       
c                       ENDIF
cC                       
cc                                   Final results
cC 
cC 
c                       q1 =w20 +w10
c                       q2 = w22 +w11 +w21 +w12
cC                       
c                       IF (q1 .EQ. 0.) q1=1.
c                       IF (q2 .EQ. 0.) q2=1.
cC                       
c                       v00=thyfin(iii,kkk,iy)
cC                       
c                       gamma00 = 1./SQRT(1.-v00**2)
cC                       
c                       vort_fin(iii,kkk,iy)=(-1./2.)*(((v20-v10)/
c     1 (q1*xsize))+((v22-v11+v21-v12)/(q2*xsize)))
cC     
c                       vort_rel_fin(iii,kkk,iy)=-(1./2.)*(gamma00*
c     1 ((v20-v10)/(q1*xsize)+(v22-v11+v21-v12)/(q2*xsize))+v00*
c     2 ((gamma20-gamma10)/(q1*xsize)+(gamma22-gamma11+gamma21-gamma12)/
c     3 (q2*xsize)))
ccC              
c               ENDIF         
c             ENDDO
c          ENDDO
c      ENDDO
      
c!  
c!     ------------------------------------------------------------------
c!                             TOTAL VORTICITY
c!
c      totclassvort=0.
c      totrelvort=0.
c!      
c      DO j=1,Nceldasy
c        DO i=1,Nceldasx
c          DO k=1,Nceldasz
c            totclassvort=totclassvort+vort_fin(i,k,j)
c            totrelvort=totrelvort+vort_rel_fin(i,k,j)
c          ENDDO
c        ENDDO
c      ENDDO    
      
      
c      mtotclassvort(neventos)=0.
c      mtotrelvort(neventos)=0.
c      mtotclassvort(neventos)=totclassvort
c      mtotrelvort(neventos)=totrelvort
      
c      IF (neventos .EQ. nteventos) THEN
      
c      OPEN(UNIT=20,FILE='plots/mtotclassvort.dat',STATUS='REPLACE')
cC
c            WRITE(20,131) mtotclassvort
cC
c      CLOSE(UNIT=20) 
      
c      OPEN(UNIT=25,FILE='plots/mtotrelvort.dat',STATUS='REPLACE')
cC
c            WRITE(25,131) mtotrelvort
cC
c      CLOSE(UNIT=25) 
      
      
c      DEALLOCATE(mtotclassvort)
c      DEALLOCATE(mtotrelvort)
      
c      ENDIF
      
      
      
cc!
c  725 CONTINUE
      
  
  
c
        IF (nteventos .GT. 1) GO TO 346
      
          DEALLOCATE(efin)
          DEALLOCATE(thyfin)
          DEALLOCATE(eadd)
          DEALLOCATE(zmomadd)
          DEALLOCATE(padd)
          DEALLOCATE(wadd)
          DEALLOCATE(ycorr)
          DEALLOCATE(press)
          DEALLOCATE(nfin)
          DEALLOCATE(Ne1)
          DEALLOCATE(Ne2)
          DEALLOCATE(Ne)
          DEALLOCATE(Ee)
          DEALLOCATE(Pze)
          DEALLOCATE(l1)
          DEALLOCATE(l2)
          DEALLOCATE(kz1min)
          DEALLOCATE(kz1max)
          DEALLOCATE(kz2min)
          DEALLOCATE(kz2max)
          DEALLOCATE(fkmin)
          DEALLOCATE(fkmax)
          DEALLOCATE(rho01)
          DEALLOCATE(rho02)
          DEALLOCATE(chy)
          DEALLOCATE(yf)
          DEALLOCATE(dens)
          DEALLOCATE(chyz)
C
  346   CONTINUE    
  
        WRITE(*,146) tnpnaverage2,tPznpnaverage2
        WRITE(*,*) '<N str>=',Nnp_str_aver
        WRITE(*,*) '<E str>=',ener_str_aver
        WRITE(*,*) '<Pz str>=',zmom_str_aver
        write(*,*)'epn (1766)= ', ener_str_aver/Nnp_str_aver
C      
C     ######################################################################
C
        CALL CPU_TIME(finish)
C      
        PRINT '("Time = ",f10.4," seconds.")',finish-start
C
c     ######################################################################
c

  999 CONTINUE
C
      
C     ------------------------- END SIMULATIONS ------------------------
c
c
      ALLOCATE(T00hydro(1:Nceldasx,1:495,1:Nceldasy))
      ALLOCATE(T0Zhydro(1:Nceldasx,1:495,1:Nceldasy))
      ALLOCATE(NZhydro(1:Nceldasx,1:495,1:Nceldasy))
      ALLOCATE(N0hydro(1:Nceldasx,1:495,1:Nceldasy))
      
      T00hydro=0.
      T0Zhydro=0.
      N0hydro=0.
      NZhydro=0.
      ehydro=0.
      vhydro=0.
      pzhydro=0.
      presshydro=0.
      nhydro=0.
                
c     Always =0 in our model    ! Volodya   
                pxhydro(i,k,j)=0.
                pyhydro(i,k,j)=0.
                u(i,k,j)=0.
                v(i,k,j)=0.
        
C     ######################################################################
C                      QUANTITIES IN HYDRO CELL
C     ######################################################################
C  
        Nnp=0.
        Enp=0.
C
C                                test VM 2
C
        IF (for100 .NE. 1) GO TO 345
c
c                     FORT.100
c
        DO 350 i=ismin,ismax
          DO 355 j=1,Nceldasy
c
            ll=1951
            kk=247 
c
  365       CONTINUE
c
            kk=kk+1
c
              DO 360 k=ll,ll+99

                T00hydro(i,kk,j) = T00hydro(i,kk,j) + eaddtot   (i,k,j) ! T00
                NZhydro (i,kk,j) = NZhydro (i,kk,j) + waddtot   (i,k,j) ! Nz
                T0Zhydro(i,kk,j) = T0Zhydro(i,kk,j) + zmomaddtot(i,k,j) ! T0Z
                N0hydro (i,kk,j) = N0hydro (i,kk,j) + paddtot   (i,k,j) ! N0
C
  360         CONTINUE
C
              Nnp = Nnp + N0hydro (i,kk,j)
              Enp = Enp + T00hydro(i,kk,j)
C
              ll=ll+100
C
C                  v,e,n de las leyes de conservación
C
              IF (N0hydro(i,kk,j) .EQ. 0. .AND. 
     1 T00hydro(i,kk,j) .EQ. 0.) THEN
C
                vhydro    (i,kk,j) =0.
                ehydro    (i,kk,j) =0.
                presshydro(i,kk,j) =0.
                pzhydro   (i,kk,j) =0.
                nhydro    (i,kk,j) =0.
C
              ENDIF
C
              IF (N0hydro(i,kk,j) .EQ. 0. .AND. 
     1 T00hydro(i,kk,j) .NE. 0.) THEN
C
                WRITE(*,*)'hydro N=0,E>0'   

                vhydro(i,kk,j) = T0Zhydro(i,kk,j)/T00hydro(i,kk,j)! definición Landau  ??????????
C      
                IF (vhydro(i,kk,j) .GE. 1.) THEN
C
                  WRITE(*,*) 'Hydro Velocidad>1',i,kk,j,vhydro(i,kk,j)

                  vhydro(i,kk,j) = 0.99
C
                ENDIF
C
                IF (vhydro(i,kk,j) .LE. -1.) THEN
C
                  WRITE(*,*) 'Hydro Velocidad<-1',i,kk,j,vhydro(i,kk,j)

                  vhydro(i,kk,j) = -0.99
C
                ENDIF
C
                pzhydro   (i,kk,j) = T0Zhydro(i,kk,j)/100.    ! T0Z/V   ?????
                ehydro    (i,kk,j) = T00hydro(i,kk,j)/100./
     1 ((1+Co2)/(1.-vhydro(i,kk,j)**2)-Co2)
                presshydro(i,kk,j) = Co2*ehydro(i,kk,j) - (1+Co2)*Bag
C
                IF (presshydro(i,kk,j).LT.0.)  presshydro(i,kk,j) = 0.
C
                nhydro(i,kk,j) = 0.
C
              ENDIF
C
              IF (N0hydro(i,kk,j) .NE. 0.) THEN
C
                vhydro(i,kk,j)=NZhydro(i,kk,j)/N0hydro(i,kk,j) ! definición Eckart 
C
                IF (vhydro(i,kk,j) .GE. 1.) THEN
C
                  WRITE(*,*) 'Hydro Velocidad>1',i,kk,j,vhydro(i,kk,j)
                  vhydro(i,kk,j)=0.99
C
                ENDIF
C
                IF (vhydro(i,kk,j) .LE. -1.) THEN
C
                  WRITE(*,*) 'Hydro Velocidad<-1',i,kk,j,vhydro(i,kk,j)
                  vhydro(i,kk,j)=-0.99
C
                ENDIF
C
                pzhydro(i,kk,j)=T0Zhydro(i,kk,j)/100.    ! T0Z/V   ?????
                ehydro(i,kk,j)=T00hydro(i,kk,j)/100./
     1 ((1+Co2)/(1.-vhydro(i,kk,j)**2)-Co2)
                presshydro(i,kk,j)=Co2*ehydro(i,kk,j)-(1+Co2)*Bag
C
                IF (presshydro(i,kk,j).LT.0.)  presshydro(i,kk,j)=0.
C
                nhydro(i,kk,j)=N0hydro(i,kk,j)*
     1 SQRT(1.-vhydro(i,kk,j)**2)/100.
C
              ENDIF
C
C                    v,e,n de las leyes de conservación --END
C
              IF (kk .LE. 495 .AND. ll .LE. Nceldasz - 100) GO TO 365
C
                jj=1950
                kk=248
C
  375           CONTINUE
C
                kk=kk-1
C
                DO 370 k=jj-99,jj
C
                  T00hydro(i,kk,j)=T00hydro(i,kk,j)+eaddtot(i,k,j)  ! T00
                  NZhydro(i,kk,j)=NZhydro(i,kk,j)+waddtot(i,k,j) ! NZ
                  T0Zhydro(i,kk,j)=T0Zhydro(i,kk,j)+zmomaddtot(i,k,j) ! T0Z
                  N0hydro(i,kk,j)=N0hydro(i,kk,j)+paddtot(i,k,j) ! N0
C
  370           CONTINUE 
C
                Nnp=Nnp+N0hydro(i,kk,j)
                Enp=Enp+T00hydro(i,kk,j)
C
C                                test VM 2
C
                jj=jj-100
C
C                    v,e,n de las leyes de conservación
C
                IF (N0hydro(i,kk,j) .EQ. 0. .AND. 
     1 T00hydro(i,kk,j) .EQ. 0.) THEN
C
                  vhydro(i,kk,j)=0.
                  ehydro(i,kk,j)=0.
                  presshydro(i,kk,j)=0.
                  pzhydro(i,kk,j)=0.
                  nhydro(i,kk,j)=0.
C
                ENDIF
C
                IF (N0hydro(i,kk,j) .EQ. 0. .AND. 
     1 T00hydro(i,kk,j) .NE. 0.) THEN
C
                  WRITE(*,*)'hydro N=0,E>0'   
                vhydro(i,kk,j)= T0Zhydro(i,kk,j)/T00hydro(i,kk,j)! definición Landau  ??????????     
C
                  IF (vhydro(i,kk,j) .GE. 1.) THEN
C
                    WRITE(*,*) 'Hydro Velocidad>1',i,kk,j,vhydro(i,kk,j)
                    vhydro(i,kk,j)=0.99
C
                  ENDIF
C
                  IF (vhydro(i,kk,j) .LE. -1.) THEN
C
                  WRITE(*,*) 'Hydro Velocidad<-1',i,kk,j,vhydro(i,kk,j)
                  vhydro(i,kk,j)=-0.99
C
                  ENDIF
C
                  pzhydro(i,kk,j)=T0Zhydro(i,kk,j)/100.    ! T0Z/V   ?????
                  ehydro(i,kk,j)=T00hydro(i,kk,j)/100./
     1 ((1+Co2)/(1.-vhydro(i,kk,j)**2)-Co2)   
                  presshydro(i,kk,j)=Co2*ehydro(i,kk,j)-(1+Co2)*Bag
C
                  IF (presshydro(i,kk,j).LT.0.)  presshydro(i,kk,j)=0.
C
                  nhydro(i,kk,j)=0.
C
                ENDIF
C
                IF (N0hydro(i,kk,j) .NE. 0.) THEN
C
                  vhydro(i,kk,j)=NZhydro(i,kk,j)/N0hydro(i,kk,j) ! definición Eckart 
C
                  IF (vhydro(i,kk,j) .GE. 1.) THEN
C
                    WRITE(*,*) 'Hydro Velocidad>1',i,kk,j,vhydro(i,kk,j)
                    vhydro(i,kk,j)=0.99
                  ENDIF
C          
                  IF (vhydro(i,kk,j) .LE. -1.) THEN
C          
                   WRITE(*,*) 'Hydro Velocidad<-1',i,kk,j,vhydro(i,kk,j)
                    vhydro(i,kk,j)=-0.99
C             
                 ENDIF
C
                 pzhydro(i,kk,j)=T0Zhydro(i,kk,j)/100.    ! T0Z/V   ?????
                 ehydro(i,kk,j)=T00hydro(i,kk,j)/100./
     1 ((1+Co2)/(1.-vhydro(i,kk,j)**2)-Co2)
                  presshydro(i,kk,j)=Co2*ehydro(i,kk,j)-(1+Co2)*Bag
C
                  IF (presshydro(i,kk,j).LT.0.)  presshydro(i,kk,j)=0.
C               
                  nhydro(i,kk,j)=N0hydro(i,kk,j)*
     1 SQRT(1.-vhydro(i,kk,j)**2)/100.
C
                ENDIF
C
C                  v,e,n de las leyes de conservación --END
c  
                IF (kk .GE. 1 .AND. jj .GE. 201) GO TO 375
C               
  355     CONTINUE    
  350   CONTINUE
c       
         Nnp=Nnp*dx*dy*dz  !*100.  is not need, ya esta sumado!!!!!  
         Enp=Enp*dx*dy*dz
c
        WRITE(*,*) 'N in hydro',Nnp
        WRITE(*,*) 'E in hydro',Enp
        write(*,*)'epn (2029)= ', Enp/Nnp

        
        
  345   CONTINUE



C
C
C
C     ######################################################################
C                        PRINT HYDRO QUANTITIES 
C     ######################################################################
C   
      ALLOCATE(epaper(1:Nceldasx,1:Nceldasz,1:Nceldasy))
      ALLOCATE(vpaper(1:Nceldasx,1:Nceldasz,1:Nceldasy))
      ALLOCATE(pzpaper(1:Nceldasx,1:Nceldasz,1:Nceldasy))
      ALLOCATE(npaper(1:Nceldasx,1:Nceldasz,1:Nceldasy))
      ALLOCATE(presspaper(1:Nceldasx,1:Nceldasz,1:Nceldasy))
            
      vpaper=0.
      epaper=0.
      presspaper=0.
      pzpaper=0.
      npaper=0.
      Nnp=0.
      Enp=0.
            
        DO 356 i=ismin,ismax
          DO 357 j=1,Nceldasy
            DO 358 k=1,Nceldasz
C
C                  v,e,n de las leyes de conservación
C
              IF (paddtot(i,k,j) .EQ. 0. .AND. 
     1 eaddtot(i,k,j) .EQ. 0.) THEN
C
                vpaper(i,k,j)=0.
                epaper(i,k,j)=0.
                pzpaper(i,k,j)=0.
                npaper(i,k,j)=0.
                presspaper(i,k,j)=0.
                Nnp=Nnp+npaper(i,k,j)
                Enp=Enp+epaper(i,k,j)
C
              ENDIF
cC
              IF (paddtot(i,k,j) .EQ. 0. .AND. 
     1 eaddtot(i,k,j) .NE. 0.) THEN
cC
                vpaper(i,k,j)= zmomaddtot(i,k,j)/eaddtot(i,k,j)! definición Landau  ??????????
C      
                IF (vpaper(i,k,j) .GE. 1.) THEN
C
                  vpaper(i,k,j)=0.99
C
                ENDIF
C
                IF (vpaper(i,k,j) .LE. -1.) THEN
C
                  vpaper(i,k,j)=-0.99
C
                ENDIF
cC
                pzpaper(i,k,j)=zmomaddtot(i,k,j)    ! T0Z/V   ?????
c                epaper(i,k,j)=eaddtot(i,k,j)/
c     1 ((1+Co2)/(1.-vpaper(i,k,j)**2)-Co2)
                epaper(i,k,j)=eaddtot(i,k,j)
                presspaper(i,k,j)=Co2*epaper(i,k,j)-(1+Co2)*Bag
C
                IF (presspaper(i,k,j).LT.0.)  presspaper(i,k,j)=0.
cC
                npaper(i,k,j)=0.
cC
              ENDIF
cC
                Nnp=Nnp+npaper(i,k,j)
                Enp=Enp+epaper(i,k,j)
               
              IF (paddtot(i,k,j) .NE. 0.) THEN
C
                vpaper(i,k,j)=waddtot(i,k,j)/paddtot(i,k,j) ! definición Eckart 
C
                IF (vpaper(i,k,j) .GE. 1.) THEN
C
                  vpaper(i,k,j)=0.99
C
                ENDIF
C
                IF (vpaper(i,k,j) .LE. -1.) THEN
C
                  vpaper(i,k,j)=-0.99
C
                ENDIF
cC
                pzpaper(i,k,j)=zmomaddtot(i,k,j)    ! T0Z/V   ?????
c                epaper(i,k,j)=eaddtot(i,k,j)/
c     1 ((1.+Co2)/(1.-vpaper(i,k,j)**2)-Co2)
                epaper(i,k,j)=eaddtot(i,k,j)
                presspaper(i,k,j)=Co2*epaper(i,k,j)-(1+Co2)*Bag
C
                IF (presspaper(i,k,j).LT.0.)  presspaper(i,k,j)=0.
cC
c                npaper(i,k,j)=paddtot(i,k,j)*
c     1 SQRT(1.-vpaper(i,k,j)**2)
                npaper(i,k,j)=paddtot(i,k,j)
cC
              ENDIF
              Nnp=Nnp+npaper(i,k,j)
              Enp=Enp+epaper(i,k,j)
cC  
  358       CONTINUE             
  357     CONTINUE    
  356   CONTINUE

        WRITE(*,*) 'N paper',Nnp*dx*dy*dz
        WRITE(*,*) 'E paper',Enp*dx*dy*dz
        write(*,*)'epn (2147)= ', Enp/Nnp

c      
c
c
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c




c

c     ##################################################################

      IF (vortin .NE. 1) GO TO 270
      
c     ##################################################################
c                            AVERAGE VORTICITY
C     ##################################################################
c
c     Vfin is the velocity in z direction of each point of the discrete
c     matter 
c     Wy =1/2 (dxvz - dzvx)   ---->      vx = 0   ===>  Wy = 1/2 (dxvz)
C     dx= partial derivative with repect to x
c     Wy= vorticity in the reaction plane [xz]
c     NOTATION 1 = previous point; 2 = next point; 0 = same point i.e: v12 =v-+
c     gamma = 1/sqrt(1-v^2)
    
C     ----------------------------------------------------------------------


      

c                             
      DO iy=1,Nceldasy
        DO iii=1,Nceldasx
            DO kkk=1,Nceldasz
               vort_fin(iii,kkk,iy)=0.
               vort_rel_fin(iii,kkk,iy)=0.
C                   
               IF (efintot(iii,kkk,iy) .NE. 0) THEN
c                    Derivatives on x direction
!
!
!      -----------------------------------------------------------------
C 
                       IF (efintot(iii-1,kkk+1,iy) .EQ. 0.) THEN
C                       
                          v12=thyfintot(iii,kkk+1,iy)
                          gamma12=1./SQRT(1.-v12**2)
                          w12 = 0.
C                       
                       ELSE
C                       
                          v12=thyfintot(iii-1,kkk+1,iy)
                          w12=1.
                          gamma12= 1./SQRT(1.-v12**2)
C                       
                       ENDIF
cC                       
                       IF (efintot(iii+1,kkk-1,iy) .EQ. 0.) THEN
C                       
                          v21=thyfintot(iii,kkk-1,iy)
                          gamma21=1./SQRT(1.-v21**2)
                          w21=0.
C                       
                       ELSE
C                       
                          v21=thyfintot(iii+1,kkk-1,iy)
                          gamma21= 1./SQRT(1.-v21**2)
                          w21=1.
C                       
                       ENDIF
cC                       
                       IF (efintot(iii+1,kkk+1,iy) .EQ. 0.) THEN
C                       
                          v22 = thyfintot(iii,kkk+1,iy)
                          gamma22=1./SQRT(1.-v22**2)
                          w22 = 0.
C                       
                       ELSE
C                       
                          v22 = thyfintot(iii+1,kkk+1,iy)
                          gamma22= 1./SQRT(1.-v22**2)
                          w22 = 1.
C                       
                       ENDIF
C                       
                       IF (efintot(iii-1,kkk-1,iy) .EQ. 0.) THEN
C                       
                           v11 = thyfintot(iii,kkk-1,iy)
                           gamma11=1./SQRT(1.-v11**2)
                           w11 = 0.
C                       
                       ELSE
C                       
                           v11 = thyfintot(iii-1,kkk-1,iy)
                           gamma11= 1./SQRT(1.-v11**2)
                           w11 = 1.
C                       
                       ENDIF
C                       
c                                  On z drection
C 
C 
                       IF (efintot(iii-1,kkk,iy) .EQ. 0.) THEN
C                       
                          v10 = thyfintot(iii,kkk,iy)
                          gamma10=1./SQRT(1.-v10**2)
                          w10 = 0.
C                       
                       ELSE
C                       
                          v10 = thyfintot(iii-1,kkk,iy)
                          gamma10= 1./SQRT(1.-v10**2)
                          w10=1.
C                       
                       ENDIF
C                       
                       IF (efintot(iii+1,kkk,iy) .EQ. 0.) THEN
C                       
                          v20 = thyfintot(iii,kkk,iy)
                          gamma20=1./SQRT(1.-v20**2)
                          w20 = 0.
C                       
                       ELSE
C                       
                          v20 = thyfintot(iii+1,kkk,iy)
                          gamma20= 1./SQRT(1.-v20**2)
                          w20=1.
C                       
                       ENDIF
C                       
C                       
c                              Boundary conditions
C 
C 
                       IF (v22 .EQ. 0. .AND. w22 .EQ. 0.) THEN
C                       
                          v22 = thyfintot(iii-1,kkk+1,iy)
                          gamma22= 1./SQRT(1.-v22**2)
                          w22=-1.
C                       
                       ENDIF
C                       
                       IF (v12 .EQ. 0. .AND. w12 .EQ. 0.) THEN
C                       
                          v12 = thyfintot(iii+1,kkk+1,iy)
                          gamma12= 1./SQRT(1.-v12**2)
                          w12=-1.
C                       
                       ENDIF
C                       
                       IF (v21 .EQ. 0. .AND. w21 .EQ. 0.) THEN
C                       
                          v21 = thyfintot(iii-1,kkk-1,iy)
                          gamma21= 1./SQRT(1.-v21**2)
                          w21=-1.
C                       
                       ENDIF
cC                       
                       IF (v11 .EQ. 0. .AND. w11 .EQ. 0.) THEN
C                       
                          v11 = thyfintot(iii+1,kkk-1,iy)
                          gamma11= 1./SQRT(1.-v11**2)
                          w11=-1.
C                       
                       ENDIF
C                       
c                                   Final results
C 
C 
                       q1 =w20 +w10
                       q2 = w22 +w11 +w21 +w12
C                       
                       IF (q1 .EQ. 0.) q1=1.
                       IF (q2 .EQ. 0.) q2=1.
C                       
                       v00=thyfintot(iii,kkk,iy)
C                       
                       gamma00 = 1./SQRT(1.-v00**2)
C                       
                       vort_fin(iii,kkk,iy)=(-1./2.)*(((v20-v10)/
     1 (q1*xsize))+((v22-v11+v21-v12)/(q2*xsize)))
C     
                       vort_rel_fin(iii,kkk,iy)=-(1./2.)*(gamma00*
     1 ((v20-v10)/(q1*xsize)+(v22-v11+v21-v12)/(q2*xsize))+v00*
     2 ((gamma20-gamma10)/(q1*xsize)+(gamma22-gamma11+gamma21-gamma12)/
     3 (q2*xsize)))
cC              
               ENDIF         
             ENDDO
          ENDDO
      ENDDO
      
!  
!     ------------------------------------------------------------------
!                             TOTAL VORTICITY
!
      totclassvort=0.
      totrelvort=0.
!      
      DO j=1,Nceldasy
        DO i=1,Nceldasx
          DO k=1,Nceldasz
            totclassvort=totclassvort+vort_fin(i,k,j)
            totrelvort=totrelvort+vort_rel_fin(i,k,j)
          ENDDO
        ENDDO
      ENDDO    
      
      WRITE(*,*) 'totclassvort', totclassvort
      WRITE(*,*) 'totrelvort', totrelvort
      
      DEALLOCATE(vort_fin)
      DEALLOCATE(vort_rel_fin)
      
c!
  270 CONTINUE
    
c
c
C 
C     ######################################################################
C                        INPUT FOR HYDRO (FORT.100)
c     ######################################################################
c      
c      ---------------------------------------------------------------------
C            We calculate the number of lines in fort.100 file
c
      i7=0
      eehydro=0.
      nnhydro=0.
      T00hydro=0.
      N0hydro=0.
C
      DO i=1,Nceldasx
        DO j=1,Nceldasy
          DO k=228,268
C
            IF (ehydro(i,k,j) .NE. 0. .OR. nhydro(i,k,j) .NE. 0. .OR. 
     1 presshydro(i,k,j) .NE. 0. .OR. vhydro(i,k,j) .NE. 0. .OR. 
     2 u(i,k,j) .NE. 0. .OR. v(i,k,j) .NE. 0. .OR. pzhydro(i,k,j) 
     3 .NE. 0. .OR. pxhydro(i,k,j) .NE. 0. .OR. pyhydro(i,k,j) 
     4 .NE. 0.) THEN
     
              i7=i7+1
       T00hydro(i,k,j)=ehydro(i,k,j)*((1+Co2)/(1.-vhydro(i,k,j)**2)
     1 -Co2)
       N0hydro(i,k,j)=nhydro(i,k,j)/SQRT(1.-vhydro(i,k,j)**2)
       
       eehydro=eehydro+T00hydro(i,k,j)*dx*dy*dz*100.
       nnhydro=nnhydro+N0hydro(i,k,j)*dx*dy*dz*100

            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
c      
      N_lines=i7
      eednn=eehydro/nnhydro
C
c     ----------------------------------------------------------------------
c
      OPEN(UNIT=100,FILE='printforinhydro/fort.100',STATUS='REPLACE')
      OPEN(UNIT=111,FILE='printforinhydro/fort_100.dat',STATUS='REPLACE'
     1 )
c
      IF (nteventos .GT. 100) THEN
      WRITE (100,311) nteventos,nnhydro,eehydro,eednn ! Dan's suggestion
      ELSE
      WRITE (100,312) nteventos,nnhydro,eehydro,eednn ! Dan's suggestion
      ENDIF
      WRITE (100,148) Tfin
      WRITE (100,136) epsilon00,b0,A1,A2,tnc
      WRITE (100,137) rap00
      WRITE (100,138) xsize,ysize,zsize*100
      WRITE (100,139) i7
      WRITE (100,141)
c
      Nnptot=0.
      Enptot=0.
C                                  test VM
      DO i=1,Nceldasx
        DO j=1,Nceldasy
          DO k=228,268
c
            IF (T00hydro(i,k,j) .NE. 0. .OR. N0hydro(i,k,j) .NE. 0. .OR. 
     1 presshydro(i,k,j) .NE. 0. .OR. vhydro(i,k,j) .NE. 0. .OR. 
     2 u(i,k,j) .NE. 0. .OR. v(i,k,j) .NE. 0. .OR. pzhydro(i,k,j) 
     3 .NE. 0. .OR. pxhydro(i,k,j) .NE. 0. .OR. pyhydro(i,k,j) 
     4 .NE. 0.) THEN
c
              WRITE(100,142) i,j,k,T00hydro(i,k,j),
     1            N0hydro(i,k,j),presshydro(i,k,j),
     2            vhydro(i,k,j),u(i,k,j),v(i,k,j),pzhydro(i,k,j),
     3            pxhydro(i,k,j),pyhydro(i,k,j)
              WRITE(111,142) i,j,k,T00hydro(i,k,j),
     1            N0hydro(i,k,j),presshydro(i,k,j),
     2            vhydro(i,k,j),u(i,k,j),v(i,k,j),pzhydro(i,k,j),
     3            pxhydro(i,k,j),pyhydro(i,k,j)
 
                  Nnptot=Nnptot+N0hydro(i,k,j)*dx*dy*dz*100.
                  Enptot=Enptot+T00hydro(i,k,j)*dx*dy*dz*100.
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(UNIT=111)
      CLOSE(UNIT=100)
C
      N_lines=i7*dx*dy*dz*100.
      WRITE(*,*) 'Nnptot in fort 100',Nnptot
      WRITE(*,*) 'Enptot in fort 100',Enptot
      WRITE(*,*) 'Reaction volume fort 100',N_lines
      write(*,*)'epn (2475)= ', Enp/Nnp
C   
c     ###################### END INPUT FOR HYDRO ###########################
c                  
      Nnptot=0.
      Enptot=0.
      i7=0
      N_lines=0
      OPEN(UNIT=115,FILE='printforinhydro/fort_115.dat',STATUS='REPLACE'
     1 )
C                                  
      DO i=1,Nceldasx
        DO j=1,Nceldasy
          DO k=1,Nceldasz
c
            IF (epaper(i,k,j) .NE. 0. .OR. npaper(i,k,j) .NE. 0. .OR. 
     1 presspaper(i,k,j) .NE. 0. .OR. vpaper(i,k,j) .NE. 0. .OR. 
     2 pzpaper(i,k,j) .NE. 0.) THEN
     
              i7=i7+1
c
              WRITE(115,213) i,j,k,epaper(i,k,j),
     1            npaper(i,k,j),presspaper(i,k,j),
     2            vpaper(i,k,j),pzpaper(i,k,j)
     
            Nnptot=Nnptot+npaper(i,k,j)*dx*dy*dz
            Enptot=Enptot+epaper(i,k,j)*dx*dy*dz
C
            ENDIF
C
          ENDDO
        ENDDO
      ENDDO
C
      CLOSE(UNIT=115)
      
      N_lines=i7*dx*dy*dz
      WRITE(*,*) 'Nnptot in fort 115',Nnptot
      WRITE(*,*) 'Enptot in fort 115',Enptot
      write(*,*)'epn (2516)= ', Enptot/Nnptot
      WRITE(*,*) 'Reaction Volume fort 115',N_lines
C
c     ##################################################################

      DEALLOCATE(efintot)
      DEALLOCATE(nfintot)
      DEALLOCATE(zmomtot)
      DEALLOCATE(thyfintot)
      DEALLOCATE(eaddtot)
      DEALLOCATE(epaper)
      DEALLOCATE(vpaper)
      DEALLOCATE(pzpaper)
      DEALLOCATE(npaper)
      DEALLOCATE(presspaper)
      DEALLOCATE(paddtot)
      DEALLOCATE(zmomaddtot)
      DEALLOCATE(waddtot)
      DEALLOCATE(presstot)
      DEALLOCATE(T00hydro)
      DEALLOCATE(T0Zhydro)
      DEALLOCATE(NZhydro)
      DEALLOCATE(N0hydro)
      DEALLOCATE(tfkmin)
      DEALLOCATE(tfkmax)
C
c
c
C     ######################################################################
C                                  FORMATS
C     ######################################################################
C
  106 FORMAT('################ SIZE OF CELLS (fm) ################')   
  107 FORMAT('dx= ',F10.6,5x,'dy= ',F10.6,5x,'dz= ',F10.6,//) 
  108 FORMAT('# INITIAL POSITION OF THE CENTRE OF EACH NUCLEUS (fm) #')
  109 FORMAT('x1= ',F6.3,5X,'y1= ',F6.3,5X,'z1= ',F6.3,/)
  111 FORMAT('x2= ',F6.3,5X,'y2= ',F6.3,5X,'z2= ',F6.3,//)
  112 FORMAT('############# TRANSVERSE PLANE DIMENSIONS #############')
  113 FORMAT('xmin= ',F10.6,5x,'xmax= ',F10.6,/,'ymin= ',F10.6,5x,
     1 'ymax= ',F10.6,//)    
  114 FORMAT('########## NUMBER OF NODES IN EACH DIRECTION ###########')
  116 FORMAT('nx= ',I4,5X,'ny= ',I4,5X,'nz= ',I6,//)   
  117 FORMAT('########## NUMBER OF CELLS IN EACH DIRECTION ###########')
  118 FORMAT('Nceldasx= ',i3,5X,'Nceldasy= ',i3,5X,'Nceldasz= ',i4,//)
  119 FORMAT('## LIMITS OF THE TRANSVERSE PLANE IN THE X-DIRECTION ###')
  121 FORMAT('ismin= ',i3,5X,'ismax= ',i3,/)
  122 FORMAT('xsmin= ',F7.3,' fm ',5X,'xsmax= ',F7.3,' fm ') 
  123 FORMAT('No. OF PARTICIPANT NUCLEONS METH.1= ',F10.5,/
     1 'ENERGY PARTICIPANT NUCLEONS METH.1= ',F14.5,2x,'GeV',/
     2 'MOMENTUM PARTICIPANT NUCLEONS METH.1= ',F14.5,2x,'GeV')    
  124 FORMAT('TOTAL BARYON CHARGE= ',F10.5,/
     1 'TOTAL ENERGY= ',F14.5,2X,'GeV',/
     2 'TOTAL MOMENTUM= ',F14.5,2X,'GeV')
  126 FORMAT('MMMT:', F16.10)
  127 FORMAT('Problem: MINT_lim < MMMT')
  128 FORMAT('error: Give me Tfin >= MMMT')
  129 FORMAT('Okey') 
  131 FORMAT(E11.4)
  132 FORMAT('Tfin>MINT - THIS EVENT IS NOT COUNTED')
  133 FORMAT('Number of simulations:',   i10)
C	formats for fort.100	
  311 FORMAT ('# Au+Au simulation GESRM smooth Nevents=',I8,
     1 ' N = ',E11.6,' E = ',E11.6,' GeV',' E/N = ',E11.6,' GeV')
  312 FORMAT ('# Au+Au simulation GESRM fluct Nevents= ',I8,
     1 ' N = ',E11.6,' E = ',E11.6,' GeV',' E/N = ',E11.6,' GeV')
  136 FORMAT('# Au+Au, epsilon0 [GeV/A]= ',e8.2,3x,'b= ',e8.2,' r1+r2''
     1 - Ap, At = ',i3,3x,i3,' ,  No. of cells in tr. diameter, ktnc= '
     2 ,i2)
  137 FORMAT('# beam y0= ',e11.5)
C  137 FORMAT('# beam y0= ',e11.5,3x,'rho0 [fm-3]= ',e11.5,3x,' e0 [GeV/
C     1fm3] = ',e11.5)
  138 FORMAT('# Cell size ','dx= ',E12.6,2x,'dy= ',E12.6,2x,'dz= ',E12.
     1 6,' [fm]')
  139 FORMAT('# All quantities are given in CALCULATIONAL frame.',' No.
     1 of lines in table= ',i7) 
  141 FORMAT('#  i   j  k       T00         ',' N0           p
     1vz          vx          vy','         zmom        xmom         ymo
     2m   '/'#               [GeV/fm3]','     [1/fm3]    [GeV/fm3]    ',
     3' [c]         [c]         [c]','      [GeV/fm3]   [GeV/fm3]    ',
     4'[GeV/fm3]')
  142 FORMAT(3i4,' ',9(E11.4,' '))
C	formats for fort.100	
  213 FORMAT(3i6,' ',9(E11.4,' '))
  144 FORMAT('No. OF PARTICIPANT NUCLEONS= ',F10.5,/
     2 'ENERGY PARTICIPANT NUCLEONS= ',F14.5,' GeV'/
     3 'MOMENTUM PARTICIPANT NUCLEONS= ',F14.5,' GeV')
  146 FORMAT('<PARTICIPANT NUCLEONS METH2>= ',F12.6,/
     1 '<MOMENTUM PARTICIPANT NUCLEONS METH2>= ',F12.6)
c  147 FORMAT('<N in hydro cells> = ',F12.6)
  148 FORMAT('# Tfin [fm/c]= ',E9.3)
c  151 FORMAT('# Initial energy per nucleon: epsilon00= ',F8.3,'[GeV]')
  151 FORMAT(48X,F10.3)
C  152 FORMAT('# Beam parameter: b0= ',F5.3)
  152 FORMAT(23X,F5.3)
C  153 FORMAT('# Atomic mass target: A1= ',I5)
  153 FORMAT(27X,I3)
C  154 FORMAT('# Atomic mass projectile: A2= ',I5)
  154 FORMAT(31X,I3)
C  156 FORMAT('# Target radius: R1= ',F6.3,'[fm]')
  156 FORMAT(27X,F6.3)
C  157 FORMAT('# Projectile radius: R2= ',F6.3,'[fm]')
  157 FORMAT(31X,F6.3)
C  158 FORMAT('# Skin depth: aws1= ',F6.3,'[fm]')
  158 FORMAT(26X,F6.3)
C  159 FORMAT('# Skin depth: aws2= ',F6.3,'[fm]')
  159 FORMAT(26X,F6.3)
C  161 FORMAT('# Nucleon mass: Mn= ',F6.3,'[GeV]')
  161 FORMAT(27X,F6.3)
C  171 FORMAT('# Nucleon radius: rn= ',F6.3,'[GeV]')
  171 FORMAT(28X,F6.3)
C  162 FORMAT('# String tension parameter: A= ',F6.3)
  162 FORMAT(32X,F6.3)
C  163 FORMAT('# Normal nuclear density: rho00= ',F6.3,'[fm-3]')
  163 FORMAT(41X,F6.3)
C  164 FORMAT('# tnc= ',I4)
  164 FORMAT(8X,I2)
C  166 FORMAT('# Number of simulations: neventos= ',I6)
  166 FORMAT(36X,I6)
C  167 FORMAT('# Vorticity calculations: vortin= ',I4)
  167 FORMAT(35X,I4)
C  168 FORMAT('# Print fort.100: for100= ',I4)
  168 FORMAT(27X,I4)
C  169 FORMAT('# Check calculations: check= ',I4)
  169 FORMAT(30X,I4)
  211 FORMAT(10X,F4.3)
  173 FORMAT('#################### INITIAL PARAMETERS #################
     1 ##')
  174 FORMAT ('#Initial energy density e00= ',F8.4,' [GeV/fm3]')
  176 FORMAT ('#Gamma factor gsmma0= ',F8.4)
  177 FORMAT ('#Initial velocity v0=',F8.4,' [c]')
  178 FORMAT ('#Initial rapidity rap00= ',F8.4)
  179 FORMAT ('#Impact parameter b= ',F8.4,' [fm]')
  181 FORMAT ('#######################################################')
  182 FORMAT ('################ INPUT FILE PARAMETERS ###############')
  183 FORMAT ('#Initial energy per nucleon: epsilon00= ',F10.3,' [GeV]')
  184 FORMAT ('#Beam parameter: b0= ',F8.4)
  186 FORMAT ('#Atomic mass target: A1= ',I4)
  187 FORMAT ('#Atomic mass projectile: A2= ',I4)
  188 FORMAT ('#Target radius: r1= ',F8.4, '[fm]')
  189 FORMAT ('#Projectile radius: r2= ',F8.4, '[fm]')
  192 FORMAT ('#WS skin depth parameter for target:aws1= ',F8.4, '[fm]')
  193 FORMAT ('#WS skin depth parameter for proj.: aws2= ',F8.4, '[fm]')
  194 FORMAT ('#Nucleon mass: Mn= ',F8.4, '[GeV]')
  196 FORMAT ('#Nucleon radius: rn= ',F8.4, '[fm]')
  198 FORMAT ('#String tension parameter: A= ',F8.4)
  199 FORMAT ('#Normal nuclear density: rho00= ',F8.4, '[GeV/fm3]')
  201 FORMAT ('#Number of hydro cells per nuclear radius: tnc= ',I3)
  202 FORMAT ('#########SIMULATION CONTROL PARAMETERS #########')
  203 FORMAT ('#Number of simulations: nteventos= ',I7)
  204 FORMAT ('#Vorticity implementation: vortin= ',I2)
  206 FORMAT ('#Print fort.100 file for hydro: fort100= ',I2)
  207 FORMAT ('#Checking: check= ',I2)
  208 FORMAT ('#######################################################')
  209 FORMAT ('################## GRID INFORMATION ##################')
C
C     ######################################################################
C
      STOP
C
      ENDPROGRAM ESRM
C
C
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
C
C  
c   
c     ===============================================================
c                       SUBROUTINES AND FUNTIONS
c     ===============================================================
c
c
c
c
      SUBROUTINE nucl_rand_dist (xab1,yab1,xab2,yab2)
c     This subroutine calculates the nucleon random distribution for 
c     each nucleous: target and projectile.
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C
C      
      REAL :: aws1,aws2,dx,dy,r1,r2,rhomax1,rhomax2,rn,x1,x2,y1,y2,z1,
     1 z2
      INTEGER :: A1,A2,Nceldasx,Nceldasy,nx,ny,check,nteventos
c
      REAL, DIMENSION(A1) :: xx1,yy1,zz1
      REAL, DIMENSION(A2) :: xx2,yy2,zz2
c
      REAL :: x10,x20,y10,y20,z10,z20,rho10,rho20,rnx10,rnx20,rny10,
     1 rny20,rnz10,rnz20,rnrho10,rnrho20,rnxa1,rnxa2,rnya1,rnya2,rnza1,
     2 rnza2,xa1,xa2,ya1,ya2,za1,za2,alphaHO,aHO
     
      INTEGER :: i,j,Npuntosal
      REAL, PARAMETER :: pi=3.1415926535897932384626
      REAL, INTENT(OUT), DIMENSION(1:197,1:1000) :: xab1,yab1
      REAL, INTENT(OUT), DIMENSION(1:197,1:1000) :: xab2,yab2

c
C
C     ---------------------------------------------------------------
C
C
      COMMON /nrd/ aws1,aws2,A1,A2,dx,dy,Nceldasx,Nceldasy,nx,ny,r1,r2,
     1 rn,x1,x2,y1,y2,z1,z2
      COMMON /cond/ check,nteventos
c
C      
C     ---------------------------------------------------------------      
c      We generate randomly the center of each nucleon using a W-S 
C                             distribution
c      
c                               TARGET
c
      CALL RANDOM_SEED ()
            
      
      DO 10 i=1,MAX(A1,A2)
C      
      IF (i .GT. A1) GO TO 25
C      
   20 CONTINUE
C
      rhomax1=1.0
      rhomax2=1.0
      CALL RANDOM_NUMBER (rnx10)
      CALL RANDOM_NUMBER (rny10)
      CALL RANDOM_NUMBER (rnz10)
      CALL RANDOM_NUMBER (rnrho10)
C      
c     Generate a random number in the interval -2r1+x1<=x10<=2r1+x1
      x10=-2.*r1+x1+4.*r1*rnx10
c     Generate a random number in the interval -2r1+y1<=y10<=2r1+y1
      y10=-2.*r1+y1+4.*r1*rny10   
c     Generate a random number in the interval -2r1+z1<=z10<=2r1+z1
      z10=-2.*r1+z1+4.*r1*rnz10
c     Generate a random number between 0 and rhomax1
      rho10=rnrho10*rhomax1
C     

      
      IF (A1 .LT. 6 .AND. A1 .GE. 2) GO TO 300
      IF (A1 .LT. 20 .AND. A1 .GE. 6) GO TO 305
      IF (A1 .GE. 20) GO TO 310
      
  300 CONTINUE
!     HOMMOGENEOUS DISTRIBUTION
!      WRITE(*,*) 'HOMOGENEOUS DISTRIBUTION'
      
      IF ((x1-x10)**2+(y1-y10)**2+(z1-z10)**2 .LE. r1**2 ) THEN
C     
         GO TO 320
C         
      ELSE
C      
         GO TO 20
C        
      ENDIF
C      
  320 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx1(i)=0.
      xx1(i)=x10
c     Center of each nucleon of the target in y-direction
      yy1(i)=0.
      yy1(i)=y10
c     Center of each nucleon of the target in z-direction
      zz1(i)=0.
      zz1(i)=z10 
  
      GOTO 25
  
      
  305 CONTINUE  
!     HARMONIC OSCILATOR MODEL (HOM)
!     HOM PARAMETERS (taken from ATOMIC DATA AND NUCLEAR DATA TABLES 36, 495-536 (1987))
 !     WRITE(*,*) 'HO DISTRIBUTION'
      IF (A1 .EQ. 6) THEN
      WRITE(*,*) 'Error: Data not avalaible for A1=6'
      STOP
      ENDIF
      IF (A1 .EQ. 7) THEN
      alphaHO=0.327
      aHO=1.770
      ENDIF
      IF (A1 .EQ. 8) THEN
      WRITE(*,*) 'Error: Data not avalaible for1 A=8'
      STOP
      ENDIF
      IF (A1 .EQ. 9) THEN
      alphaHO=0.611
      aHO=1.791
      ENDIF
      IF (A1 .EQ. 10) THEN
      alphaHO=0.837
      aHO=1.710
      ENDIF
      IF (A1 .EQ. 11) THEN
      alphaHO=0.811
      aHO=1.690
      ENDIF
      IF (A1 .EQ. 12) THEN
      WRITE(*,*) 'Error: Data not avalaible for A1=12'
      STOP
      ENDIF
      IF (A1 .EQ. 13) THEN
      alphaHO=1.403
      aHO=1.635
      ENDIF
      IF (A1 .EQ. 14) THEN
      alphaHO=1.380
      aHO=1.730
      ENDIF
      IF (A1 .EQ. 15) THEN
      alphaHO=1.250
      aHO=1.810
      ENDIF	
      IF (A1 .EQ. 16) THEN
      alphaHO=1.544
      aHO=1.833
      ENDIF
      IF (A1 .EQ. 17) THEN
      alphaHO=1.498
      aHO=1.798
      ENDIF
      IF (A1 .EQ. 18) THEN
      alphaHO=1.513
      aHO=1.841
      ENDIF
      IF (A1 .EQ. 19) THEN
      WRITE(*,*) 'Error: Data not avalaible for A1=19'
      STOP
      ENDIF
  
  
      IF (rho10 .LE. (1.+alphaHO*((x1-x10)**2+(y1-y10)**2+
     1 (z1-z10)**2)/aHO**2)*EXP(-((x1-x10)**2+(y1-y10)**2+
     2 (z1-z10)**2)/aHO**2)) THEN
C     
         GO TO 315
C         
      ELSE
C      
         GO TO 20
C        
      ENDIF
C      
  315 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx1(i)=0.
      xx1(i)=x10
c     Center of each nucleon of the target in y-direction
      yy1(i)=0.
      yy1(i)=y10
c     Center of each nucleon of the target in z-direction
      zz1(i)=0.
      zz1(i)=z10 
  
      GOTO 25
     
  310 CONTINUE
!     WOODS-SAXON DISTRIBUTION
!      WRITE(*,*) 'WS DISTRIBUTION'
  
      IF (rho10 .LE. 1./(1.+EXP((SQRT((x1-x10)**2+(y1-y10)**2+(z1-z10)
     1 **2)-r1)/aws1))) THEN
C     
         GO TO 30
C         
      ELSE
C      
         GO TO 20
C        
      ENDIF
C      
   30 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx1(i)=0.
      xx1(i)=x10
c     Center of each nucleon of the target in y-direction
      yy1(i)=0.
      yy1(i)=y10
c     Center of each nucleon of the target in z-direction
      zz1(i)=0.
      zz1(i)=z10 
C      
   25 CONTINUE
C   
c      
c                               PROJECTILE
c
C   
      IF (i .GT. A2) GO TO 10
C      
   50 CONTINUE
C      
      CALL RANDOM_NUMBER (rnx20)
      CALL RANDOM_NUMBER (rny20)
      CALL RANDOM_NUMBER (rnz20)
      CALL RANDOM_NUMBER (rnrho20)
C      
c     Generate a random number in the interval -2r1+x1<=x10<=2r1+x1s
      x20=-2.*r2+x2+4.*r2*rnx20
c     Generate a random number in the interval -2r1+y1<=y10<=2r1+y1
      y20=-2.*r2+y2+4.*r2*rny20   
c     Generate a random number in the interval -2r1+z1<=z10<=2r1+z1
      z20=-2.*r2+z2+4.*r2*rnz20
c     Generate a random number between 0 and rhomax1
      rho20=rnrho20*rhomax2
      
      
      IF (A2 .LT. 6 .AND. A2 .GE. 2) GO TO 330
      IF (A2 .LT. 20 .AND. A2 .GE. 6) GO TO 335
      IF (A2 .GT. 20) GO TO 340
      
  330 CONTINUE
!     HOMMOGENEOUS DISTRIBUTION
  
      IF ((x2-x20)**2+(y2-y20)**2+(z2-z20)**2 .LE. r2**2) THEN
C     
         GO TO 345
C         
      ELSE
C      
         GO TO 50
C        
      ENDIF
C      
  345 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx2(i)=0.
      xx2(i)=x20
c     Center of each nucleon of the target in y-direction
      yy2(i)=0.
      yy2(i)=y20
c     Center of each nucleon of the target in z-direction
      zz2(i)=0.
      zz2(i)=z20 
  
      GOTO 10
  
      
  335 CONTINUE  
!     HARMONIC OSCILATOR MODEL (HOM)
!     HOM PARAMETERS (taken from ATOMIC DATA AND NUCLEAR DATA TABLES 36, 495-536 (1987))
      IF (A2 .EQ. 6) THEN
      WRITE(*,*) 'Error: Data not avalaible for A2=6'
      STOP
      ENDIF
      IF (A2 .EQ. 7) THEN
      alphaHO=0.327
      aHO=1.770
      ENDIF
      IF (A2 .EQ. 8) THEN
      WRITE(*,*) 'Error: Data not avalaible for A2=8'
      STOP
      ENDIF
      IF (A2 .EQ. 9) THEN
      alphaHO=0.611
      aHO=1.791
      ENDIF
      IF (A2 .EQ. 10) THEN
      alphaHO=0.837
      aHO=1.710
      ENDIF
      IF (A2 .EQ. 11) THEN 
      alphaHO=0.811
      aHO=1.690
      ENDIF
      IF (A2 .EQ. 12) THEN
      WRITE(*,*) 'Error: Data not avalaible for A2=12'
      STOP
      ENDIF
      IF (A2 .EQ. 13) THEN
      alphaHO=1.403
      aHO=1.635
      ENDIF
      IF (A2 .EQ. 14) THEN
      alphaHO=1.380
      aHO=1.730
      ENDIF
      IF (A2 .EQ. 15) THEN
      alphaHO=1.250
      aHO=1.810
      ENDIF
      IF (A2 .EQ. 16) THEN
      alphaHO=1.544
      aHO=1.833
      ENDIF
      IF (A2 .EQ. 17) THEN
      alphaHO=1.498
      aHO=1.798
      ENDIF
      IF (A2 .EQ. 18) THEN
      alphaHO=1.513
      aHO=1.841
      ENDIF
      IF (A2 .EQ. 19) THEN
      WRITE(*,*) 'Error: Data not avalaible for A2=19'
      STOP
      ENDIF
  
  
      IF (rho20 .LE. (1.+alphaHO*((x2-x20)**2+(y2-y20)**2+(z2-z20)**2)/
     1 aHO**2)*EXP(-((x2-x20)**2+(y2-y20)**2+(z2-z20)**2)/aHO**2)) 
     2 THEN
C     
         GO TO 350
C         
      ELSE
C      
         GO TO 50
C        
      ENDIF
C      
  350 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx2(i)=0.
      xx2(i)=x20
c     Center of each nucleon of the target in y-direction
      yy2(i)=0.
      yy2(i)=y20
c     Center of each nucleon of the target in z-direction
      zz2(i)=0.
      zz2(i)=z20 
  
      GOTO 10
     
  340 CONTINUE
!     WOODS-SAXON DISTRIBUTION
  
      IF (rho20 .LE. 1./(1.+EXP((SQRT((x2-x20)**2+(y2-y20)**2+(z2-z20)
     1 **2)-r2)/aws2))) THEN
C     
         GO TO 355
C         
      ELSE
C      
         GO TO 50
C        
      ENDIF
C      
  355 CONTINUE 
C   
c     Center of each nucleon of the target in x-direction
      xx2(i)=0.
      xx2(i)=x20
c     Center of each nucleon of the target in y-direction
      yy2(i)=0.
      yy2(i)=y20
c     Center of each nucleon of the target in z-direction
      zz2(i)=0.
      zz2(i)=z20 
C       
C
C      
   10 CONTINUE
C      
c  
c     ---------------------------------------------------------------
C       WE CONSIDER EACH NUCLEON FORMED BY ONE THOUSAND FICTITIOUS 
c                    PARTICLES RANDOMLY DISTRIBUTED
c
c     Number of fictitious particles: Npuntosal
C
      Npuntosal=1000
c
c
      CALL RANDOM_SEED ()
c
      DO 170 i=1,MAX(A1,A2)
         DO 180 j=1,Npuntosal
C         
      IF (i .GT. A1) GO TO 191
c
c                               TARGET
c
C         
  190 CONTINUE
C  
             CALL RANDOM_NUMBER (rnxa1)
             CALL RANDOM_NUMBER (rnya1)
             CALL RANDOM_NUMBER (rnza1)
C             
c     We generate random numbers within the interval xx-rn<=x<=xx+rn
             xa1=(xx1(i)-rn)+2.*rn*rnxa1
c     We generate random numbers within the interval yy-rn<=y<=yy+rn
             ya1=(yy1(i)-rn)+2.*rn*rnya1
c     We generate random numbers within the interval zz-rn<=z<=zz+rn
             za1=(zz1(i)-rn)+2.*rn*rnza1
C             
      IF ((xa1-xx1(i))**2+(ya1-yy1(i))**2+(za1-zz1(i))**2 .LE. 
     1 rn**2) THEN
C     
       GO TO 200
C       
      ELSE
C      
       GO TO 190
C       
      ENDIF
C      
  200    CONTINUE
C  
         xab1(i,j)=0.
         xab1(i,j)=xa1
         yab1(i,j)=0.
         yab1(i,j)=ya1
C         
  191 CONTINUE
C
      IF (i .GT. A2) GO TO 180
c      
c                             PROJECTILE
c
C         
  260 CONTINUE
C  
             CALL RANDOM_NUMBER (rnxa2)
             CALL RANDOM_NUMBER (rnya2)
             CALL RANDOM_NUMBER (rnza2)
C             
c     We generate random numbers within the interval xx-rn<=x<=xx+rn
             xa2=(xx2(i)-rn)+2.*rn*rnxa2
c     We generate random numbers within the interval yy-rn<=y<=yy+rn
             ya2=(yy2(i)-rn)+2.*rn*rnya2
c     We generate random numbers within the interval zz-rn<=z<=zz+rn
             za2=(zz2(i)-rn)+2.*rn*rnza2
C             
C             
      IF ((xa2-xx2(i))**2+(ya2-yy2(i))**2+(za2-zz2(i))**2 .le. rn**2)
     1 THEN
C     
       GO TO 250
C       
      ELSE
C      
       GO TO 260
C       
      ENDIF
C      
  250    CONTINUE
C  
         xab2(i,j)=0.
         xab2(i,j)=xa2
         yab2(i,j)=0.
         yab2(i,j)=ya2
C         
  180    CONTINUE
  170 CONTINUE
C  
C  
      RETURN
C      
      END SUBROUTINE nucl_rand_dist
c
c
c
C     ---------------------------------------------------------------
C                            XPLUSC1 FUNCTION
C     ---------------------------------------------------------------
C
C
      REAL FUNCTION xplusC1 (x,x0)
C
C      
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C
C     
      REAL, DIMENSION(1:122,1:82) :: sigma_s
      REAL :: alpha,rho00,rap00
      INTEGER :: iiig,lg,k
      REAL Co2
      REAL, DIMENSION(1:122,1:82):: aa2,b1,d1,rap01,tauminus,e01	
      REAL, INTENT(IN) :: x
      REAL, INTENT(IN) :: x0
C      
C
C     ---------------------------------------------------------------
C
C      
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
      COMMON /inxplusC1/ aa2,b1,d1,rap01,tauminus,e01
C
C
C     ---------------------------------------------------------------
C
C     
      xplusC1=0.5*2.*x0-(d1(iiig,lg)/b1(iiig,lg))*x+
     1 (alpha*e01(iiig,lg)*(1.+Co2)*EXP(2.*rap01(iiig,lg))
     2 /2./sigma_s(iiig,lg)/rho00/EXP(rap00))*(d1(iiig,lg)/b1(iiig,lg)+
     3 EXP(-2.*rap01(iiig,lg)))*((1.-x/tauminus(iiig,lg))**
     4 (-2.*sigma_s(iiig,lg)*rho00*EXP(rap00)/alpha/aa2(iiig,lg))-1.)
C      
C
C     ---------------------------------------------------------------
C
C      
      RETURN
C      
      END FUNCTION xplusC1
C
C
c
c     ---------------------------------------------------------------
c                           XMINUSC2 FUNCTION
C     ---------------------------------------------------------------
C
C
      REAL FUNCTION xminusC2 (xx,xx0)
C
C      
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C
C      
      REAL, DIMENSION(1:122,1:82) :: sigma_s
      REAL :: alpha,rho00,rap00
      INTEGER :: iiig,lg,k
      REAL :: Co2
      REAL, DIMENSION(1:122,1:82) :: aa1,b2,d2,rap02,tauplus,e02	
      REAL, INTENT(IN) :: xx
      REAL, INTENT(IN) :: xx0
C
C
C     ---------------------------------------------------------------
C
C      
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
      COMMON /inxminusC2/ aa1,b2,d2,rap02,tauplus,e02
C
C
C     ---------------------------------------------------------------  
C
C
       xminusC2=-0.5*2.*xx0-(d2(iiig,lg)/b2(iiig,lg))*xx+
     1 (alpha*e02(iiig,lg)*(1.+Co2)*EXP(2.*rap02(iiig,lg))
     2 /2./sigma_s(iiig,lg)/rho00/EXP(rap00))*(d2(iiig,lg)/b2(iiig,lg)+
     3 EXP(-2.*rap02(iiig,lg)))*((1.-xx/tauplus(iiig,lg))**
     4 (-2.*sigma_s(iiig,lg)*rho00*EXP(rap00)/alpha/aa1(iiig,lg))-1.)
C
C
C     ---------------------------------------------------------------
C
C       
      RETURN
C      
      END FUNCTION xminusC2
c
C
C
c     ---------------------------------------------------------------
C                     FIND SOLUTION SUBROUTINE
c     ---------------------------------------------------------------
C
C
      SUBROUTINE find_solution (xx1,xx2,kk,xx0,tt0,MAXTT,find_sol)
C
C      
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C
C      
      REAL, DIMENSION(1:122,1:82) :: sigma_s
      REAL :: alpha,rho00,rap00
      INTEGER :: iiig,lg,k
      REAL :: Co2 
      REAL, INTENT(INOUT) :: xx1,xx2
      INTEGER, INTENT(IN) :: kk
      REAL, INTENT(IN) :: xx0
      REAL, INTENT(IN) :: tt0
      REAL, INTENT(IN) :: MAXTT
      REAL, INTENT(OUT) :: find_sol
      REAL f1,f2,f,xplusC1,xminusC2
      REAL xx
c
c
C     ---------------------------------------------------------------
C      
C
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
C      
C      
C     ---------------------------------------------------------------
C
C     
      IF (kk .EQ. 1) THEN
C      
         f1=0.5*(xplusC1(xx1,xx0)+xx1+tt0)-MAXTT
C         
      ELSEIF (kk .EQ. 2) THEN
C      
         f1=0.5*(xminusC2(xx1,xx0)+xx1+tt0)-MAXTT
C      
      ENDIF
C      
C      
C      
      IF (kk .EQ. 1) THEN
C      
         f2=0.5*(xplusC1(xx2,xx0)+xx2+tt0)-MAXTT
C      
      ELSEIF (kk .EQ. 2) THEN
C      
         f2=0.5*(xminusC2(xx2,xx0)+xx2+tt0)-MAXTT
C      
      ENDIF
C      
C      
C      
      COND1: IF (f1*f2 .GT. 0.) THEN
C      
         WRITE(*,100)
C         
         find_sol=xx1
C      
      ELSE
C      
   10 COND2: IF (xx2-xx1 .GT. 0.0001) THEN
C
                 xx=0.5*(xx1+xx2)
C                 
                IF (kk .EQ. 1) THEN
C                
                  f1=0.5*(xplusC1(xx1,xx0)+xx1+tt0)-MAXTT
C                  
                ELSEIF (kk .EQ. 2) THEN
C                
                  f1=0.5*(xminusC2(xx1,xx0)+xx1+tt0)-MAXTT
C                  
                ENDIF
C      
                IF (kk .EQ. 1) THEN
C                
                  f=0.5*(xplusC1(xx,xx0)+xx+tt0)-MAXTT
C                  
                ELSEIF (kk .EQ. 2) THEN
C                
                  f=0.5*(xminusC2(xx,xx0)+xx+tt0)-MAXTT
C                  
                ENDIF
C      
                IF (f1*f .GT. 0.) THEN
C                
                   xx1=xx
C                   
                ELSE
C                
                   xx2=xx
C                   
                ENDIF
C                
                GO TO 10
C                
           ENDIF COND2
C      
              find_sol=xx
C              
      ENDIF COND1
C
C
C     ---------------------------------------------------------------     
C
C   
  100 FORMAT('no solution')    
C
C 
C     ---------------------------------------------------------------
C
C  
      RETURN
C      
      END SUBROUTINE find_solution
c
C
C
c     ---------------------------------------------------------------
C                       TRAJECTORY111 FUNCTION
c     ---------------------------------------------------------------
C
C
      REAL FUNCTION trajectory111(T,k1)
C
C
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------    
C
C  
      REAL, DIMENSION(1:122,1:82) :: sigma_s
      REAL :: alpha,rho00,rap00
      INTEGER :: iiig,lg,k
      REAL :: Co2 
      REAL, DIMENSION(1:122,1:82) :: aa2,b1,d1,rap01,tauminus,e01
      REAL, DIMENSION(1:122,1:82) :: thy
      REAL, DIMENSION(1:123,1:83) :: MAXT
      REAL, DIMENSION(1:122,1:82) :: xturn1
      REAL, DIMENSION(1:4001,1:122,1:82) :: tC1,zC1
      REAL, DIMENSION(1:4001) :: z
      REAL, DIMENSION(1:122,1:82) :: t0
      REAL, INTENT(IN) :: T
      INTEGER, INTENT(IN) :: k1
      REAL :: xp1,zzz,find_sol,xplusC1
      REAL :: xm1,xxx
      INTEGER :: xxx1,kk
c
C
C     ---------------------------------------------------------------
C
C
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
      COMMON /inxplusC1/ aa2,b1,d1,rap01,tauminus,e01
      COMMON /intraj111/ thy,MAXT,xturn1,tC1,zC1,z,t0
c
c
C     ---------------------------------------------------------------
C
C
      xxx1=1
      xxx=0.
      kk=1
C      
C      
      COND1: IF (T .LT. MAXT(iiig,lg)) THEN
C      
      COND2:   IF (T .LT. tC1(k1,iiig,lg)) THEN
C         
             CALL find_solution(xxx,xturn1(iiig,lg),kk,
     1 z(k1),t0(iiig,lg),T,find_sol)
C
            xm1=find_sol
            xp1=xplusC1(xm1,z(k1))
            trajectory111=0.5*(xp1-t0(iiig,lg)-xm1)
C
         ELSE
C         
c              Motion in maximum extended frame 
C
            trajectory111=zC1(k1,iiig,lg)+thy(iiig,lg)*(T-
     1 tC1(k1,iiig,lg))
C     
         ENDIF COND2
C
      ELSE
C      
      zzz=zC1(k1,iiig,lg)+thy(iiig,lg)*(MAXT(iiig,lg)-
     1 tC1(k1,iiig,lg))
C      
      trajectory111=zzz-(T-MAXT(iiig,lg))
C      
      ENDIF COND1
C      
C
C     ---------------------------------------------------------------
C
C      
      RETURN
C      
      END FUNCTION trajectory111
C
c
C
c     ---------------------------------------------------------------
c                        TRAJECTORY222 FUNCTION
C     ---------------------------------------------------------------
C
C
      REAL FUNCTION trajectory222(T,k2)
C
C      
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C
C      
      REAL, DIMENSION(1:122,1:82) :: sigma_s
      REAL :: alpha,rho00,rap00
      INTEGER :: iiig,lg,k
      REAL :: Co2 
      REAL, DIMENSION(1:122,1:82) :: aa1,b2,d2,rap02,tauplus,e02
      REAL, DIMENSION(1:122,1:82) :: thy
      REAL, DIMENSION(1:123,1:83) :: MAXT
      REAL, DIMENSION(1:122,1:82) :: xturn1
      REAL, DIMENSION(1:4001,1:122,1:82) :: tC1,zC1
      REAL, DIMENSION(1:4001) :: z
      REAL, DIMENSION(1:122,1:82) :: t0
      REAL, DIMENSION(1:122,1:82) :: xturn2
      REAL, DIMENSION(1:4001,1:122,1:82) :: tC2,zC2
      REAL, INTENT(IN) :: T
      INTEGER, INTENT(IN) :: k2
      REAL :: xm2,zzz,find_sol,xminusC2
      REAL ::xp2,zeros,xxx
      INTEGER :: kk
c
C
C     ---------------------------------------------------------------
C
C
      COMMON /inxplusminus/ sigma_s,alpha,rho00,rap00,Co2,iiig,lg,k
      COMMON /inxminusC2/ aa1,b2,d2,rap02,tauplus,e02
      COMMON /intraj111/ thy,MAXT,xturn1,tC1,zC1,z,t0
      COMMON /intraj222/ xturn2,tC2,zC2
c
c
C     ---------------------------------------------------------------
C
C
      xxx=0.
      xp2=0.
      zeros=0.
      kk=2
C      
      COND1: IF (T .LT. MAXT(iiig,lg)) THEN
C      
      COND2:   IF (T .lt. tC2(k2,iiig,lg)) THEN
C         
             CALL find_solution(xxx,xturn2(iiig,lg),kk,
     1 z(k2),t0(iiig,lg),T,find_sol)
C     
            xp2=find_sol
            xm2=xminusC2(xp2,z(k2))
            trajectory222=0.5*(xp2+t0(iiig,lg)-xm2)
C            
         ELSE
C         
c              Motion in maximum extended frame 
C
            trajectory222=zC2(k2,iiig,lg)+thy(iiig,lg)*(T-
     1 tC2(k2,iiig,lg))
C     
         ENDIF COND2
C
      ELSE
C      
      zzz=zC2(k2,iiig,lg)+thy(iiig,lg)*(MAXT(iiig,lg)-
     1 tC2(k2,iiig,lg))
C      
      trajectory222=zzz+(T-MAXT(iiig,lg))
C      
      ENDIF COND1
C
C
C     ---------------------------------------------------------------
C      
C      
      RETURN
C      
      END FUNCTION trajectory222
