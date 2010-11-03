      PROGRAM MAIN
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: GFDL_OCEAN.WESTATL
C   PRGMMR: frolov           ORG: NP22         DATE: 1999-06-02
C
C ABSTRACT: This program runs Princeton Ocean Model configured for
C   west atlantic domain. The entire code is split into two main 
C   subroutines: OCEANINIT and OSEANSTEP. OCEANINIT reads all the input
C   data and initializes model variables. OCEANSTEP run the model for
C   one time step.
C
C PROGRAM HISTORY LOG:
C   --------  BLUMBERG/MELLOR - Original implementation
C   99-06-02  frolov - restructured, reconfigured and modified
C
C INPUT FILES:
C   UNIT    8    file gdem, contains ascii ocean climatalogy data
C   UNIT   10    PARAMETERS.inp file, contains configuration parameters
C   UNIT   12    gfdl_ocean_topography (ascii ocean floor topography data)
C   UNIT   14    binary ocean restart filefile,
C                name specified on the last of PARAMETERS.inp
C   UNIT   15    storm history file, file name specified in PARAMETERS.inp
C   UNIT   21    sst.gfs.dat file, contains ascii SST data
C   UNIT   22    mask.gfs.dat file, contains ascii GFS landsea mask data
C   UNIT   23    gfdl_gfs_latlon_t170 (lons & lats for GFS 512x256 T170)
C
C OUTPUT FILES:
C   UNIT   31    binary ocean surface elevation data EL.YYMMDDHH
C   UNIT   32    binary vertically averaged velocity UVA.YYMMDDHH
C   UNIT   35    binary 3-D ocean temperature T.YYMMDDHH
C   UNIT   37    binary 3-D ocean U-velocity U.YYMMDDHH
C   UNIT   38    binary 3-D ocean V-velocity V.YYMMDDHH
C   UNIT   39    binary wind stress data TXY.YYMMDDHH
C   UNIT   39    binary wind speed data TXY.YYMMDDHH
C   UNIT   90    binary ocean restart file RST.YYMMDDHH
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$
      include 'comblk.h'
      include 'comblk1.h'

c      call w3tagb('GFDL_OCEAN.WESTATL',2001,0082,0003,'NP22')

      CALL OCEANINIT
      DO IINT=1,IEND
        CALL OCEANSTEP(0)
      END DO
c---- for debug
      i0=46
      j0=52
      print *,' check before RST: (t(i0,j0,k),k=1,KB) KB=',KB
      write(6,101) (t(i0,j0,k),k=1,KB)
 101  format(10f7.2)
C--------------------------------------------------------------------
C          SAVE RESTART FILE in FRSTO 
C--------------------------------------------------------------------
      OPEN(90,FILE=FRSTO,STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(90) TIME,
     1 WUBOT,WVBOT,AAM2D,UA,UAB,VA,VAB,EL,ELB,ET,ETB,EGB,
c     2 UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,RHO,ADVUU,ADVVV,ADVUA,ADVVA,
     2 UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,ADVUU,ADVVV,ADVUA,ADVVA,
     3 KM,KH,KQ,L,Q2,Q2B,AAM,Q2L,Q2LB
      WRITE(6,'('' Restart File saved at hour:'',F8.1)') TIME*24.0
c---  for debug
      rewind(90)
      READ(90) TIME,
     1 WUBOT,WVBOT,AAM2D,UA,UAB,VA,VAB,EL,ELB,ET,ETB,EGB,
c     2 UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,RHO,ADVUU,ADVVV,ADVUA,ADVVA,
     2 UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,ADVUU,ADVVV,ADVUA,ADVVA,
     3 KM,KH,KQ,L,Q2,Q2B,AAM,Q2L,Q2LB
      WRITE(6,'('' Restart File saved at hour:'',F8.1)') TIME*24.0
c
      CLOSE(90)
c
      print *,' check after RST: (t(i0,j0,k),k=1,KB) KB=',KB
      write(6,101) (t(i0,j0,k),k=1,KB)

c      call w3tage('GFDL_OCEAN.WESTATL')

      STOP
      END
