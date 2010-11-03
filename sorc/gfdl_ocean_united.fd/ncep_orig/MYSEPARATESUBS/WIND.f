      SUBROUTINE WIND(filename,startdate)
      INCLUDE 'comblk.h'
      INTEGER PLN
      PARAMETER(PLN=100)
c
      integer hour, lat, long, mx, rmw
      integer*4 startdate,date
      integer day,month,year
      integer garb(5),Rd1(4),Rd2(4)
      character*19 name
      character*15 filename
      character*1 letter
c
      DIMENSION X(PLN),Y(PLN),TM(PLN),PRES(PLN),PRES0(PLN),
     *      RMAXa(PLN),WSPMAX(PLN),
     *      R18v(PLN,4),R26v(PLN,4),Rref18v(5),Rref26v(5),alphv(5)
      DIMENSION RAD(14),WS(14),RADM(14),WSM(14),ANGL(14)
      REAL CMP,T1,T2,F0,F1,L0,L1,REARTH,R,A7,B,E,DELP,x0,y0
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      REAL DELTAX,DELTAX1,DELTAY,DELTAY1,DXDY,julday
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX

      DATA RM,R0/60.E3,480.E3/,RHO_0/1024.E0/
      DATA RAD/0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75
     * ,8.1,10.8,13.5/
      DATA WS/0.,.1,.5,.8,.95,1.,.97,.72,.54,.44,.4,.36
     * ,.27,.23/
      DATA ANGL/0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,
     * 21.,21.,21.,21./
c     print*,'In subroutine WIND ...'
c
      WIND_SCALE=0.8
      ROA=1.28
      RMAX=50.E3
      PI=3.1415927
      E=exp(1.)
c
c----------------------- Reading message file --------------------
  17  format(A19,I6,1x,I4,1x,I3,1x,I4,1x,I3,1x,I3,3I5,
     * 1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,
     * 1x,I4,1x,I4,1x,I4)
c     print*,'reading file ',filename
      open(15,file=filename,status='old')
      end=0.
      I=0
      do while(end.eq.0)
        read(15,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
c----------- falk 09-12-05 change output
       if(MOD(IINT,24).EQ.0) then
        print*,'reading file ',filename
        write(6,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
       end if
        if(date.eq.0) goto 20
        I=I+1
c
        date=date*100+hour/100
        call date2day(year,julday,date)
c
        TM(i)=julday
        X(i)=-long/10.
        Y(i)=lat/10.
        PRES(i)=float(garb(3))
        PRES0(i)=float(garb(4))
        WSPMAX(i)=float(mx)
        RMAXa(i)=float(rmw)
        do n=1,4
          n1=n+(1-mod(n,2))*sign(2,3-n)
          R18v(i,n)=Rd1(n1)
          R26v(i,n)=Rd2(n1)
          if(wspmax(i).le.26.or.R26v(i,n).le.RMAXa(i)) R26v(i,n)=-999
          if(wspmax(i).le.18.or.R18v(i,n).le.RMAXa(i)) R18v(i,n)=-999
          if(R26v(i,n).gt.R18v(i,n)) R26v(i,n)=-999
        end do
      end do
  20  end=1.
      ipmax=I
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) then
      print*,'Number of hurricane path datapoints read: ',ipmax
c---------------- falk 01-05-04 change printing
c     print*,'tm=',(tm(i),i=1,ipmax)
      write(6,101) (tm(i),i=1,ipmax)
 101  format(10f7.2)
      print*,'year=',year
      print*,'startdate=',startdate
      end if
      close(15)
c
c--------------------- Calculating starting day -----------------
c
      call date2day(year,julday,startdate)
      do i=1,ipmax
        TM(i)=TM(i)-julday
      end do
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTERPOLATION OF HURRICANE PATH TO DETERMINE THE CURRENT POSITION

      CMP=TM(ipmax)
      if(cmp.lt.time.or.time.lt.tm(1)) then
c---------- falk 01-05-04 change printing
        print*,'NO HURRICANE PATH DATA FOR THIS TIME'
        print*,'   time=',time
        print*,'   tm(1)=',tm(1)
        print*,'   tm(ipmax)=',tm(ipmax)
        print*,'   day=',julday
        RETURN
      end if

      call verinterp(ipmax,1,TM,X,TIME,F0)
      call verinterp(ipmax,1,TM,Y,TIME,L0)
      call verinterp(ipmax,1,TM,PRES,TIME,PRES1)
      call verinterp(ipmax,1,TM,PRES0,TIME,PRES2)
      DELP=(PRES2-PRES1)*100.
      call verinterp(ipmax,1,TM,WSPMAX,TIME,WSMAX)
      WSMAX=WSMAX*WIND_SCALE
      WS18=18*WIND_SCALE
      WS26=26*WIND_SCALE
      call verinterp(ipmax,1,TM,RMAXa,TIME,RMAX)
      RMAX=RMAX*1.e3
      do n=1,4
        call interp1d(-999.0,ipmax,1,TM,R18v(1,n),TIME,Rref18v(n))
        if(Rref18v(n).ne.-999) Rref18v(n) = Rref18v(n)*1.e3
        call interp1d(-999.0,ipmax,1,TM,R26v(1,n),TIME,Rref26v(n))
        if(Rref26v(n).ne.-999) Rref26v(n) = Rref26v(n)*1.e3
        alphv(n) = (n-1)*pi/2
      end do
      do n=2,6
        n1=mod(n-1,4)+1
        nm1=mod(n-2,4)+1
        np1=mod(n,4)+1
        if(Rref18v(n1).eq.-999) then
          if(Rref18v(nm1).ne.-999) then
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=0.5*(Rref18v(nm1)+Rref18v(np1))
            else
              Rref18v(n1)=Rref18v(nm1)
            end if
          else
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=Rref18v(np1)
            else
              Rref18v(n1)=-999
            end if
          end if
        end if
        if(Rref26v(n1).eq.-999) then
          if(Rref26v(nm1).ne.-999) then
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=0.5*(Rref26v(nm1)+Rref26v(np1))
            else
              Rref26v(n1)=Rref26v(nm1)
            end if
          else
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=Rref26v(np1)
            else
              Rref26v(n1)=-999
            end if
          end if
        end if
      end do
c
      Rref18v(5) = Rref18v(1)
      Rref26v(5) = Rref26v(1)
      alphv(5) = alphv(4)+pi/2
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) then
      print*,'Time=',Time
      print*,'Current hurricane position (x,y): ',f0,l0
      print*,'WSMAX=',WSMAX,'; DELP=',DELP,'; RMAX=',RMAX
      print*,'Rref18v=',Rref18v
      print*,'Rref26v=',Rref26v
      end if
c
      x0=f0
      y0=l0
      F0=F0*2.*PI/360
      L0=L0*2.*PI/360
C--- F0,L0 ARE THE CARRENT (LONGITUDE,LATTITUDE) COORDINATES OF THE HURRICANE
C
C--- CALCULATING UTX AND UTY (HURRICANE SPEED)
c
      cmp=tm(ipmax)
      do i=1,ipmax
        if(abs(tm(i)-time).le.cmp) then
          cmp=abs(tm(i)-time)
          ii=i
        end if
      end do
c----------- falk 01-05-04 not to go out bnd
      if((tm(ii)-time).le.0.and.ii.ne.ipmax) then
        t1=tm(ii)
        t2=tm(ii+1)
        x1=x(ii)
        x2=x(ii+1)
        y1=y(ii)
        y2=y(ii+1)
      else
        t2=tm(ii)
        t1=tm(ii-1)
        x2=x(ii)
        x1=x(ii-1)
        y2=y(ii)
        y1=y(ii-1)
      end if
      deltax1=rearth*cos(l0)*(x2-x1)*2.*pi/360
      deltay1=rearth*(y2-y1)*2.*pi/360
      utx=deltax1/((t2-t1)*24.*3600.)
      uty=deltay1/((t2-t1)*24.*3600.)
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'utx,uty: ',utx,uty
C
C--- CALCULATING PARAMETERS FOR WIND PROFILE FORMULA
c
      B=WSMAX**2*E*ROA/DELP
      A7=RMAX**B
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'B= ',B
c
      do i=1,14
        RADM(I)=RMAX*RAD(I)
      end do
C
      DO 350 J=1,JM
      DO 351 I=1,IM

C  CALCULATING U-WIND STRESS FOR I,J POINT OF U-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/((IM-0.5)-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY) +
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

c      if(R.GT.RADM(14).OR.R.EQ.0) then
c        k=14
c      else
c        k=1
c        do while(R.GE.RADM(k+1))
c          k=k+1
c        end do
c      end if

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
c      WF=WND
c      WR=-WND*TAN(RANGL)
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.LE.35.) then
      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0
      else
      WUSURF(I,J)=WUSURF(I,J)-(3.3368+
     1      (WM-34.0449)**0.3)*WX/(WM*RHO_0)
      end if


C  CALCULATING V-WIND STRESS FOR I,J POINT OF V-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/((JM-0.5)-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
c------------ falk 01-05-04 check DELTAX
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.le.35.) then
      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0
      else
      WVSURF(I,J)=WVSURF(I,J)-(3.3368+
     1      (WM-34.0449)**0.3)*WY/(WM*RHO_0)
      end if

C  CALCULATING WIND STRESS FOR I,J POINT OF DEPTH GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360.
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360.
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
c      WINDX(i,j)=WX
c      WINDY(i,j)=WY
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.GT.35.) THEN
      TMAX=3.3368+(WM-34.0449)**0.3
      TAUX(I,J)=TAUX(I,J)+TMAX*WX/WM
      TAUY(I,J)=TAUY(I,J)+TMAX*WY/WM
      ELSE
      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX

C      WTSURF(I,J)=30.*WM/(4000*1024.)

      END IF

 351  CONTINUE
 350  CONTINUE

 230  FORMAT(151e16.7)
 231  FORMAT(10F8.3)
c
c------------------------- falk 01-05-04 add call avrtau
      call avrtau(x0,y0,0)
c
c----------- falk 09-12-05 change output
      if(MOD(IINT,24).EQ.0) 
     *print*,'Exiting WIND ...'
      RETURN
      END
