        
        subroutine newwind(time,uu,vv,wusurf,wvsurf,taux,tauy)
c
c Interpolate HRD 10m wind and create wind stress for ocean model
c Biju Thomas  on 12/22/2005
c

        parameter(im=254,jm=225,nhrd=33)
        parameter(ck=0.4,roa=1.28,rho_o=1024.e0)
        integer year1,year2,yeari
        integer mon1,mon2,moni,day1,day2,dayi
        integer smon,sday,syer,shor,smin,ssec
        integer dmon,dday,dyer,dhor,dmin,dsec
        integer hour1,hour2,houri,min1,min2,mini,sec1,sec2,seci
        integer month(nhrd),day(nhrd),year(nhrd)
        integer hour(nhrd),min(nhrd),sec(nhrd)
        double precision sjul,djul,jul(nhrd),jul1,jul2,juli
        character*5 hur,dhur
        character*1 dt1,dt2,dt3
        real uu(im,jm),vv(im,jm),wind(im,jm)
        real cd(im,jm),wusurf(im,jm),wvsurf(im,jm)
        real taux(im,jm),tauy(im,jm)
 

        open(11,file='parameter1.inp')
        read(11,*)hur
        read(11,*)syer
        read(11,*)smon
        read(11,*)sday
        read(11,*)shor
        
        
        smin=00
        ssec=00
        call jhour(smon,sday,syer,shor,smin,ssec,sjul,ier)
        juli=sjul+time*24.

 
 16     format(A5,I4,A1,I2,I2,A1,I2,A1,I2)
        dsec=ssec
        do k=1,nhrd
cRMY          jul(k)=sjul-4.50+(k-1)*6.
           read(11,16)dhur,dyer,dt1,dmon,dday,dt2,dhor,dt3,dmin
           call jhour(dmon,dday,dyer,dhor,dmin,dsec,djul,ier)
           jul(k)=djul
c may be this should change later ....   

cRMY          call greg(jul(k),1,2,month(k),day(k),year(k),
cRMY     &              hour(k),min(k),sec(k),ier)
        enddo

        
        jul1=jul(1)
        jul2=jul(nhrd)
        timeh=time*24.
        do k=1,nhrd
         if(jul(k).le.juli)jul1=jul(k)
         if(timeh.eq.0.25)print*,'jul(',k,') = ',jul(k)
        enddo
        do k=nhrd,1,-1
         if(jul(k).gt.juli)jul2=jul(k)
        enddo

cRMY REMOVE THIS... FOR DEBUG ONLY
c        goto 666
cccccccccccccccccccccccccccccccccc

        print*,sjul,jul1,jul2,juli
        call greg(jul1,1,2,mon1,day1,year1,hour1,min1,sec1,ier)
        call greg(jul2,1,2,mon2,day2,year2,hour2,min2,sec2,ier)
        call greg(juli,1,2,moni,dayi,yeari,houri,mini,seci,ier)

        print*
        print*,year1,mon1,day1,hour1,min1,sec1
        print*,year2,mon2,day2,hour2,min2,sec2
        print*,yeari,moni,dayi,houri,mini,seci
        print*
       call INTERP_HRD_WIND(HUR,YEAR1,MON1,DAY1,HOUR1,MIN1,
     &                        YEAR2,MON2,DAY2,HOUR2,MIN2,
     &                        YEARI,MONI,DAYI,HOURI,MINI,UU,VV)


        do 68 i=1,im
        do 68 j=1,jm
          wind(i,j)=sqrt(uu(i,j)**2+vv(i,j)**2)
          if(wind(i,j).lt.12.5)
     &     z0=(0.0281*exp(0.190669*wind(i,j)))/1000.
          if(wind(i,j).ge.12.5)
     &     z0=(0.0739793*wind(i,j)-0.61841)/1000.
           cd(i,j)=(ck/log(10./z0))**2
           wusurf(i,j)=-roa*cd(i,j)*wind(i,j)*uu(i,j)/rho_o
           wvsurf(i,j)=-roa*cd(i,j)*wind(i,j)*vv(i,j)/rho_o
           taux(i,j)=-wusurf(i,j)*rho_o
           tauy(i,j)=-wvsurf(i,j)*rho_o
68     continue

cRMY REMOVE THIS... FOR DEBUG ONLY      
c 666    continue
cccccccccccccccccccccccccccccccccc

        RETURN
        END



	subroutine INTERP_HRD_WIND(HUR,YEAR1,MON1,DAY1,HOUR1,MIN1,
     &                        YEAR2,MON2,DAY2,HOUR2,MIN2,
     &                        YEARI,MONI,DAYI,HOURI,MINI,UU,VV)
        
c----------------------------------------------------------------
c       Program INTERP_HRD_WIND read in HRD wind analysis field
c       intepolate they to a time interval of every tint hour
c       extropolate the wind profile to a large domain on every
c       5 degree bases based on the tangitial of each profile
c                Writen by: Yalin Fan, 10/05/2005
c                Graduate School of Oceanography
c                University of Rhode Island
c----------------------------------------------------------------

        parameter(m=225,n=254)
	integer i,j,k,m,n,ier,k1,itime
	integer YEAR1,YEAR2,YEARI,GN
	integer MON1,MON2,MONI,DAY1,DAY2,DAYI
	integer HOUR1,HOUR2,HOURI,MIN1,MIN2,MINI,SEC1,SEC2,SECI
	real tint,Latup,Latlow,Lonlef,Lonrt,delxy,Re,pi
	real XD(400,400),YD(400,400), X11(400,400),Y11(400,400)
        real XDU(400,400),YDU(400,400),XDV(400,400),YDV(400,400)
	real X22(400,400),Y22(400,400),X33U(400,400),Y33U(400,400)
        real X33V(400,400),Y33V(400,400),WIND33V(400,400)
	real R11(400,400),R22(400,400),R33(400,400),WIND11(400,400)
	real U11(400,400),U22(400,400),U33(400,400),WIND22(400,400)
	real V11(400,400),V22(400,400),V33(400,400),WIND33U(400,400)
	real X1(182,72),Y1(182,72),X2(182,72),Y2(182,72)
	real X3(182,72),Y3(182,72),X(200,200),Y(200,200)
	real WIND1(200,200),WIND2(200,200),WIND3(200,200)
	real Rm1(72),Rm2(72),Rm3(72),Center(4)
	real TANG_WSP(72),TANG_WSP1(72),TANG_WSP2(72)
	real WW(82),WW1(82),WW2(82),RR(81),theta(100)
	real LX1(182*72),LY1(182*72),LX2(182*72),LY2(182*72)
	real LX3(182*72),LY3(182*72),LWIND1(182*72),LWIND2(182*72)
	real LWIND3(182*72),MMR(72)
	real Cx1,Cx2,Cx3,Cy1,Cy2,Cy3
        real uu(n,m),vv(n,m)
	double precision DATE1,DATE2,DATEI	

	character*1,t1,t2,t3,t4
	character*4,c1,c2
	character*2 c3,c4,c5,c6,CN1
	character*3 CN2
	character*5 HUR
	character*20 NAME1,NAME2
	character*10 filename,COL2
	character*9 COL1

	GN=12
	delxy=1.0/GN
	Re=6371.0
	pi=3.1415927
c-------------------------------------------------------------------
c       open input file HRD_FILES, which include the latitude and 
c       longitude of the model domain (the default grid increment
c       is 1/12 degree), the time increment in hours, and the file
c       names of the HRD measurements in time increasing order. 
c-------------------------------------------------------------------

c	m=int((Latup-Latlow)*GN)+1
c	n=int((Lonrt-Lonlef)*GN)+1

        Lonrt=-50.
        Lonlef=-98.5
        Latup=47.5
        Latlow=10.       

c--------------------------- Lat, Lon of grid point in model domain
	do 12 i=1,m
	  do 11 j=1,n
	    XDV(i,j)=Lonlef+(J-1)*(Lonrt-Lonlef)/(N-1)
	    YDV(i,j)=Latlow+(I-1)*(Latup-Latlow)/((M-0.5)-1)
            XDU(i,j)=Lonlef+(J-1)*(Lonrt-Lonlef)/((N-0.5)-1)
            YDU(i,j)=Latlow+(I-1)*(Latup-Latlow)/(M-1)
11      continue
12	continue	  	 
	
c	do i=1,81
c	  RR(i)=800+40*(i-1)
c	enddo

        do 15 i=1,72
	  theta(i)=-pi/4.0+pi*(i-1.0)/36.0
 15	continue

	do i=1,200
	   do j=1,200
	      WIND1(i,j)=0.0
	      WIND2(i,j)=0.0
	      WIND3(i,j)=0.0
	      X(i,j)=0.0
	      Y(i,j)=0.0
	   enddo
	enddo
	do i=1,182
	   do j=1,72
	      X1(i,j)=0.0
	      X2(i,j)=0.0
	      X3(i,j)=0.0
	      Y1(i,j)=0.0
	      Y2(i,j)=0.0
	      Y3(i,j)=0.0
	   enddo
	enddo

16  	format(A5,I4,A1,I2,I2,A1,I2,A1,I2)
17	format(I4,I2,I2,A1,I2,I2,I2)     
118	format(A4,A4,2x,A2,A2,A2)

    	SEC1=0
601     format(i4.4)
602	format(i2.2)
603	format(i3.3)

100	continue

	write(c1,601)YEAR1
	itime=MON1*100+DAY1
	write(c2,601)itime
	write(c3,602)HOUR1
	write(c4,602)MIN1
	write(c5,602)SEC1

	filename=c2//'_'//c3//'_'//c4
	NAME1=HUR//c1//'_'//filename
	print*,NAME1
	
	SEC2=0

	write(*,118)c1,c2,c3,c4,c5

        write(c1,601)YEAR2
        itime=MON2*100+DAY2
        write(c2,601)itime
        write(c3,602)HOUR2
        write(c4,602)MIN2
	write(c5,602)SEC2

        NAME2=HUR//c1//'_'//c2//'_'//c3//'_'//c4
	print*,NAME2      

        SECI=0

	call jhour(MON1,DAY1,YEAR1,HOUR1,MIN1,SEC1,DATE1,ier)
	call jhour(MON2,DAY2,YEAR2,HOUR2,MIN2,SEC2,DATE2,ier)
        call jhour(MONI,DAYI,YEARI,HOURI,MINI,SECI,DATEI,ier)
        print*,DATE1,DATE2,DATEI
c-----------------------------------------------------------------
c       calling subroutine GET_HRD, this subroutine take HRD wind
c       field at time 1 and 2 in cartitian coordinates, turn them
c       into polar coordinates, calculate the radius of maximum 
c       wind (Rm) for every 5 degrees for each wind field, and 
c       normalize the distance from the eye by Rm. Then interpolate
c       them onto non-dimensional grid (X,Y), Center returns the 
c	center location of the hurricane
c-----------------------------------------------------------------
	call GET_HRD(NAME1,NAME2,WIND1,WIND2,X,Y,Rm1,Rm2,Center)	
c----------------------------- WIND1,WIND2 calculation OK

    	Cx1=Center(1) 
	Cy1=Center(2) 
	Cx2=Center(3) 
	Cy2=Center(4)
 
    	do 19 i=1,101
	   do 18 j=1,72
	     X1(i,j)=X(i,j)*Rm1(j)
	     Y1(i,j)=Y(i,j)*Rm1(j)
18	   continue
19	continue

c----------------------------------------------------------
c	find out the minimum distance from the boundary of 
c	the normalized domain to the center and use it to
c	define the common domain for the two wind field
c----------------------------------------------------------
	do i=1,72
	   MMR(i)=sqrt(X1(101,i)**2+Y1(101,i)**2)
	enddo
	call sort(72,MMR)
	do i=1,81
	   RR(i)=MMR(72)+40*i
	enddo
	


        do 23 j=1,72
c--------------------------------------------------------------
c	use e-folding to extrapolate the wind profile outside
c	the observation area
c---------------------------------------------------------------
	   DIST=sqrt(X1(101,j)**2+Y1(101,j)**2)
	    do 21 i=1,81
	       X1(i+101,j)=RR(i)*sin(theta(j))
	       Y1(i+101,j)=RR(i)*cos(theta(j))
21	    continue
	    do 22 k=1,81
	        WW1(k)=WIND1(101,j)*exp((DIST-RR(k))/DIST)
	       WIND1(k+101,j)=WW1(k)
22	    continue
23        continue


c-----------------------------------------------------------------
c	begin the interpolation for wind field at time between 
c	time 1 and time 2 using time increment defined in 
c	HRD_FILES
c----------------------------------------------------------------- 	
	if(DATEI .ge. DATE2) goto 999

	do 31 i=1,101
	  do 30 j=1,72
	    WIND3(i,j)=WIND1(i,j)+(DATEI-DATE1)*(WIND2(i,j)-WIND1(i,j))
     *            /(DATE2-DATE1)
	    Rm3(j)=Rm1(j)+(Rm2(j)-Rm1(j))*(DATEI-DATE1)/(DATE2-DATE1)
 30	   continue
 31	continue
	Cx3=Cx1+(Cx2-Cx1)*(DATEI-DATE1)/(DATE2-DATE1)
	Cy3=Cy1+(Cy2-Cy1)*(DATEI-DATE1)/(DATE2-DATE1)

        do 33 i=1,101
           do 32 j=1,72
              X3(i,j)=X(i,j)*Rm3(j)
              Y3(i,j)=Y(i,j)*Rm3(j)
 32	   continue
 33	continue

	do 34 j=1,72
           TANG_WSP(j)=(WIND3(91,j)-WIND3(101,j))
     *             /(sqrt(X3(101,j)**2+Y3(101,j)**2)
     *             -sqrt(X3(91,j)**2+Y3(91,j)**2))
           if(TANG_WSP(j).lt.0.0) then
              TANG_WSP(j)=-TANG_WSP(j)
           endif
           if(TANG_WSP(j).lt.0.01) then
              TANG_WSP(j)=0.01
           endif
 34	continue

         do 37 j=1,72
	    do 36 i=1,81
	       X3(i+101,j)=RR(i)*sin(theta(j))
	       Y3(i+101,j)=RR(i)*cos(theta(j))
 36	    continue
	    DIST=sqrt(X3(101,j)**2+Y3(101,j)**2)
	    do 35 k=1,81
               WW(k)=WIND3(101,j)*exp((DIST-RR(k))/DIST)
	       WIND3(k+101,j)=WW(k)
 35	    continue
 37	 continue

	do 39 i=1,m
	   do 38 j=1,n
	      X33U(i,j)=2*pi*Re*cos(YD(i,j)*pi/180)
     *                 *(XDU(i,j)-Cx3)/360.0
    	      Y33U(i,j)=2*pi*Re*(YDU(i,j)-Cy3)/360.0
	      X33V(i,j)=2*pi*Re*cos(YD(i,j)*pi/180)
     *                 *(XDV(i,j)-Cx3)/360.0
    	      Y33V(i,j)=2*pi*Re*(YDV(i,j)-Cy3)/360.0
 38	   continue
 39	continue 
 	call Rad2Rec(m,n,X3,Y3,WIND3,X33U,Y33U,WIND33U)
        call Rad2Rec(m,n,X3,Y3,WIND3,X33V,Y33V,WIND33V)

	do 41 i=1,m
	  do 40 j=1,n 
             R33(i,j)=sqrt(X33U(i,j)**2+Y33U(i,j)**2)
             U33(i,j)=-WIND33U(i,j)*Y33U(i,j)/R33(i,j)
             V33(i,j)=WIND33V(i,j)*X33V(i,j)/R33(i,j)
 40	  continue
 41	continue

c--------------------------------------------------------
c	write U, V to output file
c--------------------------------------------------------	    
 999	continue

       do 68 i=1,m
        do 68 j=1,n
          UU(j,i)=U33(i,j)
          VV(j,i)=V33(i,j)
68      continue

	return
	end


        subroutine greg(julian,ims,isym,month,day,year,hour,min,sec,ier)

c ----------------------------------------------------------------------
c       convert julian hour to gregorian date.
c
c              BY Yalin Fan 09/02/05
c
c       Based on cacm 8, (1963)  p 444 (modified).
c----------------------------------------------------------------------
        integer ims,isym
c                  -ims indicates whether integer min. and sec. values
c                    are desired ( ims=1 indicates yes, while any other
c                    value for ims gives hours and fraction of hours)
c                  -isym indicates whether a month abbreviation or an
c                    integer representation should be used ( isym=1 indicates
c                    use abbre., while any other value gives an integer.)

        integer day,year,month,hour,min,sec,ier,montna(12)
        double precision julian,hours

        data montna/'jan','feb','mar','apr','may','jun',
     1              'jul','aug','sep','oct','nov','dec'/

        if(julian.ge.0d0)go to 10
        ier=1
        write(*,5)julian
5       format(' julian hour is negative ',1pd25.17)
        return

10      ier=0
        j=idint(julian/24d0)
        hours=julian-24d0*dble(j)
        if(ims.ne.1)go to 15
        hour=dble(idint(hours))
        hours=hours-(hour)
        sec=(3600d0*hours)
        min=idint(sec/60.0d0)
        sec=sec-dble(60*min)
        go to 16
15      hour=hours
        sec=0.0d0
        min=0
16      j=j+693902
        year=(4*j-1)/146097
        j=4*j-1-146097*year
        day=j/4
        j=(4*day+3)/1461
        day=4*day+3-1461*j
        day=(day+4)/4
        month=(5*day-3)/153
        day=5*day-3-153*month
        day=(day+5)/5
        year=100*year+j
        if(month.ge.10)go to 20
        month=month+3
        go to 25
20      month=month-9
        year=year+1
25      if(isym.ne.1)return
        month=montna(month)
        return
        end


        subroutine jhour(month,day,year,hour,min,sec,julian,ier)
c-----------------------------------------------------------------
c       convert gregorian time to julian time.
c
c           By Yalin Fan 09/02/05
c-----------------------------------------------------------------

c        character*3 month,umonth,montnm(12)
        integer montda(12),month
        integer day,year,yr,cent,nk,hour,sec,min,j,munt,temp
        double precision julian,hours
c c       data montnm/'JAN','FEB','MAR','APR','MAY','JUN',
c c    1              'JUL','AUG','SEP','OCT','NOV','DEC'/
        data montda/31,29,31,30,31,30,31,31,30,31,30,31/

        ier=1
        if(year.lt.1900)go to 95

c        do 2 i=1,12
c            call toupper(month,umonth)
c            if (umonth.eq.montnm(i)) goto 3
c2       continue
c3       if (i.lt.1.or.i.gt.12) goto 105
        munt=month

c       determine whether the day of the month be plausible.
10      if(day.lt.1.or.day.gt.montda(munt))go to 110
        cent=year/100
        yr=year-100*cent
        if(munt.ne.2.or.day.ne.29)go to 20
        if(mod(yr,4).ne.0)go to 110
        if(yr.ne.0)go to 20
        if(mod(cent,4).ne.0)go to 110
c       if the hour be not integral ignore the minutes and seconds.
20      if(hour.ge.24.or.hour.lt.0)go to 115
        hours=dble(hour)
        if(min.lt.0.or.min.gt.59)go to 120
        if(sec.lt.0.0.or.sec.ge.60.0)go to 125
        hours=hours+dble(min)/60d0+dble(sec)/3600d0
c       calculate julian hour.
25      if(munt.gt.2)go to 30
        munt=munt+9
c The following line of code commented out and the next three lines
c added - KG - 2 February 2000 - a new century problem - when subtracting
c one from the year have to use the 4 digit representation instead of
c the two digit.  The problem occurs in Jan and Feb (see statement 25)
c when 9 is added to month - have to subtract one from the year.
c temp is a temporary variable.
c        yr=yr-1
        temp=year-1
        cent=temp/100
        yr=temp-100*cent
        go to 35
30      munt=munt-3
35      j=(146097*cent)/4+(1461*yr)/4+(153*munt+2)/5+day-693902
        julian=dble(j)*24d0+hours
        ier=0
        return
c       error messages
95      write(*,96)year
96      format(i12,' preceeds 1900.')
        return
105     write(*,106)month
106     format(' incorrect month abbreviation ',a5)
        return
110     write(*,111)day
111     format(' incorrect day of month specification ',i12)
        return
115     write(*,116)hour
116     format(' incorrect hour specification ',g12.3)
        return
120     write(*,121)min
121     format(' incorrect minute specification ',i12)
        return
125     write(*,126)sec
126     format(' incorrect second specification ',g12.3)
        return
        end
                                         

	
      Subroutine GET_HRD(NAME1,NAME2,WIND1itp,WIND2itp,XI,YI,Rm1,Rm2,CC)
      
      integer i,i3,j,k,k1,LON_N1,LAT_N1,LON_N2,LAT_N2,Nline
      integer Nright,Nleft,DIFN1,DIFN2,ID(20)
      real delx1,delx2,XB,XMAX,YMAX,AveRm1,AveRm2,RR1,RR2,pi,tt
      real CC(4),XMM(4),theta(73),RI
      real X1(200,200),Y1(200,200),WIND1(200,200),U1(200,200)
      real V1(200,200),X2(200,200),Y2(200,200),WIND2(200,200)
      real U2(200,200),V2(200,200),XI(200,200),YI(200,200)
      real WIND1itp(200,200),WIND2itp(200,200),y2a1(200,200)
      real y2a2(200,200),Rn(200),WINDa(200),WINDa2(200)
      real X(200,200),Y(200,200),WINDr1(200,200),WINDr2(200,200)
      real Xn1(200,200),Yn1(200,200),Xn2(200,200),Yn2(200,200)
      real XL1(200),XL2(200),YL1(200),YL2(200)
      real LXn1(5832),LXn2(5832),LYn1(5832),LYn2(5832)
      real LWINDr1(5832),LWINDr2(5832),Xnmin(72),Ynmin(72)
      real LWIND1(32400),LWIND2(32400),LWr1(81),LWr2(81),LR1(81)
      real LR2(81),LR11(81),LR22(81),Vm1(72),Vm2(72),Rm1(72),Rm2(72)

      character*3,T1,T2,T3,T4
      character*20 NAME1,NAME2
      character*30 HEADER,JUNK
      
      pi=3.1415926
c----------------------------------------------------------------
c       Read HRD wind at time 1, LAT_N1, LON_N1 are the number 
c       of observation points in latitude and longitude direction
c----------------------------------------------------------------
      open(11,file=NAME1,status='old')
      read(11,'(A30)')HEADER
      read(11,21)HEADER,delx1,JUNK
      read(11,22)HEADER,CC(1),JUNK,CC(2)
      read(11,'(A30)')HEADER
      read(11,*)LON_N1
      LAT_N1=LON_N1

      Nline=(int(LON_N1/6)+1)*4+8
      do i=1,Nline
         read(11,'(A30)')HEADER
      enddo
      
      do 13 i=1,LAT_N1
         do 12 j=1,int(LON_N1/2)
            read(11,23)T1,U1(i,2*j-1),T2,V1(i,2*j-1),
     *                 T3,U1(i,2*j),T4,V1(i,2*j)
 12      continue
         read(11,123)T1,U1(i,LON_N1),T2,V1(i,LON_N1)
 13   continue
      close(11)
      do 15 i=1,LAT_N1
         do 14 j=1,LON_N1
            WIND1(i,j)=sqrt(U1(i,j)**2+V1(i,j)**2)
 14      continue
 15   continue
      DIFN1=(LAT_N1+1)/2
      do 17 i=1,LAT_N1
         do 16 j=1,LON_N1
            X1(i,j)=(j-DIFN1)*delx1
            Y1(i,j)=(i-DIFN1)*delx1
 16      continue
 17   continue

 21   format(A6,1x,F7.5,1x,A11)
 22   format(A22,1x,F8.4,1x,A19,1x,F8.4)
 23   format(A1,1x,F12.7,A1,1x,F12.7,A2,1x,F12.7,A1,1x,F12.7)
 123  format(A1,1x,F12.7,A1,1x,F12.7)

c-------------------------------------------------------------
c       Read HRD wind at time 2, LAT_N2, LON_N2 are the obs.
c       points in latitude and longitude directions
c-------------------------------------------------------------
      open(31,file=NAME2,status='old')
      read(31,'(A30)')HEADER
      read(31,21)HEADER,delx2,JUNK
      read(31,22)HEADER,CC(3),JUNK,CC(4)
      read(31,'(A30)')HEADER
      read(31,*)LON_N2
      LAT_N2=LON_N2

      Nline=(int(LON_N2/6)+1)*4+8
      do i=1,Nline
         read(31,'(A30)')HEADER
      enddo
      
      do 33 i=1,LAT_N2
         do 32 j=1,int(LON_N2/2)
            read(31,23)T1,U2(i,2*j-1),T2,V2(i,2*j-1),
     *                 T3,U2(i,2*j),T4,V2(i,2*j)
 32      continue 
         read(31,123)T1,U2(i,LON_N2),T2,V2(i,LON_N2)
 33   continue
      close(31)

      do 35 i=1,LAT_N2
         do 34 j=1,LON_N2
            WIND2(i,j)=sqrt(U2(i,j)**2+V2(i,j)**2)
 34      continue
 35   continue


      DIFN2=(LAT_N2+1)/2
      do 37 i=1,LAT_N2
         do 36 j=1,LON_N2
            X2(i,j)=(j-DIFN2)*delx2
            Y2(i,j)=(i-DIFN2)*delx2
 36      continue
 37   continue
      XMM(1)=X1(LAT_N1,LON_N1)
      XMM(2)=Y1(LAT_N1,LON_N1)
      XMM(3)=X2(LAT_N2,LON_N2)
      XMM(4)=Y2(LAT_N2,LON_N2)
C----------------------------HRD wind field reading finished
 
c----------------------------------------------------------
c       pick up the minimum distance from the boundary to 
c       the center of the hurricane for the two observation
c       field at time 1 and 2. and use it to generate the
c       grids in radial directions for every 5 degrees
c----------------------------------------------------------
      XB=int(min(XMM(1),XMM(2),XMM(3),XMM(4)))
	
      do 38 i=1,72
	  theta(i)=-pi/4.0+pi*(i-1.0)/36.0
 38    continue
 
      DO 40 j=1,18
         XMAX=XB*tan(theta(j))
         YMAX=XB
         do 39 i=1,81
            X(i,j)=XMAX*i/81
            Y(i,j)=YMAX*i/81
 39      continue
 40   continue

      DO 42 j=19,36
         XMAX=XB
         YMAX=XB/tan(theta(j))
         do 41 i=1,81
            X(i,j)=XMAX*i/81
            Y(i,j)=YMAX*i/81
 41      continue
 42     continue
      DO 44 j=37,54
         YMAX=-XB;
         XMAX=-XB*tan(theta(j))         
         do 43 i=1,81
            X(i,j)=XMAX*i/81
            Y(i,j)=YMAX*i/81
 43      continue
 44   continue

      DO 46 j=55,72
         XMAX=-XB
         YMAX=-XB/tan(theta(j))
         do 45 i=1,81
            X(i,j)=XMAX*i/81
            Y(i,j)=YMAX*i/81
 45      continue
 46   continue

	
C--------------------------- set up X, Y OK

	do i=1,LAT_N1
	   YL1(i)=Y1(i,1)
	enddo
	do j=1,LON_N1
	   XL1(j)=X1(1,j)
	enddo	
        call splie2(YL1,XL1,WIND1,LAT_N1,LON_N1,y2a1)

	do 48 i=1,81
          do 47 j=1,72
             call splin2(YL1,XL1,WIND1,y2a1,LAT_N1,LON_N1,Y(i,j),
     *                       X(i,j),WINDr1(i,j))
47         continue
48      continue
	
	do i=1,LAT_N2
	   YL2(i)=Y2(i,1)
	enddo
	do j=1,LON_N2
	   XL2(j)=X2(1,j)
	enddo	

        call splie2(YL2,XL2,WIND2,LAT_N2,LON_N2,y2a2)
        do 88 i=1,81
          do 87 j=1,72
             call splin2(YL2,XL2,WIND2,y2a2,LAT_N2,LON_N2,Y(i,j),
     *                       X(i,j),WINDr2(i,j))
87         continue
88      continue

c------------------------------------------------------------------
c       finishing interp WIND1, WIND2 onto X,Y to WINDr1,WINDr2 
c-------------------------------------------------------------------
        do 50 j=1,72
           do i=1,81
              LWr1(i)=WINDr1(i,j)
              LR1(i)=sqrt(X(i,j)**2+Y(i,j)**2)
              LR11(i)=LR11(i)
           enddo
           call shell(81,LWr1,LR1,LR11)
           Vm1(j)=LWr1(81)
           Rm1(j)=LR1(81)
 50     continue

        do 51 j=1,72
           do i=1,81
              LWr2(i)=WINDr2(i,j)
              LR2(i)=sqrt(X(i,j)**2+Y(i,j)**2)
              LR22(i)=LR22(i)
           enddo
           call shell(81,LWr2,LR2,LR11)
           Vm2(j)=LWr2(81)
           Rm2(j)=LR2(81)
 51     continue

        AveRm1=0
        AveRm2=0
        do i=1,72
           AveRm1=AveRm1+Rm1(i)
           AveRm2=AveRm2+Rm2(i)
        enddo
        AveRm1=AveRm1/72
        AveRm2=AveRm2/72


        do 57 j=1,72
           Nright=0.0
           Nleft=0.0
           if (abs(Rm2(j)-AveRm2) .gt. 0.3*AveRm2) then
              do 52 i=j,j+5
                 if (i.gt.72)then
 		     i3=i-72
		 else
	             i3=i
		 endif
                 if (abs(Rm2(i3)-AveRm2).lt.0.3*AveRm2) then
                    Nright=i3
                    goto 53
                 endif
 52           continue
 53           continue
              do 54 i=j,j-5,-1
                 if (i.lt.1) then
		     i3=i+72
		 else
		     i3=i
		 endif
                 if (abs(Rm2(i3)-AveRm2).lt. 0.3*AveRm2) then
                    Nleft=i3
                    goto 55
                 endif
 54           continue
 55           continue

              if(Nright.ne.0.0.and.Nleft.ne.0.0)then
                 if (Nright.gt.Nleft)then
                    do i=Nleft+1,Nright-1
                       tt=(Rm2(Nright)-Rm2(Nleft))*(i-Nleft)
     *                    /(Nright-Nleft+1)+Rm2(Nleft)
                       Rm2(i)=tt
                    enddo
                 else
                    k1=0
                    do i=Nleft+1,72
                       k1=k1+1
                       ID(k1)=i
                    enddo
                    do i=1,Nright-1
                       k1=k1+1
                       ID(k1)=i
                    enddo

                    do 56 i=1,k1
                       tt=(Rm2(Nright)-Rm2(Nleft))*i/(k1+1)+Rm2(Nleft)
                       Rm2(ID(i))=tt
 56                 continue
                 endif
              endif
           endif
 57     continue
 
        do 67 j=1,72
           Nright=0.0
           Nleft=0.0
           if (abs(Rm1(j)-AveRm1) .gt. 0.3*AveRm1) then
              do 62 i=j,j+5
                 if (i.gt.72) then
		     i3=i-72
		 else
		     i3=i
		 endif
                 if (abs(Rm1(i3)-AveRm1).lt.0.3*AveRm1) then
                    Nright=i3
                    goto 63
                 endif
 62           continue
 63           continue
              do 64 i=j,j-5,-1
                 if (i.lt.1) then
		     i3=i+72
		 else
	 	     i3=i
		 endif
                 if (abs(Rm1(i3)-AveRm1).lt. 0.3*AveRm1) then
                    Nleft=i3
                    goto 65
                 endif
 64           continue
 65           continue

              if(Nright.ne.0.0.and.Nleft.ne.0.0)then
                 if (Nright.gt.Nleft)then
                    do i=Nleft+1,Nright-1
                       tt=(Rm1(Nright)-Rm1(Nleft))*(i-Nleft)
     *                    /(Nright-Nleft+1)+Rm1(Nleft)
                       Rm1(i)=tt
                    enddo
                 else
                    k1=0
                    do i=Nleft+1,72
                       k1=k1+1
                       ID(k1)=i
                    enddo
                    do i=1,Nright-1
                       k1=k1+1
                       ID(k1)=i
                    enddo

                    do 66 i=1,k1
                       tt=(Rm1(Nright)-Rm1(Nleft))*i/(k1+1)+Rm1(Nleft)
                       Rm1(ID(i))=tt
 66                 continue
                 endif
              endif
           endif
 67     continue

        do 69 i=1,81
           do 68 j=1,72
              Xn1(i,j)=X(i,j)/Rm1(j)
              Yn1(i,j)=Y(i,j)/Rm1(j)
              Xn2(i,j)=X(i,j)/Rm2(j)
              Yn2(i,j)=Y(i,j)/Rm2(j)
 68        continue
 69     continue


	do 70 j=1,72
	   RR1=sqrt(Xn1(81,j)**2+Yn1(81,j)**2)
	   RR2=sqrt(Xn2(81,j)**2+Yn2(81,j)**2)
	   if(RR1.le.RR2) then
	       Xnmin(j)=Xn1(81,j)
	       Ynmin(j)=Yn1(81,j)
	   else
	       Xnmin(j)=Xn2(81,j)
	       Ynmin(j)=Yn2(81,j)
	   endif
70	continue

        

      DO 74 j=1,72
         do 73 i=1,101
            XI(i,j)=Xnmin(j)*i/101.0
            YI(i,j)=Ynmin(j)*i/101.0
  73      continue
 74   continue

	do 77 j=1,72
	   do 75 i=1,81
	      Rn(i)=sqrt(Xn1(i,j)**2+Yn1(i,j)**2)
	      WINDa(i)=WINDr1(i,j)
75	   continue
	   call spline(Rn,WINDa,81,1E30,1E30,WINDa2)
	   do 76 k=1,101 
		RI=sqrt(XI(k,j)**2+YI(k,j)**2)
	   	call splint(Rn,WINDa,WINDa2,81,RI,WIND1itp(k,j))
76	   continue
77	continue

        do 80 j=1,72
           do 78 i=1,81
              Rn(i)=sqrt(Xn2(i,j)**2+Yn2(i,j)**2)
              WINDa(i)=WINDr2(i,j)
78         continue
           call spline(Rn,WINDa,81,1E30,1E30,WINDa2)
           do 79 k=1,101
                RI=sqrt(XI(k,j)**2+YI(k,j)**2)
                call splint(Rn,WINDa,WINDa2,81,RI,WIND2itp(k,j))
79         continue
80      continue

	return
	end




	SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
	INTEGER j,k,m,n,NN
	REAL x1a(200),x2a(200),y2a(200,200),ya(200,200)
	PARAMETER (NN=200) 
c	Maximum expected value of n and m.
c--------- USES spline
c	Given an m by n tabulated function ya(1:m,1:n), and tabulated independent variables
c	x2a(1:n), this routine constructs one-dimensional natural cubic splines of the rows of ya
c	and returns the second-derivatives in the array y2a(1:m,1:n). (The array x1a is included
c	in the argument list merely for consistency with routine splin2.)
c-------------------------------------------------------------------------------------------------------------------
	REAL y2tmp(n),ytmp(n)

	do 13 j=1,m
	do 11 k=1,n
	     ytmp(k)=ya(j,k)
11	continue
	call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp) 
c--------Values 1×1030 signal a natural spline.
	do 12 k=1,n
	     y2a(j,k)=y2tmp(k)
12	continue
13	continue

	return
	END


	SUBROUTINE spline(x,y,n,yp1,ypn,y2)
	INTEGER n,NMAX
	REAL yp1,ypn,x(200),y(200),y2(200)
	PARAMETER (NMAX=500)
c	Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
c	x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
c	function at points 1 and n, respectively, this routine returns an array y2(1:n) of
c	length n which contains the second derivatives of the interpolating function at the tabulated
c	points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set
c	the corresponding boundary condition for a natural spline, with zero second derivative on
c	that boundary.
c	Parameter: NMAX is the largest anticipated value of n.

	INTEGER i,k
	REAL p,qn,sig,un,u(NMAX)

c	The lower boundary condition is set either to be natural 
c	or else to have a specified first derivative.

	if (yp1.gt..99e30) then 
	     y2(1)=0.
	     u(1)=0.
	else 
	     y2(1)=-0.5
	     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif

c	This is the decomposition loop of the tridiagonal
c	algorithm. y2 and u are used for temporary
c	storage of the decomposed factors.

	do 11 i=2,n-1 
	    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	    p=sig*y2(i-1)+2.
	    y2(i)=(sig-1.)/p
	    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *            /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11	continue
c	The upper boundary condition is set either to be natural 
c	or else to have a specified first derivative.
	if (ypn.gt..99e30) then 
	    qn=0.
	    un=0.
	else 
	    qn=0.5
	    un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

c	This is the backsubstitution loop of the tridiagonal algorithm. 
	do 12 k=n-1,1,-1 
	     y2(k)=y2(k)*y2(k+1)+u(k)
12	continue
	
	return
	END


	SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)

	INTEGER m,n,NN
	REAL x1,x2,y,x1a(200),x2a(200)
	real y2a(200,200),ya(200,200)
	
	PARAMETER (NN=200) 
c	Maximum expected value of n and m.
c	C USES spline,splint
c	Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by that routine;
c	and given a desired interpolating point x1,x2; this routine returns an interpolated function
c	value y by bicubic spline interpolation.

	INTEGER j,k
	REAL y2tmp(200),ytmp(200),yytmp(200)

c	Perform m evaluations of the row splines constructed by splie2, 
c	using the onedimensional spline evaluator splint.
 
	do 12 j=1,m 
	    do 11 k=1,n
                  ytmp(k)=ya(j,k)
                  y2tmp(k)=y2a(j,k)
11	    continue
	    call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12	continue

	    call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp) 
	    call splint(x1a,yytmp,y2tmp,m,x1,y)

	return
	END


	SUBROUTINE splint(xa,ya,y2a,n,x,y)
	INTEGER n
	REAL x,y,xa(200),y2a(200),ya(200)

c	Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
c	xai's in order), and given the array y2a(1:n), which is the output from spline above,
c	and given a value of x, this routine returns a cubic-spline interpolated value y.

	INTEGER k,khi,klo
	REAL a,b,h

c 	We will find the right place in the table by means of bisection.
c	This is optimal if sequential calls to this routine are at random
c	values of x. If sequential calls are in order, and closely
c	spaced, one would do better to store previous values of
c	klo and khi and test if they remain appropriate on the
c	next call.

	klo=1 
	khi=n
1	if (khi-klo.gt.1) then
	     k=(khi+klo)/2
	     if(xa(k).gt.x)then
  	          khi=k
	     else
	          klo=k
	     endif
	goto 1
	endif
 
c	klo and khi now bracket the input value of x.

	h=xa(khi)-xa(klo)

c	bad xa input in splint. The xa's must be distinct.
	if (h.eq.0.)then
           print*,'n= ',n
	   print*,khi,xa(khi),klo,xa(klo)
	   print*,xa 
	   pause
	endif 
c	Cubic spline polynomial is now evaluated.
	a=(xa(khi)-x)/h 
	b=(x-xa(klo))/h
	y=a*ya(klo)+b*ya(khi)+
     *          ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

	return
	END


	SUBROUTINE shell(n,a,b,c)
	INTEGER n,tt
	REAL a(n),b(n),c(n),aa(n),bb(n),cc(n)
c--------------------------------------------------------------------------------
c	Sorts an array a(1:n) into ascending numerical order by Shell's method 
c	(diminishing increment sort). n is input; a is replaced on output by 
c	its sorted rearrangement.
c--------------------------------------------------------------------------------
	INTEGER i,j,inc,kk,i1,j1
	REAL v,v1,v3,b1,c1

	inc=1 
1 	continue
	inc=3*inc+1
	if(inc.le.n)goto 1
2	continue 
	   inc=inc/3
	   do 11 i=inc+1,n 
	      v=a(i)
	      v1=b(i)
	      v3=c(i)
	      j=i
3 	      continue
	      if(a(j-inc).gt.v)then 
	         a(j)=a(j-inc)
		 b(j)=b(j-inc)
		 c(j)=c(j-inc)
	         j=j-inc
	         if(j.le.inc)goto 4
	      goto 3
	      endif
4 	      continue
	      a(j)=v
	      b(j)=v1
	      c(j)=v3
11	    continue
	if(inc.gt.1)goto 2

	kk=1
	do 25 i=2,n
	  aa(kk)=a(i-1)
	  bb(kk)=b(i-1)
	  cc(kk)=c(i-1)
	  kk=kk+1
	  if ((a(i)-a(i-1)).lt.1E-3) goto 25

c	  aa(kk)=a(i)
c	  bb(kk)=b(i)
c	  cc(kk)=c(i)
	  kk=kk-1
	  do 18 j=2,kk
	     b1=bb(j)
	     c1=cc(j)
	     do 17 i1=j-1,1,-1
	        if(bb(i1).le.b1) goto 16
		bb(i1+1)=bb(i1)
		cc(i1+1)=cc(i1)
17	     continue
	     i1=0
16	     continue
	     bb(i1+1)=b1
	     cc(i1+1)=c1
18	   continue
	   tt=0
	   do j1=i-kk,i-1
	        tt=tt+1
		b(j1)=bb(tt)
		c(j1)=cc(tt)
	   enddo
	kk=1
25	continue

	return
	END



	subroutine Interp_linear(m,n,x,y,z,mi,ni,xi,yi,zi)

c---------------------------------------------------------------
c put all data into one column and sort them in increasing order
c---------------------------------------------------------------
	integer i,j,k,m,n,mi,ni
	real x(m,n),y(m,n),z(m,n),xi(mi,ni),yi(mi,ni),zi(mi,ni)
	real X1(m*n),Y1(m*n), Z1(m*n)	
	k=0
	do 15 i=1,m
	   do 14 j=1,n
	      k=k+1
   	      X1(k)=x(i,j)
	      Y1(k)=y(i,j)
	      Z1(k)=z(i,j)
14	   continue
15	continue
	call shell(m*n,Y1,X1,Z1)
	
	call linear(m*n,X1,Y1,Z1,mi,ni,xi,yi,zi)

c        call cubic(m*n,X1,Y1,Z1,mi,ni,xi,yi,zi);
  
	return	
	end

	subroutine linear(len,x,y,z,mi,ni,xi,yi,zi)

	integer len,i,j,mi,ni,hal,i1,i2,k
	real x(len),y(len),z(len),xi(mi,ni)
	real yi(mi,ni),zi(mi,ni),ttt
	real s1,s2,s3,s4,p1,p2,p3,p4,x1,x2,x3,x4,y1,y2,y3,y4

	hal=int(len/2)

	p1=0.0
	p2=0.0
	p3=0.0
	p4=0.0
	s1=0.0
	s2=0.0
	s3=0.0
	s4=0.0
	do 21 i=1,mi
	   do 20 j=1,ni
	     if(yi(i,j).lt.y(hal)) then
		    do 19  k=1,hal		       
		       if(yi(i,j).ge.y(k).and.yi(i,j).le.y(k+1)) then
			 do 18 i1=k,2,-1
		           if(xi(i,j).ge.x(i1-1).and.xi(i,j).le.x(i1))then
			     p1=z(i1-1)
			     p2=z(i1)
			     x1=x(i1-1)
			     x2=x(i1)
			     y1=y(i1-1)
			     y2=y(i1)
			     goto 33
			   endif
18			 continue
 33			 continue
			 do 111 i2=k+1,len
			    if(xi(i,j).ge.x(i2).and.xi(i,j).le.x(i2+1))then
			     p3=z(i2)
			     p4=z(i2+1)
			     x3=x(i2)
			     x4=x(i2+1)
			     y3=y(i2)
			     y4=y(i2+1)
			     goto 44
			   endif
 111			continue
 44			continue
		     endif
19		   continue
	      else
                do 15  k=hal,len
                   if(yi(i,j).ge.y(k).and.yi(i,j).le.y(k+1)) then
			 do 14 i1=k,2,-1
		           if(xi(i,j).ge.x(i1-1).and.xi(i,j).le.x(i1))then
			     p1=z(i1-1)
			     p2=z(i1)
			     x1=x(i1-1)
			     x2=x(i1)
			     y1=y(i1-1)
			     y2=y(i1)
			     goto 55
			   endif
 14			continue
 55			continue
			 do 13  i2=k+1,len
			    if(xi(i,j).ge.x(i2-1).and.xi(i,j).le.x(i2))then
			     p3=z(i2-1)
			     p4=z(i2)
			     x3=x(i2-1)
			     x4=x(i2)
			     y3=y(i2-1)
			     y4=y(i2)
			     goto 66
			   endif
 13  			continue
 66			continue
		     endif
 15	        continue
	      endif
	      
	      	s1=sqrt((x1-xi(i,j))**2+(y1-yi(i,j))**2)
	      	s2=sqrt((x2-xi(i,j))**2+(y2-yi(i,j))**2)
	      	s3=sqrt((x3-xi(i,j))**2+(y3-yi(i,j))**2)
	      	s4=sqrt((x4-xi(i,j))**2+(y4-yi(i,j))**2)
	      
	      zi(i,j)=(p1*s1+p2*s2+p3*s3+p4*s4)/(s1+s2+s3+s4)

20	   continue
21	continue

	return
	end


      subroutine Rad2Rec(m,n,Xd1,Yd1,WINDd1,X11,Y11,WIND11)
      integer i,j,k,k1,i1,j1,m,n,l
      real Xd1(182,72),Yd1(182,72),WINDd1(200,200),ww(200)
      real X11(400,400),Y11(400,400),WIND11(400,400)
      real theta,num,R,W1,W2,pi,Rn(200),WW2(200)

      pi=3.1415926
      do 21 i=1,m
         do 20 j=1,n
	    theta=atan(x11(i,j)/y11(i,j))
	    if (abs(y11(i,j)).lt.1E-5.and.x11(i,j).ge.0.0) then
	       theta=pi/2
	    else if(abs(y11(i,j)).lt.1E-5.and.x11(i,j).lt.0.0) then
	       theta=pi*3/2
  	    endif

	    if(y11(i,j).gt.0.0.and.theta.lt.-pi/4) then
		theta=theta+2*pi
	    else if(y11(i,j).lt.0.0) then
		theta=theta+pi
	    endif
            num=(theta+pi/4.)/(pi/36.)+1.
            k=int(num)
            do 19 i1=1,182
               Rn(i1)=sqrt(Xd1(i1,k)**2+Yd1(i1,k)**2)
               WW(i1)=WINDd1(i1,k)
 19         continue
            R=sqrt(x11(i,j)**2+y11(i,j)**2)
            call spline(Rn,WW,182,1E30,1E30,WW2)
            call splint(Rn,WW,WW2,182,R,W1)
	    if(W1.lt.0)then
	        W1=0.0
	    endif
	    k1=k+1
	    if(k1.gt.72) then
	       k1=1
	    endif
            do 18 i1=1,182
               Rn(i1)=sqrt(Xd1(i1,k1)**2+Yd1(i1,k1)**2)
               WW(i1)=WINDd1(i1,k1)
 18         continue
            call spline(Rn,WW,182,1E30,1E30,WW2)
            call splint(Rn,WW,WW2,182,R,W2)
	    if(W2.lt.0)then
	       W2=0.0
	    endif
            WIND11(i,j)=W1*(k+1-num)+W2*(num-k)
	    if(R.gt.Rn(101).and.WIND11(i,j).gt.WW(101)) then
	       WIND11(i,j)=WIND11(i-1,j-1)
	    endif
	    
 20      continue
 21   continue
	
      return
      end

      

	SUBROUTINE sort(n,arr)	
	INTEGER n,M,NSTACK
	REAL arr(n)
	PARAMETER (M=7,NSTACK=50)
c	Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
c	is input; arr is replaced on output by its sorted rearrangement.
c	Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
c	auxiliary storage.
	INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
	REAL a,temp
	jstack=0
	l=1
	ir=n
c	 Insertion sort when subarray small enough.
1	 if(ir-l.lt.M)then
	      do 12 j=l+1,ir
	            a=arr(j)
	            do 11 i=j-1,l,-1
	                  if(arr(i).le.a)goto 2
	                  arr(i+1)=arr(i)
11	            continue
	            i=l-1
2 	        arr(i+1)=a
12	     continue
	     if(jstack.eq.0)return
	     ir=istack(jstack) 
c            Pop stack and begin a new round of partitioning.
	     l=istack(jstack-1)
	     jstack=jstack-2
	else
	     k=(l+ir)/2 
c	Choose median of left, center, and right elements as partitioning
c	     element a. Also rearrange so that a(l) ¡Ü a(l+1) ¡Ü a(ir).
	     temp=arr(k)
	     arr(k)=arr(l+1)
	     arr(l+1)=temp
	     if(arr(l).gt.arr(ir))then
	          temp=arr(l)
	          arr(l)=arr(ir)
	          arr(ir)=temp
	     endif
	     if(arr(l+1).gt.arr(ir))then
	          temp=arr(l+1)
	          arr(l+1)=arr(ir)
	          arr(ir)=temp
	     endif
	     if(arr(l).gt.arr(l+1))then
	          temp=arr(l)
	          arr(l)=arr(l+1)
	          arr(l+1)=temp
	     endif
	     i=l+1 
c	     Initialize pointers for partitioning.
	     j=ir
	     a=arr(l+1) 
c            Partitioning element.
3 	     continue 
c	     Beginning of innermost loop.
	     i=i+1 
c	     Scan up to find element > a.
	     if(arr(i).lt.a)goto 3
4 	     continue
	      j=j-1 
c	      Scan down to find element < a.
	     if(arr(j).gt.a)goto 4
	     if(j.lt.i)goto 5 
c	     Pointers crossed. Exit with partitioning complete.
	     temp=arr(i) 
c	     Exchange elements.
	     arr(i)=arr(j)
	     arr(j)=temp
	     goto 3 
c	     End of innermost loop.
5	     arr(l+1)=arr(j) 
c	     Insert partitioning element.
	     arr(j)=a
	      jstack=jstack+2
c	Push pointers to larger subarray on stack, process smaller subarray immediately.
	     if(jstack.gt.NSTACK)pause 
c	      ¡¯NSTACK too small in sort¡¯
	     if(ir-i+1.ge.j-l)then
	          istack(jstack)=ir
	          istack(jstack-1)=i
	          ir=j-1
	      else
	           istack(jstack)=j-1
	           istack(jstack-1)=l
	           l=i
	      endif
	endif
	goto 1

	END
