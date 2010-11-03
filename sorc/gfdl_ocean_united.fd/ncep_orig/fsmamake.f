      program fsmamake
c---------  
c--------- read global data (percent of water over grid point),
c-------   calculate FSMA(im+1,jm+1) mask for atm from pctw
c-------   the 1-st atm point 1/12deg to west (south)
c-------              from  the 1-st T-point in ocean       A  V  A
c-------   the last atm point 1/12deg to east (north)       U  T  U
c-------              from  the last T-point in ocean       A  V  A
c-------
c-------     A-atmos points, T,U,V-ocean temp, and vel. points
c-------    FSMA=1 over water,FSMA =0 over land
c-------    FSMA(im+1,jm+1) included in common/maskat/ in comblk.h
c
      parameter (im=254,jm=225,LON=2160,LAT=1080,ima=293,jma=227)

c----------- (im,jm) size of an ocean region with
c----------- a resolution 1/6 deg
c----------- (LON,LAT) size of global pctwat which is defined with
c-----------   a step 1/6 deg for centers with 1/12deg to the right
      dimension pctw(LON,LAT),XG(LON),YG(LAT)
      dimension XA(ima),YA(jma)
      dimension FSMA(ima,jma)
c------- 01-09-04 read H and FSM and write them together with FSMA
      dimension H(im,jm),FSM(im,jm)
      real  LONGMIN,LONGMAX,LATMIN,LATMAX,MSKVAL
      character*80 FTPCTW,FTFSMA,FTHFSM,FTHFSMA
      FTPCTW='/nwprod/fix/gfdl_pctwat'
c------- 01-09-04 read H and FSM and write them together with FSMA
c     FTFSMA='/emc1/wx20af/levitus/FSMAu'
      FTHFSM='/emc1/wx20af/gfdl/curr_initdata/gfdl_Hdeepgsu'
      FTHFSMA='/emc1/wx20af/gfdl/curr_initdata/gfdl_H_FSM_FSMAu'
c------- 01-09-04 read H and FSM and write them together with FSMA
      OPEN(47,file=FTPCTW,form='unformatted')
      OPEN(67,file=FTHFSM,form='formatted')
      OPEN(66,file=FTHFSMA,form='formatted')
c
      MSKVAL=-9.
c
c---------- specify bnd of ocean region
c
      LONGMIN=-98.5
      LONGMAX=-50.
      LATMIN=10.
      LATMAX=47.5
c--------------
      write(6,106) LONGMIN,LONGMAX,LATMIN,LATMAX
 106  format('LONGMIN,LONGMAX,LATMIN,LATMAX=',4f9.3)
c--------------
      dlon=(LONGMAX-LONGMIN)/float(im-1)
      dlat=(LATMAX-LATMIN)/float(jm-1)
      write(6,107) dlon,dlat
 107  format('         dlon,dlat=',2f9.4)
c
c--------- process pctwat, define RESN,RESNO
c
      RESN=1./6.
      RESNO=0.5*RESN
c
      do j=1,LAT
       YG(j)=float(j-1)*RESN+RESNO-90.
      end do
      do i=1,LON
       XG(i)=float(i-1)*RESN+RESNO
      end do
c
c-------  read global PCTWAT data percent of water over pixel
c
      rewind(47)
      do j=1,LAT
       read(47) (pctw(i,j),i=1,LON)
      end do
      close(47)
c
c-------- Atm mask region is larger than ocean region
c--------      for all 4 sides
c
c--------- find left low corner
      is1=(360.+LONGMIN-XG(1))/RESN+1.0001
      js1=(LATMIN-YG(1))/RESN+1.0001
c--------- find right upper corner
      ie1=(360.+LONGMAX-XG(1))/RESN+2.0001
      je1=(LATMAX-YG(1))/RESN+2.0001
c     ima=ie1-is1+1
c     jma=je1-js1+1
c------------
      print *,'is1,js1,ie1,je1,ima,jma=',is1,js1,ie1,je1,ima,jma
      write(6,210) XG(is1),XG(ie1),YG(js1),YG(je1)
 210  format('XG(is1),XG(ie1),YG(js1),YG(je1)=',4f8.3)
c------------
c
      do i=1,ima
       XA(i)=LONGMIN+RESN*float(i-1)-RESNO
      end do
      do j=1,jma
       YA(j)=LATMIN+RESN*float(j-1)-RESNO
      end do
c------------------------
      write(6,202) XA(1),XA(ima),YA(1),YA(jma)
 202  format('XA(1),XA(ima),YA(1),YA(jma)=',4f9.3)
c------------------------
c--- calculate FSMA percent of water over atmos model pixels
c
      do i=1,ima
      do j=1,jma
       FSMA(i,j)=pctw(i+is1-1,j+js1-1)
       if(FSMA(i,j).gt.50.) then
        FSMA(i,j)=1.
       else
        FSMA(i,j)=0.
       end if
      end do
      end do
c
c------- 01-09-04 read H and FSM and write them together with FSMA
c
      read(67,102) H
      read(67,102) FSM
      write(66,102) H
      write(66,102) FSM
      write(66,102) FSMA
 102  format(10f7.0)
      close (66)
      close (67)
c
      ip=178
      jp=47
      call pritop(FSMA,XA,YA,ima,jma,ip-5,ip+5,jp-10,jp+10,
     *       '     FSMA(i,j)           ',-1)
      ip=188
      jp=47
      call pritop(FSMA,XA,YA,ima,jma,ip-5,ip+5,jp-10,jp+10,
     *       '     FSMA(i,j)           ',-1)
      stop
      end
c
c---------------
c
       subroutine pritop(H,XI,YI,im,jm,ibeg,iend,j1,j2,char,mig)
       dimension H(im,jm),XI(im),YI(jm)
       character*25 char
        write(6,100) char
 100    format(/10x,a25)
        write(6,105) im,jm,ibeg,iend,j1,j2
 105    format(/1x,' PRITOP',/5x,'im,jm,ibeg,iend,j1,j2=',
     *         /5x,6i6)
        write(6,103) (i,i=ibeg,iend)
        write(6,102) (XI(i),i=ibeg,iend)
 102    format(/14x,11f7.2,/)
        do j=j1,j2
         if(mig.eq.-1) then
          jc=j2-j+j1
         else
          jc=j
         endif
         write(6,101) jc,YI(jc),(H(i,jc),i=ibeg,iend)
        end do
  101   format(i7,f7.2,11f7.2)
  103   format(14x,11i7)
        return
        end
c
c--------------
c
      subroutine INTHOR(MSKVAL,tempav,tt,XON,YON,XI,YI,idim,jdim,im,jm)
      real  tempav(idim,jdim),MSKVAL
      real  XON(idim),YON(jdim)
      real  XI(im),YI(jm)
      real  tt(im,jm),ZT(4),ZTS(4)
      ip=171
      jp=138
      RESN=(XON(2)-XON(1))
c
      do i=1,im
      do j=1,jm
c--------- find left low corner
       i1=(XI(i)+360.-XON(1))/RESN+1.0001
       j1=(YI(j)-YON(1))/RESN+1.0001
c
       x1=XI(i)+360.-XON(i1)
       x2=XON(i1+1)-360.-XI(i)
       y1=YI(j)-YON(j1)
       y2=YON(j1+1)-YI(j)
       ZT(1)=tempav(i1,j1)
       ZT(2)=tempav(i1+1,j1)
       ZT(3)=tempav(i1+1,j1+1)
       ZT(4)=tempav(i1,j1+1)
C
C******* CALCULATING FORMULA *****************************
C
       DO k=1,4
        IF(ZT(k).EQ.MSKVAL) THEN
          ZTS(k)=0.
        ELSE
          ZTS(k)=1.
        END IF
       END DO
C
       tt(I,J)=ZTS(1)*ZT(1)*X2*Y2
     1        +ZTS(2)*ZT(2)*X1*Y2
     2        +ZTS(3)*ZT(3)*X1*Y1
     3        +ZTS(4)*ZT(4)*X2*Y1
       AREA=   ZTS(1)*X2*Y2
     1        +ZTS(2)*X1*Y2
     2        +ZTS(3)*X1*Y1
     3        +ZTS(4)*X2*Y1
C
       IF(AREA.lt.1.E-6) THEN
         tt(I,J)=MSKVAL
       ELSE
         tt(I,J)=tt(I,J)/AREA
       END IF      
c-----------------
       if(i.eq.ip.and.j.eq.jp) then
        print *,' INTHOR'
        print *,'ip,jp,i1,j1=',ip,jp,i1,j1
        write(6,101) XI(i),YI(j),XON(i1),YON(j1),
     *               XON(i1+1),YON(j1+1)
 101    format('XI(i),YI(j),XON(i1),YON(j1),',
     *         'XON(i1+1),YON(j1+1)',/6f7.2)
        write(6,102) x1,x2,y1,y2
 102    format('x1,x2,y1,y2=',4f7.2)
        write(6,103) ZT
 103    format(' ZT=',4f7.2)
        write(6,104) ZTS
 104    format(' ZTS=',4f7.0)
        write(6,105) tt(i,j),AREA,RESN
 105    format('tt(i,j),AREA,RESN=',f7.2,f10.5,f7.3)
       end if
c-----------------
      end do
      end do
      return
      end
