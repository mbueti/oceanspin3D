      program day_to_date
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: GFDL_DAY2DATE
C   PRGMMR: MARCHOK          ORG: NP22        DATE: 2001-03-22
C
C ABSTRACT: this program converts year and julian date input to date 
C   in YYMMDDHH format date integer
C
C PROGRAM HISTORY LOG:
C   98-06-02  frolov - original implementation
C   00-03-21  rowe   - reorganized
C
C INPUT FILES:
C   UNIT    5    standart input
C
C OUTPUT FILES:
C   UNIT    6    standart output
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C
C$$$

c     convert year and julian date input to date in YYMMDDHH format
c
c     dail - March 21, 2000

      implicit none

      integer*8 date
      integer dat2day(12),dat2dayl(12),day,month,year,year1,hour,n
      real julday,julday1

c     number of days in months of leap and regular years
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

c      CALL W3TAGB('GFDL_DAY2DATE',2001,0081,0096,'NP22')                  

      read(5,*) year,julday

      if(int(julday).gt.365+int(1./(mod(year,4)*100.+1.))) then
         julday1=julday-365-int(1./(mod(year,4)*100.+1.))
         year1=year+1
      else
         julday1=julday
         year1=year
      end if

      day=0
      n=1
      if(mod(year1,4).eq.0) then  
         do while(day+dat2dayl(n).lt.int(julday1))
            day=day+dat2dayl(n)
            n=n+1   
         end do
      else
         do while(day+dat2day(n).lt.int(julday1))
            day=day+dat2day(n)
            n=n+1   
         end do
      end if

      month=n
      day=int(julday1-day)
      hour=int((julday1-int(julday1))*24.)
      date=year1*1000000+month*10000+day*100+hour

      write(61,*) date
c
c      CALL W3TAGE('GFDL_DAY2DATE')
c
      stop
      end
