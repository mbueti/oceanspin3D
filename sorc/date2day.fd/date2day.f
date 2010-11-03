      program date_to_day
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM: GFDL_DATE2DAY
C   PRGMMR: MARCHOK          ORG: NP22        DATE: 2001-03-22
C
C ABSTRACT: this program converts a date integer in YYMMDDHH format to
C   julian day in year YY
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

c     this program converts a date integer in YYMMDDHH format to 
c     julian day in year YY
c
c     dail - march 21, 2000

      implicit none

      integer*8 date
      integer dat2day(12),dat2dayl(12),day,month,year,hour,n
      real julday

c number of days in each month in regular and leap years
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

c      CALL W3TAGB('GFDL_DATE2DAY',2001,0081,0096,'NP22')                  

      read(5,*) date

c get year and month
      year=int(date/1000000.)
      month=nint(100*(date/1000000.-int(date/1000000.)))

c calculate julian day in year.  account for leap years.
      julday=0
      if(mod(year,4).eq.0) then
         do n=1,month-1
            julday=julday+dat2dayl(n)
         end do
      else
         do n=1,month-1
            julday=julday+dat2day(n)
         end do
      end if

      julday=julday+nint(100*(date/10000.-int(date/10000.)))
      hour=date-nint(date/100.)*100
      julday=julday+float(hour)/24.

      write(61,*) julday, year

c      call w3tage('GFDL_DATE2DAY')

      stop
      end
