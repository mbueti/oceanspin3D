      subroutine verZ2SIGsp(n,ni,x,y,xi,yi,i,j)
      real x(n),y(n),xi(ni),yi(ni),R1(500)
c-----
c----- use SP to interpolate from Z-level to sigma
c-----  after interpolation t can be unstable
c----- x(l) l=1,n    Z-levels
c----- y(l) is t or s on Z-levels
c-----   
c----- xi(k)=-ZZ(k)*(H+EL) k=1,ni  Z at sigma-levels
c-----     xi(ni)=H+EL
c----- for GDEM x(1)=0 x(2)=10
c----- for model  0 < xi(1)=-ZZ(1)*(H+EL) < -5500*ZZ(1)=5
c-----   Thus   x(1) < xi(1)
c----- yi(k) k=1,ni  interpolated values of t or s to sigma-lvels
c
c-----------------
c      if(i.eq.203.and.j.eq.48) then
c       print *,'check verZ2SIGsp i,j=',i,j
c       print *,'(x(l),l=1,n)'
c       write(6,501) (x(l),l=1,n)
c       print *,'(y(l),l=1,n)'
c       write(6,502) (y(l),l=1,n)
c      end if
c-----------------
      h=xi(ni)
c-----    find l=lb that x(lb) > h
c------------- falk 04-21-02 set lb=n
      lb=n
      do l=1,n
       if(x(l).gt.h) then
        lb=l
        go to 200
       end if
      end do
 200  continue
c-----------------
c      if(i.eq.203.and.j.eq.48) then
c       print *,'check verZ2SIGsp i,j=',i,j
c       print *,' i,j,n,ni,lb,h=',i,j,n,ni,lb,h
c       print *,'(x(l),l=1,lb)'
c       write(6,501) (x(l),l=1,lb)
c       print *,'(y(l),l=1,lb)'
c       write(6,502) (y(l),l=1,lb)
c      end if
c-----------------
c-----  sharp saved additional level under bottom
c-----  thus  though Zlev(lb) > h there is information in y(lb)
c
c-----    use SP to interpolate from Z-levels where z:
c-----         0=x(1) <= z <= h < x(lb)=Zlev(lb) to
c-----    sigma-levels for z: 
c-----       0 < xi(1) < z < xi(ni-1) < h
c
c---------- falk 06-21-05 use R1 as inner param
      call SP(x,y,yi,xi,lb,ni-1,n,ni)
c
c-----    send to ni sigma level a value from a previous level
c-----    Frolov sends t to the last sigma-level from previous 
c          
      yi(ni)=yi(ni-1)
c------------- for debug
c      if(i.eq.175.and.j.eq.25.or.i.eq.100.and.j.eq.91) then
c       print *,' verZ2SIGsp'
c       print *,' i,j,n,ni,lb,h=',i,j,n,ni,lb,h
c       print *,'(x(l),l=1,lb)'
c       write(6,501) (x(l),l=1,lb)
c       print *,'(y(l),l=1,lb)'
c       write(6,502) (y(l),l=1,lb)
 501    format(10f7.0)
 502    format(10f7.2)
c       print *,'(xi(k),k=1,ni)'
c       write(6,501) (xi(k),k=1,ni)
c       print *,'(yi(k),k=1,ni)'
c       write(6,502) (yi(k),k=1,ni)
c      end if
c-------------
      return
      end
