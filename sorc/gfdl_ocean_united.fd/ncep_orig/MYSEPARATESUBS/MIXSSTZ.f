      SUBROUTINE MIXSSTZ(f,STM,H,FSM,Zlev,im,jm,nl,TBIAS)
      parameter(kk=9,jj=31,step=5.)
c---------------------------------------------------------------------
c----- this subr. was rewritten 06/17/05 by A.Falkovich
c----- to assimilate SST real data in GDEM monthly t data for Z-levels
c---------------------------------------------------------------------
c
c---- Change temperature in the layer 0-125m.
c---- Use spline to interpolate from GDEM z-levels to even step 10m
c---- Here temperature is specified at all levels 
c----                                 (under land and under bottom)
c---- Mixed layer depth is the last level n where t(1)-t(k)<0.5deg
c---- DTS=STM-t(1);  tnew(k)=t(k)+DTS for k=1,...,n
c---- if abs(DTS)>5deg correct STM
c---- if STM<t(150m) correct STM
c
      real STM(im,jm),H(im,jm),FSM(im,jm)
      real f(IM,JM,nl),Zlev(nl)
      real DTS,t(jj),z(jj),tt(kk),zz(kk),tc(jj),te(kk),dt(kk)
      logical plg1,plg2,plg3
c
c--------- STM is used after correction for TBIAS,
c--------- f temperature before correction for TBIAS
c
c
c---------------
      print *,' MIXSSTZ: im,jm,nl=',im,jm,nl
      print *,' TBIAS=',TBIAS
      print *,'kk,jj=',kk,jj
      print *,'  Zlev, k=1,nl'
c---------------
      write(6,202) (Zlev(k),k=1,nl)
 202  format(10f7.0)
c--------------- for debug
       i0=124
       j0=92
c--------------- for debug
       i1=124
       j1=88
       j2=98
c---------------
       print *,' BEFORE SST ASSIM: i1,j1,j2=',i1,j1,j2
       print *,'     f(i1,j,k), j=j1,j2   nl=',nl
       write(6,301) (j,j=j1,j2)
       do k=1,kk
        write(6,302) Zlev(k),(f(i1,j,k),j=j1,j2)
       end do
 301   format(5x,11i7)
 302   format(f8.0,11f7.2)
c---------------
c-----  remove correction for TBIAS
      do j=1,jm
      do i=1,im
         if(FSM(i,j).eq.1.) STM(i,j)=STM(i,j)+TBIAS
      end do
      end do
c----------------
      write(6,401) i0,j0,STM(i0,j0),H(i0,j0)
 401  format('i0,j0,STM(i0,j0),H(i0,j0)=',
     *            /2i7,f7.2,f7.0)
c---------------
      print *,'     STM(i1,j), j=j1,j2 '
      write(6,303) (STM(i1,j),j=j1,j2)
 303  format(8x,11f7.2)
c
c------------------ check difference between STM and f(i,j,1)
c
      difmax=0.
      difmin=0.
      do i=1,im
      do j=1,jm
       if(FSM(i,j).gt.0.5) then
        difp=STM(i,j)-f(i,j,1)
        difm=-difp
        if(difp.gt.difmax) then
         imax=i
         jmax=j
         difmax=difp
        end if
        if(difm.gt.difmin) then
         imin=i
         jmin=j
         difmin=difm
        end if
       end if
      end do
      end do
c------------------
      print *,' MIXSSTZ: max(STM-f(1)),i,j'
      print *,' difmax,idifmax,jdifmax=',difmax,imax,jmax
      print *,' MIXSSTZ: max(f(1)-STM),i,j'
      print *,' difmin,imin,jmin=',difmin,imin,jmin
c------------------ specify z and zz
      do j=1,jj
       z(j)=step*float(j-1)
      end do
      do k=1,kk
       zz(k)=Zlev(k)
      end do
c------------------
c-----------------------------   loops in i,j
      do j=1,jm
      do i=1,im
c--------------- for debug
c      plg1=i.eq.i0.and.j.eq.j0
       plg1=i.eq.i1.and.j.ge.j1.and.j.le.j2
       plg2=i.eq.imax.and.j.eq.jmax
       plg3=i.eq.imin.and.j.eq.jmin
c---------------
       if(FSM(i,j).eq.0.) go to 1000
c------------         send tt
       do k=1,kk
        tt(k)=f(i,j,k)
       end do
c------------         send SST
       SST=STM(i,j)
c------------         SP interpolate  to even step z
       call SP(zz,tt,t,z,kk,jj,kk,jj)
c------------
c      if(i.eq.i0.and.j.eq.j0) then
       if(plg1.or.plg2.or.plg3) then
c------------
        print *,'i,j,SST=',i,j,SST
        print *,'initial depth zz'
        write(6,101) zz
 101    format(10f7.0)
c------------
        print *,'initial temperature tt'
        write(6,102) tt
 102    format(10f7.2)
c------------
        print *,' depth for interpolation'
        write(6,101) z
c------------
        print *,'temperature after SPLINE interpolation'
        write(6,102) t
c------------
       end if
c
c------------- find n for the depth of mixed layer
       n=1
       do while ((t(1)-t(n+1)).le.0.5.and.n.lt.(jj))
        n=n+1
       end do
c------------- if mixed layer more than 80m send 80m (n=17)
       if(n.gt.17) n=17
c---------------------- correct SST if it deviates too much from t(1)
       if(SST-t(1).gt.5.)  SST=t(1)+5.
       if(SST-t(1).lt.-5.) SST=t(1)-5.
c---------------------- correct SST if SST<t(jj)
       if(SST.lt.t(jj)) SST=t(jj)
c----------------------
       DTS=SST-t(1)
c----------------------
       do k=1,n
          tc(k)=t(k)+DTS
       end do
c----------------------
       m=jj-n
c----------------------
       do k=n+1,jj
        tc(k)=t(k)+DTS*float(jj-k)/float(m)
       end do
c----------------------
c      if(i.eq.i0.and.j.eq.j0) then
c      if(plg2.or.plg3) then
       if(plg1.or.plg2.or.plg3) then
         print *,'i,j,SST=',i,j,SST
         print *,'MIX: n,m=',n,m
         write(6,201) SST,t(1),DTS
 201     format('SST,t(1),DTS=',3f7.2)
c------------
         print *,' depth for interpolation'
         write(6,101) z
c------------
         print *,'temperature in MIX before checking'
         write(6,102) tc
       end if
c------------- check if(tc(k).lt.t(jj))
       do k=2,jj-1
        if(tc(k).le.t(jj)) then
          n=k-1
          m=jj-n
          DTS=t(jj)-tc(n)
c         write(6,104) n,m,DTS
 104      format('n,m,DTS=',2i7,f7.2)
          go to 100
        end if
       end do
       go to 200
 100   continue
       do k=n+1,jj
         tc(k)=t(jj)-DTS*float(jj-k)/float(m)
       end do
c------------
 200   continue
c------------- check stability
       do k=2,jj
         if(tc(k-1).lt.tc(k)) tc(k)=tc(k-1)
       end do
c------------ interpolate back to GDEM Z-levels
       call SP(z,tc,te,zz,jj,kk,jj,kk)
c------------
       do k=1,kk
         dt(k)=te(k)-tt(k)
       end do
c------------
c      if(i.eq.i0.and.j.eq.j0) then
c      if(plg2.or.plg3) then
       if(plg1.or.plg2.or.plg3) then
         print *,'i,j,SST=',i,j,SST
         print *,'temperature in MIX after checking'
         write(6,102) tc
c------------
         print *,'initial depth zz'
         write(6,101) zz
c------------
         print *,'initial temperature tt'
         write(6,102) tt
c------------
         print *,'temperature after back interpolation'
         write(6,102) te
c------------
         print *,'temperature differences after MIX'
         write(6,102) dt
       end if
c---------------------- send back adjusted temperature
       do k=1,kk
         f(i,j,k)=te(k)
       end do
c------------
 1000  continue
c--------------------- end loops in i,j
      end do
      end do
c----------------------
       print *,' AFTER SST ASSIM: i1,j1,j2=',i1,j1,j2
       print *,'     f(i1,j,k), j=j1,j2   nl=',nl
       write(6,301) (j,j=j1,j2)
       do k=1,kk
        write(6,302) Zlev(k),(f(i1,j,k),j=j1,j2)
       end do
c----------------------
      return
      end
