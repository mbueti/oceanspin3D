	SUBROUTINE MIXSSTZ(f,STM,H,FSM,Zlev,im,jm,nl,TBIAS)
c---------------------------------------------------------------------
c----- this subr. was written 06/18/01 by A.Falkovich
c----- to assimilate SST real data in GDEM monthly t data for Z-levels
c---------------------------------------------------------------------
c
c---- MXD is the max depth where climatological temperature is modified ---
c----      for SST > f(1)
c---- it is computed as climat. mixed layer depth + transition depth ------
c---- transition depth is calculated on the basis of simple criteria ------
c---- the added temperature gradient should not exceed 2 deg per 50m ------
c
c----  For cases with SST<TB temperature is modified up to level 
c----   of intersection SST with initial temp. profile
c
      real STM(im,jm),H(im,jm),FSM(im,jm)
      real f(IM,JM,nl),Zlev(nl),dt(100),dt1(100)
      real MXD, DTS
c
c--------- STM is used after correction for TBIAS,
c--------- f temperature before correction for TBIAS
c
c
c     print *,'  MIXSSTZ: check STM'
c     print *,' TBIAS=',TBIAS
c     i=113
c     j=185
c     write(6,101) i,j,STM(i,j),H(i,j)
c
c----------------------
       i1=69
       j1=162
       j2=171
c      print *,' BEFORE SST ASSIM: i1,j1,j2=',i1,j1,j2
c      print *,'     f(i1,j,k), j=j1,j2   nl=',nl
c      write(6,301) (j,j=j1,j2)
c      do k=1,nl
c       write(6,302) Zlev(k),(f(i1,j,k),j=j1,j2)
c      end do
c----------------------
c-----  remove correction for TBIAS
      do j=1,jm
      do i=1,im
         if(FSM(i,j).eq.1.) STM(i,j)=STM(i,j)+TBIAS
      end do
      end do
c     print *,'     STM(i1,j), j=j1,j2 '
c     write(6,102) (STM(i1,j),j=j1,j2)
c        
c------------------  04-17-02 check STM as in MIXSST
c
      difmax=0.
      difmin=0.
      do i=1,im
      do j=1,jm
       if(FSM(i,j).gt.0.5) then
        difp=STM(i,j)-f(i,j,1)
        difm=-difp
        if(difp.gt.difmax) then
         idifmax=i
         jdifmax=j
         difmax=difp
        end if
        if(difm.gt.difmin) then
         idifmin=i
         jdifmin=j
         difmin=difm
        end if
       end if
      end do
      end do
      print *,' MIXSSTZ: max(STM-f(1)),i,j'
      print *,' difmax,idifmax,jdifmax=',difmax,idifmax,jdifmax
      print *,' MIXSSTZ: max(f(1)-STM),i,j'
      print *,' difmin,idifmin,jdifmin=',difmin,idifmin,jdifmin
c------------------
      CALL OADJUSTZ(f,im,jm,nl,0,Zlev,FSM)
c----------------------
      nlm1=nl-1
c----------------------
c----------------------          loop in i,j
      do j=1,jm
        do i=1,im
         if(FSM(i,j).eq.0.) go to 1000
         do k=1,nlm1
          dt(k)=f(i,j,k+1)-f(i,j,k)
         end do
         dt(nl)=0.
         do k=1,nl
          dt1(k)=f(i,j,k)
         end do
c------------- find n for the depth of mixed layer
         n=1
         do while ((f(i,j,1)-f(i,j,n+1)).le.0.5.and.n.lt.(nl-1))
           n=n+1
         end do
c----------------------
c        if(i.eq.113.and.(j.ge.185.and.j.le.205)) then
c         print *,' MIXSSTZ before assimilation'
c         write(6,101) i,j,STM(i,j),H(i,j)
 101      format('i,j,STM(i,j),H(i,j)=',
     *            /2i7,f7.2,f7.0)
c         print *,'     f(i,j,k), k=1,nl   nl=',nl
c         write(6,102) (f(i,j,k),k=1,nl)
c         print *,'  Zlev, k=1,nl'
c         write(6,202) (Zlev(k),k=1,nl)
 102      format(10f7.2)
 202      format(10f7.0)
c        end if
c----------------------
c
         if(STM(i,j).gt.f(i,j,1)) then
c----------------------
c          if(i.eq.56.and.j.eq.1) then
c           print *,' MIXSSTZ before assimilation'
c           write(6,101) i,j,STM(i,j),H(i,j)
c           print *,'     f(i,j,k), k=1,nl   nl=',nl
c           write(6,102) (f(i,j,k),k=1,nl)
c           print *,'  Zlev, k=1,nl'
c           write(6,202) (Zlev(k),k=1,nl)
c          end if
c----------------------
            DTS=STM(i,j)-f(i,j,1)
c----------------------
            n1=1
            MXD=Zlev(n)+DTS/2.*50.
            do while (Zlev(n1).le.MXD.and.n1.lt.nl)
              n1=n1+1
            end do
c----------------------
            tmixbold=f(i,j,n)
c
            do k=1,n
              f(i,j,k)=f(i,j,k)+DTS
            end do
c
            tmixb=f(i,j,n)
            tbotMXD=f(i,j,n1)
            dtbold=max(0.1,tmixbold-tbotMXD)
            constgrd=(tmixb-tbotMXD)/dtbold
c
            do k=n,n1-1
              f(i,j,k+1)=f(i,j,k)+constgrd*dt(k)
            end do
            do k=1,nl
              dt1(k)=f(i,j,k)-dt1(k)
            end do
c----------------------
c           if(i.eq.56.and.j.eq.1) then
c            write(6,106) STM(i,j),DTS,
c    *              tmixb,constgrd,MXD
 106         format('STM(i,j),DTS,',
     *        'tmixb,constgrd,MXD=',/5f7.2)
c            write(6,109) tmixbold, tbotMXD,constgrd
 109         format('tmixbold,tbotMXD,constgrd=',3f10.2)
c            print *,'   level of mixed layer depth n=',n
c            print *,' AFTER assimilation'
c            print *,' the last level of changes n1=',n1
c            print *,'     f(i,j,k), k=1,nl   nl=',nl
c            write(6,102) (f(i,j,k),k=1,nl)
c            print *,'  Zlev, k=1,nl'
c            write(6,202) (Zlev(k),k=1,nl)
c            print *,' dt1=f(k)-fold(k)'
c            write(6,102) (dt1(k),k=1,nl)
c           end if
c----------------------
         else
c
            DTS=STM(i,j)-f(i,j,1)
c----------------- t has to be more than -1.
            tmixb=max(-1.,f(i,j,1)+DTS)
            f(i,j,1)=tmixb
c----------------   find level where tmixb > t(k)
            n1=2
            do while (tmixb.lt.f(i,j,n1).and.n1.lt.nl)
             n1=n1+1
            end do
c------------ put t=const up to level nbot
            if(tmixb.le.f(i,j,nl)) then
             nbot=nl
            else 
             nbot=n1-1
            end if
            do k=2,nbot
             f(i,j,k)=tmixb
            end do
c             
c----------------------
c           if(i.eq.113.and.(j.ge.185.and.j.le.205)) then
c            print *,'    i,j=',i,j
c            print *,'   level of mixed layer depth n=',n
c            print *,'   level of intersection n1=',n1
c            print *,'   level of const temp nbot=',nbot
c            print *,' the last level of changes n2=',n2
c            write(6,103) STM(i,j),DTS,tmixb
 103         format('STM(i,j),DTS,tmixb',/3f10.2)
c           end if
c----------------------
            do k=1,nl
              dt1(k)=f(i,j,k)-dt1(k)
            end do
c----------------------
         end if
c----------- remove instability
          do k=2,nl
            if(f(i,j,k).gt.f(i,j,k-1)) f(i,j,k)=f(i,j,k-1)
          end do
c----------------------
c         if(i.eq.113.and.(j.ge.185.and.j.le.205)) then
c          print *,' AFTER assimilation'
c          print *,'     f(i,j,k), k=1,nl   nl=',nl
c          write(6,102) (f(i,j,k),k=1,nl)
c          print *,' dt1=f(k)-fold(k)'
c          write(6,102) (dt1(k),k=1,nl)
c         end if
c----------------------
 1000    continue
        end do
      end do
c
c----------------------
       i1=69
       j1=162
       j2=171
c      print *,' AFTER SST ASSIM: lll,i1,j1,j2=',lll,i1,j1,j2
c      print *,'     f(i1,j,k), j=j1,j2   nl=',nl
c      write(6,301) (j,j=j1,j2)
c      do k=1,nl
c       write(6,302) Zlev(k),(f(i1,j,k),j=j1,j2)
c      end do
c----------------------
         
      return
      end
