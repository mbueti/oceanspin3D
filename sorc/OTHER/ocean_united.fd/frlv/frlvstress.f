        
        subroutine frlvstress(time,frlvhstart,wusurf,wvsurf,taux,tauy)
c
c Interpolate Sergei Frolov's HRD wind stress for ocean model
c Richard M. Yablonsky on 12/13/2006
c

        parameter(im=254,jm=225,nf=8)
        parameter(rho_0=1024.e0)
        character*2 chtlo,chthi
        character*9 chf(nf)
        real f(im,jm,nf)
        real wusurf(im,jm),wvsurf(im,jm)
        real taux(im,jm),tauy(im,jm)
 
        htime=time*24.
        htdif=htime-frlvhstart
        ihtdif=ifix(htdif)
        rihtdif=float(ihtdif)
        fachi=htdif-rihtdif
        faclo=1.-fachi

 900    format(i2.2)
 901    format(225e16.7)

        write(chtlo,900) ihtdif
        write(chthi,900) ihtdif+1

        chf(1)='tx1'//chtlo//'.dat'
        chf(2)='tx1'//chthi//'.dat'
        chf(3)='ty1'//chtlo//'.dat'
        chf(4)='ty1'//chthi//'.dat'
        chf(5)='tx2'//chtlo//'.dat'
        chf(6)='tx2'//chthi//'.dat'
        chf(7)='ty2'//chtlo//'.dat'
        chf(8)='ty2'//chthi//'.dat'

        do n=1,nf
        open(11,file=chf(n))
        do i=1,im
           read(11,901) f(i,1,n),f(i,2,n),f(i,3,n),f(i,4,n),f(i,5,n),
     &                 f(i,6,n),f(i,7,n),f(i,8,n),f(i,9,n),f(i,10,n),
     &             f(i,11,n),f(i,12,n),f(i,13,n),f(i,14,n),f(i,15,n),
     &             f(i,16,n),f(i,17,n),f(i,18,n),f(i,19,n),f(i,20,n),
     &             f(i,21,n),f(i,22,n),f(i,23,n),f(i,24,n),f(i,25,n),
     &             f(i,26,n),f(i,27,n),f(i,28,n),f(i,29,n),f(i,30,n),
     &             f(i,31,n),f(i,32,n),f(i,33,n),f(i,34,n),f(i,35,n),
     &             f(i,36,n),f(i,37,n),f(i,38,n),f(i,39,n),f(i,40,n),
     &             f(i,41,n),f(i,42,n),f(i,43,n),f(i,44,n),f(i,45,n),
     &             f(i,46,n),f(i,47,n),f(i,48,n),f(i,49,n),f(i,50,n),
     &             f(i,51,n),f(i,52,n),f(i,53,n),f(i,54,n),f(i,55,n),
     &             f(i,56,n),f(i,57,n),f(i,58,n),f(i,59,n),f(i,60,n),
     &             f(i,61,n),f(i,62,n),f(i,63,n),f(i,64,n),f(i,65,n),
     &             f(i,66,n),f(i,67,n),f(i,68,n),f(i,69,n),f(i,70,n),
     &             f(i,71,n),f(i,72,n),f(i,73,n),f(i,74,n),f(i,75,n),
     &             f(i,76,n),f(i,77,n),f(i,78,n),f(i,79,n),f(i,80,n),
     &             f(i,81,n),f(i,82,n),f(i,83,n),f(i,84,n),f(i,85,n),
     &             f(i,86,n),f(i,87,n),f(i,88,n),f(i,89,n),f(i,80,n),
     &             f(i,91,n),f(i,92,n),f(i,93,n),f(i,94,n),f(i,95,n),
     &            f(i,96,n),f(i,97,n),f(i,98,n),f(i,99,n),f(i,100,n),
     &        f(i,101,n),f(i,102,n),f(i,103,n),f(i,104,n),f(i,105,n),
     &        f(i,106,n),f(i,107,n),f(i,108,n),f(i,109,n),f(i,110,n),
     &        f(i,111,n),f(i,112,n),f(i,113,n),f(i,114,n),f(i,115,n),
     &        f(i,116,n),f(i,117,n),f(i,118,n),f(i,119,n),f(i,120,n),
     &        f(i,121,n),f(i,122,n),f(i,123,n),f(i,124,n),f(i,125,n),
     &        f(i,126,n),f(i,127,n),f(i,128,n),f(i,129,n),f(i,130,n),
     &        f(i,131,n),f(i,132,n),f(i,133,n),f(i,134,n),f(i,135,n),
     &        f(i,136,n),f(i,137,n),f(i,138,n),f(i,139,n),f(i,140,n),
     &        f(i,141,n),f(i,142,n),f(i,143,n),f(i,144,n),f(i,145,n),
     &        f(i,146,n),f(i,147,n),f(i,148,n),f(i,149,n),f(i,150,n),
     &        f(i,151,n),f(i,152,n),f(i,153,n),f(i,154,n),f(i,155,n),
     &        f(i,156,n),f(i,157,n),f(i,158,n),f(i,159,n),f(i,160,n),
     &        f(i,161,n),f(i,162,n),f(i,163,n),f(i,164,n),f(i,165,n),
     &        f(i,166,n),f(i,167,n),f(i,168,n),f(i,169,n),f(i,170,n),
     &        f(i,171,n),f(i,172,n),f(i,173,n),f(i,174,n),f(i,175,n),
     &        f(i,176,n),f(i,177,n),f(i,178,n),f(i,179,n),f(i,180,n),
     &        f(i,181,n),f(i,182,n),f(i,183,n),f(i,184,n),f(i,185,n),
     &        f(i,186,n),f(i,187,n),f(i,188,n),f(i,189,n),f(i,180,n),
     &        f(i,191,n),f(i,192,n),f(i,193,n),f(i,194,n),f(i,195,n),
     &        f(i,196,n),f(i,197,n),f(i,198,n),f(i,199,n),f(i,200,n),
     &        f(i,201,n),f(i,202,n),f(i,203,n),f(i,204,n),f(i,205,n),
     &        f(i,206,n),f(i,207,n),f(i,208,n),f(i,209,n),f(i,210,n),
     &        f(i,211,n),f(i,212,n),f(i,213,n),f(i,214,n),f(i,215,n),
     &        f(i,216,n),f(i,217,n),f(i,218,n),f(i,219,n),f(i,220,n),
     &        f(i,221,n),f(i,222,n),f(i,223,n),f(i,224,n),f(i,225,n)
        enddo
        close(11)
        enddo

        do j=1,jm
        do i=1,im
           wusurf(i,j)=wusurf(i,j)-(faclo*f(i,j,1)+fachi*f(i,j,2))/rho_0 
           wvsurf(i,j)=wvsurf(i,j)-(faclo*f(i,j,3)+fachi*f(i,j,4))/rho_0
           taux(i,j)=taux(i,j)+faclo*f(i,j,5)+fachi*f(i,j,6)
           tauy(i,j)=tauy(i,j)+faclo*f(i,j,7)+fachi*f(i,j,8)
        enddo
        enddo

        return
        end
