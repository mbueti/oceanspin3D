        
        subroutine frlvstress(time,frlvhstart,wusurf,wvsurf,taux,tauy)
c
c Interpolate Sergei Frolov's HRD wind stress for ocean model
c Richard M. Yablonsky on 12/13/2006
c

        parameter(im=254,jm=225,nf=8)
        parameter(rho_0=1024.e0)
        character*2 chtlo,chthi
        character*9 chf(nf)
        real wul(im,jm),wuh(im,jm),wvl(im,jm),wvh(im,jm)
        real txl(im,jm),txh(im,jm),tyl(im,jm),tyh(im,jm)
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

        open(11,file=chf(1))
           read(11,901) wul
        close(11)
 
        open(11,file=chf(2))
           read(11,901) wuh
        close(11)

        open(11,file=chf(3))
           read(11,901) wvl
        close(11)
                                                                                                                                                             
        open(11,file=chf(4))
           read(11,901) wvh
        close(11)

        open(11,file=chf(5))
           read(11,901) txl
        close(11)
                                                                                                                                                             
        open(11,file=chf(6))
           read(11,901) txh
        close(11)
                                                                                                                                                             
        open(11,file=chf(7))
           read(11,901) tyl
        close(11)
                                                                                                                                                             
        open(11,file=chf(8))
           read(11,901) tyh
        close(11)

        do i=1,im
        do j=1,jm
           wusurf(i,j)=wusurf(i,j)-(faclo*wul(i,j)+fachi*wuh(i,j))/rho_0 
           wvsurf(i,j)=wvsurf(i,j)-(faclo*wvl(i,j)+fachi*wvh(i,j))/rho_0
           taux(i,j)=taux(i,j)+faclo*txl(i,j)+fachi*txh(i,j)
           tauy(i,j)=tauy(i,j)+faclo*tyl(i,j)+fachi*tyh(i,j)
        enddo
        enddo

        return
        end
