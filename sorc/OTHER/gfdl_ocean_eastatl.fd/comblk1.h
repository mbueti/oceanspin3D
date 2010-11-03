      PARAMETER (IMA=66,JMA=66)
      COMMON /TOCN/  TOSTEP,TOCEAN,TOHOUR,TOCAT,ISM,PRT,PRD,IOSTEP 
      REAL MINUTS1,MINUTS2
      REAL DTE2,DTI2,DAYI, julday
      REAL ISPI,ISP2I,TMP,RHO_0
      integer iend,year 
      integer*4 startdate,enddate
      CHARACTER*40 FRSTI, FRSTO
      CHARACTER*15 pathname,pathname1
c
      COMMON/connect1/UTB(IM,JM),VTB(IM,JM),UTF(IM,JM),VTF(IM,JM),
     1     ADVX(IM,JM,KB),ADVY(IM,JM,KB),ADVUA(IM,JM),ADVVA(IM,JM),
     2     TSURF(IM,JM),SSURF(IM,JM),WINDX2(IM,JM),WINDY2(IM,JM),
     2     DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM),
     3     ADVUU(IM,JM),ADVVV(IM,JM),WINDX1(IM,JM),WINDY1(IM,JM),
     4     SWRAD(IM,JM),SSTIN(IM,JM),IDAMP,ISMOTH,IHOURS,
     5     MINUTS1,MINUTS2,DTE2,DTI2,IEND,NREAD,IPRTH1,INOWINDH,
     6     DAYI,ISPI,ISP2I,FRSTI,FRSTO,SMH,MODE,TIME0,ISPLIT,IMAY
c	    
      COMMON/ATMOCN1/XSU(IM),XSV(IM),XST(IM),YSU(JM),YSV(JM),YST(JM), 
     1                               TMA(IM,JM)   
c
      COMMON/ATMOCN2/SST(IMA,JMA),XSTR(IMA,JMA),YSTR(IMA,JMA),                   
     1           HFLX(IMA,JMA),XLON(IMA),YLAT(JMA),
     2           RFSW(IMA,JMA)               
      real*4 SST,XSTR,YSTR,HFLX,XLON,YLAT,RFSW
C
      REAL LAMBDA,PHI,REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX,CMP
c
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
c
      common/misc/NBC,nbc2d(im,jm),nbc2ds(im,jm),RHO_0,pi,small,
     1            beta,ispadv,smoth,horcon,iswtch,iprtd2
c
      common/dating/julday,startdate,enddate,year,pathname,pathname1
