      REAL KM,KH,KQ,L
      PARAMETER (MMX=6,NMX=10)
c
Cfr THIS PARAMETER SPECIFIES THE LENGTH OF DATA VECTOR FOR HURRICANE PATH
      PARAMETER (IM=254,JM=225,KB=23,IMM1=IM-1,JMM1=JM-1,KBM1=KB-1)
      PARAMETER (IMM2=IM-2,JMM2=JM-2,KBM2=KB-2)
      PARAMETER (LIJ=IM*JM,LIJ2=LIJ*2,LIJK=LIJ*KB,LIJKM2=LIJ*KBM2)
      PARAMETER (LIJM1=IM*(JM-1),LIJKM1=LIJ*KBM1)
      PARAMETER (IMAO=293,JMAO=227)
      COMMON/BLKCON/
     1          H_MAX,IINT,IPRINT,DTE,DTI,TPRNU,UMOL,
     2          GRAV,TIME,RAMP,
     3          TBIAS,SBIAS
C---------------- 1-D ARRAYS --------------------------------------
      COMMON/BLK1D/
     1      DZR(KB),Z(KB),ZZ(KB),DZ(KB),DZZ(KB)
C---------------- 2-D ARRAYS --------------------------------------
      COMMON/BLK2D/H(IM,JM),DX(IM,JM),DY(IM,JM),D(IM,JM),DT(IM,JM),
     1     ART(IM,JM),ARU(IM,JM),ARV(IM,JM),CBC(IM,JM),
     2     ALON(IM,JM),ALAT(IM,JM),ANG(IM,JM),
     3     DUM(IM,JM),DVM(IM,JM),FSM(IM,JM),COR(IM,JM),
     4     WUSURF(IM,JM),WVSURF(IM,JM),WUBOT(IM,JM),WVBOT(IM,JM),
     5     WTSURF(IM,JM),WSSURF(IM,JM),TPS(IM,JM),AAM2D(IM,JM),
     6     UAF(IM,JM),UA(IM,JM),UAB(IM,JM),VAF(IM,JM),VA(IM,JM),
     7     VAB(IM,JM),ELF(IM,JM),EL(IM,JM),ELB(IM,JM),PSI(IM,JM),
     8     ETF(IM,JM),ET(IM,JM),ETB(IM,JM),FLUXUA(IM,JM),FLUXVA(IM,JM),
     9     EGF(IM,JM),EGB(IM,JM),STM(IM,JM)
     *    ,TAUX(IM,JM),TAUY(IM,JM)
     *    ,windx(IM,JM),windy(IM,JM)
C---------------- 3-D ARRAYS --------------------------------------
      COMMON/BLK3D/
     1     A(IM,JM,KB),C(IM,JM,KB),EE(IM,JM,KB),GG(IM,JM,KB),
     1     UF(IM,JM,KB),VF(IM,JM,KB),
     2     KM(IM,JM,KB),KH(IM,JM,KB),KQ(IM,JM,KB),L(IM,JM,KB),
     3     Q2(IM,JM,KB),Q2B(IM,JM,KB),AAM(IM,JM,KB),
     4     Q2L(IM,JM,KB),Q2LB(IM,JM,KB),
     5     U(IM,JM,KB),UB(IM,JM,KB),W(IM,JM,KB),
     6     V(IM,JM,KB),VB(IM,JM,KB),TMEAN(IM,JM,KB),SMEAN(IM,JM,KB),
     7     T(IM,JM,KB),TB(IM,JM,KB),TBIN(IM,JM,KB),SBIN(IM,JM,KB),
     8     S(IM,JM,KB),SB(IM,JM,KB),RMEAN1(IM+NMX-1,JM+MMX-1,KB),
     9     RHO(IM,JM,KB),DTEF(IM,JM,KB),RMEAN(IM,JM,KB)
C----------- 1 AND 2-D BOUNDARY VALUE ARRAYS ------------------------
      COMMON/BDRY/
     1     TBW(JM,KB),TBE(JM,KB),TBN(IM,KB),TBS(IM,KB),SBN(IM,KB),
     2     SBW(JM,KB),SBE(JM,KB),SBS(IM,KB),UABE(JM),UABW(JM),
     3     VABN(IM),VABS(IM),ELBN(IM),ELBE(JM),ELBW(JM),ELBS(IM),
     4     COVRHN(IM),COVRHS(IM),COVRHE(JM),COVRHW(JM)
c----------- Arrays for nudging ------------------------------------
      COMMON/NDGE/ spwghts(im,jm),tnudge(im,jm,kb),snudge(im,jm,kb)
c----------- falk 01-05-04 add parameters for avrtau and FSMA
      common/tau/ tauavr,taumax,awucon,bwucon,migtau
      common/maskat/ FSMA(imao,jmao)
