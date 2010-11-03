      SUBROUTINE OADJUST(iflag)

      include 'comblk.h'
      real tmn
      integer counter

      small=1.e-3
      counter=1
      ll=0
      DO WHILE(counter.gt.0)
      counter=0
      ll=ll+1
      do j=1,jm
         do i=1,im
            do n=1,2
               do k=n,kb-1,2
               if(t(i,j,k).lt.t(i,j,k+1).and.fsm(i,j).eq.1.) then
                   tmn=(t(i,j,k)*dz(k)+t(i,j,k+1)*dz(k+1))/
     1             (dz(k)+dz(k+1))
               t(i,j,k)=tmn
               t(i,j,k+1)=tmn-small
      if(iflag.eq.1) print*,'Convective adjustment T at (i,j,k)=',i,j,k
               counter=counter+1
               end if
               end do
            end do
         end do
      end do
      if(ll.eq.1) ii=counter
      END DO
      print*,'Adjusted in ',ii,' points after',ll,' iterations'
C
      counter=1
      DO WHILE(counter.gt.0)
      counter=0
      do j=1,jm
         do i=1,im
            do n=1,2
               do k=n,kb-1,2
               if(tb(i,j,k).lt.tb(i,j,k+1).and.fsm(i,j).eq.1.) then
                   tmn=(tb(i,j,k)*dz(k)+tb(i,j,k+1)*dz(k+1))/
     1             (dz(k)+dz(k+1))
               tb(i,j,k)=tmn
               tb(i,j,k+1)=tmn-small
c               print*,'Convective adjustment TB at (i,j,k)=',i,j,k
               counter=counter+1
               end if
               end do
            end do
         end do
      end do
      END DO
c
      return
      end
      SUBROUTINE ADVAVE(ADVUA,ADVVA,MODE)
       INCLUDE 'comblk.h'
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),CURV2D(IM,JM)
      EQUIVALENCE (TPS,CURV2D)
C---------------------------------------------------------------------
C      CALCULATE U-ADVECTION & DIFFUSION
C---------------------------------------------------------------------
C
C-------- ADVECTIVE FLUXES -------------------------------------------
C
      DO 200 J=1,JM
      DO 200 I=1,IM
 200  ADVUA(I,J)=0.
      DO 300 J=2,JM
      DO 300 I=2,IMM1
 300  FLUXUA(I,J)=.125E0*((D(I+1,J)+D(I,J))*UA(I+1,J)
     1                 +(D(I,J)+D(I-1,J))*UA(I,J))
     2                  *(UA(I+1,J)+UA(I,J))
      DO 400 J=2,JM
      DO 400 I=2,IM
 400  FLUXVA(I,J)=.125E0*((D(I,J)+D(I,J-1))*VA(I,J)
     1                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     2                    *(UA(I,J)+UA(I,J-1))
C----------- ADD VISCOUS FLUXES ---------------------------------
      DO  460 J=2,JM
      DO  460 I=2,IMM1
 460  FLUXUA(I,J)=FLUXUA(I,J)
     1         -D(I,J)*2.E0*AAM2D(I,J)*(UAB(I+1,J)-UAB(I,J))/DX(I,J)
      DO  470 J=2,JM
      DO  470 I=2,IM
      TPS(I,J)=.25E0*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     1            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     2                *((UAB(I,J)-UAB(I,J-1))
     3                /(DY(I,J)+DY(I-1,J)+DY(I,J-1)+DY(I-1,J-1))
     4                 +(VAB(I,J)-VAB(I-1,J))
     5                /(DX(I,J)+DX(I-1,J)+DX(I,J-1)+DX(I-1,J-1)) )
      FLUXUA(I,J)=FLUXUA(I,J)*DY(I,J)
      FLUXVA(I,J)=(FLUXVA(I,J)-TPS(I,J))
     1            *.25E0*(DX(I,J)+DX(I-1,J)+DX(I,J-1)+DX(I-1,J-1))
 470  CONTINUE
C----------------------------------------------------------------
      DO  480 J=2,JMM1
      DO  480 I=2,IMM1
 480  ADVUA(I,J)=FLUXUA(I,J)-FLUXUA(I-1,J)
     1           +FLUXVA(I,J+1)-FLUXVA(I,J)
C----------------------------------------------------------------
C       CALCULATE V-ADVECTION & DIFFUSION
C----------------------------------------------------------------
C
      DO 481 J=1,JM
      DO 481 I=1,IM
 481  ADVVA(I,J)=0.
C
C---------ADVECTIVE FLUXES ----------------------------
      DO 700 J=2,JM
      DO 700 I=2,IM
 700  FLUXUA(I,J)=.125E0*((D(I,J)+D(I-1,J))*UA(I,J)
     1         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     2                        (VA(I-1,J)+VA(I,J))
      DO 800 J=2,JMM1
      DO 800 I=2,IM
 800  FLUXVA(I,J)=.125E0*((D(I,J+1)+D(I,J))
     1       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     2      *(VA(I,J+1)+VA(I,J))
C------- ADD VISCOUS FLUXES -----------------------------------
      DO  860 J=2,JMM1
      DO  860 I=2,IM
 860  FLUXVA(I,J)=FLUXVA(I,J)
     1        -D(I,J)*2.E0*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/DY(I,J)
      DO  870 J=2,JM
      DO  870 I=2,IM
      FLUXVA(I,J)=FLUXVA(I,J)*DX(I,J)
 870  FLUXUA(I,J)=(FLUXUA(I,J)-TPS(I,J))
     3             *.25E0*(DY(I,J)+DY(I-1,J)+DY(I,J-1)+DY(I-1,J-1))
C---------------------------------------------------------------
      DO  880 J=2,JMM1
      DO  880 I=2,IMM1
 880  ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     1          +FLUXVA(I,J)-FLUXVA(I,J-1)
C
C---------------------------------------------------------------
      IF(MODE.NE.2) GO TO 5000
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      WUBOT(I,J)=-0.5E0*(CBC(I,J)+CBC(I-1,J))
     1     *SQRT(UAB(I,J)**2+(.25E0*(VAB(I,J)
     2     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
  100 CONTINUE
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
      WVBOT(I,J)=-0.5E0*(CBC(I,J)+CBC(I,J-1))
     1     *SQRT((.25E0*(UAB(I,J)+UAB(I+1,J)
     2     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
  102 CONTINUE
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
      CURV2D(I,J)=.25*((VA(I,J+1)+VA(I,J))*(DY(I+1,J)-DY(I-1,J))
     1               -(UA(I+1,J)+UA(I,J))*(DX(I,J+1)-DX(I,J-1)) )
     2                /(DX(I,J)*DY(I,J))
      ADVUA(I,J)=ADVUA(I,J)   
     1    -ARU(I,J)*.25*(  CURV2D(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     2          +CURV2D(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
      ADVVA(I,J)=ADVVA(I,J)
     1    +ARV(I,J)*.25*(  CURV2D(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     2          +CURV2D(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
  120 CONTINUE  
 5000 CONTINUE
C 
      RETURN
      END
      SUBROUTINE ADVCT(ADVX,ADVY)        
C======================================================================
C  This routine calculates the horizontal portions of momentum advection
C  well in advance of their use in ADVU and ADVV so that their vertical 
C  integrals (created in MAIN) may be used in the external mode calculation.
C======================================================================
       INCLUDE 'comblk.h'   
      DIMENSION ADVX(IM,JM,KB),ADVY(IM,JM,KB),
     1          XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),
     1          CURV(IM,JM,KB)
      EQUIVALENCE (A,XFLUX),(C,YFLUX),(EE,CURV)
C
C
      DO 60 K=1,KBM1
      DO 60 J=2,JMM1
      DO 60 I=2,IMM1
  60  CURV(I,J,K)=
     1      +.25E0*((V(I,J+1,K)+V(I,J,K))*(DY(I+1,J)-DY(I-1,J))
     2               -(U(I+1,J,K)+U(I,J,K))*(DX(I,J+1)-DX(I,J-1)) )
     3                /(DX(I,J)*DY(I,J))
C-----------------------------------------------------------------
C      Calculate x-component of velocity advection                 
C-----------------------------------------------------------------
      DO 59 K=1,KB
      DO 59 J=1,JM
      DO 59 I=1,IM
      ADVX(I,J,K)=0.E0
      XFLUX(I,J,K)=0.E0
  59  YFLUX(I,J,K)=0.E0
C
C******** HORIZONTAL ADVECTION FLUXES *****************************
      DO 100 K=1,KBM1
      DO 100 J=1,JM
      DO 100 I=2,IMM1
 100  XFLUX(I,J,K)=.125E0*((DT(I+1,J)+DT(I,J))*
     1              U(I+1,J,K)+(DT(I,J)+DT(I-1,J))
     2             *U(I,J,K))*(U(I+1,J,K)+U(I,J,K))
      DO 120 K=1,KBM1
      DO 120 J=2,JM
      DO 120 I=2,IM
 120  YFLUX(I,J,K)=.125E0*((DT(I,J)+DT(I,J-1))
     1        *V(I,J,K)+(DT(I-1,J)+DT(I-1,J-1))
     2        *V(I-1,J,K))*(U(I,J,K)+U(I,J-1,K))
C****** ADD HORIZONTAL DIFFUSION FLUXES ****************************
      DO  130 K=1,KBM1
      DO  130 J=2,JM
      DO  130 I=2,IMM1
      XFLUX(I,J,K)=XFLUX(I,J,K)
     1        -DT(I,J)*AAM(I,J,K)*2.E0*(UB(I+1,J,K)-UB(I,J,K))/DX(I,J)
      DTAAM=.25E0*(DT(I,J)+DT(I-1,J)+DT(I,J-1)+DT(I-1,J-1))
     1           *(AAM(I,J,K)+AAM(I-1,J,K)+AAM(I,J-1,K)+AAM(I-1,J-1,K))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     1          -DTAAM*((UB(I,J,K)-UB(I,J-1,K))
     2                    /(DY(I,J)+DY(I-1,J)+DY(I,J-1)+DY(I-1,J-1))
     3                    +(VB(I,J,K)-VB(I-1,J,K))
     4                    /(DX(I,J)+DX(I-1,J)+DX(I,J-1)+DX(I-1,J-1)))
C
      XFLUX(I,J,K)=DY(I,J)*XFLUX(I,J,K)
      YFLUX(I,J,K)=
     1     .25E0*(DX(I,J)+DX(I-1,J)+DX(I,J-1)+DX(I-1,J-1))*YFLUX(I,J,K)
  130 CONTINUE 
C
C
C******** DO HORIZ. ADVECTION *******
      DO 146 K=1,KBM1
      DO 146 J=2,JMM1
      DO 146 I=2,IMM1
 146  ADVX(I,J,K)= 
     1            +XFLUX(I,J,K)-XFLUX(I-1,J,K)
     2            +YFLUX(I,J+1,K)-YFLUX(I,J,K)
     3   -ARU(I,J)*.25*(CURV(I,J,K)*DT(I,J)*(V(I,J+1,K)+V(I,J,K))
     4             +CURV(I-1,J,K)*DT(I-1,J)*(V(I-1,J+1,K)+V(I-1,J,K)))
C
C-----------------------------------------------------------------
C      Calculate y-component of velocity advection                 
C-----------------------------------------------------------------
      DO 299 K=1,KB
      DO 299 J=1,JM
      DO 299 I=1,IM
      ADVY(I,J,K)=0.E0
      XFLUX(I,J,K)=0.E0
 299  YFLUX(I,J,K)=0.E0
C
C********** HORIZONTAL ADVECTION FLUXES **************************
      DO 300 K=1,KBM1
      DO 300 J=2,JM
      DO 300 I=2,IM
 300  XFLUX(I,J,K)=.125E0*((DT(I,J)+DT(I-1,J))*U(I,J,K)
     1                  +(DT(I,J-1)+DT(I-1,J-1))*U(I,J-1,K))
     2                        *(V(I,J,K)+V(I-1,J,K))
      DO 320 K=1,KBM1
      DO 320 J=2,JMM1
      DO 320 I=1,IM
 320  YFLUX(I,J,K)=.125E0*((DT(I,J+1)+DT(I,J))*V(I,J+1,K)
     1                  +(DT(I,J)+DT(I,J-1))*V(I,J,K))
     2                        *(V(I,J+1,K)+V(I,J,K))
C******* ADD HORIZONTAL DIFFUSION FLUXES **************************
      DO  700 K=1,KBM1
      DO  700 J=2,JMM1
      DO  700 I=2,IM
      DTAAM=.25E0*(DT(I,J)+DT(I-1,J)+DT(I,J-1)+DT(I-1,J-1))
     1           *(AAM(I,J,K)+AAM(I-1,J,K)+AAM(I,J-1,K)+AAM(I-1,J-1,K))
      XFLUX(I,J,K)=XFLUX(I,J,K)
     1     -DTAAM*((UB(I,J,K)-UB(I,J-1,K))
     2         /(DY(I,J)+DY(I-1,J)+DY(I,J-1)+DY(I-1,J-1))
     3               +(VB(I,J,K)-VB(I-1,J,K))
     4         /(DX(I,J)+DX(I-1,J)+DX(I,J-1)+DX(I-1,J-1)))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     1        -DT(I,J)*AAM(I,J,K)*2.E0*(VB(I,J+1,K)-VB(I,J,K))/DY(I,J)
C
      XFLUX(I,J,K)
     1    =.25E0*(DY(I,J)+DY(I-1,J)+DY(I,J-1)+DY(I-1,J-1))*XFLUX(I,J,K)
      YFLUX(I,J,K)=DX(I,J)*YFLUX(I,J,K)
  700 CONTINUE
C
C********** DO HORIZ. ADVECTION ************
      DO 400 K=1,KBM1
      DO 400 J=2,JMM1
      DO 400 I=2,IMM1
 400  ADVY(I,J,K)=
     1            +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     2            +YFLUX(I,J,K)-YFLUX(I,J-1,K)
     3   +ARV(I,J)*.25*(CURV(I,J,K)*DT(I,J)*(U(I+1,J,K)+U(I,J,K))
     4             +CURV(I,J-1,K)*DT(I,J-1)*(U(I+1,J-1,K)+U(I,J-1,K)))
C
      RETURN
      END
      SUBROUTINE ADVQ(QB,Q,DTI2,QF)
       INCLUDE 'comblk.h'   
      DIMENSION QB(IM,JM,KB),Q(IM,JM,KB),QF(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
C
C******* HORIZONTAL ADVECTION ************************************
      DO 110 K=2,KBM1
      DO 110 J=2,JM
      DO 110 I=2,IM
      XFLUX(I,J,K)=.125E0*(Q(I,J,K)+Q(I-1,J,K))
     1     *(DT(I,J)+DT(I-1,J))*(U(I,J,K)+U(I,J,K-1))
 110  YFLUX(I,J,K)=.125E0*(Q(I,J,K)+Q(I,J-1,K))*(DT(I,J)+DT(I,J-1))
     1        *(V(I,J,K)+V(I,J,K-1))
C******* HORIZONTAL DIFFUSION ************************************
      DO 315 K=2,KBM1
      DO 315 J=2,JM
      DO 315 I=2,IM
      XFLUX(I,J,K)=XFLUX(I,J,K)
     1    -.25*(AAM(I,J,K)+AAM(I-1,J,K)+AAM(I,J,K-1)+AAM(I-1,J,K-1))
     2     *(H(I,J)+H(I-1,J))*(QB(I,J,K)-QB(I-1,J,K))*DUM(I,J)
     3           /(DX(I,J)+DX(I-1,J))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     1    -.25*(AAM(I,J,K)+AAM(I,J-1,K)+AAM(I,J,K-1)+AAM(I,J-1,K-1))
     2     *(H(I,J)+H(I,J-1))*(QB(I,J,K)-QB(I,J-1,K))*DVM(I,J)
     3           /(DY(I,J)+DY(I,J-1))
      XFLUX(I,J,K)=.5E0*(DY(I,J)+DY(I-1,J))*XFLUX(I,J,K)
      YFLUX(I,J,K)=.5E0*(DX(I,J)+DX(I,J-1))*YFLUX(I,J,K)
  315 CONTINUE
C****** VERTICAL ADVECTION; ADD FLUX TERMS ;THEN STEP FORWARD IN TIME **
      DO 230 K=2,KBM1
      DO 230 J=2,JMM1
      DO 230 I=2,IMM1
      QF(I,J,K)=(W(I,J,K-1)*Q(I,J,K-1)-W(I,J,K+1)*Q(I,J,K+1))
     1                     /(DZ(K)+DZ(K-1))*ART(I,J)
     2                      +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     3                      +YFLUX(I,J+1,K)-YFLUX(I,J,K)
 230  QF(I,J,K)=((H(I,J)+ETB(I,J))*ART(I,J)*QB(I,J,K)-DTI2*QF(I,J,K))
     1             /((H(I,J)+ETF(I,J))*ART(I,J))
C
      RETURN
      END
      SUBROUTINE ADVT(FB,F,FMEAN,DTI2,FF)
C
C     THIS SUBROUTINE INTEGRATES CONSERVATIVE SCALAR EQUATIONS
C
       INCLUDE 'comblk.h'   
      DIMENSION FB(IM,JM,KB),F(IM,JM,KB),FF(IM,JM,KB),FMEAN(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
C
      DO 529 J=1,JM
      DO 529 I=1,IM
      F(I,J,KB)=F(I,J,KBM1)
 529  FB(I,J,KB)=FB(I,J,KBM1)
C
C******* DO ADVECTION FLUXES **************************************
      DO 530 K=1,KBM1
      DO 530 J=2,JM
      DO 530 I=2,IM
      XFLUX(I,J,K)=.25E0*((DT(I,J)+DT(I-1,J))
     2            *(F(I,J,K)+F(I-1,J,K))*U(I,J,K))
 530  YFLUX(I,J,K)=.25E0*((DT(I,J)+DT(I,J-1))
     2            *(F(I,J,K)+F(I,J-1,K))*V(I,J,K))
C******  ADD DIFFUSIVE FLUXES *************************************
C
      DO 99 K=1,KB
      DO 99 J=1,JM
      DO 99 I=1,IM
  99  FB(I,J,K)=FB(I,J,K)-FMEAN(I,J,K)
C
      DO 100 K=1,KBM1
      DO 100 J=2,JM
      DO 100 I=2,IM
      XFLUX(I,J,K)=XFLUX(I,J,K)
     1    -.5E0*(AAM(I,J,K)+AAM(I-1,J,K))*(H(I,J)+H(I-1,J))
     2    *(FB(I,J,K)-FB(I-1,J,K))*DUM(I,J)/(TPRNU*(DX(I,J)+DX(I-1,J)))
      YFLUX(I,J,K)=YFLUX(I,J,K)
     1    -.5E0*(AAM(I,J,K)+AAM(I,J-1,K))*(H(I,J)+H(I,J-1))
     2    *(FB(I,J,K)-FB(I,J-1,K))*DVM(I,J)/(TPRNU*(DY(I,J)+DY(I,J-1)))
      XFLUX(I,J,K)=.5E0*(DY(I,J)+DY(I-1,J))*XFLUX(I,J,K)
      YFLUX(I,J,K)=.5E0*(DX(I,J)+DX(I,J-1))*YFLUX(I,J,K)
  100 CONTINUE
C
      DO 101 K=1,KB
      DO 101 J=1,JM
      DO 101 I=1,IM
 101  FB(I,J,K)=FB(I,J,K)+FMEAN(I,J,K)
C
C****** DO VERTICAL ADVECTION *************************************
      DO 505 J=2,JMM1
      DO 505 I=2,IMM1
 505  FF(I,J,1)=-.5E0*DZR(1)*(F(I,J,1)+F(I,J,2))*W(I,J,2)*ART(I,J)
      DO 520 K=2,KBM1
      DO 520 J=2,JMM1
      DO 520 I=1,IM
 520  FF(I,J,K)=.5E0*DZR(K)*((F(I,J,K-1)+F(I,J,K))*W(I,J,K)
     1                  -(F(I,J,K)+F(I,J,K+1))*W(I,J,K+1))*ART(I,J)
C****** ADD NET HORIZONTAL FLUXES; THEN STEP FORWARD IN TIME **********
      DO 120 K=1,KBM1
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
      FF(I,J,K)=FF(I,J,K)
     1              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     2              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
      FF(I,J,K)=(FB(I,J,K)*(H(I,J)+ETB(I,J))*ART(I,J)-DTI2*FF(I,J,K))
     1                    /((H(I,J)+ETF(I,J))*ART(I,J))
 120  CONTINUE   
      RETURN
      END
      SUBROUTINE ADVU(DRHOX,ADVX,DTI2)
       INCLUDE 'comblk.h'   
      DIMENSION DRHOX(IM,JM,KB),ADVX(IM,JM,KB)                        
C
C  Do vertical advection
      DO 100 K=1,KB  
      DO 100 J=1,JM
      DO 100 I=1,IM
  100 UF(I,J,K)=0.
      DO 140 K=2,KBM1
      DO 140 J=1,JM
      DO 140 I=2,IM
 140  UF(I,J,K)=.25E0*(W(I,J,K)+W(I-1,J,K))*(U(I,J,K)+U(I,J,K-1))
C****COMBINE HOR. and VERT. ADVECTION with
C           -FVD + GDEG/DX + BAROCLINIC TERM **********************
      DO 150 K=1,KBM1
      DO 150 J=2,JMM1
      DO 150 I=2,IMM1
 150  UF(I,J,K)=ADVX(I,J,K)+DZR(K)*(UF(I,J,K)-UF(I,J,K+1))*ARU(I,J)
     1   -ARU(I,J)*.25*(COR(I,J)*DT(I,J)*(V(I,J+1,K)+V(I,J,K))
     2             +COR(I-1,J)*DT(I-1,J)*(V(I-1,J+1,K)+V(I-1,J,K)))
     3        +GRAV*.25E0*(DT(I,J)+DT(I-1,J))
     4        *.5*(EGF(I,J)-EGF(I-1,J)+EGB(I,J)-EGB(I-1,J))
     5        *(DY(I,J)+DY(I-1,J))
     6        +DRHOX(I,J,K)
C******* STEP FORWARD IN TIME ***********************************
      DO 190 K=1,KBM1
      DO 190 J=2,JMM1
      DO 190 I=2,IMM1
 190  UF(I,J,K)=
     1      ((H(I,J)+ETB(I,J)+H(I-1,J)+ETB(I-1,J))*ARU(I,J)*UB(I,J,K)
     2         -2.E0*DTI2*UF(I,J,K))
     3     /((H(I,J)+ETF(I,J)+H(I-1,J)+ETF(I-1,J))*ARU(I,J))
      RETURN
      END
      SUBROUTINE ADVV(DRHOY,ADVY,DTI2)
       INCLUDE 'comblk.h'   
      DIMENSION DRHOY(IM,JM,KB),ADVY(IM,JM,KB)                    
C
C  Do vertical advection
      DO 100 K=1,KB  
      DO 100 J=1,JM
      DO 100 I=1,IM
  100 VF(I,J,K)=0.
      DO 140 K=2,KBM1
      DO 140 J=2,JM
      DO 140 I=1,IM
 140  VF(I,J,K)=.25*(W(I,J,K)+W(I,J-1,K))*(V(I,J,K)+V(I,J,K-1))
C****COMBINE HOR. and VERT. ADVECTION with
C           +FUD + GDEG/DY + BAROCLINIC TERM **********************
      DO 340 K=1,KBM1
      DO 340 J=2,JMM1
      DO 340 I=2,IMM1
 340  VF(I,J,K)=ADVY(I,J,K)+DZR(K)*(VF(I,J,K)-VF(I,J,K+1))*ARV(I,J)
     1   +ARV(I,J)*.25*(COR(I,J)*DT(I,J)*(U(I+1,J,K)+U(I,J,K))
     2             +COR(I,J-1)*DT(I,J-1)*(U(I+1,J-1,K)+U(I,J-1,K)))
     3       +GRAV*.25E0*(DT(I,J)+DT(I,J-1))
     4       *.5*(EGF(I,J)-EGF(I,J-1)+EGB(I,J)-EGB(I,J-1))
     5       *(DX(I,J)+DX(I,J-1))
     6       +DRHOY(I,J,K)
C******* STEP FORWARD IN TIME ***************************************
      DO 390 K=1,KBM1
      DO 390 J=2,JMM1
      DO 390 I=2,IMM1
  390 VF(I,J,K)=
     1     ((H(I,J)+ETB(I,J)+H(I,J-1)+ETB(I,J-1))*ARV(I,J)*VB(I,J,K)
     2          -2.E0*DTI2*VF(I,J,K))
     3    /((H(I,J)+ETF(I,J)+H(I,J-1)+ETF(I,J-1))*ARV(I,J))
      RETURN
      END
      SUBROUTINE BAROPG(DRHOX,DRHOY)
       INCLUDE 'comblk.h'   
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB)
C
      DO 299 K=1,KB
      DO 299 J=1,JM
      DO 299 I=1,IM
 299  RHO(I,J,K)=RHO(I,J,K)-RMEAN(I,J,K)
C
C               X COMPONENT OF BAROCLINIC PRESSURE GRADIENT
C
      DO 300 J=2,JM
      DO 300 I=2,IM
  300 DRHOX(I,J,1)=.5E0*GRAV*(-ZZ(1))*(DT(I,J)+DT(I-1,J))
     2     *(RHO(I,J,1)-RHO(I-1,J,1))
      DO 310 K=2,KBM1
      DO 310 J=2,JM
      DO 310 I=2,IM
 310  DRHOX(I,J,K)=DRHOX(I,J,K-1)
     1      +GRAV*.25E0*(ZZ(K-1)-ZZ(K))*(DT(I,J)+DT(I-1,J))
     2      *(RHO(I,J,K)-RHO(I-1,J,K)+RHO(I,J,K-1)-RHO(I-1,J,K-1))
     3      +GRAV*.25E0*(ZZ(K-1)+ZZ(K))*(DT(I,J)-DT(I-1,J))
     4      *(RHO(I,J,K)+RHO(I-1,J,K)-RHO(I,J,K-1)-RHO(I-1,J,K-1))
C
      DO 360 K=1,KBM1
      DO 360 J=2,JM
      DO 360 I=2,IM
 360  DRHOX(I,J,K)=.25E0*(DT(I,J)+DT(I-1,J))*DRHOX(I,J,K)*DUM(I,J)
     1     *(DY(I,J)+DY(I-1,J))
C
C               Y COMPONENT OF BAROCLINIC PRESSURE GRADIENT
C
      DO 500 J=2,JM
      DO 500 I=2,IM
 500  DRHOY(I,J,1)=.5E0*GRAV*(-ZZ(1))*(DT(I,J)+DT(I,J-1))
     1     *(RHO(I,J,1)-RHO(I,J-1,1))
      DO 510 K=2,KBM1
      DO 510 J=2,JM
      DO 510 I=2,IM
 510  DRHOY(I,J,K)=DRHOY(I,J,K-1)
     1        +GRAV*.25E0*(ZZ(K-1)-ZZ(K))*(DT(I,J)+DT(I,J-1))
     2        *(RHO(I,J,K)-RHO(I,J-1,K)+RHO(I,J,K-1)-RHO(I,J-1,K-1))
     3        +GRAV*.25E0*(ZZ(K-1)+ZZ(K))*(DT(I,J)-DT(I,J-1))
     4        *(RHO(I,J,K)+RHO(I,J-1,K)-RHO(I,J,K-1)-RHO(I,J-1,K-1))
C
      DO 560 K=1,KBM1
      DO 560 J=2,JM
      DO 560 I=2,IM
 560  DRHOY(I,J,K)=.25E0*(DT(I,J)+DT(I,J-1))*DRHOY(I,J,K)*DVM(I,J)
     1       *(DX(I,J)+DX(I,J-1))
C
      DO 561 K=1,KBM1
      DO 561 J=2,JM
      DO 561 I=2,IM
      DRHOX(I,J,K)=RAMP*DRHOX(I,J,K)
 561  DRHOY(I,J,K)=RAMP*DRHOY(I,J,K)
C
      DO 571 K=1,KB
      DO 571 J=1,JM
      DO 571 I=1,IM
 571  RHO(I,J,K)=RHO(I,J,K)+RMEAN(I,J,K)
C
      RETURN
      END
      SUBROUTINE BCOND(IDX)
       INCLUDE 'comblk.h'   
C
C  Closed boundary conditions are automatically enabled through
C  specification of the masks, DUM, DVM and FSM. For problems with
C  open boundaries, dependent variables must be specified on the
C  the boundary grid points. Sample boundary conditions, left over
C  from a problem with three open boundaries are to be seen below
C  but are commented out for the closed basin application.
C
C              SPECIFICATION OF OPEN BOUNDARY CONDITIONS
C        BRACKETS OF *  INDICATE INTERIOR (NON-BOUNDARY) GRID POINTS.
C        HORIZONTAL LOCATIONS OF EL, T AND S ARE COINCIDENT.
C
C        I-1           I         I+1
C
C
C      U(JM)       EL(JM)  =.   U(JM)           BC
C                           .
C                   V(JM)   .                   BC
C                           .
C     *U(JMM1)*   *EL(JMM1)*   *U(JMM1)*
C
C                 *V(JMM1)*
C                       .                         *V(IMM1)*           V(IM)
C                       .
C                       .              *U(IMM1)* *EL(IMM1)*   U(IM)   EL(IM)
C                       .                              ". . . . . . . ."
C                       .
C                       .                         *V(IMM1)*           V(IM)
C                       .
C                       .                                     BC        BC
C
C                    *V(3)*
C
C       *U(2)*       *EL(2)*   *U(2)*
C                           .
C                     V(2)  .             BC
C                           .
C         U(1)        EL(1)=.   U(1)      BC
C --------------------------------------------------------------------
C                     V(1) IS NOT USED
C
C
        real dval(6)
        real dvale(6)
        real dvalt(6)
c  this dval is for the internal northern boundary vel-
c    2.77 days (80.) to 25. days (720.) using Tau=2*deltat*dval
        data dval /80.,120.,240.,360.,540.,720./
c--- this is for temperature
c        data dvalt /80.,120.,240.,360.,540.,720./
        data dvalt /8.,12.,24.,36.,54.,72./
c  this dval is for the external vel-
c    2.77 days (4000.) to 25. days (36000.) using Tau=2*deltat*dval
        data dvale /4000.,6000.,12000.,18000.,27000.,36000./

        damp(ad,f,fo)=(ad*f+fo)/(1.+ad)
C
      PI=3.14167E0
      GEE=9.807E0
c
      GO TO (10,20,30,40,50,60), IDX
C
 10   CONTINUE
C-----------------------------------------------------------------------
C                   EXTERNAL ELEV. B.C.'S
C-----------------------------------------------------------------------
      do j=1,jm
        elf(im,j)=elf(imm1,j)
        elf(1,j)=elf(2,j)
      end do
      do i=1,im
        elf(i,jm)=elf(i,jmm1)
        elf(i,1)=elf(i,2)
      end do

      do j=1,jm
         do i=1,im
           elf(i,j)=elf(i,j)*fsm(i,j)
         end do
      end do

      return

 20   CONTINUE
C-----------------------------------------------------------------------
C                   EXTERNAL VEL B.C.'S
C-----------------------------------------------------------------------
C                   ***** SEAWARD ********
      DO 120 J=1,JM
      UAF(IM,J)=RAMP*UABE(J)+COVRHE(J)*(EL(IMM1,J)-ELBE(J))
      UAF(2,J)=RAMP*UABW(J)-COVRHW(J)*(EL(2,J)-ELBW(J))
 120  continue
      UMID=UA(IM,1)
      VAF(IM,1)=VA(IM,1)-DTI/(DX(IM,1)+DX(IMM1,1))*
     1  (  (UMID+ABS(UMID))*(VA(IM,1)-VA(IMM1,1))+
     2  (UMID-ABS(UMID))*(0.E0 -VA(IM,1))  )
      UMID=UA(2,1)
      VAF(1,1)=VA(1,1)+DTI/(DX(1,1)+DX(2,1))*
     1  (  (UMID-ABS(UMID))*(VA(1,1)-VA(2,1))+
     2  (UMID+ABS(UMID))*(0.E0 -VA(1,1))  )
      do j=2,jmm1
        UMID=.5E0*(UA(IM,J)+UA(IM,J-1))
        VAF(IM,J)=VA(IM,J)-DTI/(DX(IM,J)+DX(IMM1,J))*
     1  (  (UMID+ABS(UMID))*(VA(IM,J)-VA(IMM1,J))+
     2  (UMID-ABS(UMID))*(0.E0 -VA(IM,J))  )
        UMID=.5E0*(UA(2,J)+UA(2,J-1))
        VAF(1,J)=VA(1,J)+DTI/(DX(1,J)+DX(2,J))*
     1  (  (UMID-ABS(UMID))*(VA(1,J)-VA(2,J))+
     2  (UMID+ABS(UMID))*(0.E0 -VA(1,J))  )
      end do
c         
C                   **** NORTH AND SOUTH ****
c
      DO 124 I=1,IMM1
      vaf(i,2)=vabs(i)-covrhs(i)*(el(i,2)-elbs(i))
  124 VAF(I,JM)=RAMP*VABN(I)+COVRHN(I)*(EL(I,JMM1)-ELBN(I))
c
      VMID=VA(1,JM)
      UAF(1,JM)=UA(1,JM)-DTI/(DY(1,JM)+DY(1,JMM1))
     1    *(   (VMID+ABS(VMID))*(UA(1,JM)-UA(1,JMM1))
     2        +(VMID-ABS(VMID))*(0.E0-UA(1,JM))  )
      VMID=VA(1,2)
      UAF(1,1)=UA(1,1)-DTI/(DY(1,1)+DY(1,2))
     1    *(   (VMID+ABS(VMID))*(UA(1,1)-0.E0)
     2        +(VMID-ABS(VMID))*(UA(1,2)-UA(1,1)) )
      DO I=2,IMM1
        VMID=.5E0*(VA(I,JM)+VA(I-1,JM))
        UAF(I,JM)=UA(I,JM)-DTI/(DY(I,JM)+DY(I,JMM1))*
     1  (   (VMID+ABS(VMID))*(UA(I,JM)-UA(I,JMM1))+
     2  (VMID-ABS(VMID))*(0.E0-UA(I,JM))  )
        VMID=.5E0*(VA(I,2)+VA(I-1,2))
        UAF(I,1)=UA(I,1)-DTI/(DY(I,1)+DY(I,2))*
     1  (   (VMID+ABS(VMID))*(UA(I,1)-0.E0)+
     2  (VMID-ABS(VMID))*(UA(I,2)-UA(I,1)) )
      end do
C
        do i=1,4
          do j=i+2,jm-i
            vaf(i+1,j)=damp(dvale(i),vaf(i+1,j),0.)
            vaf(im-i,j)=damp(dvale(i),vaf(im-i,j),0.)
          enddo
        enddo
        do i=1,4
          do j=i+2,jmm1-i
            uaf(i+2,j)=damp(dvale(i),uaf(i+2,j),uabw(j))
            uaf(im-i,j)=damp(dvale(i),uaf(im-i,j),uabe(j))
          enddo
        enddo
        do j=1,4
          do i=j+2,im-j
            uaf(i,j+1)=damp(dvale(j),uaf(i,j+1),0.)
            uaf(i,jm-j)=damp(dvale(j),uaf(i,jm-j),0.)
          enddo
        enddo
        do j=1,4
          do i=j+2,imm1-j
            vaf(i,j+2)=damp(dvale(j),vaf(i,j+2),vabs(i))
            vaf(i,jm-j)=damp(dvale(j),vaf(i,jm-j),vabn(i))
          enddo
        enddo
c    
      DO 131 J=1,JM
      DO 131 I=1,IM
      UAF(I,J)=UAF(I,J)*DUM(I,J)
 131  VAF(I,J)=VAF(I,J)*DVM(I,J)
      RETURN
C
 30   CONTINUE
C-----------------------------------------------------------------------
C                   INTERNAL VEL B.C.'S
C-----------------------------------------------------------------------
C                         **** EAST AND WEST *******
      DO K=1,KBM1
        DO J=2,JMM1
          GA=SQRT(H(IM,J)/3000.E0)
          UF(IM,J,K)=
     1    GA*(.25E0*U(IMM1,J-1,K)+.5E0*U(IMM1,J,K)+.25E0*U(IMM1,J+1,K))+
     2    (1.E0-GA)*(.25E0*U(IM,J-1,K)+.5E0*U(IM,J,K)+.25E0*U(IM,J+1,K))
        END DO
        UF(IM,JM,K)=0.E0
        UF(IM,1,K)=0.E0
        DO J=2,JM
          UMID=.5E0*(U(IM,J,K)+U(IM,J-1,K))
          VF(IM,J,K)=V(IM,J,K)-DTI/(DX(IM,J)+DX(IMM1,J))*
     1    (  (UMID+ABS(UMID))*(V(IM,J,K)-V(IMM1,J,K))+
     2    (UMID-ABS(UMID))*(0.E0 -V(IM,J,K))  )
        END DO
        UMID=U(IM,1,K)
        VF(IM,1,K)=V(IM,1,K)-DTI/(DX(IM,1)+DX(IMM1,1))*
     1  (  (UMID+ABS(UMID))*(V(IM,1,K)-V(IMM1,1,K))+
     2  (UMID-ABS(UMID))*(0.E0 -V(IM,1,K))  )
      END DO
c
      do K=1,KBM1
        do J=2,JMM1
          GA=SQRT(H(2,J)/3000.E0)
          UF(2,J,K)
     1    =GA*(.25E0*U(3,J-1,K)+.5E0*U(3,J,K)+.25E0*U(3,J+1,K))
     2    +(1.E0-GA)*(.25E0*U(2,J-1,K)+.5E0*U(2,J,K)+.25E0*U(2,J+1,K))
        end do
        UF(2,JM,K)=0.E0
        UF(2,1,K)=0.E0
        do J=2,JM
          UMID=.5E0*(U(2,J,K)+U(2,J-1,K))
          VF(1,J,K)=V(1,J,K)-DTI/(DX(1,J)+DX(2,J))*
     1    (  (UMID+ABS(UMID))*(V(1,J,K)-0.E0)+
     2    (UMID-ABS(UMID))*(V(2,J,K)-V(1,J,K))  )
        end do
        UMID=U(2,1,K)
        VF(1,1,K)=V(1,1,K)-DTI/(DX(1,1)+DX(2,1))*
     1  (  (UMID+ABS(UMID))*(V(1,1,K)-0.E0)+
     2  (UMID-ABS(UMID))*(V(2,1,K)-V(1,1,K))  )
      end do


C                       ***** NORTH AND SOUTH ***
      DO 158 K=1,KBM1
      DO 150 I=2,IMM1
      GA=SQRT(H(I,JMM1)/3000.E0)
      VF(I,JM,K)=
     1   GA*(.25E0*V(I-1,JMM1,K)+.5E0*V(I,JMM1,K)+.25E0*V(I+1,JMM1,K))
     2 + (1.E0-GA)
     3 *(.25E0*V(I-1,JM,K)+.5E0*V(I,JM,K)+.25E0*V(I+1,JM,K))
  150 CONTINUE
      DO 152 I=2,IMM1
      GA=SQRT(H(I,2)/4000.E0)
      VF(I,2,K)=GA*(.25E0*V(I-1,3,K)+.5E0*V(I,3,K)+.25E0*V(I+1,3,K))
     1    +(1.E0-GA)*(.25E0*V(I-1,2,K)+.5E0*V(I,2,K)+.25E0*V(I+1,2,K))
  152 VF(I,1,K)=VF(I,2,K)
      DO I=2,IMM1
        VMID=.5E0*(V(I,JM,K)+V(I-1,JM,K))
        UF(I,JM,K)=U(I,JM,K)-DTI/(DY(I,JM)+DY(I,JMM1))*
     1  (   (VMID+ABS(VMID))*(U(I,JM,K)-U(I,JMM1,K))+
     2  (VMID-ABS(VMID))*(0.E0-U(I,JM,K))  )
        VMID=.5E0*(V(I,2,K)+V(I-1,2,K))
        UF(I,1,K)=U(I,1,K)-DTI/(DY(I,1)+DY(I,2))*
     1  (   (VMID+ABS(VMID))*(U(I,1,K)-0.E0)+
     2  (VMID-ABS(VMID))*(U(I,2,K)-U(I,1,K)) )
      END DO
      VMID=V(1,JM,K)
      UF(1,JM,K)=U(1,JM,K)-DTI/(DY(1,JM)+DY(1,JMM1))*
     1(   (VMID+ABS(VMID))*(U(1,JM,K)-U(1,JMM1,K))+
     2(VMID-ABS(VMID))*(0.E0-U(1,JM,K))  )
      VMID=V(1,2,K)
      UF(1,1,K)=U(1,1,K)-DTI/(DY(1,1)+DY(1,2))*
     1(   (VMID+ABS(VMID))*(U(1,1,K)-0.E0)+
     2(VMID-ABS(VMID))*(U(1,2,K)-U(1,1,K)) )
          do j=1,4
            do i=j+2,imm1-j
              vf(i,jm-j,k)=damp(dval(j),vf(i,jm-j,k),v(i,jm,k))
            enddo
          enddo
          do j=1,4
            do i=j+2,imm1-j
              vf(i,j+2,k)=damp(dval(j),vf(i,j+2,k),v(i,2,k))
            enddo
          enddo
  158 CONTINUE
C                    **********************
      DO 160 K=1,KBM1
C    
      DO 160 J=1,JM
      DO 160 I=1,IM
      UF(I,J,K)=UF(I,J,K)*DUM(I,J)
      VF(I,J,K)=VF(I,J,K)*DVM(I,J)
  160 CONTINUE
C      
      RETURN
C
 40   CONTINUE
C-----------------------------------------------------------------------
C                   TEMP & SAL B.C.'S
C-----------------------------------------------------------------------
C                         ***** EAST AND WEST *****
      DO 220 K=1,KBM1
      DO 220 J=1,JM
      UF(IM,J,K)=T(IM,J,K)-DTI/(DX(IM,J)+DX(IMM1,J))
     1     *(  (U(IM,J,K)+ABS(U(IM,J,K)))*(T(IM,J,K)-T(IMM1,J,K))
     2        +(U(IM,J,K)-ABS(U(IM,J,K)))*(TBE(J,K) -T(IM,J,K))  )
      VF(IM,J,K)=S(IM,J,K)-DTI/(DX(IM,J)+DX(IMM1,J))
     1     *(  (U(IM,J,K)+ABS(U(IM,J,K)))*(S(IM,J,K)-S(IMM1,J,K))
     2        +(U(IM,J,K)-ABS(U(IM,J,K)))*(SBE(J,K) -S(IM,J,K))  )
      UF(1,J,K)=T(1,J,K)-DTI/(DX(1,J)+DX(2,J))
     1     *(  (U(2,J,K)+ABS(U(2,J,K)))*(T(1,J,K)-TBW(J,K))
     2        +(U(2,J,K)-ABS(U(2,J,K)))*(T(2,J,K)-T(1,J,K))  )
      VF(1,J,K)=S(1,J,K)-DTI/(DX(1,J)+DX(2,J))
     1     *(  (U(2,J,K)+ABS(U(2,J,K)))*(S(1,J,K)-SBW(J,K))
     2        +(U(2,J,K)-ABS(U(2,J,K)))*(S(2,J,K) -S(1,J,K))  )
  220 CONTINUE
      do k=1,kbm1
        do i=1,6
          do j=i,jm-i
            vf(im-i,j,k)=damp(dvalt(i),vf(im-i,j,k),sbin(im-i,j,k))
            uf(im-i,j,k)=damp(dvalt(i),uf(im-i,j,k),tbin(im-i,j,k))
          enddo
          do j=i,jm-i
            vf(i,j,k)=damp(dvalt(i),vf(i,j,k),sbin(i,j,k))
            uf(i,j,k)=damp(dvalt(i),uf(i,j,k),tbin(i,j,k))
          enddo
        enddo
      end do
c
      DO 230 K=1,KBM1
      DO 230 I=1,IM
      UF(I,JM,K)=T(I,JM,K)-DTI/(DY(I,JM)+DY(I,JMM1))
     1    *(   (V(I,JM,K)+ABS(V(I,JM,K)))*(T(I,JM,K)-T(I,JMM1,K))
     2        +(V(I,JM,K)-ABS(V(I,JM,K)))*(TBN(I,K)-T(I,JM,K))  )
      VF(I,JM,K)=S(I,JM,K)-DTI/(DY(I,JM)+DY(I,JMM1))
     1    *(   (V(I,JM,K)+ABS(V(I,JM,K)))*(S(I,JM,K)-S(I,JMM1,K))
     2        +(V(I,JM,K)-ABS(V(I,JM,K)))*(SBN(I,K)-S(I,JM,K))  )
      UF(I,1,K)=T(I,1,K)-DTI/(DY(I,1)+DY(I,2))
     1    *(   (V(I,2,K)+ABS(V(I,2,K)))*(T(I,1,K)-TBS(I,K))
     2        +(V(I,2,K)-ABS(V(I,2,K)))*(T(I,2,K)-T(I,1,K)) )
      VF(I,1,K)=S(I,1,K)-DTI/(DY(I,1)+DY(I,2))
     1    *(   (V(I,2,K)+ABS(V(I,2,K)))*(S(I,1,K)-SBS(I,K))
     2        +(V(I,2,K)-ABS(V(I,2,K)))*(S(I,2,K)-S(I,1,K)) )
  230 CONTINUE
      do k=1,kbm1
        do j=1,6
          do i=j,im-j
            uf(i,j,k)=damp(dvalt(j),uf(i,j,k),tbin(i,j,k))
            vf(i,j,k)=damp(dvalt(j),vf(i,j,k),sbin(i,j,k))
          enddo
          do i=j,im-j
            uf(i,jm-j,k)=damp(dvalt(j),uf(i,jm-j,k),tbin(i,jm-j,k))
            vf(i,jm-j,k)=damp(dvalt(j),vf(i,jm-j,k),sbin(i,jm-j,k))
          enddo
        enddo
      enddo
C    
      DO 240 K=1,KBM1
      DO 240 J=1,JM
      DO 240 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
 
 240  CONTINUE
      RETURN
C
C
 50   CONTINUE
C---------------VERTICAL VEL. B. C.'S --------------------------------
      DO 250 K=1,KBM1
      DO 250 J=1,JM
      DO 250 I=1,IM
      W(I,J,K)=W(I,J,K)*FSM(I,J)
 250  CONTINUE
      DO 251 K=1,KB
      DO 251 J=1,JM
 251  W(IM,J,K)=0.E0
C    
      RETURN
C
 60   CONTINUE
C---------------- Q2 AND Q2L B.C.'S -----------------------------------
      DO 300 K=1,KB
C     
      DO 295 J=1,JM
      UF(IM,J,K)=1.E-10
 295  VF(IM,J,K)=1.E-10
      DO 296 I=1,IM
      UF(I,JM,K)=1.E-10
      VF(I,JM,K)=1.E-10
      UF(I,1,K)=1.E-10
 296  VF(I,1,K)=1.E-10
      DO 297 J=1,JM
      DO 297 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 297  VF(I,J,K)=VF(I,J,K)*FSM(I,J)
 300  CONTINUE
      RETURN
      END
      SUBROUTINE BDAMP
C=========================================================
C Damp the boundary of model domain using Dr. Ginis's idea
C Written by Erxuan Fu, GSO,URI,11/3/94
C=========================================================
      INCLUDE 'comblk.h'
C
C 
      DIMENSION TC(IM),TCJ(JM),A0(IM,JM)
      DAMP(AD,F,FO)=(AD*F+FO)/(1.0+AD)
      DATA IONCE/1/,A0/LIJ*0.E0/
C do only once! 
      IF(IONCE.EQ.1) THEN
      IONCE=0
C
      PRINT*,'INSIDE (IF) in BDAMP'
C
      DO 1215 J=1,JM
1215  TCJ(J)=1.0E30
      DO 1216 I=1,IM
1216  TC(I)=1.0E30
C 
      TC(IM)=20.
      TC(IM-1)=20.
      TC(IM-2)=40.
      TC(IM-3)=60.
      TC(IM-4)=90.
      TC(IM-5)=120.
      TC(1)=20.
      TC(2)=20.
      TC(3)=40.
      TC(4)=60.
      TC(5)=90.
      TC(6)=120.
      TCJ(JM)=20.
      TCJ(JM-1)=20.
      TCJ(JM-2)=40.
      TCJ(JM-3)=60.
      TCJ(JM-4)=90.
      TCJ(JM-5)=120.
      TCJ(1)=20.
      TCJ(2)=20.
      TCJ(3)=40.
      TCJ(4)=60.
      TCJ(5)=90.
      TCJ(6)=120.
      DO 1217 I=1,IM
      DO 1217 J=1,JM
      A0(I,J) = 2.0*AMIN1(TC(I),TCJ(J))
1217  CONTINUE
      ENDIF
C
C      WRITE(6,'(''Coeficient A along J=(JM+1)/2'')')
C      WRITE(6,'(110E12.4)') (A(I,(JM+1)/2),I=1,IM)
C
      DO 1 I=1,IM
      DO 1 J=1,JM
      ELF(I,J)=DAMP(A0(I,J),ELF(I,J),0.)
      UAF(I,J)=DAMP(A0(I,J),UAF(I,J),0.)
      VAF(I,J)=DAMP(A0(I,J),VAF(I,J),0.)
  1   CONTINUE
C
      RETURN
      END
      SUBROUTINE CLIMAT2MODEL(SB1,TB2)
c---------------------------------------------------------------
c      This subroutine calculates initial fiels for model from
c      climatological temperature and salinity fields (SB1,TB2)
c      and streamfunction (phi)
c---------------------------------------------------------------
      INCLUDE 'comblk.h'
c
      real SB1(IM,JM,33),TB2(IM,JM,33)
      real DRHOX(IM,JM,KB),DRHOY(IM,JM,KB)
      real TMPT(33),TMPTI(KB),TMPS(33),TMPSI(KB),ZI(KB)
      real yvec(jm),tvec(jm)
      REAL CMP,L1,L2,LATMIN,LATMAX,LONGMIN,LONGMAX
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
      common/mean/TBM(33),SBM(33),RHM(33),DB(33)
c
      print *,'   '
      print *,' begin CLIMAT2MODEL'
c
      mid=int(im/2)
      mjd=int((40.0-latmin)*jm/(latmax-latmin))
C
c---------------- calculate SSE -------------------------------
c
      do j=1,jm
         do i=1,im
            ELB(i,j)=0.0
            etb(i,j)=elb(i,j)
            el(i,j)=elb(i,j)
            et(i,j)=elb(i,j)
            D(i,j)=H(i,j)+ELB(i,j)
            DT(i,j)=H(i,j)+ELB(i,j)
         end do
      end do
c
c-------------- INTERPOLATING INTO SIGMA COORDINATES ----------------
c
      DO J=1,JM
c        print *,' check verZ2SIGsp  j=',j
         DO I=1,IM
                 IF(FSM(I,J).EQ.1.) THEN
c         if(j.eq.48) print *,' check verZ2SIGsp i,j=',i,j
c
            DO K=1,KB
            ZI(K)=-D(I,J)*ZZ(K)
            END DO
C
            DO K=1,33
            TMPT(K)=TB2(I,J,K)
            TMPS(K)=SB1(I,J,K)
            END DO
C
c----------- 02-04-02 use SP interpolation
c
      call verZ2SIGsp(33,KB,DB,TMPT,ZI,TMPTI,i,j)
      call verZ2SIGsp(33,KB,DB,TMPS,ZI,TMPSI,i,j)
c     CALL VERINTERP(33,KB,DB,TMPT,ZI,TMPTI)
c     CALL VERINTERP(33,KB,DB,TMPS,ZI,TMPSI)
c
c----------- 02-04-02 remove instability in t
c
            DO K=2,KB
             if(TMPTI(k).gt.TMPTI(k-1)) TMPTI(k)=TMPTI(k-1)
            END DO
C
            DO K=1,KB
            TB(I,J,K)=TMPTI(k)
            SB(I,J,K)=TMPSI(k)
            END DO
C
      call verZ2SIGsp(33,KB,DB,TBM,ZI,TMPTI,i,j)
      call verZ2SIGsp(33,KB,DB,SBM,ZI,TMPSI,i,j)
c     CALL VERINTERP(33,KB,DB,TBM,ZI,TMPTI)
c     CALL VERINTERP(33,KB,DB,SBM,ZI,TMPSI)
C
            DO K=1,KB
            TMEAN(I,J,K)=TMPTI(k)
            SMEAN(I,J,K)=TMPSI(k)
            END DO
C
                 ELSE
C
            DO K=1,KB
            TB(I,J,K)=0
            SB(I,J,K)=0
            TMEAN(I,J,K)=0
            SMEAN(I,J,K)=0
            END DO
C
                 END IF 
C
         END DO
      END DO
c--------------------------------------------------------
      CALL DENS(SB,TB,RHO)
      CALL DENS(SMEAN,TMEAN,RMEAN)
      CALL BAROPG(DRHOX,DRHOY)
      do k=1,kbm1
         do j=2,jm
            do i=2,im
               drhox(i,j,k)=drhox(i,j,k)+.25E0*GRAV*(DY(I,J)+DY(I-1,J))
     1         *(D(I,J)+D(I-1,J))*(EL(I,J)-EL(I-1,J))*DUM(i,j)
               drhoy(i,j,k)=drhoy(i,j,k)+.25E0*GRAV*(DX(I,J)+DX(I,J-1))
     1         *(D(I,J)+D(I,J-1))*(EL(I,J)-EL(I,J-1))*DVM(i,j)
               drhox(1,j,k)=drhox(2,j,k)
               drhoy(1,j,k)=drhoy(2,j,k)
               drhox(i,1,k)=drhox(i,2,k)
               drhoy(i,1,k)=drhoy(i,2,k)
            end do
         end do
         drhox(1,1,k)=0.5*(drhox(1,2,k)+drhox(2,1,k))
         drhoy(1,1,k)=0.5*(drhoy(1,2,k)+drhoy(2,1,k))
      end do
      small=1.e-15
      do k=1,kbm1
         do j=2,jm-1
            do i=2,im-1
               cmp=dum(i+1,j)+dum(i,j)+dum(i,j-1)+dum(i+1,j-1)+small
               vb(i,j,k)= 1/aru(i,j)*( drhox(i+1,j,k)+drhox(i,j,k)+
     1         drhox(i,j-1,k)+drhox(i+1,j-1,k) )/cmp/d(i,j)/cor(i,j)
               cmp=dvm(i-1,j)+dvm(i,j)+dvm(i,j+1)+dvm(i-1,j+1)+small
               ub(i,j,k)=-1/arv(i,j)*( drhoy(i-1,j,k)+drhoy(i,j,k)+
     1         drhoy(i,j+1,k)+drhoy(i-1,j+1,k) )/cmp/d(i,j)/cor(i,j)
            end do
            cmp=dvm(im-1,j)+dvm(im,j)+dvm(im,j+1)+dvm(im-1,j+1)+small
            ub(im,j,k)=-1/arv(im,j)*( drhoy(im-1,j,k)+drhoy(im,j,k)+
     1      drhoy(im,j+1,k)+drhoy(im-1,j+1,k) )/cmp/d(im,j)/cor(im,j)            
            cmp=dum(im,j)+dum(im,j-1)+small
            vb(im,j,k)= 1/aru(im,j)*( drhox(im,j,k)+
     1      drhox(im,j-1,k) )/cmp/d(im,j)/cor(im,j)
         end do
         do i=2,im-1
            cmp=dum(i+1,jm)+dum(i,jm)+dum(i,jm-1)+dum(i+1,jm-1)+small
            vb(i,jm,k)= 1/aru(i,jm)*( drhox(i+1,jm,k)+drhox(i,jm,k)+
     1      drhox(i,jm-1,k)+drhox(i+1,jm-1,k) )/cmp/d(i,jm)/cor(i,jm)
            cmp=dvm(i-1,jm)+dvm(i,jm)+small
            ub(i,jm,k)=-1/arv(i,jm)*( drhoy(i-1,jm,k)+drhoy(i,jm,k)
     1       )/cmp/d(i,jm)/cor(i,jm)
         end do
c     print*,'drhoy=',cmp,h(2,jm),elb(2,jm),arv(2,jm)
         ub(1,1,k)=0.5*(ub(1,2,k)+ub(2,1,k))
         ub(1,jm,k)=0.5*(ub(1,jm-1,k)+ub(2,jm,k))
         ub(im,1,k)=0.5*(ub(im,2,k)+ub(im-1,1,k))
         ub(im,jm,k)=0.5*(ub(im,jm-1,k)+ub(im-1,jm,k))
         vb(1,1,k)=0.5*(vb(1,2,k)+vb(2,1,k))
         vb(1,jm,k)=0.5*(vb(1,jm-1,k)+vb(2,jm,k))
         vb(im,1,k)=0.5*(vb(im,2,k)+vb(im-1,1,k))
         vb(im,jm,k)=0.5*(vb(im,jm-1,k)+vb(im-1,jm,k))
         do j=1,jm
            do i=1,im
               ub(i,j,k)=ub(i,j,k)*dum(i,j)
               vb(i,j,k)=vb(i,j,k)*dvm(i,j)
            end do
         end do
      end do
c----------- falk 09/28/00 correct u,v that u=v=0 at lev kb-1
c     do k=1,kbm1
c      do j=1,jm
c       do i=1,im
c        ub(i,j,k)=ub(i,j,k)-ub(i,j,kbm1)
c        vb(i,j,k)=vb(i,j,k)-vb(i,j,kbm1)
c       end do
c      end do
c     end do
c
c---------- 01-23-02 do not use calculation of u,v below 2000m
c
      do i=1,im
      do j=1,jm
       if(D(i,j).gt.2000.) then
        do k=1,kb
         if(-D(i,j)*ZZ(k).gt.2000.) then
          k2000=k-1
          go to 987
         end if
        end do
 987    continue
        do k=k2000+1,kb
         ub(i,j,k)=ub(i,j,k2000)
         vb(i,j,k)=vb(i,j,k2000)
        end do
       end if
      end do
      end do
c
c------  Calculating vertically averaged velocities  ------------
c
      do j=1,jm
         do i=1,im
            uab(i,j)=0.
            vab(i,j)=0.
            cmp=0.
            do k=1,kb-1
               uab(i,j)=uab(i,j)+ub(i,j,k)*dz(k)
               vab(i,j)=vab(i,j)+vb(i,j,k)*dz(k)
               cmp=cmp+dz(k)
            end do
            uab(i,j)=uab(i,j)/cmp
            vab(i,j)=vab(i,j)/cmp
         end do
      end do
c------ falk 12/04/00 check south bnd
c     print *,' check south bnd (vab(i,2),i=1,im)'
c     write(6,1101) (vab(i,2),i=1,im)
 1101 format(10f7.2)
c
c-------------- Adjusting the velocities for transports ------------
c
      CALL TRANSPORTS
c-------------------------------------------------------------------
      print*,'Finished initialization ...'
c
      return
      end
      SUBROUTINE CURL(FX,FY,DX,DY,IM,JM,CF)  
      DIMENSION FX(IM,JM),FY(IM,JM),DX(IM,JM),DY(IM,JM)   
      DIMENSION CF(IM,JM) 
      DO I=3,IM-1      
      DO J=3,JM-1      
      CF(I,J)=-FX(I,J)*(DX(I,J)+DX(I-1,J))  
     1       +FX(I,J-1)*(DX(I,J-1)+DX(I-1,J-1))  
     2       +FY(I,J)*(DY(I,J)+DY(I,J-1))  
     1       -FY(I-1,J)*(DY(I-1,J)+DY(I-1,J-1))  
      AREA   =   0.25*(DX(I,J)+DX(I-1,J)  
     1                +DX(I,J-1)+DX(I-1,J-1))     
     2          *0.25*(DY(I,J)+DY(I,J-1)  
     1                +DY(I-1,J)+DY(I-1,J-1))  
      CF(I,J)=CF(I,J)/AREA    
      ENDDO       
      ENDDO       
      RETURN   
      END    
	SUBROUTINE DATE2DAY(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,hour
      real julday
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

      year=int(date/1000000.)
      month=nint(100*(date/1000000.-int(date/1000000.)))
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

      return
      end
	SUBROUTINE DAY2DATE(year,julday,date)
      integer*4 date
      integer dat2day(12),dat2dayl(12),day,month,year,year1,hour
      real julday,julday1
      real*8 tmp
      data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
      data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/

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
      hour=nint((julday1-int(julday1))*24.)
      date=year1*1000000+month*10000+day*100+hour
c
      return
      end
      SUBROUTINE DENS(SI,TI,RHOO)
       INCLUDE 'comblk.h'
      DIMENSION SI(IM,JM,KB),TI(IM,JM,KB),RHOO(IM,JM,KB)
      REAL*8  TR,SR,P,RHOR,CR
C       If using 32 bit precision, it is recommended that
C       TR, SR, P, RHOR , CR be made double precision.
C
C         THIS SUBROUTINE COMPUTES DENSITY- 1.025
C         T = POTENTIAL TEMPERATURE
C

      DO 1 K=1,KBM1
      DO 1 J=1,JM
      DO 1 I=1,IM
      TR=TI(I,J,K)+TBIAS
      SR=SI(I,J,K)+SBIAS
C            Approximate pressure in units of bars
      P=-GRAV*1.025*ZZ(K)*DT(I,J)*0.01
C
      RHOR = 999.842594 + 6.793952E-2*TR
     $        - 9.095290E-3*TR**2 + 1.001685E-4*TR**3
     $        - 1.120083E-6*TR**4 + 6.536332E-9*TR**5
C
      RHOR = RHOR + (0.824493 - 4.0899E-3*TR
     $       + 7.6438E-5*TR**2 - 8.2467E-7*TR**3
     $       + 5.3875E-9*TR**4) * SR
     $       + (-5.72466E-3 + 1.0227E-4*TR
     $       - 1.6546E-6*TR**2) * SR**1.5
     $       + 4.8314E-4 * SR**2
C
      CR=1449.1+.0821*P+4.55*TR-.045*TR**2
     $                              +1.34*(SR-35.)
      RHOR=RHOR + 1.E5*P/CR**2
     $     *(1.-2.*P/CR**2)
C
      RHOO(I,J,K)=(RHOR-1025.)*1.E-3*FSM(I,J)
    1 CONTINUE
C
      DO 3 J=1,JM
      DO 3 I=1,IM
    3 RHOO(I,J,KB)=1.
      RETURN
      END
      SUBROUTINE DENS1(DB,SI,TI)
      include 'comblk.h'
      DIMENSION SI(IM,JM,33),TI(IM,JM,33),RHOO(IM,JM,33)
      DIMENSION DB(33)
      REAL*8  TR,SR,P,RHOR,CR
C       If using 32 bit precision, it is recommended that
C       TR, SR, P, RHOR , CR be made double precision.
C
C         THIS SUBROUTINE COMPUTES DENSITY- 1.025
C         T = POTENTIAL TEMPERATURE
C
      DO 1 K=1,33
      DO 1 J=1,JM
      DO 1 I=1,IM
      TR=TI(I,J,K)
      SR=SI(I,J,K)
C            Approximate pressure in units of bars
      P=GRAV*1.025*DB(K)*0.01
C
      RHOR = 999.842594 + 6.793952E-2*TR
     $        - 9.095290E-3*TR**2 + 1.001685E-4*TR**3
     $        - 1.120083E-6*TR**4 + 6.536332E-9*TR**5
C
      RHOR = RHOR + (0.824493 - 4.0899E-3*TR
     $       + 7.6438E-5*TR**2 - 8.2467E-7*TR**3
     $       + 5.3875E-9*TR**4) * SR
     $       + (-5.72466E-3 + 1.0227E-4*TR
     $       - 1.6546E-6*TR**2) * SR**1.5
     $       + 4.8314E-4 * SR**2
C
      CR=1449.1+.0821*P+4.55*TR-.045*TR**2
     $                              +1.34*(SR-35.)
      RHOR=RHOR + 1.E5*P/CR**2
     $     *(1.-2.*P/CR**2)
C
      RHOO(I,J,K)=RHOR*FSM(I,J)
    1 CONTINUE
C
      do k=1,33
         do j=1,jm
            do i=1,im
            si(i,j,k)=rhoo(i,j,k)
            end do
         end do
      end do
c
      RETURN
      END
      SUBROUTINE DEPTH(Z,ZZ,DZ,DZZ,DZR,KB,KBM1,H_MAX)
C >>>
      DIMENSION Z(KB),ZZ(KB),DZ(KB),DZZ(KB),DZR(KB)
C *************************************************
C    SUBROUTINE DEPTH.f.ori is the original one
C *************************************************
C
      Z(1)=0.E0
      DO K=2,8
        Z(K)=Z(K-1)-10.E0/H_MAX
      ENDDO
      Z(9)=-85.E0/H_MAX
      Z(10)=-100.E0/H_MAX
      Z(11)=-120.E0/H_MAX
      Z(12)=-150.E0/H_MAX
      Z(13)=-200.E0/H_MAX
      Z(14)=-300.E0/H_MAX
      Z(15)=-450.E0/H_MAX
      Z(16)=-650.E0/H_MAX
      Z(17)=-900.E0/H_MAX
      Z(18)=-1300.E0/H_MAX
      Z(19)=-1800.E0/H_MAX
      Z(20)=-2400.E0/H_MAX
      Z(21)=-3200.E0/H_MAX
      Z(22)=-4200.E0/H_MAX
      Z(KB)=-1.E0
C
      DO K=1,KBM1
        ZZ(K)=0.5*(Z(K)+Z(K+1))
      ENDDO
      ZZ(KB)=-1.E0
C
      DO K=1,KBM1
        DZ(K)=Z(K)-Z(K+1)
        DZZ(K)=ZZ(K)-ZZ(K+1)
      ENDDO
      DZ(KB)=DZ(KBM1)
      DZZ(KB)=DZZ(KBM1)
C
      DO K=1,KB
        DZR(K)=1./DZ(K)
      ENDDO
C
      WRITE(6,70)
   70 FORMAT(/2X,'K',7X,'Z',9X,'ZZ',9X,'DZ',9X,'DZZ',9x,'H'/)
cbt      OPEN(95,FILE='depth.dat')
      DO 90 K=1,KB
      WRITE(6,80) K,Z(K),ZZ(K),DZ(K),DZZ(K)
cbt      WRITE(95,*) ZZ(K)
   80 FORMAT(' ',I5,5F10.4)
   90 CONTINUE
      CLOSE(95)
C
      RETURN
      END
      FUNCTION EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
      real Rref18,Rref26,R,RMAX,WSMAX,WS18,WS26
      real b,expwnd
c 
       r1=0.5*(Rref18+Rref26)
       WS=0.5*(WS18+WS26)
       if(Rref18.le.0.) then
       r1=Rref26
       WS=WS26
       end if
       if(Rref26.le.0.) then
       r1=Rref18
       WS=WS18
       end if
c
         if(R.GE.RMAX) then
         b=(RMAX-r1)/log(WS/WSMAX)                     
         expwnd=WSMAX*exp((RMAX-R)/b)
         else
         expwnd=R*WSMAX/RMAX
         end if
c
      return
      end
      SUBROUTINE FINDPSI
         INCLUDE 'comblk.h'    
        DO 9004 J=1,JM
        DO 9004 I=1,IM
 9004   PSI(I,J)=0.E0
        DO 9005 J=2,JM
        DO 9005 I=1,IMM1
 9005   PSI(I+1,J)=PSI(I,J)-VAB(I,J)*.25E0*(D(I,J)+D(I,J-1))
     1           *(DX(I,J)+DX(I,J-1))
c---------- falk 04-23-02 remove error
c       DO 9006 J=2,JMM1
c       DO 9006 I=1,IM
        DO 9006 J=1,JMM1
        DO 9006 I=2,IM
 9006   PSI(I,J+1)=PSI(I,J)+.25E0*UAB(I,J)*(D(I,J)+D(I-1,J))
     1           *(DY(I,J)+DY(I-1,J))
      RETURN
      END
      SUBROUTINE HORINTERP(MSKVAL,IM,JM,IM1,JM1,X,Y,XI,YI,F,FI)
      REAL X1,X2,Y1,Y2,AREA,MSKVAL
      REAL F(IM1,JM1),FI(IM,JM)
      DIMENSION X(IM1),Y(JM1),XI(IM),YI(JM),ZT(4),ZTS(4)
      INTEGER II,JJ
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  FIRST INTERPOLATING HORIZONTALLY ALONG EACH LAYER
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF(XI(IM).GT.X(IM1).OR.XI(1).LT.X(1)) THEN
      print*,'X = ',x,'; XI = ',XI
      PRINT*,'THERE IS NO DATA FOR THIS VICINITY'
      STOP
      ELSE
      IF(YI(JM).GT.Y(JM1).OR.YI(1).LT.Y(1)) THEN
      print*,'Y = ',Y,'; YI = ',YI  
      PRINT*,'THERE IS NO DATA FOR THIS VICINITY'
      STOP
      END IF
      END IF

      DO J1=1,JM
         DO I1=1,IM
C
C******** LOOKING FOR CLOSEST POINTS IN X ****************
C
      II=1
      DO WHILE (ii.le.im1.and.(X(II)-XI(I1)).LT.0)
      II=II+1
      END DO
      if(ii.eq.1) II=II+1 
C
C******** LOOKING FOR CLOSEST POINTS IN Y ****************
C
      JJ=1
      DO WHILE (jj.le.jm1.and.(Y(JJ)-YI(J1)).LT.0)
      JJ=JJ+1
      END DO
      if(jj.eq.1) jj=jj+1 
C
C*****  CALCULATING PARAMETERS FOR FORMULA  **************
C
      X1=X(II-1)
      X2=X(II)
      Y1=Y(JJ-1)
      Y2=Y(JJ)
      ZT(1)=F(II-1,JJ-1)
      ZT(2)=F(II,JJ-1)
      ZT(3)=F(II,JJ)
      ZT(4)=F(II-1,JJ)
C
C******* CALCULATING FORMULA *****************************
C
      DO I=1,4
         IF(ZT(I).EQ.MSKVAL) THEN
      ZTS(I)=0.
         ELSE
      ZTS(I)=1.
         END IF
      END DO
C
      FI(I1,J1)=ZTS(1)*ZT(1)*(X2-XI(I1))*(Y2-YI(J1))
     1         +ZTS(2)*ZT(2)*(XI(I1)-X1)*(Y2-YI(J1))
     2         +ZTS(3)*ZT(3)*(XI(I1)-X1)*(YI(J1)-Y1)
     3         +ZTS(4)*ZT(4)*(X2-XI(I1))*(YI(J1)-Y1)
      AREA=ZTS(1)*ABS((X2-XI(I1))*(Y2-YI(J1)))
     1         +ZTS(2)*ABS((XI(I1)-X1)*(Y2-YI(J1)))
     2         +ZTS(3)*ABS((XI(I1)-X1)*(YI(J1)-Y1))
     3         +ZTS(4)*ABS((X2-XI(I1))*(YI(J1)-Y1))
C
c------------ falk 04-19-02 check more carefully
c     IF(AREA.EQ.0.) THEN
      IF(AREA.LT.(0.1*ABS(X2-X1)*ABS(Y2-Y1))) THEN
      FI(I1,J1)=MSKVAL
      ELSE
      FI(I1,J1)=FI(I1,J1)/AREA
      END IF
c
      END DO
         END DO
C
      RETURN
      END
      SUBROUTINE INTERP(F,F1)

      INCLUDE 'comblk.h'
      DIMENSION F(IM,JM,KB),F1(IM,JM,KB),DB(KB-1),DB1(KB-1)
      REAL tmp(kb-1),tmpi(kb-1)

      print*,'in subroutine interp'

      DO 15 J=1,JM
      DO 15 I=1,IM
c
      IF(FSM(I,J).EQ.0.0) THEN
        DO K=1,KB
        F1(I,J,K)=-99.9990
        END DO
      ELSE
        DO K=1,KB-1
        DB(K)=-H(I,J)*ZZ(K)
        DB1(K)=-ZZ(K)*H_MAX
        tmp(k)=f(i,j,k)
        END DO
        DB1(1)=0.0
      call verinterp(kb-1,kb-1,db,tmp,db1,tmpi)
        DO K=1,KB-1
            if(DB1(K).gt.DB(KB-1)) then
        f1(i,j,k)=-99.9990
            else
        f1(i,j,k)=tmpi(k)
            end if
        END DO
        f1(i,j,kb)=f1(i,j,kb-1)
      END IF
  15  continue
c
      RETURN
      END
c--------- falk 04-15-02 include NREAD in LOADGDEM
c--------- falk 04-16-02 include sstsource in LOADGDEM
      SUBROUTINE LOADGDEM(F,F1,NREAD,sstsource)
C**********************************************************
C  THIS SUBROUTINE TEMPR READS 3-D ARRAYS OF GDEM
C  TEMPERATURE AND SALINITY FROM FILES 'gdem.#' 
C**********************************************************

      INCLUDE 'comblk.h'
      PARAMETER(IM1=171,JM1=101,nl=33)
      real X(IM1),Y(JM1),XI(IM),YI(JM)
      real TEM(IM1,JM1),TEMI(IM,JM),f(IM,JM,nl),Zlev(nl),
     *     f1(IM,JM,nl),fin(IM1,JM1,nl),fin1(IM1,JM1,nl)
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      REAL XMIN,XMAX,YMIN,YMAX
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
c------ falk 01/08/00
      CHARACTER*50 FTSHARPb, FTSHARPa
c     FTSHARP='/migr/data/wd20af/forecast_system/ocean/sharpdata'
      FTSHARPb='bsharpdata'
      FTSHARPa='asharpdata'
c------ falk 01/08/00
c     OPEN(18,FILE=FTSHARPb,STATUS='unknown',form='formatted')
c     OPEN(28,FILE=FTSHARPa,STATUS='unknown',form='formatted')
c------
c
      DATA Zlev/0,10,20,30,50,75,100,125,150,200,250,300,400
     1  ,500,600,700,800,900,1000,1100,1200,1300,1400,1500
     2  ,1750,2000,2500,3000,3500,4000,4500,5000,5500/
c
      do j=1,jm
        do i=1,im
          xi(i)=longmin+(i-1)*(longmax-longmin)/(im-1)
          yi(j)=latmin+(j-1)*(latmax-latmin)/(jm-1)
        end do
      end do
c
c--------- falk 04-15-02 include NREAD in LOADGDEM
c
      print *,'    '
      print *,' LOADGDEM         NREAD=',NREAD
      print *,' LOADGDEM         sstsource=',sstsource
c
c------ falk 05-08-02 remove NREAD=-1; read initdata (sharpened T,S)
c
       rewind(13)
       read(13) f
       read(13) f1
       close(13)
c--------- falk 04-16-02 include sstsource in LOADGDEM
      if(sstsource.eq.1.or.sstsource.eq.2) then
          print *,' before MIXSSTZ'
        call MIXSSTZ(f,STM,H,FSM,Zlev,im,jm,nl,TBIAS)
          print *,' after MIXSSTZ'
      end if
c---------
c
c--------- falk 12/06/00 put in two lines along bnd mean values
      do i=1,im
       if(FSM(i,1).gt.0.5.and.FSM(i,2).gt.0.5) then
        do k=1,nl
         if(Zlev(k).le.H(i,1).and.Zlev(k).le.H(i,2)) then
          ac=(f(i,1,k)+f(i,2,k))*0.5
          bc=(f1(i,1,k)+f1(i,2,k))*0.5
          f(i,1,k)=ac
          f(i,2,k)=ac
          f1(i,1,k)=bc
          f1(i,2,k)=bc
         end if
        end do
       end if
       if(FSM(i,jm).gt.0.5.and.FSM(i,jm-1).gt.0.5) then
        do k=1,nl
         if(Zlev(k).le.H(i,jm).and.Zlev(k).le.H(i,jm-1)) then
          ac=(f(i,jm,k)+f(i,jm-1,k))*0.5
          bc=(f1(i,jm,k)+f1(i,jm-1,k))*0.5
          f(i,jm,k)=ac
          f(i,jm-1,k)=ac
          f1(i,jm,k)=bc
          f1(i,jm-1,k)=bc
         end if
        end do
       end if
      end do
c
      do j=1,jm
       if(FSM(1,j).gt.0.5.and.FSM(2,j).gt.0.5) then
        do k=1,nl
         if(Zlev(k).le.H(1,j).and.Zlev(k).le.H(2,j)) then
          ac=(f(1,j,k)+f(2,j,k))*0.5
          bc=(f1(1,j,k)+f1(2,j,k))*0.5
          f(1,j,k)=ac
          f(2,j,k)=ac
          f1(1,j,k)=bc
          f1(2,j,k)=bc
         end if
        end do
       end if
       if(FSM(im,j).gt.0.5.and.FSM(im-1,j).gt.0.5) then
        do k=1,nl
         if(Zlev(k).le.H(im,j).and.Zlev(k).le.H(im-1,j)) then
          ac=(f(im,j,k)+f(im-1,j,k))*0.5
          bc=(f1(im,j,k)+f1(im-1,j,k))*0.5
          f(im,j,k)=ac
          f(im-1,j,k)=ac
          f1(im,j,k)=bc
          f1(im-1,j,k)=bc
         end if
        end do
       end if
      end do
c
c------ falk 01-25-02 include write after MIXSSTZ and sharp
c
c     rewind(28)
c     write(28,101) IM,JM,nl
c     write(28,100) XI
c     write(28,100) YI
c     write(28,102) Zlev
c     write(28,102) H
c     write(28,103) FSM
c     write(28,103) f
c     write(28,103) f1
c     close(28)
c------
      return
      end
      SUBROUTINE OUTPUT(startdate)
C=========================================================
C save the output for graphics/analysis
C Written by Erxuan Fu, GSO,URI,11/3/94
C=========================================================
      INCLUDE 'comblk.h'
      REAL TMP1(IM,JM),TMP2(IM,JM),TMP3(IM,JM,KB),julday
      DIMENSION TB1(IM,JM,KB),U1(IM,JM,KB),V1(IM,JM,KB)
      INTEGER year
      CHARACTER FN*15, DOUT*8
      integer*4 startdate,date
C
      CALL DATE2DAY(year,julday,startdate)
CC
       julday=julday+time+1.e-5
CC
      call day2date(year,julday,date)
      WRITE(DOUT,'(I8.8)') date
c
c---------falk 12/12/00 include for GRADS printing
c
cbt      FN = 'GRADS.'//DOUT
cbt      OPEN(91,FILE=FN,STATUS='UNKNOWN',form='unformatted')
cbt      WRITE(91) T,S,RHO,U,V,UA,VA,ELB
cbt      WRITE(6,*) 'OUTPUT GRADS FILE '
cbt      WRITE(6,*) ' --> ',FN
cbt      CLOSE(91)
C
      FN = 'EL.'//DOUT
      OPEN(31,FILE=FN,STATUS='UNKNOWN',form='unformatted')
      WRITE(31) ELB
      CLOSE(31)
      WRITE(6,*) 'OUTPUT FILES:'
      WRITE(6,*) ' --> ',FN
C
      FN = 'TXY.'//DOUT
      OPEN(39,FILE=FN,STATUS='UNKNOWN',form='unformatted')
      WRITE(39) TAUX
      WRITE(39) TAUY
      CLOSE(39)
      WRITE(6,*) ' --> ',FN
C
c      FN = 'WXY.'//DOUT
c      OPEN(39,FILE=FN,STATUS='UNKNOWN',form='unformatted')
c      WRITE(39) WINDX
c      WRITE(39) WINDY
c      CLOSE(39)
c      WRITE(6,*) ' --> ',FN
C
      FN = 'UVA.'//DOUT
      do i=1,im-1
        do j=1,jm-1
          tmp1(i,j)=(uab(i,j)+uab(i+1,j))/2.0
          tmp2(i,j)=(vab(i,j)+vab(i,j+1))/2.0
          if(i.eq.1) tmp1(1,j)=uab(2,j)
          if(j.eq.1) tmp2(i,1)=vab(i,2)
        enddo
      enddo
      do i=1,im-1
        tmp1(i,jm)=(uab(i,jm)+uab(i+1,jm))/2.0
        tmp2(i,jm)=vab(i,jm-1)
        if(i.eq.1) tmp1(1,j)=uab(2,j)
      end do
      do j=2,jm-1
        tmp1(im,j)=uab(im-1,j)
        tmp2(im,j)=(vab(im,j)+vab(im,j+1))/2.0
        if(j.eq.1) tmp2(i,1)=vab(i,2)
      end do
      tmp1(im,jm)=uab(im,jm)
      tmp2(im,jm)=vab(im,jm-1)
      OPEN(32,FILE=FN,STATUS='UNKNOWN',form='unformatted')
      WRITE(32) TMP1
      WRITE(32) TMP2
      CLOSE(32)
      WRITE(6,*) ' --> ',FN
C
      CALL INTERP(TB,TB1)
c      CALL INTERP(tnudge,TB1)
      FN = 'T.'//DOUT
      OPEN(35,FILE=FN,STATUS='UNKNOWN',form='unformatted')
      WRITE(35) TB1
      CLOSE(35)
c
c      do k=1,kb
c        DO I=1,IM
c          DO J=1,JM
c            U1(I,J,K)=RHO(I,J,K)-RMEAN(I,J,K)
c          END DO
c        END DO
c      end do
c      CALL INTERP(U1,TB1)
c      FN = 'DN.'//DOUT
c      OPEN(35,FILE=FN,STATUS='UNKNOWN',form='unformatted')
c      WRITE(35) TB1
c      CLOSE(35)
cC
c      CALL INTERP(SB,TB1)
c      FN = 'S.'//DOUT
c      OPEN(35,FILE=FN,STATUS='UNKNOWN',form='unformatted')
c      WRITE(35) TB1
c      CLOSE(35)
C      
      FN = 'U.'//DOUT
c
C      Transfer into A-grid
c
      DO K=1,KB-1
      do j=1,jm
         do i=2,im-1
         if(dum(i,j).eq.1.and.dum(i+1,j).eq.1) 
     *   tmp3(i,j,k)=0.5*(ub(i,j,k)+ub(i+1,j,k))
         if(dum(i,j).eq.0.and.dum(i+1,j).eq.0) 
     *   tmp3(i,j,k)=0.0
         if(dum(i,j).eq.1.and.dum(i+1,j).eq.0)
     *   tmp3(i,j,k)=ub(i,j,k)
         if(dum(i,j).eq.0.and.dum(i+1,j).eq.1)
     *   tmp3(i,j,k)=ub(i+1,j,k)
         end do
      end do
      do j=1,jm
      tmp3(im,j,k)=ub(im,j,k)*dum(im,j)
      tmp3(1,j,k)=ub(2,j,k)*dum(2,j)
      end do
      ENDDO
c
      CALL INTERP(TMP3,U1)
      OPEN(37,FILE=FN,STATUS='UNKNOWN',form='unformatted')
c      WRITE(37,2030) TIME,IM,JM,KB
      WRITE(37) U1
      CLOSE(37)
C
      FN = 'V.'//DOUT
c
C      Transfer into A-grid
c
      DO K=1,KB-1
      do j=2,jm-1
         do i=1,im
         if(dvm(i,j).eq.1.and.dvm(i,j+1).eq.1) 
     *   tmp3(i,j,k)=0.5*(vb(i,j,k)+vb(i,j+1,k))
         if(dvm(i,j).eq.0.and.dvm(i,j+1).eq.0) 
     *   tmp3(i,j,k)=0.0
         if(dvm(i,j).eq.1.and.dvm(i,j+1).eq.0)
     *   tmp3(i,j,k)=vb(i,j,k)
         if(dvm(i,j).eq.0.and.dvm(i,j+1).eq.1)
     *   tmp3(i,j,k)=vb(i,j+1,k)
         end do
      end do
      do i=1,im
      tmp3(i,jm,k)=vb(i,jm,k)*dvm(i,jm)
      tmp3(i,1,k)=vb(i,2,k)*dvm(i,2)
      end do
      ENDDO
c
      CALL INTERP(TMP3,V1)
      OPEN(38,FILE=FN,STATUS='UNKNOWN',form='unformatted')
c      WRITE(38,2030) TIME,IM,JM,KB
      WRITE(38) V1
      CLOSE(38)
c
c      FN = 'US.'//DOUT
c      OPEN(38,FILE=FN,STATUS='UNKNOWN')
c      WRITE(38,2030) U
c      close(38)
c      FN = 'VS.'//DOUT
c      OPEN(38,FILE=FN,STATUS='UNKNOWN')
c      WRITE(38,2030) V
c      close(38)
c      FN = 'TS.'//DOUT
c      OPEN(38,FILE=FN,STATUS='UNKNOWN')
c      WRITE(38,2030) T
c      close(38)
C
C      FN = 'W.'//DOUT
C      OPEN(39,FILE=FN,STATUS='UNKNOWN',form='unformatted')
C      WRITE(39,2030) TIME,IM,JM,KB
C      WRITE(39,2031)(((W(I,J,K)*3600.,I=1,IM),J=1,JM),K=1,KB)
C      CLOSE(39)
C      WRITE(6,*) ' --> ',FN
C
c 2030 FORMAT(F8.3,3I8)
 2030 FORMAT(121F8.3)
 2031 FORMAT(10F8.3)
C
      RETURN
      END
      SUBROUTINE PROFQ(DT2)
C
       INCLUDE 'comblk.h'    
      REAL KN,KAPPA,UTIL       
      DIMENSION GH(IM,JM,KB),SM(IM,JM,KB),SH(IM,JM,KB)
      DIMENSION PROD(IM,JM,KB),KN(IM,JM,KB),BOYGR(IM,JM,KB)
      DIMENSION DH(IM,JM),CC(IM,JM,KB)
      EQUIVALENCE (A,SM),(C,SH),(PROD,KN),(TPS,DH)
      EQUIVALENCE (DTEF,CC,GH)
      DATA A1,B1,A2,B2,C1/0.92,16.6,0.74,10.1,0.08/
      DATA E1/1.8/,E2/1.33/,E3/1.0/
      DATA KAPPA/0.40/,SQ/0.20/,SEF/1./
      DATA GEE/9.806/
      DATA SMALL/1.E-8/
C
C
      DO 50 J=1,JM
      DO 50 I=1,IM
   50 DH(I,J)=H(I,J)+ETF(I,J)
      DO 100 K=2,KBM1
      DO 100 J=1,JM
      DO 100 I=1,IM
      A(I,J,K)=-DT2*(KQ(I,J,K+1)+KQ(I,J,K)+2.*UMOL)*.5
     1    /(DZZ(K-1)*DZ(K)*DH(I,J)*DH(I,J))
      C(I,J,K)=-DT2*(KQ(I,J,K-1)+KQ(I,J,K)+2.*UMOL)*.5
     1    /(DZZ(K-1)*DZ(K-1)*DH(I,J)*DH(I,J))
  100 CONTINUE
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C        DT2*(KQ*Q2')' - Q2*(2.*DT2*DTEF+1.) = -Q2B                    *
C                                                                      *
C***********************************************************************
C------  SURFACE AND BOTTOM B.C.S ------------
      CONST1=16.6**.6666667*SEF
      DO 90 J=1,JM-1
      DO 90 I=1,IM-1
      EE(I,J,1)=0.
      GG(I,J,1)=SQRT( (.5*(WUSURF(I,J)+WUSURF(I+1,J)))**2
     1                +(.5*(WVSURF(I,J)+WVSURF(I,J+1)))**2 )*CONST1
      UF(I,J,KB)=SQRT( (.5*(WUBOT(I,J)+WUBOT(I+1,J)))**2
     1                +(.5*(WVBOT(I,J)+WVBOT(I,J+1)))**2 )*CONST1
   90 CONTINUE
      DO J=1,JM-1
        EE(IM,J,1)=0.
        GG(IM,J,1)=SQRT( WUSURF(IM,J)**2+
     1  (.5*(WVSURF(IM,J)+WVSURF(IM,J+1)))**2 )*CONST1
        UF(IM,J,KB)=SQRT( WUBOT(IM,J)**2+
     1  (.5*(WVBOT(IM,J)+WVBOT(IM,J+1)))**2 )*CONST1
      END DO
      DO I=1,IM-1
        EE(I,JM,1)=0.
        GG(I,JM,1)=SQRT( (.5*(WUSURF(I,JM)+WUSURF(I+1,JM)))**2+
     1  WVSURF(I,JM)**2 )*CONST1
        UF(I,JM,KB)=SQRT( (.5*(WUBOT(I,JM)+WUBOT(I+1,JM)))**2+
     1  WVBOT(I,JM)**2 )*CONST1
      END DO
      EE(IM,JM,1)=0.
      GG(IM,JM,1)=SQRT( WUSURF(IM,JM)**2
     1                +WVSURF(IM,JM)**2 )*CONST1
      UF(IM,JM,KB)=SQRT( WUBOT(IM,JM)**2
     1                +WVBOT(IM,JM)**2 )*CONST1
C---------------------------------------------------------------------
C   RESTORE 10 DEG, 35 PPT AND 1.025 BIAS IN T, S AND RHO
C---------------------------------------------------------------------
      DO 101 K=1,KBM1
      DO 101 J=1,JM
      DO 101 I=1,IM
      TP=T(I,J,K)+TBIAS
      SP=S(I,J,K)+SBIAS
C      Calculate pressure in units of decibars
      P=-GEE*1.025*ZZ(K)*DH(I,J)*.1 
      CC(I,J,K)=1449.1+.00821*P+4.55*TP -.045*TP**2
     $                              +1.34*(SP-35.)
      CC(I,J,K)=CC(I,J,K)/SQRT((1.-.01642*P/CC(I,J,K))
     $       *(1.-0.40*P/CC(I,J,K)**2))
  101 CONTINUE
      DO 102 K=2,KBM1
      DO 102 J=1,JM
      DO 102 I=1,IM
      Q2B(I,J,K)=ABS(Q2B(I,J,K))
      Q2LB(I,J,K)=ABS(Q2LB(I,J,K))
      BOYGR(I,J,K)=GEE*((RHO(I,J,K-1)-RHO(I,J,K))/(DZZ(K-1)*DH(I,J)))      
     1 +GEE**2*2.*1.025/(CC(I,J,K-1)**2+CC(I,J,K)**2)
  102 CONTINUE
C------ CALC. T.K.E. PRODUCTION -----------------------------------
      DO 120 K=2,KBM1
      DO 120 J=1,JMM1
      DO 120 I=1,IMM1
      PROD(I,J,K)=KM(I,J,K)*.25*SEF
     1       *( (U(I,J,K)-U(I,J,K-1)+U(I+1,J,K)-U(I+1,J,K-1))**2
     2         +(V(I,J,K)-V(I,J,K-1)+V(I,J+1,K)-V(I,J+1,K-1))**2 )
     3              /(DZZ(K-1)*DH(I,J))**2
  120 PROD(I,J,K)=PROD(I,J,K)+KH(I,J,K)*BOYGR(I,J,K)
      do k=2,kbm1
        do j=1,jmm1
          prod(im,j,k)=km(im,j,k)*sef*
     1    ( (u(im,j,k)-u(im,j,k-1))**2+ 0.25*
     2    (v(im,j,k)-v(im,j,k-1)+v(im,j+1,k)-v(im,j+1,k-1))**2 )/
     3    (dzz(k-1)*dh(i,j))**2
        end do
        do i=1,imm1
          prod(i,jm,k)=km(i,jm,k)*sef*
     1    (.25*(u(i,jm,k)-u(i,jm,k-1)+u(i+1,jm,k)-u(i+1,jm,k-1))**2
     2    +(v(i,jm,k)-v(i,jm,k-1))**2 )/(dzz(k-1)*dh(i,jm))**2
        end do
        prod(im,jm,k)=km(im,jm,k)*sef*
     1  ( (u(im,jm,k)-u(im,jm,k-1))**2+
     2  (v(im,jm,k)-v(im,jm,k-1))**2 )/(dzz(k-1)*dh(im,jm))**2
      end do
      DO 110 K=2,KBM1
      DO 110 J=1,JM
      DO 110 I=1,IM
  110 DTEF(I,J,K)=Q2B(I,J,K)*SQRT(Q2B(I,J,K))/(B1*Q2LB(I,J,K)+SMALL)
      DO 140 K=2,KBM1
      DO 140 J=1,JM
      DO 140 I=1,IM
      GG(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1))
     1    -(2.*DT2*DTEF(I,J,K)+1.) )
      EE(I,J,K)=A(I,J,K)*GG(I,J,K)
      GG(I,J,K)=(-2.*DT2*PROD(I,J,K)
     1  +C(I,J,K)*GG(I,J,K-1)-UF(I,J,K))*GG(I,J,K)
  140 CONTINUE
      DO 150 K=1,KBM1
      KI=KB-K
      DO 150 J=1,JM
      DO 150 I=1,IM
      UF(I,J,KI)=EE(I,J,KI)*UF(I,J,KI+1)+GG(I,J,KI)
  150 CONTINUE
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C        DT2(KQ*Q2L')' - Q2L*(DT2*DTEF+1.) = -Q2LB                     *
C                                                                      *
C***********************************************************************
      DO 155 J=1,JM
      DO 155 I=1,IM
      EE(I,J,1)=0.0
      GG(I,J,1)=0.0
      VF(I,J,KB)=0.     
  155 CONTINUE   
      DO 160 K=2,KBM1
      DO 160 J=1,JM
      DO 160 I=1,IM
      DTEF(I,J,K) =DTEF(I,J,K)*(1.+E2*((1./ABS(Z(K)-Z(1))+
     1    1./ABS(Z(K)-Z(KB))) *L(I,J,K)/(DH(I,J)*KAPPA))**2)
      GG(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1))
     1    -(DT2*DTEF(I,J,K)+1.))
      EE(I,J,K)=A(I,J,K)*GG(I,J,K)
      GG(I,J,K)=(DT2*(-PROD(I,J,K)
     1   *L(I,J,K)*E1)+C(I,J,K)*GG(I,J,K-1)-VF(I,J,K))*GG(I,J,K)
  160 CONTINUE
      DO 170 K=1,KBM1
      KI=KB-K
      DO 170 J=1,JM
      DO 170 I=1,IM
      VF(I,J,KI)=EE(I,J,KI)*VF(I,J,KI+1)+GG(I,J,KI)
  170 CONTINUE
      DO 180 K=2,KBM1
      DO 180 J=1,JM
      DO 180 I=1,IM
      IF(UF(I,J,K).LT.SMALL.OR.VF(I,J,K).LT.SMALL) THEN
        UF(I,J,K)=SMALL
        VF(I,J,K)=SMALL
      ENDIF
  180 CONTINUE
C***********************************************************************
C                                                                      *
C               THE FOLLOWING SECTION SOLVES FOR KM AND KH             *
C                                                                      *
C***********************************************************************
      COEF1=A2*(1.-6.*A1/B1)
      COEF2=3.*A2*B2+18.*A1*A2
      COEF3=A1*(1.-3.*C1-6.*A1/B1)
      COEF4=18.*A1*A1+9.*A1*A2
      COEF5=9.*A1*A2
C NOTE THAT SM,SH LIMIT TO INFINITY WHEN GH APPROACHES 0.0288
      DO 220 K=2,KBM1
      DO 220 J=1,JM
      DO 220 I=1,IM
      L(I,J,K)=VF(I,J,K)/UF(I,J,K)
 220  GH(I,J,K)=L(I,J,K)**2/UF(I,J,K)*BOYGR(I,J,K)              
      DO 225 J=1,JM
      DO 225 I=1,IM
      L(I,J,1)=0.   
      L(I,J,KB)=0.   
      GH(I,J,1)=0.  
 225  GH(I,J,KB)=0.  
      DO 230 K=1,KB
      DO 230 J=1,JM
      DO 230 I=1,IM
      GH(I,J,K)=MIN(GH(I,J,K),.028)
      SH(I,J,K)=COEF1/(1.-COEF2*GH(I,J,K))
      SM(I,J,K)=COEF3+SH(I,J,K)*COEF4*GH(I,J,K)
      SM(I,J,K)=SM(I,J,K)/(1.-COEF5*GH(I,J,K))
  230 CONTINUE
C
      DO 280 K=1,KB
      DO 280 J=1,JM
      DO 280 I=1,IM
      KN(I,J,K)=L(I,J,K)*SQRT(ABS(Q2(I,J,K)))
      KQ(I,J,K)=(KN(I,J,K)*.41*SM(I,J,K)+KQ(I,J,K))*.5
      KM(I,J,K)=(KN(I,J,K)*SM(I,J,K)+KM(I,J,K))*.5
      KH(I,J,K)=(KN(I,J,K)*SH(I,J,K)+KH(I,J,K))*.5
  280 CONTINUE
Cfr**********************************************************
C    LIMITATION ON VERTICAL DIFFUSION COEFFICIENT
Cfr**********************************************************
C      DO 290 J=1,JM
C      DO 290 I=1,IM
C      DO 290 K1=1,KB-2
C      K=KB-K1
C      UTIL=0.25*DZ(K)*(H(I,J)+ETF(I,J))**2/DTI
C      IF(UTIL.LE.(KH(I,J,K+1)+UMOL)/DZZ(K)) 
C     1    KH(I,J,K+1)=0.1*UTIL*DZZ(K)-UMOL
C      IF((KH(I,J,K)+UMOL).GT.
C     1    (UTIL-(UMOL+KH(I,J,K+1))/DZZ(K))*DZZ(K-1)) KH(I,J,K)
C     2    =(UTIL-(UMOL+KH(I,J,K+1))/DZZ(K))*DZZ(K-1)-UMOL
C      IF(KH(I,J,K+1).LE.0) KH(I,J,K+1)=1.E-12
C      IF(KH(I,J,K).LE.0) KH(I,J,K)=1.E-12
C      IF(KH(I,J,K-1).LE.0) KH(I,J,K-1)=1.E-12
C  290 CONTINUE
Cfr**********************************************************

      RETURN   
      END
      SUBROUTINE PROFT(F,WFSURF,SWRAD,FSURF,NBC2D,DT2)
C
	INCLUDE 'comblk.h'
      DIMENSION F(IM,JM,KB),WFSURF(IM,JM),FSURF(IM,JM),DH(IM,JM)
      DIMENSION SWRAD(IM,JM),RAD(IM,JM,KB),TR(5),EXTC(5),NBC2D(im,jm)
      EQUIVALENCE (TPS,DH)
C
      NTP=2
C       NTP         =     1      2       3       4       5
C   JERLOV TYPE     =     I      IA      IB      II      III
      DATA TR/          .32  ,  .31  ,  .29  ,  .26  ,  .24   /
      DATA EXTC/       .037  , .042  , .056  , .073  , .127   /
C
      UMOLPR=UMOL
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C         DTI2*(KH*F')'-F=-FB                                          *
C                                                                      *
C***********************************************************************
      DO 10 J=1,JM
      DO 10 I=1,IM
      DH(I,J)=H(I,J)+ETF(I,J)
   10 CONTINUE
      DO 20 K=2,KBM1
      DO 20 J=1,JM
      DO 20 I=1,IM
      A(I,J,K-1)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K-1)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
      C(I,J,K)=-DT2*(KH(I,J,K)+UMOLPR)/(DZ(K)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
   20 CONTINUE
C-----------------------------------------------------------------------
C   NBC=1: SURF. B.C. IS WFSURF. NO SW RADIATIVE PENETRATION.
C   NBC=2; SURF. B.C. IS WFSURF+SWRAD*(1.-TR). TR*SWRAD PENETRATES WATER COLUMN
C   NBC=3; SURF. B.C. IS TSURF. NO SW RADIATIVE PENETRATION.
C   NBC=4; SURF. B.C. IS TSURF. TR*SWRAD PENETRATES WATER COLUMN
C
C NOTE THAT WTSURF (=WFSURF) AND SWRAD ARE NEG. VALUES WHEN WATER COLUMN IS
C            WARMING.
C-----------------------------------------------------------------------
C------------------------------------------------------------------
C     Penetrative Radiation Calculation. At the bottom any 
C     unattenuated radiation is deposited in the bottom layer.
C------------------------------------------------------------------
        DO 512 K=1,KB
        DO 512 J=1,JM
        DO 512 I=1,IM
  512   RAD(I,J,K)=0.
      DO K=1,KBM1
        DO J=1,JM
          DO I=1,IM
            IF(NBC2D(i,j).EQ.2.OR.NBC2D(i,j).EQ.4.) THEN   
              RAD(I,J,K)=SWRAD(I,J)*TR(NTP)*
     1        EXP(EXTC(NTP)*Z(K)*DH(I,J))
            ENDIF
          END DO
        END DO
      END DO
      DO J=1,JM
        DO I=1,IM
          GO TO (50,51,52,52), NBC2D(i,j)
   50     CONTINUE
          EE(I,J,1)=A(I,J,1)/(A(I,J,1)-1.0)
          GG(I,J,1)=-DT2*WFSURF(I,J)/(-DZ(1)*DH(I,J))-F(I,J,1)
          GG(I,J,1)=GG(I,J,1)/(A(I,J,1)-1.0)
          GO TO 53
C
   51     CONTINUE
          EE(I,J,1)=A(I,J,1)/(A(I,J,1)-1.0)
          GG(I,J,1)=DT2*(WFSURF(I,J)+SWRAD(I,J)*(1.0-TR(NTP))
     2    +RAD(I,J,1)-RAD(I,J,2))
     1    /(DZ(1)*DH(I,J))-F(I,J,1)
          GG(I,J,1)=GG(I,J,1)/(A(I,J,1)-1.0)
          GO TO 53
C
   52     CONTINUE
          EE(I,J,1)=0.
          GG(I,J,1)=FSURF(I,J)
c
   53     CONTINUE
        END DO
      END DO
C
C----------------------------------------------------------------------
      DO 101 K=2,KBM2
      DO 101 J=1,JM
      DO 101 I=1,IM
      GG(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1))-1.)
      EE(I,J,K)=A(I,J,K)*GG(I,J,K)
      GG(I,J,K)=(C(I,J,K)*GG(I,J,K-1)-F(I,J,K)
     1     +DT2*(RAD(I,J,K)-RAD(I,J,K+1))/(DH(I,J)*DZ(K)))*GG(I,J,K)
  101 CONTINUE
C-----  BOTTOM ADIABATIC B.C. ------------------------------------------
      DO 102 J=1,JM
      DO 102 I=1,IM
  102 F(I,J,KBM1)=((C(I,J,KBM1)*GG(I,J,KBM2)-F(I,J,KBM1)
     1         +DT2*(RAD(I,J,KBM1)-RAD(I,J,KB))/(DH(I,J)*DZ(K)))
     2          /(C(I,J,KBM1)*(1.-EE(I,J,KBM2))-1.))
C----------------------------------------------------------------------
      DO 105 K=2,KBM1
      KI=KB-K
      DO 105 J=1,JM
      DO 105 I=1,IM
      F(I,J,KI)=(EE(I,J,KI)*F(I,J,KI+1)+GG(I,J,KI))
  105 CONTINUE
C
      RETURN
      END
      SUBROUTINE PROFU(DT2)
       INCLUDE 'comblk.h'    
      DIMENSION DH(IM,JM)
      DATA UMOLB/1.36E-6 /
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C         DT2*(KM*U')'-U=-UB                                           *
C                                                                      *
C***********************************************************************
C
      DO 84 J=1,JM
      DO 84 I=1,IM
84    DH(I,J)=1.0
 
      DO 85 J=2,JM
      DO 85 I=2,IM
   85 DH(I,J)=.5E0*(H(I,J)+ETF(I,J)+H(I-1,J)+ETF(I-1,J))
      DO 90 K=1,KB
      DO 90 J=2,JM
      DO 90 I=2,IM
   90 C(I,J,K)=(KM(I,J,K)+KM(I-1,J,K))*.5E0
      DO 100 K=2,KBM1
      DO 100 J=1,JM
      DO 100 I=1,IM
      A(I,J,K-1)=-DT2*(C(I,J,K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
      C(I,J,K)=-DT2*(C(I,J,K)+UMOL  )/(DZ(K)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
  100 CONTINUE
C
      DO 1011 J=1,JM
      DO 1011 I=1,IM
      EE(I,J,1)=A(I,J,1)/(A(I,J,1)-1.E0)
      GG(I,J,1)=(-DT2*WUSURF(I,J)/(-DZ(1)*DH(I,J))-UF(I,J,1))
     1   /(A(I,J,1)-1.E0)
1011  CONTINUE
      DO 101 K=2,KBM2
      DO 101 J=1,JM
      DO 101 I=1,IM
      GG(I,J,K)=1.E0/(A(I,J,K)+C(I,J,K)*(1.-EE(I,J,K-1))-1.)
      EE(I,J,K)=A(I,J,K)*GG(I,J,K)
      GG(I,J,K)=(C(I,J,K)*GG(I,J,K-1)-UF(I,J,K))*GG(I,J,K)
  101 CONTINUE
      DO 102 J=2,JMm1
      DO 102 I=2,IMm1
      TPS(I,J)=0.5E0*(CBC(I,J)+CBC(I-1,J))
     1     *SQRT(UB(I,J,KBM1)**2+(.25E0*(VB(I,J,KBM1)
     2     +VB(I,J+1,KBM1)+VB(I-1,J,KBM1)+VB(I-1,J+1,KBM1)))**2)
      UF(I,J,KBM1)=(C(I,J,KBM1)*GG(I,J,KBM2)-UF(I,J,KBM1))/(TPS(I,J)
     1 *DT2/(-DZ(KBM1)*DH(I,J))-1.E0-(EE(I,J,KBM2)-1.E0)*C(I,J,KBM1))
  102 UF(I,J,KBM1)=UF(I,J,KBM1)*DUM(I,J)
      DO 103 K=2,KBM1
      KI=KB-K
      DO 103 J=2,JMM1
      DO 103 I=2,IMM1
      UF(I,J,KI)=(EE(I,J,KI)*UF(I,J,KI+1)+GG(I,J,KI))*DUM(I,J)
  103 CONTINUE
C
      DO 104 J=2,JMM1
      DO 104 I=2,IMM1
104   WUBOT(I,J)=-TPS(I,J)*UF(I,J,KBM1)
      RETURN
      END
      SUBROUTINE PROFV(DT2)
       INCLUDE 'comblk.h'    
      DIMENSION DH(IM,JM)
      DATA UMOLB/1.36E-6/
C***********************************************************************
C                                                                      *
C        THE FOLLOWING SECTION SOLVES THE EQUATION                     *
C         DT2*(KM*U')'-U=-UB                                           *
C                                                                      *
C***********************************************************************
C
      DO 84 J=1,JM
      DO 84 I=1,IM
84    DH(I,J)=1.
      DO 85 J=2,JM
      DO 85 I=2,IM
   85 DH(I,J)=.5E0*(H(I,J)+ETF(I,J)+H(I,J-1)+ETF(I,J-1))
      DO 90 K=1,KB
      DO 90 J=2,JM
      DO 90 I=2,IM
   90 C(I,J,K)=(KM(I,J,K)+KM(I,J-1,K))*.5E0
      DO 100 K=2,KBM1
C
      DO 100 J=1,JM
      DO 100 I=1,IM
      A(I,J,K-1)=-DT2*(C(I,J,K)+UMOL  )/(DZ(K-1)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
      C(I,J,K)=-DT2*(C(I,J,K)+UMOL  )/(DZ(K)*DZZ(K-1)*DH(I,J)
     1     *DH(I,J))
  100 CONTINUE
C
      DO 1001 J=1,JM
      DO 1001 I=1,IM
      EE(I,J,1)=A(I,J,1)/(A(I,J,1)-1.E0)
      GG(I,J,1)=(-DT2*WVSURF(I,J)/(-DZ(1)*DH(I,J))-VF(I,J,1))
     1   /(A(I,J,1)-1.E0)
1001  CONTINUE
      DO 101 K=2,KBM2
      DO 101 J=1,JM
      DO 101 I=1,IM
      GG(I,J,K)=1.E0/(A(I,J,K)+C(I,J,K)*(1.E0-EE(I,J,K-1))-1.E0)
      EE(I,J,K)=A(I,J,K)*GG(I,J,K)
      GG(I,J,K)=(C(I,J,K)*GG(I,J,K-1)-VF(I,J,K))*GG(I,J,K)
  101 CONTINUE
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
      TPS(I,J)=0.5E0*(CBC(I,J)+CBC(I,J-1))
     1     *SQRT((.25E0*(UB(I,J,KBM1)+UB(I+1,J,KBM1)
     2     +UB(I,J-1,KBM1)+UB(I+1,J-1,KBM1)))**2+VB(I,J,KBM1)**2)
      VF(I,J,KBM1)=(C(I,J,KBM1)*GG(I,J,KBM2)-VF(I,J,KBM1))/(TPS(I,J)
     1  *DT2/(-DZ(KBM1)*DH(I,J))-1.E0-(EE(I,J,KBM2)-1.E0)*C(I,J,KBM1))
  102 VF(I,J,KBM1)=VF(I,J,KBM1)*DVM(I,J)
      DO 103 K=2,KBM1
      KI=KB-K
      DO 103 J=2,JMM1
      DO 103 I=2,IMM1
      VF(I,J,KI)=(EE(I,J,KI)*VF(I,J,KI+1)+GG(I,J,KI))*DVM(I,J)
  103 CONTINUE
C
      DO 104 J=2,JMM1
      DO 104 I=2,IMM1
104   WVBOT(I,J)=-TPS(I,J)*VF(I,J,KBM1)
      RETURN
      END
c
c-------------------
c
c--------- falk 07-06-01 use SERFTEMPR which reads GFS sst

      SUBROUTINE SERFTEMPR(IM1,JM1)
C Changes are made on 11/01/00 to read GFS global sst
C Biju Thomas

      INCLUDE 'comblk.h'
c     PARAMETER(IM1=512,JM1=256)
      DIMENSION XI(IM),YI(JM),X(IM1),Y(JM1)
      DIMENSION XGFS(IM1),YGFS(JM1)
      DIMENSION TEM(IM1,JM1),TB1(IM,JM),ZT(4),ZTS(4)
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      REAL XMIN,YMIN,MSK(IM1,JM1)
      REAL CMP
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
c
CCC      XMIN=255.5-360
CCC      read(45,*)XMIN
CCC      close(45)
	XMIN=-180.0
c
C***********************************************************
C   SUBROUTINE SURFTEMPR ASSIGN THE ARRAY OF SEA-SURFACE 
C   TEMPERATURE. THE DATA ARE STORED IN THE FILE 'surftem.dat'
C***********************************************************
c
CC       OPEN(89,FILE=
CC     *'sst.gfs.dat',status='old',FORM='UNFORMATTED')
CC      OPEN(99,FILE=
CC     *'mask.gfs.dat',status='old',FORM='UNFORMATTED')
CC      OPEN(98,FILE=
CC     *'lonlat_gfs_t170',status='old',FORM='FORMATTED')
C
CCC      read(8,1900) cmp
c
c-------------- falk 03-26-02 
      print *,'SERFTEMPR: dim gfs SST data IM1,JM1=',IM1,JM1
      cmp=0.
      READ(21) TEM
      READ(22) MSK
      rewind(21)
      rewind(22)
C      CLOSE(21)
C      CLOSE(22)
C      CLOSE(8)
      READ(23,124)IM3,JM3
      READ(23,123)XGFS,YGFS
 123  FORMAT(1X,10E10.4)
 124  FORMAT(1X,2I5)
      rewind(23)
 100  FORMAT(x,10F7.2)
 1900 FORMAT(x,F7.2)
  101 FORMAT(75e16.5)
c
c      OPEN(7,FILE='serftem.dat',FORM='FORMATTED',recl=1600)
c      WRITE(7,101) TEM
c      CLOSE(7)
c      print*,'Finished writinf serftem.dat'
cc
c      OPEN(7,FILE='serfmsk.dat',FORM='FORMATTED',recl=1600)
c      WRITE(7,101) MSK
c      CLOSE(7)
c      print*,'Finished writinf serfmsk.dat'
c
C*********************************************************
C  EXTRAPOLATING ORIGINAL DATA TO FILL MISSING DATA POINTS
C*********************************************************
      DO N1=1,10
      DO 9 J=2,JM1-1
      DO 9 I=2,IM1-1
         IF(MSK(I,J).EQ.1) THEN
      TEM(I,J)=0
      CMP=0
      ZT(1)=TEM(I+1,J)
      ZTS(1)=1-MSK(I+1,J)
      ZT(2)=TEM(I-1,J)
      ZTS(2)=1-MSK(I-1,J)
      ZT(3)=TEM(I,J+1)
      ZTS(3)=1-MSK(I,J+1)
      ZT(4)=TEM(I,J-1)
      ZTS(4)=1-MSK(I,J-1)
        DO N=1,4
        TEM(I,J)=TEM(I,J)+ZTS(N)*ZT(N)
        CMP=CMP+ZTS(N)
        END DO
      IF(CMP.EQ.0) THEN
        TEM(I,J)=0.0
      ELSE
        TEM(I,J)=TEM(I,J)/CMP
        MSK(I,J)=0.0
      END IF
         END IF
   9  CONTINUE
      END DO
c
C----------- INTERPOLATING HORIZONTALLY ALONG EACH LAYER ---------
c
      DO J=1,JM
         DO I=1,IM
         XI(I)=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
         YI(J)=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
         END DO
      END DO
c
      DO i=1,IM1
c        x(i)=xmin+(i-1)*0.703125
         x(i)=XGFS(I)
      END DO
         DO J=1,JM1
         y(j)=YGFS(J)
         END DO

c
      CALL HORINTERP(0.0,IM,JM,IM1,JM1,x,y,xi,yi,TEM,TB1)
c
      DO J=1,JM
         DO I=1,IM
         STM(I,J)=TB1(I,J)-273.19-TBIAS
         STM(I,J)=STM(I,J)*FSM(I,J)
         end do
      end do
c
      RETURN
      END
c
c----------
c
      SUBROUTINE SMOOTHING(F,NI,NJ,SM)
C*****************************************
C Smoothing(SM>0) and de-smoothing(SM<0).
C This code was developed by Dr. Ginis,
C rewrited by Erxuan Fu. 12/1/94,GSO,URI.
C*****************************************
      INCLUDE 'comblk.h'
      DIMENSION F(NI,NJ),F1(500,500)
C
C**********************************************
C X-SMOOTHING
C**********************************************
      DO 112 J=2,NJ-1
      DO 112 I=2,NI-1
  112 F1(I,J)=F(I,J)
      DO 111 J=2,NJ-1
      DO 111 I=2,NI-1
      IF(FSM(I,J).NE.0.) THEN
      IF(FSM(I-1,J).EQ.0) THEN
      IF(FSM(I+1,J).EQ.0) THEN
      F1(I,J)=F(I,J)
      ELSE
      F1(I,J)=F(I,J)+SM*(F(I+1,J)+F(I,J)-2.*F(I,J))
      END IF
      ELSE
      IF(FSM(I+1,J).EQ.0) THEN
      F1(I,J)=F(I,J)+SM*(F(I,J)+F(I-1,J)-2.*F(I,J))
      ELSE
      F1(I,J)=F(I,J)+SM*(F(I+1,J)+F(I-1,J)-2.*F(I,J))
      END IF
      END IF
      ELSE
      F1(I,J)=F(I,J)
      END IF
  111 CONTINUE
      DO 97 I=2,NI-1
      DO 97 J=2,NJ-1
   97 F(I,J)=F1(I,J)
C
C**********************************************
C Y-SMOOTHING
C**********************************************
      DO 114 I=2,NI-1
      DO 114 J=2,NJ-1
      IF(FSM(I,J).NE.0.) THEN
      IF(FSM(I,J-1).EQ.0) THEN
      IF(FSM(I,J+1).EQ.0) THEN
      F1(I,J)=F(I,J)
      ELSE
      F1(I,J)=F(I,J)+SM*(F(I,J+1)+F(I,J)-2.*F(I,J))
      END IF
      ELSE
      IF(FSM(I,J+1).EQ.0) THEN
      F1(I,J)=F(I,J)+SM*(F(I,J)+F(I,J-1)-2.*F(I,J))
      ELSE
      F1(I,J)=F(I,J)+SM*(F(I,J+1)+F(I,J-1)-2.*F(I,J))
      END IF
      END IF
      ELSE
      F1(I,J)=F(I,J)
      END IF
  114 CONTINUE
C 
      DO 113 J=2,NJ-1
      DO 113 I=2,NI-1
  113 F(I,J)=F1(I,J)
C
      RETURN
      END          
c-------- falk 04-15-02 include NREAD in TEMPR
c-------- falk 04-16-02 include sstsource in TEMPR
      SUBROUTINE TEMPR(NREAD,sstsource)

      INCLUDE 'comblk.h'
      CHARACTER*20 FILENAME
      DIMENSION ZI(KB),ZI1(KB)
      DIMENSION TMPT(33),TMPTI(KB),TMPS(33),TMPSI(KB)
      DIMENSION TEMI(IM,JM),TB1(IM,JM,33),SB1(IM,JM,33)
      REAL CMP,L1,L2,LATMIN,LATMAX,LONGMIN,LONGMAX
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
      common/mean/TBM(33),SBM(33),RHM(33),DB(33)
C
C**********************************************************
C  SUBROUTINE TEMPR CREATES ARRAY OF INITIAL TEMPERATURE 
C  INTERPOLATING LEVITUS DATA.  IT ALSO CALCULATES 
C  MEAN TEMPERATURE AND DENSITY FIELDS.
C**********************************************************
C
      DATA DB/0,10,20,30,50,75,100,125,150,200,250,300,400
     1  ,500,600,700,800,900,1000,1100,1200,1300,1400,1500
     2  ,1750,2000,2500,3000,3500,4000,4500,5000,5500/
  102 FORMAT(251F8.3)
c
      PRINT*,'IN SUBROUTINE TEMPR'
cC
c      CALL LOADLEV(1,TB1)
cc
c      CALL LOADLEV(2,SB1)
cc      read(61) TB1,SB1
c
c-------- falk 04-15-02 include NREAD in TEMPR and LOADGDEM
      CALL LOADGDEM(TB1,SB1,NREAD,sstsource)
c      do k=1,33
c         do j=1,jm
c            do i=1,im
c            sb1(i,j,k)=35.
c            end do
c         end do
c      end do      
c
C*********************************************************
C   INCORPORATING SEA-SURFACE TEMPERATURE
C*********************************************************
c      DO 21 J=1,JM
c      DO 21 I=1,IM
c      TB1(I,J,1)=STM(I,J)+TBIAS
c      DO K=2,33
c      IF(TB1(I,J,K).GT.TB1(I,J,1)) TB1(I,J,K)=TB1(I,J,1)
cC      IF(H(I,J).LE.50.AND.H(I,J).GT.1) TB1(I,J,K)=TB1(I,J,1)
c      END DO
c  21  CONTINUE

c      OPEN(7,FILE='controllevtem.dat',FORM='FORMATTED')
c      WRITE(7,102) TB1
c      CLOSE(7)

c      print*,'Finished writing temperature test file'
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    CALCULATING HORIZONTAL MEAN TEMPERATURE AND DENSITY
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DO 88 K=1,33
      DO 88 J=1,JM
      DO 88 I=1,IM
      SB1(I,J,K)=SB1(I,J,K)-SBIAS
  88  TB1(I,J,K)=TB1(I,J,K)-TBIAS

      DO K=1,33
      CMP=0
      CMP1=0
      TBM(K)=0
      SBM(K)=0
      DO 15 J=1,JM
      DO 15 I=1,IM
      IF(TB1(I,J,K).LE.-99) THEN
      CMP=CMP-1
      ELSE
      TBM(K)=TBM(K)+TB1(I,J,K)*FSM(I,J)
      END IF
      IF(SB1(I,J,K).LE.-99) THEN
      CMP1=CMP1-1
      ELSE
      SBM(K)=SBM(K)+SB1(I,J,K)*FSM(I,J)
      END IF
  15  CONTINUE
      DO J=1,JM
      DO I=1,IM
      CMP=CMP+FSM(I,J)
      CMP1=CMP1+FSM(I,J)
      END DO
      END DO
      TBM(K)=TBM(K)/CMP
      SBM(K)=SBM(K)/CMP1
      END DO
c
      CALL CLIMAT2MODEL(SB1,TB1)
c
      RETURN
      END

     
      SUBROUTINE TPROF(H0,H1,T0,T1)
C=========================================================
C Define temperature profile in a form:
C       T(h)=A0/(h-B0).
C Parameters:H0,H1>0 in ocean here! Variable:DD < 0 in ocean!
C H0=20m,H1=760m,T0=29,T1=5.8 
C=========================================================
C                     T1         T0
C       -----------|----------|---
C               |             *
C          H0   -             *
C               |           *
C               |         *
C               |       *
C               |     *
C               |    *
C               |   *
C          H1   -  *
C===========================================================
      INCLUDE 'comblk.h'
      REAL H0,H1,T0,T1,A0,B0,DD
C
      B0=(-H0*T0+H1*T1)/(T0-T1)
      A0=(-H1-B0)*T1
      DO 100 I=1,IM
      DO 100 J=1,JM
      DO 200 K=1,KB
      DD=ZZ(K)*H(I,J)
      IF(ABS(DD).LE.H0) THEN
      TB(I,J,K)=T0
      ELSE
      IF(ABS(DD).GT.H1) THEN
      TB(I,J,K)=T1
      ELSE
      TB(I,J,K)=A0/(DD-B0)       
      END IF
      END IF
      TB(I,J,K)=TB(I,J,K)-TBIAS
  200 CONTINUE
C     TB(I,J,KB)=TB(I,J,KBM1)
  100 CONTINUE

      print*,tbias
      OPEN(25,FILE='controltem.dat')
      WRITE(25,*) TB
      CLOSE(25)
C
C	PRINT*,'Temperature profile at center of domain'
C	PRINT*,'           K      H              TB'
C	DO K=1,KBM1
C       PRINT*,K,H(IM/2,JM/2)*ZZ(K),TB(IM/2,JM/2,K)+TBIAS
C      ENDDO
C
C
C      DO 150 K=1,KB
C
      RETURN
      END
      SUBROUTINE VERINTERP(n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi 
c   interpolating between the n values y at the points x
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp
      integer n,ni,ii

            do i=1,ni
      if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *   (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
        if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
        else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
        end if
      else
        do j=1,n-1
        tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
        end do
c
         do j=1,n-1
           if(tmp(j).le.0) then
             ii=j
           end if
         end do
c
        yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))
     1        /abs(x(ii+1)-x(ii))
      end if
            end do

      return
      end
      SUBROUTINE INTERP1d(mask,n,ni,x,y,xi,yi)
c--------------------------------------------------------------
c   This subroutine determines ni values yi at the points xi 
c   interpolating between the n values y at the points x
c   values equal to mask are ignored
c--------------------------------------------------------------
      real x(n),y(n),xi(ni),yi(ni),tmp(500)
      real cmp,mask
      integer n,ni,ii
c
      do i=1,ni
        if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.
     *  (xi(i).lt.x(1).and.xi(i).lt.x(n))) then
          if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
            if(x(n).gt.x(1)) then
              yi(i)=y(n)
            else
              yi(i)=y(1)
            end if
          else
            if(x(n).gt.x(1)) then
              yi(i)=y(1)
            else
              yi(i)=y(n)
            end if
          end if
        else
          do j=1,n-1
            tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
          end do
          do j=1,n-1
            if(tmp(j).le.0) ii=j
          end do
          if(y(ii).eq.mask.or.y(ii+1).eq.mask) then
            yi(i)=mask
          else
            yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))/
     1      abs(x(ii+1)-x(ii))
          end if
        end if
      end do
c
      return
      end
      SUBROUTINE VERTVL(DTI2)
       INCLUDE 'comblk.h'    
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
C
C
C CALCULATE NEW VERTICAL VELOCITY
C
C REESTABLISH BOUNDARY CONDITIONS
      DO 100 K=1,KBM1
      DO 100 J=1,JM
      DO 100 I=2,IM
 100  XFLUX(I,J,K)
     1 =.25E0*(DY(I,J)+DY(I-1,J))*(DT(I,J)+DT(I-1,J))*U(I,J,K)
      DO 120 K=1,KBM1
      DO 120 J=2,JM
      DO 120 I=1,IM
  120 YFLUX(I,J,K)
     1  =.25E0*(DX(I,J)+DX(I,J-1))*(DT(I,J)+DT(I,J-1))*V(I,J,K)
C
      DO 125 J=1,JM
      DO 125 I=1,IM
125   W(I,J,1)=0.E0
      DO 710 K=1,KBM1
      DO 710 J=2,JMM1
      DO 710 I=2,IMM1
 710  W(I,J,K+1)=W(I,J,K)
     1    +DZ(K)*((XFLUX(I+1,J,K)-XFLUX(I,J,K)
     2            +YFLUX(I,J+1,K)-YFLUX(I,J,K))/(DX(I,J)*DY(I,J))
     3                        +(ETF(I,J)-ETB(I,J))/DTI2 )
      RETURN
      END
      SUBROUTINE VORT(ADVUA,ADVVA,ADVUU,ADVVV,TRNU,TRNV,IIM,JJM)  
       INCLUDE 'comblk.h'     
      REAL JBARX,JBARY    
      DIMENSION ADVUA(IIM,JJM),ADVVA(IIM,JJM),ADVUU(IIM,JJM),
     1    ADVVV(IIM,JJM),TRNU(IIM,JJM),TRNV(IIM,JJM)
      DIMENSION TXSURF(IM,JM),TYSURF(IM,JM),CTSURF(IM,JM)   
      DIMENSION TXBOT(IM,JM),TYBOT(IM,JM),CTBOT(IM,JM)   
      DIMENSION PVFX(IM,JM),PVFY(IM,JM),CPVF(IM,JM)  
      DIMENSION JBARX(IM,JM),JBARY(IM,JM),CJBAR(IM,JM)  
      DIMENSION ADVX(IM,JM),ADVY(IM,JM),CADV(IM,JM)   
      DIMENSION TENX(IM,JM),TENY(IM,JM),CTEN(IM,JM)   
      DIMENSION ELGX(IM,JM),TOTX(IM,JM)  
      DIMENSION CTOT(IM,JM)                           
      alpha =0.225        
      DO I=2,IM     
      DO J=2,JM     
      DMX=0.5*(D(I,J)+D(I-1,J))  
      DMY=0.5*(D(I,J)+D(I,J-1))  
      TXSURF(I,J)=WUSURF(I,J)/DMX*DUM(I,J)                         
      TYSURF(I,J)=WVSURF(I,J)/DMY*DVM(I,J)                     
      TXBOT(I,J)=-WUBOT(I,J)/DMX*DUM(I,J)                 
      TYBOT(I,J)=-WVBOT(I,J)/DMY*DVM(I,J)                 
      PVFY(I,J)=+.25*(  COR(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     2           +COR(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     3       /DMY*DVM(I,J)      
      PVFX(I,J)=-.25*(  COR(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     2           +COR(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     3       /DMX*DUM(I,J)     
      JBARX(I,J)=TRNU(I,J)/(ARU(I,J)*DMX)*DUM(I,J)       
      JBARY(I,J)=TRNV(I,J)/(ARV(I,J)*DMY)*DVM(I,J)     
      ADVX(I,J)=(ADVUU(I,J)+ADVUA(I,J))/(ARU(I,J)*DMX)*DUM(I,J)  
      ADVY(I,J)=(ADVVV(I,J)+ADVVA(I,J))/(ARV(I,J)*DMY)*DVM(I,J)       
      TENX(I,J)=(UAF(I,J)*(H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))
     1  -UAB(I,J)*(H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J)))/(4.*DTE*DMX)       
     2   *DUM(I,J)  
      TENY(I,J)=(VAF(I,J)*(H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))
     1  -VAB(I,J)*(H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1)))/(4.*DTE*DMY)  
     2   *DVM(I,J)   
      ELGX(I,J)=.25E0*GRAV*(DY(I,J)+DY(I-1,J))*(D(I,J)+D(I-1,J))
     4             *( (1.E0-2.E0*ALPHA)*(EL(I,J)-EL(I-1,J))
     4            +ALPHA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
      ELGX(I,J)=ELGX(I,J)/(ARU(I,J)*DMX)*DUM(I,J)  
      TOTX(I,J)=TXSURF(I,J)+TXBOT(I,J)+PVFX(I,J)+JBARX(I,J)
     1   +ADVX(I,J)+TENX(I,J)+ELGX(I,J)      
      ENDDO     
      ENDDO     
      CALL CURL(TXSURF,TYSURF,DX,DY,IM,JM,CTSURF)  
      CALL CURL(TXBOT,TYBOT,DX,DY,IM,JM,CTBOT)  
      CALL CURL(PVFX,PVFY,DX,DY,IM,JM,CPVF)   
      CALL CURL(JBARX,JBARY,DX,DY,IM,JM,CJBAR)   
      CALL CURL(TENX ,TENY ,DX,DY,IM,JM,CTEN )   
      CALL CURL(ADVX,ADVY,DX,DY,IM,JM,CADV)    
      DO I=2,IM   
      DO J=2,JM   
      CTOT(I,J)=CTSURF(I,J)+CTBOT(I,J)+CPVF(I,J)+CJBAR(I,J)
     1   +CTEN(I,J)+CADV(I,J)
      ENDDO       
      ENDDO       
      RETURN   
      END  
      SUBROUTINE WIND(filename,startdate)
      INCLUDE 'comblk.h'
      INTEGER PLN
      PARAMETER(PLN=100)
c
      integer hour, lat, long, mx, rmw
      integer*4 startdate,date
      integer day,month,year
      integer garb(5),Rd1(4),Rd2(4)
      character*19 name
      character*15 filename
      character*1 letter
c
      DIMENSION X(PLN),Y(PLN),TM(PLN),PRES(PLN),PRES0(PLN),
     *      RMAXa(PLN),WSPMAX(PLN),
     *      R18v(PLN,4),R26v(PLN,4),Rref18v(5),Rref26v(5),alphv(5)
      DIMENSION RAD(14),WS(14),RADM(14),WSM(14),ANGL(14)
      REAL CMP,T1,T2,F0,F1,L0,L1,REARTH,R,A7,B,E,DELP,x0,y0
      REAL LATMIN,LATMAX,LONGMIN,LONGMAX
      REAL DELTAX,DELTAX1,DELTAY,DELTAY1,DXDY,julday
      COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX

      DATA RM,R0/60.E3,480.E3/,RHO_0/1024.E0/
      DATA RAD/0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75
     * ,8.1,10.8,13.5/
      DATA WS/0.,.1,.5,.8,.95,1.,.97,.72,.54,.44,.4,.36
     * ,.27,.23/
      DATA ANGL/0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,
     * 21.,21.,21.,21./
      print*,'In subroutine WIND ...'
c
      WIND_SCALE=0.8
      ROA=1.28
      RMAX=50.E3
      PI=3.1415927
      E=exp(1.)
c
c----------------------- Reading message file --------------------
  17  format(A19,I6,1x,I4,1x,I3,1x,I4,1x,I3,1x,I3,3I5,
     * 1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,
     * 1x,I4,1x,I4,1x,I4)
      print*,'reading file ',filename
      open(15,file=filename,status='old')
      end=0.
      I=0
      do while(end.eq.0)
        read(15,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
        write(6,17) name,date,hour,lat,long,garb,mx,rmw,Rd1,Rd2
        if(date.eq.0) goto 20
        I=I+1
c
        date=date*100+hour/100
        call date2day(year,julday,date)
c
        TM(i)=julday
        X(i)=-long/10.
        Y(i)=lat/10.
        PRES(i)=float(garb(3))
        PRES0(i)=float(garb(4))
        WSPMAX(i)=float(mx)
        RMAXa(i)=float(rmw)
        do n=1,4
          n1=n+(1-mod(n,2))*sign(2,3-n)
          R18v(i,n)=Rd1(n1)
          R26v(i,n)=Rd2(n1)
          if(wspmax(i).le.26.or.R26v(i,n).le.RMAXa(i)) R26v(i,n)=-999
          if(wspmax(i).le.18.or.R18v(i,n).le.RMAXa(i)) R18v(i,n)=-999
          if(R26v(i,n).gt.R18v(i,n)) R26v(i,n)=-999
        end do
      end do
  20  end=1.
      ipmax=I
      print*,'Number of hurricane path datapoints read: ',ipmax
      print*,'tm=',(tm(i),i=1,ipmax)
      print*,'year=',year
      print*,'startdate=',startdate
      close(15)
c
c--------------------- Calculating starting day -----------------
c
      call date2day(year,julday,startdate)
      do i=1,ipmax
        TM(i)=TM(i)-julday
      end do
C++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTERPOLATION OF HURRICANE PATH TO DETERMINE THE CURRENT POSITION

      CMP=TM(ipmax)
      if(cmp.lt.time.or.time.lt.tm(1)) then
        print*,'NO HURRICANE PATH DATA FOR THIS TIME'
        print*,'   time=',time
        print*,'   tm(1)=',tm(1)
        print*,'   tm(ipmax)=',tm(ipmax)
        print*,'   day=',julday
        RETURN
      end if

      call verinterp(ipmax,1,TM,X,TIME,F0)
      call verinterp(ipmax,1,TM,Y,TIME,L0)
      call verinterp(ipmax,1,TM,PRES,TIME,PRES1)
      call verinterp(ipmax,1,TM,PRES0,TIME,PRES2)
      DELP=(PRES2-PRES1)*100.
      call verinterp(ipmax,1,TM,WSPMAX,TIME,WSMAX)
      WSMAX=WSMAX*WIND_SCALE
      WS18=18*WIND_SCALE
      WS26=26*WIND_SCALE
      call verinterp(ipmax,1,TM,RMAXa,TIME,RMAX)
      RMAX=RMAX*1.e3
      do n=1,4
        call interp1d(-999.0,ipmax,1,TM,R18v(1,n),TIME,Rref18v(n))
        if(Rref18v(n).ne.-999) Rref18v(n) = Rref18v(n)*1.e3
        call interp1d(-999.0,ipmax,1,TM,R26v(1,n),TIME,Rref26v(n))
        if(Rref26v(n).ne.-999) Rref26v(n) = Rref26v(n)*1.e3
        alphv(n) = (n-1)*pi/2
      end do
      do n=2,6
        n1=mod(n-1,4)+1
        nm1=mod(n-2,4)+1
        np1=mod(n,4)+1
        if(Rref18v(n1).eq.-999) then
          if(Rref18v(nm1).ne.-999) then
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=0.5*(Rref18v(nm1)+Rref18v(np1))
            else
              Rref18v(n1)=Rref18v(nm1)
            end if
          else
            if(Rref18v(np1).ne.-999) then
              Rref18v(n1)=Rref18v(np1)
            else
              Rref18v(n1)=-999
            end if
          end if
        end if
        if(Rref26v(n1).eq.-999) then
          if(Rref26v(nm1).ne.-999) then
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=0.5*(Rref26v(nm1)+Rref26v(np1))
            else
              Rref26v(n1)=Rref26v(nm1)
            end if
          else
            if(Rref26v(np1).ne.-999) then
              Rref26v(n1)=Rref26v(np1)
            else
              Rref26v(n1)=-999
            end if
          end if
        end if
      end do
c
      Rref18v(5) = Rref18v(1)
      Rref26v(5) = Rref26v(1)
      alphv(5) = alphv(4)+pi/2
c
      print*,'Time=',Time
      print*,'Current hurricane position (x,y): ',f0,l0
      print*,'WSMAX=',WSMAX,'; DELP=',DELP,'; RMAX=',RMAX
      print*,'Rref18v=',Rref18v
      print*,'Rref26v=',Rref26v
c
      x0=f0
      y0=l0
      F0=F0*2.*PI/360
      L0=L0*2.*PI/360
C--- F0,L0 ARE THE CARRENT (LONGITUDE,LATTITUDE) COORDINATES OF THE HURRICANE
C
C--- CALCULATING UTX AND UTY (HURRICANE SPEED)
c
      cmp=tm(ipmax)
      do i=1,ipmax
        if(abs(tm(i)-time).le.cmp) then
          cmp=abs(tm(i)-time)
          ii=i
        end if
      end do
      if((tm(ii)-time) .le. 0. .and. ii .ne. ipmax) then
        t1=tm(ii)
        t2=tm(ii+1)
        x1=x(ii)
        x2=x(ii+1)
        y1=y(ii)
        y2=y(ii+1)
      else
        t2=tm(ii)
        t1=tm(ii-1)
        x2=x(ii)
        x1=x(ii-1)
        y2=y(ii)
        y1=y(ii-1)
      end if
      deltax1=rearth*cos(l0)*(x2-x1)*2.*pi/360
      deltay1=rearth*(y2-y1)*2.*pi/360
      utx=deltax1/((t2-t1)*24.*3600.)
      uty=deltay1/((t2-t1)*24.*3600.)
c
      print*,'utx,uty: ',utx,uty
C
C--- CALCULATING PARAMETERS FOR WIND PROFILE FORMULA
c
      B=WSMAX**2*E*ROA/DELP
      A7=RMAX**B
      print*,'B= ',B
c
      do i=1,14
        RADM(I)=RMAX*RAD(I)
      end do
C
      DO 350 J=1,JM
      DO 351 I=1,IM

C  CALCULATING U-WIND STRESS FOR I,J POINT OF U-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/((IM-0.5)-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY) +
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

c      if(R.GT.RADM(14).OR.R.EQ.0) then
c        k=14
c      else
c        k=1
c        do while(R.GE.RADM(k+1))
c          k=k+1
c        end do
c      end if

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
c      WF=WND
c      WR=-WND*TAN(RANGL)
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.LE.35.) then
      WUSURF(I,J)=WUSURF(I,J)-CD*ROA*WM*WX/RHO_0
      else
      WUSURF(I,J)=WUSURF(I,J)-(3.3368+
     1      (WM-34.0449)**0.3)*WX/(WM*RHO_0)
      end if


C  CALCULATING V-WIND STRESS FOR I,J POINT OF V-VELOCITY GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/((JM-0.5)-1)
      L1=L1*2.*PI/360
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.le.35.) then
      WVSURF(I,J)=WVSURF(I,J)-CD*ROA*WM*WY/RHO_0
      else
      WVSURF(I,J)=WVSURF(I,J)-(3.3368+
     1      (WM-34.0449)**0.3)*WY/(WM*RHO_0)
      end if

C  CALCULATING WIND STRESS FOR I,J POINT OF DEPTH GRID
      F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      F1=F1*2.*PI/360.
      L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      L1=L1*2.*PI/360.
      DELTAX=REARTH*COS(L0)*(F1-F0)
      if (deltax .eq. 0.) deltax=1.e-8
      DELTAY=REARTH*(L1-L0)
      DXDY=DELTAX*DELTAY
      R=SQRT(DELTAX**2+DELTAY**2)
      alpha=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+
     1      (1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
      if(alpha.ge.pi/4) then
        alpha = alpha - pi/4
      else
        alpha = alpha - pi/4 + 2*pi
      end if
      call verinterp(5,1,alphv,Rref18v,alpha,Rref18)
      call verinterp(5,1,alphv,Rref26v,alpha,Rref26)

      call verinterp(14,1,radm,angl,R,RANGL)

C    CALCULATING WIND SPEED 

      if(Rref18.le.0.and.Rref26.le.0) then
        WND=SQRT(A7*B*DELP*EXP(-A7/R**B)/(ROA*R**B)+R**2*
     1  COR(I,J)**2/4.)-R*COR(I,J)/2.
        UTXa=UTX/2.
        UTYa=UTY/2.
      else
        WND=EXPWND(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
        UTXa=0.
        UTYa=0.
      end if
c      RANGL=ANGL(K)*PI/180.
      RANGL=RANGL*PI/180.
      WF=WND*COS(RANGL)
      WR=-WND*SIN(RANGL)
      WX=WR*(DELTAX)/R-WF*(DELTAY)/R+UTXa
      WY=WF*(DELTAX)/R+WR*(DELTAY)/R+UTYa
c      WINDX(i,j)=WX
c      WINDY(i,j)=WY
      WM=SQRT(WX**2+WY**2)
      IF(WM.LT.10.) CD=1.14*1.E-3
      IF(WM.GE.10.) CD=(0.49+.065*WM)*1.E-3
      IF(WM.GT.35.) THEN
      TMAX=3.3368+(WM-34.0449)**0.3
      TAUX(I,J)=TAUX(I,J)+TMAX*WX/WM
      TAUY(I,J)=TAUY(I,J)+TMAX*WY/WM
      ELSE
      TAUY(I,J)=TAUY(I,J)+CD*ROA*WM*WY
      TAUX(I,J)=TAUX(I,J)+CD*ROA*WM*WX

C      WTSURF(I,J)=30.*WM/(4000*1024.)

      END IF

 351  CONTINUE
 350  CONTINUE

 230  FORMAT(151e16.7)
 231  FORMAT(10F8.3)
      print*,'Exiting WIND ...'
      RETURN
      END
	subroutine pdens(rhor)

      include 'comblk.h'

	real cr(im,jm,kb),tr(im,jm,kb),sr(im,jm,kb),p(im,jm,kb)
	real rhor(im,jm,kb)

	equivalence (a,cr),(ee,tr),(uf,sr),(c,p)

c  this subroutine computes density-1.025
c  t = potential temperature

	do k=1,kbm1
	  do j=1,jm
	    do i=1,im
	      tr(i,j,k)=t(i,j,k)+tbias
	      sr(i,j,k)=s(i,j,k)+sbias

c  approximate pressure in units of bars

	      p(i,j,k)=-grav*1.025*zz(k)*dt(i,j)*0.01

	      rhor(i,j,k) = 999.842594 + 6.793952e-2*tr(i,j,k)
     #- 9.095290e-3*tr(i,j,k)**2 + 1.001685e-4*tr(i,j,k)**3
     #- 1.120083e-6*tr(i,j,k)**4 + 6.536332e-9*tr(i,j,k)**5

	      rhor(i,j,k)=rhor(i,j,k) + (0.824493 - 4.0899e-3*tr(i,j,k) 
     #+ 7.6438e-5*tr(i,j,k)**2 - 8.2467e-7*tr(i,j,k)**3
     #+ 5.3875e-9*tr(i,j,k)**4) * sr(i,j,k)
     #+ (-5.72466e-3 + 1.0227e-4*tr(i,j,k)
     #- 1.6546e-6*tr(i,j,k)**2) * sr(i,j,k)**1.5
     #+ 4.8314e-4 * sr(i,j,k)**2

c  at this point, rhor has potential density, referenced to the surface
c  this next section takes it down to pressure p

c	      cr(i,j,k)
c     #=1449.1+.0821*p(i,j,k)+4.55*tr(i,j,k)-.045*tr(i,j,k)**2
c     #+1.34*(sr(i,j,k)-35.)
c	      rhor(i,j,k)=rhor(i,j,k) + 1.e5*p(i,j,k)/cr(i,j,k)**2 
c     #*(1.-2.*p(i,j,k)/cr(i,j,k)**2)

	      rhor(i,j,k)=(rhor(i,j,k)-1000.)*fsm(i,j)
	    enddo
	  enddo
	enddo

	do j=1,jm
	  do i=1,im
	    rhor(i,j,kb)=0.e0
	  enddo
	enddo

	return
	end
	subroutine rikh

c  this subroutine calculates a Richardson number dependent kh (and km)

      include 'comblk.h'
        real sigthet(im,jm,kb)

	dimension shear2(im,jm,kb),uz(im,jm,kb),vz(im,jm,kb)
     #,rz(im,jm,kb),dh(im,jm)

	equivalence (uz,a),(vz,c),(dh,tps)


        do j=1,jm
          do i=1,im
            dh(i,j)=h(i,j)+etf(i,j)
          enddo
        enddo

        call pdens(sigthet)

	do k=2,kbm1
	  do j=1,jm
	    do i=1,im
	      tmp=1./(dzz(k-1)*dh(i,j))
	      uz(i,j,k)=(u(i,j,k-1)-u(i,j,k))*tmp
	      vz(i,j,k)=(v(i,j,k-1)-v(i,j,k))*tmp
	      rz(i,j,k)=
     #-grav*(sigthet(i,j,k-1)-sigthet(i,j,k))*tmp*.001
c  *.001 is for 1/rho0 = 1/1000.
	    enddo
	  enddo
	enddo

	do k=2,kbm1
	  do j=2,jmm1
	    do i=2,imm1
	      shear2(i,j,k)=
     #(.5*(uz(i+1,j,k)+uz(i,j,k)))**2
     #+(.5*(vz(i,j+1,k)+vz(i,j,k)))**2
	    enddo
	  enddo
	enddo

	do k=2,kbm1
	  do j=2,jmm1
	    do i=2,imm1
	      tmp=rz(i,j,k)/(shear2(i,j,k)+1.e-8)  !=richardson#
	      tmp2=-zz(k)*dh(i,j)
	      if(tmp2.gt.500.) tmp2=500.
	      if(tmp.lt.1.) then
	        kh(i,j,k)=sqrt((1.-tmp)*shear2(i,j,k))*tmp2*tmp2
	      else
	        kh(i,j,k)=0.
	      endif
	    enddo
	    kh(1,j,k)=kh(2,j,k)
	    kh(im,j,k)=kh(imm1,j,k)
	  enddo
	  do i=1,im
	    kh(i,1,k)=kh(i,2,k)
	    kh(i,jm,k)=kh(i,jmm1,k)
	  enddo
	enddo

	do j=1,jm
	  do i=1,im
	    kh(i,j,1)=kh(i,j,2)
	    kh(i,j,kb)=kh(i,j,kbm1)
	  enddo
	enddo

	do k=1,kb
	  do j=1,jm
	    do i=1,im
	      km(i,j,k)=kh(i,j,k)
	    enddo
	  enddo
	enddo

	return
	end

	SUBROUTINE ATMOS2OCEAN
c-----------------------------------------------------------------------
c   This subroutine interpolates atmospheric fluxes onto the ocean grid
c-----------------------------------------------------------------------
      INCLUDE 'comblk.h'
      INCLUDE 'comblk1.h'
C
       IINT = IOSTEP
C
C--- CONVERT THE HURRICANE WIND STRESS AND HEAT FLUX  TO SI UNITS
C
      do J=1,JMA
        do I=1,IMA
          XSTR(I,J)=-(0.1*XSTR(I,J))/RHO_0
          YSTR(I,J)=-(0.1*YSTR(I,J))/RHO_0
          HFLX(I,J)=HFLX(I,J)*1.0E-3/(RHO_0*4200)
          RFSW(I,J)=RFSW(I,J)*1.0E-3/(RHO_0*4200)
        end do
      end do
c
C   PASSING THE WIND STRESS AND THE NET HEAT FLUX FROM
C      THE ATMOSPHERE MODEL TO THE OCEAN MODEL
C
      DO 222 I=1,IM
      XSU(I)=(LONGMIN+360.)*PI/180+(LONGMAX-LONGMIN)*PI*
     1       (I-1-0.5)/(180*(IM-1))
      XSV(I)=(LONGMIN+360.)*PI/180+(LONGMAX-LONGMIN)*PI*
     1       (I-1)/(180*(IM-1))
 222  XST(I)=XSV(I)
      DO 223 J=1,JM
      YSV(J)=LATMIN*PI/180+(LATMAX-LATMIN)*PI*
     1       (J-1-0.5)/(180*(JM-1))
      YSU(J)=LATMIN*PI/180+(LATMAX-LATMIN)*PI*
     1       (J-1)/(180*(JM-1))
 223  YST(J)=YSU(J)
C
      PRINT*, 'XLON ', XLON(1),XLON(IMA)
      PRINT*, 'YLAT ', YLAT(1),YLAT(JMA)
      PRINT*, 'XST ',  XST(1) ,XST(IM)
      PRINT*, 'YST ',  YST(1) ,YST(JM)
C      
      SATMMAX=0.
      DO I=1,IMA
      DO J=1,JMA
      SATMMOD=SQRT(XSTR(I,J)**2+YSTR(I,J)**2)
      IF (SATMMOD.GT.SATMMAX) SATMMAX=SATMMOD
      ENDDO
      ENDDO
      PRINT*,'ATM MODEL: MAXIMUM WIND STRESS ', SATMMAX
c
      CALL INTERP1(XSTR,WUSURF,IMA,JMA,IM,JM,XLON,YLAT,XSU,YSU)
      CALL INTERP1(YSTR,WVSURF,IMA,JMA,IM,JM,XLON,YLAT,XSV,YSV)
      CALL INTERP1(XSTR,TAUX,IMA,JMA,IM,JM,XLON,YLAT,XST,YST)
      CALL INTERP1(YSTR,TAUY,IMA,JMA,IM,JM,XLON,YLAT,XST,YST)
      CALL INTERP1(HFLX,WTSURF,IMA,JMA,IM,JM,XLON,YLAT,XST,YST)
      CALL INTERP1(RFSW,SWRAD,IMA,JMA,IM,JM,XLON,YLAT,XST,YST)
c
      STRMAX=0.
      SHFMAX=1000.
      DO I=1,IM
      DO J=1,JM
      STRMOD=SQRT(WUSURF(I,J)**2+WVSURF(I,J)**2)
      IF (STRMOD.GT.STRMAX) STRMAX=STRMOD
      IF (WTSURF(I,J).LT.SHFMAX) SHFMAX=WTSURF(I,J)
C      WTSURF(I,J)=0
      ENDDO
      ENDDO
      PRINT*,'OCEAN MODEL: MAXIMUM WIND STRESS ', STRMAX
      PRINT*,'OCEAN MODEL: MINIMUM HEAT FLUX   ', SHFMAX
c
c      DO 85 J=2,JMM1
c      DO 85 I=2,IMM1
c   85 SWRAD(I,J)=0.E0
c
      do j=1,jm
        do i=1,im
          taux(i,j)=-taux(i,j)*rho_0
          tauy(i,j)=-tauy(i,j)*rho_0
        enddo
      enddo
c
      return
      end
	function avrsst(xln, ylt)

      include 'comblk.h'
      include 'comblk1.h'
      real*8 avrsst
      real xln,ylt,tavr,counter,x1,y1,x0,y0,r,deltax,deltay

      tavr=0.0
      counter=0.0
      do j=1,jm
        do i=1,im
          x1=( LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1) )*pi/180.
          y1=( LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1) )*pi/180.
          x0=xln*pi/180.
          y0=ylt*pi/180.
          DELTAX=REARTH*COS(y0)*(x1-x0)
          DELTAY=REARTH*(y1-y0)
          R=SQRT(DELTAX**2+DELTAY**2)
          if(r.lt.100.e3) then
            tavr=tavr+(T(i,j,1)-SSTIN(i,j))*fsm(i,j)
            counter=counter+fsm(i,j)
          end if
        end do
      end do
      if(counter.gt.0) then
        avrsst=tavr/counter
      else
        avrsst=0.0
      end if

      print*,'avrsst: xln=',xln,', ylt=',ylt,', avrsst=',avrsst

      return
      end

      SUBROUTINE INTERP1(OLD,NEW,IMX,JMX,INMAX,JNMAX,XO,YO,XN,YN)   
      REAL NEW     
      DIMENSION OLD(IMX,JMX),     NEW(INMAX,JNMAX),        IDI(INMAX),
     1          IDJ(JNMAX),       WTI(INMAX),              WTJ(JNMAX),
     2          XO(IMX),     YO(JMX),    XN(INMAX),        YN(JNMAX) 
      CALL WGHTS(IMX,JMX,INMAX,JNMAX,XO,YO,XN,YN,WTI,WTJ,IDI,IDJ)  
      DO 500 J=1,JNMAX
      DO 500 I=1,INMAX
* THIS CONDITION IS FOR Y AXES DIRECTED NORTHWARD     
      IF (YN(J).GE.YO(JMX).OR.YN(J).LE.YO(1)) GO TO 150 
* THIS CONDITION IS FOR X AXES DIRECTED WESTWARD
      IF (XN(I).GE.XO(IMX).OR.XN(I).LE.XO(1)) GO TO 150   
       
      WJ1=1.-WTJ(J)      
      WJ=WTJ(J)     
      IJ=IDJ(J)     
      IJ1=IDJ(J)+1      
  
      II=IDI(I)       
      II1=IDI(I)+1   
      NEW(I,J)=WJ1*(1.-WTI(I))*OLD(II,IJ)+ 
     1         WJ1*WTI(I)*OLD(II1,IJ)+            
     2         WJ*(1.-WTI(I))*OLD(II,IJ1)+     
     3         WJ*WTI(I)*OLD(II1,IJ1)
      GO TO 500
 150  NEW(I,J)=0.     
  500 CONTINUE  
      RETURN     
      END
      SUBROUTINE INTERP2(OLD,NEW,IMX,JMX,INMAX,JNMAX,XO,YO,XN,YN)   
      REAL NEW     
      DIMENSION OLD(IMX,JMX),     NEW(INMAX,JNMAX),        IDI(INMAX),
     1          IDJ(JNMAX),       WTI(INMAX),              WTJ(JNMAX),
     2          XO(IMX),     YO(JMX),    XN(INMAX),        YN(JNMAX)
C
      REAL FLAG(4)
C 
      CALL WGHTS(IMX,JMX,INMAX,JNMAX,XO,YO,XN,YN,WTI,WTJ,IDI,IDJ)  
      DO 500 J=1,JNMAX
      DO 500 I=1,INMAX
* THIS CONDITION IS FOR Y AXES DIRECTED NORTHWARD     
      IF (YN(J).GE.YO(JMX).OR.YN(J).LE.YO(1)) GO TO 150 
* THIS CONDITION IS FOR X AXES DIRECTED WESTWARD
      IF (XN(I).GE.XO(IMX).OR.XN(I).LE.XO(1)) GO TO 150   
C
      IJ=IDJ(J)     
      IJ1=IDJ(J)+1        
      II=IDI(I)       
      II1=IDI(I)+1   
C
      WJ1=(YO(IJ1)-YO(IJ))-WTJ(J)      
      WJ=WTJ(J)
C
      if(OLD(II,IJ).eq.-99.) then
      flag(1)=0.
      else
      flag(1)=1.
      endif     
      if(OLD(II1,IJ).eq.-99.) then
      flag(2)=0.
      else
      flag(2)=1. 
      endif
      if(OLD(II,IJ1).eq.-99.) then
      flag(3)=0.
      else
      flag(3)=1. 
      endif
      if(OLD(II1,IJ1).eq.-99.) then
      flag(4)=0.
      else
      flag(4)=1. 
      endif
C
      s=(flag(1)*WJ1*((XO(II1)-XO(II))-WTI(I))+
     5         flag(2)*WJ1*WTI(I)+
     6         flag(3)*WJ*((XO(II1)-XO(II))-WTI(I))+
     7         flag(4)*WJ*WTI(I))
      if(s.eq.0.) then
      NEW(I,J)=29.3
      else
      NEW(I,J)=(flag(1)*WJ1*((XO(II1)-XO(II))-WTI(I))*OLD(II,IJ)+ 
     1         flag(2)*WJ1*WTI(I)*OLD(II1,IJ)+            
     2         flag(3)*WJ*((XO(II1)-XO(II))-WTI(I))*OLD(II,IJ1)+     
     3         flag(4)*WJ*WTI(I)*OLD(II1,IJ1))/s
      end if
      GO TO 500
 150  NEW(I,J)=29.3    
  500 CONTINUE  
      RETURN     
      END 
      SUBROUTINE OCEANINIT
C
         INCLUDE 'comblk.h'    
         INCLUDE 'comblk1.h'    
c-------- 03-08-01 include for reading smoothed topogr.
      CHARACTER*55  FTTOPOG
c
      REARTH=6371.E3
      LATMIN=10.
      LATMAX=47.5
      LONGMIN=-98.5
      LONGMAX=-50.
c
      PI=3.1415927E0
      RAMP=1.E0
      SMALL=1.E-10
      TIME0=0.E0
      BETA=1.98E-11
      GRAV=9.806E0
      UMOL=2.E-5
C--- specify bathymetery
       H_MAX=5500.E0
C
      IINT=0
      TIME=0.
c
c  Reading gfs data dimensions Biju Thomas on 01/27/02
      READ(23,224)IM1,JM1
      rewind(23)
 224  FORMAT(1X,2I5)
c
c---------- falk 03-26-02 
      print *,' read gfs data dimensions IM1,JM1=',IM1,JM1
c---------------------
c
C--- Model contral parameters----------------------------
C--- see file PARAMETERS.inp for meanings of parameters
191   format(a15)
c---------- falk 07-06-01 use description of file in script
C     OPEN(10,FILE='PARAMETERS.inp',STATUS='OLD')
      READ(10,*) MODE
      READ(10,*) NBC
      READ(10,*) TIME0
      READ(10,*) IMAY
      READ(10,*) sstsource
      READ(10,*) startdate
      READ(10,191) pathname
      READ(10,191) pathname1
      READ(10,*) DTI
      READ(10,*) ISPLIT
      READ(10,*) NREAD
      READ(10,*) IPRTH1
      READ(10,*) INOWINDH,VC
      READ(10,*) IDAMP,ISMOTH,SMH
      READ(10,*) wndg
      READ(10,*) IHOURS
      IF(NREAD.GT.0) READ(10,'(A)') FRSTI
      rewind(10)
C     CLOSE(10)
      IF(IDAMP.EQ.1) PRINT*,'APPLY BOUNDARY DAMPING'
      IF(ISMOTH.EQ.1) PRINT*,'APPLY U,V,T,S SMOOTHING. SMH= ',SMH
C
      ISPADV=5
      TPRNU=1.   
      SMOTH=0.1
      HORCON=0.1
      TBIAS=10.
      SBIAS=35.
      RHO_0=1024.e0
      ISWTCH=100000
      IPRTD2=1
c
      do j=1,jm
        do i=1,im
          nbc2d(i,j)=nbc
          nbc2ds(i,j)=1
        end do
      end do
C--------------------------------------------------------------------------
C     READ NAMELIST (DISCONNECTED HERE) FOR PROBLEM PARAMETERS: 
C       MODE = 2; 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
C              3; 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
C              4; 3-D CALCULATION WITH T AND S HELD FIXED
C       NBC = TYPE OF THERMO. B.C.; SEE SUBROUTINE PROFT
C       DTI = INTERNAL TIME STEP
C       DTE = EXTERNAL TIME STEP
C       ISPLIT = DTI/DTE
C       NREAD=0, NO RESTART INPUT FILE; NREAD=1, RESTART
C       IPRTD1 = PRINT INTERVAL IN DAYS; AFTER ISWTCH DAYS, PRINT
C                  INTERVAL = IPRTD2
C       IDAYS = LENGTH OF RUN IN DAYS
C       ISPADV = STEP INTERVAL WHERE EXTERNAL MODE ADVECTIVE TERMS ARE
C                  NOT UPDATED
C       HORCON = CONSTANT IN SMAGORINSKY HORIZONTAL DIFFUSIVITY
C       TPNU = HORIZONTAL DIFFUSIVITY PRANDTL NUMBER
C       SMOTH = CONSTANT IN TIME SMOOTHER TO PREVENT SOLUTION SPLITTING
C
C     READ(5,PRMTR)
C--------------------------------------------------------------------------
c
      if(IMAY.eq.0)then
      print*
      print*,'RUN WITHOUT MELLOR & YAMADA TURBULENT MIXING'
      print*
      endif
      if(IMAY.eq.1)then
      print*
      print*,'RUN WITH MELLOR & YAMADA TURBULENT MIXING'
      print*
      endif
c
      DTE=DTI/FLOAT(ISPLIT)
      DTE2=DTE*2
      DTI2=DTI*2
      IPRINT=IPRTH1*3600/INT(DTI)
      ISWTCH=ISWTCH*24*3600/INT(DTI)   
      IEND=IHOURS*3600/INT(DTI)
      WRITE(6,7030) MODE,DTE,DTI,IEND,ISPLIT,ISPADV,IPRINT,SMOTH,
     1   HORCON,TPRNU   
 7030 FORMAT(//,' MODE =',I3,
     1 '  DTE =',F7.2,'   DTI =',F9.1,'     IEND =',I6,
     2  '   ISPLIT =',I6,'    ISPADV =',I6,'   IPRINT =',I6,/,
     3  '   SMOTH =',F7.2,'   HORCON =',F7.3,'   TPRNU =',F7.3,/)
      WRITE(6,'('' NREAD ='',I5,//)') NREAD
C
C----------------------------------------------------------------------
C             ESTABLISH PROBLEM CHARACTERISTICS
C          ****** ALL UNITS IN M.K.S. SYSTEM ******
C      F,BLANK AND B REFERS TO FORWARD,CENTRAL AND BACKWARD TIME LEVELS.
C----------------------------------------------------------------------
C
      DAYI=1.E0/86400.E0
      ISPI=1.E0/FLOAT(ISPLIT)
      ISP2I=1.E0/(2.E0*FLOAT(ISPLIT))
C
C     SET UP  GRID AND INITIAL CONDITIONS FOR STAND-ALONE TEST RUN
      CALL DEPTH(Z,ZZ,DZ,DZZ,DZR,KB,KBM1,H_MAX)
      print*,'out of depth...'
      
      do k=1,3
         do j=1,jm
            do i=1,im
            KH(I,J,K)=2.E-2
            KM(I,J,K)=KH(I,J,K)
            end do
         end do
      end do
      do k=4,kb
         do j=1,jm
            do i=1,im
            KH(I,J,K)=2.E-5
            KM(I,J,K)=KH(I,J,K)
            end do
         end do
      end do
      DO 11 K=1,KB
      DO 11 J=1,JM
      DO 11 I=1,IM
      Q2B(I,J,K)=1.E-8   
      Q2LB(I,J,K)=1.E-8   
  11  AAM(I,J,K)=500.
      DO 12 J=1,JM
      DO 12 I=1,IM
      AAM2D(I,J)=AAM(I,J,1)    
  12  CONTINUE
C
C--- ASSIGNING SPACING AND CORIOLIS
C--- STEP IN LONGITUDE(PHI) IS (LONGMAX-LONGMIN)/(IM-1)
C--- STEP IN LATTITUDE(LAMBDA) IS (LATMAX-LATMIN)/(JM-1)
C
      DO 14 J=1,JM
      DO 14 I=1,IM
      PHI=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
      PHI=2.*3.1415927*PHI/360.
      LAMBDA=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
      LAMBDA=2.*3.1415927*LAMBDA/360
      COR(I,J)=4.*3.1415927/(24.*60.*60.)*SIN(LAMBDA)
      DX(I,J)=REARTH*COS(LAMBDA)*(LONGMAX-LONGMIN)/(IM-1)
      DX(I,J)=DX(I,J)*2.*3.1415927/360.
      DY(I,J)=REARTH*(LATMAX-LATMIN)/(JM-1)
      DY(I,J)= DY(I,J)*2.*3.1415927/360.
      ART(I,J)=DX(I,J)*DY(I,J)
  14  CONTINUE
C
c----------- falk 07-06-01 use description of file in a script
c     FTTOPOG='/migr/data/wd20af/forecast_system/ocean/Hdeepgs'
c     OPEN(66,FILE=FTTOPOG,STATUS='unknown',form='formatted')
c
c----------- read smoothed topography and mask
c
      read(66,102) H
      read(66,102) FSM
 102  format(10f7.0)
c     close(66)
c
c--------- check reading
c
      write(6,141) H(im,jm),FSM(im,jm)
 141  format(' check reading H, FSM: H(im,jm),FSM(im,jm)=',2f7.0)
c
c--------- 03-08-01 define DUM,DVM (it was done before in TOPOGR)
c
      DO J=1,JM
        DO I=1,IM
          IF(H(I,J).gt.1.E0.and.FSM(i,j).eq.1.0) then
c           FSM(I,J)=1.E0
            DUM(I,J)=1.E0
            DVM(I,J)=1.E0
          else
c           FSM(I,J)=0.E0
            DUM(I,J)=0.E0
            DVM(I,J)=0.E0
          end if
        end do
      end do
      DO J=1,JM
        DO I=1,IM
          IF(H(I,J).ne.1.0.and.FSM(i,j).eq.0.0) then
            print*,'set h=1 at (i,j)=',i,j,H(i,j),FSM(i,j)
            H(I,J)=1.0
          end if
        end do
      end do
c
      DO J=1,JMM1
         DO I=1,IM
         IF(FSM(I,J).EQ.0.E0.AND.FSM(I,J+1).NE.0.E0) DVM(I,J+1)=0.E0
         end do
      end do
      DO J=1,JM
         DO I=1,IMM1
         IF(FSM(I,J).EQ.0.E0.AND.FSM(I+1,J).NE.0.E0) DUM(I+1,J)=0.E0
         end do
      end do
C
      DO J=1,JM
        DUM(1,J)=0.E0
        DVM(1,J)=0.E0
      end do
c
c--------------- 03-11-02 check north bnd
c
c     print *,' check north bnd  (dvm(i,jm),i=1,im)'
c     write(6,102) (dvm(i,jm),i=1,im)
c     print *,' check north bnd  (dvm(i,jm-1),i=1,im)'
c     write(6,102) (dvm(i,jm-1),i=1,im)
c
c--------- 03-08-01 end changes
C
      DO 13 J=2,JM
      DO 13 I=2,IM
      ARU(I,J)=.25*(DX(I,J)+DX(I-1,J))*(DY(I,J)+DY(I-1,J))
  13  ARV(I,J)=.25*(DX(I,J)+DX(I,J-1))*(DY(I,J)+DY(I,J-1))
      DO 121 J=1,JM
      ARU(1,J)=ARU(2,J)
      ARV(1,J)=ARV(2,J)
 121  continue
      DO 122 I=1,IM
      ARU(I,1)=ARU(I,2)
      ARV(I,1)=ARV(I,2)
 122  continue
C
C The following temperature field will be barotropic only 
C for the flat bottomed case.
      DO 15 K=1,KBM1
      DO 15 J=1,JM
      DO 15 I=1,IM
      SB(I,J,K)=35.0-SBIAS  
  15  CONTINUE
C
c     OPEN(13,FILE='initdata',STATUS='OLD',FORM='UNFORMATTED')
c
      if(sstsource.eq.1.or.sstsource.eq.2) then
        print *,' before SERFTEMPR'
c-----      falk 02-28-02 include IM1,JM1 in SERFTEMPR
        CALL SERFTEMPR(IM1,JM1)
c       CALL MIXSST
        print *,'after SERFTEMPR IM1,JM1=',IM1,JM1
      end if
c
      print *,' before TEMPR      NREAD=',NREAD
      CALL TEMPR(NREAD,sstsource)
c
c--------- falk 04-15-02 
c     print *,'(TB(im,174,k),k=1,kb) east bnd GS center before MIXSST'
c     write(6,201) (TB(im,174,k),k=1,kb)
c---------
c
      CALL DENS(SB,TB,RHO)
      CALL DENS(SMEAN,TMEAN,RMEAN)
C
      do k=1,kb
        do j=1,jm
          do i=1,im
            TMEAN(I,J,K)=TB(I,J,K)
            SMEAN(I,J,K)=SB(I,J,K)
          end do
        end do
      end do
c
      DO 16 J=1,JM
      DO 16 I=1,IM
  16  SSURF(I,J)=SB(I,J,1)
      DO 20 K=1,KBM1
  20  DZR(K)=1./DZ(K)
      PERIOD=DAYI*(2.E0*PI)/COR(IM/2,JM/2)

C A very empirical specification of the roughness parameter follows
      DO 45 J=1,JM
      DO 45 I=1,IM
      Z0B=.01
      CBCMIN=.0025E0
  45  CBC(I,J)=MAX(CBCMIN,.16E0/LOG((ZZ(KBM1)-Z(KB))*H(I,J)
     1        /Z0B)**2)*FSM(I,J)
C
      IP=IM
      ISK=1
C
        DO 47 J=1,JM
        DO 47 I=1,IM
   47   TPS(I,J)=0.5E0/SQRT(1.E0/DX(I,J)**2+1.E0/DY(I,J)**2)
     1         /SQRT(GRAV*H(I,J))*FSM(I,J)
C
      DO 48 I=1,IM
      COVRHS(I)=.1E0*SQRT(GRAV/H(I,1))
   48 COVRHN(I)=.1E0*SQRT(GRAV/H(I,JMM2))
      DO 49 J=1,JM
      COVRHW(J)=.1E0*SQRT(GRAV/H(1,J))
   49 COVRHE(J)=.1E0*SQRT(GRAV/H(IMM1,J))
      DO 50 J=1,JM
      DO 50 I=1,IM
      L(I,J,1)=0.E0
   50 L(I,J,KB)=0.E0
c
      do j=1,jm
         uabe(j)=UAB(IM,J)
         uabw(j)=UAB(2,J)
         elbe(j)=ELB(IM,J)
         elbw(j)=ELB(2,J)
      end do
      do i=1,im
         vabn(i)=vab(i,jm)
         vabs(i)=vab(i,2)
         elbn(i)=elb(i,jm)
         elbs(i)=elb(i,2)
      end do
c
c     print *,'Western bound UABW i=2,j=1,jm jm=',jm
c     write(6,201) (uabw(j),j=1,jm)
c     print *,'Southern bound VABS i=1,im,j=2 im=',im
c     write(6,201) (vabs(i),i=1,im)
c     print *,'Eastern bound UABE i=im,j=1,jm jm=',jm
c     write(6,201) (uabe(j),j=1,jm)
c     print *,'Nothern bound VABN i=1,im,j=jm im=',im
c     write(6,201) (vabn(i),i=1,im)
c     print *,'Western bound ELBW i=2,j=1,jm jm=',jm
c     write(6,201) (elbw(j),j=1,jm)
c     print *,'Southern bound ELBS i=1,im,j=2 im=',im
c     write(6,201) (elbs(i),i=1,im)
c     print *,'Eastern bound ELBE i=im,j=1,jm jm=',jm
c     write(6,201) (elbe(j),j=1,jm)
c     print *,'Nothern bound ELBN i=1,im,j=jm im=',im
c     write(6,201) (elbn(i),i=1,im)
 201  format(10f7.3)
c    
      DO 51 J=1,JM
      DO 51 I=1,IM
      UA(I,J)=0.0
      VA(I,J)=0.0
      UAB(I,J)=0.0
      VAB(I,J)=0.0
      ELB(I,J)=0.0
      ETB(I,J)=0.0
      EL(I,J)=ELB(I,J)
      ET(I,J)=ETB(I,J)
      ETF(I,J)=ET(I,J)
      D(I,J)=H(I,J)+EL(I,J)
   51 DT(I,J)=H(I,J)+ET(I,J)
CCCC      print*,'uabs = ',uabs
c     print*,'elbn = ',elbn
c     print*,'vabn = ',vabn
c     print*,'elbs = ',elbs
      DO 52 K=1,KB
      DO 52 J=1,JM
      DO 52 I=1,IM
      L(I,J,K)=Q2LB(I,J,K)/(Q2B(I,J,K)+SMALL)
      KQ(I,J,K)=0.2E0*L(I,J,K)*SQRT(Q2B(I,J,K))
      Q2(I,J,K)=Q2B(I,J,K)
      Q2L(I,J,K)=Q2LB(I,J,K)
      T(I,J,K)=TB(I,J,K)
 52   S(I,J,K)=SB(I,J,K)
c
	do j=1,jm
	do i=1,im
	   T(i,j,kb)=T(i,j,kb-1)
	   S(i,j,kb)=S(i,j,kb-1)
	   TB(i,j,kb)=TB(i,j,kb-1)
	   SB(i,j,kb)=SB(i,j,kb-1)
	end do
       end do
c
      call oadjust(0)
c
      do j=1,jm
         do i=1,im
            TSURF(i,j)=T(i,j,1)
            SSURF(i,j)=S(i,j,1)
            wtsurf(i,j)=T(i,j,1)
         end do
      end do
c
      DO K=1,KB
        DO J=1,JM
          DO I=1,IM
c--------- falk 12/13/00 put ub=0,vb=0
c--------- they were calculated for elev=0 in CLIMAT2MODEL
c           UB(I,J,K)=0.
c           VB(I,J,K)=0.
c---------
            TBW(j,k)=T(1,j,k)
            TBE(j,k)=T(im,j,k)
            TBN(i,k)=T(i,jm,k)
            TBS(i,k)=T(i,1,k)
            S(I,J,K)=SB(I,J,K)
            SBW(j,k)=S(1,j,k)
            SBE(j,k)=S(im,j,k)
            SBN(i,k)=S(i,jm,k)
            SBS(i,k)=S(i,1,k)
          end do
        end do
      end do
c--------- falk 04-15-02 
c     print *,'(TBE(174,k),k=1,kb) east bnd GS center after MIXSST'
c     write(6,201) (TBE(174,k),k=1,kb)
c---------
      DO K=1,KB
        DO J=1,JM
          DO I=1,IM
            U(I,J,K)=UB(I,J,K)
            V(I,J,K)=VB(I,J,K)
          end do
        end do
      end do
Cfr****************************************************************
Cfr  MEMORISING OF INITIAL TEMPERATURE FIELD
      DO 33 K=1,KB
      DO 33 J=1,JM
      DO 33 I=1,IM
      TBIN(I,J,K)=TB(I,J,K)*FSM(I,J)
      SBIN(I,J,K)=SB(I,J,K)*FSM(I,J)
  33  continue
c
Cfr****************************************************************
C
      CALL BAROPG(DRHOX,DRHOY)
C
C  THE FOLLOWING DATA ARE NEEDED FOR A SEAMLESS RESTART.
C
c-------- falk 07-06-01 use description of file in script
      IF(NREAD.GT.0) THEN
C       OPEN(70,FILE=FRSTI,STATUS='OLD',FORM='UNFORMATTED')
        DO N=1,NREAD
c         READ(70) TIME0,
          READ(14) TIME0,
     1    WUBOT,WVBOT,AAM2D,UA,UAB,VA,VAB,EL,ELB,ET,ETB,EGB,
C     2 UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,RHO,ADVUU,ADVVV,ADVUA,ADVVA,
     2    UTB,VTB,U,UB,W,V,VB,T,TB,S,SB,ADVUU,ADVVV,ADVUA,ADVVA,
     3    KM,KH,KQ,L,Q2,Q2B,AAM,Q2L,Q2LB
        END DO
C       CLOSE(10)
        WRITE(6,*) TIME0*24.0, ' <= ',FRSTI
C
        if(sstsource.eq.1) then
	  do j=1,jm
	    do i=1,im
	      T(i,j,kb)=T(i,j,kb-1)
	      S(i,j,kb)=S(i,j,kb-1)
	      TB(i,j,kb)=TB(i,j,kb-1)
	      SB(i,j,kb)=SB(i,j,kb-1)
	    end do
          end do
c
c--- falk 06-19-01 comment call to 2 subr.
c---      use MIXSSTZ after reading sharpened T,S (read(13))
c
c         CALL SERFTEMPR(IM1,JM1)
c         CALL MIXSST
c
          do j=1,jm
            do i=1,im
              TSURF(i,j)=T(i,j,1)
              SSURF(i,j)=S(i,j,1)
              wtsurf(i,j)=T(i,j,1)
            end do
          end do
c
        end if
        CALL DENS(S,T,RHO)
      ENDIF

      call date2day(year,julday,startdate)
      print*,'startdate=',startdate,';  julday=',julday
      julday=julday+ihours/24.
      call day2date(year,julday,enddate)
      WRITE(FRSTO,'(''RST.'',I8.8)') enddate
      WRITE(6,'(/'' * New ReStart File Name : '',A/)') FRSTO
C
      DO 81 J=1,JM
      DO 81 I=1,IM
      D(I,J)=H(I,J)+EL(I,J)
   81 DT(I,J)=H(I,J)+ET(I,J)
c
      do j=1,jm
        do i=1,im
          SSTIN(i,j)=T(i,j,1)
        end do
      end do
c
      TIME0=0
      TIME=TIME0
      CALL OUTPUT(startdate)
      print *,' OCEANINIT end'
c
      return
      end

      SUBROUTINE OCEANSTEP(icoupling)
C----------------------------------------------------------------------
C     THIS SUBROUTINE DOES ONE STEP IN TIME OF THE OCEAN MODEL
C----------------------------------------------------------------------
C
      INCLUDE 'comblk.h'
      INCLUDE 'comblk1.h'
c
      TIME=DAYI*DTI*FLOAT(IINT)+TIME0
      RAMP=1.E0
C
      IF(ICOUPLING.EQ.1) THEN
        call ATMOS2OCEAN
      ELSE
        do j=2,jmm1
          do i=2,imm1
            wtsurf(i,j)=0.e0
            swrad(i,j)=0.e0
          end do
        end do
c
        do j=1,jm
          do i=1,im
            taux(i,j)=0.e0
            tauy(i,j)=0.e0
            wusurf(i,j)=0.e0
            wvsurf(i,j)=0.e0
          enddo
        enddo
c
        if((time*24.0).lt.float(inowindh)) then
          call wind(pathname,startdate)
          if(pathname1.ne.pathname) then
            call wind(pathname1,startdate)
          end if
c Changes Biju Thomas to slowly increase WUSURF and WVSURF on 03/24/2003
          WSURFC=1.
c         if(time*24..lt.6.) WSURFC=1.-(6.-time*24.)/6.
          if(time*24..lt.12.) WSURFC=1.-(12.-time*24.)/12.
c          print*,'watch',time*24.,WSURFC
          do i=1,im
          do j=1,jm
c          print*,'watch',time*24.,WSURFC,WUSURF(i,j),WSURFC*WUSURF(i,j)
            WUSURF(i,j)=WSURFC*WUSURF(i,j)
            WVSURF(i,j)=WSURFC*WVSURF(i,j)
          end do
          end do
c  -------------------------------03/24/2003----------------------------
        else
          do j=1,jm
            do i=1,im
              taux(i,j)=0.e0
              tauy(i,j)=0.e0
              wusurf(i,j)=0.e0
              wvsurf(i,j)=0.e0
            enddo
          enddo
        endif
c
        do j=1,jm
          do i=1,im
            if(taux(i,j).gt.0.1.and.nbc2d(i,j).eq.3) nbc2d(i,j)=1
          end do
        end do
      END IF
c
      IF(MODE.NE.2) THEN
        CALL ADVCT(ADVX,ADVY)
        CALL BAROPG(DRHOX,DRHOY)
C**********************************************************************
C     HOR VISC = HORCON*DX*DY*SQRT((DU/DX)**2+(DV/DY)**2
C                                 +.5*(DU/DY+DV/DX)**2)
C**********************************************************************
C  If MODE.EQ.2 then initial values of AAM2D are used. If one wishes 
C  Smagorinsky lateral viscosity and diffusion for an external mode 
C  calculation, then appropiate code can be adapted from that below
C  and installed after s.n 102 and before s.n. 5000 in subroutine ADVAVE.
        DO 95 K=1,KBM1
        DO 95 J=2,JMM1
        DO 95 I=2,IMM1
        AAM(I,J,K)=HORCON*DX(I,J)*DY(I,J)
     1          *SQRT( ((U(I+1,J,K)-U(I,J,K))/DX(I,J))**2
     2                 +((V(I,J+1,K)-V(I,J,K))/DY(I,J))**2
     3 +.5E0*(.25E0*(U(I,J+1,K)+U(I+1,J+1,K)-U(I,J-1,K)-U(I+1,J-1,K))
     4              /DY(I,J)
     5       +.25E0*(V(I+1,J,K)+V(I+1,J+1,K)-V(I-1,J,K)-V(I-1,J+1,K))
     6              /DX(I,J)) **2)
   95   CONTINUE
C --  Form vertical averages of 3-D fields for use in external mode --
        DO 96 J=1,JM
        DO 96 I=1,IM
        ADVUU(I,J)=0.
        ADVVV(I,J)=0.
        TRNU(I,J)=0.
        TRNV(I,J)=0.
        AAM2D(I,J)=0.
   96   CONTINUE    
        DO 100 K=1,KBM1
        DO 100 J=1,JM
        DO 100 I=1,IM
        ADVUU(I,J)=ADVUU(I,J)+ADVX(I,J,K)*DZ(K)
        ADVVV(I,J)=ADVVV(I,J)+ADVY(I,J,K)*DZ(K)
        TRNU(I,J)=TRNU(I,J)+DRHOX(I,J,K)*DZ(K)
        TRNV(I,J)=TRNV(I,J)+DRHOY(I,J,K)*DZ(K)
        AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
  100   CONTINUE
C ----------------------------------------------------------------------
        CALL ADVAVE(ADVUA,ADVVA,MODE)
        DO 87 J=1,JM
        DO 87 I=1,IM
        ADVUU(I,J)=ADVUU(I,J)-ADVUA(I,J)
   87   ADVVV(I,J)=ADVVV(I,J)-ADVVA(I,J)
C
      ENDIF
C------------------------------------------------------------------------
      DO 399 J=1,JM
      DO 399 I=1,IM
  399 EGF(I,J)=EL(I,J)*ISPI
      DO 400 J=2,JM
      DO 400 I=2,IM
      UTF(I,J)=UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
 400  VTF(I,J)=VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
C********** BEGIN EXTERNAL MODE ***************************************
                      DO 8000 IEXT=1,ISPLIT
      DO 405 J=2,JM
      DO 405 I=2,IM
      FLUXUA(I,J)=.25E0*(D(I,J)+D(I-1,J))*(DY(I,J)+DY(I-1,J))*UA(I,J)
 405  FLUXVA(I,J)=.25E0*(D(I,J)+D(I,J-1))*(DX(I,J)+DX(I,J-1))*VA(I,J)
C
      DO 410 J=2,JMM1
      DO 410 I=2,IMM1
  410 ELF(I,J)=ELB(I,J)
     1    -DTE2*(FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))
     2                    /ART(I,J)
C
      CALL BCOND(1)
C
      IF(MOD(IEXT,ISPADV).EQ.0) CALL ADVAVE(ADVUA,ADVVA,MODE)
C  Note that ALPHA = 0. is perfectly acceptable. The value, ALPHA = .225
C  permits a longer time step.
      ALPHA=0.225     
      DO 420 J=2,JMM1
      DO 420 I=2,IMM1
  420 UAF(I,J)=ADVUU(I,J)+ADVUA(I,J)   
     1    -ARU(I,J)*.25*(  COR(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     2              +COR(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     3         +.25E0*GRAV*(DY(I,J)+DY(I-1,J))*(D(I,J)+D(I-1,J))
     4             *( (1.E0-2.E0*ALPHA)*(EL(I,J)-EL(I-1,J))
     4            +ALPHA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
     5               +TRNU(I,J)   
     6      +ARU(I,J)*( WUSURF(I,J)-WUBOT(I,J)   )
      DO 425 J=2,JMM1
      DO 425 I=2,IMM1
  425 UAF(I,J)=
     1         ((H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J))*ARU(I,J)*UAB(I,J)
     2                -4.E0*DTE*UAF(I,J))
     3        /((H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))*ARU(I,J))
      DO 430 J=2,JMM1
      DO 430 I=2,IMM1
 430  VAF(I,J)=ADVVV(I,J)+ADVVA(I,J)  
     1    +ARV(I,J)*.25*(  COR(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     2               +COR(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     3         +.25E0*GRAV*(DX(I,J)+DX(I,J-1))*(D(I,J)+D(I,J-1))
     4             *( (1.E0-2.E0*ALPHA)*(EL(I,J)-EL(I,J-1))
     4            +ALPHA*(ELB(I,J)-ELB(I,J-1)+ELF(I,J)-ELF(I,J-1)) )
     5               +TRNV(I,J)  
     6    + ARV(I,J)*( WVSURF(I,J)-WVBOT(I,J)   )
      DO 435 J=2,JMM1
      DO 435 I=2,IMM1
  435 VAF(I,J)=
     1        ((H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1))*VAB(I,J)*ARV(I,J)
     2              -4.E0*DTE*VAF(I,J))
     3       /((H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))*ARV(I,J))
      CALL BCOND(2)
C
      IF(IEXT.LT.(ISPLIT-2)) GO TO 440
      IF(IEXT.EQ.(ISPLIT-2))THEN
        DO 431 J=1,JM
        DO 431 I=1,IM
  431   ETF(I,J)=.25*SMOTH*ELF(I,J)
       ENDIF
      IF(IEXT.EQ.(ISPLIT-1)) THEN
        DO 432 J=1,JM
        DO 432 I=1,IM
  432   ETF(I,J)=ETF(I,J)+.5*(1.-.5*SMOTH)*ELF(I,J)
       ENDIF
      IF(IEXT.EQ.(ISPLIT-0)) THEN
        DO 433 J=1,JM
        DO 433 I=1,IM
  433   ETF(I,J)=(ETF(I,J)+.5*ELF(I,J))*FSM(I,J)
       ENDIF
 440  CONTINUE
C
C  TEST FOR CFL VIOLATION. IF SO, PRINT AND STOP
Cfu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      UAMAX = 0.E0
      VAMAX = 0.E0
      DO J = 1, JM
        DO I = 1, IM
          IF(UAMAX.LT.ABS(UAF(I,J))) THEN
            UAMAX = ABS(UAF(I,J))
            IUMAX = I
            JUMAX = J
          ENDIF
          IF(VAMAX.LT.ABS(VAF(I,J))) THEN
            VAMAX = ABS(VAF(I,J))
            IVMAX = I
            JVMAX = J
          ENDIF
        END DO
      END DO
      UVAMAX = AMAX1(UAMAX,VAMAX)
      IF(UVAMAX.GT.100.E0) GO TO 9001
Cfu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Cfu      VAMAX=0.       
Cfu      DO 442 J=1,JM
Cfu      DO 442 I=1,IM
Cfu      IF(ABS(VAF(I,J)).GE.VAMAX)VAMAX=ABS(VAF(I,J))   
Cfu  442 CONTINUE
Cfu      IF(VAMAX.GT.100.) GO TO 9001
C    
C     IF(MOD(IINT,IPRINT).EQ.0..AND.IEXT.EQ.ISPLIT/2)   
C    1    CALL VORT(ADVUA,ADVVA,ADVUU,ADVVV,TRNU,TRNV,IM,JM)  
C
Cfu Boundary-Damping of ELF, UAF, VAF
      IF(IDAMP.EQ.1) CALL BDAMP
C       APPLY FILTER TO REMOVE TIME SPLIT. RESET TIME SEQUENCE.
      DO 445 J=1,JM
      DO 445 I=1,IM
      UA(I,J)=UA(I,J)+.5E0*SMOTH*(UAB(I,J)-2.E0*UA(I,J)+UAF(I,J))
      VA(I,J)=VA(I,J)+.5E0*SMOTH*(VAB(I,J)-2.E0*VA(I,J)+VAF(I,J))
      EL(I,J)=EL(I,J)+.5E0*SMOTH*(ELB(I,J)-2.E0*EL(I,J)+ELF(I,J))
      ELB(I,J)=EL(I,J)
      EL(I,J)=ELF(I,J)
      D(I,J)=H(I,J)+EL(I,J)
      UAB(I,J)=UA(I,J)
      UA(I,J)=UAF(I,J)
      VAB(I,J)=VA(I,J)
      VA(I,J)=VAF(I,J)
  445 CONTINUE

c      kbb=1
c      kbb1=jm-70
c      i=monitor(im,jm,kbb,EL,'Surface elevation',
c     *' flat')
C
      IF(IEXT.EQ.ISPLIT) GO TO 8000
      DO 450 J=2,JM
      DO 450 I=2,IM
      EGF(I,J)=EGF(I,J)+EL(I,J)*ISPI
      UTF(I,J)=UTF(I,J)+UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
 450  VTF(I,J)=VTF(I,J)+VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
C
C
 8000                    CONTINUE
C---------------------------------------------------------------------
C          END EXTERNAL (2-D) MODE CALCULATION
C     AND CONTINUE WITH INTERNAL (3-D) MODE CALCULATION
C---------------------------------------------------------------------
      IF(IINT.EQ.1.AND.TIME0.EQ.0.) GO TO 8200
      IF(MODE.EQ.2) GO TO 8200
C---------------------------------------------------------------------
C      ADJUST U(Z) AND V(Z) SUCH THAT
C      VERTICAL AVERAGE OF (U,V) = (UA,VA)
C---------------------------------------------------------------------
C
      DO 299 J=1,JM
      DO 299 I=1,IM
  299 TPS(I,J)=0.E0
      DO 300 K=1,KBM1
      DO 300 J=2,JMM1
      DO 300 I=2,IMM1
 300  TPS(I,J)=TPS(I,J)+U(I,J,K)*DZ(K)
      DO 302 K=1,KBM1
      DO 302 J=2,JMM1
      DO 302 I=2,IMM1
 302  U(I,J,K)=(U(I,J,K)-TPS(I,J))+
     1     (UTB(I,J)+UTF(I,J))/(DT(I,J)+DT(I-1,J))
      DO 303 J=1,JM
      DO 303 I=1,IM
 303  TPS(I,J)=0.
      DO 304 K=1,KBM1
      DO 304 J=2,JMM1
      DO 304 I=2,IMM1
 304  TPS(I,J)=TPS(I,J)+V(I,J,K)*DZ(K)
      DO 306 K=1,KBM1
      DO 306 J=2,JMM1
      DO 306 I=2,IMM1
 306  V(I,J,K)=(V(I,J,K)-TPS(I,J))+
     2     (VTB(I,J)+VTF(I,J))/(DT(I,J)+DT(I,J-1))
C
C----------------------------------------------------------------
C     VERTVL INPUT = U,V,DT(=H+ET),ETF,ETB; OUTPUT = W
C----------------------------------------------------------------
      CALL VERTVL(DTI2)
      CALL BCOND(5)
C
C
      DO 307 K=1,KB
      DO 307 J=1,JM
      DO 307 I=1,IM
      UF(I,J,K)=0.E0
  307 VF(I,J,K)=0.E0
C
      print*,'IMAY=',IMAY
      if(IMAY.EQ.1)THEN
C
C----------------------------------------------------------------
C     COMPUTE Q2F AN Q2LF USING UF AND VF AS TEMPORARY VARIABLES
C----------------------------------------------------------------
      CALL ADVQ(Q2B,Q2,DTI2,UF)
      CALL ADVQ(Q2LB,Q2L,DTI2,VF)
      CALL PROFQ(DTI2)
      CALL BCOND(6)
      DO 310 K=1,KB
      DO 310 J=1,JM
      DO 310 I=1,IM
      Q2(I,J,K)=Q2(I,J,K)
     1    +.5*SMOTH*(UF(I,J,K)+Q2B(I,J,K)-2.*Q2(I,J,K))
      Q2L(I,J,K)=Q2L(I,J,K)
     1    +.5*SMOTH*(VF(I,J,K)+Q2LB(I,J,K)-2.*Q2L(I,J,K))
      Q2B(I,J,K)=Q2(I,J,K)
      Q2(I,J,K)=UF(I,J,K)
      Q2LB(I,J,K)=Q2L(I,J,K)
      Q2L(I,J,K)=VF(I,J,K)
310   CONTINUE
C
       ELSE
       call rikh
       ENDIF
c
C----------------------------------------------------------------
C     COMPUTE TF AN SF USING UF AND VF AS TEMPORARY VARIABLES
C----------------------------------------------------------------
      IF(MODE.EQ.4) GO TO 360
C  In this calculation, S = constant. However, we calculate
C  S to demonstrate that the same subroutines can be used to 
C  calculate multiple scalars.
      CALL ADVT(TB,T,TMEAN,DTI2,UF)
      CALL ADVT(SB,S,SMEAN,DTI2,VF)
      CALL PROFT(UF,WTSURF,SWRAD,TSURF,NBC2D,DTI2)
      CALL PROFT(VF,WSSURF,SWRAD,SSURF,NBC2Ds,DTI2)
      CALL BCOND(4)
C
Cfu-------------------------------------------------
Cfr*******************************************************
Cfr  ELIMINATING TEMPERATURE INCREASE
Cfr      DO 352 K=1,KB
Cfr      DO 352 J=1,JM
Cfr      DO 352 I=1,IM
Cfr      IF(UF(I,J,K).GT.19.01) THEN
Cfr      UF(I,J,K)=19
Cfr      end if
Cfr  352 continue
Cfr*******************************************************
      IF((AMOD(FLOAT(IINT),SMH).EQ.0.0).AND.(ISMOTH.EQ.1)) THEN
      PRINT*,'DO SMOOTHING OF T,S AT TIME= ',TIME*24.0

      DO 354 K=1,KB
      DO 354 J=1,JM
      DO 354 I=1,IM
  354 UF(I,J,K)=UF(I,J,K)-TBIN(I,J,K)

      DO K=1,KB
        CALL SMOOTHING(UF(1,1,K),IM,JM,0.25)
        CALL SMOOTHING(UF(1,1,K),IM,JM,-0.27)
        CALL SMOOTHING(VF(1,1,K),IM,JM,0.25)
        CALL SMOOTHING(VF(1,1,K),IM,JM,-0.27)
      ENDDO

      DO 353 K=1,KB
      DO 353 J=1,JM
      DO 353 I=1,IM
  353 UF(I,J,K)=UF(I,J,K)+TBIN(I,J,K)

      ENDIF
Cfu------------------------------------------------
C
Cfr*************************************************************
Cfr   ELIMINATING TEMPERATURE GROWTH
Cfr*************************************************************
      DO 358 J=2,JM-1
      DO 358 I=2,IM-1
      IF(UF(I,J,1).GT.T(I,J,1)) THEN
      IF(UF(I,J,1).GE.T(I-1,J,1).AND.
     1   UF(I,J,1).GE.T(I+1,J,1).AND.
     2   UF(I,J,1).GE.T(I,J-1,1).AND.
     3   UF(I,J,1).GE.T(I,J+1,1).AND.
C
     4   UF(I,J,1).GE.T(I,J,2)) UF(I,J,1)=T(I,J,1)
      END IF
  358 CONTINUE


      DO 359 K=2,KB-1
      DO 359 J=2,JM-1
      DO 359 I=2,IM-1
      IF(UF(I,J,K).GT.T(I,J,K)) THEN
      IF(UF(I,J,K).GE.T(I,J,K-1).AND.
     1   UF(I,J,K).GE.T(I-1,J,K).AND.
     2   UF(I,J,K).GE.T(I+1,J,K).AND.
     3   UF(I,J,K).GE.T(I,J-1,K).AND.
     4   UF(I,J,K).GE.T(I,J+1,K).AND.
C
     5   UF(I,J,K).GE.T(I,J,K+1))  UF(I,J,K)=T(I,J,K)
      END IF
  359 CONTINUE

      DO 357 J=2,JM-1
      DO 357 I=2,IM-1
      IF(UF(I,J,KB).GT.T(I,J,KB)) THEN
      IF(UF(I,J,KB).GE.T(I-1,J,KB).AND.
     1   UF(I,J,KB).GE.T(I+1,J,KB).AND.
     2   UF(I,J,KB).GE.T(I,J-1,KB).AND.
     3   UF(I,J,KB).GE.T(I,J+1,KB).AND.
C
     4   UF(I,J,KB).GE.T(I,J,KB-1)) UF(I,J,KB)=T(I,J,KB)
      END IF
  357 CONTINUE
c
c-----------------  NUDGING  ----------------------------------
c
c      CALL NUDGE(2,wndg)
c
Cfr*************************************************************

      DO 355 K=1,KB
      DO 355 J=1,JM
      DO 355 I=1,IM
      T(I,J,K)=T(I,J,K)+.5*SMOTH*(UF(I,J,K)
     1  +TB(I,J,K)-2.E0*T(I,J,K))
      S(I,J,K)=S(I,J,K)+.5*SMOTH*(VF(I,J,K)
     1  +SB(I,J,K)-2.E0*S(I,J,K))
      TB(I,J,K)=T(I,J,K)
      T(I,J,K)=UF(I,J,K)
      SB(I,J,K)=S(I,J,K)
  355 S(I,J,K)=VF(I,J,K)
c      SB(I,J,K)=SB(I,J,K)
c  355 S(I,J,K)=S(I,J,K)
      print*,iint

c      i=monitor(im,jm,kb,T,'Temperature',' flat')
Cfr*************************************************************
Cfr   ELIMINATING TEMPERATURE GROWTH
Cfr*************************************************************
C      DO 358 J=2,JM-1
C      DO 358 I=2,IM-1
C      IF(T(I,J,1).GT.TB(I,J,1)) THEN
C      IF(T(I,J,1).GE.TB(I-1,J,1).AND.
C     1   T(I,J,1).GE.TB(I+1,J,1).AND.
C     2   T(I,J,1).GE.TB(I,J-1,1).AND.
C     3   T(I,J,1).GE.TB(I,J+1,1)) T(I,J,1)=TB(I,J,1)
C      END IF
C  358 CONTINUE
C
C      DO 359 K=2,KB
C      DO 359 J=2,JM-1
C      DO 359 I=2,IM-1
C      IF(T(I,J,K).GT.TB(I,J,K)) THEN
C      IF(T(I,J,K).GE.TB(I,J,K-1).AND.
C     1   T(I,J,K).GE.TB(I-1,J,K).AND.
C     2   T(I,J,K).GE.TB(I+1,J,K).AND.
C     3   T(I,J,K).GE.TB(I,J-1,K).AND.
C     4   T(I,J,K).GE.TB(I,J+1,K)) T(I,J,K)=TB(I,J,K)
C      END IF
C  359 CONTINUE
Cfr*************************************************************
C
      CALL DENS(S,T,RHO)
C
  360 CONTINUE
C
C----------------------------------------------------------------
C     COMPUTE UF AND VF
C----------------------------------------------------------------
      CALL ADVU(DRHOX,ADVX,DTI2)
      CALL ADVV(DRHOY,ADVY,DTI2)
      CALL PROFU(DTI2)
      CALL PROFV(DTI2)
      CALL BCOND(3)
C
Cfu-------------------------------------------------
      IF((AMOD(float(IINT),SMH*3.).EQ.0.0).AND.(ISMOTH.EQ.1)) THEN
        PRINT*,'DO SMOOTHING OF U,V AT TIME= ',TIME*24.0
        DO K=1,KBM1
          CALL SMOOTHING(UF(1,1,K),IM,JM,0.25)
          CALL SMOOTHING(UF(1,1,K),IM,JM,-0.25)
          CALL SMOOTHING(VF(1,1,K),IM,JM,0.25)
          CALL SMOOTHING(VF(1,1,K),IM,JM,-0.25)
        ENDDO
      ENDIF
Cfu-------------------------------------------------
C
      DO 369 J=1,JM
      DO 369 I=1,IM
  369 TPS(I,J)=0.E0
      DO 370 K=1,KBM1
      DO 370 J=1,JM
      DO 370 I=1,IM
  370 TPS(I,J)=TPS(I,J)+(UF(I,J,K)+UB(I,J,K)-2.E0*U(I,J,K))*DZ(K)
      DO 372 K=1,KBM1
      DO 372 J=1,JM
      DO 372 I=1,IM
  372 U(I,J,K)=U(I,J,K)+.5*SMOTH*(UF(I,J,K)+UB(I,J,K)-2.E0*U(I,J,K)
     1        -TPS(I,J))
      DO 3721 J=1,JM
      DO 3721 I=1,IM
 3721 TPS(I,J)=0.E0
      DO 374 K=1,KBM1
      DO 374 J=1,JM
      DO 374 I=1,IM
  374 TPS(I,J)=TPS(I,J)+(VF(I,J,K)+VB(I,J,K)-2.E0*V(I,J,K))*DZ(K)
      DO 376 K=1,KBM1
      DO 376 J=1,JM
      DO 376 I=1,IM
  376 V(I,J,K)=V(I,J,K)+.5*SMOTH*(VF(I,J,K)+VB(I,J,K)-2.E0*V(I,J,K)
     1        -TPS(I,J))
      DO 377 K=1,KB
      DO 377 J=1,JM
      DO 377 I=1,IM
      UB(I,J,K)=U(I,J,K)
      U(I,J,K)=UF(I,J,K)
      VB(I,J,K)=V(I,J,K)
  377 V(I,J,K)=VF(I,J,K)
C
 8200 CONTINUE
C
      DO 8210 J=1,JM
      DO 8210 I=1,IM
      EGB(I,J)=EGF(I,J)
      ETB(I,J)=ET(I,J)
      ET(I,J)=ETF(I,J)
      DT(I,J)=H(I,J)+ET(I,J)
      UTB(I,J)=UTF(I,J)
 8210 VTB(I,J)=VTF(I,J)

C
C------------- CREATING SST FOR THE ATMOSPHERE MODEL ------------------
C
      IF(ICOUPLING.EQ.1) THEN
        TDIFMAX=0.
        do j=1,jm
          do i=1,im
            TMA(I,J)=T(I,J,1)+TBIAS+273.15
            if(t(i,j,1).eq.0) tma(i,j)=-99.
            TDIF=ABS(T(I,J,1)-TB(I,J,1))
            IF (TDIF. GT. TDIFMAX) TDIFMAX=TDIF
          end do
        end do
        PRINT*, 'OCEAN MODEL:  MAX SST DECREASE ', TDIFMAX 
c--- passing sst from the ocean model to the atmosphere model
C
        CALL INTERP2(TMA,SST,IM,JM,IMA,JMA,XST,YST,XLON,YLAT)
      END IF
c
C--------------------------------------------------------------------
C           BEGIN PRINT SECTION
C--------------------------------------------------------------------
c
      IF(MOD(IINT,IPRINT).NE.0.) GO TO 7000
C
C--- MODIFYING SURFACE ELEVATION DATA FOR PLOTTING
      do j=1,jm
        do i=1,im
          elb(i,j)=elb(i,j)+(1.-fsm(i,j))*999
        end do
      end do
C
      CALL OUTPUT(startdate)
C
        IF(IINT.GE.ISWTCH) IPRINT=IPRTD2*24*3600/INT(DTI)
        WRITE(6,'(/,'' TIME ='',F10.2,'' HOUR   IINT ='',I8,
     1  ''   IEXT ='',I8,''   IPRINT ='',I5,//)')
     2  TIME*24.E0,IINT,IEXT,IPRINT
c---------- falk 04-23-02 check 
        print *,' check before FINDPSI'
        CALL FINDPSI
c---------- falk 04-23-02 check 
        print *,' check after FINDPSI'
c---------- falk 04-23-02 comment 2 next lines
c       IF(MODE.NE.2) THEN
c       ENDIF
C
 9001 CONTINUE
C
      VTOT=0.E0
      TTOT=0.E0
      STOT=0.E0
      DO 8888 K=1,KBM1
      DO 8888 J=1,JM
      DO 8888 I=1,IM
      DVTOT=DX(I,J)*DY(I,J)*DT(I,J)*DZ(K)*FSM(I,J)
      VTOT=VTOT+DVTOT
      TTOT=TTOT+TB(I,J,K)*DVTOT
      STOT=STOT+SB(I,J,K)*DVTOT
 8888 CONTINUE
c---------- falk 04-23-02 chrck
      print *,' VTOT=',VTOT
      TTOT=TTOT/VTOT
      STOT=STOT/VTOT
      WRITE(6,'(''  VTOT ='',E20.8,''    TTOT ='',E20.8,
     1      ''   STOT ='',E20.8)') VTOT,TTOT,STOT
      IF (UVAMAX.GT.100.) THEN
      CALL OUTPUT(startdate) 
      WRITE(6,'(//////////////////)')
      WRITE(6,'(''************************************************'')')
      WRITE(6,'(''************ ABNORMAL JOB END ******************'')')
      WRITE(6,'(''************* USER TERMINATED ******************'')')
      WRITE(6,'(''************************************************'')')
      WRITE(6,'(''UAMAX= '',E12.3,'' I,J,H ='',2I4,F7.0)')
     *            UAMAX,IUMAX,JUMAX,H(IUMAX,JUMAX)
      WRITE(6,'(''VAMAX= '',E12.3,'' I,J,H ='',2I4,F7.0)')
     *            VAMAX,IVMAX,JVMAX,H(IVMAX,JVMAX)
      WRITE(6,'('' IINT,IEXT ='',2I6)') IINT,IEXT
      STOP
      ENDIF
 7000 CONTINUE
C
C----------------------- END PRINT SECTION -----------------------------
C
      return
      END
      SUBROUTINE WGHTS(IMX,JMX,INMAX,JNMAX,XO,YO,XN,YN,WTI,WTJ,IDI,IDJ)  
      DIMENSION XO(IMX),  YO(JMX),  XN(INMAX),  YN(JNMAX),  WTI(INMAX), 
     ,         WTJ(JNMAX), IDI(INMAX), IDJ(JNMAX)       
C     1ST SEARCH FOR Y PARTIAL WEIGHTS        
      JMXM1=JMX-1      
      IMXM1=IMX-1     
      DO 200 JN=1,JNMAX     
      DO 100 JO=1,JMXM1
      IF((YN(JN).GE.YO(JO)).AND.(YN(JN).LT.YO(JO+1)))  GO TO 150    
  100 CONTINUE
      GO TO 200
  150 WTJ(JN)=(YN(JN)-YO(JO))
      IDJ(JN)=JO      
  200 CONTINUE      
C     NOW SEARCH FOR X PARTIAL WEIGHTS    
      DO 400 IN=1,INMAX     
      DO 300 IO=1,IMXM1       
      IF((XN(IN).GE.XO(IO)).AND.(XN(IN).LT.XO(IO+1)))  GO TO 350  
  300 CONTINUE
       GO TO 400
  350 WTI(IN)=(XN(IN)-XO(IO))   
      IDI(IN)=IO      
  400 CONTINUE    
      RETURN 
      END  
c
c-------------------
c
      SUBROUTINE TRANSPORTS
c---------------------------------------------------------------
C    THIS SUBROUTINE CORRECTS THE BOUNDARY VELOCITIES 
C    BY MATCHING THE TRANSPORTS TO THE SPECIFIED VALUES
c    The integration along the boundary is done counterclockwise:
c    western boundary from north to south - southern boundary -
c    eastern boundary - northern boundary
c    The boundary integral of the added barotropic velocity usually
c    is nonzero, i.e., the surface elevation jump is created in the
c    north-west corner of the domain. This jump is being eliminated
c    by adding variable velocity correlated with depth. These additions
c    are chosen to have minimum possible relative velocity deviation.
c    The additional velocity has the form of alpha*D(i,j)/1000+const.
C---------------------------------------------------------------
      INCLUDE 'comblk.h'
      INCLUDE 'comblk1.h'
c
c-------------- falk 12/06/00 change west boundary
      parameter(nrs=1,nrn=1,nrw=1,nre=7,
     *          je1=24,je2=48,je3=174,je4=184,je5=187,je6=198)
c---------- j1=24(13.85N) j2=48(17.87) j3=174(38.96) j4=184(40.63)
c---------- j5=187(41.14) j6=198(42.98)
c     parameter(nrs=1,nrn=1,nrw=3,nre=5)
      real level1,level2,levnm,level
      real zi(kb-1),zi1(kb-1),tmpsi(kb-1),tmpti(kb-1)
      real trnsw(nrw),trnsw_inc(nrw)
      real trnss(nrs),trnss_inc(nrs)
      real trnse(nre),trnse_inc(nre)
      real trnsn(nrn),trnsn_inc(nrn)
      real secarw(nrw),secar2w(nrw),secvelw(nrw),constw(nrw)
      real secars(nrs),secar2s(nrs),secvels(nrs),consts(nrs)
      real secare(nre),secar2e(nre),secvele(nre),conste(nre)
      real secarn(nrn),secar2n(nrn),secveln(nrn),constn(nrn)
      real ptsw(nrw),ptss(nrs),ptse(nre),ptsn(nrn)
      real wgcf(nrw+nrs+nre+nrn),mxd(nrw+nrs+nre+nrn)
      real alpha(nrw+nrs+nre+nrn)
      integer indw(nrw+1),inds(nrs+1),inde(nre+1),indn(nrn+1)
      real secar,const,dEL1,dEL2
      common/indexes/ind1,ind2
      PI=3.1415927
c
c---------- Specifying transports ------------------------------
c
c----- Southern boundary:
      data trnss/0.0e6/,ptss/0.0/
c----- Northern boundary:
      data trnsn/0.0e6/,ptsn/0.0/
c----- Western boundary:
c------------ 01-09-02 change western bnd for united1 domain
      data trnsw/0.0e6/,ptsw/0.0/
c----- Eastern boundary:
c------------ 01-09-02 change eastern bnd for united1 domain
      data trnse/0.e6,-60.0e6,0.0e6,90.0e6,0.0e6,-30.0e6,0.0e6/,
     *      ptse/13.85,17.85,38.96,40.63,41.14,43.0,0.0/
c
c----------  calculating vectors of indexes  -------------------
c
      inds(1)=1
      if(nrs.gt.1) then
      do n=2,nrs
        inds(n)=nint((ptss(n-1)-longmin)/(longmax-longmin)*(im-1))+1
        alonind=longmin+(inds(n)-1)*(longmax-longmin)/(im-1)
c       write(6,202) n,alonind
 202    format('             n,alonind=',i7,f7.2)
      end do
      end if
      inds(nrs+1)=im
c
      indn(1)=1
      if(nrn.gt.1) then
      do n=2,nrn
        indn(n)=nint((ptsn(n-1)-longmin)/(longmax-longmin)*(im-1))+1
        alonind=longmin+(indn(n)-1)*(longmax-longmin)/(im-1)
c       write(6,202) n,alonind
      end do
      end if
      indn(nrn+1)=im
c
      indw(1)=1
      if(nrw.gt.1) then
      do n=2,nrw
        indw(n)=nint((ptsw(n-1)-latmin)/(latmax-latmin)*(jm-1))+1
        alatind=latmin+(indw(n)-1)*(latmax-latmin)/(jm-1)
c       write(6,203) n,alatind
 203    format('             n,alatind=',i7,f7.2)
      end do
      end if
      indw(nrw+1)=jm
c
      inde(1)=1
c     print *,'     eastern bnd'
      do n=2,nre
        inde(n)=nint((ptse(n-1)-latmin)/(latmax-latmin)*(jm-1))+1
        alatind=latmin+(inde(n)-1)*(latmax-latmin)/(jm-1)
c       write(6,203) n,alatind
      end do
c
      inde(nre+1)=jm
      print *,'   '
      print *,' TRANSPORTS, check indexes'
      print *,' inds=',inds
      print *,' indn=',indn
      print *,' indw=',indw
      print *,' inde=',inde
c---------- put el=0.
      do j=1,jm
        do i=1,im
          D(i,j)=D(i,j)*dum(i,j)*dvm(i,j)-0.0001*(dum(i,j)*dvm(i,j)-1)
          el(i,j)=0.
        end do
      end do
c
c----------------  SOUTHERN  BOUNDARY  ------------------------------
c----   specify barotr. velocity as a first-guess
      do i=1,im
       VAB(i,2)=0.
      end do
c
c     print*,'SOUTHERN BOUNDARY: inds=',inds
c
c----------------  EASTERN  BOUNDARY  -------------------------------
c----   specify barotr. velocity as a first-guess
      do j=1,jm
       UAB(im,j)=0.
      end do
      do j=je1,je2
       UAB(im,j)=-sin(PI*(j-je1)/(je2-je1))
      end do
      do j=je3,je4
       UAB(im,j)=sin(PI*(j-je3)/(je4-je3))
      end do
      do j=je5,je6
       UAB(im,j)=-sin(PI*(j-je5)/(je6-je5))
      end do
c
 524  FORMAT('EASTERN BOUNDARY, REGION ',I2,':')
      do n=1,nre
        print*,'      '
        write(6,524) n
        trnse_inc(n)=0.0
        do j=inde(n)+1,inde(n+1)
          trnse_inc(n)=trnse_inc(n)+D(im,j)*DY(im,j)*UAB(im,j)
        end do
        if(abs(trnse_inc(n)).gt.1.) then
         cmp=trnse(n)/trnse_inc(n)
         PRINT*,'Multiply VELOCITY to: ',cmp,' M/S...'
c        print *,'trnse(n)=',trnse(n)
c        print *,'trnse_inc(n)=',trnse_inc(n)
         do j=inde(n)+1,inde(n+1)
          uab(im,j)=uab(im,j)*cmp
         end do
        end if
      end do
c
c     print *,'je1,je2,je3,je4,je5,je6=',je1,je2,je3,je4,je5,je6
c     print*,'inde=',inde
c---------- falk 09/28/00 print uab and el for eastern bnd
c     print *,'Eastern bound first-guess UAB after transp correction'
      do n=1,nre
       jbeg=inde(n)
       jend=inde(n+1)
c      print *,'n,jbeg,jend=',n,jbeg,jend
c      print *,'(uab(im,j),j=jbeg,jend)'
c      write(6,201) (uab(im,j),j=jbeg,jend)
 201   format(10f7.3)
      end do
c
c----------------  NORTHERN  BOUNDARY  -------------------------------
c
c----   specify barotr. velocity as a first-guess
      do i=1,im
       VAB(i,jm)=0.
      end do
c     print*,' NORTHERN BND: indn=',indn
c
c----------------  WESTERN  BOUNDARY  -------------------------------
c
c----   specify barotr. velocity as a first-guess
      do j=1,jm
       UAB(2,j)=0.
      end do
c     print*,'WESTERN BND: indw=',indw
c
      CALL BAROPG(DRHOX,DRHOY)
c--- Form vertical averages of 3-D fields for use in external mode ---
      do j=1,jm
        do i=1,im
          TRNU(i,j)=0.
          TRNV(i,j)=0.
        end do
      end do
      do k=1,KBM1
        do j=1,jm
          do i=1,im
            TRNU(i,j)=TRNU(i,j)+DRHOX(i,j,k)*DZ(k)
            TRNV(i,j)=TRNV(i,j)+DRHOY(i,j,k)*DZ(k)
          end do
        end do
      end do
c
      do j=1,jm
        do i=1,im
          D(i,j)=D(i,j)*dum(i,j)*dvm(i,j)-0.0001*(dum(i,j)*dvm(i,j)-1)
          uab(i,j)=uab(i,j)*dum(i,j)
          vab(i,j)=vab(i,j)*dvm(i,j)
        end do
      end do
c
c-----  Calculate the coefficient which relates alpha to dEL --------
c
c----------------  WESTERN  BOUNDARY  -------------------------------
c
      do n=nrw,1,-1
        secarw(n)=0.0
        secar2w(n)=0.0
        do j=indw(n)+1,indw(n+1)
          secar=.25E0*(D(2,j)+D(1,j))*(DY(2,j)+DY(1,j))
          secarw(n)=secarw(n)+secar
          secar2w(n)=secar2w(n)+secar*D(2,j)/1000
        end do
        const=-secar2w(n)/secarw(n)
        constw(n)=const
        wgcf(nrw-n+1)=0.0
        mxd(nrw-n+1)=0.0
        do j=indw(n+1),indw(n)+1,-1
          wgcf(nrw-n+1) = wgcf(nrw-n+1)
     1    +ARV(2,j)/((DX(2,j)+DX(2,j-1))*(D(2,j)+D(2,j-1)))
     2    *( COR(2,j)*D(2,j)*2*(D(2,j)/1000+const)
     3      +COR(2,j-1)*D(2,j-1)*2*(D(2,j-1)/1000+const)
     4       )/grav*dum(2,j)
          mxd(nrw-n+1) = max(mxd(nrw-n+1),D(2,j)/1000)
        end do
      end do
c
c----------------  SOUTHERN  BOUNDARY  ------------------------------
c
      do n=1,nrs
        secars(n)=0.0
        secar2s(n)=0.0
        do i=inds(n)+1,inds(n+1)
          secar=.25E0*(D(I,2)+D(I,1))*(DX(I,2)+DX(I,1))
          secars(n)=secars(n)+secar
          secar2s(n)=secar2s(n)+secar*D(i,2)/1000
        end do
        const=-secar2s(n)/secars(n)
        consts(n)=const
        wgcf(nrw+n)=0.0
        mxd(nrw+n)=0.0
        do i=inds(n)+1,inds(n+1)
          wgcf(nrw+n) = wgcf(nrw+n)
     1    +ARU(i,2)/((DY(i,2)+DY(i-1,2))*(D(i,2)+D(i-1,2)))
     2    *( COR(i,2)*D(i,2)*2*(D(i,2)/1000+const)
     3      +COR(i-1,2)*D(i-1,2)*2*(D(i-1,2)/1000+const)
     4       )/grav*dvm(i,2)
          mxd(nrw+n) = max(mxd(nrw+n),D(i,2)/1000)
         end do
      end do
c
c----------------  EASTERN  BOUNDARY  -------------------------------
c
      do n=1,nre
        secare(n)=0.0
        secar2e(n)=0.0
        do j=inde(n)+1,inde(n+1)
          secar=.25E0*(D(im,J)+D(im-1,J))*(DY(im,J)+DY(im-1,J))
          secare(n)=secare(n)+secar
          secar2e(n)=secar2e(n)+secar*D(im,j)/1000
        end do
        const=-secar2e(n)/secare(n)
        conste(n)=const
        wgcf(nrw+nrs+n)=0.0
        mxd(nrw+nrs+n)=0.0
        do j=inde(n)+1,inde(n+1)
          wgcf(nrw+nrs+n) = wgcf(nrw+nrs+n)-ARV(im-1,j)
     1    /((DX(im-1,j)+DX(im-1,j-1))*(D(im-1,j)+D(im-1,j-1)))
     2    *( COR(im-1,j)*D(im-1,j)*2*(D(im,j)/1000+const)
     3      +COR(im-1,j-1)*D(im-1,j-1)*2*(D(im,j-1)/1000+const)
     4       )/grav*dum(im,j)
          mxd(nrw+nrs+n) = max(mxd(nrw+nrs+n),D(im,j)/1000)
        end do
      end do
c
c----------------  NORTHERN  BOUNDARY  -------------------------------
c
      do n=nrn,1,-1
        secarn(n)=0.0
        secar2n(n)=0.0
        do i=indn(n)+1,indn(n+1)
          secar=.25E0*(D(I,jm)+D(I,jm-1))*(DX(I,jm)+DX(I,jm-1))
          secarn(n)=secarn(n)+secar
          secar2n(n)=secar2n(n)+secar*D(i,jm)/1000
        end do
        const=-secar2n(n)/secarn(n)
        constn(n)=const
        wgcf(nrw+nrs+nre+nrn-n+1)=0.0
        mxd(nrw+nrs+nre+nrn-n+1)=0.0
        do i=indn(n+1),indn(n)+1,-1
         wgcf(nrw+nrs+nre+nrn-n+1)=wgcf(nrw+nrs+nre+nrn-n+1)-ARU(i,jm-1)
     1    /((DY(i,jm-1)+DY(i-1,jm-1))*(D(i,jm-1)+D(i-1,jm-1)))
     2    *( COR(i,jm-1)*D(i,jm-1)*2*(D(i,jm)/1000+const)
     3      +COR(i-1,jm-1)*D(i-1,jm-1)*2*(D(i-1,jm)/1000+const)
     4       )/grav*dvm(i,jm)
          mxd(nrw+nrs+nre+nrn-n+1)=
     1    max(mxd(nrw+nrs+nre+nrn-n+1),D(i-1,jm)/1000)
        end do
      end do
c
      do nnn=1,3
c
c----------  Calculating the first guess of transport ----------
c
      do j=1,jm
        do i=1,im
          D(i,j)=D(i,j)*dum(i,j)*dvm(i,j)-0.0001*(dum(i,j)*dvm(i,j)-1)
          uab(i,j)=uab(i,j)*dum(i,j)
          vab(i,j)=vab(i,j)*dvm(i,j)
        end do
      end do
      EL(1,jm) = 0.0
      EL(2,jm) = 0.0
c-------- falk 07-20-01 correct calculation of elevation
c     EL(3,jm) = 0.0
c
c----------------  WESTERN  BOUNDARY  -------------------------------
c
      print *,'   Western boundary, iteration nnn=',nnn
      do n=nrw,1,-1
        trnsw_inc(n)=0.0
        secarw(n)=0.0
        do j=indw(n)+1,indw(n+1)
          secar=.25E0*(D(2,j)+D(1,j))*(DY(2,j)+DY(1,j))
          trnsw_inc(n)=trnsw_inc(n)+secar*UAB(2,j)
          secarw(n)=secarw(n)+secar
        end do
c--- Calculate new velocity
        if(secarw(n).gt.0) then
          secvelw(n)=(trnsw(n)-trnsw_inc(n))/secarw(n)
        else
          secvelw(n)=0.0
        end if
         print*,'n,secvelw(n)=',n,secvelw(n)
         print *,'trnsw(n)=',trnsw(n)
         print *,'trnsw_inc(n)=',trnsw_inc(n)
        do j=indw(n+1),indw(n)+1,-1
          uab(2,j)=uab(2,j)+secvelw(n)*dum(2,j)
          uab(1,j)=uab(2,j)
          uab(3,j)=uab(2,j)
        end do
        do j=indw(n+1),indw(n)+1,-1
          EL(2,j-1)=EL(2,j)
     1    +ARV(2,j)/((DX(2,j)+DX(2,j-1))*(D(2,j)+D(2,j-1)))
     2    *( COR(2,j)*D(2,j)*(UAB(3,j)+UAB(2,j))
     3      +COR(2,j-1)*D(2,j-1)*(UAB(3,j-1)+UAB(2,j-1))
     4      +4*TRNV(2,j)/ARV(2,j) )/grav
          EL(1,j-1)=EL(2,j-1)
c-------- falk 07-20-01 correct calculation of elevation
c         EL(3,j-1)=EL(2,j-1)
        end do
      end do
c     print *,'   Western boundary, iteration nnn=',nnn
      print *,' EL(2,jm-1)=',EL(2,jm-1)
      print *,' EL(2,2)=',EL(2,2)
      do n=1,nrw
       jbeg=indw(n)
       jend=indw(n+1)
c      print *,'n,jbeg,jend=',n,jbeg,jend
c      print *,'(uab(2,j),j=jbeg,jend)'
c      print *,'(el(2,j),j=jbeg,jend)'
c      write(6,201) (el(2,j),j=jbeg,jend)
c      print *,'(dum(2,j),j=jbeg,jend)'
c      write(6,200) (dum(2,j),j=jbeg,jend)
 200   format(10f7.0)
      end do
c
c----------------  SOUTHERN  BOUNDARY  ------------------------------
c
      print *,'   Southern boundary, iteration nnn=',nnn
      do n=1,nrs
        trnss_inc(n)=0.0
        secars(n)=0.0
        do i=inds(n)+1,inds(n+1)
          secar=.25E0*(D(I,2)+D(I,1))*(DX(I,2)+DX(I,1))
          trnss_inc(n)=trnss_inc(n)+secar*VAB(I,2)
          secars(n)=secars(n)+secar
        end do
c--- Calculate new velocity 
        if(secars(n).gt.0) then
          secvels(n)=(trnss(n)-trnss_inc(n))/secars(n)
        else
          secvels(n)=0.0
        end if
        print*,'n, secvels(n)=',n,secvels(n)
        print *,'trnss(n)=',trnss(n)
        print *,'trnss_inc(n)=',trnss_inc(n)
        do i=inds(n)+1,inds(n+1)
          vab(i,2)=vab(i,2)+secvels(n)*dvm(i,2)
          vab(i,1)=vab(i,2)
          vab(i,3)=vab(i,2)
        end do
        do i=inds(n)+1,inds(n+1)
c-------- falk 07-20-01 correct calculation of elevation
         if(i.eq.2) go to 101
          EL(i,2)=EL(i-1,2) 
     1    +ARU(i,2)/((DY(i,2)+DY(i-1,2))*(D(i,2)+D(i-1,2)))
     2    *( COR(i,2)*D(i,2)*(VAB(i,3)+VAB(i,2))
     3      +COR(i-1,2)*D(i-1,2)*(VAB(i-1,3)+VAB(i-1,2))
     4      -4*TRNU(i,2)/ARU(i,2) )/grav
 101     continue
          EL(i,1)=EL(i,2)
c-------- falk 07-20-01 correct calculation of elevation
c         EL(i,3)=EL(i,2)
         end do
      end do
c---------------- use ient as IEND is used in comblk1.h 
c     print *,'   Southern boundary, iteration nnn=',nnn
      print *,' EL(2,2)=',EL(2,2)
      print *,' EL(im-1,2)=',EL(im-1,2)
      do n=1,nrs
       ibeg=inds(n)
       ient=inds(n+1)
c      print *,'n,ibeg,ient=',n,ibeg,ient
c      print *,'(vab(i,2),i=ibeg,ient)'
c      write(6,201) (vab(i,2),i=ibeg,ient)
c      print *,'n,ibeg,ient=',n,ibeg,ient
c      print *,'(el(i,2),i=ibeg,ient)'
c      write(6,201) (el(i,2),i=ibeg,ient)
c      print *,'(dvm(i,2),i=ibeg,ient)'
c      write(6,200) (dvm(i,2),i=ibeg,ient)
      end do
c
c----------------  EASTERN  BOUNDARY  -------------------------------
c
      print *,'   Eastern boundary, iteration nnn=',nnn
      do n=1,nre
        trnse_inc(n)=0.0
        secare(n)=0.0
        do j=inde(n)+1,inde(n+1)
          secar=.25E0*(D(im,J)+D(im-1,J))*(DY(im,J)+DY(im-1,J))
          trnse_inc(n)=trnse_inc(n)+secar*UAB(im,j)
          secare(n)=secare(n)+secar
        end do
c--- Calculate new velocity
        if(secare(n).gt.0) then
          secvele(n)=(trnse(n)-trnse_inc(n))/secare(n)
        else
          secvele(n)=0.0
        end if
        print*,'n,secvele(n)=',n,secvele(n)
        print *,'trnse(n)=',trnse(n)
        print *,'trnse_inc(n)=',trnse_inc(n)
        do j=inde(n)+1,inde(n+1)
          uab(im,j)=uab(im,j)+secvele(n)*dum(im,j)
          uab(im-1,j)=uab(im,j)
        end do
        do j=inde(n)+1,inde(n+1)
c-------- falk 07-20-01 correct calculation of elevation
         if(j.eq.2) go to 102
          EL(im-1,j)=EL(im-1,j-1)-ARV(im-1,j)
     1    /((DX(im-1,j)+DX(im-1,j-1))*(D(im-1,j)+D(im-1,j-1)))
     2    *( COR(im-1,j)*D(im-1,j)*(UAB(im,j)+UAB(im-1,j))
     3      +COR(im-1,j-1)*D(im-1,j-1)*(UAB(im,j-1)+UAB(im-1,j-1))
     4      +4*TRNV(im-1,j)/ARV(im-1,j) )/grav
 102     continue
          EL(im,j)=EL(im-1,j)
        end do
      end do
c     print *,'   Eastern boundary, iteration nnn=',nnn
      print *,' EL(im-1,2)=',EL(im-1,2)
      print *,' EL(im-1,jm-1)=',EL(im-1,jm-1)
      do n=1,nre
       jbeg=inde(n)
       jend=inde(n+1)
c      print *,'n,jbeg,jend=',n,jbeg,jend
c      print *,'(uab(im,j),j=jbeg,jend)'
c      write(6,201) (uab(im,j),j=jbeg,jend)
c      print *,'n,jbeg,jend=',n,jbeg,jend
c      print *,'(el(im,j),j=jbeg,jend)'
c      write(6,201) (el(im,j),j=jbeg,jend)
c      print *,'(dum(im,j),j=jbeg,jend)'
c      write(6,200) (dum(im,j),j=jbeg,jend)
      end do
c
c----------------  NORTHERN  BOUNDARY  -------------------------------
c
      print *,'   Northern boundary iteration nnn=',nnn
      do n=nrn,1,-1
        trnsn_inc(n)=0.0
        secarn(n)=0.0
        do i=indn(n)+1,indn(n+1)
          secar=.25E0*(D(I,jm)+D(I,jm-1))*(DX(I,jm)+DX(I,jm-1))
          trnsn_inc(n)=trnsn_inc(n)+secar*VAB(I,jm)
          secarn(n)=secarn(n)+secar
        end do
c--- Calculate new velocity
        if(secarn(n).gt.0) then
          secveln(n)=(trnsn(n)-trnsn_inc(n))/secarn(n)
        else
          secveln(n)=0.0
        end if
        print*,'n, secveln(n)=',n,secveln(n)
        print *,'trnsn(n)=',trnsn(n)
        print *,'trnsn_inc(n)=',trnsn_inc(n)
        do i=indn(n)+1,indn(n+1)
          vab(i,jm)=vab(i,jm)+secveln(n)*dvm(i,jm)
          vab(i,jm-1)=vab(i,jm)
        end do
        do i=indn(n+1),indn(n)+1,-1
c-------- falk 07-20-01 correct calculation of elevation
         if(i.eq.im) go to 103
          EL(i-1,jm-1)=EL(i,jm-1)-ARU(i,jm-1)
     1    /((DY(i,jm-1)+DY(i-1,jm-1))*(D(i,jm-1)+D(i-1,jm-1)))
     2    *( COR(i,jm-1)*D(i,jm-1)*(VAB(i,jm)+VAB(i,jm-1))
     3      +COR(i-1,jm-1)*D(i-1,jm-1)*(VAB(i-1,jm)+VAB(i-1,jm-1))
     4      -4*TRNU(i,jm-1)/ARU(i,jm-1) )/grav
 103     continue
          EL(i-1,jm)=EL(i-1,jm-1)
c
        end do
      end do
c     print *,'   Northern boundary iteration nnn=',nnn
      print *,' EL(im-1,jm-1)=',EL(im-1,jm-1)
      print *,' EL(2,jm-1)=',EL(2,jm-1)
      do n=1,nrn
       ibeg=indn(n)
       ient=indn(n+1)
c      print *,'n,ibeg,ient=',n,ibeg,ient
c      print *,'(vab(i,jm),i=ibeg,ient)'
c      write(6,201) (vab(i,jm),i=ibeg,ient)
c      print *,'n,ibeg,ient=',n,ibeg,ient
c      print *,'(el(i,jm),i=ibeg,ient)'
c      write(6,201) (el(i,jm),i=ibeg,ient)
c      print *,'(dvm(i,jm),i=ibeg,ient)'
c      write(6,200) (dvm(i,jm),i=ibeg,ient)
      end do
c
c----------- 03-10-02 check north bnd
c
      TRNUmax=0.
      ARUmax=0.
       DYmax=0.
        Dmax=0.
      do i=1,im
       if(abs(TRNU(i,jm-1)).gt.TRNUmax) TRNUmax=abs(TRNU(i,jm-1))
       if(abs(ARU(i,jm-1)).gt.ARUmax) ARUmax=abs(ARU(i,jm-1))
       if(abs(DY(i,jm-1)).gt.DYmax) DYmax=abs(DY(i,jm-1))
       if(abs(D(i,jm-1)).gt.Dmax) Dmax=abs(D(i,jm-1))
      end do      
c     print *,' check north bnd TRNUmax=',TRNUmax
c     print *,'(TRNU(i,jm-1)/TRNUmax,i=1,im)'
c     write(6,201) (TRNU(i,jm-1)/TRNUmax,i=1,im)
c
c     print *,' check north bnd ARUmax=',ARUmax
c     print *,'(ARU(i,jm-1)/ARUmax,i=1,im)'
c     write(6,201) (ARU(i,jm-1)/ARUmax,i=1,im)
c
c     print *,' check north bnd DYmax=',DYmax
c     print *,'(DY(i,jm-1)/DYmax,i=1,im)'
c     write(6,201) (DY(i,jm-1)/DYmax,i=1,im)
c
c     print *,' check north bnd Dmax=',Dmax
c     print *,'(D(i,jm-1)/Dmax,i=1,im)'
c     write(6,201) (D(i,jm-1)/Dmax,i=1,im)
c
c     print *,'(vab(i,jm-1),i=1,im)'
c     write(6,201) (vab(i,jm-1),i=1,im) 
c     print *,'(COR(i,jm-1)*1.E4,i=1,im)'
c     write(6,201) (COR(i,jm-1)*1.E4,i=1,im) 
c
c     print *,'(dvm(i,jm-1),i=1,im)'
c     write(6,200) (dvm(i,jm-1),i=1,im) 
c     print *,'(FSM(i,jm-1),i=1,im)'
c     write(6,200) (FSM(i,jm-1),i=1,im) 
c     print *,'(FSM(i,jm),i=1,im)'
c     write(6,200) (FSM(i,jm),i=1,im) 
c     print *,'(H(i,jm),i=1,im)'
c     write(6,200) (H(i,jm),i=1,im) 
c     print *,'(H(i,jm-1),i=1,im)'
c     write(6,200) (H(i,jm-1),i=1,im) 
c
c     stop
c
      print *,'  iteration    nnn=',nnn
      print*,'total SSE mismatch is ',EL(1,jm)
      cmp=0.0
      do n=1,nrw+nrs+nre+nrn
        if(mxd(n).gt.0.011) then
          cmp=cmp+abs(wgcf(n))/mxd(n)
        end if
      end do
      const=-EL(1,jm)/cmp
      do n=1,nrw+nrs+nre+nrn
        if(mxd(n).gt.0.011) then
          alpha(n)=const*abs(wgcf(n))/wgcf(n)/mxd(n)
c      print*,'dEL(',n,')=',const*abs(wgcf(n))/mxd(n)
        else
          alpha(n)=0.0
        end if
      end do
c---------- 06-01-01 put additional printing
c     print *,'          nnn=',nnn
c     print *,'    cmp,const'
c     write(6,205) cmp,const
 205  format(5f10.4)
c     print *,'    wgcf'
c     write(6,205) wgcf
c     print *,'    mxd'
c     write(6,205) mxd
c     print *,'    alpha'
c     write(6,205) alpha
c     print *,'    constw'
c     write(6,205) constw
c     print *,'    consts'
c     write(6,205) consts
c     print *,'    conste'
c     write(6,205) conste
c     print *,'    constn'
c     write(6,205) constn
c
      do n=nrw,1,-1
        do j=indw(n+1),indw(n)+1,-1
          uab(2,j)=uab(2,j)
     1    +alpha(nrw-n+1)*(D(2,j)/1000+constw(n))*dum(2,j)
          uab(1,j)=uab(2,j)
          uab(3,j)=uab(2,j)
        end do
      end do
      do n=1,nrs
        do i=inds(n)+1,inds(n+1)
          vab(i,2)=vab(i,2)
     1    +alpha(nrw+n)*(D(i,2)/1000+consts(n))*dvm(i,2)
          vab(i,1)=vab(i,2)
          vab(i,3)=vab(i,2)
        end do
      end do
      do n=1,nre
        do j=inde(n)+1,inde(n+1)
          uab(im,j)=uab(im,j)
     1    +alpha(nrw+nrs+n)*(D(im,j)/1000+conste(n))*dum(im,j)
          uab(im-1,j)=uab(im,j)
        end do
      end do
      do n=nrn,1,-1
        do i=indn(n)+1,indn(n+1)
          vab(i,jm)=vab(i,jm)
     1    +alpha(nrw+nrs+nre+nrn-n+1)*(D(i,jm)/1000+constn(n))*dvm(i,jm)
          vab(i,jm-1)=vab(i,jm)
        end do
      end do

      end do
c
      do j=1,jm
        do i=1,im
          D(i,j)=H(i,j)
        end do
      end do
c
      do j=1,jm
        do i=1,im
           elb(i,j)=el(i,j)
        end do
      end do
c
      return
      end
c
c------------
c
	SUBROUTINE MIXSST
c---------------------------------------------------------------------
c  this subroutine 'mixes' NMC temperature down
c---------------------------------------------------------------------
c---- MXD is the max depth where climatological temperature is modified ---
c---- it is computed as climat. mixed layer depth + transition depth ------
c---- transition depth id calculated on the basis of simple criteria ------
c---- the added temperature gradient should not exceed 2 deg per 50m ------
c
      include 'comblk.h'
      real MXD, DTS
c
c------------ falk 03-26-02 check STM
c
      difmax=0.
      difmin=0.
      do i=1,im
      do j=1,jm
       if(FSM(i,j).gt.0.5) then
        difp=STM(i,j)-TB(i,j,1)
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
      print *,' MIXSST: max(STM-TB(1)),i,j'
      print *,' difmax,idifmax,jdifmax=',difmax,idifmax,jdifmax
      print *,' MIXSST: max(TB(1)-STM),i,j'
      print *,' difmin,idifmin,jdifmin=',difmin,idifmin,jdifmin
c
      call oadjust(0)
      do j=1,jm
        do i=1,im
          if(STM(i,j).gt.TB(i,j,1)) then
            n=1
            do while ((TB(i,j,1)-TB(i,j,n+1)).le.0.5.and.n.lt.kbm1)
              n=n+1
            end do
            n1=1
            MXD=-D(I,J)*ZZ(n)+(STM(i,j)-TB(i,j,1))/2.*50.
            DTS=STM(i,j)-TB(i,j,1)
            do while (-D(I,J)*ZZ(n1).le.MXD.and.n1.lt.kbm1)
              n1=n1+1
            end do
            do k=1,kb
              if(k.lt.n) then
                TB(i,j,k)=TB(i,j,k)+DTS
                T(i,j,k)=T(i,j,k)+DTS
              else
                if(k.lt.n1) TB(i,j,k)=TB(i,j,k)+float(n1-k)/(n1-n)*DTS
                if(k.lt.n1) T(i,j,k)=T(i,j,k)+float(n1-k)*DTS/(n1-n)
              end if
            end do
          else
            TB(i,j,1)=STM(i,j)
            T(i,j,1)=STM(i,j)
          end if
          do k=2,kb
            if(T(i,j,k).gt.T(i,j,k-1)) T(i,j,k)=T(i,j,k-1)
            if(TB(i,j,k).gt.TB(i,j,k-1)) TB(i,j,k)=TB(i,j,k-1)
          end do
        end do
      end do
      call oadjust(0)
c
      return
      end
c
c-------------
c
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
       i1=155
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
c        
c-------------- for debug
c     STM(100,50)=f(100,50,1)+10.
c
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
c        if(i.eq.199.and.(j.ge.185.and.j.le.205)) then
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
c          if(i.eq.idifmax.and.j.eq.jdifmax) then
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
c           if(i.eq.idifmax.and.j.eq.jdifmax) then
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
c           if(i.eq.idifmin.and.j.eq.jdifmin) then
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
c         if(i.eq.idifmax.and.j.eq.jdifmax) then
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
       i1=155
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
c
      SUBROUTINE OADJUSTZ(t,im,jm,nl,iflag,Zlev,FSM)

      real t(im,jm,nl),Zlev(nl),dz(100),FSM(im,jm)
c
      real tmn
      integer counter

      do k=1,nl-1
       dz(k)=Zlev(k+1)-Zlev(k)
      end do
       dz(nl)=0.

      small=1.e-3
      counter=1
      ll=0
      DO WHILE(counter.gt.0)
      counter=0
      ll=ll+1
      do j=1,jm
         do i=1,im
            do n=1,2
               do k=n,nl-1,2
               if(t(i,j,k).lt.t(i,j,k+1).and.fsm(i,j).eq.1.) then
                   tmn=(t(i,j,k)*dz(k)+t(i,j,k+1)*dz(k+1))/
     1             (dz(k)+dz(k+1))
               t(i,j,k)=tmn
               t(i,j,k+1)=tmn-small
      if(iflag.eq.1) print*,'Convective adjustment T at (i,j,k)=',i,j,k
               counter=counter+1
               end if
               end do
            end do
         end do
      end do
      if(ll.eq.1) ii=counter
      END DO
      print*,'OADJUSTZ in ',ii,' points after',ll,' iterations'
c
      return
      end
c
c--------------
c
       subroutine smtobot(t,H,Zlev,im,jm,lm,i1,i2,j1,j2,hmax)
       parameter(mm=500)
       dimension t(im,jm,lm),H(im,jm),Zlev(lm)
       dimension t1(mm),h1(mm)
       print *,'smtobot: i1,i2,j1,j2,kf=',i1,i2,j1,j2,kf
c------------  removes waves with n=f(kf)*step to bottom
c------------  when temp. is defined over ocean up to zlev(lm)
c------------  H=1 over land and H>=11. over ocean
c------- kf=1      removes sinusoidal waves 2h
c------- kf=2      removes sinusoidal waves 3h
c------- kf=3 or 4 removes sinusoidal waves 4h
c------- kf=5      removes sinusoidal waves 5h
c------- kf=6      removes sinusoidal waves 6h
c------- kf=7      removes sinusoidal waves 7h
c------- kf=9      removes sinusoidal waves 8h
c------- kf=10 or 11 removes sinusoidal waves 9h
c
c-------- falk 07-23-01 not to smooth 2 upper levels
c      do l=1,lm
       do l=3,lm
        if(l.le.3) then
         kf=1
        else
         if(l.le.5) then
          kf=2
c         kf=11
         else
          if(l.le.7) then
           kf=5
c          kf=11
          else
           if(l.le.9) then
c           kf=11
            kf=7
c           kf=5
           else
            kf=11
c           kf=5
           end if
          end if
         end if
        end if
c       print *,'  l,kf=',l,kf
       zl=Zlev(l)
       do j=j1,j2
        do i=i1,i2
         h1(i)=H(i,j)
         if(h1(i).gt.1.1) h1(i)=hmax
         t1(i)=t(i,j,l)
        end do
c---------------
c       if(l.eq.7.and.j.eq.80) then
c        print *,'smtobot in y dir, l,j,i1,i2=',l,j,i1,i2
c        print *,'h1(i),i=i1,i2'
c        write(6,101) (h1(i),i=i1,i2)
 101     format(10f7.0)
c        print *,'before smoothing t1(i),i=i1,i2'
c        write(6,102) (t1(i),i=i1,i2)
 102     format(10f7.2)
c       end if
c---------------
        call smooth1D(h1,t1,mm,i1,i2,zl,kf)
c---------------
c       if(l.eq.7.and.j.eq.80) then
c        print *,'after  smoothing t1(i),i=i1,i2'
c        write(6,102) (t1(i),i=i1,i2)
c       end if
c---------------
        do i=i1,i2
         t(i,j,l)=t1(i)
        end do
       end do
c
       do i=i1,i2
        do j=j1,j2
         h1(j)=H(i,j)
         if(h1(j).gt.1.1) h1(j)=hmax
         t1(j)=t(i,j,l)
        end do
c---------------
c       if(l.eq.7.and.i.eq.40) then
c        print *,'smtobot in x dir., l,i,j1,j2=',l,i,j1,j2
c        print *,'h1(j),j=j1,j2'
c        write(6,101) (h1(j),j=j1,j2)
c        print *,'before smoothing t1(j),j=j1,j2'
c        write(6,102) (t1(j),j=j1,j2)
c       end if
c---------------
        call smooth1D(h1,t1,mm,j1,j2,zl,kf)
c---------------
c       if(l.eq.7.and.i.eq.40) then
c        print *,'after  smoothing t1(j),j=j1,j2'
c        write(6,102) (t1(j),j=j1,j2)
c       end if
c---------------
        do j=j1,j2
         t(i,j,l)=t1(j)
        end do
       end do
       end do
       return
       end
c
c-------------
c
       subroutine smooth1D(h1,t1,mm,mbeg,mend,zl,kf)
       parameter(mmm=500)
       dimension t1(mm),h1(mm)
       dimension tr(mmm),mbr(mmm),mer(mmm)
c----------      find intervals with data (over bottom) for file 
c----------      t1(m) from m=mbeg to m=mend 
c----------      mbr(k)<=m<= mer(k)  k-interval kk is a number of
c                                              intervals
         mig=1
         k=1
         do m=mbeg,mend
c-------- falk 12/06/00 check if a point over land
          if(mig.eq.1.and.h1(m).lt.1.1) go to 100
          if(mig.eq.1.and.h1(m).gt.zl) then
            mbr(k)=m
            mig=2
            go to 100
          end if
c-------- falk 12/06/00 check if a point over land
          if(mig.eq.2.and.(h1(m).lt.zl.or.h1(m).lt.1.1)) then
c         if(mig.eq.2.and.h1(m).lt.zl) then
            mer(k)=m-1
            k=k+1
            mig=1
            go to 100
          end if
          if(mig.eq.2.and.m.eq.mend) then
            mer(k)=m
            k=k+1
          end if
 100      continue
         end do
         kk=k-1
c
         do k=1,kk
          m1=mbr(k)
          m2=mer(k)
          jj=m2-m1+1
          if(jj.gt.2) then
c-------    smooth t1 in points near land and islands 
           a=(t1(m1)+t1(m1+1)+t1(m1+2))/3.
           b=(t1(m2)+t1(m2-1)+t1(m2-2))/3.
           if(k.ne.1)  t1(m1)=a
           if(k.ne.kk) t1(m2)=b
           do m=1,jj
            i=m+m1-1
            tr(m)=t1(i)
           end do
c-------    smooth tr in points 1 and jj (near land or bnd of array)
c          a=(tr(1)+tr(2)+tr(3))/3.
c          b=(tr(jj)+tr(jj-1)+tr(jj-2))/3.
c          tr(1)=a
c          tr(jj)=b
c
c-----  filtr waves with n=f(kf)*step 
c
c           call filtr1D(tr,mmm,1,jj,kf)
           do m=1,jj
            i=m+m1-1
            t1(i)=tr(m)
           end do
          end if
         end do
       return
       end
c
c--------------
c
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
      call SP(x,y,yi,R1,xi,lb,ni-1,n,ni)
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
c
c
c-------------------
c
C  *** SPLINE INTERP. BOUND. COND.:SECOND DERIVATIVES IS EQUEL TO 0
C  ***GIVEN Y=Y(X),X(K),K=1,KK CALCULATE R(Z) IN POINTS Z(J),
C  ***J=1,JJ  R1=DR/DZ: KK and JJ number of points processed
C  ***                  KKM<500, KKM is a size of X and Y
C  ***                  JJM<500, JJM is a size of R,R1,Z
C
      SUBROUTINE SP(X,Y,R,R1,Z,KK,JJ,KKM,JJM)
      DIMENSION X(KKM),Y(KKM),R(JJM),R1(JJM),Z(JJM),
     *P(500),PM(500),H(500),A(500),B(500),C(500),D(500),
     *V(500),Q(500)
      DO 200 K=2,KK
  200 H(K)=X(K)-X(K-1)
      PM(1)=0.
      PM(KK)=0.
      K1=KK-1
      DO 201 K=2,K1
      C1=H(K)
      C2=H(K+1)
      B(K)=2.
      C(K)=C2/(C1+C2)
      A(K)=1-C(K)
  201 D(K)=6.*((Y(K+1)-Y(K))/C2-(Y(K)-Y(K-1))/C1)/(C1+C2)
C
      Q(1)=0.
      V(1)=0.
      DO 202 K=2,K1
      P(K)=A(K)*Q(K-1)+B(K)
      Q(K)=-C(K)/P(K)
      V(K)=(D(K)-A(K)*V(K-1))/P(K)
  202 CONTINUE
C
      DO 203 K=2,K1
      I=KK+1-K
  203 PM(I)=Q(I)*PM(I+1)+V(I)
      K=2
C
      DO 204 J=1,JJ
      E=Z(J)
      F=ABS(E-X(K-1))
      G=ABS(X(K)-X(K-1))*1.E-6
      IF(F.LT.G) GOTO 207
  205 IF((E-X(K-1))*(E-X(K))) 207,207,206
  206 IF(K.EQ.KK) GOTO 207
      K=K+1
      GOTO 205
C
  207 A1=PM(K-1)
      A2=PM(K)
      B1=Y(K-1)
      B2=Y(K)
      C1=X(K)-E
      C2=E-X(K-1)
      C3=H(K)
C
      R(J)=(A1*(C1**3)+A2*(C2**3)+(B1*6.-A1*C3*C3)
     **C1+(B2*6.-A2*C3*C3)*C2)/(6.*C3)
      R1(J)=(-A1*C1*C1*0.5+A2*C2*C2*0.5+B2-B1)/C3-
     *(A2-A1)*C3/6.
  204 CONTINUE
      RETURN
      END
