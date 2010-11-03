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
