************************************************************************

      DOUBLE PRECISION FUNCTION FK(IMODEL,H,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M,KS,KR,KK,LAMBDA
      INTEGER PPAR
      DOUBLE PRECISION H,PAR(10)


      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)
      KS=PAR(5)
      BPAR=PAR(6)

      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN               ! VG AND MODIFIED VG
C       BPAR=.5D0
        PPAR=2
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
          QK=QS
          KK=KS
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
          QK=PAR(9)
          KK=PAR(10)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        HMIN=-1.D300**(1.D0/N)/DMAX1(ALFA,1.D0)
        HH=DMAX1(DBLE(H),HMIN)
        QEES=DMIN1((QS-QA)/(QM-QA),.999999999999999D0)
        QEEK=DMIN1((QK-QA)/(QM-QA),QEES)
        HS=-1.D0/ALFA*(QEES**(-1.D0/M)-1.D0)**(1.D0/N)
        HK=-1.D0/ALFA*(QEEK**(-1.D0/M)-1.D0)**(1.D0/N)
        IF(DBLE(H).LT.HK) THEN
          QEE=(1.D0+(-ALFA*HH)**N)**(-M)
          QE =(QM-QA)/(QS-QA)*QEE
          QEK=(QM-QA)/(QS-QA)*QEEK
          FFQ =1.D0-(1.D0-QEE **(1.D0/M))**M
          FFQK=1.D0-(1.D0-QEEK**(1.D0/M))**M
          IF(FFQ.LE.0.D0) FFQ=M*QEE**(1.D0/M)
          KR=(QE/QEK)**BPAR*(FFQ/FFQK)**PPAR*KK/KS
          FK=DMAX1(KS*KR,1.D-37)
        END IF
        IF(H.GE.HK.AND.H.LT.HS) THEN
          KR=(1.D0-KK/KS)/(HS-HK)*(H-HS)+1.D0
          FK=KS*KR
        END IF
        IF(H.GE.HS) FK=KS
      ELSE IF(IMODEL.EQ.2) THEN                  ! BROOKS AND CORES
C       BPAR=1.D0
        LAMBDA=2.D0  !  !=2 FOR MUALEM MODEL, =1.5 FOR BURDINE MODEL
        HS=-1.D0/ALFA
        IF(H.LT.HS) THEN
          KR=1.D0/(-ALFA*H)**(N*(BPAR+LAMBDA)+2.D0)
          FK=DMAX1(KS*KR,1.D-37)
        ELSE
          FK=KS
        END IF
      ELSE IF(IMODEL.EQ.4) THEN                  ! LOG-NORMAL MODEL
        HS=0.D0
        IF(H.LT.HS) THEN
          QEE=QNORM(DLOG(-H/ALFA)/N)
          T=QNORM(DLOG(-H/ALFA)/N+N)
          KR=QEE**BPAR*T*T
          FK=DMAX1(KS*KR,1.D-37)
        ELSE
          FK=KS
        END IF
      END IF

      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FKQ(IMODEL,TH,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M,KS,KR,KK
      INTEGER PPAR
      DOUBLE PRECISION TH,PAR(10)

      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)
      KS=PAR(5)
      BPAR=PAR(6)
      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN               ! VG AND MODIFIED VG
        PPAR=2
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
          QK=QS
          KK=KS
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
          QK=PAR(9)
          KK=PAR(10)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        QEES=DMIN1((QS-QA)/(QM-QA),.999999999999999D0)
        QEEK=DMIN1((QK-QA)/(QM-QA),QEES)
        IF(DBLE(TH).LT.QK) THEN
          QEE=(DBLE(TH)-QA)/(QM-QA)
          QE =(QM-QA)/(QS-QA)*QEE
          QEK=(QM-QA)/(QS-QA)*QEEK
          FFQ =1.D0-(1.D0-QEE **(1.D0/M))**M
          FFQK=1.D0-(1.D0-QEEK**(1.D0/M))**M
          IF(FFQ.LE.0.D0) FFQ=M*QEE**(1.D0/M)
          KR=(QE/QEK)**BPAR*(FFQ/FFQK)**PPAR*KK/KS
          FKQ=DMAX1(KS*KR,1.D-37)
        END IF
        IF(DBLE(TH).GE.QS) FKQ=KS
      END IF

      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FC1(IMODEL,H,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M
      DOUBLE PRECISION H,PAR(10)

      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)
      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        HMIN=-1.D300**(1.D0/N)/DMAX1(ALFA,1.D0)
        HH=DMAX1(H,HMIN)
        QEES=DMIN1((QS-QA)/(QM-QA),.999999999999999D0)
        HS=-1.D0/ALFA*(QEES**(-1.D0/M)-1.D0)**(1.D0/N)
        IF(H.LT.HS) THEN
          C1=(1.D0+(-ALFA*HH)**N)**(-M-1.D0)
          C2=(QM-QA)*M*N*(ALFA**N)*(-HH)**(N-1.D0)*C1
          FC1=MAX(C2,1.D-37)
          RETURN
        ELSE
          FC1=0.0
        END IF
      ELSE IF(IMODEL.EQ.2) THEN
        HS=-1.D0/ALFA
        IF(H.LT.HS) THEN
          C2=(QS-QR)*N*ALFA**(-N)*(-H)**(-N-1.D0)
          FC1=MAX(C2,1.D-37)
        ELSE
          FC1=0.0
        END IF
      ELSE IF(IMODEL.EQ.4) THEN
        HS=0.D0
        IF(H.LT.HS) THEN
          T=DEXP(-1.D0*(DLOG(-H/ALFA))**2.D0/(2.D0*N**2.D0))
          C2=(QS-QR)/(2.D0*3.141592654)**0.5D0/N/(-H)*T
          FC1=MAX(C2,1.D-37)
        ELSE
          FC1=0.0
        END IF
      END IF


      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FQ(IMODEL,H,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M
      DOUBLE PRECISION H,PAR(10)

      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)
      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        HMIN=-1.D300**(1.D0/N)/DMAX1(ALFA,1.D0)
        HH=DMAX1(H,HMIN)
        QEES=DMIN1((QS-QA)/(QM-QA),.999999999999999D0)
        HS=-1.D0/ALFA*(QEES**(-1.D0/M)-1.D0)**(1.D0/N)
        IF(DBLE(H).LT.HS) THEN
          QEE=(1.D0+(-ALFA*HH)**N)**(-M)
          FQ=DMAX1(QA+(QM-QA)*QEE,1.D-37)
          RETURN
        ELSE
          FQ=QS
        END IF
      ELSE IF(IMODEL.EQ.2) THEN
        HS=-1.D0/ALFA
        IF(H.LT.HS) THEN
          QEE=(-ALFA*H)**(-N)
          FQ=DMAX1(QR+(QS-QR)*QEE,1.D-37)
        ELSE
          FQ=QS
        END IF
      ELSE IF(IMODEL.EQ.4) THEN
        HS=0.D0
        IF(H.LT.HS) THEN
          QEE=QNORM(DLOG(-H/ALFA)/N)
          FQ=DMAX1(QR+(QS-QR)*QEE,1.D-37)
        ELSE
          FQ=QS
        END IF
      END IF

      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FH(IMODEL,QE,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M
      DOUBLE PRECISION QE,PAR(10)

      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)

      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        HMIN=-1.D300**(1.D0/N)/MAX(ALFA,1.D0)
        QEEM=(1.D0+(-ALFA*HMIN)**N)**(-M)
        QEE=DMIN1(DMAX1(QE*(QS-QA)/(QM-QA),QEEM),.999999999999999D0)

        FH=MAX(-1.D0/ALFA*(QEE**(-1.D0/M)-1.D0)**(1.D0/N),-1.D37)
      ELSE IF(IMODEL.EQ.2) THEN
        FH=MAX(-1.D0/ALFA*QE**(-1.D0/N),-1.D37)
      ELSE IF(IMODEL.EQ.4) THEN
        IF(QE.GT.0.9999) THEN
          FH=0.0
        ELSE IF(QE.LT.0.00001) THEN
          FH=-1.E+8
        ELSE
          Y=QE*2.D0
          IF(Y.LT.1.) P=DSQRT(-DLOG(Y/2.D0))
          IF(Y.GE.1.) P=DSQRT(-DLOG(1-Y/2.D0))
          X=P-(1.881796+0.9425908*P+0.0546028*P**3)/
     1        (1.+2.356868*P+0.3087091*P**2+0.0937563*P**3
     2         +0.021914*P**4)
          IF(Y.GE.1.) X=-X
          FH=-ALFA*EXP(SQRT(2.)*N*X)
        END IF
      END IF
      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FS(IMODEL,H,PAR)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION N,M
      DOUBLE PRECISION H,PAR(10)

      QR=PAR(1)
      QS=PAR(2)
      ALFA=PAR(3)
      N=PAR(4)
      IF(IMODEL.LE.1.OR.IMODEL.EQ.3) THEN
        IF(IMODEL.EQ.0.OR.IMODEL.EQ.3) THEN
          QM=QS
          QA=QR
        ELSE IF(IMODEL.EQ.1) THEN
          QM=PAR(7)
          QA=PAR(8)
        END IF
        IF(IMODEL.EQ.3) QM=PAR(7)
        M=1.D0-1.D0/N
        QEES=DMIN1((QS-QA)/(QM-QA),.999999999999999D0)
        HS=-1.D0/ALFA*(QEES**(-1.D0/M)-1.D0)**(1.D0/N)
        IF(H.LT.HS) THEN
          HMIN=-1.D300**(1./N)/MAX(ALFA,1.D0)
          HH=DMAX1(H,HMIN)
          QEE=(1.D0+(-ALFA*HH)**N)**(-M)
          QE=QEE*(QM-QA)/(QS-QA)
          FS=MAX(QE,1.D-37)
        ELSE
          FS=1.0
        END IF
      ELSE IF(IMODEL.EQ.2) THEN
        HS=-1.D0/ALFA
        IF(H.LT.HS) THEN
          QE=(-ALFA*H)**(-N)
          FS=DMAX1(QE,1.D-37)
        ELSE
          FS=1.0
        END IF
      ELSE IF(IMODEL.EQ.4) THEN
        HS=0.D0
        IF(H.LT.HS) THEN
          QEE=QNORM(DLOG(-H/ALFA)/N)
          FS=DMAX1(QEE,1.D-37)
        ELSE
          FS=1.0
        END IF
      END IF

      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION QNORM(X)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      Z=DABS(X/2.**0.5)
      T=1./(1.+0.5*Z)
      ERFC=T*DEXP(-Z*Z-1.26551223+T*(1.00002368+T*(0.37409196+
     !     T*(0.09678418+T*(-0.18628806+T*(0.27886807+T*(-1.13520398+
     !     T*(1.48851587+T*(-0.82215223+T*0.17087277)))))))))
      IF(X.LT.0.) ERFC=2.-ERFC
      QNORM=ERFC/2.

      RETURN
      END

************************************************************************

      DOUBLE PRECISION FUNCTION FALFA(TPOT,H,P0,P1,P2H,P2L,P3,R2H,R2L)

      DOUBLE PRECISION:: TPOT,H,P0,P1,P2H,P2L,P3,R2H,R2L
      IF(TPOT.LT.R2L) P2=P2L
      IF(TPOT.GT.R2H) P2=P2H
      IF((TPOT.GE.R2L).AND.(TPOT.LE.R2H))
     !  P2=P2H+(R2H-TPOT)/(R2H-R2L)*(P2L-P2H)
      FALFA=0.0D0
      IF((H.GT.P3).AND.(H.LT.P2)) FALFA=(H-P3)/(P2-P3)
      IF((H.GE.P2).AND.(H.LE.P1)) FALFA=1.0
      IF((H.GT.P1).AND.(H.LT.P0)) FALFA=(H-P0)/(P1-P0)

      RETURN
      END
