 SUBROUTINE GWF2UNSF1PROFUP1(PID,HDRY)
    ! *************************************************
    ! modified by AdamS, according to the method proposed by Sahila Beegum
    ! update pressure profile in the unsaturated zone
    ! so that the initial bottom flux in the new Modflow time step
    ! matches the average bottom flux from the previous time step
    ! in order to eliminate spurious flux oscillations
    ! profile is calculated using the formulation of the average
    ! internodal conductivity consistent with Hydrus approximation
    ! ***************************************************

    USE GWFUNSFMODULE, ONLY: HHOLD,HHNEW,HTEMP,HB,X,XSURF,NUMNP,NUMNPD,&
        THOLD,THNEW,CON,PAR,MATNUM,HTAB,CONTAB,NTAB,VBOT,THORG,THMOD,xmasschange

    IMPLICIT NONE
    INTEGER,INTENT(IN):: PID
    DOUBLE PRECISION,INTENT(IN):: HDRY ! critical value of pressure for dry soil
    DOUBLE PRECISION, ALLOCATABLE:: HH1(:),HP(:),HZ(:),ZZ(:),ZLAY(:),QQ1(:),QQ2(:),QQ3(:),&
        KPREV(:),KNEXT(:),KNODE(:)
    DOUBLE PRECISION:: QFLUX,KK1,KK2,HC1,HC2,DZ12,DZ,ZGWT,ZMAX,ZMIN,FF1,FF2,KS,HK,&
        DH,ZC,QFLUX1,MAXQ,MINQ,Dx
     INTEGER:: J,INODE,INEXT,ISTEP,IDGW,NLAY,MAT1,MAT2,FLAG,IPROF,K,NODEGW
    INTEGER,ALLOCATABLE:: MATLAY(:)
    !INTEGER,PARAMETER:: NSTEPS=100,NMAX=10000

    ALLOCATE(HH1(NUMNPD),HP(NUMNPD))
    ALLOCATE(QQ1(NUMNP(PID)),QQ2(NUMNP(PID)),QQ3(NUMNP(PID)),&
        KPREV(NUMNP(PID)),KNODE(NUMNP(PID)))
    !HHOLD(1,PID)=HB(PID)
    ! this is the thetaold that should be used for modifying the concentration
    
    DO 22 J=1,NUMNP(PID)
       
        THORG(J,PID)=THOLD(J,PID)
        
22  continue

    HH1=HHOLD(:,PID)
    xmasschange(pid)=0
    
     
    DO J=1,NUMNP(PID)-1
        DZ=X(J+1,PID)-X(J,PID)
        KPREV(J)=(CON(J,PID)+CON(J+1,PID))/2.
        QQ1(J)=-(CON(J,PID)+CON(J+1,PID))/2.*&
            &((HH1(J+1)-HH1(J))/DZ+1.)
    END DO
    write(100,*)' Head_not_adjusted '
    write(100,'(1es12.3)') hh1
    !  DZ=X(2,PID)-X(1,PID)
    !  QFLUX=-(CON(1,PID)+CON(2,PID))/2.*&
    !            &((HH1(2)-HH1(1))/DZ+1.)
    QFLUX = VBOT(PID)
    PRINT *, 'QFLUX=',QFLUX

    HP=HH1
    HP(1)=HB(PID)
    KNODE(1)=FKNODE(HP(1),MATNUM(J,PID))
    IF (HP(1)<0.0) THEN
        PRINT *,'Negative water pressure at the bottom of profile ',PID
        STOP
    ENDIF

    ! compute steady state profile
    DO J=1,NUMNP(PID)-1
        DZ=X(J+1,PID)-X(J,PID)
        CALL SOLVEPRESSURE(QFLUX,HP(J),HC2,KNODE(J),KK2,&
            MATNUM(J,PID), MATNUM(J+1,PID),DZ,HDRY)
        HP(J+1) = HC2
        KNODE(J+1) = KK2
        QQ2(J)=-(KNODE(J)+KNODE(J+1))/2.*&
            &((HP(J+1)-HP(J))/DZ+1.)

    ENDDO
    ! write(100,*)'hp,hhold,qq1,qq2'
    !  write(100,'(4es12.3)')HP,hhold,qq1,qq2
    ! update the bottom part of Hydrus profile
    HHOLD(1,PID)=HP(1)
    DO J=1,NUMNP(PID)-1
        HHOLD(J+1,PID)=HP(J+1)

        ! compute flux in the Hydrus transient profile
        !DZ=X(J+1,PID)-X(J,PID)
        !QFLUX1=-(CON(J,PID)+CON(J+1,PID))/2.*&
        !        &((HH1(J+1)-HH1(J))/DZ+1.)
        ! the flux
        MAXQ=MAX(ABS(QQ2(J)),ABS(QQ1(J)))
        MINQ=MIN(ABS(QQ2(J)),ABS(QQ1(J)))
        !IF (QFLUX1*QFLUX<0.0) EXIT
        IF (ABS(QQ1(J)-QQ2(J))>(1.0d-12+1.0D-3*MAXQ))write(100,*)'the limit',j

        IF (ABS(QQ1(J)-QQ2(J))>(1.0d-12+1.0D-3*MAXQ)) EXIT
    ENDDO

    !HHOLD(:,PID)=HP(:)

    HHNEW(:,PID)=HHOLD(:,PID)  ! this is the modified head
    write (100,*)'Head_adjusted'
    write(100, '(1es12.3)') hHnew
    HTEMP(:,PID)=HHOLD(:,PID)
    CALL GWF2UNSF1SETMAT(PID) ! this takes HHNEW and gives the valus of thnew
    THOLD(:,PID)=THNEW(:,PID)!  this is the modified thta
    ! this is the value of thnew that should be used for modifying the pressure head profile.
  
    DO  23 J=1,NUMNP(PID)
        
        THMOD(J,PID)=THOLD(J,PID)
23  continue
    
 
    CONTAINS

    DOUBLE PRECISION FUNCTION FKNODE(HH,M)
    ! computes hydraulic conductivity as a function of water pressure
    USE GWFUNSFMODULE, ONLY: NTAB,IMODEL,NUMNP,HTAB,CONTAB,MATNUM,PAR,HSAT

    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)::HH ! water pressure
    INTEGER,INTENT(IN)::M ! material index
    !	EXTERNAL FK

    DOUBLE PRECISION::ALH1,DLH,DH,FK
    INTEGER::IT

    ALH1 = DLOG10(-HTAB(1))
    DLH  =(DLOG10(-HTAB(NTAB))-ALH1)/(NTAB-1)
    IF(HH.GE.HSAT(M)) THEN
        FKNODE=PAR(5,M)
    ELSE IF(HH.GT.HTAB(NTAB) .AND. HH.LE.HTAB(1)) THEN
        IT=INT((DLOG10(-HH)-ALH1)/DLH)+1
        DH=(HH-HTAB(IT))/(HTAB(IT+1)-HTAB(IT))
        FKNODE=CONTAB(IT,M)+(CONTAB(IT+1,M)-CONTAB(IT,M))*DH
    ELSE
        FKNODE=FK(IMODEL,HH,PAR(:,M))
    END IF
    END FUNCTION


    DOUBLE PRECISION FUNCTION FNODE(QQ,H1,H2,K1,K2,IMAT1,IMAT2,DZ)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN):: QQ,H1,H2,K1,DZ
    DOUBLE PRECISION,INTENT(OUT):: K2
    INTEGER,INTENT(IN):: IMAT1,IMAT2

    DOUBLE PRECISION:: KAV !,FKNODE

    K2=FKNODE(H2,IMAT2)
    KAV=0.5*(K1+K2)
    FNODE=QQ+KAV/DZ*(H2-H1+DZ)

    END FUNCTION

    SUBROUTINE SOLVEPRESSURE(Q,H1,H2,K1,K2,IMAT1,IMAT2,DZ,HDRY)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN):: Q,H1,K1,DZ,HDRY
    DOUBLE PRECISION,INTENT(OUT):: H2,K2
    INTEGER,INTENT(IN):: IMAT1,IMAT2

    DOUBLE PRECISION:: HA,HB,HC,FA,FB,FC !,FNODE,FKNODE
    INTEGER:: ICOUNT

    IF (ABS(Q)<1.0d-12) THEN
        ! very small flux, pressure distribution approximately hydrostatic
        H2=H1-DZ
        K2=FKNODE(H2,IMAT2)
        RETURN
    ENDIF
    IF ((Q>0).AND.(H1<=HDRY)) THEN
        H2=HDRY
        K2=FKNODE(H2,IMAT2)
        RETURN
    ENDIF
    ICOUNT=0
    IF (Q>0) THEN
        ! upward flux, root interval < (H1-DZ)
        HA=H1-2*DZ
        FA=FNODE(Q,H1,HA,K1,K2,IMAT1,IMAT2,DZ)
        HB=H1-0.9*DZ
        FB=FNODE(Q,H1,HB,K1,K2,IMAT1,IMAT2,DZ)
        DO WHILE (FA*FB>0)
            !PRINT *, HA,HB,FA,FB
            HB=HA
            FB=FA
            HA=HA-DZ
            IF (HA<=HDRY) THEN
                H2=HDRY
                K2=FKNODE(H2,IMAT2)
                RETURN
            ENDIF
            FA=FNODE(Q,H1,HA,K1,K2,IMAT1,IMAT2,DZ)
            ICOUNT=ICOUNT+1
            !IF (ICOUNT==10000) THEN
            !	STOP 'Cannot find root, upward flow'
            !ENDIF
        ENDDO
    ELSE
        ! downward flux, root interval > (H1-DZ)
        HA=H1-1.1*DZ
        FA=FNODE(Q,H1,HA,K1,K2,IMAT1,IMAT2,DZ)
        HB=H1
        FB=FNODE(Q,H1,HB,K1,K2,IMAT1,IMAT2,DZ)
        DO WHILE (FA*FB>0)
            HA=HB
            FA=FB
            HB=HB+DZ
            FB=FNODE(Q,H1,HB,K1,K2,IMAT1,IMAT2,DZ)
            ICOUNT=ICOUNT+1
            IF (ICOUNT==1000) THEN
                STOP 'Cannot find root, downward flow'
            ENDIF
        ENDDO
    ENDIF

    DO
        HC=(HA*FB-HB*FA)/(FB-FA)
        FC=FNODE(Q,H1,HC,K1,K2,IMAT1,IMAT2,DZ)
        !PRINT *,HC,FC
        IF ((ABS(FC/Q)<1.0d-4).OR.(ABS(FC)<1.0D-12)) THEN
            EXIT
        ELSE IF (FA*FC<0) THEN
            HB=HC
            FB=FC
        ELSE
            HA=HC
            FA=FC
        ENDIF
    END DO
    H2=HC

    END SUBROUTINE


    END SUBROUTINE
