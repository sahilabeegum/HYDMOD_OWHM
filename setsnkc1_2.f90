    SUBROUTINE GWF2UNSF1SETSNK1(PID)
        ! ******************************************************************
        ! SET SINKS
        ! WARNING: THETA SHOULD BE LINKED BEFORE DOING THIS
        ! (USE THETA=>THOLD OR THETA=>THNEW)
        ! SETS: SINK,HROOT,VROOT
        ! ******************************************************************
        ! CHECK: TPOT -> RROOT (SEE LINES 575,2097 OF THE ORIGINAL CODE)
        USE GWFUNSFMODULE, ONLY:NMAT
        USE GWFUNSFMODULE, ONLY:NUMNP,HROOT,VROOT,RROOT,P0,POPTM,P2H,&
            &P2L,P3,R2H,R2L,DT
        USE GWFUNSFMODULE, ONLY:MATNUM,X,SINK,HHNEW,BETA,THETA,PAR
        IMPLICIT NONE
        ! VARIABLES
        INTEGER,INTENT(IN) :: PID
        INTEGER :: I,M,N
        DOUBLE PRECISION :: AROOT,DXM,ALFA,FALFA
        ! BODY
        N=NMAT
        VROOT(PID)=0.
        HROOT(PID)=0.
        AROOT=0.
        DO I=2,NUMNP(PID)
            IF (BETA(I,PID).GT.0) THEN
                IF (I.EQ.NUMNP(PID)) THEN
                    DXM=(X(I,PID)-X(I-1,PID))/2.
                ELSE
                    DXM=(X(I+1,PID)-X(I-1,PID))/2.
                END IF
                M=MATNUM(I,PID)
                ALFA=FALFA(RROOT(PID),HHNEW(I,PID),P0(PID),POPTM(PID,M),&
                    &P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID))
                SINK(I,PID)=ALFA*BETA(I,PID)*RROOT(PID)
                IF(THETA(I,PID)-0.00025 .LT. PAR(1,M)) SINK(I,PID)=0.
                SINK(I,PID)=DMIN1(SINK(I,PID),0.5*(THETA(I,PID)-PAR(1,M))/&
                    &DT(PID))
                VROOT(PID)=VROOT(PID)+SINK(I,PID)*DXM
                HROOT(PID)=HROOT(PID)+HHNEW(I,PID)*DXM
                AROOT=AROOT+DXM
            ELSE
                SINK(I,PID)=0.
            END IF
            !PRINT *, PID, I, SINK(I,PID), BETA(I,PID)
        END DO
        !READ(*,*)
        IF(AROOT.GT.0.001) HROOT(PID)=HROOT(PID)/AROOT
        RETURN
        END SUBROUTINE
  
    SUBROUTINE GWF2UNSF1SETSNKC1(PID)
        ! ******************************************************************
        ! SET SINKS
        ! WARNING: THETA SHOULD BE LINKED BEFORE DOING THIS
        ! (USE THETA=>THOLD OR THETA=>THNEW)
        ! SETS: SINK,HROOT,VROOT
        ! ******************************************************************
        ! CHECK: TPOT -> RROOT (SEE LINES 575,2097 OF THE ORIGINAL CODE)
        USE GWFUNSFMODULE, ONLY:NMAT
        USE GWFUNSFMODULE, ONLY:NUMNP,HROOT,VROOT,RROOT,P0,POPTM,P2H,&
            &P2L,P3,R2H,R2L,DT
        USE GWFUNSFMODULE, ONLY:MATNUM,X,SINK,HHNEW,BETA,THETA,PAR
        IMPLICIT NONE
        ! VARIABLES
        INTEGER,INTENT(IN) :: PID
        INTEGER :: I,M,N
        DOUBLE PRECISION :: AROOT,DXM,ALFA,FALFA
        ! BODY
        N=NMAT
        VROOT(PID)=0.
        HROOT(PID)=0.
        AROOT=0.
        DO I=2,NUMNP(PID)
            IF (BETA(I,PID).GT.0) THEN
                IF (I.EQ.NUMNP(PID)) THEN
                    DXM=(X(I,PID)-X(I-1,PID))/2.
                ELSE
                    DXM=(X(I+1,PID)-X(I-1,PID))/2.
                END IF
                M=MATNUM(I,PID)
                ALFA=FALFA(RROOT(PID),HHNEW(I,PID),P0(PID),POPTM(PID,M),&
                    &P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID))
                SINK(I,PID)=ALFA*BETA(I,PID)*RROOT(PID)
                IF(THETA(I,PID)-0.00025 .LT. PAR(1,M)) SINK(I,PID)=0.
                SINK(I,PID)=DMIN1(SINK(I,PID),0.5*(THETA(I,PID)-PAR(1,M))/&
                    &DT(PID))
                VROOT(PID)=VROOT(PID)+SINK(I,PID)*DXM
                HROOT(PID)=HROOT(PID)+HHNEW(I,PID)*DXM
                AROOT=AROOT+DXM
            ELSE
                SINK(I,PID)=0.
            END IF
            !PRINT *, PID, I, SINK(I,PID), BETA(I,PID)
        END DO
        !READ(*,*)
        IF(AROOT.GT.0.001) HROOT(PID)=HROOT(PID)/AROOT
        RETURN
    END SUBROUTINE
    
        SUBROUTINE GWF2UNSF1SETSNK2(PID)
        ! ******************************************************************
        ! SET SINKS
        ! WARNING: THETA SHOULD BE LINKED BEFORE DOING THIS
        ! (USE THETA=>THOLD OR THETA=>THNEW)
        ! SETS: SINK,HROOT,VROOT
        ! ******************************************************************
        ! CHECK: TPOT -> RROOT (SEE LINES 575,2097 OF THE ORIGINAL CODE)
        USE GWFUNSFMODULE, ONLY:NMAT
        USE GWFUNSFMODULE, ONLY:NUMNP,HROOT,VROOT,RROOT,P0,POPTM,P2H,&
            &P2L,P3,R2H,R2L,DT
        USE GWFUNSFMODULE, ONLY:MATNUM,X,SINK,HHNEW,BETA,THETA,PAR
        IMPLICIT NONE
        ! VARIABLES
        INTEGER,INTENT(IN) :: PID
        INTEGER :: I,M,N
        DOUBLE PRECISION :: AROOT,DXM,ALFA,FALFA
        ! BODY
        N=NMAT
        VROOT(PID)=0.
        HROOT(PID)=0.
        AROOT=0.
        DO I=2,NUMNP(PID)
            IF (BETA(I,PID).GT.0) THEN
                IF (I.EQ.NUMNP(PID)) THEN
                    DXM=(X(I,PID)-X(I-1,PID))/2.
                ELSE
                    DXM=(X(I+1,PID)-X(I-1,PID))/2.
                END IF
                M=MATNUM(I,PID)
                ALFA=FALFA(RROOT(PID),HHNEW(I,PID),P0(PID),POPTM(PID,M),&
                    &P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID))
                SINK(I,PID)=ALFA*BETA(I,PID)*RROOT(PID)
                IF(THETA(I,PID)-0.00025 .LT. PAR(1,M)) SINK(I,PID)=0.
                SINK(I,PID)=DMIN1(SINK(I,PID),0.5*(THETA(I,PID)-PAR(1,M))/&
                    &DT(PID))
                VROOT(PID)=VROOT(PID)+SINK(I,PID)*DXM
                HROOT(PID)=HROOT(PID)+HHNEW(I,PID)*DXM
                AROOT=AROOT+DXM
            ELSE
                SINK(I,PID)=0.
            END IF
            !PRINT *, PID, I, SINK(I,PID), BETA(I,PID)
        END DO
        !READ(*,*)
        IF(AROOT.GT.0.001) HROOT(PID)=HROOT(PID)/AROOT
        RETURN
        END SUBROUTINE
   
    SUBROUTINE GWF2UNSF1SETSNKC2(PID)
        ! ******************************************************************
        ! SET SINKS
        ! WARNING: THETA SHOULD BE LINKED BEFORE DOING THIS
        ! (USE THETA=>THOLD OR THETA=>THNEW)
        ! SETS: SINK,HROOT,VROOT
        ! ******************************************************************
        ! CHECK: TPOT -> RROOT (SEE LINES 575,2097 OF THE ORIGINAL CODE)
        USE GWFUNSFMODULE, ONLY:NMAT
        USE GWFUNSFMODULE, ONLY:NUMNP,HROOT,VROOT,RROOT,P0,POPTM,P2H,&
            &P2L,P3,R2H,R2L,DT
        USE GWFUNSFMODULE, ONLY:MATNUM,X,SINK,HHNEW,BETA,THETA,PAR
        IMPLICIT NONE
        ! VARIABLES
        INTEGER,INTENT(IN) :: PID
        INTEGER :: I,M,N
        DOUBLE PRECISION :: AROOT,DXM,ALFA,FALFA
        ! BODY
        N=NMAT
        VROOT(PID)=0.
        HROOT(PID)=0.
        AROOT=0.
        DO I=2,NUMNP(PID)
            IF (BETA(I,PID).GT.0) THEN
                IF (I.EQ.NUMNP(PID)) THEN
                    DXM=(X(I,PID)-X(I-1,PID))/2.
                ELSE
                    DXM=(X(I+1,PID)-X(I-1,PID))/2.
                END IF
                M=MATNUM(I,PID)
                ALFA=FALFA(RROOT(PID),HHNEW(I,PID),P0(PID),POPTM(PID,M),&
                    &P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID))
                SINK(I,PID)=ALFA*BETA(I,PID)*RROOT(PID)
                IF(THETA(I,PID)-0.00025 .LT. PAR(1,M)) SINK(I,PID)=0.
                SINK(I,PID)=DMIN1(SINK(I,PID),0.5*(THETA(I,PID)-PAR(1,M))/&
                    &DT(PID))
                VROOT(PID)=VROOT(PID)+SINK(I,PID)*DXM
                HROOT(PID)=HROOT(PID)+HHNEW(I,PID)*DXM
                AROOT=AROOT+DXM
            ELSE
                SINK(I,PID)=0.
            END IF
            !PRINT *, PID, I, SINK(I,PID), BETA(I,PID)
        END DO
        !READ(*,*)
        IF(AROOT.GT.0.001) HROOT(PID)=HROOT(PID)/AROOT
        RETURN
        END SUBROUTINE