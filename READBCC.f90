
    SUBROUTINE GWF2UNSF1READBCc(PID,IN)
    !     ******************************************************************
    !	Modified by AdamS
    !     Read and set initial values of the Boundary Conditions
    !      Sets: MAXAL,HCRITS,TTATM,PRECI,RSO,RR,HCA,TATM,PREC,RSOIL
    !            RTOP,LMINSTEP,RROOT,HBOT
    !     ******************************************************************
    USE GLOBAL, ONLY: IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY:TATM,PREC,RSOIL,RR,HCA,RB,HB,HT,RTOP,&
        &RROOT,RBOT,HCRITA,HBOT,HTOP,GWL0L,&
        &INFTOP,INFBOT,KODTOP,LMINSTEP,&
        &HCRITS,TTATM,PRECI,RSO,MAXAL,&
        &HHOLD
    USE GWFUNSFMODULE, ONLY:CB,CT,CBOT,CTOP,NS
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: PID,IN
    CHARACTER(LEN=100) :: LINE
    INTEGER :: I,JJ
    DOUBLE PRECISION :: RTOPOLD
    !      Body
    GWL0L=0.0
    RB=0.0
    READ (IN,*) MAXAL(PID),HCRITS(PID)
    CALL UHRCOM (IN,IOUT,LINE)
    BACKSPACE IN
    DO I=1,MAXAL(PID)
        !     CHECK: These values are read, but only the first of each series
        !            is used
        READ(IN,*) TTATM(PID,I),PRECI(PID,I),RSO(PID,I),RR(PID,I),&
            &HCA(PID,I),(cT(pid,jj,i),cB(PID,jj,i),jj=1,NS)
        write(Iout,*) TTATM(PID,I),PRECI(PID,I),RSO(PID,I),RR(PID,I),&
            &HCA(PID,I),(cT(pid,jj,i),cB(PID,jj,i),jj=1,NS)
        HT(PID,I)=0.
    END DO
    HB(PID)=HHOLD(1,PID)
    TATM(PID)=TTATM(PID,1)
    PREC(PID)=PRECI(PID,1)
    RSOIL(PID)=RSO(PID,1)

    ! Top of the profile
    IF(INFTOP(PID).GE.0)THEN
        RTOPOLD = RTOP(PID)
        HCRITA(PID)=-ABS(HCA(PID,1))
        RTOP(PID)=ABS(RSO(PID,1))-ABS(PRECI(PID,1))
        IF(ABS(RTOPOLD-RTOP(PID)).GT.ABS(RTOP(PID))*0.2 .AND. &
            &RTOP(PID).LT.0.) THEN
        LMINSTEP(PID)=0
        END IF
        IF(KODTOP(PID).EQ.3) THEN
            HTOP(PID)=HT(PID,1)
            IF(ABS(HTOP(PID)-HT(PID,1)).GT.ABS(HTOP(PID))*0.2) THEN
                LMINSTEP(PID)=0
            END IF
        END IF
        RROOT(PID)=ABS(RR(PID,1))
    END IF

    ! Bottom of the profile
    IF(INFBOT(PID).GE.0) THEN
        IF(ABS(RBOT(PID)-RB).GT.ABS(RBOT(PID))*0.2) LMINSTEP=0
        RBOT(PID)=RB
        IF(ABS(HBOT(PID)-HB(PID)-GWL0L).GT.ABS(HBOT(PID))*0.2) THEN
            LMINSTEP=0
        END IF
        HBOT(PID)=HB(PID)+GWL0L
    END IF
  
        DO   JJ=1,NS
            CTOP(pid,JJ)=CT(PID,JJ,1)
            CBOT(pid,JJ)=CB(PID,JJ,1)
             write(iout,*)'ctop',CTOP(pid,JJ)
            IF (KODTOP(piD).EQ.-4) THEN
                IF (PREC(PID)-RSOIL(PID).GT.0.) THEN
                    CTOP(pid,JJ)=PREC(PID)/(PREC(PID)-RSOIL(PID))*CT(PID,JJ,1)
                ELSE
                    CTOP(pid,JJ)=0.
                END IF
            end if
            write (iout,*)'ctop',ctop(pid,JJ)
        END DO
        CBOT(pid,JJ)=0.
    RETURN
    END SUBROUTINE