
    
    SUBROUTINE GWF2UNSF1SETBCC1(PID)
    !     ******************************************************************
    !     SET Boundary Conditions (Modified version used in AD)
    !      Sets: MAXAL,HCRITS,TTATM,PRECI,RSO,RR,HCA,TATM,PREC,RSOIL
    !            RTOP,LMINSTEP,RROOT,HBOT
    !     ******************************************************************
    USE GLOBAL, ONLY: IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY:TATM,PREC,RSOIL,RR,HCA,RB,HB,HT,RTOP,&
        &RROOT,RBOT,HCRITA,HBOT,HTOP,GWL0L,&
        &INFTOP,INFBOT,KODTOP,LMINSTEP,&
        &TTATM,PRECI,RSO,MAXAL,CTOP ,CBOT,NS,CT,CB
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: PID
    INTEGER :: I,JJ
    DOUBLE PRECISION :: RTOPOLD=0.
    !      Body
    DO I=1,MAXAL(PID)-1
        !     CHECK: The Loop appears unnecessary, since it only enters the
        !            following IF statement once before exiting
        IF(TATM(PID).EQ.TTATM(PID,I)) THEN
            TATM(PID)=TTATM(PID,I+1)
            PREC(PID)=PRECI(PID,I+1)
            ! modified by AdamS, apparent typo
            RSOIL(PID)=RSO(PID,I+1)
            !        Top of the profile
            IF(INFTOP(PID).GE.0)THEN
                RTOPOLD = RTOP(PID)
                HCRITA(PID)=-ABS(HCA(PID,I+1))
                RTOP(PID)=ABS(RSOIL(PID))-ABS(PREC(PID))
                IF(ABS(RTOPOLD-RTOP(PID)).GT.ABS(RTOP(PID))*0.2 .AND. &
                    &RTOP(PID).LT.0.) THEN
                LMINSTEP(PID)=0
                END IF
                IF(KODTOP(PID).EQ.3) THEN
                    HTOP(PID)=HT(PID,I+1)
                    IF(ABS(HTOP(PID)-HT(PID,I+1)).GT.ABS(HTOP(PID))*0.2) THEN
                        LMINSTEP(PID)=0
                    END IF
                END IF
                RROOT(PID)=ABS(RR(PID,I+1))
            END IF
            !        Bottom of the profile
            IF(INFBOT(PID).GE.0) THEN
                IF(ABS(RBOT(PID)-RB).GT.ABS(RBOT(PID))*0.2) LMINSTEP(PID)=0
                RBOT(PID)=RB
                IF(ABS(HBOT(PID)-HB(PID)-GWL0L).GT.ABS(HBOT(PID))*0.2) THEN
                    LMINSTEP(PID)=0
                END IF
                !         HBOT(PID)=HB(PID)+GWL0L
                
                
            END IF
            
            DO 2 JJ=1,NS
				CTOP(PID, JJ)=CT(PID,JJ,I+1)
				CBOT(PID, JJ)=CB(PID,JJ,I+1)
				IF (KODTOP(PID).EQ.-4) THEN
				IF (PREC(PID)-RSOIL(PID).GT.0.) THEN
					CTOP(PID,JJ)=PREC(PID)/(PREC(PID)-RSOIL(PID))*CT(PID,JJ,I+1)	
				ELSE
					CTOP(PID,JJ)=0.
				END IF	 
				end if
					CBOT(PID,JJ)=0.

2			CONTINUE
            
            
            EXIT
        END IF
    
    END DO
   
    RETURN
    
    END SUBROUTINE