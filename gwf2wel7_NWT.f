      MODULE GWFWELMODULE
        USE TABLEFILE_INTERFACE,ONLY:TABFILETYPE3IDX                    !seb
        INTEGER,SAVE,POINTER  ::NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL
        INTEGER,SAVE,POINTER  ::NPWEL,IWELPB,NNPWEL,IRDPSI
        CHARACTER(LEN=16),SAVE,DIMENSION(:),  POINTER,CONTIGUOUS::WELAUX
        REAL,             SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS::WELL
        REAL,             SAVE,               POINTER     ::PSIRAMP
        INTEGER,          SAVE,               POINTER     ::IUNITRAMP
        TYPE(TABFILETYPE3IDX),SAVE,           POINTER     ::WELTABFILE  !seb ADDED TABFILE SUPPORT FROM TABFILE_INTERFACE
      TYPE GWFWELTYPE
        INTEGER,POINTER  ::NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL
        INTEGER,POINTER  ::NPWEL,IWELPB,NNPWEL,IRDPSI
        CHARACTER(LEN=16), DIMENSION(:),   POINTER,CONTIGUOUS::WELAUX
        REAL,              DIMENSION(:,:), POINTER,CONTIGUOUS::WELL
        REAL,                              POINTER           ::PSIRAMP
        INTEGER,                           POINTER           ::IUNITRAMP
        TYPE(TABFILETYPE3IDX),             POINTER         :: WELTABFILE!seb ADDED TABFILE SUPPORT FROM TABFILE_INTERFACE
      END TYPE
      TYPE(GWFWELTYPE), SAVE:: GWFWELDAT(10)
      END MODULE GWFWELMODULE


      SUBROUTINE GWF2WEL7AR(IN,IUNITNWT,IGRID)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR WELL PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GLOBAL,       ONLY:LSTCHK
      USE GWFWELMODULE, ONLY:NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL,NPWEL,
     1                       IWELPB,NNPWEL,WELAUX,WELL,PSIRAMP,IUNITRAMP
     2                       ,WELTABFILE                                !seb
      USE TABLEFILE_INTERFACE,ONLY: TABFILEPARSE,TABFILELINKS
C
      CHARACTER*700 LINE
C     ------------------------------------------------------------------
      ALLOCATE(NWELLS,MXWELL,NWELVL,IWELCB,IPRWEL)
      ALLOCATE(NPWEL,IWELPB,NNPWEL,PSIRAMP,IUNITRAMP)
C
C1------IDENTIFY PACKAGE AND INITIALIZE NWELLS.
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,1)IN
      ENDIF
    1 FORMAT(1X,/1X,'WEL -- WELL PACKAGE, VERSION 1.0.6, 12/05/2012',
     1' INPUT READ FROM UNIT ',I4)
      NWELLS=0
      NNPWEL=0
      PSIRAMP = 0.1
      IUNITRAMP=IOUT
C
C2------READ MAXIMUM NUMBER OF WELLS AND UNIT OR FLAG FOR
C2------CELL-BY-CELL FLOW TERMS.
      LLOC=1
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPWEL,MXPW)
      !
      ALLOCATE(WELTABFILE)                                             !seb ADDED TABFILE SUPPORT FROM TABFILE_INTERFACE
      CALL TABFILEPARSE(IN,IOUT,LINE,WELTABFILE)
      CALL TABFILELINKS(IN,IOUT,LINE,WELTABFILE)
      !
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') MXACTW,IWELCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTW,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWELCB,R,IOUT,IN)
      END IF
      i = 0
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,3) MXACTW
      ENDIF
    3 FORMAT(1X,'MAXIMUM OF ',I15,' ACTIVE WELLS AT ONE TIME')
      IF(LSTCHK(3)) THEN
        IF(IWELCB.LT.0) WRITE(IOUT,7)
      ENDIF
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(LSTCHK(3)) THEN
        IF(IWELCB.GT.0) WRITE(IOUT,8) IWELCB
      ENDIF
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
    9 FORMAT(1X,'NEGATIVE PUMPING RATES WILL BE REDUCED IF HEAD '/
     +       ' FALLS WITHIN THE INTERVAL PHIRAMP TIMES THE CELL '/
     +       ' THICKNESS. THE VALUE SPECIFIED FOR PHIRAMP IS ',E12.5,/
     +       ' WELLS WITH REDUCED PUMPING WILL BE '
     +       'REPORTED TO FILE UNIT NUMBER',I5)
C
C3------READ AUXILIARY VARIABLES AND PRINT FLAG.
      ALLOCATE(WELAUX(20))
      NAUX=0
      IPRWEL=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.20) THEN
            NAUX=NAUX+1
            WELAUX(NAUX)=LINE(ISTART:ISTOP)
            IF(LSTCHK(3)) THEN
              WRITE(IOUT,12) WELAUX(NAUX)
            ENDIF
   12       FORMAT(1X,'AUXILIARY WELL VARIABLE: ',A)
         END IF
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,13)
         ENDIF
   13    FORMAT(1X,'LISTS OF WELL CELLS WILL NOT BE PRINTED')
         IPRWEL = 0
         GO TO 10
      END IF
! Check keyword for specifying PSI (NWT).
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPP,MXVL)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SPECIFY') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,PSIRAMP,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUNITRAMP,R,IOUT,IN)
         IF(PSIRAMP.LT.1.0E-5) PSIRAMP=1.0E-5
         IF ( IUNITRAMP.EQ.0 ) IUNITRAMP = IOUT
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,*)
         ENDIF
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,9) PSIRAMP,IUNITRAMP
         ENDIF
      ELSE
         BACKSPACE IN
         IF ( IUNITNWT.GT.0 )THEN
           IUNITRAMP = IOUT
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,*)' PHIRAMP WILL BE SET TO A DEFAULT VALUE OF 0.05'
      ENDIF
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,*) ' WELLS WITH REDUCED PUMPING WILL BE '
     +                      ,'REPORTED TO THE MAIN LISTING FILE'
      ENDIF
         END IF
      END IF
!
C3A-----THERE ARE FOUR INPUT VALUES PLUS ONE LOCATION FOR
C3A-----CELL-BY-CELL FLOW.
      NWELVL=5+NAUX
C
C4------ALLOCATE SPACE FOR THE WELL DATA.
      IWELPB=MXACTW+1
      MXWELL=MXACTW+MXPW
      IF(MXACTW.LT.1) THEN
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,17)
         ENDIF
   17    FORMAT(1X,
     1'Deactivating the Well Package because MXACTW=0')
         IN=0
      END IF
      ALLOCATE (WELL(NWELVL,MXWELL))
C
C5------READ NAMED PARAMETERS.
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,18) NPWEL
      ENDIF
   18 FORMAT(1X,//1X,I5,' Well parameters')
      IF(NPWEL.GT.0) THEN
        LSTSUM=IWELPB
        DO 120 K=1,NPWEL
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXWELL,IN,IOUT,IP,'WEL','Q',1,
     &                   NUMINST)
          NLST=LSTSUM-LSTBEG
          IF(NUMINST.EQ.0) THEN
C5A-----READ PARAMETER WITHOUT INSTANCES.
            CALL ULSTRD(NLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.  LAYER   ROW   COL   STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
          ELSE
C5B-----READ INSTANCES.
            NINLST=NLST/NUMINST
            DO 110 I=1,NUMINST
            CALL UINSRP(I,IN,IOUT,IP,IPRWEL)
            CALL ULSTRD(NINLST,WELL,LSTBEG,NWELVL,MXWELL,1,IN,
     &        IOUT,'WELL NO.  LAYER   ROW   COL   STRESS FACTOR',
     &        WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
            LSTBEG=LSTBEG+NINLST
  110       CONTINUE
          END IF
  120   CONTINUE
      END IF
C
C6------RETURN
      CALL SGWF2WEL7PSV(IGRID)
      RETURN
      END
      SUBROUTINE GWF2WEL7RP(IN,IGRID)
C     ******************************************************************
C     READ WELL DATA FOR A STRESS PERIOD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GLOBAL,      ONLY:LSTCHK
      USE GWFWELMODULE, ONLY:NWELLS,MXWELL,NWELVL,IPRWEL,NPWEL,
     1                       IWELPB,NNPWEL,WELAUX,WELL
C
      CHARACTER*6 CWELL
C     ------------------------------------------------------------------
      CALL SGWF2WEL7PNT(IGRID)
C
C1----READ NUMBER OF WELLS (OR FLAG SAYING REUSE WELL DATA).
C1----AND NUMBER OF PARAMETERS
      IF(NPWEL.GT.0) THEN
        IF(IFREFM.EQ.0) THEN
           READ(IN,'(2I10)') ITMP,NP
        ELSE
           READ(IN,*) ITMP,NP
        END IF
      ELSE
         NP=0
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10)') ITMP
         ELSE
            READ(IN,*) ITMP
         END IF
      END IF
C
C------Calculate some constants.
      NAUX=NWELVL-5
      
      IOUTU = IOUT
      IF (IPRWEL.EQ.0) IOUTU=-IOUTU
C
C1A-----IF ITMP LESS THAN ZERO REUSE NON-PARAMETER DATA. PRINT MESSAGE.
C1A-----IF ITMP=>0, SET NUMBER OF NON-PARAMETER WELLS EQUAL TO ITMP.
      IF(ITMP.LT.0) THEN
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,6)
         ENDIF
    6    FORMAT(1X,/
     1    1X,'REUSING NON-PARAMETER WELLS FROM LAST STRESS PERIOD')
      ELSE
         NNPWEL=ITMP
      END IF
C
C1B-----IF THERE ARE NEW NON-PARAMETER WELLS, READ THEM.
      MXACTW=IWELPB-1
      IF(ITMP.GT.0) THEN
         IF(NNPWEL.GT.MXACTW) THEN
            IF(LSTCHK(1)) THEN
              WRITE(IOUT,99) NNPWEL,MXACTW
            ENDIF
   99       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE WELLS (',I6,
     1                     ') IS GREATER THAN MXACTW(',I6,')')
            CALL USTOP(' ')
         END IF
         CALL ULSTRD(NNPWEL,WELL,1,NWELVL,MXWELL,1,IN,IOUT,
     1            'WELL NO.  LAYER   ROW   COL   STRESS RATE',
     2             WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
      END IF
      NWELLS=NNPWEL
C
C1C-----IF THERE ARE ACTIVE WELL PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('Q')
      NREAD=NWELVL-1
      IF(NP.GT.0) THEN
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'WEL',IOUTU,'Q',WELL,NWELVL,MXWELL,NREAD,
     1                MXACTW,NWELLS,4,4,
     2            'WELL NO.  LAYER   ROW   COL   STRESS RATE',
     3            WELAUX,20,NAUX)
   30    CONTINUE
      END IF
C
C3------PRINT NUMBER OF WELLS IN CURRENT STRESS PERIOD.
      CWELL=' WELLS'
      IF(NWELLS.EQ.1) CWELL=' WELL '
      IF(LSTCHK(3)) THEN
        WRITE(IOUT,101) NWELLS,CWELL
      ENDIF
  101 FORMAT(1X,/1X,I6,A)
C
C6------RETURN
      RETURN
      END
      !
      SUBROUTINE GWF2WEL7AD(KSTP,IGRID)
C     ******************************************************************
C     IF TABFILES ARE PRESENT APPLY THEM DEVELOPED BY SCOTT E. BOYCE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFWELMODULE, ONLY:WELL,NWELLS,WELTABFILE
      USE TABLEFILE_INTERFACE,ONLY:TABFILEPACKINDEX,TABFILEUPDATE
      IMPLICIT NONE
      INTEGER::KSTP,IGRID
C Local
      INTEGER,DIMENSION(:,:),ALLOCATABLE:: LRC
C     ------------------------------------------------------------------ 
C
      CALL SGWF2GHB7PNT(IGRID)
C     
      IF(WELTABFILE%NTAB.EQ.0) RETURN   ! NO TABFILES FOR WEL
      !
      IF(KSTP.EQ.1) THEN                                                !NEW STRESS PERIOD, SO THERE MAYBE CHANGES TO WELL ORDER OR COUNT. UPDATE TABFILE INDICIES
         ALLOCATE(LRC(3,NWELLS))
         LRC=INT(WELL(1:3,1:NWELLS))                                    !GET LAY ROW COL INDEX
         !
         CALL TABFILEPACKINDEX(WELTABFILE,INDEX3LST=LRC)     
         !
         DEALLOCATE(LRC)
      END IF
      !
      CALL TABFILEUPDATE( WELTABFILE,'WEL', WELL(4,1:NWELLS) )          !UPDATES THE PUMPING RATES IN WELL
      !
C
C4------RETURN
      RETURN
      END SUBROUTINE
      !
      SUBROUTINE GWF2WEL7FM(Iunitnwt, IGRID)
C     ******************************************************************
C     SUBTRACT Q FROM RHS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IBOUND,RHS,HCOF,LBOTM,BOTM,HNEW,IOUT
      USE GWFWELMODULE, ONLY:NWELLS,WELL,PSIRAMP
      USE GWFNWTMODULE, ONLY: A, IA, Heps, Icell
      USE GWFUPWMODULE, ONLY: LAYTYPUPW
!External function interface
      INTERFACE 
        FUNCTION SMOOTH3(H,T,B,dQ,ISS)
        DOUBLE PRECISION SMOOTH3
        DOUBLE PRECISION, INTENT(IN) :: H
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION, INTENT(IN) :: B
        DOUBLE PRECISION, INTENT(OUT) :: dQ
        INTEGER, INTENT(IN) :: ISS
        END FUNCTION SMOOTH3
      END INTERFACE
!
      DOUBLE PRECISION Qp,Hh,Ttop,Bbot,dQp
      INTEGER Iunitnwt,ISS
C     ------------------------------------------------------------------
      CALL SGWF2WEL7PNT(IGRID)
      ZERO=0.0D0
      Qp = 0.0
C
C1------IF NUMBER OF WELLS <= 0 THEN RETURN.
      IF(NWELLS.LE.0) RETURN
C
C2------PROCESS EACH WELL IN THE WELL LIST.
      DO 100 L=1,NWELLS
      IR=WELL(2,L)
      IC=WELL(3,L)
      IL=WELL(1,L)
      Q=WELL(4,L)
C
C2A-----IF THE CELL IS INACTIVE THEN BYPASS PROCESSING.
      IF(IBOUND(IC,IR,IL).LE.0) GO TO 100
C
C2B-----IF THE CELL IS VARIABLE HEAD THEN SUBTRACT Q FROM
C       THE RHS ACCUMULATOR.
      IF ( Q .LT. ZERO .AND. IUNITNWT.NE.0 ) THEN
        IF ( LAYTYPUPW(il).GT.0 ) THEN
          Hh = HNEW(ic,ir,il)
          bbot = Botm(IC, IR, Lbotm(IL))
          ttop = Botm(IC, IR, Lbotm(IL)-1)
          Qp = Q*smooth3(Hh,Ttop,Bbot,dQp,1)
          RHS(IC,IR,IL)=RHS(IC,IR,IL)-Qp
! Derivative for RHS
          ij = Icell(IC,IR,IL)
          A(IA(ij)) = A(IA(ij)) + dQp
        ELSE
          RHS(IC,IR,IL)=RHS(IC,IR,IL)-Q
        END IF
      ELSE
        RHS(IC,IR,IL)=RHS(IC,IR,IL)-Q
      END IF
  100 CONTINUE
C
C3------RETURN
      RETURN
      END
      SUBROUTINE GWF2WEL7BD(KSTP,KPER,Iunitnwt,IGRID)
C     ******************************************************************
C     CALCULATE VOLUMETRIC BUDGET FOR WELLS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF,BOTM,LBOTM,
     1                      HNEW
      USE GLOBAL,      ONLY:LSTCHK
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,
     1                      VBVL,VBNM
      USE GWFWELMODULE,ONLY:NWELLS,IWELCB,WELL,NWELVL,WELAUX,PSIRAMP,
     1                      IUNITRAMP,IPRWEL
      USE GWFUPWMODULE, ONLY: LAYTYPUPW
!External function interface
      INTERFACE 
        FUNCTION SMOOTH3(H,T,B,dQ,ISS)
        DOUBLE PRECISION SMOOTH3
        DOUBLE PRECISION, INTENT(IN) :: H
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION, INTENT(IN) :: B
        DOUBLE PRECISION, INTENT(OUT) :: dQ
        INTEGER, INTENT(IN) :: ISS
        END FUNCTION SMOOTH3
      END INTERFACE
      CHARACTER*16 TEXT
      DOUBLE PRECISION RATIN,RATOUT,QQ,QSAVE
      double precision Qp,Hh,Ttop,Bbot,dQp
      INTEGER Iunitnwt, iw1, ISS
      DATA TEXT /'           WELLS'/
C     ------------------------------------------------------------------
      CALL SGWF2WEL7PNT(IGRID)
C
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      ZERO=0.
      RATIN=ZERO
      RATOUT=ZERO
      IBD=0
      Qp = 1.0D0               !RICH CHANGED FROM 0.0
      IF(IWELCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
      IF(IWELCB.GT.0) IBD=ICBCFL
      IBDLBL=0
      iw1 = 1
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) THEN
         NAUX=NWELVL-5
         IF(IAUXSV.EQ.0) NAUX=0
         CALL UBDSV4(KSTP,KPER,TEXT,NAUX,WELAUX,IWELCB,NCOL,NROW,NLAY,
     1          NWELLS,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      END IF
C
C3------CLEAR THE BUFFER.
      DO 50 IL=1,NLAY
      DO 50 IR=1,NROW
      DO 50 IC=1,NCOL
      BUFF(IC,IR,IL)=ZERO
50    CONTINUE
C
C4------IF THERE ARE NO WELLS, DO NOT ACCUMULATE FLOW.
      IF(NWELLS.EQ.0) GO TO 200
C
C5------LOOP THROUGH EACH WELL CALCULATING FLOW.
      DO 100 L=1,NWELLS
C
C5A-----GET LAYER, ROW & COLUMN OF CELL CONTAINING WELL.
      IR=WELL(2,L)
      IC=WELL(3,L)
      IL=WELL(1,L)
      Q=ZERO
      QSAVE = ZERO                                                      !seb ADDED TO KEEP BOOK KEEPING CLEAN
      bbot = Botm(IC, IR, Lbotm(IL))                                    !seb MOVED TO KEEP BOOK KEEPING CLEAN
      ttop = Botm(IC, IR, Lbotm(IL)-1)                                  !seb MOVED TO KEEP BOOK KEEPING CLEAN
      Hh = HNEW(ic,ir,il)                                               !seb MOVED TO KEEP BOOK KEEPING CLEAN
C
C5B-----IF THE CELL IS NO-FLOW OR CONSTANT HEAD, IGNORE IT.
      IF(IBOUND(IC,IR,IL).LE.0)GO TO 99
C
C5C-----GET FLOW RATE FROM WELL LIST.
      Q=WELL(4,L)
      QSAVE = Q
      IF ( Q.LT.zero  .AND. Iunitnwt.NE.0) THEN
        IF ( LAYTYPUPW(il).GT.0 ) THEN
          Qp = smooth3(Hh,Ttop,Bbot,dQp,1)
          Q = Q*Qp
        END IF
      END IF
      QQ=Q
C
C5D-----PRINT FLOW RATE IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(LSTCHK(3)) THEN
           IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
         ENDIF
   61    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
         IF(LSTCHK(3)) THEN
           WRITE(IOUT,62) L,IL,IR,IC,Q
         ENDIF
   62    FORMAT(1X,'WELL ',I6,'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '   RATE ',1PG15.6)
         IBDLBL=1
      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
      BUFF(IC,IR,IL)=BUFF(IC,IR,IL)+Q
C
C5F-----SEE IF FLOW IS POSITIVE OR NEGATIVE.
      IF(Q.GE.ZERO) THEN
C
C5G-----FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
        RATIN=RATIN+QQ
      ELSE
C
C5H-----FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
        RATOUT=RATOUT-QQ
      END IF
C
C5I-----IF SAVING CELL-BY-CELL FLOWS IN A LIST, WRITE FLOW.  ALSO
C5I-----COPY FLOW TO WELL LIST.
   99 IF(IBD.EQ.2) CALL UBDSVB(IWELCB,NCOL,NROW,IC,IR,IL,Q,
     1                  WELL(:,L),NWELVL,NAUX,5,IBOUND,NLAY)
      WELL(NWELVL,L)=Q
! write wells with reduced pumping
      IF ( Qp.LT.0.9999D0 .AND. Iunitnwt.NE.0 .AND. 
     +     IPRWEL.NE.0 ) THEN
        IF ( iw1.EQ.1 ) THEN
          WRITE(IUNITRAMP,*)
          WRITE(IUNITRAMP,300)KPER,KSTP
          WRITE(IUNITRAMP,400)
        END IF
ccrth        WRITE(IUNITRAMP,500)IL,IR,IC,QSAVE,Q,FLOAT(hh),FLOAT(bbot)
        WRITE(IUNITRAMP,500)IL,IR,IC,QSAVE,Q,hh,bbot
        iw1 = iw1 + 1
      END IF
  300 FORMAT(' WELLS WITH REDUCED PUMPING FOR STRESS PERIOD ',I5,
     1      ' TIME STEP ',I5)
  400 FORMAT('   LAY   ROW   COL         APPL.Q          ACT.Q',
     1       '        GW-HEAD       CELL-BOT')
  500 FORMAT(3I6,4E15.6)

  100 CONTINUE
      IF (iw1.GT.1 )WRITE(IUNITRAMP,*)
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IWELCB,BUFF,NCOL,NROW,
     1                          NLAY,IOUT)
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(4,MSUM)=ROUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUM=MSUM+1
C
C9------RETURN
      RETURN
      END
      DOUBLE PRECISION FUNCTION smooth3(H,T,B,dQ,ISS)
C     ******************************************************************
C     SMOOTHLY REDUCES PUMPING TO ZERO FOR DEWATERED CONDITIONS
C     ******************************************************************
! h is the depth 
! dC is the derivative of well conductance with respect to well head
      USE GWFWELMODULE,ONLY:PSIRAMP
      USE FMPMODULE, ONLY:PSIRAMPF,SATTHK               !added variable for connection between Smoothing and FMP single-aquifer wells by rth
      IMPLICIT NONE
      DOUBLE PRECISION s, aa, ad, bb, x, y
      DOUBLE PRECISION cof1, cof2, cof3, Qp
      DOUBLE PRECISION, INTENT(IN) :: H
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, INTENT(IN) :: B
      DOUBLE PRECISION, INTENT(OUT) :: dQ
      INTEGER, INTENT(IN) ::  ISS
      smooth3 = 0.0D0
      IF(ISS.eq.1)s = PSIRAMP
      IF(ISS.eq.2)s = PSIRAMPF                    !added variable for connection between Smoothing and FMP single-aquifer wells by rth
      s = s*(T-B)   ! puming rate begins to be ramped down.
      IF(ISS.eq.2.and.SATTHK.lt.s)s=SATTHK
      x = (H-B)
      IF ( x.LT.0.0D0 ) THEN
        Qp = 0.0D0
        dQ = 0.0D0
      ELSEIF ( x-s.GT.-1.0e-14 ) THEN
        Qp = 1.0D0
        dQ = 0.0D0
      ELSE
        aa = -6.0d0/(s**3.0d0)
        bb = -6.0d0/(s**2.0d0)
        cof1 = x**2.0D0
        cof2 = -(2.0D0*x)/(s**3.0D0)
        cof3 = 3.0D0/(s**2.0D0)
        Qp = cof1*(cof2+cof3)
        dQ = (aa*x**2.0D0-bb*x)
      END IF
      smooth3 = Qp
      END FUNCTION smooth3
C
      SUBROUTINE GWF2WEL7DA(IGRID)
C  Deallocate WEL MEMORY
      USE GWFWELMODULE
C
        DEALLOCATE(GWFWELDAT(IGRID)%PSIRAMP) 
        DEALLOCATE(GWFWELDAT(IGRID)%IUNITRAMP) 
        DEALLOCATE(GWFWELDAT(IGRID)%NWELLS)
        DEALLOCATE(GWFWELDAT(IGRID)%MXWELL)
        DEALLOCATE(GWFWELDAT(IGRID)%NWELVL)
        DEALLOCATE(GWFWELDAT(IGRID)%IWELCB)
        DEALLOCATE(GWFWELDAT(IGRID)%IPRWEL)
        DEALLOCATE(GWFWELDAT(IGRID)%NPWEL)
        DEALLOCATE(GWFWELDAT(IGRID)%IWELPB)
        DEALLOCATE(GWFWELDAT(IGRID)%NNPWEL)
        DEALLOCATE(GWFWELDAT(IGRID)%WELAUX)
        DEALLOCATE(GWFWELDAT(IGRID)%WELL)
        DEALLOCATE(GWFWELDAT(IGRID)%WELTABFILE)
C
C NULLIFY THE LOCAL POINTERS
      IF(IGRID.EQ.1)THEN
          PSIRAMP   =>NULL()
          IUNITRAMP =>NULL()
          NWELLS    =>NULL()
          MXWELL    =>NULL()
          NWELVL    =>NULL()
          IWELCB    =>NULL()
          IPRWEL    =>NULL()
          NPWEL     =>NULL()
          IWELPB    =>NULL()
          NNPWEL    =>NULL()
          WELAUX    =>NULL()
          WELL      =>NULL()
          WELTABFILE=>NULL()
      END IF
      RETURN
      END
      SUBROUTINE SGWF2WEL7PNT(IGRID)
C  Change WEL data to a different grid.
      USE GWFWELMODULE
C
        PSIRAMP=>GWFWELDAT(IGRID)%PSIRAMP
        IUNITRAMP=>GWFWELDAT(IGRID)%IUNITRAMP 
        NWELLS=>GWFWELDAT(IGRID)%NWELLS
        MXWELL=>GWFWELDAT(IGRID)%MXWELL
        NWELVL=>GWFWELDAT(IGRID)%NWELVL
        IWELCB=>GWFWELDAT(IGRID)%IWELCB
        IPRWEL=>GWFWELDAT(IGRID)%IPRWEL
        NPWEL=>GWFWELDAT(IGRID)%NPWEL
        IWELPB=>GWFWELDAT(IGRID)%IWELPB
        NNPWEL=>GWFWELDAT(IGRID)%NNPWEL
        WELAUX=>GWFWELDAT(IGRID)%WELAUX
        WELL=>GWFWELDAT(IGRID)%WELL
        WELTABFILE=>GWFWELDAT(IGRID)%WELTABFILE                         !seb
C
      RETURN
      END
      SUBROUTINE SGWF2WEL7PSV(IGRID)
C  Save WEL data for a grid.
      USE GWFWELMODULE
C 
        GWFWELDAT(IGRID)%PSIRAMP=>PSIRAMP
        GWFWELDAT(IGRID)%IUNITRAMP=>IUNITRAMP
        GWFWELDAT(IGRID)%NWELLS=>NWELLS
        GWFWELDAT(IGRID)%MXWELL=>MXWELL
        GWFWELDAT(IGRID)%NWELVL=>NWELVL
        GWFWELDAT(IGRID)%IWELCB=>IWELCB
        GWFWELDAT(IGRID)%IPRWEL=>IPRWEL
        GWFWELDAT(IGRID)%NPWEL=>NPWEL
        GWFWELDAT(IGRID)%IWELPB=>IWELPB
        GWFWELDAT(IGRID)%NNPWEL=>NNPWEL
        GWFWELDAT(IGRID)%WELAUX=>WELAUX
        GWFWELDAT(IGRID)%WELL=>WELL
        GWFWELDAT(IGRID)%WELTABFILE=>WELTABFILE                         !seb
C
      RETURN
      END
