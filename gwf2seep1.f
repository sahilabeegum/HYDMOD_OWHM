!     ******************************************************************
!     Author : Gleb Goussarov
!     Code based on gwf2drn7_NWT.f, in line with the general style of
!     of the MF OWHM project. Note however, that some f77 features have
!     been replaced by f90 features.
!     Subroutines in this package use "IMPLICIT NONE", but variables do
!     follow the naming conventions of f77
!
!     THIS PACKAGE IS LINKED TO THE SUBSIDENCE PACKAGE
!     ******************************************************************
!
      MODULE GWFSEEPMODULE

!     ******************************************************************
!     Module Definition
!     ******************************************************************
       INTEGER, SAVE, POINTER :: NSEEP,ISPCB,NSPVL
       INTEGER, SAVE, POINTER :: MXSEEP,IPRSP
       INTEGER, SAVE, POINTER :: NPSEEP,ISPPB,NNPSP
       INTEGER, SAVE, POINTER :: IUSUB
!      NSEEP : amouNt of SEEPage cells
!      ISPCB : SeePage output unit
!      NSPVL : amouNt of SeePage parameter VaLues
!      MXSEEP: MaX number of SEEPage cells
!      IPRSP : PRint SeePage cells ? (1-yes, 0-no)
!      NPSEEP: amount of Named Parameters
!      ISPPB : maximum amount of active seepage cells + 1
!      NNPSP : Number of Non-Parameter SeePage cells (stress period)
!      IUSUB : Use the SUBsidence package ? (1-yes, 0-no)
       REAL,    DIMENSION(:,:), POINTER, CONTIGUOUS :: SEEP
!      SEEP  : List of SEEPage cells
       CHARACTER(LEN=16),SAVE, DIMENSION(:),POINTER,CONTIGUOUS :: SPAUX
!      SPAUX : AUXiliary data for SeePage cells
       TYPE GWFSEEPTYPE
        INTEGER,POINTER :: NSEEP,ISPCB,NSPVL,MXSEEP,IPRSP
        INTEGER,POINTER :: NPSEEP,ISPPB,NNPSP
        INTEGER,POINTER :: IUSUB
        CHARACTER(LEN=16),DIMENSION(:),  POINTER,CONTIGUOUS::SPAUX
        REAL,             DIMENSION(:,:),POINTER,CONTIGUOUS::SEEP
       END TYPE GWFSEEPTYPE
       TYPE(GWFSEEPTYPE),SAVE :: GWFSEEPDAT(10)
      END MODULE
!
      SUBROUTINE GWF2SEEP1AR(IN,IGRID)
!     ******************************************************************
!     Allocate array storage and Read parameter definitions
!     ******************************************************************
       USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
       USE GLOBAL,       ONLY:LSTCHK
       USE GWFSEEPMODULE,ONLY:NSEEP,MXSEEP,NSPVL,ISPCB,IPRSP,NPSEEP,ISPPB,NNPSP,SPAUX,SEEP,IUSUB
       IMPLICIT NONE
!      Parameters
       CHARACTER(*),PARAMETER :: FMT1 = "(1X,1X, &
       'SEEP1 -- SEEPAGE PACKAGE VERSION 1, 19/8/2015' &
       'INPUT READ FROM UNIT ',I3)"
       CHARACTER(*),PARAMETER :: FMT2 = "(1X,&
       &'MAXIMUM OF ',I6,' ACTIVE SEEPAGE CELLS AT ONE TIME')"
       CHARACTER(*),PARAMETER :: FMT3 = "(1X,&
       &'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')"
       CHARACTER(*),PARAMETER :: FMT4 = "(1X,&
       &'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)"
       CHARACTER(*),PARAMETER :: FMT5 = "(1X,&
       &'AUXILIARY SEEPAGE VARIABLE:',A)"
       CHARACTER(*),PARAMETER :: FMT6 = "(1X,&
       &'LIST OF SEEPAGE CELLS WILL NOT BE PRINTED')"
       CHARACTER(*),PARAMETER :: FMT7 = "(1X,//1X,I5,&
       &' SEEPAGE PARAMETERS')"
       CHARACTER(*),PARAMETER :: HEADER = &
       &'SEEP  NO.  LAYER   ROW   COL     SEEP  EL.  STRESS FACTOR'
!
       CHARACTER(LEN=700) :: LINE
       INTEGER :: IN,IGRID
       INTEGER :: ISTART,ISTOP,N,K,I,IP
       INTEGER :: MXPS,MXACTS,LLOC,NAUX,LSTSUM,LSTBEG,NUMNST,NLST,NINLST
       REAL    :: R
       INTEGER :: CTN
!
       ALLOCATE(NSEEP,MXSEEP,NSPVL,ISPCB,IPRSP,IUSUB)
       ALLOCATE(NPSEEP,ISPPB,NNPSP)
!
       IF(LSTCHK(3)) THEN
        WRITE(IOUT,FMT1) IN
       END IF
!
       NSEEP=0
       NNPSP=0
!      Read the maximum number of seepage cells and I/O unit or flag
!      for Cell-by-cell flow terms.
       CALL URDCOM(IN,IOUT,LINE)
       CALL UPARLSTAL(IN,IOUT,LINE,NPSEEP,MXPS)
!      Read MXACTS - MaXimum number of ACTive Seepage cells
!      Read ISPCB  - SeePage data output unit
       IF(IFREFM.EQ.0)THEN
        READ(LINE,"(2I10)") MXACTS,ISPCB
        LLOC = 21
       ELSE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTS,R,IOUT,IN)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ISPCB ,R,IOUT,IN)
       END IF
       IF(LSTCHK(3))THEN
        WRITE(IOUT,FMT2) MXACTS
        IF(ISPCB.LT.0) WRITE(IOUT,FMT3)
        IF(ISPCB.GT.0) WRITE(IOUT,FMT4) ISPCB
       END IF
!
       ALLOCATE(SPAUX(20))
       NAUX=0
       IPRSP=1
       IUSUB=0
!      Read names of the auxiliary values for seepage
!      Check whether things should be printed or not
       CTN = 1
       DO WHILE(CTN.EQ.1)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
        SELECT CASE(LINE(ISTART:ISTOP))
         CASE ('AUX','AUXILIARY')
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
          IF(NAUX.LT.20) THEN
           NAUX = NAUX+1
           SPAUX(NAUX)=LINE(ISTART:ISTOP)
           IF(LSTCHK(3)) WRITE(IOUT,FMT5) SPAUX(NAUX)
          END IF
         CASE ('CBC','CBCALLOCATE')
!         Left for compatibility with older files, but has no effect
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         CASE ('NOPRINT')
          IF(LSTCHK(3)) WRITE(IOUT,FMT6)
          IPRSP = 0
         CASE ('USESUB')
          IUSUB = 1
         CASE DEFAULT
          CTN = 0
        END SELECT
       END DO
!      In addition to auxiliary values, each cell contains
!      5 input data values and 1 location for cell-by-cell flow
       NSPVL=6+NAUX
!
       ISPPB=MXACTS+1
       MXSEEP=MXACTS+MXPS
       ALLOCATE(SEEP(NSPVL,MXSEEP))
!     ------------------------------------------------------------------
!     Read Named Parameters
       IF(LSTCHK(3)) THEN
        WRITE(IOUT,FMT7) NPSEEP
       END IF
       IF(NPSEEP.GT.0) THEN
        LSTSUM=ISPPB
        DO K=1,NPSEEP
         LSTBEG=LSTSUM
         CALL UPARLSTRP(LSTSUM,MXSEEP,IN,IOUT,IP,'SEEP','SEEP',1,&
                       &NUMNST)
         NLST=LSTSUM-LSTBEG
         IF(NUMNST.EQ.0) THEN
!         Read parameter without instances
          CALL ULSTRD(NLST,SEEP,LSTBEG,NSPVL,MXSEEP,1,IN,IOUT,HEADER,&
                       &SPAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRSP)
         ELSE
!         Read instances
          NINLST=NLST/NUMNST
          DO I=1,NUMNST
           CALL UINSRP(I,IN,IOUT,IP,IPRSP)
           CALL ULSTRD(NLST,SEEP,LSTBEG,NSPVL,MXSEEP,1,IN,IOUT,HEADER,&
                       &SPAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRSP)
           LSTBEG=LSTBEG+NINLST
          END DO
         END IF
        END DO
       END IF
!
       CALL SGWF2SEEP1PSV(IGRID)
       RETURN
      END SUBROUTINE
!
      SUBROUTINE GWF2SEEP1RP(IN,IGRID,IUNIT)
!     ******************************************************************
!     Read seepage Parameters ( Head, Conductance, Bottom Elevation)
!     ******************************************************************
       USE GLOBAL,    ONLY: IOUT,NCOL,NROW,NLAY,IFREFM,BOTM,LBOTM,IBOUND
       USE GLOBAL,    ONLY: LSTCHK
       USE GWFSEEPMODULE, ONLY: NSEEP,MXSEEP,NSPVL,IPRSP,NPSEEP,ISPPB,&
                               &NNPSP,SPAUX,SEEP
!
       IMPLICIT NONE
!      Parameters
       CHARACTER(*),PARAMETER :: FMT1 = "(1X,/1X,&
       &'REUSING NON-PARAMETER SEEPAGE CELL FROM LAST STRESS PERIOD')"
       CHARACTER(*),PARAMETER :: FMT2 = "(1X,/1X,'THE NUMBER OF ACTIVE&
       & SEEPAGE CELLS (',I6,') IS GREATER THAN MXACTS(',I6,')')"
       CHARACTER(*),PARAMETER :: FMT3 = "('DRAIN SET TO BELOW CELL&
       & BOTTOM. MODEL STOPPING. CELL WITH ERROR (IC,IR,IL): ',3I5)"
       CHARACTER(*),PARAMETER :: HEADER =&
       &'SEEP  NO.  LAYER   ROW   COL     SEEP  EL.  CONDUCTANCE'
!      Variables
       INTEGER :: IN,IGRID,IUNIT
       CHARACTER(LEN=5):: CSEEP = "SEEP"
       INTEGER :: ITMP,NP,IOUTU,MXACTS,NAUX,NREAD,N,L,IC,IR,IL
       REAL :: BOT,EL
!      Read ITMP (Number of seepage cells or flag to reuse data) and
!      number of parameters
       CALL SGWF2SEEP1PNT(IGRID)
       IF(NPSEEP.GT.0) THEN
        IF(IFREFM.EQ.0) THEN
         READ(IN,'(2I10)') ITMP,NP
        ELSE
         READ(IN,*)ITMP,NP
        END IF
       ELSE
        NP=0
        IF(IFREFM.EQ.0) THEN
         READ(IN,'(I10)') ITMP
        ELSE
         READ(IN,*)ITMP
        END IF
       END IF

       NAUX = NSPVL-6
       IOUTU= IOUT
       IF(IPRSP.EQ.0) IOUTU=-IOUT
!
       IF(ITMP.LT.0) THEN
        IF(LSTCHK(3)) WRITE(IOUT,FMT1)
       ELSE
        NNPSP=ITMP
       END IF
!
       MXACTS=ISPPB-1
       IF(ITMP.GT.0) THEN
        IF(NNPSP.GT.MXACTS) THEN
         IF(LSTCHK(3)) WRITE(IOUT,FMT2)NNPSP,MXACTS
         CALL USTOP(' ')
        END IF
        CALL ULSTRD(NNPSP,SEEP,1,NSPVL,MXSEEP,1,IN,IOUT,HEADER,&
                   &SPAUX,20,IFREFM,NCOL,NROW,NLAY,5,5,IPRSP)
       END IF
       NSEEP = NNPSP

       CALL PRESET('SEEP')
       IF(NP.GT.0)THEN
        NREAD = NSPVL - 1
        DO N=1,NP
         CALL UPARLSTSUB(IN,'SEEP',IOUTU,'SEEP',SEEP,NSPVL,MXSEEP,&
                        &NREAD,HEADER,SPAUX,20,NAUX)
        END DO
       END IF
       IF(LSTCHK(3)) WRITE (IOUT,"(1X,/1X,I6,' SEEPAGE CELLS')")NSEEP
       DO L=1,NSEEP
        IL = SEEP(1,L)
        IR = SEEP(2,L)
        IC = SEEP(3,L)
        IF(IBOUND(IC,IR,IL).GT.0)THEN
         BOT = BOTM(IC,IR,LBOTM(IL))
         EL = SEEP(4,L)
         IF(EL.LT.BOT .AND. IUNIT.GT.0)THEN
          IF(LSTCHK(1)) WRITE(IOUT,FMT3)IC,IR,IL
         END IF
        END IF
       END DO
       RETURN
      END SUBROUTINE GWF2SEEP1RP
!
      SUBROUTINE GWF2SEEP1AD(KSTP,IGRID)
!     *****************************************************************
!     Apply the effects of Subsidence on seepage flow (ADvance)
!     *****************************************************************
       USE GWFSEEPMODULE, ONLY: NSEEP,SEEP,IUSUB
       USE GLOBAL,        ONLY: IBOUND,SUBLNK
       USE GWFSUBMODULE,  ONLY: DVZ
       IMPLICIT NONE
!
       INTEGER :: KSTP,IGRID
       INTEGER :: L,IL,IR,IC
!      If the seepage package does not need to interact with the SUB
!      package, we can just skip this step
       IF(IUSUB.eq.1)THEN
!       However, we still need to make sure the SUB package is linked
        IF(SUBLNK) THEN
         DO L=1,NSEEP
          IL = SEEP(1,L)
          IR = SEEP(2,L)
          IC = SEEP(3,L)
!
          IF (IBOUND(IC,IR,IL).LE.0) CYCLE
          SEEP(4,L)=SEEP(4,L) - DVZ(IC,IR,IL)
         END DO
        END IF
       END IF
       RETURN
      END SUBROUTINE
!
      SUBROUTINE GWF2SEEP1FM(IGRID)
!     *****************************************************************
!     add seepage flow to the source term (ForMulate)
!     *****************************************************************
       USE GLOBAL,       ONLY:IBOUND,HNEW,RHS,HCOF
       USE GWFSEEPMODULE,ONLY:NSEEP,SEEP
!
       IMPLICIT NONE
!
       INTEGER :: IGRID
       INTEGER,PARAMETER :: IROW = 2,ICOL = 3, ILAY = 1
       INTEGER,PARAMETER :: IEL = 4
       DOUBLE PRECISION  :: EL,C
       INTEGER :: IL,IR,IC,L
!
       CALL SGWF2SEEP1PNT(IGRID)
!
!      If NSEEP<=0, there is no seepage
       IF(NSEEP.LE.0) RETURN
!
       DO L=1,NSEEP
!       Get the coordinates of the seepage cell
        IL = SEEP(ILAY,L)
        IR = SEEP(IROW,L)
        IC = SEEP(ICOL,L)
!       External cells should be skipped
        IF(IBOUND(IC,IR,IL).LE.0) CYCLE
!       Obtain seepage data from internal cells
        EL=SEEP(IEL,L)
!       If the head is lower than seepage head skip this cell
        IF(HNEW(IC,IR,IL).LE.EL) CYCLE
!       Add terms to RHS and HCOF(see page 2-12 of the doc for meaning)
        C=SEEP(5,L)
        HCOF(IC,IR,IL)=HCOF(IC,IR,IL)-C
        RHS(IC,IR,IL)=RHS(IC,IR,IL)-C*EL
       END DO
       RETURN
      END SUBROUTINE
!
      SUBROUTINE GWF2SEEP1BD(KSTP,KPER,IGRID)
!     *****************************************************************
!     calculate volumetric BuDget for seepage
!     *****************************************************************
       USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,HNEW,BUFF
       USE GLOBAL,       ONLY:LSTCHK
       USE GWFBASMODULE, ONLY:ICBCFL,IAUXSV,DELT,PERTIM,TOTIM,VBVL,&
                             &VBNM,MSUM

       USE GWFSEEPMODULE,ONLY:NSEEP,ISPCB,SEEP,NSPVL,SPAUX
!
       IMPLICIT NONE
!      Parameters
       REAL,PARAMETER    :: ZERO = 0.
       CHARACTER(LEN=*),PARAMETER :: fmt1 = "(1X,/1X,A,'   PERIOD ',&
       &I4,'   STEP ',I3)"
       CHARACTER(LEN=*),PARAMETER :: fmt2 = "(1X,'SEEP. ',I6,&
       &'   LAYER ',I3,'   ROW ',I5,'   COL ',I5,'   RATE ',1PG15.6)"
!
       INTEGER :: KSTP,KPER,IGRID
       CHARACTER(LEN=16) :: TEXT
       DOUBLE PRECISION  :: HHNEW,CEL,RATOUT,Q,QQ,C,CC,EL,EEL
       INTEGER :: IBD,IBDLBL,IC,IR,IL,L,NAUX
!
       DATA TEXT /'      DTM SEEPAGE'/
!
       CALL SGWF2SEEP1PNT(IGRID)
!      Initialize Cell-by-cell flow term flag (IBD) and
!      accumulator (RATOUT)
       RATOUT=ZERO
       IBD=0
       IF(ISPCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
       IF(ISPCB.GT.0) IBD=ICBCFL
       IBDLBL=0
!      If Cell-by-cell flows are saved as a list, write header
       IF(IBD.EQ.2) THEN
!        NAUX=NSPVL-5-ISPAL
        NAUX = NSPVL-6
        IF(IAUXSV.EQ.0) NAUX=0
        CALL UBDSV4(KSTP,KPER,TEXT,NAUX,SPAUX,ISPCB,NCOL,NROW,NLAY,&
                   &NSEEP,IOUT,DELT,PERTIM,TOTIM,IBOUND)
       END IF
!      Clear the Buffer
       DO IL=1,NLAY
        DO IR=1,NROW
         DO IC=1,NCOL
          BUFF(IC,IR,IL)=ZERO
         END DO
        END DO
       END DO
!      Check if there are actually seepage cells
       IF(NSEEP.GT.0) THEN
        DO L=1,NSEEP
!        Get cell coordinates
         IL=SEEP(1,L)
         IR=SEEP(2,L)
         IC=SEEP(3,L)
         Q=ZERO
!        If cell is no-flow or constant head, ignore it
         IF(IBOUND(IC,IR,IL).GT.0) THEN
!         Get seepage parameters from the seepage list
          EL=SEEP(4,L)
          EEL=EL
          C=SEEP(5,L)
          HHNEW=HNEW(IC,IR,IL)
!         If the head is higher than seepage, calculate Q=C*(EL-HHNEW)
!         Substract Q from RATOUT
          IF(HHNEW.GT.EEL) THEN
           CC=C
           CEL=C*EL
           QQ=CEL - CC*HHNEW
           Q=QQ
           RATOUT=RATOUT-QQ
          END IF
!         Print Individual rates if requested (ISPCB<0)
          IF(IBD.LT.0) THEN
           IF(LSTCHK(3)) THEN
            IF(IBDLBL.EQ.0) WRITE(IOUT,fmt1) TEXT,KPER,KSTP
            WRITE(IOUT,fmt2) L,IL,IR,IC,Q
           END IF
           IBDLBL=1
          END IF
!         Add Q to the buffer.
          BUFF(IC,IR,IL)=BUFF(IC,IR,IL)+Q
!         If Saving Cell-by-cell flows in a list, write flow.
!         If Returning the flow in the SEEP array, copy it there.
         END IF
         IF(IBD.EQ.2) CALL UBDSVB(ISPCB,NCOL,NROW,IC,IR,IL,Q,&
                       &SEEP(1,L),NSPVL,NAUX,6,IBOUND,NLAY)
!         IF(ISPAL.NE.0) SEEP(NSPVL,L)=Q
         SEEP(NSPVL,L)=Q
        END DO
!       To save Cell-by-cell flow as a 3D array, use UBUDSV
        IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,ISPCB,BUFF,NCOL,NROW,&
                               &NLAY,IOUT)
       END IF
!      Prepare the data (rates,volumes and labels) for printing
       VBVL(3,MSUM) = ZERO
       VBVL(4,MSUM) = RATOUT
       VBVL(2,MSUM) = VBVL(2,MSUM)+RATOUT*DELT
       VBNM(MSUM) = TEXT
       MSUM = MSUM+1
      END SUBROUTINE
!
      SUBROUTINE GWF2SEEP1DA(IGRID)
!     *****************************************************************
!     DeAllocate memory for the SEEP module
!     *****************************************************************
       USE GWFSEEPMODULE
       IMPLICIT NONE
       INTEGER :: IGRID
!
       DEALLOCATE(GWFSEEPDAT(IGRID)%NSEEP )
       DEALLOCATE(GWFSEEPDAT(IGRID)%MXSEEP)
       DEALLOCATE(GWFSEEPDAT(IGRID)%NSPVL )
       DEALLOCATE(GWFSEEPDAT(IGRID)%ISPCB )
       DEALLOCATE(GWFSEEPDAT(IGRID)%IPRSP )
       DEALLOCATE(GWFSEEPDAT(IGRID)%NPSEEP)
       DEALLOCATE(GWFSEEPDAT(IGRID)%ISPPB )
       DEALLOCATE(GWFSEEPDAT(IGRID)%NNPSP )
       DEALLOCATE(GWFSEEPDAT(IGRID)%SPAUX )
       DEALLOCATE(GWFSEEPDAT(IGRID)%SEEP  )
       DEALLOCATE(GWFSEEPDAT(IGRID)%IUSUB )
!      There was something here in DRN that seemed useless to me
      END SUBROUTINE
!
      SUBROUTINE SGWF2SEEP1PNT(IGRID)
!     *****************************************************************
!     change seepage data to a different grid
!     *****************************************************************
       USE GWFSEEPMODULE
       IMPLICIT NONE
       INTEGER :: IGRID
!
       NSEEP =>GWFSEEPDAT(IGRID)%NSEEP
       MXSEEP=>GWFSEEPDAT(IGRID)%MXSEEP
       NSPVL =>GWFSEEPDAT(IGRID)%NSPVL
       ISPCB =>GWFSEEPDAT(IGRID)%ISPCB
       IPRSP =>GWFSEEPDAT(IGRID)%IPRSP
       NPSEEP=>GWFSEEPDAT(IGRID)%NPSEEP
       ISPPB =>GWFSEEPDAT(IGRID)%ISPPB
       NNPSP =>GWFSEEPDAT(IGRID)%NNPSP
       SPAUX =>GWFSEEPDAT(IGRID)%SPAUX
       SEEP  =>GWFSEEPDAT(IGRID)%SEEP
       IUSUB =>GWFSEEPDAT(IGRID)%IUSUB
      END SUBROUTINE
!
      SUBROUTINE SGWF2SEEP1PSV(IGRID)
!     *****************************************************************
!     SaVe seepage data for a grid
!     *****************************************************************
       USE GWFSEEPMODULE
       IMPLICIT NONE
       INTEGER :: IGRID
!
       GWFSEEPDAT(IGRID)%NSEEP =>NSEEP
       GWFSEEPDAT(IGRID)%MXSEEP=>MXSEEP
       GWFSEEPDAT(IGRID)%NSPVL =>NSPVL
       GWFSEEPDAT(IGRID)%ISPCB =>ISPCB
       GWFSEEPDAT(IGRID)%IPRSP =>IPRSP
       GWFSEEPDAT(IGRID)%NPSEEP=>NPSEEP
       GWFSEEPDAT(IGRID)%ISPPB =>ISPPB
       GWFSEEPDAT(IGRID)%NNPSP =>NNPSP
       GWFSEEPDAT(IGRID)%SPAUX =>SPAUX
       GWFSEEPDAT(IGRID)%SEEP  =>SEEP
       GWFSEEPDAT(IGRID)%IUSUB =>IUSUB
!
      END SUBROUTINE
