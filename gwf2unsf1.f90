!     ******************************************************************
!     Unsaturated Flow (UNSF) Package
!
!     Author : Gleb Goussarov
!     Code based on gwf1unsf.f by SEO et al., in line with the general
!     style of the MF OWHM project. Note however, that some f77 features
!     have been replaced by f90 features.
!     Subroutines in this package use "IMPLICIT NONE"
!     ******************************************************************
MODULE GWFUNSFMODULE
!     ******************************************************************
!     Module Definition
!     ******************************************************************
  INTEGER,PARAMETER:: ICOUP = 1 ! Adam Szymkiewicz
! ICOUP = 0 no flow bottom boundary condition
! ICOUP = 1 head boundary condition (from MODFLOW)
!  INTEGER,PARAMETER:: IPROU = 97 ! unit for writing profiles
  INTEGER,PARAMETER:: IPROUH = 191 ! unit for writing pressure profiles
  INTEGER,PARAMETER:: IPROUW = 192 ! unit for water content profiles	  
  INTEGER,PARAMETER:: ITWTU = 193 ! unit for writing water table
  INTEGER,PARAMETER:: ITFLU = 194 ! unit for writing bottom flux
  INTEGER,PARAMETER:: ITLOU = 96
       LOGICAL,SAVE,POINTER :: LPRINT,LPTIME,LOUTF,LHEAD
       INTEGER,SAVE,POINTER :: NPUNSF,NUNSFOP,IUNSFCB,IUNSFPR,NIZ,NUMNPD
       INTEGER,SAVE,POINTER :: IMODEL,ITMIN,ITMAX,MAXIT,PROPR,PROINF
       INTEGER,SAVE,POINTER :: MAXATM
!      NPUNSF : Number of Profiles for UNSat Flow
!      NUNSFOP: UNSat Flow OPtion
!      IUNSFCB: UNSat Flow output unit
!      IUNSFPR:
!      NIZ    : Cell/Profile association
!                 0: All cells are associated with the same profile
!                 2: Each zone has its own profile
!      CHECK:     1: Something else
!      NUMPD  :
!      IMODEL :
!      ITMIN  :
!      ITMAX  :
!      MAXIT  :
!      PROPR  :
!      PROINF :
       DOUBLE PRECISION,SAVE,POINTER    :: TOLTH,TOLH
       DOUBLE PRECISION,SAVE,POINTER    :: TPRINT
       DOUBLE PRECISION,SAVE,POINTER    :: XCONV,GWL0L,RB,TOTIMOLD,DMUL,DMUL2
!      TOLH   :
!      TOLTH  :
!      TPRINT :
!      RHENTRY:
!      XCONV  :
!      DMUL   ,DMUL2   :
       INTEGER,SAVE,POINTER :: NMAT,NTAB
!      NMAT   :
!      NTAB   :
       DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS :: HTAB
       DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: PAR
!      ------------------------ SIZE = NPUNSF --------------------------
       INTEGER,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::&
          &IZ,PLEVEL,ALEVEL,TLEVEL,NUMNP,&
          &ITCUM,ITER,LMINSTEP,&
          &KODTOP,KODBOT,INFTOP,INFBOT,KTOLD,KBOLD,NOLAY,NOMAT,&
          &SINKF,SHORTO,ATMBC,&
          &LWAT,WLAYER,LINITW,&
          &MAXAL,PERIMP
!      IZ    : Zone ID for each profile
!      PLEVEL,ALEVEL,TLEVEL:
!      NUMNP :
!      ITCUM :
!      ITER  :
!      LMINSTEP:
!      KODTOP,KODBOT: Type of boundary condition for Top/Bottom
!                     > 0 : Dirichlet BC
!                     ==0 : Can Change from Neumann BC to Dirichlet BC
!                     ==-1: Atmospheric BC
!                     < -1: Neumann BC
!      INFTOP,INFBOT: Time dependency of the boundary condition
!                     >=0 : BC is time-dependent
!                     < 0 : BC is time-independent
!      KTOLD ,KBOLD : Old values of KODTOP and KODBOT
!      NOLAY ,NOMAT :
!      SINKF :
!      SHORTO:
!      ATMBC :
!      LWAT  :
!      WLAYER:
!      LINITW:
!      MAXAL :
!      PERIMP:
       DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::&
          &HBOT,VBOT,HTOP,HROOT,VROOT,WCUMT,WCUMA,WVOLI,SUMVBOT,&
          &RTOP,RBOT,RROOT,HCRITS,HCRITA,XSURF,&
          &TINIT,TMAX,TOLD,T,DT,DTMIN,DTMAX,DTOPT,DTINIT,DTOLD,&
          &P0,P2H,P2L,P3,R2H,R2L,POPTM,&
          &TATM,PREC,RSOIL,&
          &HB,ZONEFLUX ! modifed by Adam Szymkiewicz
!      HBOT  ,VBOT  :
!      HROOT ,VROOT :
!      WCUMT ,VCUMA :
!      WVOLI :
!      SUMVBOT:
!      XSURF :
!      TINIT ,TMAX  ,TOLD  ,T     :
!      DT,DTMIN,DTMAX,DTOPT,DTINIT:
!      TATM  :
!      PREC  :
!      RSOIL :
!     CHECK: This variables is not very useful. Consider removing it.
!      HB    : Same as HBOT, but used in other places
!      -----------------------------------------------------------------
!      -------------------- SIZE = NPUNSF*NUMNPD -----------------------
       INTEGER,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &MATNUM,LAYNUM,WATIN
!      MATNUM:
!      LAYNUM:
!      WATIN :
       DOUBLE PRECISION,   SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &SINK,X,HHNEW,HHOLD,HTEMP,BETA,THOLD,THNEW,CON,CAP
!      SINK  :
!      X     :
!      HHNEW ,HHOLD ,HTEMP : Head values for HYDRUS
!      BETA  :
!      THOLD ,THNEW :
!      CON   :
!      CAP   :
       DOUBLE PRECISION,   SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::THETA
!      THETA : temporary pointer to either THOLD or THNEW
!      -----------------------------------------------------------------
!      -------------------- SIZE = NPUNSF*MAXATM -----------------------
       DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &TTATM,PRECI,RSO,HT,RR,HCA
!      TTATM :
!      PRECI :
!      RSO   :
!      HT    :
!     CHECK: HB: For some reason, this variable is never used as a 2D
!                aray in the original code. So it was moved to the
!                " SIZE = NPUNSF " section
!      RR    :
!      HCA   :
!      -----------------------------------------------------------------
!      ---------------------- SIZE = NTAB*NMAT -------------------------
       DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: CONTAB,CAPTAB,THETAB
	   DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::HSAT
	   ! HSAT modified to 1D array by Adam Szymkiewicz
!      CONTAB:
!      CAPTAB:
!      THETAB:
!      -----------------------------------------------------------------
       DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: UNSFLUX,UNSFLUXV
       DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: CUMQ
       DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS :: PRTIME
!      UNSFLUX :
!      UNSFLUXV:
!      PRTIME  :
!      CUMQ    :
       INTEGER,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: KZON
!      KZON  : ZONE ID associated with each cell
	   
TYPE GWFUNSFTYPE
  LOGICAL,POINTER :: LPRINT,LPTIME,LOUTF,LHEAD
  INTEGER,POINTER :: &
          &NPUNSF,NUNSFOP,IUNSFCB,IUNSFPR,NIZ,NUMNPD,&
          &IMODEL,ITMIN,ITMAX,MAXIT,PROPR,PROINF,&
          &MAXATM,&
          &NMAT,NTAB
  DOUBLE PRECISION   ,POINTER :: &
          &TOLTH,TOLH,&
          &TPRINT,&
          &XCONV,GWL0L,RB,TOTIMOLD,DMUL,DMUL2
  DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS :: HTAB
  DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: PAR
  INTEGER,DIMENSION(:),POINTER,CONTIGUOUS ::&
          &IZ,PLEVEL,ALEVEL,TLEVEL,NUMNP,&
          &ITCUM,ITER,LMINSTEP,&
          &KODTOP,KODBOT,INFTOP,INFBOT,KTOLD,KBOLD,NOLAY,NOMAT,&
          &SINKF,SHORTO,ATMBC,&
          &LWAT,WLAYER,LINITW,&
          &MAXAL,PERIMP
  DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS ::&
          &HBOT,VBOT,HTOP,HROOT,VROOT,WCUMT,WCUMA,WVOLI,SUMVBOT,&
          &RTOP,RBOT,RROOT,HCRITS,HCRITA,XSURF,&
          &TINIT,TMAX,TOLD,T,DT,DTMIN,DTMAX,DTOPT,DTINIT,DTOLD,&
          &P0,P2H,P2L,P3,R2H,R2L,POPTM,&
          &TATM,PREC,RSOIL,&
          &HB,ZONEFLUX ! modifed by Adam Szymkiewicz
  INTEGER,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &MATNUM,LAYNUM,WATIN
  DOUBLE PRECISION,   DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &SINK,X,HHNEW,HHOLD,HTEMP,BETA,THOLD,THNEW,CON,CAP
  DOUBLE PRECISION,   DIMENSION(:,:),POINTER,CONTIGUOUS ::&
          &TTATM,PRECI,RSO,HT,RR,HCA
  DOUBLE PRECISION,   DIMENSION(:,:),POINTER,CONTIGUOUS :: CONTAB,CAPTAB,THETAB
  DOUBLE PRECISION,   DIMENSION(:),POINTER,CONTIGUOUS ::HSAT
  DOUBLE PRECISION,   DIMENSION(:,:),POINTER,CONTIGUOUS :: UNSFLUX,UNSFLUXV
  DOUBLE PRECISION,   DIMENSION(:,:),POINTER,CONTIGUOUS :: CUMQ
  DOUBLE PRECISION,   DIMENSION(:),POINTER,CONTIGUOUS :: PRTIME
  INTEGER,DIMENSION(:,:),POINTER,CONTIGUOUS :: KZON
END TYPE

TYPE(GWFUNSFTYPE),SAVE :: GWFUNSFDAT(10)

END MODULE

	
!     ******************************************************************
SUBROUTINE GWF2UNSF1AR(IN,IGRID,IEVT,INRCH)
!     ******************************************************************
!     Allocate and Read subroutine
!     ******************************************************************
       USE GLOBAL, ONLY: IOUT,NCOL,NROW
       USE GLOBAL, ONLY: LSTCHK
       USE GWFUNSFMODULE
       IMPLICIT NONE
!      Declarations
       INTEGER,INTENT(IN) :: IN,IGRID,IEVT,INRCH
       CHARACTER(LEN=*),PARAMETER :: FMT1 = "(1X,1X,'UNSF -- UNSAT FLOW&
       & PACKAGE,VERSION 1, 22/08/2015','INPUT READ FROM UNIT ',I4)"
       CHARACTER(LEN=*),PARAMETER :: FMT2 = "(1X,I5,' HYDRUS PROFILES&
       & SELECTED')"
       CHARACTER(LEN=*),PARAMETER :: FMT3 = "(1X,'OPTION 1 -- UNSAT&
       & FLUX APPLIES TO THE TOP LAYER')"
       CHARACTER(LEN=*),PARAMETER :: FMT4 = "(1X,'OPTION 2 -- UNSAT&
       & FLUX APPLIES TO HIGHEST ACTIVE CELL IN EACH VERTICAL COLUMN')"
       CHARACTER(LEN=*),PARAMETER :: FMT5 = "(1X,'CELL-BY-CELL FLOWS&
       & WILL BE SAVED ON UNIT ',I4)"
       CHARACTER(LEN=*),PARAMETER :: ERRM1 = '*** ERROR: NO HYDRUS&
       & PROFILE HAS BEEN FOUND !'
       CHARACTER(LEN=*),PARAMETER :: ERRM2 = '*** ERROR: ILLEGAL UNSAT&
       & OPTION CODE. SIMULATION ABORTING'
       CHARACTER(LEN=*),PARAMETER :: ERRM3 = '*** ERROR: CHECK NUMEBER&
       & OF SOIL MATERIALS.'
       INTEGER,PARAMETER :: NMATD=50
       CHARACTER(LEN=700) :: LINE
       INTEGER :: MAXNP
       INTEGER :: LLOC,ISTART,ISTOP
       REAL :: R
!      Body
!      Allocate variables of known size
       CALL GWF2UNSF1PREAL()
!      Identify package
       WRITE(IOUT,FMT1) IN
!      Read Comments and print them, return the first non-comment line
       CALL UHRCOM(IN,IOUT,LINE)
!      Read Profile information
       LLOC = 1
       CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPUNSF ,R,IOUT,IN)
       CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUNSFOP,R,IOUT,IN)
       CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUNSFCB,R,IOUT,IN)
       CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUNSFPR,R,IOUT,IN)
!      Print profile information
       IF(NPUNSF.GT.0) THEN
        IF(LSTCHK(3)) WRITE (IOUT,FMT2) NPUNSF
       ELSE
        NPUNSF=0
        IF(LSTCHK(1)) WRITE(IOUT,'(1X,A)') ERRM1
        CALL USTOP(' ')
       END IF
!      Check if the unsat option is legal
       SELECT CASE(NUNSFOP)
        CASE(1)
         IF(LSTCHK(3)) WRITE (IOUT,FMT3)
        CASE(2)
         IF(LSTCHK(3)) WRITE (IOUT,FMT4)
        CASE DEFAULT
         IF(LSTCHK(1)) WRITE (IOUT,'(1X,A)')ERRM2
       END SELECT
!      Print unit or flag for Cell-by-Cell Flow terms(IUNSFCB)
       IF(IUNSFCB.GT.0) THEN
        IF(LSTCHK(3)) WRITE(IOUT,FMT5) IUNSFCB
       END IF
!      Read reference data to allocate space
       CALL UHRCOM(IN,IOUT,LINE)
       BACKSPACE IN
       LLOC = 1
!     CHECK:  Consider using URWORD(...)
       READ(IN,*) NMAT,MAXNP,MAXATM
       NTAB = 100  ! modified by Adam Szymkiewicz, 26.05.2016
       NUMNPD=MAXNP !+1 modified by Adam Szymkiewicz, 28.05.2016
       IF(NMAT.GT.NMATD)THEN
        IF(LSTCHK(1)) WRITE(IOUT,'(1X,A)') ERRM3
        CALL USTOP(' ')
       END IF
!      Allocate space (Previously IR table)
       ALLOCATE (INFTOP(NPUNSF))
       ALLOCATE (INFBOT(NPUNSF))
       ALLOCATE (KODTOP(NPUNSF))
       ALLOCATE (KODBOT(NPUNSF))
       ALLOCATE (LWAT(NPUNSF))
       ALLOCATE (SINKF(NPUNSF))
       ALLOCATE (WLAYER(NPUNSF))
       ALLOCATE (LMINSTEP(NPUNSF))
       ALLOCATE (ATMBC(NPUNSF))
       ALLOCATE (LINITW(NPUNSF))
       ALLOCATE (NOLAY(NPUNSF))
       ALLOCATE (NOMAT(NPUNSF))
       ALLOCATE (PERIMP(NPUNSF))
       ALLOCATE (KTOLD(NPUNSF))
       ALLOCATE (KBOLD(NPUNSF))
       ALLOCATE (ITCUM(NPUNSF))
       ALLOCATE (TLEVEL(NPUNSF))
       ALLOCATE (ITER(NPUNSF))
       ALLOCATE (NUMNP(NPUNSF))
       ALLOCATE (PLEVEL(NPUNSF))
       ALLOCATE (IZ(NPUNSF))
       ALLOCATE (ALEVEL(NPUNSF))
       ALLOCATE (MAXAL(NPUNSF))
       ALLOCATE (SHORTO(NPUNSF))
!      2D
       ALLOCATE (MATNUM(NUMNPD,NPUNSF))
       ALLOCATE (LAYNUM(NUMNPD,NPUNSF))
       ALLOCATE (KZON(NCOL,NROW))
!      3D
!     UNUSED:  ALLOCATE (IBO(NCOL,NROW,NLAY))
!      Allocate space(Previously RX table)
       ALLOCATE (HTEMP(NUMNPD,NPUNSF))
       ALLOCATE (HHOLD(NUMNPD,NPUNSF))
       ALLOCATE (HHNEW(NUMNPD,NPUNSF))
       ALLOCATE (THETA(NUMNPD,NPUNSF))
       ALLOCATE (THOLD(NUMNPD,NPUNSF))
       ALLOCATE (WATIN(NUMNPD,NPUNSF))
       ALLOCATE (THNEW(NUMNPD,NPUNSF))
       ALLOCATE (SINK (NUMNPD,NPUNSF))
       ALLOCATE (CON  (NUMNPD,NPUNSF))
       ALLOCATE (CAP  (NUMNPD,NPUNSF))
       ALLOCATE (X    (NUMNPD,NPUNSF))
       ALLOCATE (BETA (NUMNPD,NPUNSF))
!
       ALLOCATE (TTATM(NPUNSF,MAXATM))
       ALLOCATE (PRECI(NPUNSF,MAXATM))
       ALLOCATE (RSO  (NPUNSF,MAXATM))
       ALLOCATE (RR   (NPUNSF,MAXATM))
       ALLOCATE (HCA  (NPUNSF,MAXATM))
!       ALLOCATE (HB   (NPUNSF,MAXATM))
       ALLOCATE (HB   (NPUNSF))
       ALLOCATE (HT   (NPUNSF,MAXATM))
!
       ALLOCATE (UNSFLUX(NCOL,NROW))
       ALLOCATE (UNSFLUXV(NCOL,NROW))
!
       ALLOCATE (CAPTAB(NTAB,NMAT))
       ALLOCATE (THETAB(NTAB,NMAT))
       ALLOCATE (CONTAB(NTAB,NMAT))
!
       ALLOCATE (HSAT(NMAT))
       ALLOCATE (HTAB(NTAB))
!
       ALLOCATE (PAR(10,NMAT))
!
!     UNUSED:  ALLOCATE (ADEPTH(NPUNSF))
       ALLOCATE (HCRITS (NPUNSF))
       ALLOCATE (HBOT  (NPUNSF))
       ALLOCATE (HTOP  (NPUNSF))
       ALLOCATE (WCUMT (NPUNSF))
       ALLOCATE (WCUMA (NPUNSF))
       ALLOCATE (WVOLI (NPUNSF))
       ALLOCATE (HROOT (NPUNSF))
       ALLOCATE (VROOT (NPUNSF))
!
       ALLOCATE (T     (NPUNSF))
       ALLOCATE (TOLD  (NPUNSF))
       ALLOCATE (TINIT (NPUNSF))
       ALLOCATE (DT    (NPUNSF))
       ALLOCATE (DTOLD (NPUNSF))
       ALLOCATE (DTOPT (NPUNSF))
       ALLOCATE (DTMIN (NPUNSF))
       ALLOCATE (DTINIT(NPUNSF))
       ALLOCATE (TMAX  (NPUNSF))
!     UNUSED:  ALLOCATE (TFIX  (NPUNSF))
       ALLOCATE (DTMAX (NPUNSF))
!     UNUSED:  ALLOCATE (DTMAXW(NPUNSF))
       ALLOCATE (TATM  (NPUNSF))
!
       ALLOCATE (XSURF (NPUNSF))
       ALLOCATE (RTOP  (NPUNSF))
       ALLOCATE (HCRITA(NPUNSF))
       ALLOCATE (PREC  (NPUNSF))
       ALLOCATE (RSOIL (NPUNSF))
       ALLOCATE (RROOT (NPUNSF))
       ALLOCATE (RBOT  (NPUNSF))
       ALLOCATE (P0    (NPUNSF))
       ALLOCATE (POPTM (NPUNSF))
       ALLOCATE (P2H   (NPUNSF))
       ALLOCATE (P2L   (NPUNSF))
       ALLOCATE (P3    (NPUNSF))
       ALLOCATE (R2H   (NPUNSF))
       ALLOCATE (R2L   (NPUNSF))
       ALLOCATE (VBOT  (NPUNSF))
       ALLOCATE (SUMVBOT(NPUNSF))
	   ALLOCATE (ZONEFLUX(NPUNSF))
       ALLOCATE (CUMQ  (12,NPUNSF))
!     ------------------------------------------------------------------
       CALL GWF2UNSF1FOPEN()
       CALL GWF2UNSF1RPP(IN,IEVT,INRCH)
       CALL GWF2UNSF1PSV(IGRID)
	END SUBROUTINE
	
SUBROUTINE GWF2UNSF1RPP(IN,INEVT,INRCH)
!     ******************************************************************
!     Read and Prepare Parameters
!     ******************************************************************
  USE GLOBAL, ONLY: IOUT,NPER,PERLEN
  USE GWFUNSFMODULE,ONLY:DTINIT,DT,PROPR,PROINF,LINITW,SINKF,THOLD
  USE GWFUNSFMODULE,ONLY:INFTOP,INFBOT,ATMBC,TATM,TMAX,PRTIME,TINIT
  USE GWFUNSFMODULE,ONLY:NPUNSF,THETA,LPRINT,LPTIME,LOUTF,&
                             &LHEAD !,IPROU
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: IN,INEVT,INRCH
  CHARACTER(LEN=100) :: LINE
  INTEGER :: I,K
!     CHECK: In the original code, these variables are arrays, however
!            they are not used, so I have replaced them with DOUBLE PRECISIONs
  DOUBLE PRECISION :: ADEPTH, AHEAD
!     CHECK: ZONARR is a variable present in the original code, at
!            line 439, but not DOUBLE PRECISIONly used for anything - removed
!      Body
!      The following line is required to make certain subroutines work
!      correctly.
       THETA=>THOLD
!     CHECK:
!      Initialize
  CALL GWF2UNSF1INI()
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
  CALL GWF2UNSF1INIT(IN,INEVT,INRCH)
!      Read Time information
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
!      Note: Requires KPER, which is 1 at this point
  CALL GWF2UNSF1TMIN(IN,1)
  DTINIT=DT
!      Read Material information
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
  CALL GWF2UNSF1MATIN(IN)
!      Generate Hydraulic Properties
  CALL GWF2UNSF1GENMAT()
!      Read profile printing information
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
  READ(IN,*) PROPR
  IF (PROPR.GT.0) THEN
    LPRINT=.TRUE. ! modified by Adam Szymkiewicz
    ALLOCATE(PRTIME(PROPR))
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    READ(IN,*) PROINF
	IF (PROINF.GT.0) THEN ! modified by Adam Szymkiewicz
      LHEAD=.FALSE.
	ELSE
	  LHEAD=.TRUE.
	END IF
	CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
!     CHECK: The original code use a temporary variable here
!            It was removed, since it isn't modified during READ
!            See line 505 of the original code for more information
    READ(IN,*) (PRTIME(I),I=1,PROPR)
  ELSE ! modified by Adam Szymkiewicz, print only the final profile
	PROPR=1
    LPRINT=.TRUE.
    ALLOCATE(PRTIME(1))
	PRTIME(1)=0.0
	DO I=1, NPER
	  PRTIME(1)=PRTIME(1)+PERLEN(I)
	END DO
    !PRTIME=>NULL()
  END IF
!      For each profile ...
  DO I=1,NPUNSF
!       Read basic information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    CALL GWF2UNSF1BASINF(I,IN)
!       Read nodal point information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    CALL GWF2UNSF1NODINF(I,IN)
!     CHECK: GWD has been modified to have a different behavior if
!            KPER=0
!            The original code contained a separate function GWDEPTH1,
!            which was used here exclusively. The only differences
!            between GWDEPTH and GWDEPTH1 are as follows:
!             The first line ( T=T+TD ) does not appear in GWDEPTH1
!             No assignments are made to HBOT(HBOT=HB) in GWDEPTH1
!            GWDEPTH  is located at line 1686 in the original code
!            GWDEPTH1 is located at line 2778 in the original code
    CALL GWF2UNSF1GWD(I,0,0,ADEPTH,AHEAD)
!       Read initial condition given in terms of water content
    IF(LINITW(I).GE.0)THEN
      CALL GWF2UNSF1INITW(I)
    END IF
!       Determine nodal values of hydraulic properties
    THETA=>THOLD
    CALL GWF2UNSF1SETMAT(I)
		!DO K = 1, 34
		!  PRINT *, 'TH ', THOLD(K,1), THOLD(K,2)
		!ENDDO
		!READ(*,*)
!       Read root water uptake information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    IF (SINKF(I).GE.0) THEN
      CALL GWF2UNSF1SINKIN(I,IN)
	END IF
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
!       Read atmospheric information
    IF (INFTOP(I).GE.0 .OR. INFBOT(I).GE.0 .OR. ATMBC(I).GE.0) THEN
      CALL GWF2UNSF1SETBC(I,IN)
    ELSE
	  TATM(I)=TMAX(I)
	END IF
!       Calculate Root water extraction rate
    IF (SINKF(I).GE.0) THEN
      THETA=>THOLD
	  CALL GWF2UNSF1SETSNK(I)
	END IF
    CALL GWF2UNSF1SUBREG(I)
  END DO

	   
!      Print information if requested
!       IF(PROPR.GT.0) THEN
!  DO K=1,PROPR
!         IF(ABS(PRTIME(K)-TINIT(1)).LT.0.005)THEN
!          IF(LOUTF)THEN
!           LPTIME=.TRUE.
!           IF(LPRINT)WRITE(IPROU,*) ' Depth (vs) Heads/ Water content'
!          END IF
  CALL GWF2UNSF1NODOUT1(TINIT(1),LHEAD)
!          EXIT
!         END IF
!        END DO
!       END IF
   RETURN
END SUBROUTINE
	
SUBROUTINE GWF2UNSF1AD(IGRID,KSTP,KPER)
!     ******************************************************************
!     calculate unsat flow (ADvance subroutine)
!     ******************************************************************
  USE GWFBASMODULE,  ONLY: TOTIM,DELT
  USE GWFUNSFMODULE, ONLY: NPUNSF,TMAX,THNEW,THETA,SINKF,PLEVEL,SUMVBOT,ZONEFLUX,PERIMP
  USE GWFUNSFMODULE, ONLY: T,TATM,DT,TOLD,TLEVEL,DTOLD
  USE GWFUNSFMODULE, ONLY: INFTOP,INFBOT,ATMBC,ALEVEL
  USE GWFUNSFMODULE, ONLY: KTOLD,KBOLD,KODBOT,KODTOP
  USE GWFUNSFMODULE, ONLY: LPTIME,LPRINT,LHEAD
  USE GWFUNSFMODULE, ONLY: PROPR,PRTIME
  USE GWFUNSFMODULE, ONLY: HHOLD,HHNEW,HTEMP,HB,X,XSURF,NUMNP,THOLD,THNEW

  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: IGRID,KSTP,KPER
  DOUBLE PRECISION :: ADEPTH,AHEAD
  INTEGER :: I,K,SKIP,INOD
!      Body
  CALL GWF2UNSF1PNT(IGRID)
!      The following line is required to make some subroutines work
!      correctly
  THETA=>THNEW
!     CHECK: Nonsensical line in the original code at this point
!            See line 687 of the original code for more information
       !SKIP=0
  !SUMVBOT=0.0  ! Adam Szymkiewicz
  DO I=1,NPUNSF  ! loop over soil profiles
        !IF(SKIP.EQ.0) THEN
    TMAX(I)=TOTIM
    IF(KSTP.EQ.1 .AND. KPER.EQ.1)THEN ! modified by Adam Szymkiewicz
      CONTINUE
!     CHECK: Why copy?
!          IBOUND0(1:NCOL,1:NROW,1:NLAY) = IBOUND(1:NCOL,1:NROW,1:NLAY)
	ELSE
      CALL GWF2UNSF1GWD(I,KSTP,KPER,ADEPTH,AHEAD)
	END IF

!	DO INOD=1,NUMNP(I)
!	  PRINT *,HHOLD(INOD,I), X(INOD,I), XSURF(I)
!	END DO	
!	READ(*,*)
	
	HHOLD(1,I)=HB(I)  
	DO INOD=1,NUMNP(I)-1
	  IF ((X(INOD+1,I)-X(1,I)).GT.HB(I)) EXIT
	  HHOLD(INOD+1,I)=HHOLD(INOD,I)-X(INOD+1,I)+X(INOD,I)
	END DO

	HHNEW(:,I)=HHOLD(:,I)
	HTEMP(:,I)=HHOLD(:,I)
	CALL GWF2UNSF1SETMAT(I)
	THOLD=THNEW
!	END IF
        !END IF
        !SKIP=0
	DO  ! loop over time steps for a single profile
	  !PRINT *, 'T(I)=', T(I), ' DT(I)=', DT(I),' TMAX(I)=',TMAX(I)
	  !READ(*,*)
		  ! solve unsaturated flow for a single time step in the profile
	  !PRINT *, 'T(I)=',T(I), 'DT(I)=',DT(I)
!	  DO INOD=1,NUMNP(I)
!	    PRINT *,I,INOD,HHOLD(INOD,I)
!	  END DO	
!	  READ(*,*)
		
	  CALL GWF2UNSF1WF(I)
!       Calculate time-averaged VBOT
      CALL GWF2UNSF1TOTF(I)
!       Root zone calculations
      IF(SINKF(I).GE.0) CALL GWF2UNSF1SETSNK(I)
!      CALL GWF2UNSF1TLINF(I)
      PLEVEL(I)=PLEVEL(I)+1
!       Read time-dependent boundary condition
      IF(ABS(T(I)-TATM(I)).LE.0.001*DT(I)) THEN
!     CHECK: Redundant IF statement at line 771 in the original code
        IF(INFTOP(I).GE.0 .OR. INFBOT(I).GE.0 .OR. ATMBC(I).GE.0) THEN
          CALL GWF2UNSF1SETBC1(I)
          ALEVEL(I)=ALEVEL(I)+1
        END IF
      END IF
!       Time Governing
!      TOLD(I)=T(I)
!	  DTOLD(I)=DT(I)
!      KTOLD(I)=KODTOP(I)
!      KBOLD(I)=KODBOT(I)
	  IF(ABS(T(I)-TMAX(I)).LE.0.001*DT(I)) THEN
        EXIT
      ELSE
        TOLD(I)=T(I)
		DTOLD(I)=DT(I)
        KTOLD(I)=KODTOP(I)
        KBOLD(I)=KODBOT(I)
        CALL GWF2UNSF1TMCONT(I)
        T(I)=T(I)+DT(I)
        TLEVEL(I)=TLEVEL(I)+1
        IF(TLEVEL(I).GT.1000000) TLEVEL(I)=2
	  END IF
      CALL GWF2UNSF1UPDATE(I)
        !SKIP=1
	END DO ! loop over time steps in the profile	  
	ZONEFLUX(I)=SUMVBOT(I)/DELT*((100.-PERIMP(I))/100.)
	!PRINT *, I, ZONEFLUX(I)
  END DO ! loop over soil profiles

  CALL GWF2UNSF1TLSUM
  
  DO K=1,PROPR
    IF (ABS(PRTIME(K)-T(1)).LT.DT(1)) THEN
      CALL GWF2UNSF1NODOUT1(T(1),LHEAD)
      EXIT
	END IF
  END DO
  
  !  LPTIME=.FALSE.
	   !LPRINT=.TRUE. ! modified by Adam Szymkiewicz
!  IF(LPRINT) THEN
	    !PRINT *, 'T(1)=', T(1)
		!READ(*,*)
!  DO K=1,PROPR
!      IF(ABS(PRTIME(K)-T(1))/PRTIME(K).LT.0.05) THEN
!        LPTIME=.TRUE.
!      END IF
!    END DO
!  END IF
!  IF(LPTIME) THEN
!     CALL GWF2UNSF1NODOUT1(T(1),LHEAD)
  !END IF
  RETURN
END SUBROUTINE
	   
SUBROUTINE GWF2UNSF1FM(IGRID)
!     ******************************************************************
!     ForMulate subroutine
!     ******************************************************************
  USE GLOBAL, ONLY: NROW,NCOL,NLAY, RHS,DELC,DELR,IBOUND
  USE GWFBASMODULE, ONLY: DELT,TOTIM
  USE GWFUNSFMODULE, ONLY: NPUNSF,NUNSFOP,UNSFLUX,SUMVBOT,UNSFLUXV,ZONEFLUX
  USE GWFUNSFMODULE, ONLY: NIZ,PERIMP,IZ,KZON,TOTIMOLD
  IMPLICIT NONE
!      Unknown variables

  INTEGER,INTENT(IN) :: IGRID
  INTEGER ICOUNT
!      Variables
  INTEGER :: I,J,K,IC,IR,IL
!      Body
  CALL GWF2UNSF1PNT(IGRID)
  TOTIMOLD=TOTIM
  ICOUNT=0
!  DO K=1,NPUNSF
!	ZONEFLUX(K)=SUMVBOT(K)/DELT*((100.-PERIMP(K))/100.)
!  END DO
!      Apply SUMVBOT to each cell
  IF(NIZ.EQ.0) THEN
    UNSFLUX(1:NCOL,1:NROW) = ZONEFLUX(1) !SUMVBOT(1)/DELT*((100.-PERIMP(1))/100.)
  ELSE IF (NIZ.EQ.2) THEN
    ICOUNT=1
    DO I=1,NROW
      DO J=1,NCOL
        DO K=1,NPUNSF
          UNSFLUX(J,I) = ZONEFLUX(K) !SUMVBOT(K)/DELT*((100.-PERIMP(K))/100.)
          ICOUNT = ICOUNT+1
!     CHECK: The original code had a line here that I did not understand
        END DO
      END DO
    END DO
  ELSE
    DO K=1,NPUNSF
      DO I=1,NROW
        DO J=1,NCOL
          IF(KZON(J,I).EQ.IZ(K)) THEN
            UNSFLUX(J,I) = ZONEFLUX(K) !SUMVBOT(K)/DELT*((100.-PERIMP(K))/100.)
          ENDIF
        END DO
      END DO
    END DO
  END IF
!      Multiply bottom flux by cell area to get volumetric rate
	   !PRINT *, DELC
	   !READ(*,*)
	   !PRINT *, DELR
	   !READ(*,*)
  DO IR=1,NROW
    DO IC=1,NCOL
      UNSFLUXV(IC,IR)=UNSFLUX(IC,IR)*DELC(IR)*DELR(IC) ! modified by Adam Szymkiewicz
	  !PRINT *, IC,IR,UNSFLUX(IC,IR),UNSFLUXV(IC,IR)
	END DO
  END DO
!      Process each horizontal cell location
!      For NUNSFOP=1, the leyer is 1
!      For NUNSFOP=2, find the uppermost active cell
  DO IR=1,NROW
    DO IC=1,NCOL
      IF(NUNSFOP.EQ.1) THEN
        IL=1
      ELSE
        DO IL=1,NLAY
          IF(IBOUND(IC,IR,IL).NE.0) EXIT
        END DO
        IF(IL.GT.NLAY) IL=1
      END IF
!        Ignore external cells
      IF(IBOUND(IC,IR,IL).GT.0) THEN
        RHS(IC,IR,IL) = RHS(IC,IR,IL)+UNSFLUXV(IC,IR) ! Adam Szymkiewicz
!		PRINT *, IC,IR,'UNSF = ',UNSFLUXV(IC,IR)
		  !READ(*,*)
      END IF
    END DO
  END DO
 ! READ(*,*)
  RETURN
END SUBROUTINE
	
SUBROUTINE GWF2UNSF1BD(IGRID,KSTP,KPER)
!     ******************************************************************
!     Calculate volumetric BuDget for unsaturated flow
!     ******************************************************************
  USE GLOBAL, ONLY:IOUT,NLAY,NROW,NCOL,BUFF,IBOUND
  USE GLOBAL, ONLY:LSTCHK
  USE GWFBASMODULE, ONLY:MSUM,ICBCFL,DELT,PERTIM,TOTIM,VBVL,VBNM
  USE GWFUNSFMODULE, ONLY:IUNSFCB,NUNSFOP,UNSFLUXV,IUNSFPR,NPUNSF,&
                              &UNSFLUX,KZON,ZONEFLUX
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: IGRID,KSTP,KPER
  DOUBLE PRECISION,PARAMETER :: ZERO=0.
  CHARACTER(LEN=*),PARAMETER :: FMT1 = "('1',/2X,'FLUX AT END OF&
       & TIME STEP',1X,I4,' IN STRESS PERIOD ',I4/2X,75('-'))"
  DOUBLE PRECISION :: RATIN,RATOUT,Q,QQ
  INTEGER,DIMENSION(3,NPUNSF+1)::UZO
  INTEGER,DIMENSION(NROW,NCOL)::IFLUX
  INTEGER :: IBD,IL,IR,IC,J,I,N,K,KK,ZCOUNT
  INTEGER :: UHFUAC
  CHARACTER(LEN=16) :: TEXT
!      Body
  CALL GWF2UNSF1PNT(IGRID)
  DATA TEXT /'      UNSAT FLOW'/
!      Clear the rate accumulator
  RATIN = ZERO
  RATOUT = ZERO
!      Set Cell-by-cell save flag (IBD) and clear the buffer
  IBD=0
  IF(IUNSFCB.GT.0) IBD=ICBCFL
  DO IL=1,NLAY
    DO IR=1,NROW
      DO IC=1,NCOL
        BUFF(IC,IR,IL)=ZERO
      END DO
    END DO
  END DO
!      Process each horizontal cell location
!      For NUNSFOP=1, the layer is 1
!      For NUNSFOP=2, find the uppermost active cell
  DO IR=1,NROW
    DO IC=1,NCOL
      IF(NUNSFOP.EQ.1) THEN
        IL=1
      ELSE
        IL=UHFUAC(IC,IR,NLAY)
        IFLUX(IC,IR)=IL
      END IF
!        Ignore external cells
      IF(IBOUND(IC,IR,IL).GT.0) THEN
        Q = UNSFLUXV(IC,IR)
        QQ= Q
        IF(Q.GT.0) THEN
          RATOUT = RATOUT+QQ
          BUFF(IC,IR,IL) = SNGL(-Q)
        ELSE IF(Q.LT.0) THEN
		  RATIN = RATIN-QQ
          BUFF(IC,IR,IL) = SNGL(-Q)
		END IF
	  END IF
	END DO
  END DO
!      If cell by cell flows should be saved, call the appropriate
!      utility module to write them
  IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IUNSFCB,BUFF,NCOL,NROW,&
                               &NLAY,IOUT)
  IF(IBD.EQ.2) CALL UBDSV3(KSTP,KPER,TEXT,IUNSFCB,BUFF,IFLUX,&
                       &NUNSFOP,NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,&
                       &IBOUND)
!      Move total unsat flow rate into VBVL for printing
  VBVL(3,MSUM) = SNGL(RATIN)
  VBVL(4,MSUM) = SNGL(RATOUT)
!      Move unsat flow for this time step to the VBVL accumulator
  VBVL(1,MSUM) = VBVL(1,MSUM) + SNGL(RATIN*DELT)
  VBVL(2,MSUM) = VBVL(2,MSUM) + SNGL(RATOUT*DELT)
!      Move the budget term label to VBNM for printing
  VBNM(MSUM) = TEXT
!      Increment the budget counter
  MSUM=MSUM+1
!     CHECK: Modified by BLe
!      Print flux for unsat at each cell
!      ZCOUNT set to zero to avoid the 'maybe-uninitialized' warning
  ZCOUNT=0
  IF(LSTCHK(3))THEN
!    DO I=1,NROW
!      DO J=1,NCOL
!        IF(I.EQ.1.AND.J.EQ.1) THEN
!          UZO(1,1)=KZON(J,I)
!          UZO(2,1)=J
!          UZO(3,1)=I
!          ZCOUNT=1
!        ELSE IF (ALL(UZO(1,1:ZCOUNT).NE.KZON(J,I))) THEN
!          ZCOUNT=ZCOUNT+1
!          UZO(1,ZCOUNT)=KZON(J,I)
!          UZO(2,ZCOUNT)=J
!          UZO(3,ZCOUNT)=I
!        END IF
!      END DO
!    END DO
    IF(IUNSFPR.GE.0) THEN
      WRITE (IOUT,FMT1) KSTP,KPER
      DO N=1,NPUNSF  !+1 ! modified by Adam Szymkiewicz
!        K=UZO(2,N)
!        KK=UZO(3,N)
        WRITE (IOUT,'(1X,I3,2X,ES13.6)') N, ZONEFLUX(N) !UZO(1,N) !,UNSFLUX(K,KK) to check
      END DO
    END IF
  END IF
  RETURN
END SUBROUTINE
	
SUBROUTINE GWF2UNSF1INI()
!     ******************************************************************
!     INI... subroutine, called by RPP
!      Sets ITCUM,TLEVEL,ALEVEL,PLEVEL,HROOT,VROOT,ITER,WCUMT,WCUMA,
!           LMINSTEP,SUMVBOT,SINK and CUMQ
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: NPUNSF,TOTIMOLD,LPRINT,LOUTF
       USE GWFUNSFMODULE, ONLY: NUMNPD
       USE GWFUNSFMODULE, ONLY: TLEVEL,ALEVEL,PLEVEL
       USE GWFUNSFMODULE, ONLY: ITCUM,WCUMT,WCUMA,CUMQ
       USE GWFUNSFMODULE, ONLY: HROOT,VROOT,ITER,LMINSTEP,SUMVBOT
       USE GWFUNSFMODULE, ONLY: SINK
       IMPLICIT NONE
!      Variables
       INTEGER :: I,L,N
       DO I=1,NPUNSF
        ITCUM(I)=0
        TLEVEL(I)=1! (changed 06/09/13) TLEVEL(I)=1
        ALEVEL(I)=1
        PLEVEL(I)=1! (changed 30/09/13) TLEVEL(I)=1
        HROOT(I)=0.
        VROOT(I)=0.
        ITER(I)=0
        WCUMT(I)=0.
        WCUMA(I)=0.
        LMINSTEP(I)=0
        SUMVBOT(I)=0.0
       END DO
       DO I=1,NPUNSF
        DO N=1,NUMNPD
         SINK(N,I)=0.
        END DO
       END DO
       DO I=1,NPUNSF
        DO L=1,12
         CUMQ(L,I)=0.
        END DO
       END DO
       LPRINT=LOUTF
       TOTIMOLD=0.0
	END SUBROUTINE
	
SUBROUTINE GWF2UNSF1INIT(IN,INEVT,INRCH)
!     ******************************************************************
!     INIT... subroutine, called by RP
!      Sets NIZ,KZON,IZ,MAXIT,TOLTH,TOLH
!     ******************************************************************
  USE GLOBAL, ONLY:IOUT,NROW,NCOL
  USE GLOBAL, ONLY:LSTCHK
  USE PARAMMODULE, ONLY: IZON,NZONAR,ZONNAM
!     CHECK: NZONAR replaced NZN in MF2005
  USE GWFUNSFMODULE, ONLY: IZ,KZON,NIZ
  USE GWFUNSFMODULE, ONLY: MAXIT,TOLTH,TOLH
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: IN,INEVT,INRCH
  CHARACTER(LEN=*),PARAMETER :: WARN1 =&
       & '----- WARNING - EVAPORATION PACKAGE EXISTS! ----'
  CHARACTER(LEN=*),PARAMETER :: WARN2 =&
       & '----- WARNING - RECHARGE PACKAGE EXISTS! -----'
  CHARACTER(LEN=*),PARAMETER :: FMT1 ="(' ABS ',A,' TOLERANCE IN&
       & UNSAT FLOW REGION = ',F8.5)"
  CHARACTER(LEN=100):: LINE
  CHARACTER(LEN=10) :: ZONARR,ZONAME,ZONARR0
  INTEGER :: I,J,K,NTOT
  LOGICAL:: LZFOUND
!      Body
!      Warn the user if either the evaporation package or the recharge
!      package is active
  IF(LSTCHK(2))THEN
!     CHECK: Consider using UPKGS (see utl7ext.f90)
    IF (INEVT.GT.0) WRITE(IOUT,'(//A/)') WARN1
    IF (INRCH.GT.0) WRITE(IOUT,'(//A/)') WARN2
  END IF
!      Read Zone Array
  IF(LSTCHK(3)) WRITE(IOUT,'(//A/)') 'UNSF -- READ ZONE ARRAY'
  READ(IN,'(A)') ZONARR
  CALL UPCASE(ZONARR)
  DO I=1, 10  ! modified by Adam Szymkiewicz, get rid of trailing tabs
	ZONARR0(I:I)=CHAR(32)
  END DO
  J=1
  DO I=1, LEN(ZONARR)
	IF ((ICHAR(ZONARR(I:I)).GE.32).AND.(ICHAR(ZONARR(I:I)).LE.126)) THEN
	  ZONARR0(J:J)=ZONARR(I:I)
	  J=J+1
	END IF
  END DO
  ZONARR=ZONARR0
  
  IF (ZONARR.EQ.'ALL') THEN
    IF (LSTCHK(3)) WRITE(IOUT,*)' ALL CELLS ARE ASSOCIATED WITH ONE PROFILE'
	NIZ=0
  ELSE IF(ZONARR.EQ.'EACH') THEN
    NIZ=2
    NTOT=NROW*NCOL
    DO I=1,NTOT
      IZ(I)=I
    END DO
  ELSE
    NIZ=1
!       Read zone names to find ZONEARR and copy the
!       correct values into KZON
	LZFOUND=.FALSE.
    DO K=1,NZONAR
      ZONAME=ZONNAM(K)
      CALL UPCASE(ZONAME)
      IF (TRIM(ZONAME).EQ.TRIM(ZONARR)) THEN
		LZFOUND=.TRUE.
        IF(LSTCHK(3)) WRITE(IOUT,*) 'ZONE ARRAY NAME = ',ZONARR
		DO I=1,NROW
          DO J=1,NCOL
            KZON(J,I)=IZON(J,I,K)
          END DO
        END DO
	  END IF
	END DO
	IF (.NOT.(LZFOUND)) STOP 'ZONE ARRAY NOT DEFINED'
!    CALL UHRCOM(IN,IOUT,LINE)
!    BACKSPACE IN
  END IF
!     CHECK: Consider using URWORD
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
  READ (IN,*) MAXIT,TOLTH,TOLH
  IF (LSTCHK(3)) THEN
    WRITE(IOUT,*)
    WRITE(IOUT,*) 'MAXIMUM ITERATION = ',MAXIT
    WRITE(IOUT,FMT1) 'WATER CONTENT',TOLTH
    WRITE(IOUT,FMT1) 'PRESSURE HEAD',TOLH
  END IF
  RETURN
END SUBROUTINE
  
SUBROUTINE GWF2UNSF1TMIN(IN,KPER)
!     ******************************************************************
!     TiMe INistialization subroutine
!      Sets TINIT,TMAX,TOLD,DT,DTMIN,DTMAX,DTOPT,ITMIN,ITMAX,TPRINT
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:TINIT,TMAX,TOLD,DT,DTMIN,DTMAX,DTOPT,T
       USE GWFUNSFMODULE, ONLY:ITMIN,ITMAX,TPRINT
       USE GWFUNSFMODULE, ONLY:NPUNSF,DMUL,DMUL2
       USE GLOBAL, ONLY:NSTP,PERLEN,NPER,TSMULT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: IN,KPER
       INTEGER :: I,J
!     CHECK: DMUL,DMUL2
       DOUBLE PRECISION :: DDT,DDTMIN,DDTMAX,TTMAX
!      Body
       TTMAX=PERLEN(1)/FLOAT(NSTP(KPER))
       IF(TSMULT(1).NE.1.) THEN
!     CHECK: Should this DOUBLE PRECISIONly be TSMULT(1) ?
        TTMAX=PERLEN(1)*(1.-TSMULT(1))/(1.-(TSMULT(1))**NSTP(KPER))
       END IF
       READ(IN,*) DDT,DDTMIN,DDTMAX,DMUL,DMUL2,ITMIN,ITMAX
       DO I=1,NPUNSF
        TINIT(I)=0.0
        TMAX(I)=TTMAX
        DT(I)=DDT
        DTMIN(I)=DDTMIN
        DTMAX(I)=DDTMAX
        DTOPT(I)=DT(I)
        TOLD(I)=TINIT(I)
        T(I)=TINIT(I)+DT(I)
       END DO
       TPRINT=0.0
       DO J=1,NPER
        TPRINT=TPRINT+PERLEN(J)
       END DO
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1MATIN(IN)
!     ******************************************************************
!     Read MATerial INformation
!      Sets values for HTAB(1),HTAB(NTAB),IMODEL and PAR
!     ******************************************************************
       USE GLOBAL, ONLY: IOUT,LSTCHK,LENUNI
       USE GWFUNSFMODULE, ONLY: NMAT,XCONV,IUNSFPR,IMODEL,NTAB
       USE GWFUNSFMODULE, ONLY: PAR,HTAB,IMODEL
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: IN
       DOUBLE PRECISION,PARAMETER :: EPS = -0.00001
       INTEGER,PARAMETER :: NPAR=6
       CHARACTER(LEN=*),PARAMETER :: FMT1 = &
         &"(/'MATNUM, PARAM. ARRAY:'//'   MAT    THR    THS    ',&
            &'  ALPHA          N           KS        L'/)"
       CHARACTER(LEN=*),PARAMETER :: FMT2 = &
                                   &'(I5,2X,2F7.4,3E12.3,4F7.3,2E12.3)'
       CHARACTER(LEN=100) :: LINE
       INTEGER :: I,M
       DOUBLE PRECISION :: HTAB1,HTABN,RHENTRY
!      Body
       READ(IN,*) HTAB1,HTABN,IMODEL
       XCONV = 1.0D0
!      Convert unit (LENUNI==3 -> cm , LENUNI==1 -> ft , LENUNI==2 -> m)
       IF(LENUNI.EQ.3) XCONV=100.
       IF(LENUNI.EQ.1) XCONV=1./0.3048
       RHENTRY = 0.2*XCONV
!      MODEL TYPES: 0 - Van Genuchten
!                   1 - Van Genuchten with Air entry value of 2cm
!                   2 - Brooks and Corey
       HTAB1=-DMIN1(ABS(HTAB1),ABS(HTABN))
       HTABN=-DMAX1(ABS(HTAB1),ABS(HTABN))
       IF((HTAB1.GT.EPS .AND. HTABN.GT.EPS) .OR. HTAB1.EQ.HTABN) THEN
        HTAB1 = -0.0001*XCONV
        HTABN = -100.  *XCONV
       END IF
!      ------
       IF(IUNSFPR.EQ.0) THEN
        IF(IMODEL.EQ.2 .OR. IMODEL.EQ.1 .OR. IMODEL.EQ.0) THEN
         IF(LSTCHK(3)) WRITE(IOUT,FMT1)
        END IF
       END IF
!      ------
       CALL UHRCOM(IN,IOUT,LINE)
       BACKSPACE IN
       DO M=1,NMAT
        READ(IN,*) (PAR(I,M), I=1,NPAR)
        IF(IMODEL.EQ.1)THEN
         PAR(7,M) = PAR(1,M)+(PAR(2,M)-PAR(1,M))*&
                  &(1.+(PAR(3,M)*RHENTRY)**PAR(4,M)**(1.-1./PAR(4,M)))
        END IF
        IF(IUNSFPR.EQ.0) THEN
         IF(LSTCHK(3)) WRITE(IOUT,FMT2) M,(PAR(I,M),I=1,NPAR)
        END IF
       END DO
       HTAB(1)   =HTAB1
       HTAB(NTAB)=HTABN
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1GENMAT()
!     ******************************************************************
!     GENerate hydraulic (MATerial) properties
!      Sets values for HSAT,CONTAB,CAPTAB,THETAB,QE,A10H,A10H
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:NTAB,NMAT,IMODEL,PAR,HSAT,HTAB
       USE GWFUNSFMODULE, ONLY:CONTAB,CAPTAB,THETAB
       IMPLICIT NONE
!      Variables
       INTEGER :: M,J,L
       DOUBLE PRECISION :: HTAB1,HTABN,ALH,DLH,FH,QE,A10H,A10K,FK,FC1,FQ,FS
!      Body
       HTAB1=HTAB(1)
       HTABN=HTAB(NTAB)
       DLH = (DLOG10(-HTABN)-DLOG10(-HTAB1))/(NTAB-1)
       DO J=1,NTAB
        ALH=DLOG10(-HTAB1)+(J-1)*DLH
        HTAB(J)=-10**ALH
       END DO
       DO M=1,NMAT
!     CHECK: This subroutine only affects the first profile.
! modified by Adam Szymkiewicz, 23.03.2016
        HSAT(M) = FH(IMODEL,1.0d0,PAR(:,M))
        DO L=1,NTAB
         CONTAB(L,M)=FK(IMODEL,HTAB(L),PAR(:,M))
         CAPTAB(L,M)=FC1(IMODEL,HTAB(L),PAR(:,M))
         THETAB(L,M)=FQ(IMODEL,HTAB(L),PAR(:,M))
         QE         =FS(IMODEL,HTAB(L),PAR(:,M))
         A10H=DLOG10(DMAX1(-HTAB(L),1d-30))
         A10K=DLOG10(CONTAB(L,M))
        END DO
	   END DO
	   !DO J=1, NTAB
	   !  PRINT *, HTAB(J), THETAB(J,1), CONTAB (J,1)
	   !END DO
	   !PRINT *, FQ(IMODEL,-7.99D0,PAR(:,2))
	   !DO M=1,NMAT
	!	   PRINT *, HSAT(M), HSAT(M)
	   !END DO
	   !READ(*,*)
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1BASINF(PID,IN)
!     ******************************************************************
!     get BASic INFormation for selected profile (PID)
!     ******************************************************************
       USE GLOBAL, ONLY: IOUT,LSTCHK
       USE GWFUNSFMODULE, ONLY: INFTOP,INFBOT,KODTOP,LWAT,SINKF,SHORTO,&
                               &ATMBC,WLAYER,LINITW,NOLAY,NOMAT,PERIMP,&
                               &RTOP,RROOT,RBOT,HCRITA,KTOLD,&
                               &KBOLD,KODBOT,NIZ,IZ,ICOUP
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID,IN
       CHARACTER(LEN=*),PARAMETER :: FMT1 = &
        &"('INITIAL CONDITION IS GIVEN IN TERMS OF THE ',A)"
       CHARACTER(LEN=*),PARAMETER :: FMT2 = &
        &"('NUMBER OF ',A,' IN PROFILE ',I4,' = ',I4)"
!       CHARACTER(LEN=17)  :: EXTD
       CHARACTER(LEN=100) :: LINE
!      Body
!     CHECK: Several lines moved to MATIN, since RHENTRY is used there
!            See line 2497 of the original code for more information
       WRITE(IOUT,*)
       WRITE(IOUT,*) '=============================================='
       WRITE(IOUT,*) '             PROFILE ',PID
       WRITE(IOUT,*) '=============================================='
!      Initialize some values
       IF(NIZ.EQ.1) READ(IN,*) IZ(PID)
	   PRINT *, PID, IZ(PID)
       CALL UHRCOM(IN,IOUT,LINE)
       BACKSPACE IN

       READ(IN,*) SINKF(PID),WLAYER(PID),LINITW(PID),NOMAT(PID) !, &
               !  &PERIMP(PID)
       PERIMP(PID) = 0 ! modified by Adam Szymkiewicz
	   ATMBC(PID) =0
       LWAT(PID)  =0
       SHORTO(PID)=-1
       INFTOP(PID)=1
       IF (ICOUP.LE.0) THEN
	     INFBOT(PID)=-1 !modified by Adam Szymkiewicz, impermeable bottom
         KODBOT(PID)=-3  !modified by Adam Szymkiewicz
	   ELSE
	     INFBOT(PID)=1 !modified by Adam Szymkiewicz, variable head at the bottom
         KODBOT(PID)=3 !modified by Adam Szymkiewicz
	   END IF
 	   KODTOP(PID)=-4
       NOLAY(PID) =1
!     CHECK: Removed lines that had no effect.
!            See line 2522 of the original code for more information
       RROOT(PID) =0.0
       RTOP(PID)  =0.0
       RBOT(PID)  =0.0
!     CHECK: Moved input code from the bottom of the subroutine to here.
!            Most of them were deemed unnecessary and removed.
!            See line 2630 of the original code for more information
       HCRITA(PID)= -ABS(1.E+10)
       IF(WLAYER(PID).GE.0) KODTOP(PID)=-IABS(KODTOP(PID))
       KTOLD(PID) = KODTOP(PID)
       KBOLD(PID) = KODBOT(PID)
!      Write the Retrieved information
!     CHECK: Since most of this information is defined directly instead
!            of being read from IN, most of it was not outputted in the
!            original code. Nevertheless, the code is necessary code was
!            present in the form of comments between lines 2531 and 2605
!            Refer to these lines when extending this code.
       IF(LSTCHK(3))THEN
        IF(SINKF(PID).GE.0) THEN
         WRITE(IOUT,*) 'WATER EXTRACTION FROM THE ROOT ZONE OCCURS.'
        END IF
        IF(LINITW(PID).GE.0) THEN
         WRITE(IOUT,FMT1) 'WATER CONTENT'
        ELSE
         WRITE(IOUT,FMT1) 'PRESSURE HEAD'
        END IF
        IF(WLAYER(PID).GE.0) THEN
         WRITE(IOUT,*) 'WATER CAN ACCUMULATE AT THE SURFACE WITH ZERO',&
                      &' SURFACE RUNOFF.'
        END IF
        WRITE(IOUT,*)
        WRITE(IOUT,FMT2) 'SOIL MATERIALS',PID,NOMAT(PID)
!     CHECK: Commented in the original code, because it is set to 1
!            See line 2627 of the original code for more information
        WRITE(IOUT,FMT2) 'SUBREGIONS',PID,NOLAY(PID)
       END IF
       RETURN
	END SUBROUTINE
	
SUBROUTINE GWF2UNSF1NODINF(PID,IN)
!     ******************************************************************
!     get NODal INFormation for selected profile (PID)
!     ******************************************************************
  USE GLOBAL, ONLY: IOUT,LSTCHK
  USE GWFUNSFMODULE, ONLY: NUMNPD,IUNSFPR
  USE GWFUNSFMODULE, ONLY: NUMNP,HTOP,HBOT,XSURF
  USE GWFUNSFMODULE, ONLY: X,HHNEW,HHOLD,MATNUM,HTEMP,LAYNUM,BETA
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: PID,IN
  CHARACTER(LEN=*),PARAMETER::FMT1="(/'NODAL POINT INFORMATION'//&
       &'NODE      X         HOLD    MATN  BETA'/)"
  CHARACTER(LEN=*),PARAMETER::FMT2="(I4,2F11.3,I5,F8.3,E12.4)"
  CHARACTER(LEN=100) :: LINE
  INTEGER :: I,J,M,N,NOLD
  DOUBLE PRECISION :: H,B,X1,DX,SBETA,SHOLD
!      Body
  WRITE(IOUT,'(//A/)') 'READING NODAL INFORMATION'
  READ(IN,*) NUMNP(PID)
  IF (NUMNP(PID).GT.NUMNPD) THEN
    IF (LSTCHK(1))WRITE(IOUT,*)'*** ERROR: CHECK THE NUMBER OF NODES'
!     CHECK: Replaced STOP by USTOP to better conform to MODFLOW
    CALL USTOP(' ')
  END IF
  CALL UHRCOM(IN,IOUT,LINE)
  BACKSPACE IN
  J=NUMNP(PID)+1
!     CHECK: Added the following line to prevent ambiguity
  NOLD=0
  DO WHILE(J.GT.1)
    J=J-1
    READ(IN,*) N,X1,H,M,B
    N=NUMNP(PID)-N+1
    X(N,PID)=X1
    HHOLD(N,PID)=H
    MATNUM(N,PID)=M
	BETA(N,PID)=B
!       LAYNUM is no longer an input
    LAYNUM(N,PID)=1
        IF(J-N.LT.0)THEN
         IF(LSTCHK(1))WRITE(IOUT,*)'ERROR IN NODINF AT NODE =', N
         CALL USTOP(' ')
        ELSE IF(J-N.GT.0) THEN
         DX=X(NOLD,PID)-X(N,PID)
         SHOLD=(HHOLD(N,PID)-HHOLD(N,PID))/DX
         SBETA=(BETA(NOLD,PID)-BETA(N,PID))/DX
         DO I=NOLD-1,N+1,-1
          DX=X(NOLD,PID)-X(N,PID)
          HHOLD(I,PID)=HHOLD(NOLD,PID)-SHOLD*DX
          BETA (I,PID)=BETA (NOLD,PID)-SBETA*DX
          MATNUM(I,PID)=MATNUM(I+1,PID)
          LAYNUM(I,PID)=LAYNUM(I+1,PID)
         END DO
         J=N
        END IF
        J=N
        NOLD=N
       END DO
       SBETA=BETA(NUMNP(PID),PID) * &
                           &(X(NUMNP(PID),PID)-X(NUMNP(PID)-1,PID)) /2.
       DO I=2,NUMNP(PID)-1
        SBETA=SBETA+BETA(I,PID)*(X(I+1,PID)-X(I-1,PID))/2.
       END DO
       DO I=2,NUMNP(PID)
        IF(SBETA.GT.0) THEN
         BETA(I,PID)=BETA(I,PID)/SBETA
        ELSE
         BETA(I,PID)=0.
        END IF
       END DO
       XSURF(PID)=X(NUMNP(PID),PID)
!      Print Nodal Information
       IF(IUNSFPR.EQ.0 .AND. LSTCHK(3)) THEN
        WRITE(IOUT,FMT1)
        DO N=NUMNP(PID),1,-1
         WRITE(IOUT,FMT2) NUMNP(PID)-N+1,X(N,PID),HHOLD(N,PID),&
                         &MATNUM(N,PID),BETA(N,PID)
         HHNEW(N,PID) = HHOLD(N,PID)
         HTEMP(N,PID) = HHOLD(N,PID)
        END DO
       ELSE
        DO N=NUMNP(PID),1,-1
         HHNEW(N,PID) = HHOLD(N,PID)
         HTEMP(N,PID) = HHOLD(N,PID)
        END DO
       END IF
  HBOT(PID)=HHNEW(1,PID)
  HTOP(PID)=HHNEW(NUMNP(PID),PID)
!  PRINT *, HBOT(PID),HTOP(PID)
!  READ(*,*)
  IF (LSTCHK(3)) THEN
    WRITE(IOUT,*)
  END IF
  RETURN
END SUBROUTINE

SUBROUTINE GWF2UNSF1INITW(PID)
!     ******************************************************************
!     INITial condition given in terms of Water content
!      Sets: HHNEW,HHOLD,HTEMP,HBOT,HTOP for a specific profile
!     ******************************************************************
       USE GLOBAL, ONLY:IOUT,LSTCHK
       USE GWFUNSFMODULE, ONLY: IMODEL
       USE GWFUNSFMODULE, ONLY: NUMNP,MATNUM,HHNEW,HHOLD,HTEMP,PAR
       USE GWFUNSFMODULE, ONLY: HTOP,HBOT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I,M
       DOUBLE PRECISION :: QE,FH
!      Body
       DO I=1,NUMNP(PID)
        M = MATNUM(I,PID)
        QE= MIN((HHNEW(I,PID)-PAR(1,M))/(PAR(2,M)-PAR(1,M)),1.)
        IF(QE.LT.0) THEN
         IF(LSTCHK(1))THEN
          WRITE(IOUT,*)'*** ERROR: INITIAL WATER CONTENT CONDITION ',&
                      &'IS LOWER THAN QR !'
         END IF
         CALL USTOP(' ')
        END IF
        HHNEW(I,PID)=FH(IMODEL,QE,PAR(:,M))
        HHOLD(I,PID)=HHNEW(I,PID)
        HTEMP(I,PID)=HHNEW(I,PID)
       END DO
       HBOT(PID)=HHNEW(1,PID)
       HTOP(PID)=HHNEW(NUMNP(PID),PID)
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SINKIN(PID,IN)
!     ******************************************************************
!     Read SINK INformation
!      Sets: P0,P2H,P2L,P3,R2H,R2L,POPTM
!     ******************************************************************
       USE GLOBAL, ONLY: IOUT,LSTCHK
       USE GWFUNSFMODULE, ONLY: P0,P2H,P2L,P3,R2H,R2L,POPTM
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID,IN
!      Body
       IF (LSTCHK(3)) THEN
        WRITE(IOUT,'(//A/)') 'READING SINK INFORMATION'
       END IF
       READ(IN,*) P0(PID),P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID),&
                 &POPTM(PID)
       P0(PID) =-ABS(P0(PID))
       P2L(PID)=-ABS(P2L(PID))
       P2H(PID)=-ABS(P2H(PID))
       P3(PID) =-ABS(P3(PID))
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SETBC(PID,IN)
!     ******************************************************************
!     SET Boundary Conditions
!      Sets: MAXAL,HCRITS,TTATM,PRECI,RSO,RR,HCA,TATM,PREC,RSOIL
!            RTOP,LMINSTEP,RROOT,HBOT
!     ******************************************************************
       USE GLOBAL, ONLY: IOUT,LSTCHK
       USE GWFUNSFMODULE, ONLY:TATM,PREC,RSOIL,RR,HCA,RB,HB,HT,RTOP,&
                              &RROOT,RBOT,HCRITA,HBOT,HTOP,GWL0L,&
                              &INFTOP,INFBOT,KODTOP,LMINSTEP,&
                              &HCRITS,TTATM,PRECI,RSO,MAXAL,&
                              &HHOLD
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID,IN
       CHARACTER(LEN=100) :: LINE
       INTEGER :: I
       DOUBLE PRECISION :: RTOPOLD
!      Body
       GWL0L=0.0
       RB=0.0
       READ (IN,*)MAXAL(PID),HCRITS(PID)
       CALL UHRCOM(IN,IOUT,LINE)
       BACKSPACE IN
       DO I=1,MAXAL(PID)
!     CHECK: These values are read, but only the first of each series
!            is used
        READ(IN,*)TTATM(PID,I),PRECI(PID,I),RSO(PID,I),RR(PID,I),&
                 &HCA(PID,I)
        HT(PID,I)=0.
       END DO
       HB(PID)=HHOLD(1,PID)
       TATM(PID)=TTATM(PID,1)
       PREC(PID)=PRECI(PID,1)
       RSOIL(PID)=RSO(PID,1)
!      Top of the profile
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
!      Bottom of the profile
       IF(INFBOT(PID).GE.0) THEN
        IF(ABS(RBOT(PID)-RB).GT.ABS(RBOT(PID))*0.2) LMINSTEP=0
        RBOT(PID)=RB
        IF(ABS(HBOT(PID)-HB(PID)-GWL0L).GT.ABS(HBOT(PID))*0.2) THEN
         LMINSTEP=0
        END IF
        HBOT(PID)=HB(PID)+GWL0L
       END IF
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SETBC1(PID)
!     ******************************************************************
!     SET Boundary Conditions (Modified version used in AD)
!      Sets: MAXAL,HCRITS,TTATM,PRECI,RSO,RR,HCA,TATM,PREC,RSOIL
!            RTOP,LMINSTEP,RROOT,HBOT
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:TATM,PREC,RSOIL,RR,HCA,RB,HB,HT,RTOP,&
                              &RROOT,RBOT,HCRITA,HBOT,HTOP,GWL0L,&
                              &INFTOP,INFBOT,KODTOP,LMINSTEP,&
                              &TTATM,PRECI,RSO,MAXAL
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I
       DOUBLE PRECISION :: RTOPOLD=0.
!      Body
       DO I=1,MAXAL(PID)-1
!     CHECK: The Loop appears unnecessary, since it only enters the
!            following IF statement once before exiting
        IF(TATM(PID).EQ.TTATM(PID,I)) THEN
         TATM(PID)=TTATM(PID,I+1)
         PREC(PID)=PRECI(PID,I+1)
         RSOIL(PID)=RSO(PID,I+1) ! modified by Adam Szymkiewicz
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
         EXIT
        END IF
       END DO
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SETSNK(PID)
!     ******************************************************************
!     SET SiNKs
!      WARNING: THETA should be linked before doing this
!       (use THETA=>THOLD or THETA=>THNEW)
!      Sets: SINK,HROOT,VROOT
!     ******************************************************************
!     CHECK: TPOT -> RROOT (see lines 575,2097 of the original code)
       USE GWFUNSFMODULE, ONLY:NMAT
       USE GWFUNSFMODULE, ONLY:NUMNP,HROOT,VROOT,RROOT,P0,POPTM,P2H,&
                              &P2L,P3,R2H,R2L,DT
       USE GWFUNSFMODULE, ONLY:MATNUM,X,SINK,HHNEW,BETA,THETA,PAR
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I,M,N
       DOUBLE PRECISION :: AROOT,DXM,ALFA,FALFA
!      Body
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
         ALFA=FALFA(RROOT(PID),HHNEW(I,PID),P0(PID),POPTM(PID),&
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
       END DO
       IF(AROOT.GT.0.001)HROOT(PID)=HROOT(PID)/AROOT
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SUBREG(PID)
!     CHECK: This function barely serves any purpose, as most of it is
!            commented in the original code. The main purpose of this
!            subroutine appears to have been printing (before that was
!            removed). The non-commented parts have been left, despite
!            serving no purpose.
!            See line 1879 of the original code for more information
!     ******************************************************************
!      Sets: WATIN, if LWAT(PID)>=0 and PLEVEL(PID)==0
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: NUMNP,NOLAY,HHNEW,THETA,THOLD,&
                               &X,MATNUM,LAYNUM,DT,CON,PLEVEL,&
                               &WATIN,LWAT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I,J,N,MI,MJ,LAY
       DOUBLE PRECISION    :: ATOT,VOLUME,CHANGE,HTOT,DELTW,HE,VNEWI,VOLDI,DX,DX1,&
                 &V1,DXN,VN,DH1,DHN
!     CHECK: Consider removing the hard limit on array size
       DOUBLE PRECISION,DIMENSION(10) :: AREA,SUBVOL,SUBCHA,HMEAN
!      Body
       N=NUMNP(PID)
       ATOT=0.
       HE=0.
       IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
        VOLUME=0.
        CHANGE=0.
        HTOT  =0.
        DELTW =0.
       END IF
!     CHECK: NOLAY(PID) has been set to 1 in GWF2UNSF1BASINF, consider
!            removing the loop if there are no plans to restore this
!            functionality
       DO LAY=1,NOLAY(PID)
        AREA(LAY)=0.
        IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
         SUBVOL(LAY)=0.
         SUBCHA(LAY)=0.
         HMEAN (LAY)=0.
        END IF
       END DO
       DO I=N-1,1,-1
        J=I+1
        MI=MATNUM(I,PID)
        MJ=MATNUM(J,PID)
        LAY=LAYNUM(J,PID)
        DX=X(J,PID)-X(I,PID)
        AREA(LAY)=AREA(LAY)+DX
        ATOT=ATOT+DX
        IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
         HE=(HHNEW(I,PID)+HHNEW(J,PID))/2.
         VNEWI=DX*(THETA(I,PID)+THETA(J,PID))/2.
         VOLDI=DX*(THOLD(I,PID)+THOLD(J,PID))/2
         VOLUME=VOLUME+VNEWI
         CHANGE=CHANGE+(VNEWI-VOLDI)/DT(PID)
         SUBCHA(LAY)=SUBCHA(LAY)+(VNEWI-VOLDI)/DT(PID)
         SUBVOL(LAY)=SUBVOL(LAY)+VNEWI
         HTOT=HTOT+HE*DX
         HMEAN(LAY)=HMEAN(LAY)+HE*DX
        END IF
        IF (LWAT(PID).GE.0) THEN
         IF (PLEVEL(PID).EQ.0) THEN
          WATIN(I,PID)=INT(VNEWI)
         ELSE
          DELTW=DELTW+ABS(WATIN(I,PID)-VNEWI)
         END IF
        END IF
       END DO
       DO LAY=1,NOLAY(PID)
        IF (AREA(LAY).GT.0) THEN
         IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
          HMEAN(LAY)=HMEAN(LAY)/AREA(LAY)
         END IF
        END IF
       END DO
!     CHECK: Verify if this is correct as operator priority is easy
!            to misunderstand without parentheses.
!            (A .AND. B .OR. C) <=> ((A.AND.B) .OR. C)
!            See line 1940 of the original code for more information
       IF(LWAT(PID).GE.0 .AND. ATOT.GT.0. .OR. PLEVEL(PID).EQ.0) THEN
        HTOT=HTOT/ATOT
       END IF
       DX1=X(2,PID)-X(1,PID)
       DH1=HHNEW(2,PID)-HHNEW(1,PID)
       V1=-(CON(1,PID)+CON(2,PID))/2.*(DH1/DX1+1.)
       DXN=X(N,PID)-X(N-1,PID)
       DHN=HHNEW(N,PID)-HHNEW(N-1,PID)
       VN=-(CON(N,PID)+CON(N-1,PID))/2.*(DHN/DXN+1.)
!     CHECK: The original code has 54 lines of commented code here
!            See line 1946 of the original code for more information
       RETURN
	END SUBROUTINE
	
SUBROUTINE GWF2UNSF1GWD(PID,KSTP,KPER,ADEPTH,AHEAD)
!     ******************************************************************
!     Calculate average GW Depth and average Head
!     Input KPER=0 for the initialization
!     ******************************************************************
  USE GLOBAL, ONLY: IOUT,BOTM,NCOL,NROW,NLAY,HNEW,IBOUND
  USE GLOBAL, ONLY: LSTCHK
  USE GWFBASMODULE, ONLY: HDRY
  USE GWFUNSFMODULE, ONLY: KZON,NIZ,IZ,HBOT,HB,X,XSURF,T,DT
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN) :: KPER,KSTP,PID
  DOUBLE PRECISION,INTENT(OUT) :: ADEPTH,AHEAD
!      min depth
  DOUBLE PRECISION,PARAMETER :: MNDPTH = 0.005
  CHARACTER(LEN=*),PARAMETER :: FMT1= "('PRESSURE HEAD AT THE&
       & BOTTOM OF PROFILE = ',F13.5,4X,'FOR STRESS PERIOD : ',I5,4X,&
       &'TIME STEP: ',I5)"
  CHARACTER(LEN=*),PARAMETER :: FMT2= "('PRESSURE HEAD AT THE&
       & BOTTOM OF PROFILE = ',F13.5,4X,'AT COLUMN',I5,4X,'ROW',I5,&
       &'FOR STRESS PERIOD : ',I5,4X,'TIME STEP: ',I5)"
  CHARACTER(LEN=*),PARAMETER :: FMT3= "('PRESSURE HEAD AT THE&
       & BOTTOM OF PROFILE = ',F13.5,4X,' FOR ZONE NUMBER : ',I3,4X,&
       &'FOR STRESS PERIOD : ',I5,4X,'TIME STEP: ',I5)"
  CHARACTER(LEN=*),PARAMETER :: FMT4= "&
             &('ADEPTH=',F13.2,4X, ' STDDEPTH=',F13.2,4X,&
             &' SKWDEPTH=',F13.2,4X,'AHEAD=',F13.2,4X,&
             &' STDHEAD=',F13.2,4X, ' SKWHEAD=',F13.2,4X,&
             &' FOR ZONE ',I3,4X)"
  INTEGER :: I,J,K,ICOUNT
  DOUBLE PRECISION :: SUM1,SUM2,SUM3,SUM4,SUM5,SUM6,X1
  DOUBLE PRECISION :: STDDEPTH,STDHEAD,SKWDEPTH,SKWHEAD
  DOUBLE PRECISION,DIMENSION(NCOL,NROW):: TTOP,NHED,DEPTH
!      Body
       !IF(KPER.GT.0)T(PID) = T(PID)+DT(PID) !Adam Szymkiewicz
!     CHECK: This should be correct
  X1=X(1,PID)
!      Read top elevation
  DO I=1,NROW
    DO J=1,NCOL
      DO K=1,NLAY
        IF(IBOUND(J,I,K).NE.0)THEN
          TTOP(J,I)=BOTM(J,I,K-1)
          EXIT
        END IF
      END DO
    END DO
  END DO
!      Read head and calculate GW depth
!     CHECK: the LAYER now varies slowly as opposed to the original code
!            See line 1728 of the original code for more information
  DO K=NLAY,1,-1 ! check - Adam Szymkiewicz
    DO I=1,NROW
      DO J=1,NCOL
        IF(HNEW(J,I,K).EQ.HDRY) THEN
          IF(K.EQ.1) THEN
            NHED(J,I) = BOTM(J,I,K)
          END IF
        ELSE
           NHED(J,I) = HNEW(J,I,K)
        END IF
      END DO
    END DO
  END DO
  DEPTH = TTOP - NHED
  !PRINT *, KZON
  !READ(*,*)
!      Calculate average GW depth and head for the zone specified by PID
  SUM1=0.0
  SUM2=0.0
  SUM3=0.0
  SUM4=0.0
  SUM5=0.0
  SUM6=0.0
  ICOUNT=0
  IF(NIZ.EQ.0) THEN
    DO I=1,NROW
      DO J=1,NCOL
        SUM1=SUM1+DEPTH(J,I)
        ICOUNT=ICOUNT+1
      END DO
    END DO
    ADEPTH=SUM1/ICOUNT
    HB(PID)=XSURF(PID)-ADEPTH-X1
    IF(KPER.GT.0) HBOT(PID)=HB(PID)
	IF(LSTCHK(3)) WRITE (IOUT,FMT1) HB(PID),KPER,KSTP
  ELSE IF(NIZ.EQ.2) THEN
!     CHECK: removed the loop, as it was deemed unnecessary since the
!            correct indexes can be found mathematically
!            see line 1787 of the original code for more information
    I = IZ(PID)/NROW
    J = MOD(IZ(PID),NCOL)
    HB(PID) = XSURF(PID)-DEPTH(J,I)-X1
    IF(KPER.GT.0)HBOT(PID)=HB(PID)
    IF(LSTCHK(3)) WRITE(IOUT,FMT2) HB(PID),J,I,KPER,KSTP
  ELSE
    DO I=1,NROW
      DO J=1,NCOL
        IF(KZON(J,I).EQ.IZ(PID)) THEN
!     CHECK: Removed DEPTHTMP to make the program shorter and clearer
!            see line 1805 of the original code for more information
          IF(DEPTH(J,I).LT.0) THEN
            SUM1 = SUM1+MNDPTH
          ELSE
            SUM1 = SUM1+DEPTH(J,I)
          END IF
          SUM2=SUM2+NHED(J,I)
          ICOUNT=ICOUNT+1
		END IF
	  END DO
	END DO
    ADEPTH=SUM1/ICOUNT
    AHEAD=SUM2/ICOUNT
    HB(PID)= XSURF(PID)-ADEPTH-X1
		!PRINT *, HB(PID)
		!READ(*,*)
    IF(KPER.GT.0) HBOT(PID)=HB(PID)
!     CHECK: Added by BLe on 11/12/12
!              to account for seepage and compute std and skewness head
	IF(LSTCHK(3))THEN
      DO I=1,NROW
        DO J=1,NCOL
          IF(KZON(J,I).EQ.IZ(PID)) THEN
            SUM3=SUM3+DEPTH(J,I)**2.-(ADEPTH**2.)
            SUM4=SUM4+NHED(J,I)**2.-AHEAD**2.
            SUM5=SUM5+(DEPTH(J,I)-ADEPTH)**3.
            SUM6=SUM6+(NHED(J,I)-AHEAD)**3.
          END IF
        END DO
      END DO
      STDDEPTH=SUM3/ICOUNT
      STDHEAD =SUM4/ICOUNT
      SKWDEPTH=(SUM5/ICOUNT)/(STDDEPTH)**1.5
      SKWHEAD =(SUM6/ICOUNT)/(STDHEAD)**1.5
      WRITE(IOUT,FMT3)HB(PID),IZ(PID),KPER,KSTP
!      WRITE(IOUT,FMT4)ADEPTH,STDDEPTH,SKWDEPTH,AHEAD,STDHEAD,&
!                        &SKWHEAD,IZ(PID)
    END IF
  END IF
  
!  PRINT *, PID, ICOUNT, HB(PID), HBOT(PID)
!  PRINT *, XSURF(PID), ADEPTH, X1, IZ(PID)
!  READ(*,*)
  RETURN
END SUBROUTINE
	
SUBROUTINE GWF2UNSF1WF(PID)
!     ******************************************************************
!     Water Flow subroutine
!     ******************************************************************
  USE GLOBAL, ONLY:IOUT,LSTCHK
  USE GWFUNSFMODULE, ONLY: NUMNP,HTEMP,HHNEW,ITER,ITCUM,MATNUM,&
                               &HSAT,MAXIT,TOLH,TOLTH,THNEW,CAP,DT,&
                               &HBOT,DTMIN,HHOLD,KODTOP,KODBOT,KTOLD,&
                               &KBOLD,DTOPT,T,WLAYER,HCRITS,HTOP,VBOT,&
                               &CON,TOLD,HCRITA,SINK,PAR,X
  IMPLICIT NONE
!      Variables
  INTEGER,INTENT(IN)::PID
  DOUBLE PRECISION,PARAMETER::RMAX=1.E+10
  DOUBLE PRECISION,PARAMETER::RMIN=1.D-100
  LOGICAL :: CONVGF,ITCRIT,RESTART
  INTEGER :: I,M
  DOUBLE PRECISION :: EPSH,EPSTH,TH,DX
  DOUBLE PRECISION :: PB,RB,SB,PT,RT,ST
  DOUBLE PRECISION,DIMENSION(NUMNP(PID)) :: P,R,S
!      Body
!      CONTINUE ! to check -  Adam Szymkiewicz
  ITER(PID)=0
  CONVGF=.TRUE.
  RESTART=.TRUE.
  DO WHILE(RESTART) ! possible repetition of the same step
    RESTART=.FALSE.
!       Generate terms of matrix equation and solve by Gauss elimination
    CALL GWF2UNSF1SETMAT(PID)
    CALL GWF2UNSF1RESET(PID,P,R,S,PB,RB,SB,PT,RT,ST)
    CALL GWF2UNSF1SHIFT(PID)
    DO I=1,NUMNP(PID)
      HTEMP(I,PID)=HHNEW(I,PID)
!      CHECK: Commented lines were present here in the original code
	END DO
    CALL GWF2UNSF1GAUSS(PID,P,R,S,PB,RB,SB,PT,RT,ST,RMIN)
    DO I=1,NUMNP(PID)
      IF(ABS(HHNEW(I,PID)).GT.RMAX)THEN
        HHNEW(I,PID)=SIGN(RMAX,HHNEW(I,PID))
      END IF
      IF(ABS(KODTOP(PID)).EQ.4 .AND. HHNEW(I,PID).LT.HCRITA(PID)) THEN
		IF(I.EQ.NUMNP(PID)) THEN
          HHNEW(I,PID)=HCRITA(PID)
        ELSE IF(I.GT.NUMNP(PID)*9/10 .AND. SINK(I,PID).LE.0.) THEN
          HHNEW(I,PID)=HCRITA(PID)
        END IF
	  END IF
	END DO
    ITER(PID)=ITER(PID)+1
    ITCUM(PID)=ITCUM(PID)+1
!       Test for convergence
    ITCRIT=.TRUE.
    DO I=1,NUMNP(PID)
      M = MATNUM(I,PID)
      EPSTH = 0.
      EPSH  = 0.
      IF(HTEMP(I,PID).LT.HSAT(M) .AND. &
          &HHNEW(I,PID).LT.HSAT(M)) THEN
        TH=THNEW(I,PID)+CAP(I,PID)*(HHNEW(I,PID)-HTEMP(I,PID))/&
                                                    &(PAR(2,M)-PAR(1,M))
        EPSTH=ABS(THNEW(I,PID)-TH)
	  ELSE
		EPSH = ABS(HHNEW(I,PID)-HTEMP(I,PID)) ! modified by Adam Szymkiewicz
      END IF
	  IF(EPSTH.GT.TOLTH.OR.EPSH.GT.TOLH.OR. &
          &ABS(HHNEW(I,PID)).GT.RMAX*0.999) THEN
        ITCRIT=.FALSE.
        IF(ABS(HHNEW(I,PID)).GT.RMAX*0.999) ITER(PID)=MAXIT
        EXIT
	  END IF
	END DO
    IF(.NOT.ITCRIT.OR.ITER(PID).LE.1) THEN ! no convergence or the first iteration
	  IF(ITER(PID).LT.MAXIT) THEN 
	  ! max. number of iterations not exceeded, do the next iteration
        RESTART=.TRUE.
	  ELSE IF(DT(PID).LE.DTMIN(PID)) THEN
	  ! no convergence for minimum time step, stop
        CONVGF=.FALSE.
        IF(LSTCHK(3)) THEN
          WRITE(IOUT,*) 'with HBOT = ', HBOT(PID)
          WRITE(IOUT,*) ' THE NUMERICAL SOLUTION HAS NOT CONVERGED ! '
		END IF
		!PRINT *,'Numerical solution has not converged for minimum value of Dt'
        RETURN
	  ELSE
	  ! max. number of iterations exceeded, re-try the time step with reduced dt
        DO I=1,NUMNP(PID)
          HHNEW(I,PID)=HHOLD(I,PID)
          HTEMP(I,PID)=HHOLD(I,PID)
        END DO
        KODTOP(PID)=KTOLD(PID)
        KODBOT(PID)=KBOLD(PID)
        DT(PID)=DMAX1(DT(PID)/3,DTMIN(PID))
        DTOPT(PID)=DT(PID)
        T=TOLD+DT
        RESTART=.TRUE.
        ITER(PID)=0
        CONVGF=.TRUE.
		!PRINT *, 'Retrying time step with reduced Dt'
      END IF
	END IF
  END DO
!      ------
  IF(ITCRIT) THEN
    DO I=1,NUMNP(PID)
      M=MATNUM(I,PID)
      THNEW(I,PID)=THNEW(I,PID)+CAP(I,PID)*(HHNEW(I,PID)-&
                                                          &HTEMP(I,PID))
    END DO
  END IF
  CALL GWF2UNSF1SETMAT(PID) ! modified by Adam Szymkiewicz
  IF(WLAYER(PID).GE.0) THEN
    IF (HHNEW(NUMNP(PID),PID).GT. HCRITS(PID)) THEN
      KODTOP(PID)=4
      HTOP(PID)=HCRITS(PID)
	END IF
  END IF
!  DO I=1,NUMNP(PID)-2
!    IF (HHNEW(I,PID).GE.0.0D0 .AND. HHNEW(I+1,PID).LT.0.0D0) THEN
!      DX=X(I+2,PID)-X(I+1,PID)
!      VBOT(PID)=-(CON(I+2,PID)+CON(I+1,PID))/2.*&
!            &((HHNEW(I+2,PID)-HHNEW(I+1,PID))/DX+1.)
!	  EXIT
!	END IF
!  END DO
  DX=X(2,PID)-X(1,PID)
       VBOT(PID)=-(CON(1,PID)+CON(2,PID))/2.*&
            &((HHNEW(2,PID)-HHNEW(1,PID))/DX+1.)
!  DO I=1,NUMNP(PID)
!    PRINT *, I, 'H= ', HHNEW(I,PID), 'TH= ', THNEW(I,PID)
!  END DO
!  PRINT *, 'VBOT=',VBOT(PID)
!  READ(*,*)
       RETURN
END SUBROUTINE
	
SUBROUTINE GWF2UNSF1SETMAT(PID)
!     ******************************************************************
!     SET MATerial subroutine
!      WARNING: THETA should be linked before doing this
!       (use THETA=>THOLD or THETA=>THNEW)
!      Sets CON,CAP and THETA for a specific profile
!     CHECK: Depending on where this is called from, THETA is either
!            THOLD or THNEW
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: NTAB,IMODEL
       USE GWFUNSFMODULE, ONLY: NUMNP,HTAB
       USE GWFUNSFMODULE, ONLY: CONTAB,CAPTAB,HHNEW,MATNUM,PAR,CON,CAP
       USE GWFUNSFMODULE, ONLY: HSAT,HTEMP,THETAB,THETA
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER I,M,IT
       DOUBLE PRECISION :: ALH1,DLH,HI1,HI2,HIM,CONI,DH,CAPI,THEI,FK,FC1,FQ
!      Body
       ALH1 = DLOG10(-HTAB(1))
       DLH  =(DLOG10(-HTAB(NTAB))-ALH1)/(NTAB-1)
       DO I=1,NUMNP(PID)
        M = MATNUM(I,PID)
        HI1 = DMIN1(HSAT(M),HTEMP(I,PID))
        HI2 = DMIN1(HSAT(M),HHNEW(I,PID))
        HIM = 0.1*HI1+0.9*HI2
		!PRINT *, HIM, M
        IF(HI1.GE.HSAT(M) .AND. HI2.GE.HSAT(M)) THEN
         CONI=PAR(5,M)
        ELSE IF(HIM.GT.HTAB(NTAB) .AND. HIM.LE.HTAB(1)) THEN
         IT=INT((DLOG10(-HIM)-ALH1)/DLH)+1
         DH=(HIM-HTAB(IT))/(HTAB(IT+1)-HTAB(IT))
         CONI=CONTAB(IT,M)+(CONTAB(IT+1,M)-CONTAB(IT,M))*DH
        ELSE
         CONI=FK(IMODEL,HIM,PAR(:,M))
        END IF
        IF(HIM.GE.HSAT(M)) THEN
         CAPI=0.
         THEI=PAR(2,M)
        ELSE IF(HIM.GE.HTAB(NTAB).AND.HIM.LE.HTAB(1)) THEN
         IT=INT((DLOG10(-HIM)-ALH1)/DLH)+1
         DH=(HIM-HTAB(IT))/(HTAB(IT+1)-HTAB(IT))
         CAPI=CAPTAB(IT,M)+(CAPTAB(IT+1,M)-CAPTAB(IT,M))*DH
         THEI=THETAB(IT,M)+(THETAB(IT+1,M)-THETAB(IT,M))*DH
        ELSE
         CAPI=FC1(IMODEL,HIM,PAR(:,M))
         THEI=FQ(IMODEL,HIM,PAR(:,M))
        END IF
        CON(I,PID)=CONI
        CAP(I,PID)=CAPI
        THETA(I,PID)=PAR(1,M)+(THEI-PAR(1,M))
		!PRINT *, CON(I,PID), THETA(I,PID)
	   END DO
	   !READ(*,*)
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1RESET(PID,P,R,S,PB,RB,SB,PT,RT,ST)
!     ******************************************************************
!     RESET P,R,S,PB,RB,SB,PT,RT,ST
!     CHECK: Variables PB,RB,SB,PT,RT and ST are actually the first and
!            last values of the P,R and S arrays. As such, consider
!            removing them.
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: NUMNP,RTOP,RBOT,DT,X,HHOLD,CON,CAP,&
                               &WLAYER,HHNEW,SINK,THNEW,THOLD
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       DOUBLE PRECISION,INTENT(OUT) :: PB,RB,SB,PT,RT,ST
       DOUBLE PRECISION,DIMENSION(NUMNP(PID)),INTENT(OUT) :: P,R,S
       INTEGER :: N,I
       DOUBLE PRECISION :: COSALF,DXB,DX,CONB,DXA,CONA
       DOUBLE PRECISION :: B,F2,A2,A3
!      Body
       N=NUMNP(PID)
!      Finite differences
       COSALF=1.
!      Bottom BC
       DXB=X(2,PID)-X(1,PID)
       DX =DXB/2.
       CONB=(CON(1,PID)+CON(2,PID))/2.
       B  =CONB*COSALF
       S(1)=-CONB/DXB
       F2=CAP(1,PID)*DX/DT(PID)
       RB=CONB/DXB+F2
       SB=-CONB/DXB
       PB=B-SINK(1,PID)*DX+F2*HHNEW(1,PID)-(THNEW(1,PID)-THOLD(1,PID))*&
                                                   &DX/DT(PID)+RBOT(PID)
       DO I=2,N-1
        DXA=X(I  ,PID)-X(I-1,PID)
        DXB=X(I+1,PID)-X(I  ,PID)
        DX =(DXA+DXB)/2.
        CONA = (CON(I,PID)+CON(I-1,PID))/2.
        CONB = (CON(I,PID)+CON(I+1,PID))/2.
        B=(CONA-CONB)*COSALF
        A2=CONA/DXA+CONB/DXB
        A3=-CONB/DXB
        F2=CAP(I,PID)*DX/DT(PID)
        R(I)=A2+F2
        P(I)=F2*HHNEW(I,PID)-(THNEW(I,PID)-THOLD(I,PID))*DX/DT(PID)-B-&
                                                        &SINK(I,PID)*DX
        S(I)=A3
       END DO
!      Top BC
       DXA=X(N,PID)-X(N-1,PID)
       DX=DXA/2.
       CONA=(CON(N,PID)+CON(N-1,PID))/2.
       B=CONA*COSALF
       F2=CAP(N,PID)*DX/DT(PID)
       RT=CONA/DXA+F2
       ST=-CONA/DXA
       PT=F2*HHNEW(N,PID)-(THNEW(N,PID)-THOLD(N,PID))*DX/DT(PID)-&
                                                      &SINK(N,PID)*DX-B
       PT=PT-RTOP(PID)
      IF(WLAYER(PID).GE.0) THEN
        IF(HHNEW(N,PID).GT.0.) THEN
          RT=RT+1./DT(PID)
          PT=PT+DMAX1(HHOLD(N,PID),0.)/DT(PID)
        ELSE
          PT=PT+DMAX1(HHOLD(N,PID),0.)/DT(PID)
        END IF
      END IF
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1SHIFT(PID)
!     ******************************************************************
!     Apply Atmospheric BC
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:NUMNP,KODTOP,KODBOT,RTOP,HTOP,HCRITA,&
                              &WLAYER,CON,HHNEW,X,INFTOP,THNEW,THOLD,&
                              &SINK,DT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: M,N
       DOUBLE PRECISION    :: DX,VTOP
!      Body
       N=NUMNP(PID)
       IF(INFTOP(PID).GE.0 .AND. ABS(KODTOP(PID)).EQ.4 .OR. (&
         & ABS(KODBOT(PID)).EQ.1 .AND. RTOP(PID).GT.0.)) THEN
        IF(KODTOP(PID).GT.0)THEN
         M=N-1
         DX=X(N,PID)-X(M,PID)
         VTOP=-(CON(N,PID)+CON(M,PID))/2.&
              &*((HHNEW(N,PID)-HHNEW(M,PID))/DX+1.)&
              &-(THNEW(N,PID)-THOLD(N,PID))*DX/2./DT(PID)&
              &-SINK(N,PID)*DX/2.
         IF(ABS(VTOP).GT.ABS(RTOP(PID)).OR.VTOP*RTOP(PID).LE.0) THEN
          IF(ABS(KODTOP(PID)).EQ.4) KODTOP(PID)=-4
          IF(ABS(KODTOP(PID)).EQ.1) KODTOP(PID)=-1
         END IF
        ELSE
         IF(WLAYER(PID).LT.0) THEN
          IF(HHNEW(N,PID).GT.0) THEN
           IF(ABS(KODTOP(PID)).EQ.4) KODTOP=4
           IF(ABS(KODTOP(PID)).EQ.1) KODTOP=1
           HTOP(PID)=0.
          END IF
         END IF
         IF(HHNEW(N,PID).LE.HCRITA(PID))THEN
          IF(ABS(KODTOP(PID)).EQ.4) KODTOP=4
          IF(ABS(KODTOP(PID)).EQ.1) KODTOP=1
          HTOP(PID)=HCRITA(PID)
         END IF
        END IF
       END IF
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1GAUSS(PID,P,R,S,PB,RB,SB,PT,RT,ST,RMIN)
!     ******************************************************************
!     GAUSS elimination
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: NUMNP,KODTOP,KODBOT,HHNEW,HTOP,HBOT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       DOUBLE PRECISION,INTENT(IN) :: RMIN
       DOUBLE PRECISION,INTENT(INOUT) :: PB,RB,SB,PT,RT,ST
       DOUBLE PRECISION,DIMENSION(NUMNP(PID)),INTENT(INOUT) :: P,R,S
       INTEGER :: I,N
!      Body
       N=NUMNP(PID)
!      Forward
       IF(KODBOT(PID).GE.0) THEN
        P(2)=P(2)-S(1)*HBOT(PID)
       ELSE
        IF(DABS(RB).LT.RMIN) RB=RMIN
        P(2)=P(2)-PB*S(1)/RB
        R(2)=R(2)-SB*S(1)/RB
       END IF
       DO I=3,N-1
        IF(DABS(R(I-1)).LT.RMIN) R(I-1)=RMIN
        P(I)=P(I)-P(I-1)*S(I-1)/R(I-1)
        R(I)=R(I)-S(I-1)*S(I-1)/R(I-1)
       END DO
       IF(KODTOP(PID).GT.0) THEN
        P(N-1)=P(N-1)-S(N-1)*HTOP(PID)
       ELSE
        IF(DABS(R(N-1)).LT.RMIN) R(N-1)=RMIN
        P(N)=PT-P(N-1)*ST/R(N-1)
        R(N)=RT-S(N-1)*ST/R(N-1)
       END IF
!      Back
       IF(DABS(R(N-1)).LT.RMIN) R(N-1)=RMIN
       IF(KODTOP(PID).GT.0) THEN
        HHNEW(N,PID)=HTOP(PID)
        HHNEW(N-1,PID)=P(N-1)/R(N-1)
       ELSE
        HHNEW(N,PID)=P(N)/R(N)
        HHNEW(N-1,PID)=(P(N-1)-S(N-1)*HHNEW(N,PID))/R(N-1)
       END IF
       DO I=N-2,2,-1
        IF(DABS(R(I)).LT.RMIN) R(I)=RMIN
        HHNEW(I,PID)=(P(I)-S(I)*HHNEW(I+1,PID))/R(I)
       END DO
       IF(KODBOT(PID).GE.0) THEN
        HHNEW(1,PID)=HBOT(PID)
       ELSE
        IF(DABS(RB).LT.RMIN) RB=RMIN
        HHNEW(1,PID)=(PB-SB*HHNEW(2,PID))/RB
       END IF
      RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1TOTF(PID)
!     ******************************************************************
!     Calculates time-averaged VBOT
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY: VBOT,DT,SUMVBOT,TOLD,T,TOTIMOLD,HHNEW
       IMPLICIT NONE
       INTEGER,INTENT(IN)::PID
	   ! modified by Adam Szymkiewicz
       IF(TOLD(PID).LT.TOTIMOLD) THEN
        SUMVBOT(PID)=0.0
		!PRINT *, 'SUMVBOT SET TO 0.0'
		!READ(*,*)
       ELSE
        SUMVBOT(PID)=SUMVBOT(PID)+VBOT(PID)*DT(PID)
		!PRINT *, PID, SUMVBOT(PID), HHNEW(1,PID), HHNEW(2,PID)
		!PRINT *, T(PID), DT(PID), TOLD(PID)
		!READ(*,*)
       END IF
       RETURN
	END SUBROUTINE

SUBROUTINE GWF2UNSF1TLSUM
!     ******************************************************************
!     CHECK:
!     ******************************************************************
  USE GWFUNSFMODULE, ONLY: T,HB,ZONEFLUX,ITWTU,ITFLU,NPUNSF
  IMPLICIT NONE
!      Variables
  INTEGER :: K
!      Body

  WRITE(ITWTU,'(30ES12.3)') T(1),(HB(K),K=1,NPUNSF)
  WRITE(ITFLU,'(30ES12.3)') T(1),(ZONEFLUX(K),K=1,NPUNSF)

  RETURN
END SUBROUTINE
	
	SUBROUTINE GWF2UNSF1TLINF(PID)
!     ******************************************************************
!     CHECK:
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:NUMNP,CON,X,T,DT,TLEVEL,&
                              &RTOP,RROOT,VROOT,HHNEW,CUMQ,&
                              &LWAT,WCUMT,WCUMA,THNEW,&
                              &THOLD,SINK,RSOIL,PREC,&
                              &LPRINT,LPTIME,ITLOU
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I,J,N,M
       DOUBLE PRECISION :: DX,DXN,DX1,VTOP,VBOT,RINFIL,REVAP,VOLUME
!      Body
       N=NUMNP(PID)
       M=N-1
       DXN=X(N,PID)-X(M,PID)
       VTOP=-(CON(N,PID)+CON(M,PID))/2.*((HHNEW(N,PID)-HHNEW(M,PID))&
                                                            &/DXN+1.)-&
          &(THNEW(N,PID)-THOLD(N,PID))*DXN/2./DT(PID)-SINK(N,PID)*DXN/2.
       DX1=X(2,PID)-X(1,PID)
       VBOT=-(CON(2,PID)+CON(1,PID))/2.*((HHNEW(2,PID)-HHNEW(1,PID))&
                                                            &/DXN+1.)-&
          &(THNEW(1,PID)-THOLD(1,PID))*DXN/2./DT(PID)-SINK(1,PID)*DXN/2.
       RINFIL=0.
       REVAP=0.
       IF(VTOP.LT.0. .AND. PREC(PID).GT.0) RINFIL=-VTOP+RSOIL(PID)
       IF(VTOP.GE.0. .AND. PREC(PID).GT.0) RINFIL=PREC(PID)
       IF(VTOP.GT.0.)                      REVAP=VTOP+PREC(PID)
       IF(VTOP.LE.0. .AND. RSOIL(PID).GT.0 .AND. PREC(PID).GT.0) THEN
        REVAP=RSOIL(PID)
       END IF
       CUMQ(1,PID) =CUMQ(1,PID)+RTOP(PID) *DT(PID)
       CUMQ(2,PID) =CUMQ(2,PID)+RROOT(PID)*DT(PID)
       CUMQ(3,PID) =CUMQ(3,PID)+VTOP      *DT(PID)
       CUMQ(4,PID) =CUMQ(4,PID)+VROOT(PID)*DT(PID)
       CUMQ(5,PID) =CUMQ(5,PID)+VBOT      *DT(PID)
       CUMQ(6,PID) =0.
       CUMQ(7,PID) =CUMQ(7,PID)+RINFIL     *DT(PID)
       CUMQ(8,PID) =CUMQ(8,PID)+REVAP      *DT(PID)
       CUMQ(9,PID) =CUMQ(9,PID) +PREC(PID) *DT(PID)
       CUMQ(10,PID)=CUMQ(10,PID)+RSOIL(PID)*DT(PID)
       CUMQ(11,PID)=0.
       WCUMT(PID)=WCUMT(PID)+(VBOT-VTOP-VROOT(PID))*DT(PID)
       WCUMA(PID)=WCUMA(PID)+(ABS(VBOT)+ABS(VTOP)+ABS(VROOT(PID)))*&
                                                            &DT(PID)
       VOLUME=0.
       DO I=N-1,1,-1 ! modfied by Adam Szymkiewicz 31.05.2016
        J=I+1
        DX=X(J,PID)-X(I,PID)
        VOLUME=VOLUME+DX*(THNEW(I,PID)+THNEW(J,PID))/2.
       END DO
       IF(ABS(NINT(T(PID))/T(PID)-1).LT.0.00001) THEN
        IF(LPRINT.AND.LPTIME) THEN
         IF(LWAT(PID).GE.0 .OR. TLEVEL(PID).EQ.1) THEN
          IF(T(PID).LT.9999999.) THEN
           WRITE(ITLOU,"(I3,F13.4,11E13.5,2E13.5,4E13.5,I7)") &
            &PID,T,RTOP(PID),RROOT(PID),VTOP,VROOT(PID),VBOT,&
            &(CUMQ(I,PID),I=1,5),HHNEW(1,PID),CUMQ(6,PID),VOLUME,&
            &CUMQ(7,PID),CUMQ(8,PID),TLEVEL(PID)
          ELSE
           WRITE(ITLOU,"(I3,E14.8,11E13.5,2E13.5,4E13.5,I7)") &
            &PID,T,RTOP(PID),RROOT(PID),VTOP,VROOT(PID),VBOT,&
            &(CUMQ(I,PID),I=1,5),HHNEW(1,PID),CUMQ(6,PID),VOLUME,&
            &CUMQ(7,PID),CUMQ(8,PID),TLEVEL(PID)
          END IF
         END IF
        END IF
       END IF
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1TMCONT(PID)
!     ******************************************************************
!     Determines DT
!     ******************************************************************
       USE GLOBAL, ONLY: LSTCHK,IOUT
       USE GWFUNSFMODULE, ONLY:DT,DTMAX,DTOPT,DMUL,DMUL2,DTMIN,ITER,&
                              &TPRINT,TATM,T,TMAX,ITMIN,ITMAX,LMINSTEP,&
                              &DTINIT
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       DOUBLE PRECISION :: DTMX,TFIX
!      Body
       IF(LMINSTEP(PID).GE.0) THEN
        DTMX=DMIN1(DTMAX(PID),DTINIT(PID),DTOPT(PID))
        DTOPT(PID)=DTMX
        LMINSTEP(PID)=-1
       ELSE
        DTMX=DTMAX(PID)
       END IF
       TFIX=DMIN1(TPRINT,TATM(PID),TMAX(PID))
	   !PRINT *, 'T(PID)=',T(PID), 'TFIX=',TFIX
	   !PRINT *, 'DTMAX=',DTMX,'DTOPT=',DTOPT(PID)
       IF(ITER(PID).LE.ITMIN .AND. (TFIX-T(PID)).GE.DMUL*DTOPT(PID))THEN
        DTOPT(PID)=DMIN1(DTMX,DMUL*DTOPT(PID))
       END IF
       IF(ITER(PID).GE.ITMAX) THEN
        DTOPT(PID)=DMAX1(DTMIN(PID),DMUL2*DTOPT(PID))
	   END IF
	   !PRINT *, 'DTOPT2=',DTOPT(PID)
       DT(PID)=DMIN1(DTOPT(PID),TFIX-T(PID))
	   !PRINT *, 'DT1=',DT(PID)
       DT(PID)=DMIN1((TFIX-T(PID))/DNINT((TFIX-T(PID))/DT(PID)),DTMX)
       !PRINT *, 'DT2=',DT(PID)
	   IF(DT(PID).LE.0)THEN
        IF(LSTCHK(1)) THEN
         WRITE(IOUT,"(/1X,'DT<0 ! wrong time discr MODFLOW/HYDRUS')")
        END IF
        CALL USTOP(' ')
       END IF
       IF((TFIX-T(PID)).NE.DT(PID).AND.DT(PID).GT.(TFIX-T(PID)/2.))THEN
        DT(PID)=(TFIX-T(PID))/2.
	   END IF
	   !READ(*,*)
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1UPDATE(PID)
!     ******************************************************************
!     UPDATEs head values
!     ******************************************************************
       USE GWFUNSFMODULE, ONLY:NUMNP,LWAT,DT,DTOLD,HTEMP,HHNEW,HHOLD,&
                              &THOLD,THNEW
       IMPLICIT NONE
!      Variables
       INTEGER,INTENT(IN) :: PID
       INTEGER :: I
!      Body
       DO I=1,NUMNP(PID)
        IF(LWAT(PID).GE.0)THEN
         IF(HHNEW(I,PID).LT.0. .AND. HHOLD(I,PID).LT.0.)THEN
          HTEMP(I,PID)=HHNEW(I,PID)+(HHNEW(I,PID)-HHOLD(I,PID))*&
                                                     &DT(PID)/DTOLD(PID)
         ELSE
          HTEMP(I,PID)=HHNEW(I,PID)
         END IF
         HHOLD(I,PID)=HHNEW(I,PID)
         HHNEW(I,PID)=HTEMP(I,PID)
         THOLD(I,PID)=THNEW(I,PID)
        END IF
       END DO
       RETURN
	END SUBROUTINE
	
SUBROUTINE GWF2UNSF1NODOUT1(TVAL,LHEAD)
!     ******************************************************************
!     Print node information for the current time step into a specific
!     unit.
!     CHECK: NODOUT exists in the original code, but the line where it
!            is called has been commented
!     CHECK: This subroutine's goal is clearly to print data into a file
!            However, modifies HHOLD, which should not happen, even if
!            the values stored in HHOLD do not matter at the moment
!     ******************************************************************
  USE GWFUNSFMODULE, ONLY:NUMNP,HHNEW,THOLD,X,XSURF,NUMNPD,NPUNSF,&
                              &LOUTF,IPROUH,IPROUW !,IPROU
  IMPLICIT NONE
!      Variables
  LOGICAL,INTENT(IN) :: LHEAD
  DOUBLE PRECISION, INTENT(IN) :: TVAL
  CHARACTER(LEN=*),PARAMETER:: FMT1="(//'PROFILE AT TIME:',F14.3/)"
  CHARACTER(LEN=*),PARAMETER:: FMT2="(//'PROFILE AT TIME:',E15.8/)"
  CHARACTER(LEN=*),PARAMETER:: FMT3=&
        &"('Node',100(I3, '_Depth[L]',4X, I3, '_Head[L]',4X))"
  CHARACTER(LEN=*),PARAMETER:: FMT4=&
        &"('Node',100(I3, '_Depth[L]',2X, I3, '_WC[-]',3X))"
  CHARACTER(LEN=*),PARAMETER:: FMT5=&
        &"(1X,I3,4X,100(F8.3,5X,F11.3,6X))"
  CHARACTER(LEN=*),PARAMETER:: FMT6=&
        &"(I3,3X,100(F8.3,6X,F6.4,6X))"
  INTEGER :: I,J,K,PID
  DOUBLE PRECISION,DIMENSION(NUMNPD,NPUNSF)::DEPTHS,VALUES
!      Body
!      Prepare values
  DO PID=1,NPUNSF
		   ! modified by Adam Szymkiewicz, 28.05.2016
    DEPTHS(1:NUMNP(PID),PID)=X(1:NUMNP(PID),PID)-XSURF(PID)
    IF (NUMNP(PID).LT.NUMNPD) THEN
      DEPTHS((NUMNP(PID)+1):NUMNPD,PID) = 999999
      HHNEW ((NUMNP(PID)+1):NUMNPD,PID) = 999999
      THOLD ((NUMNP(PID)+1):NUMNPD,PID) = 999999
    END IF  
  END DO
!      Print Values
  WRITE(IPROUH,FMT1) TVAL
  WRITE(IPROUH,FMT3)(I,I,I=1,NPUNSF)
  DO J=NUMNPD,1,-1
    WRITE(IPROUH,FMT5) NUMNPD-J+1,(DEPTHS(J,K),HHNEW(J,K),&
                         &K=1,NPUNSF)
  END DO

  WRITE(IPROUW,FMT1) TVAL
  WRITE(IPROUW,FMT4)(I,I,I=1,NPUNSF)
  DO J=NUMNPD,1,-1
    WRITE(IPROUW,FMT6) NUMNPD-J+1,(DEPTHS(J,K),THOLD(J,K),&
                         &K=1,NPUNSF)
  END DO

  !  IF(LOUTF)THEN
!        IF(TVAL.LT.99999999.) THEN
!         WRITE(IPROU,FMT1) TVAL
!        ELSE
!         WRITE(IPROU,FMT2) TVAL
!        END IF
!        IF(LHEAD) THEN
!         WRITE(IPROU,FMT3)(I,I,I=1,NPUNSF)
!        ELSE
!         WRITE(IPROU,FMT4)(I,I,I=1,NPUNSF)
!        END IF
!        IF (LHEAD) THEN
!         VALUES=HHNEW
!         DO J=NUMNPD,1,-1
!          WRITE(IPROU,FMT5)NUMNPD-J+1,(DEPTHS(J,K),VALUES(J,K),&
!                         &K=1,NPUNSF)
!         END DO
!        ELSE
!         VALUES=THOLD
!         DO J=NUMNPD,1,-1
!          WRITE(IPROU,FMT6)NUMNPD-J+1,(DEPTHS(J,K),VALUES(J,K),&
!                         &K=1,NPUNSF)
!         END DO
!        END IF
!       END IF
  RETURN
END SUBROUTINE

SUBROUTINE GWF2UNSF1FOPEN()
!     ******************************************************************
!     Opens Hydrus-specific output files
!     ******************************************************************
  USE GLOBAL, ONLY:IOUT,LSTCHK
!       USE GWFUNSFMODULE, ONLY: IPROU,ITLOU,LOUTF
  USE GWFUNSFMODULE, ONLY: IPROUH,IPROUW,ITWTU,ITFLU,LOUTF
  IMPLICIT NONE
  LOGICAL O1,O2
!       INQUIRE(UNIT=IPROU,OPENED=O1) ! modified by Adam Szymkiewicz
!       INQUIRE(UNIT=ITLOU,OPENED=O2)
  INQUIRE(UNIT=IPROUH,OPENED=O1)
  INQUIRE(UNIT=IPROUW,OPENED=O2)
  IF (.NOT.O1) OPEN (UNIT=IPROUH, FILE="prof_h.out")
  IF (.NOT.O2) OPEN (UNIT=IPROUW, FILE="prof_w.out")
  INQUIRE(UNIT=IPROUH,OPENED=O1)
  INQUIRE(UNIT=IPROUW,OPENED=O2)
  IF (.NOT.(O1 .AND. O2)) STOP "COULD NOT OPEN 'prof_h.out' AND 'prof_w.out' "

  INQUIRE(UNIT=ITWTU,OPENED=O1)
  INQUIRE(UNIT=ITFLU,OPENED=O2)
  IF (.NOT.O1) OPEN (UNIT=ITWTU, FILE="tlev_h.out")
  IF (.NOT.O2) OPEN (UNIT=ITFLU, FILE="tlev_q.out")
  INQUIRE(UNIT=ITWTU,OPENED=O1)
  INQUIRE(UNIT=ITFLU,OPENED=O2)
  IF (.NOT.(O1 .AND. O2)) STOP "COULD NOT OPEN 'tlev_h.out' AND 'tlev_q.out' "
  
  !       IF(O1 .AND. O2) THEN
!        LOUTF=.TRUE.
!       ELSE
!        IF(LSTCHK(2)) WRITE(IOUT,*) "COULD NOT OPEN 'prof_h.out' AND&
!                                    & 'prof_w.out' "
!        LOUTF=.FALSE.
!       END IF
  RETURN
END SUBROUTINE
 
SUBROUTINE GWF2UNSF1FCLOSE()
!     ******************************************************************
!     CLOSES Hydrus-specific output files
!     ******************************************************************
  USE GWFUNSFMODULE, ONLY: IPROUH,IPROUW,ITWTU,ITFLU ! IPROU,ITLOU
  IMPLICIT NONE
  LOGICAL O1,O2
  INQUIRE(UNIT=IPROUH,OPENED=O1)
  INQUIRE(UNIT=IPROUW,OPENED=O2)
  IF(O1)CLOSE(UNIT=IPROUH)
  IF(O2)CLOSE(UNIT=IPROUW)
  INQUIRE(UNIT=ITWTU,OPENED=O1)
  INQUIRE(UNIT=ITFLU,OPENED=O2)
  IF(O1)CLOSE(UNIT=ITWTU)
  IF(O2)CLOSE(UNIT=ITFLU)
END SUBROUTINE

	!     ******************************************************************
!     *                DEALLOCATION AND DATA TRANSFER                  *
!     ******************************************************************
      SUBROUTINE GWF2UNSF1PSV(IGRID)
!     ******************************************************************
!     Saves the current data to target Grid
!     ******************************************************************
       USE GWFUNSFMODULE
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IGRID
!      Body
       GWFUNSFDAT(IGRID)%LPRINT=>LPRINT
       GWFUNSFDAT(IGRID)%LPTIME=>LPTIME
       GWFUNSFDAT(IGRID)%LOUTF=>LOUTF
       GWFUNSFDAT(IGRID)%LHEAD=>LHEAD
       GWFUNSFDAT(IGRID)%NPUNSF=>NPUNSF
       GWFUNSFDAT(IGRID)%NUNSFOP=>NUNSFOP
       GWFUNSFDAT(IGRID)%IUNSFCB=>IUNSFCB
       GWFUNSFDAT(IGRID)%IUNSFPR=>IUNSFPR
       GWFUNSFDAT(IGRID)%NIZ=>NIZ
       GWFUNSFDAT(IGRID)%NUMNPD=>NUMNPD
       GWFUNSFDAT(IGRID)%IMODEL=>IMODEL
       GWFUNSFDAT(IGRID)%ITMIN=>ITMIN
       GWFUNSFDAT(IGRID)%ITMAX=>ITMAX
       GWFUNSFDAT(IGRID)%MAXIT=>MAXIT
       GWFUNSFDAT(IGRID)%PROPR=>PROPR
       GWFUNSFDAT(IGRID)%PROINF=>PROINF
       GWFUNSFDAT(IGRID)%MAXATM=>MAXATM
       GWFUNSFDAT(IGRID)%TOLTH=>TOLTH
       GWFUNSFDAT(IGRID)%TOLH=>TOLH
       GWFUNSFDAT(IGRID)%TPRINT=>TPRINT
       GWFUNSFDAT(IGRID)%XCONV=>XCONV
       GWFUNSFDAT(IGRID)%GWL0L=>GWL0L
       GWFUNSFDAT(IGRID)%RB=>RB
       GWFUNSFDAT(IGRID)%TOTIMOLD=>TOTIMOLD
       GWFUNSFDAT(IGRID)%DMUL=>DMUL
       GWFUNSFDAT(IGRID)%DMUL2=>DMUL2
       GWFUNSFDAT(IGRID)%NMAT=>NMAT
       GWFUNSFDAT(IGRID)%NTAB=>NTAB
       GWFUNSFDAT(IGRID)%HTAB=>HTAB
       GWFUNSFDAT(IGRID)%PAR=>PAR
       GWFUNSFDAT(IGRID)%IZ=>IZ
       GWFUNSFDAT(IGRID)%PLEVEL=>PLEVEL
       GWFUNSFDAT(IGRID)%ALEVEL=>ALEVEL
       GWFUNSFDAT(IGRID)%TLEVEL=>TLEVEL
       GWFUNSFDAT(IGRID)%NUMNP=>NUMNP
       GWFUNSFDAT(IGRID)%ITCUM=>ITCUM
       GWFUNSFDAT(IGRID)%ITER=>ITER
       GWFUNSFDAT(IGRID)%LMINSTEP=>LMINSTEP
       GWFUNSFDAT(IGRID)%KODTOP=>KODTOP
       GWFUNSFDAT(IGRID)%KODBOT=>KODBOT
       GWFUNSFDAT(IGRID)%INFTOP=>INFTOP
       GWFUNSFDAT(IGRID)%INFBOT=>INFBOT
       GWFUNSFDAT(IGRID)%KTOLD=>KTOLD
       GWFUNSFDAT(IGRID)%KBOLD=>KBOLD
       GWFUNSFDAT(IGRID)%NOLAY=>NOLAY
       GWFUNSFDAT(IGRID)%NOMAT=>NOMAT
       GWFUNSFDAT(IGRID)%SINKF=>SINKF
       GWFUNSFDAT(IGRID)%SHORTO=>SHORTO
       GWFUNSFDAT(IGRID)%ATMBC=>ATMBC
       GWFUNSFDAT(IGRID)%LWAT=>LWAT
       GWFUNSFDAT(IGRID)%WLAYER=>WLAYER
       GWFUNSFDAT(IGRID)%LINITW=>LINITW
       GWFUNSFDAT(IGRID)%MAXAL=>MAXAL
       GWFUNSFDAT(IGRID)%PERIMP=>PERIMP
       GWFUNSFDAT(IGRID)%HBOT=>HBOT
       GWFUNSFDAT(IGRID)%VBOT=>VBOT
       GWFUNSFDAT(IGRID)%HTOP=>HTOP
       GWFUNSFDAT(IGRID)%HROOT=>HROOT
       GWFUNSFDAT(IGRID)%VROOT=>VROOT
       GWFUNSFDAT(IGRID)%WCUMT=>WCUMT
       GWFUNSFDAT(IGRID)%WCUMA=>WCUMA
       GWFUNSFDAT(IGRID)%WVOLI=>WVOLI
       GWFUNSFDAT(IGRID)%SUMVBOT=>SUMVBOT
	   GWFUNSFDAT(IGRID)%ZONEFLUX=>ZONEFLUX
       GWFUNSFDAT(IGRID)%RTOP=>RTOP
       GWFUNSFDAT(IGRID)%RBOT=>RBOT
       GWFUNSFDAT(IGRID)%RROOT=>RROOT
       GWFUNSFDAT(IGRID)%HCRITS=>HCRITS
       GWFUNSFDAT(IGRID)%HCRITA=>HCRITA
       GWFUNSFDAT(IGRID)%XSURF=>XSURF
       GWFUNSFDAT(IGRID)%TINIT=>TINIT
       GWFUNSFDAT(IGRID)%TMAX=>TMAX
       GWFUNSFDAT(IGRID)%TOLD=>TOLD
       GWFUNSFDAT(IGRID)%T=>T
       GWFUNSFDAT(IGRID)%DT=>DT
       GWFUNSFDAT(IGRID)%DTMIN=>DTMIN
       GWFUNSFDAT(IGRID)%DTMAX=>DTMAX
       GWFUNSFDAT(IGRID)%DTOPT=>DTOPT
       GWFUNSFDAT(IGRID)%DTINIT=>DTINIT
       GWFUNSFDAT(IGRID)%DTOLD=>DTOLD
       GWFUNSFDAT(IGRID)%P0=>P0
       GWFUNSFDAT(IGRID)%P2H=>P2H
       GWFUNSFDAT(IGRID)%P2L=>P2L
       GWFUNSFDAT(IGRID)%P3=>P3
       GWFUNSFDAT(IGRID)%R2H=>R2H
       GWFUNSFDAT(IGRID)%R2L=>R2L
       GWFUNSFDAT(IGRID)%POPTM=>POPTM
       GWFUNSFDAT(IGRID)%TATM=>TATM
       GWFUNSFDAT(IGRID)%PREC=>PREC
       GWFUNSFDAT(IGRID)%RSOIL=>RSOIL
       GWFUNSFDAT(IGRID)%HB=>HB
       GWFUNSFDAT(IGRID)%MATNUM=>MATNUM
       GWFUNSFDAT(IGRID)%LAYNUM=>LAYNUM
       GWFUNSFDAT(IGRID)%WATIN=>WATIN
       GWFUNSFDAT(IGRID)%SINK=>SINK
       GWFUNSFDAT(IGRID)%X=>X
       GWFUNSFDAT(IGRID)%HHNEW=>HHNEW
       GWFUNSFDAT(IGRID)%HHOLD=>HHOLD
       GWFUNSFDAT(IGRID)%HTEMP=>HTEMP
       GWFUNSFDAT(IGRID)%BETA=>BETA
       GWFUNSFDAT(IGRID)%THOLD=>THOLD
       GWFUNSFDAT(IGRID)%THNEW=>THNEW
       GWFUNSFDAT(IGRID)%CON=>CON
       GWFUNSFDAT(IGRID)%CAP=>CAP
       GWFUNSFDAT(IGRID)%TTATM=>TTATM
       GWFUNSFDAT(IGRID)%PRECI=>PRECI
       GWFUNSFDAT(IGRID)%RSO=>RSO
       GWFUNSFDAT(IGRID)%HT=>HT
       GWFUNSFDAT(IGRID)%RR=>RR
       GWFUNSFDAT(IGRID)%HCA=>HCA
       GWFUNSFDAT(IGRID)%CONTAB=>CONTAB
       GWFUNSFDAT(IGRID)%CAPTAB=>CAPTAB
       GWFUNSFDAT(IGRID)%THETAB=>THETAB
       GWFUNSFDAT(IGRID)%HSAT=>HSAT
       GWFUNSFDAT(IGRID)%UNSFLUX=>UNSFLUX
       GWFUNSFDAT(IGRID)%UNSFLUXV=>UNSFLUXV
       GWFUNSFDAT(IGRID)%CUMQ=>CUMQ
       GWFUNSFDAT(IGRID)%PRTIME=>PRTIME
       GWFUNSFDAT(IGRID)%KZON=>KZON
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1PNT(IGRID)
!     ******************************************************************
!     Changes the current data to point to data stored in target grid
!     ******************************************************************
       USE GWFUNSFMODULE
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IGRID
!      Body
       LPRINT=>GWFUNSFDAT(IGRID)%LPRINT
       LPTIME=>GWFUNSFDAT(IGRID)%LPTIME
       LOUTF=>GWFUNSFDAT(IGRID)%LOUTF
       LHEAD=>GWFUNSFDAT(IGRID)%LHEAD
       NPUNSF=>GWFUNSFDAT(IGRID)%NPUNSF
       NUNSFOP=>GWFUNSFDAT(IGRID)%NUNSFOP
       IUNSFCB=>GWFUNSFDAT(IGRID)%IUNSFCB
       IUNSFPR=>GWFUNSFDAT(IGRID)%IUNSFPR
       NIZ=>GWFUNSFDAT(IGRID)%NIZ
       NUMNPD=>GWFUNSFDAT(IGRID)%NUMNPD
       IMODEL=>GWFUNSFDAT(IGRID)%IMODEL
       ITMIN=>GWFUNSFDAT(IGRID)%ITMIN
       ITMAX=>GWFUNSFDAT(IGRID)%ITMAX
       MAXIT=>GWFUNSFDAT(IGRID)%MAXIT
       PROPR=>GWFUNSFDAT(IGRID)%PROPR
       PROINF=>GWFUNSFDAT(IGRID)%PROINF
       MAXATM=>GWFUNSFDAT(IGRID)%MAXATM
       TOLTH=>GWFUNSFDAT(IGRID)%TOLTH
       TOLH=>GWFUNSFDAT(IGRID)%TOLH
       TPRINT=>GWFUNSFDAT(IGRID)%TPRINT
       XCONV=>GWFUNSFDAT(IGRID)%XCONV
       GWL0L=>GWFUNSFDAT(IGRID)%GWL0L
       RB=>GWFUNSFDAT(IGRID)%RB
       TOTIMOLD=>GWFUNSFDAT(IGRID)%TOTIMOLD
       DMUL=>GWFUNSFDAT(IGRID)%DMUL
       DMUL2=>GWFUNSFDAT(IGRID)%DMUL2
       NMAT=>GWFUNSFDAT(IGRID)%NMAT
       NTAB=>GWFUNSFDAT(IGRID)%NTAB
       HTAB=>GWFUNSFDAT(IGRID)%HTAB
       PAR=>GWFUNSFDAT(IGRID)%PAR
       IZ=>GWFUNSFDAT(IGRID)%IZ
       PLEVEL=>GWFUNSFDAT(IGRID)%PLEVEL
       ALEVEL=>GWFUNSFDAT(IGRID)%ALEVEL
       TLEVEL=>GWFUNSFDAT(IGRID)%TLEVEL
       NUMNP=>GWFUNSFDAT(IGRID)%NUMNP
       ITCUM=>GWFUNSFDAT(IGRID)%ITCUM
       ITER=>GWFUNSFDAT(IGRID)%ITER
       LMINSTEP=>GWFUNSFDAT(IGRID)%LMINSTEP
       KODTOP=>GWFUNSFDAT(IGRID)%KODTOP
       KODBOT=>GWFUNSFDAT(IGRID)%KODBOT
       INFTOP=>GWFUNSFDAT(IGRID)%INFTOP
       INFBOT=>GWFUNSFDAT(IGRID)%INFBOT
       KTOLD=>GWFUNSFDAT(IGRID)%KTOLD
       KBOLD=>GWFUNSFDAT(IGRID)%KBOLD
       NOLAY=>GWFUNSFDAT(IGRID)%NOLAY
       NOMAT=>GWFUNSFDAT(IGRID)%NOMAT
       SINKF=>GWFUNSFDAT(IGRID)%SINKF
       SHORTO=>GWFUNSFDAT(IGRID)%SHORTO
       ATMBC=>GWFUNSFDAT(IGRID)%ATMBC
       LWAT=>GWFUNSFDAT(IGRID)%LWAT
       WLAYER=>GWFUNSFDAT(IGRID)%WLAYER
       LINITW=>GWFUNSFDAT(IGRID)%LINITW
       MAXAL=>GWFUNSFDAT(IGRID)%MAXAL
       PERIMP=>GWFUNSFDAT(IGRID)%PERIMP
       HBOT=>GWFUNSFDAT(IGRID)%HBOT
       VBOT=>GWFUNSFDAT(IGRID)%VBOT
       HTOP=>GWFUNSFDAT(IGRID)%HTOP
       HROOT=>GWFUNSFDAT(IGRID)%HROOT
       VROOT=>GWFUNSFDAT(IGRID)%VROOT
       WCUMT=>GWFUNSFDAT(IGRID)%WCUMT
       WCUMA=>GWFUNSFDAT(IGRID)%WCUMA
       WVOLI=>GWFUNSFDAT(IGRID)%WVOLI
       SUMVBOT=>GWFUNSFDAT(IGRID)%SUMVBOT
	   ZONEFLUX=>GWFUNSFDAT(IGRID)%ZONEFLUX
       RTOP=>GWFUNSFDAT(IGRID)%RTOP
       RBOT=>GWFUNSFDAT(IGRID)%RBOT
       RROOT=>GWFUNSFDAT(IGRID)%RROOT
       HCRITS=>GWFUNSFDAT(IGRID)%HCRITS
       HCRITA=>GWFUNSFDAT(IGRID)%HCRITA
       XSURF=>GWFUNSFDAT(IGRID)%XSURF
       TINIT=>GWFUNSFDAT(IGRID)%TINIT
       TMAX=>GWFUNSFDAT(IGRID)%TMAX
       TOLD=>GWFUNSFDAT(IGRID)%TOLD
       T=>GWFUNSFDAT(IGRID)%T
       DT=>GWFUNSFDAT(IGRID)%DT
       DTMIN=>GWFUNSFDAT(IGRID)%DTMIN
       DTMAX=>GWFUNSFDAT(IGRID)%DTMAX
       DTOPT=>GWFUNSFDAT(IGRID)%DTOPT
       DTINIT=>GWFUNSFDAT(IGRID)%DTINIT
       DTOLD=>GWFUNSFDAT(IGRID)%DTOLD
       P0=>GWFUNSFDAT(IGRID)%P0
       P2H=>GWFUNSFDAT(IGRID)%P2H
       P2L=>GWFUNSFDAT(IGRID)%P2L
       P3=>GWFUNSFDAT(IGRID)%P3
       R2H=>GWFUNSFDAT(IGRID)%R2H
       R2L=>GWFUNSFDAT(IGRID)%R2L
       POPTM=>GWFUNSFDAT(IGRID)%POPTM
       TATM=>GWFUNSFDAT(IGRID)%TATM
       PREC=>GWFUNSFDAT(IGRID)%PREC
       RSOIL=>GWFUNSFDAT(IGRID)%RSOIL
       HB=>GWFUNSFDAT(IGRID)%HB
       MATNUM=>GWFUNSFDAT(IGRID)%MATNUM
       LAYNUM=>GWFUNSFDAT(IGRID)%LAYNUM
       WATIN=>GWFUNSFDAT(IGRID)%WATIN
       SINK=>GWFUNSFDAT(IGRID)%SINK
       X=>GWFUNSFDAT(IGRID)%X
       HHNEW=>GWFUNSFDAT(IGRID)%HHNEW
       HHOLD=>GWFUNSFDAT(IGRID)%HHOLD
       HTEMP=>GWFUNSFDAT(IGRID)%HTEMP
       BETA=>GWFUNSFDAT(IGRID)%BETA
       THOLD=>GWFUNSFDAT(IGRID)%THOLD
       THNEW=>GWFUNSFDAT(IGRID)%THNEW
       CON=>GWFUNSFDAT(IGRID)%CON
       CAP=>GWFUNSFDAT(IGRID)%CAP
       TTATM=>GWFUNSFDAT(IGRID)%TTATM
       PRECI=>GWFUNSFDAT(IGRID)%PRECI
       RSO=>GWFUNSFDAT(IGRID)%RSO
       HT=>GWFUNSFDAT(IGRID)%HT
       RR=>GWFUNSFDAT(IGRID)%RR
       HCA=>GWFUNSFDAT(IGRID)%HCA
       CONTAB=>GWFUNSFDAT(IGRID)%CONTAB
       CAPTAB=>GWFUNSFDAT(IGRID)%CAPTAB
       THETAB=>GWFUNSFDAT(IGRID)%THETAB
       HSAT=>GWFUNSFDAT(IGRID)%HSAT
       UNSFLUX=>GWFUNSFDAT(IGRID)%UNSFLUX
       UNSFLUXV=>GWFUNSFDAT(IGRID)%UNSFLUXV
       CUMQ=>GWFUNSFDAT(IGRID)%CUMQ
       PRTIME=>GWFUNSFDAT(IGRID)%PRTIME
       KZON=>GWFUNSFDAT(IGRID)%KZON
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1DA(IGRID)
!     ******************************************************************
!     Deallocates the data
!     ******************************************************************
       USE GWFUNSFMODULE
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: IGRID
!      Body
       DEALLOCATE(GWFUNSFDAT(IGRID)%LPRINT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LPTIME)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LOUTF)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LHEAD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NPUNSF)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NUNSFOP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%IUNSFCB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%IUNSFPR)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NIZ)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NUMNPD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%IMODEL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ITMIN)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ITMAX)
       DEALLOCATE(GWFUNSFDAT(IGRID)%MAXIT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PROPR)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PROINF)
       DEALLOCATE(GWFUNSFDAT(IGRID)%MAXATM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TOLTH)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TOLH)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TPRINT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%XCONV)
       DEALLOCATE(GWFUNSFDAT(IGRID)%GWL0L)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TOTIMOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DMUL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DMUL2)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NMAT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NTAB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HTAB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PAR)
       DEALLOCATE(GWFUNSFDAT(IGRID)%IZ)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PLEVEL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ALEVEL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TLEVEL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NUMNP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ITCUM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ITER)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LMINSTEP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%KODTOP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%KODBOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%INFTOP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%INFBOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%KTOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%KBOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NOLAY)
       DEALLOCATE(GWFUNSFDAT(IGRID)%NOMAT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%SINKF)
       DEALLOCATE(GWFUNSFDAT(IGRID)%SHORTO)
       DEALLOCATE(GWFUNSFDAT(IGRID)%ATMBC)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LWAT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%WLAYER)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LINITW)
       DEALLOCATE(GWFUNSFDAT(IGRID)%MAXAL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PERIMP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HBOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%VBOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HTOP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HROOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%VROOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%WCUMT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%WCUMA)
       DEALLOCATE(GWFUNSFDAT(IGRID)%WVOLI)
       DEALLOCATE(GWFUNSFDAT(IGRID)%SUMVBOT)
	   DEALLOCATE(GWFUNSFDAT(IGRID)%ZONEFLUX)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RTOP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RBOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RROOT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HCRITS)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HCRITA)
       DEALLOCATE(GWFUNSFDAT(IGRID)%XSURF)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TINIT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TMAX)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%T)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DTMIN)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DTMAX)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DTOPT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DTINIT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%DTOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%P0)
       DEALLOCATE(GWFUNSFDAT(IGRID)%P2H)
       DEALLOCATE(GWFUNSFDAT(IGRID)%P2L)
       DEALLOCATE(GWFUNSFDAT(IGRID)%P3)
       DEALLOCATE(GWFUNSFDAT(IGRID)%R2H)
       DEALLOCATE(GWFUNSFDAT(IGRID)%R2L)
       DEALLOCATE(GWFUNSFDAT(IGRID)%POPTM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TATM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PREC)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RSOIL)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%MATNUM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%LAYNUM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%WATIN)
       DEALLOCATE(GWFUNSFDAT(IGRID)%SINK)
       DEALLOCATE(GWFUNSFDAT(IGRID)%X)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HHNEW)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HHOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HTEMP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%BETA)
       DEALLOCATE(GWFUNSFDAT(IGRID)%THOLD)
       DEALLOCATE(GWFUNSFDAT(IGRID)%THNEW)
       DEALLOCATE(GWFUNSFDAT(IGRID)%CON)
       DEALLOCATE(GWFUNSFDAT(IGRID)%CAP)
       DEALLOCATE(GWFUNSFDAT(IGRID)%TTATM)
       DEALLOCATE(GWFUNSFDAT(IGRID)%PRECI)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RSO)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%RR)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HCA)
       DEALLOCATE(GWFUNSFDAT(IGRID)%CONTAB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%CAPTAB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%THETAB)
       DEALLOCATE(GWFUNSFDAT(IGRID)%HSAT)
       DEALLOCATE(GWFUNSFDAT(IGRID)%UNSFLUX)
       DEALLOCATE(GWFUNSFDAT(IGRID)%UNSFLUXV)
       DEALLOCATE(GWFUNSFDAT(IGRID)%CUMQ)
       IF(ASSOCIATED(PRTIME))DEALLOCATE(GWFUNSFDAT(IGRID)%PRTIME)
       DEALLOCATE(GWFUNSFDAT(IGRID)%KZON)
       RETURN
      END SUBROUTINE
      SUBROUTINE GWF2UNSF1PREAL()
!     ******************************************************************
!     Allocates elements of known size
!     ******************************************************************
       USE GWFUNSFMODULE
       ALLOCATE(LPRINT)
       ALLOCATE(LPTIME)
       ALLOCATE(LOUTF)
       ALLOCATE(LHEAD)
       ALLOCATE(NPUNSF)
       ALLOCATE(NUNSFOP)
       ALLOCATE(IUNSFCB)
       ALLOCATE(IUNSFPR)
       ALLOCATE(NIZ)
       ALLOCATE(NUMNPD)
       ALLOCATE(IMODEL)
       ALLOCATE(ITMIN)
       ALLOCATE(ITMAX)
       ALLOCATE(MAXIT)
       ALLOCATE(PROPR)
       ALLOCATE(PROINF)
       ALLOCATE(MAXATM)
       ALLOCATE(TOLTH)
       ALLOCATE(TOLH)
       ALLOCATE(TPRINT)
       ALLOCATE(XCONV)
       ALLOCATE(GWL0L)
       ALLOCATE(RB)
       ALLOCATE(TOTIMOLD)
       ALLOCATE(DMUL)
       ALLOCATE(DMUL2)
       ALLOCATE(NMAT)
       ALLOCATE(NTAB)
       RETURN
      END SUBROUTINE
