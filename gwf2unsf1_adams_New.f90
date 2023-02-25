    !     ******************************************************************
    !     Unsaturated Flow (UNSF) Package
    !     Author : Gleb Goussarov
    !     Code based on gwf1unsf.dtf by SEO et al., in line with the general
    !     style of the MF OWH project. Note however, that some f77 features
    !     have been replaced by f90 features.
    !     Subroutines in this package use "IMPLICIT NONE"
    !     **************************************************************
    MODULE GWFUNSFMODULE
    !     ******************************************************************
    !     Module Definition
    !     ******************************************************************
    !     INTEGER,PARAMETER:: IPROU = 97 ! unit for writing profiles, AdamS: changed to 2 separate files
    INTEGER,PARAMETER:: IPROUH = 191    !AdamS: unit for writing pressure profiles
    INTEGER,PARAMETER:: IPROUW = 192    !AdamS: unit for water content profiles
    INTEGER,PARAMETER:: ITWTU = 193     !AdamS: unit for writing water table at each Modflow step
    INTEGER,PARAMETER:: ITFLU = 194     !AdamS: unit for writing bottom flux at each Modflow step
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
    integer:: iBact,switch
    !      NMAT   :
    !      NTAB   :
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS :: HTAB
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: PAR
    !      ------------------------ SIZE = NPUNSF --------------------------
    INTEGER,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::&
        &IZ,PLEVEL,ALEVEL,TLEVEL,NUMNP,&
        &ITCUM,ITER,LMINSTEP,&
        &KODTOP,KODBOT,INFTOP,INFBOT,KTOLD,KBOLD,NOLAY,NOMAT,&
        &SHORTO,ATMBC,&
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
        &P0,P2H,P2L,P3,R2H,R2L,&
        &TATM,PREC,RSOIL,&
        &HB,ZONEFLUX,time
    !      ZONEFLUX : added by AdamS, contains averaged bottom fluxes for each profile, used for Modflow coupling
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
        &MATNUM,LAYNUM
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: POPTM
    !      MATNUM:
    !      LAYNUM:
    !      WATIN :
    DOUBLE PRECISION,   SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
        &SINK,X,HHNEW,HHOLD,HTEMP,BETA,THOLD,THNEW,CON,CAP,thorg,thmod,hh
    
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
    !       ------------------------------------------------------------------
    !      Defining all variable and parameters for Solute transport: Modified by Sahila
    LOGICAL :: LCHEM,lUpW,lTDep ,lArtD,lScreen ,lMassIni, &
        &lTort,lBact,lFiltr,iVer,lDualNEq, lEqInit, lvar,&
        &lmoistdep,iequil,lequil!,Ibact
    LOGICAL ::  lActRSU , Lconv,lDensity,  &
        &LNEQUIL,  lVapor
    LOGICAL,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::lLinear,lMobIm
    integer,save,dimension (:),pointer,CONTIGUOUS ::sinkf
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::cprev
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::s
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
        &cNew,cvBot,cvCh0,cvCh1,cvChIm,cvChR,cvTop,&
        &g0,g1,q0,q1,Retard, SinkIm,sSink,STrans,TempN,TempO,B,&
        &ThNIm,ThOIm,wc  ,Disp ,thnewim,tholdim, cPrevO,cTemp,D,E,F,&
        &p,r,q,cRootMax
    INTEGER,SAVE,POINTER :: NS
    INTEGER::  npar2,ierr,INONEQUL,num
    INTEGER,PARAMETER:: NSD=6
    INTEGER,PARAMETER:: NMATD=20
    INTEGER,PARAMETER:: Nlevel=2
    REAL::epsi,tPulse, cTolA,cTolR,PeCr,cAtm,cmin,iconctype,maxitc
    real::OmegaS,OmegaW, rkm,ikod
    DOUBLE PRECISION,SAVE,DIMENSION(:,:,:,:),POINTER,CONTIGUOUS :: DMOIST
    DOUBLE PRECISION,SAVE,DIMENSION(:,:,:),POINTER,CONTIGUOUS ::concnew1,conc,sorb,sorb2,&
        &concmod,concnew,CT,CB,concold
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS ::CHpar,Wdep,CTOP,CBOT
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: CumCh,SorbN,SorbN2,TEMP,concf
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: WATIN
    DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER,CONTIGUOUS :: VOld,VNew,vmid,vc
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::tdep,cctop,ccbot
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS :: C
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS ::ths,thsat,N
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS :: xmasschange
    DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER,CONTIGUOUS:: MATNUM_S
    DOUBLE PRECISION,save,DIMENSION(:),POINTER,CONTIGUOUS:: VO_S,VN_S,x_s
    real::DTOLD_S,T_S,TLEVEL_S,Courant,peclet
    double precision dtMaxC,thimobo, smax2,tconv, thimob,smax,thw
    real::theta0,theta1,theta2,theta3,ReacMin1,DT_S ,DTMIN_S,rbot_s,courmax
    real::idualpor ,itort,spot
    integer ktopch,kbotch
    integer:: atmbc_s, LWAT_S,sinkf_s ,level,js,jreact,M,iter_s,imoistdep,iterc
    DOUBLE PRECISION,SAVE,POINTER::b1
    DOUBLE PRECISION  thwo,thg,smax1i,smaxij,Smax2i,smax2j,thetai,psi1i,&
        & xnu,xku,v,taug,ro,henry,frac,fexp,dw,dsconc,dg,cmid,smid,xks,gams1oi,&
        &xmuso,xmus,xmulo,xmul,xkso,ssorb2,ssorb,sconcs,sconcos,sconco,sconc,&
        &rkd1,rka2o,rka2,rka1o,rka1,omegao,omega,gamso,gams1o,gams1,gams,&
        &gaml1,gaml,dsconcs,cc,xksp,sconcps,sconcp,rkd2o,rkd2,rkd1o,psi2o,psi2,&
        & psi1o,psi1,gams1pi,gaml1pi,flmacro,f_em,dmobi,fn,f1,e1,dn,d1,bn,alf&
        &,vj,thj,ss1,ss2,gamli,gaml1i,gamloi,gaml1oi,gamsi,gams1i,gamsoi, henryj&
        &,gamg1i,smax1,ipsi2,ipsi1,dc,aa,dp,alfa1,alfa2,xnup,henryp,gamlo,gaml1o,&
        &gaml1p,gams1p,gamg1p,tto,xnuo,fexpo,henryo,dsurf,dsurft,psi2i,psi2j,&
        &DERK,DX,CG1,CG,DCONCS,DRETARDS,DRETARD,CONCP,DCONC,smax2o,smax1o,ss,&
        &psi1j,smax1j
    DOUBLE PRECISION,SAVE, DIMENSION(:),POINTER,CONTIGUOUS:: THETA_S,VELOC
    DOUBLE PRECISION,save, DIMENSION(:),POINTER,CONTIGUOUS:: wvol
    DOUBLE PRECISION,save, DIMENSION(:),POINTER,CONTIGUOUS:: area,subvol
    DOUBLE PRECISION,save, DIMENSION(:),POINTER,CONTIGUOUS:: ccuma,cvoli,solin
    DOUBLE PRECISION,save, DIMENSION(:),POINTER,CONTIGUOUS:: ccumt ,cRunOff,cgwl
    DOUBLE PRECISION,save, DIMENSION(:),POINTER,CONTIGUOUS:: solfluxwat
    DOUBLE PRECISION VOLUME ,initime
    !      ------------------------------------------------------------------


    TYPE GWFUNSFTYPE
        LOGICAL,POINTER :: LPRINT,LPTIME,LOUTF,LHEAD
        INTEGER,POINTER :: &
            &NPUNSF,NUNSFOP,IUNSFCB,IUNSFPR,NIZ,NUMNPD,&
            &IMODEL,ITMIN,ITMAX,MAXIT,PROPR,PROINF,&
            &MAXATM,&
            &NMAT,NTAB,ns
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
            &SHORTO,ATMBC,&
            &LWAT,WLAYER,LINITW,&
            &MAXAL,PERIMP,sinkf
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS ::&
            &HBOT,VBOT,HTOP,HROOT,VROOT,WCUMT,WCUMA,WVOLI,SUMVBOT,&
            &RTOP,RBOT,RROOT,HCRITS,HCRITA,XSURF,&
            &TINIT,TMAX,TOLD,T,DT,DTMIN,DTMAX,DTOPT,DTINIT,DTOLD,&
            &P0,P2H,P2L,P3,R2H,R2L,&
            &TATM,PREC,RSOIL,&
            &HB,ZONEFLUX,time ! modifed by Adam Szymkiewicz
        INTEGER,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
            &MATNUM,LAYNUM
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS ::POPTM
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
            &SINK,X,HHNEW,HHOLD,HTEMP,BETA,THOLD,THNEW,CON,CAP,thorg,thmod,hh
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS ::&
            &TTATM,PRECI,RSO,HT,RR,HCA
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: CONTAB,CAPTAB,THETAB
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS ::HSAT
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: UNSFLUX,UNSFLUXV
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: CUMQ
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS :: PRTIME
        INTEGER,DIMENSION(:,:),POINTER,CONTIGUOUS :: KZON

        !--------------------------------------------------------------------
        !DEFINING THE PARAMETERS AND VARIABLES FOR SOLUTE TRANSPORT

        LOGICAL, DIMENSION(:),POINTER,CONTIGUOUS ::lLinear,lMobIm
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS ::cprev

        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS ::s
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS ::cNew,cvBot,cvCh0,&
            &cvCh1,cvChIm,cvChR,cvTop,g0,g1,q0,q1,Retard, SinkIm,sSink,STrans,&
            &TempN,TempO,B,ThNIm,ThOIm,wc  ,Disp ,thnewim,tholdim,cPrevO,cTemp,&
            &D,E,F,p,r,q,cRootMax
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS ::CHpar,Wdep
        INTEGER :: MaxItC, npar2,ierr
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS :: C!,vo,vn
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS :: ths,thsat,n,xmasschange
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS :: tdep,cctop,ccbot
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: CumCh, CTOP,CBOT
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: watin
        DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS :: VOld,VNew ,vmid,vc,SorbN,SorbN2,TEMP,concf
        DOUBLE PRECISION,DIMENSION(:,:,:,:),POINTER,CONTIGUOUS :: DMOIST
        DOUBLE PRECISION,DIMENSION(:,:,:),POINTER,CONTIGUOUS ::concnew1,conc,sorb,sorb2,concmod,concnew,CT,CB,concold

        !Chnaging the dimension for variables for including the solute transport subroutine
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: MATNUM_S
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: VO_S,VN_S,x_s
        DOUBLE PRECISION,POINTER:: b1
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: THETA_S,VELOC
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: wvol
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: area,subvol
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: ccuma,cvoli,solin
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: ccumt ,cRunOff,cgwl
        DOUBLE PRECISION,DIMENSION(:),POINTER,CONTIGUOUS:: solfluxwat
        !--------------------------------------------------------------------
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
    !INTEGER,PARAMETER :: NMATD=50
    CHARACTER(LEN=700) :: LINE
    INTEGER :: MAXNP
    INTEGER :: LLOC,ISTART,ISTOP
    REAL :: O
    !      Body
    !      Allocate variables of known size
    CALL GWF2UNSF1PREAL()
    !      Identify package
    WRITE(IOUT,FMT1) IN
    !      Read Comments and print them, return the first non-comment line
    !      Reads the logical variable for solute tranport and the no of solute; modified by Sahila
    !      --------------------------------------------------------
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    lChem=.false.
    READ(IN,*) lChem! Switch for  solute transport read from the first line in the .uns file
    IF (LCHEM)  then
        WRITE(IOUT,*) 'The Solute transport simultaion is considered '
    ELSE
        WRITE(IOUT,*) 'The Solute transport simulation is not considered '
    end if

    !     -----------------------------------------------------------
    CALL UHRCOM(IN,IOUT,LINE)

    !      Read Profile information
    LLOC = 1
    CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPUNSF ,o,IOUT,IN)
    CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUNSFOP,o,IOUT,IN)
    CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUNSFCB,o,IOUT,IN)
    CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUNSFPR,o,IOUT,IN)
    !      Print profile information
    IF(NPUNSF.GT.0) THEN
        IF(LSTCHK(3)) WRITE (IOUT,FMT2) NPUNSF
    ELSE
        NPUNSF=0
        WRITE(IOUT,'(1X,A)') ERRM1
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
    WRITE(IOUT,*) NMAT,MAXNP,MAXATM
    ! modified by AdamS, uninitalized NTAB produced runtime errors
    NTAB = 100
    NUMNPD=MAXNP
    IF(NMAT.GT.NMATD)THEN
        WRITE(IOUT,'(1X,A)') ERRM3
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
    ALLOCATE (HH (NUMNPD,NPUNSF))
    ALLOCATE (THETA(NUMNPD,NPUNSF))
    ALLOCATE (THOLD(NUMNPD,NPUNSF))
    ALLOCATE (THORG(NUMNPD,NPUNSF))
    ALLOCATE (THMOD(NUMNPD,NPUNSF))
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
    ALLOCATE (TIME(NPUNSF))
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
    ALLOCATE (POPTM (NPUNSF,nmat))
    ALLOCATE (P2H   (NPUNSF))
    ALLOCATE (P2L   (NPUNSF))
    ALLOCATE (P3    (NPUNSF))
    ALLOCATE (R2H   (NPUNSF))
    ALLOCATE (R2L   (NPUNSF))
    ALLOCATE (VBOT  (NPUNSF))
    ALLOCATE (SUMVBOT(NPUNSF))
    ALLOCATE (ZONEFLUX(NPUNSF))
    ALLOCATE (CUMQ  (12,NPUNSF))
    allocate (wc(numnpd,NPUNSF))

    !     ------------------------------------------------------------------
    !Allocation of Solute parameters and variables
    ALLOCATE (LMOBIM(NMAT))
    ALLOCATE (LLINEAR(NSD))
    ALLOCATE (CHpar(NSD*16+4,NMat))
    ALLOCATE (Wdep(2+NMatD,NSD*9))
    ALLOCATE (ctop(npunsf,nsd))
    ALLOCATE (CBOT(npunsf,nSD))
    ALLOCATE (cctop(nsd))
    ALLOCATE (CcBOT(nSD))
    ALLOCATE (TDEP(NSD*16+4))
    allocate (Cumch(10,nsd))
    allocate (conc(NSD,NUMNPD,NPUNSF))
    allocate (concnew1(NSD,NUMNPD,NPUNSF))
    allocate (concold(NSD,NUMNPD,NPUNSF))
    allocate (concnew(nsd,numnpd,NPUNSF))
    allocate (concmod(nsd,numnpd,NPUNSF))! added for modifying the concentration
    allocate (sorb(NSD,numnpd,NPUNSF))
    allocate (sorb2(NSD,numnpd,NPUNSF))
    Allocate (C(numnpd))
    allocate (S(numnpd))
    ALLOCATE (CT(NPUNSF,NSD,MAXATM))
    ALLOCATE (cB(NPUNSF,NSD,MAXATM))
    ALLOCATE (VOld(NUMNPD,NPUNSF))
    ALLOCATE (VNew(NUMNPD,NPUNSF))
    allocate (vmid(numnpd,npunsf))
    allocate (vc(numnpd,npunsf))
    ALLOCATE (B(NUMNPD,npunsf))
    ALLOCATE (cPrevO(NUMNPD,npunsf))
    ALLOCATE (cPrev(NUMNPD))

    ALLOCATE (cTemp(NUMNPD,npunsf))
    ALLOCATE (F	(NUMNPD,npunsf))
    ALLOCATE (g0(NUMNPD,npunsf))
    ALLOCATE (g1(NUMNPD,npunsf))
    ALLOCATE (q0(NUMNPD,npunsf))
    ALLOCATE (q1(NUMNPD,npunsf))
    ALLOCATE (SorbN(NUMNPD,npunsf))
    ALLOCATE (SorbN2(NUMNPD,npunsf))
    ALLOCATE (sSink(NUMNPD,npunsf))
    ALLOCATE (STrans(NUMNPD,npunsf))
    ALLOCATE (TempO(NUMNPD,npunsf))
    ALLOCATE (cRootMax(npunsf,NSD))
    ALLOCATE (cvBot(NSD,npunsf))
    ALLOCATE (cvCh0(NSD,npunsf))
    ALLOCATE (cvCh1(NSD,npunsf))
    ALLOCATE (cvChIm(NSD,npunsf))
    ALLOCATE (cvChR(NSD,npunsf))
    ALLOCATE (cvTop(NSD,npunsf))
    ALLOCATE (DMOIST(NMatD,NSD,13,6))
    ALLOCATE (THSAT(NMATD))
    allocate (ths(nmatd))
    ALLOCATE (MATNUM_S(NUMNPD))
    ALLOCATE (VN_S(NUMNPD))
    ALLOCATE (VO_S(NUMNPD))
    allocate (x_s(numnpd))
    ALLOCATE (THETA_S(NUMNPD))
    ALLOCATE (VELOC(NUMNPD))
    ALLOCATE (TEMP(NUMNPD,NPUNSF))
    Allocate (concf(numnpd,NPUNSF))
    allocate (p(numnpd,NPUNSF))
    allocate (q(numnpd,NPUNSF))
    allocate (r(numnpd,NPUNSF))
    allocate (s(numnpd))
    allocate (thoim(numnpd,NPUNSF))
    allocate (thnim(numnpd,NPUNSF))
    allocate (thnewim(numnpd,NPUNSF))
    allocate (tholdim(numnpd,NPUNSF))
    allocate (d(numnpd,NPUNSF))
    allocate (disp(numnpd,NPUNSF))
    allocate (E(numnpd,NPUNSF))
    allocate (n(numnpd))
    allocate (xmasschange(npunsf))
    allocate (retard(numnpd,NPUNSF))
    allocate (sinkim(numnpd,NPUNSF))
    allocate (tempn(numnpd,NPUNSF))
    allocate (cnew(numnpd,NPUNSF))
    allocate (area(10))
    allocate (subvol(10))
    allocate (ccuma(nsd))
    allocate (cgwl(nsd))
    allocate (ccumt(nsd))
    allocate (cRunOff(nsd))
    allocate (cvoli(nsd))
    allocate (solin(numnpd))
    allocate (solfluxwat(npunsf))
    CALL GWF2UNSF1FOPEN()
    CALL GWF2UNSF1RPP(IN,IEVT,INRCH)
    CALL GWF2UNSF1PSV(IGRID)
    !-----------------------------------------------------------------------

    END SUBROUTINE

    SUBROUTINE GWF2UNSF1RPP(IN,INEVT,INRCH)
    !******************************************************************
    !Read and Prepare Parameters
    !******************************************************************
    USE GLOBAL, ONLY: IOUT,NPER,PERLEN
    USE GWFUNSFMODULE,ONLY:DTINIT,DT,PROPR,PROINF,LINITW,SINKF,THOLD
    USE GWFUNSFMODULE,ONLY:INFTOP,INFBOT,ATMBC,TATM,TMAX,PRTIME,TINIT
    USE GWFUNSFMODULE,ONLY:NPUNSF,THETA,LPRINT,LPTIME,LOUTF,&
        &LHEAD!,IPROU
    USE GWFUNSFMODULE
    !---------------------------------------------------------------------
    !Solute variables: Added by Sahila
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !Variables
    INTEGER,INTENT(IN) :: IN,INEVT,INRCH
    CHARACTER(LEN=100) :: LINE
    INTEGER :: I,K
    !CHECK: In the original code, these variables are arrays, however
    !they are not used, so I have replaced them with DOUBLE PRECISIONs
    DOUBLE PRECISION :: ADEPTH, AHEAD
    !CHECK: ZONARR is a variable present in the original code, at
    !line 439, but not DOUBLE PRECISIONly used for anything - removed
    !Body
    !The following line is required to make certain subroutines work
    !correctly.

    !if (lchem) then
    !    open ( unit = 101, file = 'HYD_Conc_GWT.OUT' )
    !    open ( unit = 103, file = 'HYD_ConFlux_GWT.OUT' )
    !end if

    !open(unit=98,FILE='HYD_ReadMe.out')
    !CALL HYD_README()

    THETA=>THOLD
    !CHECK:
    !Initialize
    CALL GWF2UNSF1INI()
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    CALL GWF2UNSF1INIT(IN,INEVT,INRCH)
    !Read Time information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    !Note: Requires KPER, which is 1 at this point
    CALL GWF2UNSF1TMIN(IN,1)
    DTINIT=DT
    !Read Material information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    CALL GWF2UNSF1MATIN(IN)
    !Generate Hydraulic Properties
    CALL GWF2UNSF1GENMAT()
    !Read profile printing information
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    !----------------------------------------------------------------------
    !Read the Input Paramters for the solute transport simulation:Added by Sahila

    if(lChem) then
        call ChemIn(IN)
    end if
    !-------------------------------------------------------------------
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN

    READ(IN,*) PROPR
    IF (PROPR.GT.0) THEN
        LPRINT=.TRUE. ! modified by AdamS
        ALLOCATE(PRTIME(PROPR))
        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        READ(IN,*) PROINF
        IF (PROINF.GT.0) THEN ! modified by AdamS
            LHEAD=.FALSE.
        ELSE
            LHEAD=.TRUE.
        END IF
        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        !CHECK: The original code use a temporary variable here
        !It was removed, since it isn't modified during READ
        !See line 505 of the original code for more information

        READ(IN,*) (PRTIME(I),I=1,PROPR)
    ELSE ! modified by AdamS, print only the final profile
        PROPR=1
        LPRINT=.TRUE.
        ALLOCATE(PRTIME(1))
        PRTIME(1)=0.0
        DO I=1, NPER
            PRTIME(1)=PRTIME(1)+PERLEN(I)
        END DO

        !PRTIME=>NULL()
    END IF
    !----------------------------------------------------------------------
    !For each profile ...
    DO I=1,NPUNSF
        !Read basic information
        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        CALL GWF2UNSF1BASINF(I,IN)
        !Read nodal point information
        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        !-----------------------------------------------------------------
        !Including nodal solute information: Added Sahila
        IF(.NOT.LCHEM) THEN

            CALL GWF2UNSF1NODINF(I,IN)
        ELSE

            CALL GWF2UNSF1NODINFC(I,IN)

        END IF
        !------------------------------------------------------------------
        !CALL GWF2UNSF1NODINF(I,IN)
        !CHECK: GWD has been modified to have a different behavior if
        !KPER=0
        !The original code contained a separate function GWDEPTH1,
        !which was used here exclusively. The only differences
        !between GWDEPTH and GWDEPTH1 are as follows:
        !The first line ( T=T+TD ) does not appear in GWDEPTH1
        !No assignments are made to HBOT(HBOT=HB) in GWDEPTH1
        !GWDEPTH  is located at line 1686 in the original code
        !GWDEPTH1 is located at line 2778 in the original code
        CALL GWF2UNSF1GWD(I,0,0,ADEPTH,AHEAD)
        !Read initial condition given in terms of water content
        IF(LINITW(I).GE.0)THEN
            CALL GWF2UNSF1INITW(I)
        END IF
        !Determine nodal values of hydraulic properties
        THETA=>THOLD
        CALL GWF2UNSF1SETMAT(I)
        !Read root water uptake information
        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        !----------------------------------------------------------
        ! Sink information for Solute transport: Added by Sahila
        IF (SINKF(I).GE.0.and.lchem) then
            CALL GWF2UNSF1SINKINc(I,IN)
        else if (SINKF(I).GE.0.and..not.lChem)then
            CALL GWF2UNSF1SINKIN(I,IN)
        END IF
        !----------------------------------------------------------

        !IF (SINKF(I).GE.0) THEN
        !    CALL GWF2UNSF1SINKIN(I,IN)
        ! END IF
        !CALL UHRCOM(IN,IOUT,LINE)
        !BACKSPACE IN
        !       Read atmospheric information
        !-------------------------------------------------------------------

        CALL UHRCOM(IN,IOUT,LINE)
        BACKSPACE IN
        IF (.not.lchem) then
            IF(INFTOP(I).GE.0.OR.INFBOT(I).GE.0.OR.ATMBC(I).GE.0) THEN
                CALL GWF2UNSF1READBC(I,in)
            ELSE
                TATM(I)=TMAX(I)
            END IF
        else
            IF(INFTOP(I).GE.0.OR.INFBOT(I).GE.0.OR.ATMBC(I).GE.0) THEN
                CALL GWF2UNSF1READBCC(I,in)
            ELSE
                TATM(I)=TMAX(I)
            END IF
        end if
        !-----------------------------------------------------------------------------
        !IF (INFTOP(I).GE.0 .OR. INFBOT(I).GE.0 .OR. ATMBC(I).GE.0) THEN
        !    CALL GWF2UNSF1READBC(I,IN)
        !ELSE
        !    TATM(I)=TMAX(I)
        !END IF
        !       Calculate Root water extraction rate
        THETA=>THOLD
        IF (SINKF(I).GE.0)then
            if (.not.lchem)then
                CALL GWF2UNSF1SETSNK1(I)
                if (lchem)then
                    call gwf2unsf1setsnkc1(I)
                END IF
            end if
        end if
        ! plevel(I)=0

        ! t(i)=0
        !  CALL GWF2UNSF1WVOLI(I)
        !--------------------------------------------------------------------
        !To find the Vold: Added by  Sahila
        !thnew=thold
        THETA=>THOLD
        CALL GWF2UNSF1SETMAT(I)
        thold=theta

        CALL GWF2UNSF1VELOC(I)
        vnew=vold
        thnew=thold

        !----------------------------------------------------------------------
    END DO

    ! Output style changed by AdamS
    !      Print information if requested
    !       IF(PROPR.GT.0) THEN
    !  DO K=1,PROPR
    !         IF(ABS(PRTIME(K)-TINIT(1)).LT.0.005)THEN
    !          IF(LOUTF)THEN
    !           LPTIME=.TRUE.
    !           IF(LPRINT)WRITE(IPROU,!) ' Depth (vs) Heads/ Water content'
    !          END IF
    !CALL GWF2UNSF1NODOUT1(TINIT(1),LHEAD)
    !          EXIT
    !         END IF
    !        END DO
    !       END IF
    write(77,"(a12,3x,a4,8x,a7,7x,a9,5x,a9,5x,a13)") "Profile NO:","Time","Conc_WT","Conc_flux","Time step", "Stress period"
    write(77,"(4x,a3,8x,a3,9x,a7,7x,a9,5x,a7,5x,a7)") "[-]" ,"[T]" ,"[M/L3]", "[M/L2/T]" ,"[-]","[-]"

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
    ! modified by Adam Szymkiewicz, get rid of trailing tabs
    ! when the input file is imported from Excel
    DO I=1, 10
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
    SUBROUTINE GWF2UNSF1VELOC(PID)
    !USE GLOBAL, ONLY:IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY:  NUMNP, CON, HHNEW,DT,THNEW,THOLD,SINK,X,VOld,&
        &numnpd,npunsf
    !IMPLICIT NONE
    INTEGER,INTENT(IN) :: PID
    INTEGER :: I,M
    DOUBLE PRECISION ::DXN,DXA,DXB,DX1,VA,VB,xx(numnpd,npunsf)
    
    CHARACTER(LEN=*),PARAMETER::FMT2="(I4,E24.4, F8.3,E12.4)"

    do i=1,numnp(pid)
        xx(i,pid)=x(i,pid)
    end do


    M=(NUMNP(PID)-1)
    DXN=Xx(NUMNP(PID),PID)-Xx(M,PID)
    DX1=Xx(2,PID)-Xx(1,PID)

    FRE=1.
    GRAV=1.

    DO 11 I=2,(NUMNP(PID)-1)
        DXA=  Xx(I+1,PID)-Xx(I,PID)
        DXB=  Xx(I,PID)-Xx(I-1,PID)
        VA=-(CON(I,PID)+CON(I+1,PID))/2.*((HHNEW(I+1,PID)-HHNEW(I,PID))/DXA+1.)
        VB=-(CON(I,PID)+CON(I-1,PID))/2.*((HHNEW(I,PID)-HHNEW(I-1,PID))/DXB+1.)
        vold(I,PID)=(VA*DXB+VB*DXA)/(DXA+DXB)

11  CONTINUE

   ! vold(Numnp(pid),pid)=-(Con(NUMNP(PID),PID)+Con(M,PID))/2.*((hhNew(NUMNP(PID),PID)-hhNew(M,PID))/dxN+1*1)-&
   !     &     dxN/2.*(1*(Thold(Numnp(pid),pid)-ThOld(Numnp(pid),pid))/dt(pid)+Sink(Numnp(pid),pid))

   ! vold(1,pid)=-(Con(1,pid)+Con(2,pid))/2.*((((hhNew(2,pid)-hhNew(1,pid))/dx1)+1))+&
   !     &      dx1/2.*(fre*(Thold(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))


    vold(Numnp(pid),pid)=-(Con(Numnp(pid),pid)+Con(M,pid))/2.*((hhNew(Numnp(pid),pid)-hhNew(M,pid))/dxN+1.)-dxN/2.*(1.*(ThOld(Numnp(pid),pid)-ThOld(Numnp(pid),pid))/dt(pid)+Sink(Numnp(pid),pid))
    vold(1,pid)=-(Con(1,PID)+Con(2,pid))/2.*((hhNew(2,pid)-hhNew(1,pid))/dx1+1.)+dx1/2.*(1.*(ThOld(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))
12  continue
    return
    END SUBROUTINE


    SUBROUTINE GWF2UNSF1VELOCAD(PID)
    !USE GLOBAL, ONLY:IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY:  NUMNP,x, CON, HHNEW,DT,THNEW,THOLD,SINK,numnpd,npunsf,VNew,vbot
    !IMPLICIT NONE
    INTEGER,INTENT(IN) :: PID
    INTEGER :: I,M
    double precision  ::DXN,DXA,DXB,DX1,VA,VB,fre,grav,xx(numnpd,npunsf)
    CHARACTER(LEN=*),PARAMETER::FMT2="(I4,E24.4, F8.3,E12.4)"
    
     
     
    do i=1,numnp(pid)
        j=numnp(pid)
        xx(i,pid)=x(i,pid)
    end do

    M=NUMNP(PID)-1
    DXN=xx(NUMNP(PID),PID)-Xx(M,PID)
    DX1=Xx(2,PID)-Xx(1,PID)

    DO 11 I=2,(NUMNP(PID)-1)
        DXA=  Xx(I+1,PID)-Xx(I,PID)
        DXB=  xX(I,PID)-Xx(I-1,PID)
        
        VA=-(CON(i,PID)+CON(i+1,PID))/2.* ((Hhnew(i+1,pid)-HHnew(i,pid))/Dxa+1.)
       
       vB=-(Con(i,pid)+Con(i-1,pid))/2.*((hhNew(i,pid)-hhNew(i-1,pid))/dxB+1.)
        
        VNew(I,PID)=(VA*DXB+VB*DXA)/(DXA+DXB)
         
          
11  CONTINUE

    vnew(numnp(pid),pid)=-(Con(NUMNP(PID),pid)+Con(m,pid))/2.*((hhNew(NUMNP(PID),pid)-hhNew(m,pid))/dx1+1.)+&
        &      dx1/2.*(1.*(Thnew(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))

    vnew(1,pid)=-(Con(1,pid)+Con(2,pid))/2.*((hhNew(2,pid)-hhNew(1,pid))/dx1+1.)+&
        &      dx1/2.*(1.0*(ThNew(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))

    vbot(pid)= vnew(1,pid)

    
    
    return
    END SUBROUTINE
    SUBROUTINE GWF2UNSF1SOLCONC(PID,kstp,kper)
    !USE GLOBAL, ONLY:IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY:  conc,tmax,solfluxwat,VC, NUMNP,x, CON, HHNEW,DT,THNEW,THOLD,SINK,numnpd,npunsf,VNew,vbot
    !IMPLICIT NONE
    INTEGER,INTENT(IN) :: PID, KPER,KSTP
    INTEGER :: I,M,k
    DOUBLE PRECISION ::DXN,DXA,DXB,DX1,VA,VB,vd,xx(numnpd,npunsf)
    CHARACTER(LEN=*),PARAMETER::FMT2="(I4,E24.4, F8.3,E12.4)"

   
    do i=1,numnp(pid)
        j=numnp(pid)
        xx(i,pid)=x(i,pid)
    end do


    M=NUMNP(PID)-1
    DXN=xx(NUMNP(PID),PID)-Xx(M,PID)
    DX1=Xx(2,PID)-Xx(1,PID)


    DO 11 I=2,(NUMNP(PID)-1)
        DXA=  Xx(I+1,PID)-Xx(I,PID)
        DXB=  xX(I,PID)-Xx(I-1,PID)
        VA=-(CON(I,PID)+CON(I+1,PID))/2.*((HHNEW(I+1,PID)-HHNEW(I,PID))/DXA+1.)
        VB=-(CON(I,PID)+CON(I-1,PID))/2.*((HHNEW(I,PID)-HHNEW(I-1,PID))/DXB+1.)
        VC(I,PID)=(VA*DXB+VB*DXA)/(DXA+DXB)

11  CONTINUE

    vC(numnp(pid),pid)=-(Con(NUMNP(PID),pid)+Con(m,pid))/2.*((hhNew(NUMNP(PID),pid)-hhNew(m,pid))/dx1+1.*1.)+&
        &      dx1/2.*(1.*(Thnew(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))
    vC(1,pid)=-(Con(1,pid)+Con(2,pid))/2.*((hhNew(2,pid)-hhNew(1,pid))/dx1+1.*1.)+&
        &      dx1/2.*(1.*(ThNew(1,pid)-ThOld(1,pid))/dt(pid)+Sink(1,pid))

    return

    END SUBROUTINE




    SUBROUTINE GWF2UNSF1PROFUP(PID,HDRY)
    ! *************************************************
    ! modified by AdamS, according to the method proposed by Sahila Beegum
    ! update pressure profile in the unsaturated zone
    ! so that the initial bottom flux in the new Modflow time step
    ! matches the average bottom flux from the previous time step
    ! in order to eliminate spurious flux oscillations
    ! ***************************************************

    USE GWFUNSFMODULE, ONLY: HHOLD,HHNEW,HTEMP,HB,X,XSURF,NUMNP,NUMNPD,&
        THOLD,THNEW,CON,PAR,MATNUM,HTAB,CONTAB,NTAB,VBOT

    IMPLICIT NONE
    INTEGER,INTENT(IN):: PID
    DOUBLE PRECISION,INTENT(IN):: HDRY ! critical value of pressure for dry soil
    DOUBLE PRECISION, ALLOCATABLE:: HH1(:),HP(:),HZ(:),ZZ(:),ZLAY(:),QQ1(:),QQ2(:),QQ3(:),&
        KPREV(:),KNEXT(:)
    DOUBLE PRECISION:: QFLUX,KK1,KK2,HC1,HC2,DZ12,DZ,ZGWT,ZMAX,ZMIN,FF1,FF2,KS,HK,&
        DH,ZC,QFLUX1,MAXQ,MINQ
    INTEGER:: J,INODE,INEXT,ISTEP,IDGW,NLAY,MAT1,MAT2,FLAG,IPROF,K,NODEGW
    INTEGER,ALLOCATABLE:: MATLAY(:)
    INTEGER,PARAMETER:: NSTEPS=100,NMAX=10000

    ALLOCATE(HH1(NUMNPD),HP(NUMNPD),HZ(NMAX),ZZ(NMAX),MATLAY(NUMNP(PID)),ZLAY(NUMNP(PID)))
    ALLOCATE(QQ1(NUMNP(PID)),QQ2(NUMNP(PID)),QQ3(NUMNP(PID)),KPREV(NUMNP(PID)),KNEXT(NUMNP(PID)))
    !  HHOLD(1,PID)=HB(PID)

    HH1=HHOLD(:,PID)

    DO J=1,NUMNP(PID)-1
        DZ=X(J+1,PID)-X(J,PID)
        KPREV(J)=(CON(J,PID)+CON(J+1,PID))/2.
        QQ1(J)=-(CON(J,PID)+CON(J+1,PID))/2.*&
            &((HH1(J+1)-HH1(J))/DZ+1.)
    
    END DO
     write (98,*)'----------------------'
    !  DZ=X(2,PID)-X(1,PID)
    !  QFLUX=-(CON(1,PID)+CON(2,PID))/2.*&
    !            &((HH1(2)-HH1(1))/DZ+1.)
    QFLUX = VBOT(PID)
    PRINT *, 'QFLUX=',QFLUX

    HP=HH1
    HP(1)=HB(PID)
    IF (HP(1)<0.0) THEN
        PRINT *,'Negative water pressure at the bottom of profile ',PID
        STOP
    ENDIF

    ! identify material layers
    NLAY=1
    MATLAY(1)=MATNUM(1,PID)
    DO J=1, NUMNP(PID)-1
        IF (MATNUM(J,PID)/=MATNUM(J+1,PID)) THEN
            ZLAY(NLAY)=0.5*X(J,PID)+X(J+1,PID)
            NLAY=NLAY+1
            MATLAY(NLAY)=MATNUM(J+1,PID)
        ENDIF
    ENDDO
    ZLAY(NLAY)=X(NUMNP(PID),PID)

    ! compute profile below groundwater table, find GWT position
    DO INODE=1,NUMNP(PID)-1
        MAT1=MATNUM(INODE,PID)
        KK1=PAR(5,MAT1)
        DZ12=0.5*(X(INODE+1,PID)-X(INODE,PID))
        HC1=HP(INODE)-(QFLUX/KK1+1.)*DZ12
        IF (HC1<0.0) THEN
            ZGWT=X(INODE,PID)+HP(INODE)/(1.+QFLUX/KK1)
            EXIT
        ENDIF
        MAT2=MATNUM(INODE+1,PID)
        KK2=PAR(5,MAT2)
        HC2=HC1-(QFLUX/KK2+1.)*DZ12
        IF (HC2<0.0) THEN
            ZGWT=X(INODE,PID)+DZ12+HC1/(1.+QFLUX/KK2)
            EXIT
        ELSE
            HP(INODE+1)=HC2
            QFLUX1=-KK1*((HP(INODE+1)-HP(INODE))/(2.*DZ12)+1.)
            !PRINT *, X(INODE+1,PID), HP(INODE+1), QFLUX1
            !READ(*,*)
        ENDIF
    END DO

    NODEGW=INODE+1 ! save the first node above groundwater table

    DO J=1,NLAY
        IF (ZGWT<ZLAY(J)) THEN
            IDGW=J
            EXIT
        ENDIF
    ENDDO

    PRINT *, 'IDGW= ',IDGW, ' ZGWT= ',ZGWT
    !  PRINT *, ZLAY

    ! compute profile above the groundwater table

    INODE=1
    HZ(INODE)=0
    ZZ(INODE)=ZGWT

    DO INEXT=IDGW,NLAY !loop over soil layers

        MAT1=MATLAY(INEXT)
        ZMAX=ZLAY(INEXT)
        KK1 = FINTERPOL(HTAB,CONTAB(:,MAT1),HZ(INODE),NTAB)
        PRINT *,'HZ=',HZ(INODE),' KK=',KK1
        FF1 = -1.0/(1.0+(QFLUX/KK1))
        KS=PAR(5,MAT1)

        FLAG=0;
        IF (QFLUX<-KS) THEN
            IF (HZ(INODE)>=-1e-5) THEN ! saturated profile
                FLAG=2
            ELSE ! unsat-sat transition possible
                FLAG=1;
                DH=(-1e-5-HZ(INODE))/NSTEPS
            ENDIF
        ELSEIF (QFLUX<0) THEN! unsaturated profile tending towards constant h
            HK=FINTERPOL(CONTAB(:,MAT1),HTAB,ABS(QFLUX),NTAB)
            PRINT *,'Q= ',QFLUX,' HK=',HK
            DH=0.99999*(HK-HZ(INODE))/NSTEPS
            FLAG=0
        ELSE ! evaporation
            IF (ABS((HZ(INODE)-HDRY)/HDRY)>0.0001) THEN
                FLAG=-1;
                DH=-0.01*(X(NUMNP(PID),PID)-X(1,PID))
            ELSE  ! dry profile
                FLAG=-2;
            ENDIF
        ENDIF

        IF (FLAG==-2) THEN
            ! pressure at critical dry value
            HZ(INODE+1)=HDRY
            ZZ(INODE+1)=ZMAX
            INODE=INODE+1
            EXIT

        ELSE IF (FLAG<2) THEN
            ISTEP=0;
            DO WHILE (.TRUE.)
                IF (INODE==NMAX-1) THEN
                    PRINT *, 'Max. number of points in profile reconstruction exceeded'
                    STOP
                ENDIF
                ISTEP=ISTEP+1
                HC1=HZ(INODE)+DH
                KK2 = FINTERPOL(HTAB,CONTAB(:,MAT1),HC1,NTAB)
                FF2 = -1.0/(1.0+(QFLUX/KK2))
                DZ = (FF1+FF2)/2.*DH
                ZC=ZZ(INODE)+DZ
                IF (ZC>ZMAX) THEN ! z above the top of layer
                    ZZ(INODE+1)=ZMAX
                    HZ(INODE+1)=HZ(INODE)+(HC1-HZ(INODE))*(ZMAX-ZZ(INODE))/(ZC-ZZ(INODE))
                    INODE=INODE+1
                    EXIT
                ELSEIF ((FLAG==-1).AND.((HZ(INODE)+DH)<HDRY)) THEN
                    ! critical dry pressure exceeded in evaporation
                    HZ(INODE+1)=HDRY
                    ZZ(INODE+1)=ZMAX
                    INODE=INODE+1
                    EXIT
                ELSE
                    ZZ(INODE+1)=ZC
                    HZ(INODE+1)=HC1
                    FF1=FF2
                    INODE=INODE+1
                    IF (FLAG==-1) DH=1.05*DH ! increase DH for evaporation
                ENDIF
                IF ((ISTEP==NSTEPS).AND.(FLAG>=0)) EXIT
            ENDDO
        ENDIF

        IF ((FLAG>0).AND.(ZZ(INODE)<ZMAX)) THEN
            ! integrate over the saturated part (perched water)
            DZ=ZMAX-ZZ(INODE)
            DH=-DZ*(1+QFLUX/KS)
            ZZ(INODE+1)=ZMAX
            HZ(INODE+1)=HZ(INODE)+DH
            INODE=INODE+1;
        ENDIF

    ENDDO

    ! find the nodal pressure values by interpolation

    IPROF=1
    DO J=NODEGW,NUMNP(PID)
        DO K=IPROF,INODE-1
            IF ((X(J,PID)>=ZZ(K)).AND.(X(J,PID)<=ZZ(K+1))) THEN
                HP(J)=HZ(K)+(HZ(K+1)-HZ(K))*(X(J,PID)-ZZ(K))/(ZZ(K+1)-ZZ(K))
                IPROF=K
                EXIT
            ENDIF
        ENDDO
    ENDDO

    DO J=1,NUMNP(PID)-1
        DZ=X(J+1,PID)-X(J,PID)
        KK1 = FINTERPOL(HTAB,CONTAB(:,MATNUM(J,PID)),HP(J),NTAB)
        KK2 = FINTERPOL(HTAB,CONTAB(:,MATNUM(J+1,PID)),HP(J+1),NTAB)
        KPREV(J)=(KK1+KK2)/2.
        QQ2(J)=-(KK1+KK2)/2.*&
            &((HP(J+1)-HP(J))/DZ+1.)
    END DO

    ! update the bottom part of Hydrus profile
    HHOLD(1,PID)=HP(1)
    DO J=1,NUMNP(PID)-1
        HHOLD(J+1,PID)=HP(J+1)
        ! compute flux in the Hydrus transient profile
        DZ=X(J+1,PID)-X(J,PID)
        QFLUX1=-(CON(J,PID)+CON(J+1,PID))/2.*&
            &((HH1(J+1)-HH1(J))/DZ+1.)
        ! the flux
        MAXQ=MAX(ABS(QFLUX),ABS(QFLUX1))
        MINQ=MIN(ABS(QFLUX),ABS(QFLUX1))
        IF (QFLUX1*QFLUX<0.0) EXIT
        IF (ABS((QFLUX1-QFLUX)/MAXQ)>0.001*MAXQ) EXIT
    ENDDO

    !HHOLD(:,PID)=HP(:)
    HHNEW(:,PID)=HHOLD(:,PID)
    HTEMP(:,PID)=HHOLD(:,PID)
    CALL GWF2UNSF1SETMAT(PID)
    THOLD(:,PID)=THNEW(:,PID)

    DEALLOCATE(HH1,HP,HZ,ZZ,MATLAY,ZLAY)

    CONTAINS
    DOUBLE PRECISION FUNCTION FINTERPOL(XT,YT,XP,NTAB)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN):: XT(:),YT(:),XP
    INTEGER,INTENT(IN):: NTAB
    INTEGER:: I

    !	FINTERPOL=1.
    IF (XP>=XT(1)) THEN
        FINTERPOL=YT(1)
    ELSE IF (XP<=XT(NTAB)) THEN
        FINTERPOL=YT(NTAB)
    ELSE
        DO I=1,NTAB-1
            IF ((XP<=XT(I)).AND.(XP>=XT(I+1))) THEN
                FINTERPOL=YT(I)+(YT(I+1)-YT(I))*(XP-XT(I))/(XT(I+1)-XT(I))
                EXIT
            ENDIF
        ENDDO
    ENDIF
    END FUNCTION


    END SUBROUTINE

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
    
    QFLUX = VBOT(PID)
    PRINT *, 'QFLUX1=',QFLUX

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
        IF (ABS(QQ1(J)-QQ2(J))>(1.0d-12+.1*MAXQ)) EXIT
    ENDDO

    !HHOLD(:,PID)=HP(:)

    HHNEW(:,PID)=HHOLD(:,PID)  ! this is the modified head
    HTEMP(:,PID)=HHOLD(:,PID)
    CALL GWF2UNSF1SETMAT(PID) ! this takes HHNEW and gives the valus of thnew
    THOLD(:,PID)=THNEW(:,PID)!  this is the modified thEta
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
           PRINT *, 'icount',icount
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
           ! PRINT *, 'icount',icount
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


    SUBROUTINE GWF2UNSF1AD(IGRID,KSTP,KPER)
    !     ******************************************************************
    !     calculate unsat flow (ADvance subroutine)
    !     ******************************************************************
    USE GWFBASMODULE,  ONLY: TOTIM,DELT
    USE GWFUNSFMODULE, ONLY: NPUNSF,TMAX, THETA,SINKF,PLEVEL,VBOT,SUMVBOT,ZONEFLUX,PERIMP
    USE GWFUNSFMODULE, ONLY: T,TATM,DT,TOLD,TLEVEL,DTOLD,DTMIN
    USE GWFUNSFMODULE, ONLY: INFTOP,INFBOT,ATMBC,ALEVEL,HCA
    USE GWFUNSFMODULE, ONLY: KTOLD,KBOLD,KODBOT,KODTOP
    USE GWFUNSFMODULE, ONLY: LPTIME,LPRINT,LHEAD,solfluxwat,VC
    USE GWFUNSFMODULE, ONLY: PROPR,PRTIME,vmid,con,epsi
    USE GWFUNSFMODULE, ONLY: HHOLD,HHNEW,HTEMP,HB,X,XSURF,NUMNP,THOLD,THNEW
    USE GWFUNSFMODULE, ONLY: LCHEM,KBOTCH,KTOPCH ,ths,tholdim,switch
    USE GWFUNSFMODULE, ONLY: num,ikod,s,r,p,q,vnew,vold,thnewim,nmat,beta,matnum,rbot,lwat
    USE GWFUNSFMODULE, ONLY: ths, sink,tlevel,ns,cmid,smid,tpulse,nsd,llinear,retard,thsat,solfluxwat
    USE GWFUNSFMODULE, ONLY: numnp,numnpd,conc,hh,wvoli,xmasschange,WCUMT,WCUMA,THORG,volume

    use msflib
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: IGRID,KSTP,KPER
    DOUBLE PRECISION :: ADEPTH,AHEAD,Rtime1,RTime,RTime2 ,dxa,dxb,va,vb
    INTEGER :: I,K,SKIP,INOD,nlay,l,Tswitch,m
    integer*2 i2,iYear,iMonth,iDay,iHours,iMins,iSecs,i100th
    !      Body
    CALL GWF2UNSF1PNT(IGRID)
    !      The following line is required to make some subroutines work
    !      correctly
    THETA=>THNEW
    !     CHECK: Nonsensical line in the original code at this point
    !            See line 687 of the original code for more information
    !SKIP=0

    Tswitch=1
    DO I=1,NPUNSF  ! loop over soil profiles
        !IF(SKIP.EQ.0) THEN

        TMAX(I)=TOTIM
        IF (KSTP.EQ.1 .AND. KPER.EQ.1)   plevel(I)=1
        IF (KSTP.EQ.1 .AND. KPER.EQ.1) THEN
            CONTINUE
            !     CHECK: Why copy?
            !          IBOUND0(1:NCOL,1:NROW,1:NLAY) = IBOUND(1:NCOL,1:NROW,1:NLAY)
        ELSE
            CALL GWF2UNSF1GWD(I,KSTP,KPER,ADEPTH,AHEAD)
        END if
        ! modified by AdamS - adjustment of the pressure profile
        ! according to Sahila's method (steady state profile)


        CALL GWF2UNSF1PROFUP1(I,-ABS(HCA(I,ALEVEL(I))))
        !change the wvoli to adjust for the water flow balance
        ! masschange(i)=xMassChange(i)
        plevel(I)=0
        volume=0.

        CALL GWF2UNSF1SETMAT(I)
        THOLD(:,I)=THNEW(:,I)
        ! vold(:,I)=vnew(:,I)
        SUMVBOT(I)=0.

        !start calculating the time required for simulation of water flow and Solute transport using HYDRUS
        !runstart time
        call getdat(iYear,iMonth,iDay)
        call gettim(iHours,iMins,iSecs,i100th)
        Rtime1=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)
        switch=0

        !  initime=t(pid)
        ! modified by AdamS - the loop below was missing in the previous version

        DO  ! loop over time steps for a single profile
            ! solve unsaturated flow for a single time step in the profile


            CALL GWF2UNSF1WF(I)

            !Calculate time-averaged VBOT
            !CALL GWF2UNSF1TOTF(I)
            !----------------------------------------------------------------
            CALL GWF2UNSF1VELOCAD(I)
            !----------------------------------------------------------------

            SUMVBOT(I)=SUMVBOT(I)+VBOT(I)*DT(I)
            WRITE ( 97,'(4(a7,es22.8))') 'T= ', T(I), ' DT= ',DT(I),' VBOT= ', VBOT(I),'SUMV= ',SUMVBOT(I)

            IF(SINKF(I).GE.0..and..not.lChem)then
                THETA=>THNEW
                CALL GWF2UNSF1SETSNK2(I)
            else IF(SINKF(I).GE.0..and.lChem)then
                CALL gwf2unsf1SETSNKC2(I)

            ENDIF


            !Solute transport
            !---------------------------------------------------------------------------
            IF (lChem) then
                CALL GWF2UNSF1SOLUTE(i,kstp,kper,thold,thnew,vold,vnew,htemp)
            End if

            switch=1
            PLEVEL(I)=PLEVEL(I)+1
            !   Read time-dependent boundary condition
            IF(ABS(T(I)-TATM(I)).LE.0.001*DT(I)) THEN
                !   CHECK: Redundant IF statement at line 771 in the original code
                ! IF(INFTOP(I).GE.0 .OR. INFBOT(I).GE.0 .OR. ATMBC(I).GE.0..and..not.lChem) THEN
                IF(.not.lChem) THEN
                    CALL GWF2UNSF1SETBC1(I)
                ELSE
                    CALL GWF2UNSF1SETBCC1(I)
                END IF
                ALEVEL(I)=ALEVEL(I)+1
            END IF


            ! modified by AdamS
            HTEMP(:,I)=HHNEW(:,I)
            HHOLD(:,I)=HHNEW(:,I)
            THOLD(:,I)=THNEW(:,I)
            VOLD(:,I)=VNEW(:,I)

            TOLD(I)=T(I)
            DTOLD(I)=DT(I)
            KTOLD(I)=KODTOP(I)
            KBOLD(I)=KODBOT(I)


            IF(ABS(T(I)-TMAX(I)).LE.0.00001*DT(I)) THEN
                T(I)=T(I)+DT(I)
                TLEVEL(I)=TLEVEL(I)+1
                IF(TLEVEL(I).GT.1000000) TLEVEL(I)=2
                EXIT
            ELSE
                CALL GWF2UNSF1TMCONT(I)
                T(I)=T(I)+DT(I)
                TLEVEL(I)=TLEVEL(I)+1
                IF(TLEVEL(I).GT.1000000) TLEVEL(I)=2
            END IF
        END DO   ! loop over time steps in the profile

        call GWF2UNSF1SOLCONC(i,kstp,kper )

        do 35 k=1,(numnp(i)-1)
            do 45 l=1,ns
                if (hhnew(numnp(i),i).ge.0)then
                    solfluxwat(i)=vc(numnp(i),i)*conc(l,numnp(i),i)
                    write(77,"(I4,5X,f10.2,5X, E10.3,5X,E10.3,7x, I6,7x,I6)") i,Tmax(i),conc(l,numnp(i),i),solfluxwat(i),kstp,kper
                else if (hhnew(K,i).ge.0.and.hhnew(K+1,i).le.0)then
                    solfluxwat(i)=(conc(l,K,i)-((conc(l,K+1,i)-conc(l,K,i))*((hhnew(k,i))&
                        & /(hhnew(k+1,i)-hhnew(k,i)))))*((VC(K,i)-((VC(K+1,i)-VC(K,i))*((hhnew(k,i))/(hhnew(k+1,i)-hhnew(k,i))))))
                    write(77,"(I4,5X,f10.2,5X, E10.3,5X,E10.3,7x, I6,7x,I6)") i,Tmax(i),(conc(l,K,i)&
                        &-((conc(l,K+1,i)-conc(l,K,i))*((hhnew(k,i))/(hhnew(K+1,i)-hhnew(k,i))))), solfluxwat(i),kstp,kper
                end if
45          continue
35      continue

        ! run end time
        call getdat(iYear,iMonth,iDay)
        call gettim(iHours,iMins,iSecs,i100th)
        Rtime2=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)
        Call runtime(rtime1,rtime2,KPER, KSTP,I)

        ! AdamS: introduced ZONEFLUX variable facilitate calculations and output

        ZONEFLUX(I)=SUMVBOT(I)/DELT*((100.-PERIMP(I))/100.)
        VOLD(:,i)=VNEW(:,i)
        Tswitch=0
    END DO ! loop over soil profiles


    ! AdamS: modified output

    CALL GWF2UNSF1TLSUM

    !DO K=1,PROPR
    !    IF (ABS(PRTIME(K)-T(1)).LT.DT(1)) THEN
    !        CALL GWF2UNSF1NODOUT1(T(1),LHEAD)
    !        EXIT
    !    END IF
    !END DO

    !  LPTIME=.FALSE.
    !LPRINT=.TRUE.
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
    ! AdamS: the line below was apparently omitted when porting the code
    ! it is crucial for the calculation of exchange flux
    TOTIMOLD=TOTIM
    ICOUNT=0
    !      Apply SUMVBOT to each cell
    ! AdamS: use ZONEFLUX instead of SUMVBOT
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
    DO IR=1,NROW
        DO IC=1,NCOL
            ! modified by AdamS: DELC depends on IR, DELR on IC
            UNSFLUXV(IC,IR)=UNSFLUX(IC,IR)*DELC(IR)*DELR(IC)
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
                RHS(IC,IR,IL) = RHS(IC,IR,IL)+UNSFLUXV(IC,IR)
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
        ! modified by AdamS : use ZONEFLUX to get fluxes for soil profiles
        ! use of UZO caused runtime errors
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
            DO N=1,NPUNSF
                !        K=UZO(2,N)
                !        KK=UZO(3,N)
                WRITE (IOUT,'(1X,I3,2X,ES13.6)') N, ZONEFLUX(N) !UZO(1,N) !,UNSFLUX(K,KK) to check
            END DO
        END IF
    END IF
    RETURN
    END SUBROUTINE
    !
    ! SUBROUTINE GWF2UNSF1INIT(IN,INEVT,INRCH)
    !     ******************************************************************
    !     INIT... subroutine, called by RP
    !      Sets NIZ,KZON,IZ,MAXIT,TOLTH,TOLH
    !     ******************************************************************
    !  USE GLOBAL, ONLY:IOUT,NROW,NCOL
    !  USE GLOBAL, ONLY:LSTCHK
    !  USE PARAMMODULE, ONLY: IZON,NZONAR,ZONNAM
    !     CHECK: NZONAR replaced NZN in MF2005
    !  USE GWFUNSFMODULE, ONLY: IZ,KZON,NIZ
    !  USE GWFUNSFMODULE, ONLY: MAXIT,TOLTH,TOLH
    !  IMPLICIT NONE
    !      Variables
    !  INTEGER,INTENT(IN) :: IN,INEVT,INRCH
    !  CHARACTER(LEN=*),PARAMETER :: WARN1 =&
    !       & '----- WARNING - EVAPORATION PACKAGE EXISTS! ----'
    !  CHARACTER(LEN=*),PARAMETER :: WARN2 =&
    !       & '----- WARNING - RECHARGE PACKAGE EXISTS! -----'
    !  CHARACTER(LEN=*),PARAMETER :: FMT1 ="(' ABS ',A,' TOLERANCE IN&
    !       & UNSAT FLOW REGION = ',F8.5)"
    !  CHARACTER(LEN=100):: LINE
    !  CHARACTER(LEN=10) :: ZONARR,ZONAME,ZONARR0
    !  INTEGER :: I,J,K,NTOT
    !  LOGICAL:: LZFOUND
    !      Body
    !      Warn the user if either the evaporation package or the recharge
    !      package is active
    !  IF(LSTCHK(2))THEN
    !     CHECK: Consider using UPKGS (see utl7ext.f90)
    !    IF (INEVT.GT.0) WRITE(IOUT,'(//A/)') WARN1
    !    IF (INRCH.GT.0) WRITE(IOUT,'(//A/)') WARN2
    !  END IF
    !      Read Zone Array
    !  IF(LSTCHK(3)) WRITE(IOUT,'(//A/)') 'UNSF -- READ ZONE ARRAY'
    !  READ(IN,'(A)') ZONARR
    !  CALL UPCASE(ZONARR)
    ! modified by Adam Szymkiewicz, get rid of trailing tabs
    ! when the input file is imported from Excel
    !  DO I=1, 10
    !	ZONARR0(I:I)=CHAR(32)
    !  END DO
    !  J=1
    !  DO I=1, LEN(ZONARR)
    !	IF ((ICHAR(ZONARR(I:I)).GE.32).AND.(ICHAR(ZONARR(I:I)).LE.126)) THEN
    !	  ZONARR0(J:J)=ZONARR(I:I)
    !	  J=J+1
    !	END IF
    !  END DO
    !  ZONARR=ZONARR0
    !
    !  IF (ZONARR.EQ.'ALL') THEN
    !    IF (LSTCHK(3)) WRITE(IOUT,*)' ALL CELLS ARE ASSOCIATED WITH ONE PROFILE'
    !	NIZ=0
    !  ELSE IF(ZONARR.EQ.'EACH') THEN
    !    NIZ=2
    !    NTOT=NROW*NCOL
    !    DO I=1,NTOT
    !      IZ(I)=I
    !    END DO
    !  ELSE
    !    NIZ=1
    !       Read zone names to find ZONEARR and copy the
    !       correct values into KZON
    !	LZFOUND=.FALSE.
    !    DO K=1,NZONAR
    !      ZONAME=ZONNAM(K)
    !      CALL UPCASE(ZONAME)
    !      IF (TRIM(ZONAME).EQ.TRIM(ZONARR)) THEN
    !		LZFOUND=.TRUE.
    !        IF(LSTCHK(3)) WRITE(IOUT,*) 'ZONE ARRAY NAME = ',ZONARR
    !		DO I=1,NROW
    !          DO J=1,NCOL
    !            KZON(J,I)=IZON(J,I,K)
    !          END DO
    !        END DO
    !	  END IF
    !	END DO
    !	IF (.NOT.(LZFOUND)) STOP 'ZONE ARRAY NOT DEFINED'
    !    CALL UHRCOM(IN,IOUT,LINE)
    !    BACKSPACE IN
    !  END IF
    !     CHECK: Consider using URWORD
    !  CALL UHRCOM(IN,IOUT,LINE)
    !  BACKSPACE IN
    !  READ (IN,*) MAXIT,TOLTH,TOLH
    !  IF (LSTCHK(3)) THEN
    !    WRITE(IOUT,*)
    !    WRITE(IOUT,*) 'MAXIMUM ITERATION = ',MAXIT
    !    WRITE(IOUT,FMT1) 'WATER CONTENT',TOLTH
    !    WRITE(IOUT,FMT1) 'PRESSURE HEAD',TOLH
    !  END IF
    !  RETURN
    !END SUBROUTINE


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
        HSAT(M) = FH(IMODEL,1.0d0,PAR(:,M))
        ! modified by AdamS, compilation errors in calls to FH, FK,FC1,FQ,FS
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
    !  WRITE(*,'(5es10.3)') HTAB(J), THETAB(J,1), CONTAB (J,1), CAPTAB(J,1)
    !END DO
    !PRINT *, FQ(IMODEL,-7.99D0,PAR(:,2))
    !DO M=1,NMAT
    !   PRINT *, HSAT(M)
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
        &KBOLD,KODBOT,NIZ,IZ
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
    !	   PRINT *, PID, IZ(PID)
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    ! AdamS : variable PERIMP, which I did not need, was set to 0
    READ(IN,*) SINKF(PID),WLAYER(PID),LINITW(PID),NOMAT(PID) !, &
        !  &PERIMP(PID)
        PERIMP(PID) = 0
    ATMBC(PID) =0
    LWAT(PID)  =0
    SHORTO(PID)=-1
    INFTOP(PID)=1
    INFBOT(PID)=1
    KODBOT(PID)=3
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
        WRITE(IOUT,*)'*** ERROR: CHECK THE NUMBER OF NODES'
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
        IF(J-N.LT.0) THEN
            WRITE(IOUT,*)'ERROR IN NODINF AT NODE =', N
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
            !IF(LSTCHK(1))THEN
            WRITE(IOUT,*)'*** ERROR: INITIAL WATER CONTENT CONDITION ',&
                &'IS LOWER THAN QR !'
            !END IF
            CALL USTOP(' ')
        END IF
        ! modified by AdamS, compilation error in call to FH
        HHNEW(I,PID)=FH(IMODEL,QE,PAR(:,M))
        HHOLD(I,PID)=HHNEW(I,PID)
        HTEMP(I,PID)=HHNEW(I,PID)
    END DO
    HBOT(PID)=HHNEW(1,PID)
    HTOP(PID)=HHNEW(NUMNP(PID),PID)
    RETURN
    END SUBROUTINE


    SUBROUTINE GWF2UNSF1SINKIN(PID,IN)
    ! ******************************************************************
    ! READ SINK INFORMATION
    ! SETS: P0,P2H,P2L,P3,R2H,R2L,POPTM
    ! ******************************************************************
    USE GLOBAL, ONLY: IOUT,LSTCHK
    USE GWFUNSFMODULE, ONLY: P0,P2H,P2L,P3,R2H,R2L,POPTM,nmat
    IMPLICIT NONE
    Character *500 Line
    ! VARIABLES
    INTEGER,INTENT(IN) :: PID,IN
    integer::ii
    ! BODY
    IF (LSTCHK(3)) THEN
        WRITE(IOUT,'(//A/)') 'READING SINK INFORMATION'
    END IF
    READ(IN,*) P0(PID),P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID)
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    READ(IN,*)(POPTM(PID,ii),ii=1,nmat)
    P0(PID) =-ABS(P0(PID))
    P2L(PID)=-ABS(P2L(PID))
    P2H(PID)=-ABS(P2H(PID))
    P3(PID) =-ABS(P3(PID))
    RETURN
    END SUBROUTINE

    SUBROUTINE GWF2UNSF1READBC(PID,IN)
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
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: PID,IN
    CHARACTER(LEN=100) :: LINE
    INTEGER :: I
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
            &HCA(PID,I)
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
            EXIT
        END IF
    END DO
    RETURN
    END SUBROUTINE


    !SUBROUTINE GWF2UNSF1SUBREG(PID)
    !!     CHECK: This function barely serves any purpose, as most of it is
    !!            commented in the original code. The main purpose of this
    !!            subroutine appears to have been printing (before that was
    !!            removed). The non-commented parts have been left, despite
    !!            serving no purpose.
    !!            See line 1879 of the original code for more information
    !!     ******************************************************************
    !!      Sets: WATIN, if LWAT(PID)>=0 and PLEVEL(PID)==0
    !!     ******************************************************************
    !USE GWFUNSFMODULE, ONLY: NUMNP,NOLAY,HHNEW,THETA,THOLD,&
    !    &X,MATNUM,LAYNUM,DT,CON,PLEVEL,&
    !    &WATIN,LWAT
    !IMPLICIT NONE
    !!      Variables
    !INTEGER,INTENT(IN) :: PID
    !INTEGER :: I,J,N,MI,MJ,LAY
    !DOUBLE PRECISION    :: ATOT,VOLUME,CHANGE,HTOT,DELTW,HE,VNEWI,VOLDI,DX,DX1,&
    !    &V1,DXN,VN,DH1,DHN
    !!     CHECK: Consider removing the hard limit on array size
    !DOUBLE PRECISION,DIMENSION(10) :: AREA,SUBVOL,SUBCHA,HMEAN
    !!      Body
    !N=NUMNP(PID)
    !ATOT=0.
    !HE=0.
    !IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
    !    VOLUME=0.
    !    CHANGE=0.
    !    HTOT  =0.
    !    DELTW =0.
    !END IF
    !!     CHECK: NOLAY(PID) has been set to 1 in GWF2UNSF1BASINF, consider
    !!            removing the loop if there are no plans to restore this
    !!            functionality
    !DO LAY=1,NOLAY(PID)
    !    AREA(LAY)=0.
    !    IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
    !        SUBVOL(LAY)=0.
    !        SUBCHA(LAY)=0.
    !        HMEAN (LAY)=0.
    !    END IF
    !END DO
    !DO I=N-1,1,-1
    !    J=I+1
    !    MI=MATNUM(I,PID)
    !    MJ=MATNUM(J,PID)
    !    LAY=LAYNUM(J,PID)
    !    DX=X(J,PID)-X(I,PID)
    !    AREA(LAY)=AREA(LAY)+DX
    !    ATOT=ATOT+DX
    !    IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
    !        HE=(HHNEW(I,PID)+HHNEW(J,PID))/2.
    !        VNEWI=DX*(THETA(I,PID)+THETA(J,PID))/2.
    !        VOLDI=DX*(THOLD(I,PID)+THOLD(J,PID))/2
    !        VOLUME=VOLUME+VNEWI
    !        CHANGE=CHANGE+(VNEWI-VOLDI)/DT(PID)
    !        SUBCHA(LAY)=SUBCHA(LAY)+(VNEWI-VOLDI)/DT(PID)
    !        SUBVOL(LAY)=SUBVOL(LAY)+VNEWI
    !        HTOT=HTOT+HE*DX
    !        HMEAN(LAY)=HMEAN(LAY)+HE*DX
    !    END IF
    !    IF (LWAT(PID).GE.0) THEN
    !        IF (PLEVEL(PID).EQ.0) THEN
    !            WATIN(I,PID)=INT(VNEWI)
    !        ELSE
    !            DELTW=DELTW+ABS(WATIN(I,PID)-VNEWI)
    !        END IF
    !    END IF
    !END DO
    !DO LAY=1,NOLAY(PID)
    !    IF (AREA(LAY).GT.0) THEN
    !        IF(LWAT(PID).GE.0 .OR. PLEVEL(PID).EQ.0) THEN
    !            HMEAN(LAY)=HMEAN(LAY)/AREA(LAY)
    !        END IF
    !    END IF
    !END DO
    !!     CHECK: Verify if this is correct as operator priority is easy
    !!            to misunderstand without parentheses.
    !!            (A .AND. B .OR. C) <=> ((A.AND.B) .OR. C)
    !!            See line 1940 of the original code for more information
    !IF(LWAT(PID).GE.0 .AND. ATOT.GT.0. .OR. PLEVEL(PID).EQ.0) THEN
    !    HTOT=HTOT/ATOT
    !END IF
    !DX1=X(2,PID)-X(1,PID)
    !DH1=HHNEW(2,PID)-HHNEW(1,PID)
    !V1=-(CON(1,PID)+CON(2,PID))/2.*(DH1/DX1+1.)
    !DXN=X(N,PID)-X(N-1,PID)
    !DHN=HHNEW(N,PID)-HHNEW(N-1,PID)
    !VN=-(CON(N,PID)+CON(N-1,PID))/2.*(DHN/DXN+1.)
    !!     CHECK: The original code has 54 lines of commented code here
    !!            See line 1946 of the original code for more information
    !RETURN
    !END SUBROUTINE

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
    ! AdamS: seems that without the line below the coupling works better
    !  IF(KPER.GT.0) T(PID) = T(PID)+DT(PID)
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
        &CON,TOLD,HCRITA,SINK,PAR,X,THETA
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN)::PID
    DOUBLE PRECISION,PARAMETER::RMAX=1.E+10
    DOUBLE PRECISION,PARAMETER::RMIN=1.D-100
    LOGICAL :: CONVGF,ITCRIT,RESTART,WTFOUND
    INTEGER :: I,M
    DOUBLE PRECISION :: EPSH,EPSTH,TH,DX
    DOUBLE PRECISION :: PB,RB,SB,PT,RT,ST,QQ1
    DOUBLE PRECISION,DIMENSION(NUMNP(PID)) :: P,R,S
    !      Body
    !      CONTINUE ! to check -  Adam Szymkiewicz

    !DX=X(2,PID)-X(1,PID)
    !QQ1=-(CON(1,PID)+CON(2,PID))/2.*&
    !          &((HHNEW(2,PID)-HHNEW(1,PID))/DX+1.)
    !PRINT *,HHNEW(1,PID),HHNEW(2,PID)

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
                ! modified by AdamS, missing index (PID) inserted
                T(PID)=TOLD(PID)+DT(PID)
                RESTART=.TRUE.
                ! modified by AdamS, missing index (PID) inserted
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
    THETA=>THNEW
    CALL GWF2UNSF1SETMAT(PID)
    THNEW=THETA! modified by AdamS, check if necessary
    IF(WLAYER(PID).GE.0) THEN
        IF (HHNEW(NUMNP(PID),PID).GT. HCRITS(PID)) THEN
            KODTOP(PID)=4
            HTOP(PID)=HCRITS(PID)
        END IF
    END IF

    ! modified by AdamS, the exchange flux is now computed above the water table
    ! in order to avoid artificial oscillations
    !  WTFOUND=.FALSE.
    !  DO I=1,NUMNP(PID)-3
    !    IF (HHNEW(I,PID).GE.0.0D0 .AND. HHNEW(I+1,PID).LE.0.0D0) THEN
    !        WTFOUND=.TRUE.
    !		DX=X(I+3,PID)-X(I+2,PID)
    !      VBOT(PID)=-(CON(I+3,PID)+CON(I+2,PID))/2.*&
    !            &((HHNEW(I+3,PID)-HHNEW(I+2,PID))/DX+1.)
    !	  EXIT
    !	END IF
    !  END DO
    !  IF (WTFOUND==.FALSE.) THEN
    DX=X(2,PID)-X(1,PID)
    VBOT(PID)=-(CON(1,PID)+CON(2,PID))/2.*&
        &((HHNEW(2,PID)-HHNEW(1,PID))/DX+1.)

    !  ENDIF

    DO I=1,NUMNP(PID)-1
        DX=X(I+1,PID)-X(I,PID)
        QQ1=-(CON(I+1,PID)+CON(I,PID))/2.*&
            &((HHNEW(I+1,PID)-HHNEW(I,PID))/DX+1.)

    ENDDO
    !	READ(*,*)
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
            ! modified by AdamS, compilation error in call to FK
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
            ! modified by AdamS, compilation error in call to FC1,FQ
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

    IF(TOLD(PID).LT.TOTIMOLD) THEN
        SUMVBOT(PID)=0.0
    ELSE
        SUMVBOT(PID)=SUMVBOT(PID)+VBOT(PID)*DT(PID)

    END IF
    RETURN
    END SUBROUTINE

    SUBROUTINE GWF2UNSF1TLSUM
    !     ******************************************************************
    ! AdamS: I added this subroutine to write bottom pressure heads and exchange fluxes
    ! for each Hydrus profile
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

    SUBROUTINE GWF2UNSF1TLINF(PID,THN,THO )
    !     ******************************************************************
    !     CHECK:
    !     ******************************************************************
    USE GWFUNSFMODULE, ONLY:NUMNP,CON,X,T,DT,TLEVEL,&
        &RTOP,RROOT,VROOT,HHNEW,CUMQ,&
        &LWAT,WCUMT,WCUMA,THNEW,&
        &THOLD,SINK,RSOIL,PREC,THN,THO,NPUNSF,NUMNPD,&
        &LPRINT,LPTIME,ITLOU,cumch,cvbot,cvtop,cvchr,cvch1,cvch0&
        &,ccuma,ccumt,cRunOff,ns,lchem,conc,ctop,cvch1,cvchim,cgwl,VOLUME,hh
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: PID
    INTEGER :: I,J,N,M,js,ibreak
    DOUBLE PRECISION :: DX,DXN,DX1,VTOP,VBOT,RINFIL,REVAP ,dgwl,vrunoff,&
        & ThN(Numnpd,npunsf) ,ThO(Numnpd,npunsf)
    !      Body

    WCUMT(PID)=0.
    WCUMA(PID)=0.

    N=NUMNP(PID)
    M=N-1
    DXN=X(N,PID)-X(M,PID)
    VTOP=-(CON(N,PID)+CON(M,PID))/2.*((HHNEW(N,PID)-HHNEW(M,PID))&
        &/DXN+1.)-&
        &(THN(N,PID)-THO(N,PID))*DXN/2./DT(PID)-SINK(N,PID)*DXN/2.
    DX1=X(2,PID)-X(1,PID)
    VBOT=-(CON(2,PID)+CON(1,PID))/2.*((HHNEW(2,PID)-HHNEW(1,PID))&
        &/DXN+1.)-&
        &(THN (1,PID)-THO (1,PID))*DXN/2./DT(PID)-SINK(1,PID)*DXN/2.

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

    ! AdamS : the loop should start from N-1, otherwise X array range is exceeded
    DO I=N-1,1,-1
        J=I+1
        DX=X(J,PID)-X(I,PID)
        VOLUME=VOLUME+DX*(THN (I,PID)+THN (J,PID))/2.
    END DO


    if(lChem) then
        do 11 jS=1,NS
            CumCh(1,jS)=CumCh(1,jS)-cvTop(jS,pid)*dt(pid)
            CumCh(2,jS)=CumCh(2,jS)+cvBot(jS,pid)*dt(pid)
            CumCh(3,jS)=CumCh(3,jS)+cvCh0(jS,pid)*dt(pid)
            CumCh(4,jS)=CumCh(4,jS)+cvCh1(jS,pid)*dt(pid)
            CumCh(5,jS)=CumCh(5,jS)+cvChR(jS,pid)*dt(pid)
            CumCh(6,jS)=CumCh(6,jS)+cvChIm(jS,pid)*dt(pid)
            cCumT(jS)=cCumT(jS)+(cvTop(jS,pid)-cvBot(jS,pid)-cvCh0(jS,pid)-cvCh1(jS,pid)+&
                &                        cvChR(jS,pid))*dt(pid)
            cCumA(jS)=cCumA(jS)+(abs(cvBot(jS,pid))+abs(cvTop(jS,pid))+&
                &                         abs(cvCh0(jS,pid))+abs(cvCh1(jS,pid))+&
                &                         abs(cvChR(jS,pid)))*dt(pid)


            cRunOff(jS)=vRunOff*cTop(pid,jS)
            CumCh(10,jS)=CumCh(10,jS)+cRunOff(jS)*dt(pid)

            !      Average GWL concentration
            cGWL(jS)=0.
            dGWL=0
            iBreak=0
            do 14 i=1,N-1
                j=i+1
                dx=x(j,pid)-x(i,pid)
                if(hhNew(j,pid).gt.0..and.iBreak.eq.0) then
                    cGWL(jS)=cGWL(jS)+(Conc(jS,i,PID)+Conc(jS,j,PID))/2.*dx
                    dGWL=dGWL+dx
                else
                    iBreak=1
                end if
14          continue
            if(dGWL.gt.0.) cGWL(jS)=cGWL(jS)/dGWL

11      continue

    end if


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
    IF(ITER(PID).LE.ITMIN .AND. (TFIX-T(PID)).GE.DMUL*DTOPT(PID))THEN
        DTOPT(PID)=DMIN1(DTMX,DMUL*DTOPT(PID))
    END IF
    IF(ITER(PID).GE.ITMAX) THEN
        DTOPT(PID)=DMAX1(DTMIN(PID),DMUL2*DTOPT(PID))
    END IF
    DT(PID)=DMIN1(DTOPT(PID),TFIX-T(PID))
    DT(PID)=DMIN1((TFIX-T(PID))/DNINT((TFIX-T(PID))/DT(PID)),DTMX)
    IF(DT(PID).LE.0)THEN
        !IF(LSTCHK(1)) THEN
        WRITE(IOUT,"(/1X,'DT<0 ! wrong time discr MODFLOW/HYDRUS')")
        !END IF
        CALL USTOP(' ')
    END IF
    IF((TFIX-T(PID)).NE.DT(PID).AND.DT(PID).GT.(TFIX-T(PID)/2.))THEN
        DT(PID)=(TFIX-T(PID))/2.
    END IF
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
            ! IF(HHNEW(I,PID).LT.0. .AND. HHOLD(I,PID).LT.0.)THEN
            !  HTEMP(I,PID)=HHNEW(I,PID)+(HHNEW(I,PID)-HHOLD(I,PID))*&
            !                                             &DT(PID)/DTOLD(PID)
            ! ELSE
            HTEMP(I,PID)=HHNEW(I,PID)
            !END IF
            HHOLD(I,PID)=HHNEW(I,PID)
            !HHNEW(I,PID)=HTEMP(I,PID)
            THOLD(I,PID)=THNEW(I,PID)
        END IF
    END DO
    RETURN
    END SUBROUTINE

    !SUBROUTINE GWF2UNSF1NODOUT1(TVAL,LHEAD)
    !!     ******************************************************************
    !!     Print node information for the current time step into a specific
    !!     unit.
    !!     CHECK: NODOUT exists in the original code, but the line where it
    !!            is called has been commented
    !!     CHECK: This subroutine's goal is clearly to print data into a file
    !!            However, modifies HHOLD, which should not happen, even if
    !!            the values stored in HHOLD do not matter at the moment
    !!     ******************************************************************
    !USE GWFUNSFMODULE, ONLY:NUMNP,HHNEW,THOLD,X,XSURF,NUMNPD,NPUNSF,&
    !    &LOUTF,IPROUH,IPROUW !,IPROU
    !IMPLICIT NONE
    !!      Variables
    !LOGICAL,INTENT(IN) :: LHEAD
    !DOUBLE PRECISION, INTENT(IN) :: TVAL
    !CHARACTER(LEN=*),PARAMETER:: FMT1="(//'PROFILE AT TIME:',F14.3/)"
    !CHARACTER(LEN=*),PARAMETER:: FMT2="(//'PROFILE AT TIME:',E15.8/)"
    !CHARACTER(LEN=*),PARAMETER:: FMT3=&
    !    &"('Node',100(I3, '_Depth[L]',4X, I3, '_Head[L]',4X))"
    !CHARACTER(LEN=*),PARAMETER:: FMT4=&
    !    &"('Node',100(I3, '_Depth[L]',2X, I3, '_WC[-]',3X))"
    !CHARACTER(LEN=*),PARAMETER:: FMT5=&
    !    &"(1X,I3,4X,100(F8.3,5X,F11.3,6X))"
    !CHARACTER(LEN=*),PARAMETER:: FMT6=&
    !    &"(I3,3X,100(F8.3,6X,F6.4,6X))"
    !INTEGER :: I,J,K,PID
    !DOUBLE PRECISION,DIMENSION(NUMNPD,NPUNSF)::DEPTHS,VALUES
    !!      Body
    !!      Prepare values
    !
    !! modified by AdamS: there are 2 files for soil profiles
    !! one for pressure heads, the other one for water content
    !
    !DO PID=1,NPUNSF
    !
    !    DEPTHS(1:NUMNP(PID),PID)=X(1:NUMNP(PID),PID)-XSURF(PID)
    !    IF (NUMNP(PID).LT.NUMNPD) THEN
    !        DEPTHS((NUMNP(PID)+1):NUMNPD,PID) = 999999
    !        HHNEW ((NUMNP(PID)+1):NUMNPD,PID) = 999999
    !        THOLD ((NUMNP(PID)+1):NUMNPD,PID) = 999999
    !    END IF
    !END DO
    !!      Print Values
    !WRITE(IPROUH,FMT1) TVAL
    !WRITE(IPROUH,FMT3)(I,I,I=1,NPUNSF)
    !DO J=NUMNPD,1,-1
    !    WRITE(IPROUH,FMT5) NUMNPD-J+1,(DEPTHS(J,K),HHNEW(J,K),&
    !        &K=1,NPUNSF)
    !END DO
    !
    !WRITE(IPROUW,FMT1) TVAL
    !WRITE(IPROUW,FMT4)(I,I,I=1,NPUNSF)
    !DO J=NUMNPD,1,-1
    !    WRITE(IPROUW,FMT6) NUMNPD-J+1,(DEPTHS(J,K),THOLD(J,K),&
    !        &K=1,NPUNSF)
    !END DO
    !
    !!  IF(LOUTF)THEN
    !!        IF(TVAL.LT.99999999.) THEN
    !!         WRITE(IPROU,FMT1) TVAL
    !!        ELSE
    !!         WRITE(IPROU,FMT2) TVAL
    !!        END IF
    !!        IF(LHEAD) THEN
    !!         WRITE(IPROU,FMT3)(I,I,I=1,NPUNSF)
    !!        ELSE
    !!         WRITE(IPROU,FMT4)(I,I,I=1,NPUNSF)
    !!        END IF
    !!        IF (LHEAD) THEN
    !!         VALUES=HHNEW
    !!         DO J=NUMNPD,1,-1
    !!          WRITE(IPROU,FMT5)NUMNPD-J+1,(DEPTHS(J,K),VALUES(J,K),&
    !!                         &K=1,NPUNSF)
    !!         END DO
    !!        ELSE
    !!         VALUES=THOLD
    !!         DO J=NUMNPD,1,-1
    !!          WRITE(IPROU,FMT6)NUMNPD-J+1,(DEPTHS(J,K),VALUES(J,K),&
    !!                         &K=1,NPUNSF)
    !!         END DO
    !!        END IF
    !!       END IF
    !RETURN
    !END SUBROUTINE

    SUBROUTINE GWF2UNSF1FOPEN()
    !     ******************************************************************
    !     Opens Hydrus-specific output files
    ! AdamS : changed the number of output files
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
    ! AdamS : changed the number of output files
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
    !----------------------------------------------------------------
    !Solute variables and parmeters: Modified by Sahila
    !    GWFUNSFDAT(IGRID)%NS=>NS

    !--------------------------------------------------------------
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
    !----------------------------------------------------------------
    !Solute variables and parmeters: Modifief by Sahila
    !    NS=>GWFUNSFDAT(IGRID)%NS
    !----------------------------------------------------------------
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

    !      ------------------------------------------------
    !      DEALLOCATE: SOLUTE TRANSPORT
    !    DEALLOCATE(GWFUNSFDAT(IGRID)%NS)

    !      ------------------------------------------------

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
    !      ------------------------------------------------
    !      ALLOCATE: SOLUTE TRANSPORT Modified by Sahila
    ALLOCATE(NS)
    !      ------------------------------------------------
    RETURN
    END SUBROUTINE

    double precision function RTime(iMonth,iDay,iHours,iMins,iSecs,&
        &                                i100th)

    integer*2 iMonth,iDay,iHours,iMins,iSecs,i100th

    if(iMonth.eq.1.or.iMonth.eq.3.or.iMonth.eq.5.or.iMonth.eq.7.or.&
        &   iMonth.eq.8.or.iMonth.eq.10.or.iMonth.eq.12) then
    NoDay=31
    else if(iMonth.eq.4.or.iMonth.eq.6.or.iMonth.eq.9.or.iMonth.eq.11) &
        &                                                then
    NoDay=30
    else if(iMonth.eq.2) then
        NoDay=28
    end if
    nMonth=NoDay*24.*60.*60.
    RTime=nMonth+iDay*24.*60.*60.+iHours*60.*60.+iMins*60.+iSecs+&
        &      i100th/100.

    return
    end

    subroutine runtime(time1,time2,kper,kstp,pid)
    double precision:: totaltime, time2,time1
    double precision, save :: sum
    INTEGER,INTENT(IN) ::kSTp,pid,kper
    Totaltime=time2-time1
    sum =sum +totaltime

    end subroutine

    SUBROUTINE HYD_README()
    !   logical: Lchem


    write (98,100)



    if (lchem) write (98,*)  'MODEL TYPE:WATER FLOW AND CHEMICAL TRANSPORT IS CONSIDERED'

    if (.not.lchem) write (98,*)  'MODEL TYPE:ONLY WATER FLOW IS CONSIDERED'

    write (98,110)
    write (98,120)

    IF (lchem) write (98,130)


100 Format('              HYDRUS PACKAGE FOR MODFLOW'/&
        &	   '======================================================'/&
        &       '                    version 1.2        '/&
        &       '                Updated:July 08,2008   '/&
        &       '------------------------------------------------------'//&
        &       '                       Authors		  '/&
        &	   '             Navin Twarakavi & J.Simunek'/&
        &	   '     visit: www.pc-progress.cz for more information'/&
        &   	   '------------------------------------------------------'//)


110 Format (//'The readme file gives an explanation of the various &
        & output files generated in this run.'//)

120 Format (/ &
        &  'HYDRUS_TInfo.OUT:This file has information on the ground water&
        &  table depth and bottom fluxes at  various times (time steps/ and&
        &  stress periods).In the case of chemical transport, the &
        &  concentrations at the water table as well as the flux at the &
        &  water table is also outputted' &
        &//'HYDRUS_profile.OUT:This file is generated if PROPR>0. This file &
        & has water flow and solute transport information for the hydrus&
        & profiles at selected print times.'// )

130 Format ('HYD_Conc_GWT.OUT: This file gives the average &
        & concentrations at the water table for each stress period. One may &
        & use these grids as the source function in MT3DMS for performing &
        & saturated contamination transport simulations'//&
        & 'HYD_ConFlux_GWT.OUT: This file gives the average concentration &
        & fluxes for the solutes at the water table for each stress period. &
        &One may  use these grids as the source function in MT3DMS for &
        & saturated contamination transport simulations'//)
    end  subroutine

