! CODE DEVELOPED BY SCOTT E BOYCE
!                   CONTACT <seboyce@usgs.gov> or <Boyce@engineer.com>
!
!
! MODULE XY_GRID_COORDINATE_INTERFACE
!  MODULE CREATES A DATA TYPE THAT HOLDES SPACIAL X, Y CARTESIAN COORDINATE INFORMATION.
!  IT IS MODEL/GRID INDEPENDENT AND JUST REQUIRES AN ORIGIN X,Y POINT, NROW, NCOL, AND THEIR LENGTHS/WIDTHS
!  VERSION 1.1 [6/01/2014] NEW OPTIONS FOR POINT OF ORIGIN, NEW DEFAULT POINT OF ORIGIN LOCATION [LOWER LEFT CORNER IS (0,0) TO COFORM WITH A NORMAL AXIS],
!                          ADDED NEW ARRAYS THAT STORE THE MODEL CELL CORNERS X,Y LOCATION,
!                          REMOVED DEPENDENCE ON BAS MODULE GLOBAL
!  VERSION 1.0 [8/14/2013] ORIGINAL VERSION THAT CALCULATES THE X, Y CARTESIAN COORDINATE OF THE AREAL MODEL CELL CENTERS (ROW/COL)

      MODULE XY_GRID_COORDINATE_INTERFACE
      IMPLICIT NONE
      PRIVATE
      PUBLIC:: XY_GRID_COODINATES
      !
      TYPE XY_GRID_COODINATES
         INTEGER::NROW,NCOL                                             ! CURRENT MODEL GRID # OF ROWS AND COLUMNS (COPIES OF WHAT IS STORED IN GLOBAL MODULE, THIS PREVENTS ANY LINK TO GLOBAL)
         DOUBLE PRECISION:: ROTDEG, ROTRAD, MAXDEL                      ! ROTDEG=GRID ROTATION IN DEGREES, ROTRAD=IN RADIANS, MAXDEL=LARGEST CELL WIDTH OR LENGTH MAX([DELR(:),DELC(:)]); MAXDEL IS USED TO DETRMINE A MINIMUM RADIUS AROUND A SEARCH POINT FOR LOCATING A REQUESTED XY COORDINATE.
         DOUBLE PRECISION,DIMENSION(:,:)  ,ALLOCATABLE::XCENT,YCENT     ! XCENT(NCOL,NROW)=CELL CENTER X COORDINATE,      YCENT(NCOL,NROW)=CELL CENTER Y COORDINATE
         DOUBLE PRECISION,DIMENSION(:,:)  ,ALLOCATABLE::XCORN,YCORN     ! XCORN(0:NCOL,0:NROW)= CELL CORNER X COORDINATE, YCORN(0:NCOL,0:NROW)=CELL CORNER X COORDINATE; WHERE XCORN(0,0) IS THE X-COORDINATE OF THE OUTERMOST CORNER [WITH RESPECT TO THE MODEL GRID] OF ROW 1 COL 1 AND XCORN(1,1) IS THE INNER MOST CORNER OF ROW 1 COL 1
         !
         CONTAINS
         !
         PROCEDURE, PASS(XY):: BUILD => BUILD_XYCOORD_FROM_INITIAL_POINT!CALL XY%BUILD(Xin,Yin,ROT,LLCOODRINATE,CORNERCOORD,DELR,DELC) TO BUILD XY COORDINATE SYSTEM
!         PROCEDURE, PASS(XY):: XY2RC => XYCOORD2CELL                    !CALL XY%XY2RC(X,Y,ROW,COL) RECIEVES X AND Y COORDINATE AND RETURNS THE ROW AND COL THAT THE POINT RESIDES IN
         PROCEDURE, PASS(XY):: PRINT => PRINT_XY_GRID                   !CALL XY%PRINT(IOUT) PRINTS XCENT,YCENT,XCORN,YCORN TO FILE WITH UNIT # IOUT.
         !FINAL:: DEALLOCATE_XY_GRID_COODINATES                          !DEALLOCATE(XY) PROVIDES CLEAN DEALLOCATION OF DATA TYPE XY_GRID_COODINATES
         PROCEDURE, PASS(XY):: CLEAR => DEALLOCATE_XY_GRID_COODINATES   ! REPLACES FINAL BECAUSE GCC 4.7.1 DOES NOT SUPPORT FINALIZATION
      END TYPE
      !
      INTERFACE XY_GRID_COODINATES
        MODULE PROCEDURE ALLOCATE_XY_DTYPE                              !CONSTRUCTOR FUNCTION THAT RETURNS AN ALLOCATED XY_GRID_COODINATES POINTER. INVOKED AS: XY => XY_GRID_COODINATES(NROW,NCOL)
      END INTERFACE
      !
      CONTAINS
      !
      !################################################################
      !
      FUNCTION ALLOCATE_XY_DTYPE(NROW,NCOL) RESULT(XY)
      TYPE(XY_GRID_COODINATES),POINTER:: XY
      INTEGER,INTENT(IN)::NROW,NCOL
      !
      XY=>NULL()
      ALLOCATE(XY)
      !
      XY%NROW=NROW
      XY%NCOL=NCOL
      !
      ALLOCATE( XY%XCENT(  NCOL,  NROW), XY%YCENT(  NCOL,  NROW) )
      ALLOCATE( XY%XCORN(0:NCOL,0:NROW), XY%YCORN(0:NCOL,0:NROW) )
      !
      END FUNCTION
      !
      !################################################################
      SUBROUTINE DEALLOCATE_XY_GRID_COODINATES(XY)
      CLASS(XY_GRID_COODINATES):: XY !changed
      !
      IF(ALLOCATED(XY%XCENT))  DEALLOCATE( XY%XCENT, XY%YCENT )
      IF(ALLOCATED(XY%XCORN))  DEALLOCATE( XY%XCORN, XY%YCORN )
      !
      END SUBROUTINE
      !
      !################################################################
      !
      SUBROUTINE BUILD_XYCOORD_FROM_INITIAL_POINT(XY,Xin,Yin,ROT,
     +                          LLCOODRINATE,CORNERCOORD,DELR,DELC)
      !SUBROUTINE TO PROCESS CELL CENTER COORDINATE FOR ENTIRE MODEL
      !NEW FEATURES INCLUDE PROCESSING CELL CORNERS AND CONVERSION OF MAIN ROUTINES TO TYPE-BOUND PROCEDURES FOR INCLUSION IN NEW DATA TYPE XY_GRID_COODINATES
      !CORDINATE SYSTEM IS CARTESIAN WITH A GRID ROTATION THAT STARTS AT 0 DEGREES ALONG THE POSITIVE X-AXIS AND INCREASES COUNTER CLOCWISE (POLAR ANGLE)
      !DEVELOPED BY SCOTT E BOYCE 8/14/2013
      CLASS(XY_GRID_COODINATES):: XY
      DOUBLE PRECISION, INTENT(IN)::Xin,Yin,ROT                         !Xin=X POINT OF REFERENCE,Yin=Y POINT OF REFERENCE,ROT=ROTATION OF MODEL GRID IN DEGREES
      LOGICAL,          INTENT(IN)::LLCOODRINATE,CORNERCOORD            !LLCOODRINATE=TRUE MEANSPOINT OF REFFERENCE IS AT ROW 1 COL NCOL, FALSE MEANS THAT ITS AT ROW 1 COL 1; CORNERCOORD=TRUE MEAND REFERENCE IS ON THE OUTER MOST CORNER, FALSE MEANS THAT IT IS A CELL CENTER
      REAL,DIMENSION(:),INTENT(IN)::DELR                                !DEFINED BY BAS/GLOBAL AND REPRESENTS THE MODEL SPACING. THIS CAN NOT BE CALLED FROM USE "GLOBAL, ONLY" BECAUSE IT WOULD CREATE A CIRCULAR DEPENDENCY.
      REAL,DIMENSION(:),INTENT(IN)::DELC                                !DEFINED BY BAS/GLOBAL AND REPRESENTS THE MODEL SPACING. THIS CAN NOT BE CALLED FROM USE "GLOBAL, ONLY" BECAUSE IT WOULD CREATE A CIRCULAR DEPENDENCY.
      DOUBLE PRECISION:: DISTR,DISTC,SINROT,COSROT,SIN90ROT,COS90ROT
      DOUBLE PRECISION:: X,Y,XCOL1,YCOL1                                !X,Y ARE COORD OF ROW 1, COL 1 CELL and XCOL1, YCOL1 HOLDS CELL CENTER COL 1 WHEN PROCESSING A ROW
      DOUBLE PRECISION,PARAMETER::EPS=1D-5                              !TOLERANCE FOR APPROXIMATING SINE/COSINE TO EXACT VALUE WHEN NEARBY 0, 90, 180, and 270 degrees
      DOUBLE PRECISION,PARAMETER::RAD= 1.74532925199432958D-2           !CONVERSION FACTOR FROM DEGREES TO RADIANS --  3.1415926535897932D0/180D0 = 0.0174532925199432958D0
!      DOUBLE PRECISION,PARAMETER::ANG= 57.295779513082321D0          !CONVERSION FACTOR FROM RADIANS TO DEGREES --  180D0/3.1415926535897932D0
      !DOUBLE PRECISION,DIMENSION(:,:),POINTER,CONTIGUOUS::
!     +                                         XCORD,YCORD,XPCORD,YPCORD
      DOUBLE PRECISION,PARAMETER::MINBOT=1D-20
      DOUBLE PRECISION,DIMENSION(4)::XP,YP,SLOPE,RISE,RUNN,SGN
      DOUBLE PRECISION,POINTER::MX
      DOUBLE PRECISION,DIMENSION(XY%NROW)::DC
      DOUBLE PRECISION,DIMENSION(XY%NCOL)::DR
      INTEGER::I,J,NROW,NCOL
      !
      NROW=XY%NROW
      NCOL=XY%NCOL
      DR=DBLE(DELR)                                                     !SETTING LOCAL DP NAMES TO GLOBAL DELR AND DELC WHICH ARE RENAMED AS DR AND DC, RESPECTIVELY
      DC=DBLE(DELC)
      !
      XY%ROTDEG=ROT
      XY%ROTRAD=ROT*RAD
      !MX    =>XY%MAXDEL
      !XCORD =>XY%XCENT                                              !SET UP LOCAL POINTERS
      !YCORD =>XY%YCENT
      !XPCORD=>XY%XCORN
      !YPCORD=>XY%YCORN
      !
      IF(ABS(ROT)<EPS .OR. ABS(360D0-ROT)<EPS)THEN                      !TO PREVENT NUMERICAL ERRORS SET OBVIOUS SINE/COSINE VALUES
        SINROT  =0D0
        COSROT  =1D0
        SIN90ROT=-1D0
        COS90ROT=0D0
      ELSEIF(ABS(90D0 -ROT)<EPS)THEN
        SINROT  =1D0
        COSROT  =0D0
        SIN90ROT=0D0
        COS90ROT=1D0
      ELSEIF(ABS(180D0-ROT)<EPS)THEN
        SINROT  =0D0
        COSROT  =-1D0
        SIN90ROT=1D0
        COS90ROT=0D0
      ELSEIF(ABS(270D0-ROT)<EPS)THEN
        SINROT  =-1D0
        COSROT  =0D0
        SIN90ROT=0D0
        COS90ROT=-1D0
      ELSE
        SINROT  =SIN(ROT*RAD)
        COSROT  =COS(ROT*RAD)
        SIN90ROT=SIN((ROT-90D0)*RAD)
        COS90ROT=COS((ROT-90D0)*RAD)
      END IF
      !
      IF      (LLCOODRINATE       .AND.       CORNERCOORD) THEN
         !
         DISTC=0.5D0*DC(NROW)                                           !DISTANCE TO CENTER OF ROW NROW AND SIDE OF COL 1
         X=Xin - (DISTC*COS90ROT)                                       !X COORDINATE OF ROW NROW, COL 1
         Y=Yin - (DISTC*SIN90ROT)                                       !Y COORDINATE OF ROW NROW, COL 1
         !
         DISTR=0.5D0*DR(1)                                              !DISTANCE TO CENTER OF ROW NROW, COL 1
         X=X + (DISTR*COSROT)                                           !X COORDINATE OF ROW NROW, COL 1
         Y=Y + (DISTR*SINROT)                                           !Y COORDINATE OF ROW NROW, COL 1
         !
         DISTC=0.5D0*(DC(1)+DC(NROW)) + SUM(DC(2:NROW-1))               !DISTANCE TO ROW 1, COL 1
         X=X - (DISTC*COS90ROT)                                         !X COORDINATE OF ROW 1, COL 1
         Y=Y - (DISTC*SIN90ROT)                                         !Y COORDINATE OF ROW 1, COL 1
         !
      ELSE IF (LLCOODRINATE       .AND. .NOT. CORNERCOORD) THEN
         !
         DISTC=0.5D0*(DC(1)+DC(NROW)) + SUM(DC(2:NROW-1))               !DISTANCE TO ROW 1, COL 1
         X=Xin - (DISTC*COS90ROT)                                       !X COORDINATE OF ROW 1, COL 1
         Y=Yin - (DISTC*SIN90ROT)                                       !Y COORDINATE OF ROW 1, COL 1
         !
      ELSE IF (.NOT. LLCOODRINATE .AND.       CORNERCOORD) THEN
         !
         DISTC=0.5D0*DC(NROW)                                           !DISTANCE TO CENTER OF ROW NROW AND SIDE OF COL 1
         X=Xin + (DISTC*COS90ROT)                                       !X COORDINATE OF ROW NROW, COL 1
         Y=Yin + (DISTC*SIN90ROT)                                       !Y COORDINATE OF ROW NROW, COL 1
         !
         DISTR=0.5D0*DR(1)                                              !DISTANCE TO CENTER OF ROW NROW, COL 1
         X=X + (DISTR*COSROT)                                           !X COORDINATE OF ROW NROW, COL 1
         Y=Y + (DISTR*SINROT)                                           !Y COORDINATE OF ROW NROW, COL 1
         !
      ELSE
        X=Xin
        Y=Yin
      END IF
      !IF(LLCOODRINATE) THEN                                             !USING LOWER LEFT CORNER AS CORNER, SO CALCULATE ROW 1, COL 1 COORD
      !   DISTC=0.5D0*(DC(1)+DC(NROW)) + SUM(DC(2:NROW-1))               !DISTANCE TO ROW 1, COL 1
      !   X=Xin - (DISTC*COS90ROT)                                       !X COORDINATE OF ROW 1, COL 1
      !   Y=Yin - (DISTC*SIN90ROT)                                       !Y COORDINATE OF ROW 1, COL 1
      !ELSE
      !  X=Xin
      !  Y=Yin
      !END IF
      !
      XCOL1=X                                                           !STARTING X COORDINATE IS ROW 1, COL 1
      YCOL1=Y                                                           !STARTING Y COORDINATE IS ROW 1, COL 1
      !
      DISTC=0D0
      DO I=1,NROW
        IF(I>1) THEN                                                    !MOVE DOWN 1 ROW AND PROCESS ALL COLUMNS
          DISTC=DISTC + 0.5D0*(DC(I-1)+DC(I))                           !DISTANCE TO NEXT ROW
          XCOL1=X + (DISTC*COS90ROT)                                    !NEXT ROW'S COLUMN 1's X COORDINATE
          YCOL1=Y + (DISTC*SIN90ROT)                                    !NEXT ROW'S COLUMN 1's Y COORDINATE
        END IF
        DISTR=0D0
        XY%XCENT (1,I)=XCOL1
        XY%YCENT (1,I)=YCOL1
        DO J=2,NCOL
          DISTR=DISTR + 0.5D0*(DR(J-1)+DR(J))
          XY%XCENT (J,I)=XCOL1 + (DISTR*COSROT)
          XY%YCENT (J,I)=YCOL1 + (DISTR*SINROT)
        END DO
      END DO
      !
      !BUILD CORNER POINT ARRAYS
      DO I=1,NROW                                                       !SOLVES FOR UPPER LEFT CORNER FOR EVERY CELL
         DISTC=0.5D0*DC(I)
      DO J=1,NCOL
         DISTR=0.5D0*DR(J)
         XY%XCORN(J-1,I-1)=XY%XCENT(J,I)-(DISTC*COS90ROT)-(DISTR*COSROT)
         XY%YCORN(J-1,I-1)=XY%YCENT(J,I)-(DISTC*SIN90ROT)-(DISTR*SINROT)
        END DO
      END DO
      !
      I=NROW                                                            !SOLVE FOR THE BOTTOM LEFT CORNER OF THE LAST ROW
      DISTC=0.5D0*DC(I)
      DO J=1,NCOL
         DISTR=0.5D0*DR(J)
         XY%XCORN(J-1,I)=XY%XCENT(J,I)+(DISTC*COS90ROT) - (DISTR*COSROT)
         XY%YCORN(J-1,I)=XY%YCENT(J,I)+(DISTC*SIN90ROT) - (DISTR*SINROT)
      END DO
      !
      J=NCOL                                                            !SOLVE FOR THE UPPER RIGHT OF THE LAST COLUMN
      DISTR=0.5D0*DR(J)
      DO I=1,NROW
         DISTC=0.5D0*DC(I)
         XY%XCORN(J,I-1)=XY%XCENT(J,I)-(DISTC*COS90ROT) + (DISTR*COSROT)
         XY%YCORN(J,I-1)=XY%YCENT(J,I)-(DISTC*SIN90ROT) + (DISTR*SINROT)
      END DO
      !
      I=NROW                                                            !SOLVE FOR BOTTOM RIGHT CORNER OF MODEL
      J=NCOL
      XY%XCORN(J,I)=XY%XCENT(J,I) + (DISTC*COS90ROT) + (DISTR*COSROT)   !NOTE DISTR WAS SOLVED BEFORE AT I=NROW
      XY%YCORN(J,I)=XY%YCENT(J,I) + (DISTC*SIN90ROT) + (DISTR*SINROT)
      !
      XY%MAXDEL=MAX(MAXVAL(DR),MAXVAL(DC))
      !
      !NULLIFY( XCORD, YCORD, XPCORD, YPCORD, MX )
      !DEALLOCATE(DR,DC)
      !
      END SUBROUTINE
      !
      !################################################################
      !
      SUBROUTINE PRINT_XY_GRID(XY,IOUT)
            !PRINT OUR GRID IF REQUESTED
      CLASS (XY_GRID_COODINATES):: XY
      INTEGER,INTENT(IN)::IOUT
      INTEGER::I
      CHARACTER(8):: NUM

      WRITE(NUM,'(F8.4)') XY%ROTDEG
      WRITE(IOUT,'(/,/ 4A,/,3A)')'X, Y CELL CENTER COORDINATES ',
     + 'WITH ROTATION = ',TRIM(ADJUSTL(NUM)), '﻿º',
     + 'CORDINATE SYSTEM IS CARTESIAN WITH A GRID ROTATION THAT ',
     + 'STARTS AT 0﻿º ALONG THE POSITIVE X-AXIS AND INCREASES ',
     + 'COUNTER CLOCKWISE (POLAR ANGLE)'
      WRITE(IOUT,'(/A)')'X-COORDINATES'
      DO I=1,XY%NROW
        WRITE(IOUT,'(*(G20.10))') XY%XCENT(:,I)
      END DO
      !
      WRITE(IOUT,'(/A)')'Y-COORDINATES'
      DO I=1,XY%NROW
        WRITE(IOUT,'(*(G20.10))') XY%YCENT(:,I)
      END DO
      !
      !
      WRITE(IOUT,'(/,/ 2A)')'X, Y COORDINATES OF THE MODEL CELL ',
     +'CORNERS.'
      WRITE(IOUT,'(/A)')'X-COORDINATES'
      DO I=0,XY%NROW
        WRITE(IOUT,'(*(G20.10))') XY%XCORN(:,I)
      END DO
      !
      WRITE(IOUT,'(/A)')'Y-COORDINATES'
      DO I=0,XY%NROW
        WRITE(IOUT,'(*(G20.10))') XY%YCORN(:,I)
      END DO
      !
      END SUBROUTINE
      !
      !################################################################
      !
      END MODULE XY_GRID_COORDINATE_INTERFACE
