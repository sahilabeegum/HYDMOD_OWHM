!     ******************************************************************
!     Utility Package for HYDRUS
!
!     Author : Gleb Goussarov
!     Code based on gwf1unsf.f by SEO et al., in line
!     with the general style of the MF OWHM project.
!     Note however, that some f77 features have been replaced by
!     f90 features.
!     Subroutines in this package use "IMPLICIT NONE"
!     ******************************************************************
!
      SUBROUTINE UHRCOM(IN,IOUT,LINE)
!     ******************************************************************
!     Utility for Hydrus: Read COMments (and print them)
!     Arguments:
!      IN    : ID of the input unit
!      IOUT  : ID of the output unit
!      LINE  : contains the first non-comment line
!     Note:
!      This subroutine differs from URDCOM in that also allows blank
!      lines and that it does not remove trailing blanks
!     ******************************************************************
       USE GLOBAL,  ONLY:LSTCHK
       IMPLICIT NONE
!      Declarations
       INTEGER,INTENT(IN) :: IN,IOUT
       CHARACTER(LEN=*),INTENT(INOUT) :: LINE
       INTEGER :: CTN,L
!      Body
       CTN = 1
       DO WHILE(CTN.EQ.1)
        READ(IN,'(A)') LINE
        L = LEN(LINE)
        IF(LINE(1:1).EQ.'#' .OR. LINE(1:L).EQ.' ')THEN
         IF(LSTCHK(3)) WRITE(IOUT,'(A)') LINE  !(1:L) modified by Adam Szymkiewicz, 28.05.2016
        ELSE
         CTN=0
        END IF
       END DO
       RETURN
      END SUBROUTINE
!
      INTEGER FUNCTION UHFUAC(IC,IR,NLAY)
!     ******************************************************************
!     Utility for Hydrus: Find Upper Active Cell
!     Arguments:
!      IC    : target column
!      IR    : target row
!      NLAY  : amount of layers
!     ******************************************************************
       USE GLOBAL, ONLY: IBOUND
       IMPLICIT NONE
       INTEGER :: IC,IR,NLAY
!      Body
       UHFUAC=1
       DO UHFUAC=1,NLAY
        IF(IBOUND(IC,IR,UHFUAC).NE.0)RETURN
       END DO
       IF(UHFUAC.GT.NLAY)UHFUAC=1
       RETURN
      END FUNCTION
