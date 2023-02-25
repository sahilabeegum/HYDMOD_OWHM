!     ******************************************************************
!     Author : Gleb Goussarov
!     Extension package for UTL7
!
!     ******************************************************************

      INTEGER FUNCTION UPKGS(PKGNAM)
!     ******************************************************************
!     Utility : get PacKaGe State
!     This function retrieves the value stored in the IUNIT cell
!     associated with the CUNIT cell containing PKGNAM.
!     While this is slower than directly giving the IUNIT ID, it avoids
!     passing arguments unnecessarily, when it comes to initialization
!     Arguments:
!      PKGNAM: Name of the package
!     Return Value:
!      (INTEGER) 0 means that the package was either not found or is
!               inactive
!     ******************************************************************
       USE GLOBAL, ONLY: NIUNIT
       USE GLOBAL, ONLY: IUNIT,CUNIT
       IMPLICIT NONE
!      Variables
       CHARACTER(LEN=4),INTENT(IN) :: PKGNAM
       INTEGER :: I
!      Body
       DO I=1,NIUNIT
        IF(CUNIT(I).EQ.PKGNAM)THEN
         UPKGS=IUNIT(I)
         RETURN
        END IF
       END DO
       UPKGS=0
       RETURN
      END FUNCTION
