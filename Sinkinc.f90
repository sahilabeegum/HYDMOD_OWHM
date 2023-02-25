
    subroutine GWF2UNSF1SINKINc(pid,in)
     USE GLOBAL, ONLY: IOUT,LSTCHK
     USE GWFUNSFMODULE, ONLY: ns,P0,P2H,P2L,P3,R2H,R2L,POPTM,CROOTMAX,nmat
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: PID,IN
    integer:: ii
    
	Character *500 Line
    WRITE(IOUT,100) 'READING SINK INFORMATION'
100 FORMAT(//A/)
   
   READ(IN,*) P0(PID),P2H(PID),P2L(PID),P3(PID),R2H(PID),R2L(PID)
    CALL UHRCOM(IN,IOUT,LINE)
    BACKSPACE IN
    READ(IN,*)(POPTM(PID,ii),ii=1,nmat)
    P0(PID) =-ABS(P0(PID))
    P2L(PID)=-ABS(P2L(PID))
    P2H(PID)=-ABS(P2H(PID))
    P3(PID) =-ABS(P3(PID)) 
    
        
    !READ(IN,*) P0,P2H,P2L,P3,R2H,R2L
	!CALL UHRCOM(IN,IOUT,LINE)
	!BACKSPACE IN 
	!Read (IN,*) (POPTM(pid,ii), ii=1,NMAT)
	CALL UHRCOM(IN,IOUT,LINE)
	BACKSPACE IN 
	Read (IN,*) (CRootMax(pid,ii), ii=1,NS)
! DEFAULT: IMOSINK=0, OMEGAC=1, ISOLRED=.FALSE.
! Here we are only considering the simple feddes model with no solute stress or compensated root water uptake
      
    P0(PID) =-ABS(P0(PID))
    P2L(PID)=-ABS(P2L(PID))
    P2H(PID)=-ABS(P2H(PID))
    P3(PID) =-ABS(P3(PID))
    
      RETURN
    end subroutine 
    
    
!    subroutine GWF2UNSF1SINKINc(pid,in)
!     USE GLOBAL, ONLY: IOUT,LSTCHK
!     USE GWFUNSFMODULE, ONLY: ns,P0,P2H,P2L,P3,R2H,R2L,POPTM,CROOTMAX,nmat
!    IMPLICIT NONE
!    INTEGER,INTENT(IN) :: PID,IN
!    integer:: ii
!    
!	Character *500 Line
!    WRITE(IOUT,100) 'READING SINK INFORMATION'
!100   FORMAT(//A/)
!   READ(IN,*) P0,P2H,P2L,P3,R2H,R2L
!	CALL UHRCOM(IN,IOUT,LINE)
!	BACKSPACE IN 
!	Read (IN,*) (POPTM(pid,ii), ii=1,NMAT)
!	CALL UHRCOM(IN,IOUT,LINE)
!	BACKSPACE IN 
!	Read (IN,*) (CRootMax(ii), ii=1,NS)
!! DEFAULT: IMOSINK=0, OMEGAC=1, ISOLRED=.FALSE.
!! Here we are only considering the simple feddes model with no solute stress or compensated root water uptake
!      P0 =-ABS(P0)
!      P2L=-ABS(P2L)
!      P2H=-ABS(P2H)
!      P3 =-ABS(P3)
!      RETURN
!    end subroutine 
!    
    