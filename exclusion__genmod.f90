        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 16 13:35:29 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXCLUSION__genmod
          INTERFACE 
            SUBROUTINE EXCLUSION(NUMNP,NMAT,NSD,PAR,CHPAR,THNEW,VNEW,   &
     &THOLD,VOLD)
              INTEGER(KIND=4) :: NSD
              INTEGER(KIND=4) :: NMAT
              INTEGER(KIND=4) :: NUMNP
              REAL(KIND=4) :: PAR(11,NMAT)
              REAL(KIND=4) :: CHPAR(NSD*16+4,NMAT)
              REAL(KIND=4) :: THNEW(NUMNP)
              REAL(KIND=4) :: VNEW(NUMNP)
              REAL(KIND=4) :: THOLD(NUMNP)
              REAL(KIND=4) :: VOLD(NUMNP)
            END SUBROUTINE EXCLUSION
          END INTERFACE 
        END MODULE EXCLUSION__genmod
