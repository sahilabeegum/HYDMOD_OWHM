        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 13 10:56:38 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXPANDILUT__genmod
          INTERFACE 
            SUBROUTINE EXPANDILUT(IWK,ISZ,JLU,ALU)
              INTEGER(KIND=4), INTENT(INOUT) :: IWK
              INTEGER(KIND=4), INTENT(IN) :: ISZ
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(INOUT) :: JLU(:)
              REAL(KIND=8) ,ALLOCATABLE, INTENT(INOUT) :: ALU(:)
            END SUBROUTINE EXPANDILUT
          END INTERFACE 
        END MODULE EXPANDILUT__genmod
