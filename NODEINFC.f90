      subroutine GWF2UNSF1NODINFC(PID,IN)
      USE GLOBAL, ONLY: IOUT,LSTCHK
        USE GWFUNSFMODULE, ONLY: NUMNPD,IUNSFPR
        USE GWFUNSFMODULE, ONLY: NUMNP,HTOP,HBOT,XSURF
        USE GWFUNSFMODULE, ONLY: X,HHNEW,HHOLD,MATNUM,HTEMP,LAYNUM,BETA
        USE GWFUNSFMODULE, ONLY: CONC,NSD,NS,LEQUIL,LDUALNEQ,SORB2, LCHEM,LEQUIL,LBACT,LDUALNEQ,SORB,TEMPO,TEMPN
        
        IMPLICIT NONE
        ! VARIABLES
        character*30:: Text1,Text2,Text3
        INTEGER,INTENT(IN) :: PID,IN
        lOGICAL :: lTEMP
        CHARACTER(LEN=*),PARAMETER::FMT1="(/'NODAL POINT INFORMATION'//&
            &'NODE X HOLD MATN BETA'/)"
        CHARACTER(LEN=*),PARAMETER::FMT2="(I4,2F11.3,I5,F8.3,E12.4)"
        CHARACTER(LEN=100) :: LINE
        INTEGER :: I,J,M,N,II
        DOUBLE PRECISION :: H,B,X1,DX,SHOLD,nold,IERR,SBETA
        DOUBLE PRECISION:: SConc(5),SSorb(5),C(5),S(5)
        lTEMP=.false.
        
      !Read nodal point information
      WRITE(IOUT,'(//A/)') 'READING NODAL INFORMATION'
        READ(IN,*) NUMNP(PID)
        IF (NUMNP(PID).GT.NUMNPD) THEN
            WRITE(IOUT,*)'*** ERROR: CHECK THE NUMBER OF NODES'
            ! CHECK: REPLACED STOP BY USTOP TO BETTER CONFORM TO MODFLOW
           CALL USTOP(' ')
        END IF
       CALL UHRCOM(IN,IOUT,LINE)
       BACKSPACE IN
      j=NumNP(PID)+1
      NOLD=0
      DO WHILE (J.GT.1)
          j=j-1
          write(iout,*)'j',j
            CALL UHRCOM(IN,IOUT,LINE)
           BACKSPACE IN
          if(lEquil) then
              read(IN,*,err=901) n,x1,h,M,B,(C(ii),ii=1,NS)
          else
              read(IN,*,err=901) n,x1,h,M,B,(C(ii),ii=1,NS),(S(ii),ii=1,NS)
          end if
          WRITE(IOUT,*)'SOLUTE',n,x1,h,M,B,(C(ii),ii=1,NS),(S(ii),ii=1,NS)
          n=NumNP(PID)-n+1
          x(n,PID)=x1
          HhOld(n,PID)=h
          MatNum(n,PID)=M
          LayNum(n,PID)=1
          Beta(n,PID)=B
          !TempO(n,PID)=Te

          do 1 ii=1,NS
              if(lChem) then
                  Conc(ii,n,PID)=C(ii)
                  write(iout,*)'conc:',Conc(ii,n,PID)
                  if(.not.lEquil)       Sorb (ii,n,PID)=S(ii)
                  if(lBact.or.lDualNEq) Sorb2(ii,n,PID)=0.
                  write(iout,*)'sorb_conc:',Sorb (ii,n,PID)
              end if
1         continue
          LAYNUM(N,PID)=1
          if(j-n.LT.0) THEN
              write(iout,*)'ERROR in NodInf at node =', n
              CALL Ustop(' ')
          ELSE IF (J-N.GT.0) THEN
              dx=x(nOld,PID)-x(n,PID)
              ShOld=(HhOld(nOld,PID)-HhOld(n,PID))/dx
              SBeta=(Beta(nOld,PID)-Beta(n,PID))/dx
              !STemp=(TempO(nOld)-TempO(n))/dx   !what is nold? and why we need STemp
              if(lChem) then
                  do 15 ii=1,NS
                      SConc(ii)=(Conc(ii,nOld,PID)-Conc(ii,n,PID))/dx
                      SSorb(ii)=(Sorb(ii,nOld,PID)-Sorb(ii,n,PID))/dx
15                continue
              end if
              do  i=nOld-1,n+1,-1
                  dx=x(nOld,PID)-x(i,PID)
                  HhOld(i,PID)=HhOld(nOld,PID)-ShOld*dx
                  Beta(i,PID)=Beta(nOld,PID)-SBeta*dx
                  !TempO(i)=TempO(nOld)-STemp*dx
                  MatNum(i,PID)=MatNum(i+1,PID)
                  LayNum(i,PID)=LayNum(i+1,PID)
                  if(lChem) then
                      do 16 ii=1,NS
                          Conc(ii,i,PID)=Conc(ii,nOld,PID)-SConc(ii)*dx
                          Sorb(ii,i,PID)=Sorb(ii,nOld,PID)-SSorb(ii)*dx
                          if(lBact.or.lDualNEq)    Sorb2(ii,i,PID)=Sorb2(ii,nOld,PID)-SSorb(ii)*dx
16                    continue
                  end if

              end do
              j=n
          END IF
          J=N
          nOld=n
      END DO
SBETA=BETA(NUMNP(PID),PID) * &
        &(X(NUMNP(PID),PID)-X(NUMNP(PID)-1,PID)) /2.
    DO I=2,NUMNP(PID)-1
        SBETA=SBETA+BETA(I,PID)*(X(I+1,PID)-X(I-1,PID))/2.
    END DO
      !SBeta=0.
!      if(Beta(NumNP,PID).gt.0.) SBeta=Beta(NumNP(PID),PID)*(x(NumNP(PID),PID)-x(NumNP(PID)-1,PID))/2.
!      do 19 i=2,NumNP(PID)-1
!        if(Beta(i,PID).gt.0.) SBeta=SBeta+Beta(i,PID)*(x(i+1,PID)-x(i-1,PID))/2.
!19    continue
      do 20 i=2,NumNP(PID)
        if(SBeta.gt.0.) then
          Beta(i,PID)=Beta(i,PID)/SBeta
        else
          Beta(i,PID)=0.
        end if
20    continue
      xSurf=x(NumNP(PID),PID)

!     Print nodal information
!      write(50,110,err=902)
      do 21 n=NumNP(PID),1,-1
        if(lEquil) then
          write(IOUT,120,err=902) NumNP(PID)-n+1,x(n,PID),HhOld(n,PID),MatNum(n,PID),&
     &                          Beta(n,PID),(Conc(ii,n,PID),ii=1,NS)
        else
          write(IOUT,120,err=902) NumNP(PID)-n+1,x(n,PID),HhOld(n,PID),MatNum(n,PID),&
     &                          Beta(n,PID),(Conc(ii,n,PID),ii=1,NS),&
     &                          (Sorb(ii,n,PID),ii=1,NS)
        end if
        hHNew(n,PID) =hHOld(n,PID)
        hTemp(n,PID)=hHOld(n,PID)
        TempN(n,PID)=TempO(n,PID)
21    continue
!      write(50,'(''end'')',err=902)
      hBot(PID)=HhNew(1,PID)
      hTop(PID)=HhNew(NumNP(PID),PID)
!      write(50,130,err=902) NS
      return
!     Error when reading from an input file
901   ierr=1
      return
!     Error when writing into an output file
902   ierr=2
      return

110   format (/'Nodal point information'//&
     &'node      x         hOld    MatN LayN  Beta      Ah       AK ',&
     &'     ATh     Temp    Conc(1...NS)         Sorb(1...NS)'/)
120   format (i4,2f11.3,2i5,f8.3,3f9.3,f8.2,10e12.4,10e12.4)
130   format (/' Number of species in the chain : ',i3)
140   format (///16x,10(15x,a5,i3,')', 7x))
150   format (///16x,10(15x,a5,i3,')',18x))
160   format (///16x,10(15x,a5,i3,')',29x))
170   format (///16x,10(15x,a5,i3,')',40x))
180   format (///16x,10(15x,a5,i3,')',51x))
190   format (///16x,10(15x,a5,i3,')',62x))
260   format (///16x,10(15x,a5,i3,')',73x))
261   format (///14x,10(15x,a5,i3,')',84x))
262   format (///14x,10(15x,a5,i3,')',95x))
263   format (///14x,10(15x,a5,i3,')',106x))
264   format (///14x,10(15x,a5,i3,')',117x))
265   format (///14x,10(15x,a5,i3,')',128x))
266   format (///14x,10(15x,a5,i3,')',139x))
200   format (/'         time     ',10(a29,     2x))
210   format (/'         time     ',10(a29, a12,2x))
220   format (/'         time     ',10(a29,2a12,2x))
230   format (/'         time     ',10(a29,3a12,2x))
240   format (/'         time     ',10(a29,4a12,2x))
250   format (/'         time     ',10(a29,5a12,2x))
270   format (/'         time     ',10(a29,6a12,2x))
271   format (/'       time     ',10(a29,7a12,2x))
272   format (/'       time     ',10(a29,8a12,2x))
273   format (/'       time     ',10(a29,9a12,2x))
274   format (/'       time     ',10(a29,10a12,2x))
275   format (/'       time     ',10(a29,11a12,2x))
276   format (/'       time     ',10(a29,12a12,2x))
301   format (///16x,10(9x,a5,i3,')', 8x))
302   format (/'         time     ',10(a23,3x))
      end
