SUBROUTINE GWF2UNSF1AD(IGRID,KSTP,KPER)
    !     ******************************************************************
    !     calculate unsat flow (ADvance subroutine)
    !     ******************************************************************
    USE GWFBASMODULE,  ONLY: TOTIM,DELT
    USE GWFUNSFMODULE, ONLY: NPUNSF,TMAX, THETA,SINKF,PLEVEL,VBOT,SUMVBOT,ZONEFLUX,PERIMP
    USE GWFUNSFMODULE, ONLY: T,TATM,DT,TOLD,TLEVEL,DTOLD,DTMIN
    USE GWFUNSFMODULE, ONLY: INFTOP,INFBOT,ATMBC,ALEVEL,HCA
    USE GWFUNSFMODULE, ONLY: KTOLD,KBOLD,KODBOT,KODTOP
    USE GWFUNSFMODULE, ONLY: LPTIME,LPRINT,LHEAD
    USE GWFUNSFMODULE, ONLY: PROPR,PRTIME,vmid,con,epsi
    USE GWFUNSFMODULE, ONLY: HHOLD,HHNEW,HTEMP,HB,X,XSURF,NUMNP,THOLD,THNEW
    USE GWFUNSFMODULE, ONLY:LCHEM,KBOTCH,KTOPCH ,ths,tholdim,switch
    USE GWFUNSFMODULE, ONLY: num,ikod,s,r,p,q,vnew,vold,thnewim,nmat,beta,matnum,rbot,lwat
    USE GWFUNSFMODULE, ONLY:ths, sink,tlevel,ns,cmid,tpulse,nsd,llinear,retard,thsat,solfluxwat
    USE GWFUNSFMODULE, ONLY:  numnpd,numnpd,conc ,hh,wvoli,xmasschange,WCUMT,WCUMA,THORG,volume

    use msflib
    IMPLICIT NONE
    !      Variables
    INTEGER,INTENT(IN) :: IGRID,KSTP,KPER
    DOUBLE PRECISION :: ADEPTH,AHEAD,Rtime1,RTime,RTime2 ,dxa,dxb,va,vb
    INTEGER :: I,K,SKIP,INOD,nlay,l,Tswitch
    integer*2 i2,iYear,iMonth,iDay,iHours,iMins,iSecs,i100th
    !      Body
    CALL GWF2UNSF1PNT(IGRID)
    !      The following line is required to make some subroutines work
    !      correctly
    THETA=>THNEW
    !     CHECK: Nonsensical line in the original code at this point
    !            See line 687 of the original code for more information
    !SKIP=0
    !write(92,*) 'Time step:' ,KSTP, 'Stress Period:',kper

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
        WRITE(100,*)T(i)
        CALL GWF2UNSF1PROFUP1(I,-ABS(HCA(I,ALEVEL(I))))
        !change the wvoli to adjust for the water flow balance
        ! masschange(i)=xMassChange(i)
        plevel(I)=0
        volume=0.

        CALL GWF2UNSF1SETMAT(I)
        THOLD(:,I)=THNEW(:,I)
        !  vold(:,I)=vnew(:,I)
        SUMVBOT(I)=0.

        !start calculating the time required for simulation of water flow and Solute transport using HYDRUS
        !runstart time
        call getdat(iYear,iMonth,iDay)
        call gettim(iHours,iMins,iSecs,i100th)
        Rtime1=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)
        switch=0

        !  initime=t(pid)
        ! modified by AdamS - the loop below was missing in the previous version
         write (97,*) '-----------------------------------------------------------------------------'

        DO  ! loop over time steps for a single profile
            ! solve unsaturated flow for a single time step in the profile


            CALL GWF2UNSF1WF(I)

            !Calculate time-averaged VBOT
            !CALL GWF2UNSF1TOTF(I)
            !----------------------------------------------------------------
            CALL GWF2UNSF1VELOCAD(I)
            !----------------------------------------------------------------
            !  write (96,*)vbot(i)*dt
            SUMVBOT(I)=SUMVBOT(I)+VBOT(I)*DT(I)
            WRITE ( 97,'(4(a7,es22.8))') 'T= ', T(I), ' DT= ',DT(I),' VBOT= ', VBOT(I),'SUMV= ',SUMVBOT(I)
            !
            !Root zone calculations
            IF (SINKF(I).GE.0) THEN
                THETA=>THNEW
                CALL GWF2UNSF1SETSNK(I)
            ENDIF
            !Solute transport
            !---------------------------------------------------------------------------
            IF (lChem) then

                CALL GWF2UNSF1SOLUTE(i,kstp,KPER,thold,thnew,vold,vnew,p,r,s,q,htemp,tholdim,&
                    & thnewim,ikod,tpulse,retard)

            End if
            switch=1

            PLEVEL(I)=PLEVEL(I)+1
            !   Read time-dependent boundary condition
            IF(ABS(T(I)-TATM(I)).LE.0.001*DT(I)) THEN
                !   CHECK: Redundant IF statement at line 771 in the original code
                IF(INFTOP(I).GE.0 .OR. INFBOT(I).GE.0 .OR. ATMBC(I).GE.0..and..not.lChem) THEN
                    CALL GWF2UNSF1SETBC1(I)
                ELSE
                    CALL GWF2UNSF1SETBCC1(I)
                    ALEVEL(I)=ALEVEL(I)+1
                END IF
            END IF


            ! modified by AdamS
            HTEMP(:,I)=HHNEW(:,I)
            HHOLD(:,I)=HHNEW(:,I)
            THOLD(:,I)=THNEW(:,I)

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

        write(95,*) 'sumvbot',SUMVBOT(I) ,vbot(I)
        ! write (96,*)vbot(i)*DT(I)
        ! if (tswitch.eq.1)  write(9999,"(f10.2)")

        do 33 k=1,(numnp(i)-1)
            do 44 l=1,ns
                if (hhnew(k,I).ge.0.and. hhnew(k+1,i).le.0)then
                    DXA=  X(k+1,i)-X(k,i)
                    DXB=  X(k,i)-x(k-1,i)
                    VA=-(CON(k,i)+CON(k+1,i))/2.*((HHNEW(k+1,i)-HHNEW(k,i))/DXA+1)
                    VB=-(CON(k,i)+CON(k-1,i))/2.*((HHNEW(k,i)-HHNEW(k-1,i))/DXB+1)
                    Vmid(k,i)=(VA*DXB+VB*DXA)/(DXA+DXB)
                    solfluxwat(I)=(conc(l,k,i)-((conc(l,k,i)-conc(l,k+1,i))*((hhnew(k,i))/(hhnew(k,i)-hhnew(k+1,i)))))*(vmid(k,i))
                    write(92,"(I4,5X,f10.2,5X, E10.3,5X,E10.3,7x, I6,7x,I6)") i,Tmax(i), conc(l,k,i)-((conc(l,k,i)-conc(l,k+1,i))*((hhnew(k,i))/(hhnew(k,i)-hhnew(k+1,i)))) , solfluxwat(i),kstp,kper
                end if
44          continue
33      continue


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

    DO K=1,PROPR
        IF (ABS(PRTIME(K)-T(1)).LT.DT(1)) THEN
            CALL GWF2UNSF1NODOUT1(T(1),LHEAD)
            EXIT
        END IF
    END DO

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