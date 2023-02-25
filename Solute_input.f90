    !Input information: Solute transport
    !The following is taken from the original hydrus-1d source code and
    !modified here
    SUBROUTINE ChemIn(IN)
    USE GLOBAL, ONLY: IOUT,LSTCHK,LENUNI
    USE GWFUNSFMODULE, ONLY: lUpW,lTDep,lLinear,lEquil,lArtD,lScreen,lTort,&
        &lMobIm,lBact,lFiltr,lMoistDep,lDualNEq,lMassIni,&
        &lEqInit,lVar,ChPar,TDep,CcTop,CcBot,CumCh,WDep,Par,epsi,lUpW,lArtD,&
        &lTDep,cTolA,cTolR,MaxItC,PeCr,NS,lTort,npar2,tpulse,nmat,ktopch,kbotch,&
        &ltort,imoistdep,itort,imodel,ibact,iNonEqul,ns
    implicit none
    INTEGER,INTENT(IN) :: IN
    CHARACTER(LEN=100) :: LINE
    DOUBLE PRECISION :: FQ,par1(10)
    INTEGER :: M,J,jj,jjj,I,ierr
    !       lUpW:.true. if upstream weighing formulation is to be used.
    !       false. if the original Galerkin formulation is to be used.
    !       lTDep:.true. if at least one transport or reaction coefficient (ChPar) is temperature
    !       dependent..false. otherwise. If lTDep=.true., then all values of ChPar(i,M) should be specified at a
    !       reference temperature Tr=20oC.
    !       lLinear(NSD)
    !       lEquil
    !       lArtD.true. if artificial dispersion is to be added in order to fulfill the PeCr stability
    !       criterion (see Section 8.4.5). .false. otherwise.
    !       lScreen
    !       lTort: .true. if tortuosity factor [Millington and Quirk, 1961] is to be used.
    !       .false. if tortuosity factor is assumed to be equal to one.
    !       lMobIm(NMat)
    !       lBact,
    !       lFiltr
    !       lMoistDep
    !       lDualNEq
    !       lMassIni
    !       lEqInit
    !       lVar
    !  dimension ChPar(NSD*16+4,NMat),TDep(NSD*16+4),cTop(NSD),cBot(NSD),
    !         CumCh(10,NSD),WDep(2+NMatD,NSD*9),Par(11,NMat)
    ! Following are the parameters required to include all the solute transport
    !related details; all these will be the for the whole of the HPM and not defined
    !for individual hydrus-1d profile
    !ChPar(NSD*16+4,NMat)
    !TDep(NSD*16+4)
    !cTop(NSD)
    !cBot(NSD)
    !CumCh(10,NSD)
    !WDep(2+NMatD,NSD*9)
    !Par(11,NMat)

    write(IOUT,110,err=902)
    CALL UHRCOM(IN,IOUt,LINE)
    BACKSPACE IN
    ! if(iVer.le.2) then
    !   read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,&
    !&                     PeCr,NS,lTort

    ! here we are using the version greater than 2 hence we need both
    ! details of ibact and lfilt
    ! else
    read(IN,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,&
        &                     PeCr,NS,lTort,iBact,lFiltr
    ! Reading all the details in first line of the solute data

    lBact=.false.
    if(iBact.eq.1) lBact=.true. !conversion of 1 to true
    !   end if
    ! if(iVer.eq.4) then! no need to specify the version
    !   read(30,*,err=901)
    CALL UHRCOM(IN,IOUt,LINE)
    BACKSPACE IN
    read(IN,*,err=901) iNonEqul,lMoistDep,lDualNEq,lMassIni,lEqInit,&
        &                     lVar
    ! we are reading all the values till the tortousity
    if(lMoistDep) iMoistDep=1 !converting the true to 1
    if(lVar) iTort=1 !converting the true to 1
    ! end if
    PeCr=amax1(PeCr,0.1) !limiting the lower limit of pecr as 0.1
    if(lUpW) then
        write(IOUT,120,err=902)!upstream function is used
    else
        write(IOUT,130,err=902)!upstream function not used
        if(lArtD) write(IOUT,140,err=902) PeCr !tortousity added is the pecr is greater that pecr
    end if
    write(IOUT,150,err=902) lTDep,iMoistDep,cTolA,cTolR,MaxItC

    !lEquil=.true.
    do 11 M=1,NMat! number of material
        CALL UHRCOM(IN,IOUt,LINE)
        BACKSPACE IN
        read(IN,*,err=901) (ChPar(j,M),j=1,4)!bulkdensity, long. disp,frac,mobile immobile
        write(IOUT,160,err=902) M,(ChPar(j,M),j=1,4)
        if(ChPar(3,M).lt.1..or.ChPar(4,M).gt.0..or.lBact) lEquil=.false.
        lMobIm(M)=.false. !checking if the model is an equilibrium model or physical or chemi
        if(.not.lBact.and.ChPar(4,M).gt.0.) lMobIm(M)=.true.
        if(.not.lEquil.and.ChPar(1,M).eq.0.) goto 903
11  continue
    do 13 jj=1,NS! number of solute
        jjj=(jj-1)*16
        write(IOUT,170,err=902) jj
        CALL UHRCOM(IN,IOUt,LINE)
        BACKSPACE IN
        read(IN,*,err=901) (ChPar(jjj+j,1),j=5,6)!diffusion in water and gas! DifW       DifG
        write(IOUT,180,err=902) (ChPar(jjj+j,1),j=5,6)
        lLinear(jj)=.true. !For number of solute

        do 12 M=1,NMat !number of material
            ChPar(jjj+5,M)=ChPar(jjj+5,1)
            ChPar(jjj+6,M)=ChPar(jjj+6,1)
            CALL UHRCOM(IN,IOUt,LINE)
            BACKSPACE IN
            read(IN,*,err=901) (ChPar(jjj+j,M),j=7,20)!        Ks          Nu        Beta       Henry       SnkL1       SnkS1       SnkG1       SnkL1'      SnkS1'      SnkG1'      SnkL0       SnkS0       SnkG0        Alfa
            write(IOUT,190,err=902) M,(ChPar(jjj+j,M),j=7,20)
            if     (abs(ChPar(jjj+8,M)-0.0).gt.1.e-12) then
                write(IOUT,200,err=902) M
            else if(abs(ChPar(jjj+9,M)-1.0).gt.0.001) then
                write(IOUT,210,err=902) M
            else
                write(IOUT,220,err=902) M
            end if
            if(.not.lEquil) then
                if(lMobIm(M)) then !imobim is for number of material
                    write(IOUT,222,err=902)
                else
                    write(IOUT,224,err=902)
                end if
            end if
            if(abs(ChPar(jjj+8,M)-0.0).gt.1.e-12.or.&
                &      abs(ChPar(jjj+9,M)-1.0).gt.0.001) lLinear(jj)=.false.
            if(lBact.and.(ChPar(jjj+18,M).gt.0..or.ChPar(jjj+15,M).gt.0.))&
                &                                         lLinear(jj)=.false.! reading all the information regarding the adsorption coeff and checking ifthe model is linear/ non linear./ equili/ non-equil. etc
12      continue ! no. of material

13  continue! no. of solute


    do 14 jj=1,NS*16+4
        TDep(jj)=0.
        if(jj.le.NS*9) WDep(1,jj)=1.
        if(jj.le.NS*9) WDep(2,jj)=0.
14  continue


    do 16 jj=1,NS
        do 15 i=1,10
            CumCh(i,jj)=0.
15      continue
        CALL UHRCOM(IN,IOUt,LINE)
        BACKSPACE IN
        
        if(lTDep) then ! if temperture dep, parameters are read for the number of solute
            jjj=(jj-1)*16
            read(IN,*,err=901) (TDep(jjj+j),j=5,6)
            CALL UHRCOM(IN,IOUt,LINE)
            BACKSPACE IN
            write(iout,*)line!   DifW       DifG                n-th solute
            read(IN,*,err=901) (TDep(jjj+j),j=7,20)!Ks          Nu        Beta       Henry       SnkL1       SnkS1       SnkG1       SnkL1'      SnkS1'      SnkG1'      SnkL0       SnkS0       SnkG0        Alfa
        end if
16  continue

    ! if water content dep, paramters are read for no. of solute and material

    do 19 jj=1,NS

        if(iMoistDep.eq.1) then
            if(jj.eq.1)then
                CALL UHRCOM(IN,IOUt,LINE)
                BACKSPACE IN
                read(IN,*,err=901) nPar2!no. of parameters for water content dependence
                write(iout,*)'No. of parameters for water content dependence:', npar2
            end if
            jjj=(jj-1)*9
            CALL UHRCOM(IN,IOUt,LINE)
            BACKSPACE IN
            read(IN,*,err=901) (WDep(1,jjj+j),j=1,9)
            write(iout,*) (WDep(1,jjj+j),j=1,9)
            CALL UHRCOM(IN,IOUt,LINE)
            BACKSPACE IN
            read(IN,*,err=901) (WDep(2,jjj+j),j=1,9)
            write(iout,*) (WDep(2,jjj+j),j=1,9)
            do 18 M=1,NMat
                do 17 j=1,9
                    par1=Par(1,M)
                    WDep(2+M,jjj+j)=FQ(iModel,WDep(2,jjj+j),Par1)
17              continue
18          continue
        end if
19  continue

    CALL UHRCOM(IN,IOUt,LINE)
    BACKSPACE IN
    read(IN,*,err=901) kTopCh,(CcTop(jj),jj=1,NS),kBotCh,&
        &                   (CcBot(jj),jj=1,NS)
    ! Boundary condition at the top, concentration of the incoming fluid, similarly for the bottom
    ! if(kTopCh.eq.-2) then !ours will be always -1
    !   CALL UHRCOM(IN,IOUt,LINE)
    !BACKSPACE IN
    !   read(IN,*,err=901) dSurf,cAtm
    ! end if
    write(IOUT,230,err=902) kTopCh,(CcTop(jj),jj=1,NS)
    write(IOUT,240,err=902) kBotCh,(CcBot(jj),jj=1,NS)
    CALL UHRCOM(IN,IOUt,LINE)
    BACKSPACE IN
    read(IN,*,err=901) tPulse !pulse duration of the input of contaminant
    write(IOUT,250,err=902) tPulse
    return

    !     Error when reading from an input file
901 ierr=1
    return
    !     Error when writing into an output file
902 ierr=2
    return
    !     Bulk Density is equal to zero
903 ierr=3
    return

110 format(//' Solute transport information'/1X,28('='))
120 format(/' Upstream weighting finite-element method')
130 format(/' Galerkin finite-element method')
140 format (/' Artificial dispersion is added when Peclet number is',&
        &         ' higher than',f10.3)
150 format(//' lTDep     lWDep     cTolA     cTolR   MaxItC'&
        &        /l3,6x,i3,e13.3,f10.4,i7/&
        &        //' Mat.     Bulk.D.    DispL    Fraction  Immobile WC')
160 format(i3,f13.4,3f10.4)
170 format(/'    Dif.w.      Dif.g.   ',50('-'),' (',i2,'.solute)')
180 format(2e12.4/' Mat.     KS         Nu         Beta      Henry&
        &  SinkL1     SinkS1     SinkG1     SinkL1`    SinkS1`    SinkG1`&
        &  SinkL0     SinkS0     SinkG0      Alfa')
190 format(i4,14e11.4)
200 format(/' Langmuir nonlinear adsorption isotherm for material ',&
        &       i2)
210 format(/' Freundlich nonlinear adsorption isotherm for material ',&
        &       i2)
220 format(/' No adsorption or linear adsorp. isotherm for material ',&
        &       i2)
222 format(/' Physical non-equilibrium solute transport with mobile an&
        &d imobile water.')
224 format(/' Chemical non-equilibrium solute transport with kinetic a&
        &nd equilibrium sorption sites.')
230 format(/' kTopCh      cTop(1...NS)'/i4,7x,20e10.3)
240 format(/' kBotCh      cBot(1...NS)'/i4,7x,20e10.3)
250 format(/' tPulse =   ',f15.3)
    end
