


      Subroutine GWF2UNSF1SOLUTE(pid,kstp,kper,tho,thn,vo,vn,vcorr)

      use GWFUNSFMODULE,only: js,NMat,NS,NSD,x,dt,t, ChPar,
     !thNew,vOld, vNew,Disp,epsi, cBot,ctop,cTolR,IterC,TempN,
     !p,r,s,q ,cvTop,cvBot,cvCh0,conc,thnim,g0,retard,tconv,
     !cvCh1,Peclet,Courant,dtMaxC,TempO,Lnequil,LdualNeq,switch,
     !cNew,cPrevO,cTemp,TDep,thS,cTolA,iconctype, thoim,thsat,
     !MaxItC,htemp,Sorb,SorbN,kbotch,ktopch,Lbact,Lfiltr,tMAX,
     !PeCr,q0,q1,dSurf,cAtm,Sink,cRootMax,sSink,concnew1,
     !cvChR,lMobIm,cvChIm,TLevel,Sorb2,SorbN2,hhnew,con,tpulse,
     !dtMin,dtOpt,lWat,ThOldIm,ThNewIm,LLINEAR,solfluxwat,vmid,
     !SinkIm,STrans,xConv, lVapor,rBot,concmod,dsurft,nsd,idualpor,
     !iMoistDep,NMatD,DMoist,WDep,ik od,Beta, concnew,concold,
     !AtmBC,SinkF,lActRSU,OmegaS,OmegaW,SPot,THMOD,THORG,wc,
     !rKM,cMin,lDensity,numnp,numnpd,level,npunsf,matnum,par,
     !HHNEW,b,d,e,f,lupw, Lequil,Lartd,ltort,g1,thOld,iter

      INTEGER,INTENT(IN) :: PID, KPER,KSTP
      INTEGER ::  jjj
      double precision  tho(numnpd,npunsf),thn(numnpd,npunsf),
     !tr,rr,concc,PARAM(10),Vn(numnpd,npunsf),vo(numnpd,npunsf),bn,dn,
     !fn,e1,f1,d1,dxa,dxb,va,vb,vcorr(numnpd,npunsf)

      n=numnp(pid)
      do 1 i=1,nmat
          thsat(i)=par(2,i)
   1  end do
      
      iconctype=0
      alf=1.-epsi
      IterC=1.
      NLevel=2
      Peclet=0.
      Courant=0.
      dtMaxC=1.e+30
      rMin=1.e-30
      Tr=293.15
      Rr=8.314
!     Sequential first order decay goes into equilibrium phase (lNEquil=.false.) or nonequilbrium phase (lNEquil=.true.) 
      lNEquil=.false.
      do   i=1,n
          tempo(i,pid)=20
          tempn(i,pid)=20
      end do

   10 continue

*     Loop on species in the chain

      do 18 jS=1,NS
          iter(pid)=0
          jjj=(jS-1)*16
          cvTop(jS,pid)=0.
          cvBot(jS,pid)=0.
          cvCh0(jS,pid)=0.
          cvCh1(jS,pid)=0.
          cvChR(jS,pid)=0.
          cvChIm(jS,pid)=0.
          if(t(pid)-tPulse.gt.dtMin(pid).and..not.AtmBC(pid)) then
              cTop(pid,jS)=0.
              cBot(pid,jS)=0.
          end if
          if(kBotCh.lt.0) then
              if(vO(1,pid).ge.0.) cvBot(jS,pid)=alf*cBot(pid,jS)  *vO(1,pid)
              if(vO(1,pid).lt.0.) cvBot(jS,pid)=alf*Conc(jS,1,pid)*vO(1,pid)
              if(lVapor.and.rBot(pid).eq.0.) cvBot(jS,pid)=0.
          else if(kBotCh.eq.0) then
              cvBot(jS,pid)=alf*Conc(jS,1,pid)*vO(1,pid)

          end if
          if(kTopCh.lt.0..and.TLevel(pid).ne.1) then
              if(vO(N,pid).lt.0.) cvTop(jS,pid)=alf*cTop(pid,jS)*vO(N,pid)

          end if
          if(kTopCh.eq.-2) then
              M=MatNum(N,pid)
              Tr=293.15
              RR=8.314! R changed to RR since R is already in use in gwf2unsf
              TT=(TempO(N,pid)+273.15-Tr)/RR/(TempO(N,pid)+273.15)/Tr
              Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
              Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
              dSurfT=dSurf*exp(TDep(jjj+9)*TT)
              cvTop(jS,pid)=cvTop(jS,pid)+alf*Dg/dSurfT*Henry*Conc(jS,N,pid)-
     !        Dg/dSurfT*cAtm
          end if
          if(.not.lLinear(jS)) then
              do 11 i=1,N
                  cNew(i,PID)=Conc(jS,i,pid)
                  if(.not.lEquil)       SorbN(i,pid) =Sorb(jS,i,pid)
                  if(lBact.or.lDualNEq) SorbN2(i,pid)=Sorb2(jS,i,pid)
11            continue
          end if


*       Root Solute Uptake
          if(SinkF(pid))
     !    call SetSSnk(pid,sSink)

*       Iterative loop for a nonlinear adsorption isotherm
 12       iter(pid)=iter(pid)+1
          if(.not.lLinear(jS)) then
              do 13 i=1,N
                  cTemp(i,pid)=cNew(i,PID)
13            continue
          end if

!       To construct the matrix equation
          do 15 Level=1,NLevel

*         Calculate the dispersion coefficients, retardation factors, source/
*         decay coefficients, Peclet and Courant numbers, upstream weighting
*         factorscoeff 

          call Coeff(pid,Pecl,Cour,dtMxC,vo,vn,tho,thn,vcorr,retard)

          Peclet=amax1(Peclet,Pecl)
          Courant=amax1(Courant,Cour)
          dtMaxC=amin1(dtMaxC,dtMxC)

!*        Set up the matrix equation


          call MatSet(pid, jS,dSurfT,vo,vn,tho,thn,retard,ALF)

          do 14 i=1,N
              if(Level.eq.1) vO(i,pid)=vO(i,pid)+vCorr(i,pid)
              if(Level.eq.2) vN(i,pid)=vN(i,pid)+vCorr(i,pid)
 14       continue

!*        Calculate mass-transfer fluxes at the beginning of the time interval
          if(Level.eq.1.and.iter(pid).eq.1)
     !    call MassTran(pid,jS, TempO,alf,ThO,vO)
15        continue
!
!*       Solve matrix equation
          call BanSol(pid,B,D,E,F)

!*       Test for convergence for nonlinear problem
          lConv=.true.

          do 16 i=1,N

              if((NS.gt.1.and.iter(pid).eq.1).or.lDensity) cPrevO(i,pid)=
     !         Conc(jS,i,pid)
              if(lLinear(jS)) then
                  Conc(jS,i,pid)=amax1(sngl(F(i,pid)),0.)
                  if(Conc(jS,i,pid).lt.1.e-30.and.Conc(jS,i,pid).gt.0.)
     1            Conc(jS,i,pid)=0.
              else
                  cNew(i,PID)=sngl(F(i,pid))
                  if(cNew(i,PID).lt.1.0e-30) cNew(i,PID)=0.
                  if(abs(cNew(i,PID)-cTemp(i,pid)).gt.cTolA+cTolR*
     !            Conc(jS,i,pid))
     !            lConv=.false.
              end if

16        continue

          if(.not.lLinear(jS)) then
              if(.not.lConv) then
                  if(iter(pid).lt.MaxItC) then
                      goto 12
                  else if(dt(pid).gt.dtMin(pid).and..not.lWat(pid)) then
c              ierr=1
                      dtOld=dt(pid)
                      dt(pid)=amax1(dt(pid)/3.,dtMin(pid))
                      dtOpt=dt(pid)
                      t=t-dtOld+dt(pid)
                      goto 10
                  else
                      ierr=1
                  end if
              end if
              do 17 i=1,N
                  Conc(jS,i,pid)=cNew(i,pid)
                  if(.not.lEquil)       Sorb(jS,i,pid) =SorbN(i,pid)
                  if(lBact.or.lDualNEq) Sorb2(jS,i,pid)=SorbN2(i,pid)
17            continue
          end if

!------------------------------------------------------------------------          
!*       Calculate sorbed concentration for linear noneq. adsorption or 
!*       concentration in the imobile water.
          if(.not.lEquil.and.lLinear(jS))

     !    call SorbConc(PID,js,TempN,thN, vN,ThNIm)

!*       Calculate mass-transfer fluxes at the end of the time interval

          call MassTran(pid,jS,TempN,epsi,ThN,vN)


*       Set up mass fluxes
          if(kTopCh.lt.0) then
              if(TLevel(pid).ne.1) then
                  if(vN(N,pid).lt.0.) cvTop(jS,pid)=cvTop(jS,pid)+epsi*
     Q            vN(N,pid)*cTop(pid,jS)
              else
                  if(vN(N,pid).lt.0.) cvTop(jS,pid)=cvTop(jS,pid)+
     !            vN(N,pid)*cTop(pid,jS)
              end if
          else
              cvTop(jS,pid)=FN-BN*Conc(jS,N-1,pid)-DN*Conc(jS,N,pid)
          end if
          if(kTopCh.eq.-2) then
              M=MatNum(N,pid)
              TT=(TempN(N,pid)+273.15-Tr)/rr/(TempN(N,pid)+273.15)/Tr
              Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
              Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
              cvTop(jS,pid)=cvTop(jS,pid)+epsi*Dg/dSurfT*Henry*Conc(jS,N,pid)-
     !        Dg/dSurfT*cAtm
          end if

          if(kBotCh.lt.0) then
              if(vN(1,pid).ge.0.) cvBot(jS,pid)=cvBot(jS,pid)+epsi*
     !        cBot(pid,jS)  *vN(1,pid)
              if(vN(1,pid).lt.0.) cvBot(jS,pid)=cvBot(jS,pid)+epsi*
     !        Conc(jS,1,pid)*vN(1,pid)
              if(lVapor.and.rBot(pid).eq.0.) cvBot(jS,pid)=0.
          else if(kBotCh.eq.0) then
              cvBot(jS,pid)=cvBot(jS,pid)+epsi*Conc(jS,1,pid)*vN(1,pid)
          else
              cvBot(jS,pid)=D1*Conc(jS,1,pid)+E1*Conc(jS,2,pid)-F1
          end if
          IterC=max0(IterC,iter(pid))
          if(abs(cvTop(jS,pid)).lt.rMin) cvTop(jS,pid)=0.
          if(abs(cvBot(jS,pid)).lt.rMin) cvBot(jS,pid)=0.

1001      format (2e12.3)
1002      format (e12.3)
18    continue
      
          
*     Calculate flux concentrations
      if(iConcType.eq.2)
         

     !call FluxConc(pid, x,vN,thN,TempN,cNew,ThNIm,1)

      return
      end



      SUBROUTINE CONCMODIFY(I,PID,JS,M,JJJ,CC,param)
      USE GWFUNSFMODULE,ONLY: CHPAR,THORG,THMOD ,concnew
      INTEGER,INTENT(IN) :: PID
      double precision PARAM(10)

      SSS=param(7)* CC+(  param(1)*param(3)*CC**
     !param(5) )/(1.+param(4)*CC**param(5))
      concnew(js,i,pid) =CINIT(sss,PARAM,10)
      RETURN
      END SUBROUTINE

!      ************************************************************************

*     Evaluate Liquid concentration from the total solute mass

      real function cInit(xCONC,ParAM,NPar)

      double precision ParAM(NPar)
      x1=1.e-10
      x2=1.e+10
      call ZBRAK1(X1,X2,XB1,XB2,xCONC,ParAM,NPar)
      cInit=ZBRENT1(XB1,XB2,xCONC,ParAM,NPar)

      return
      end

*************************************************************************

      real function SolMass(ConcC,xCONC,ParAM,NPar)

*     Calculate total solute mass for concentration Conc

      double precision ParAM(NPar)

      ro=ParAM(1)
      frac=ParAM(2)
      xKs=ParAM(3)
      xNu=ParAM(4)
      fExp=ParAM(5)
      !   xKH=Par(6)
      !   Theta=Par(7)
      ThetaA=ParAM(6)

      yCONC=ThetaA*ConcC+ ro*xKs*ConcC**fExp /(1.+xNu*ConCc**fExp)


      SolMass=xCONC-yCONC

      return
      end

************************************************************************

*     Bracketing of the root, Numerical recepies (345)


      subroutine ZBRAK1(X1,X2,XB1,XB2,xMass,param,NPar)

      double precision param(NPar)

      NBB=1
      NB=1000

      dlh=(alog10(X2)-alog10(X1))/(NB-1)
      FP=SolMass(X1,xMass,param,NPar)
      do 11 i=1,NB
          dx2=alog10(X1)+(i)*dlh
          X2=10**dx2
          FC=SolMass(X2,xMass,param,NPar)
          if(FC*FP.lt.0.) then
              NBB=NBB+1
              XB1=X1
              XB2=X2
              return
          end if
          FP=FC
          X1=X2
          if(NBB.eq.NB) return
11    continue

      return
      end

************************************************************************

*     Brent method of finding root that lies between x1 and x2, 
*     Numerical recepies (354)

      real function ZBRENT1(X1,X2,xMass,ParAM,NPar)

      parameter (ITMAX=100,EPS=3.E-8,TOL=1.e-6)
      double precision ParAM(NPar)
     !! double precision solmass
      A=X1
      B=X2
      FA=SolMass(A,xMass,ParAM,NPar)

      FB=SolMass(B,xMass,ParAM,NPar)
      IF(FB*FA.GT.0.) PAUSE 'Root must be bracketed for ZBRENT1.FB*FA<0'
      FC=FB
      DO 11 ITER=1,ITMAX
          IF(FB*FC.GT.0.) THEN
              C=A
              FC=FA
              D=B-A
              E=D
          ENDIF
          IF(ABS(FC).LT.ABS(FB)) THEN
              A=B
              B=C
              C=A
              FA=FB
              FB=FC
              FC=FA
          ENDIF
          TOL1=2.*EPS*ABS(B)+0.5*TOL
          XM=.5*(C-B)
          IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
              ZBRENT1=B
              RETURN
          ENDIF
          IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
              S=FB/FA
              IF(A.EQ.C) THEN
                  P=2.*XM*S
                  Q=1.-S
              ELSE
                  Q=FA/FC
                  R=FB/FC
                  P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
                  Q=(Q-1.)*(R-1.)*(S-1.)
              ENDIF
              IF(P.GT.0.) Q=-Q
              P=ABS(P)
              IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                  E=D
                  D=P/Q
              ELSE
                  D=XM
                  E=D
              ENDIF
          ELSE
              D=XM
              E=D
          ENDIF
          A=B
          FA=FB
          IF(ABS(D) .GT. TOL1) THEN
              B=B+D
          ELSE
              B=B+SIGN(TOL1,XM)
          ENDIF
          FB=SolMass(B,xMass,ParAM,NPar)
11    CONTINUE
      PAUSE 'ZBRENT1 exceeding maximum iterations.'
      ZBRENT1=B

      RETURN
      END function

      subroutine coeff (pid,Peclet,Courant,dtMaxC,vo,vn,tho,thn,vcorr,retard)

      use GWFUNSFMODULE,only:iter,js, nmat,NumNP,numnpd,npunsf,NSD,nmatd,
     !idualpor,lupw,llinear, Lequil,lartd,ltort,Lbact,Lfiltr,Lnequil,ldualneq,
     !thoim, thnim,thimobo,thimob,thwo,thw,v,SS1,SS2,thg,dw,dg,frac,xnu,xks,ro,
     !dsconc,fexp,gamli,gaml1i,gamloi,gaml1oi,gamsi,gams1i,gamsoi,gams1oi,g0,x,
     !conc,tempo,vj,thj,tempn, tdep,cnew,sorb,sorbn,ssink,sorb2,sorbn2,disp,
     !xconv,DMOIST,Wdep,dt,q1,q0,tconv,cprevo,g1,matnum,thsat,wc,chpar,
     !sinkim,henryj,omegas,omega,omegao,f_em,gamg1i,rka2,rkd2,smax1,rka1,rkd1,
     !ipsi2,ipsi1,dc,aa,psi1,psi2,dp,alfa1,alfa2,xksp,xnup,henryp,lmobim,
     !gaml1p,gaml1pi ,gams1pi,gams1p,gamg1p,tto,xkso,xnuo,fexpo,henryo,gamlo,
     !DERK,DX,DMOBI,CG1,CG,SSORB2,FLMACRO,SSORB,TAUG,SS,CC,DRETARDS,smax2,
     !SCONCPS,DCONCS,DSCONCS,SCONCOS,SCONCS,DRETARD,SCONCO,SCONC,SCONCP,DCONC,
     !smax2o,smax1o,rka2o,rka1o ,psi1o,psi2o ,xmuso,xmus,xmulo,xmul,rkd1o,rkd2o
      INTEGER,INTENT(IN) :: PID
      Real ::Peclet,Courant,dtmaxc
      INTEGER:: jjj,j,m
      double precision:: Vn(numnpd,npunsf),vo(numnpd,npunsf ),
     ! tho(numnpd,npunsf),thn(numnpd,npunsf),vcorr(numnpd,npunsf),
     !retard(numnpd,npunsf)

*     Inicialization
      jjj=(jS-1)*16
      if(jS.gt.1) jj1=jjj-16

      Peclet=0.
      Courant=0.
      CourMax=1.
      dtMaxC=1.e+30
      Tr=293.15
      rr=8.314

      do 11 i=numnp(pid),1,-1
          j=i+1
          k=i-1
          M=MatNum(i,pid)
          if(Level.eq.NLevel) then
              ThW=ThN(i,pid)
              ThWO=ThO(i,pid)
              ThG=amax1(0.,thSat(M)-ThW)
              if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
                  ThImob=ChPar(4,M)
                  ThImobO=ThImob
                  ThW=amax1(ThW-ThImob,0.001)
              end if
              if(iDualPor.gt.0) then
                  ThImob=ThNIm(i,pid)
                  ThImobO=ThOIm(i,pid)
              end if
              v=vN(i,pid)
              if(i.ne.numnp(pid)) then
                  vj=vN(j,pid)
                  Thj=ThN(j,pid)
                  if(lMobIm(M).and.iDualPor.eq.0.or.lBact)
     !            Thj=amax1(Thj-ThImob,0.001)
              end if
              TT=(TempN(i,pid)+273.15-Tr)/rr/(TempN(i,pid)+273.15)/Tr
              if(jS.gt.1) cPrev=Conc(jS-1,i,pid)
          else
              ThW=ThO(i,pid)

              ThG=amax1(0.,thSat(M)-ThW)
              if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
                  ThImob=ChPar(4,M)
                  ThW=amax1(ThW-ThImob,0.001)
              end if
              if(iDualPor.gt.0) then
                  ThImob=ThOIm(i,pid)
              end if
              v=vO(i,pid)
              if(i.ne.numnp(pid)) then
                  vj=vO(j,pid)
                  Thj=ThO(j,pid)
                  if(lMobIm(M).and.iDualPor.eq.0.or.lBact)
     !            Thj=amax1(Thj-ThImob,0.001)
              end if
              TT=(TempO(i,pid)+273.15-Tr)/rr/(TempO(i,pid)+273.15)/Tr
              if(jS.gt.1) cPrev=cPrevO(i,pid)
          end if


*       Temperature dependence
          f1=1.
          ro   =ChPar(1,     M)*exp(TDep(1)     *TT)
          Frac =ChPar(3,     M)*exp(TDep(3)     *TT)
          Dw   =ChPar(jjj+ 5,M)* exp(TDep(jjj+ 5)*TT)
          Dg   =ChPar(jjj+ 6,M)* exp(TDep(jjj+ 6)*TT)
          xKs  =ChPar(jjj+ 7,M)* exp(TDep(jjj+ 7)*TT)
          xNu  =ChPar(jjj+ 8,M)* exp(TDep(jjj+ 8)*TT)
          fExp =ChPar(jjj+ 9,M)
          Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,1,DMoist,1,WDep,ThW,iMoistDep)
          GamL  =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImob,iMoistDep)
          GamLi =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1    ! reaction in the immobile phase
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,2,DMoist,2,WDep,ThW,iMoistDep)
          GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImob,iMoistDep)
          GamSi =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,3,DMoist,3,WDep,ThW,iMoistDep)
          GamG  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)*f1
          GamGi=GamG
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,4,DMoist,4,WDep,ThW,iMoistDep)
          GamL1 =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImob,iMoistDep)
          GamL1i=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,5,DMoist,5,WDep,ThW,iMoistDep)
          GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImob,iMoistDep)
          GamS1i=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,6,DMoist,6,WDep,ThW,iMoistDep)
          GamG1 =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,7,DMoist,7,WDep,ThW,iMoistDep)
          xMuL  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,8,DMoist,8,WDep,ThW,iMoistDep)
          xMuS  =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,9,DMoist,9,WDep,ThW,iMoistDep)
          xMuG  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)*f1
          Omega =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
          f_em=1.
          if(lDualNEq) f_em  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
          if(lDualNEq) OmegaS=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
          if(lBact) then
              Dg=0.
              GamG =0.
              GamGi=0.
              GamL1=0.
              GamL1i=0.
              GamS1=0.
              GamS1i=0.
              GamG1=0.
              GamG1i=0.
              xMuL =0.
              xMuS =0.
              xMuG =0.
              Omega=0.
              SMax2 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)
              rKa2  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
              rKd2  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
              SMax1 =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TT)
              rKa1  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
              rKd1  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
              iPsi1=0
              iPsi2=0
              if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,M))
              if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,M))
              if(iPsi1.eq.0.and.SMax1.gt.0.) iPsi1=1
              if(iPsi2.eq.0.and.SMax2.gt.0.) iPsi2=1
              if(iPsi1.ge.3.or.iPsi2.ge.3)
     !        Dc=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
              if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,M)
              if(Level.eq.NLevel) then
                  ss1=SorbN(i,pid)
                  ss2=SorbN2(i,pid)
              else
                  ss1=Sorb(jS,i,pid)
                  ss2=Sorb2(jS,i,pid)
              end if
              psi1=1.
              psi2=1.
              if(iPsi1.gt.0) call Blocking(iPsi1,SMax1,psi1,x(i,pid),ss1,dc,aa)
              if(iPsi2.gt.0) call Blocking(iPsi2,SMax2,psi2,x(i,pid),ss2,dc,aa)

*         recalculate ka1 and ka2 based on filtration theory
              if(lFiltr) then
                  GamG =0.
                  GamL1=0.
                  Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
                  Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
                  Alfa1=rKa1
                  Alfa2=rKa2
                  call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,ThW,v,TempN(i,pid),
     !            xConv,tConv)
              end if
          end if

          if(jS.gt.1) then
              xKsP  =ChPar(jj1+ 7,M)*exp(TDep(jj1+ 7)*TT)
              xNuP  =ChPar(jj1+ 8,M)*exp(TDep(jj1+ 8)*TT)
              fExpP =ChPar(jj1+ 9,M) !*exp(TDep(jj1+ 9)*TT)
              HenryP=ChPar(jj1+10,M)*exp(TDep(jj1+10)*TT)
              f1=1.
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS-1,4,DMoist,4,WDep,ThW,iMoistDep)
              GamL1P =ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS-1,12,DMoist,4,WDep,ThImob,iMoistDep)
              GamL1Pi=ChPar(jj1+14,M)*exp(TDep(jj1+14)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS-1,5,DMoist,5,WDep,ThW,iMoistDep)
              GamS1P =ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS-1,13,DMoist,5,WDep,ThImob,iMoistDep)
              GamS1Pi=ChPar(jj1+15,M)*exp(TDep(jj1+15)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS-1,6,DMoist,6,WDep,ThW,iMoistDep)
              GamG1P=ChPar(jj1+16,M)*exp(TDep(jj1+16)*TT)*f1
              if(lBact) then
                  GamL1P=0.
                  GamL1Pi=0.
                  GamS1P=0.
                  GamS1Pi=0.
                  GamG1P=0.
              end if
          end if
          if(Level.eq.NLevel) then
              TTO=(TempO(i,pid)+273.15-Tr)/rr/(TempO(i,pid)+273.15)/Tr
              xKsO  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TTO)
              xNuO  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TTO)
              fExpO =ChPar(jjj+ 9,M)
              HenryO=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTO)
              f1=1.
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,1,DMoist,1,WDep,ThO(i,pid),iMoistDep)
              GamLO =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImobO,iMoistDep)
              GamLOi=ChPar(jjj+11,M)*exp(TDep(jjj+11)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,2,DMoist,2,WDep,ThO(i,pid),iMoistDep)
              GamSO =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImobO,iMoistDep)
              GamSOi=ChPar(jjj+12,M)*exp(TDep(jjj+12)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,4,DMoist,4,WDep,ThO(i,pid),iMoistDep)
              GamL1O =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImobO,iMoistDep)
              GamL1Oi=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,5,DMoist,5,WDep,ThO(i,pid),iMoistDep)
              GamS1O=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImobO,iMoistDep)
              GamS1Oi=ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,7,DMoist,7,WDep,ThO(i,pid),iMoistDep)
              xMuLO =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TTO)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,8,DMoist,8,WDep,ThO(i,pid),iMoistDep)
              xMuSO =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TTO)*f1
              OmegaO=ChPar(jjj+20,M)*exp(TDep(jjj+20)*TTO)
              dKs   =(xKs  -  xKsO)/dt(pid)
              dNu   =(xNu  -  xNuO)/dt(pid)
              ddExp =(fExp - fExpO)/dt(pid)
              dHenry=(Henry-HenryO)/dt(pid)
              if(i.ne.1)     TTi=(TempN(k,pid)+273.15-Tr)/rr/
     !        (TempN(k,pid)+273.15)/Tr
              if(i.ne.numnp(pid)) TTj=(TempN(j,pid)+273.15-Tr)/rr/
     !        (TempN(j,pid)+273.15)/Tr
              if(lBact) then
                  GamS1O=0.
                  GamS1Oi=0.
                  xMuLO =0.
                  xMuSO =0.
                  xMuGO =0.
                  OmegaO=0.
                  SMax2O =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TTO)
                  rKa2O  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TTO)
                  rKd2O  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TTO)
                  SMax1O =ChPar(jjj+18,M)*exp(TDep(jjj+18)*TTO)
                  rKa1O  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TTO)
                  rKd1O  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TTO)
                  iPsi1=0
                  iPsi2=0
                  if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,M))
                  if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,M))
                  if(iPsi1.eq.0.and.SMax1O.gt.0.) iPsi1=1
                  if(iPsi2.eq.0.and.SMax2O.gt.0.) iPsi2=1
                  if(iPsi1.ge.3.or.iPsi2.ge.3)
     !            Dc=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
                  if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,M)
                  psi1O=1.
                  psi2O=1.
                  if(iPsi1.gt.0)
     !            call Blocking(iPsi1,SMax1O,psi1O,x(i,pid),
     !            Sorb(jS,i,pid),dc,aa)
                  if(iPsi2.gt.0)
     !            call Blocking(iPsi2,SMax2O,psi2O,x(i,pid),
     !            Sorb2(jS,i,pid),dc,aa)
                  if(lFiltr) then
                      GamL1O=0.
                      Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TTO)
                      Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TTO)
                      Alfa1=rKa1O
                      Alfa2=rKa2O
                      call Deposit(rKa1O,rKa2O,Dc,Dp,Alfa1,Alfa2,ThWO,vO(i,pid),
     !                TempO(i,pid),xConv,tConv)
                  end if
              end if
          else
              TTN=(TempN(i,pid)+273.15-Tr)/rr/(TempN(i,pid)+273.15)/Tr
              xKsN  =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TTN)
              xNuN  =ChPar(jjj+ 8,M)*exp(TDep(jjj+ 8)*TTN)
              fExpN =ChPar(jjj+ 9,M) !*exp(TDep(jjj+ 9)*TTN)
              HenryN=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TTN)
              dKs   =(xKsN  -  xKs)/dt(pid)
              dNu   =(xNuN  -  xNu)/dt(pid)
              ddExp =(fExpN - fExp)/dt(pid)
              dHenry=(HenryN-Henry)/dt(pid)
              if(i.ne.1)     TTi=(TempO(k,pid)+273.15-Tr)/rr/
     !        (TempO(k,pid)+273.15)/Tr
              if(i.ne.numnp(pid)) TTj=(TempO(j,pid)+273.15-Tr)/rr/
     !        (TempO(j,pid)+273.15)/Tr
          end if
          if(i.ne.1) Henryi=ChPar(jjj+10,MatNum(k,pid))*exp(TDep(jjj+10)*TTi)
          if(i.ne.numnp(pid))
     !    Henryj=ChPar(jjj+10,MatNum(j,pid))*exp(TDep(jjj+10)*TTj)

          dSConc=1.
          dConc=1.
          SConcP=1.
          SConc=1.
          SConcO=1.
          dRetard=0.

          SConcS=1.
          SConcOS=1.
          dSConcS=1.
          dConcS=1.
          SConcPS=1.
          dRetardS=0.

*       Effects of nonlinear adsorption
          if(.not.lLinear(jS)) then
              cc=Conc(jS,i,pid)
              cMid=(Conc(jS,i,pid)+cNew(i,PID))/2.
              if(Level.eq.NLevel) cc=cNew(i,PID)
              if(cc.gt.0.) then
                  dSConc=fExp*cc**(fExp-1.)/(1.+xNu*cc**fExp)**2
                  SConc =     cc**(fExp-1.)/(1.+xNu*cc**fExp)
              end if
              if(cMid.gt.0.) then
                  dConc=fExp*cMid**(fExp-1.)/(1.+xNu*cMid**fExp)**2
                  dRetard=cMid**fExp/(1.+xNu*cMid**fExp)*dKs-
     !            xKs*cMid**(2.*fExp)/(1.+xNu*cMid**fExp)**2*dNu+
     !            xKs*(alog(cMid))*cMid**fExp/(1.+xNu*cMid**fExp)**2*ddExp
              end if
              if(Level.eq.NLevel.and..not.lEquil.and.Conc(jS,i,pid).gt.0.)
     !        SConcO=Conc(jS,i,PID)**(fExpO-1.)/(1.+xNuO*Conc(jS,i,pid)**fExpO)
              if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
                  ss=Sorb(jS,i,pid)
                  sMid=(Sorb(jS,i,pid)+SorbN(i,pid))/2.
                  if(Level.eq.NLevel) ss=SorbN(i,pid)
                  if(ss.gt.0.) then
                      dSConcS=fExp*ss**(fExp-1.)/(1.+xNu*ss**fExp)**2
                      SConcS =     ss**(fExp-1.)/(1.+xNu*ss**fExp)
                  end if
                  if(sMid.gt.0.) then
                      dConcS=fExp*sMid**(fExp-1.)/(1.+xNu*sMid**fExp)**2
                      dRetardS=sMid**fExp/(1.+xNu*sMid**fExp)*dKs-
     !                xKs*sMid**(2.*fExp)/(1.+xNu*sMid**fExp)**2*dNu+
     !                xKs*alog(sMid)*sMid**fExp/(1.+xNu*sMid**fExp)**2*ddExp
                  end if
                  if(Level.eq.NLevel.and..not.lEquil.and.Sorb(jS,i,pid).gt.0.)
     !            SConcOS=Sorb(jS,i,pid)**(fExpO-1.)/(1.+xNuO*Sorb(jS,i,pid)
     !            **fExpO)
              end if
          else
              if(Conc(jS,i,pid).gt.0.) dRetard=Conc(jS,i,pid)*dKs
              if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
                  if(Sorb(jS,i,pid).gt.0) dRetardS=Sorb(jS,i,pid)*dKs
              end if
          end if
          if(jS.gt.1) then
              if(.not.lLinear(jS-1)) then
                  if(cPrev.gt.0.)
     !            SConcP=cPrev**(fExpP-1.)/(1.+xNuP*cPrev**fExpP)
                  if(Sorb(jS-1,i,pid).gt.0.)
     !            SConcPS=Sorb(jS-1,i,pid)**(fExpP-1.)/
     !            (1.+xNuP*Sorb(jS-1,i,pid)**fExpP)
              end if
          end if
!
!*       Calculate the retardation factors
      Retard(i,pid)=(ro*Frac*f_em*xKs*dConc+ThG*Henry)/ThW+1.

     
!*       Calculate the dispersion coefficients
      call Disper(pid,i,retard)

*       Calculate the adsorbed concentration on kinetic sites or
*       the concentration in an imobile zone, before solving matrix equation
      if(.not.lEquil)
     !call NEquil(pid,I,GamLi,GamL1i,GamLOi,GamL1Oi,GamSi,GamS1i,GamSOi,
     !GamS1Oi)
      !
*       Calculate zero-order coefficient g0


      g0(i,pid)=xMuL*ThW+Frac*f_em*ro*xMuS+ThG*xMuG-sSink(i,pid)

      q0(i,pid)=xMuL*ThW+          ro*xMuS+ThG*xMuG
      if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
              g0(i,pid)=g0(i,pid)+Omega*SSorb
              if(iDualPor.gt.0.and.SinkIm(i,pid).le.0) g0(i,pid)=g0(i,pid)-FlMacro
              if(lDualNEq) g0(i,pid)=g0(i,pid)+OmegaS*ro*SSorb2
          else if(.not.lBact) then
              g0(i,pid)=g0(i,pid)+Omega*ro*SSorb
          else if(lBact) then
              g0(i,pid)=g0(i,pid)+rKd1*ro*SSorb+rKd2*ro*SSorb2
          end if
      end if
      if(jS.gt.1) then
          cG=cPrev*(GamL1P*ThW+ro*Frac*f_em*xKsP*GamS1P*SConcP+
     !    ThG*HenryP*GamG1P)
          cG1=cG
          if(.not.lEquil) then
              if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
                  aa=Sorb(jS-1,i,pid)*(ThImob*GamL1Pi+
     !            (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
                  if(.not.lNEquil) cG=cG+aa
                  cG1=cG1+aa
                  if(lDualNEq) then
                      aa=GamS1Pi*ro*Sorb2(jS-1,i,pid)
                      if(.not.lNEquil) cG=cG+aa
                      cG1=cG1+aa
                  end if
              else if(.not.lBact) then
                  aa=GamS1Pi*ro*Sorb(jS-1,i,pid)
                  if(.not.lNEquil) cG=cG+aa
                  cG1=cG1+aa
              else if(lBact) then
                  write(*,*) 'Attachment/dettachment model is implemented
     !            only for one solute'
                  write(*,*)'Press Enter to continue'
                  read(*,*)
                  stop
              end if
          end if
          g0(i,pid)=g0(i,pid)+cG
          q0(i,pid)=q0(i,pid)+cG1
      end if
      if(cMid.gt.0.) g0(i,pid)=g0(i,pid)-ro*Frac*f_em*dRetard

*       Calculate first-order coefficient g1
      g1(i,pid)=-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-
     !(GamG+GamG1)*ThG*Henry
c        if(Level.eq.NLevel) g1(i)=g1(i)-ThG*dHenry-Henry*(ThWO-ThW)/dt(pid)
      if(.not.lEquil) then
      if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
          g1(i,pid)=g1(i,pid)-Omega
          if(iDualPor.gt.0.and.SinkIm(i,pid).gt.0) g1(i,pid)=g1(i,pid)-
     !    SinkIm(i,pid)
          if(Level.eq.NLevel.and.lLinear(jS))
     !    g1(i,pid)=g1(i,pid)+Omega*dt(pid)*Omega/DMobI
          if(lDualNEq) then
              g1(i,pid)=g1(i,pid)-OmegaS*ro*Frac*(1.-f_em)*SConc*xKs
              if(Level.eq.NLevel.and.lLinear(jS)) g1(i,pid)=g1(i,pid)+OmegaS*ro*
     !        (dt(pid)*OmegaS*Frac*(1.-f_em)*xKs/
     !        (2.+dt(pid)*(OmegaS+GamSi+GamS1i)))
          end if
      else if(.not.lBact) then                             ! two-site sorption model
          g1(i,pid)=g1(i,pid)-Omega*ro*(1.-Frac)*SConc*xKs
          if(Level.eq.NLevel.and.lLinear(jS)) g1(i,pid)=g1(i,pid)+Omega*ro*
     !    (dt(pid)*Omega*(1.-Frac)*xKs/(2.+dt(pid)*(Omega+GamSi+GamS1i)))
      else if(lBact) then                                  ! filtration model
          g1(i,pid)=g1(i,pid)-ThW*(rKa1*psi1+rKa2*psi2)
          if(Level.eq.NLevel.and.lLinear(jS)) g1(i,pid)=g1(i,pid)+dt(pid)*ThW*
     !    (rKd1*rKa1/(2.+dt(pid)*(rKd1+GamSi))+
     !    rKd2*rKa2/(2.+dt(pid)*(rKd2+GamSi)))
      end if
      end if
      q1(i,pid)=(-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-
     !(GamG+GamG1)*ThG*Henry)*Conc(jS,i,pid)
      if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
              q1(i,pid)=q1(i,pid)-Sorb(jS,i,pid)*(ThImob*(GamLi+GamL1i)+
     !        (1.-Frac)*ro*xKs*SConcS*(GamSi+GamS1i))
              if(lDualNEq) q1(i,pid)=q1(i,pid)-(GamSi+GamS1i)*ro*Sorb2(jS,i,pid)
          else if(.not.lBact) then
              q1(i,pid)=q1(i,pid)-(GamSi+GamS1i)*ro*Sorb(jS,i,pid)
          else if(lBact) then
              q1(i,pid)=q1(i,pid)-ro*(GamSi+GamS1i)*(Sorb(jS,i,pid)+
     !        Sorb2(jS,i,pid))
          end if
      end if

*       Velocity corrections
      if(i.eq.1) then
          dx=x(2,pid)-x(1,pid)
          derK=(Henryj-Henry)/dx
      else if(i.eq.numnp(pid)) then
          dx=x(numnp(pid),pid)-x(numnp(pid)-1,pid)
          derK=(Henry-Henryi)/dx
      else
          dx=(x(j,pid)-x(k,pid))/2.
          derK=(Henryj-Henryi)/dx
      end if
      vCorr(i,pid)=ThG*Dg*TauG*derK
      if(Level.eq.1)      vO(i,pid)=vO(i,pid)-vCorr(i,pid)
      if(Level.eq.NLevel) vN(i,pid)=vN(i,pid)-vCorr(i,pid)
!

      call PeCour(pid, i,j)
11    continue
      return
      end


      subroutine MatSet(pid,js,dSurf,vo,vn,tho,thn,retard,ALF)

      use GWFUNSFMODULE, only:iDualPor,lVapor,rBot,lBact,NS,
     !cBot, x,Conc,Disp,matnum,wc,kBotCh,kTopCh,dt,epsi,
     !g0,g1,B,D,E,F,E1,D1,F1,BN,DN,FN,NMat,ChPar,ctop,Level,
     !TempO,TempN,TDep,cAtm,  lMobIm,numnp,numnpd,npunsf,nsd

      INTEGER,INTENT(IN) :: PID
      double precision:: dsurf,Retard(numnpd,NPUNSF), Vn(numnpd,npunsf),
     ! vo(numnpd,npunsf ),tho(numnpd,npunsf),thn(numnpd,npunsf)
     
      N=numnp(pid)
      do 10 i=1,N
          M=MatNum(i,pid)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
              ThImob=ChPar(4,M)
              if(ThImob.gt.thO(i,pid)) write(*,*) "Warning !!! ThImob > Theta"
              thN(i,pid)=max(thN(i,pid)-ThImob,0.001)
              thO(i,pid)=max(thO(i,pid)-ThImob,0.001)
          end if
10    continue

!Lower boundary condition
      b1=x(2,pid)-x(1,pid)
      if(Level.eq.1) then
          F1=           Conc(jS,1,pid)*
     !    (b1/2./dt(pid)*thO(1,pid)*Retard(1,pid)+
     !    alf*(-(thO(1,pid)*Disp(1,pid)+thO(2,pid)*Disp(2,pid))/b1/2.-
     !    ((2.+3.*wc(1,pid))*vO(1,pid)+vO(2,pid))/6.+
     !    b1/12.*(3.*g1(1,pid)+g1(2,pid))))+
     !    Conc(jS,2,pid)*
     !    alf*((thO(1,pid)*Disp(1,pid)+thO(2,pid)*Disp(2,pid))/b1/2.-
     !    (vO(1,pid)+(2.-3.*wc(1,pid))*vO(2,pid))/6.+b1/12.*
     !    (g1(1,pid)+g1(2,pid)))+alf*b1/6.*(2.*g0(1,pid)+g0(2,pid))


*       3. type  BC
          if(kBotCh.eq.-1) F(1,pid)=F1+alf*cBot(jS,pid)*vO(1,pid)
      else
          E1=epsi*(-(thN(1,pid)*Disp(1,pid)+thN(2,pid)*Disp(2,pid))/b1/2.+
     !    (vN(1,pid)+(2.-3.*wc(1,pid))*vN(2,pid))/6.-b1/12.*
     !    (g1(1,pid)+g1(2,pid)))
          D1=b1/2./dt(pid)*thN(1,pid)*Retard(1,pid)+
     !    epsi*((thN(1,pid)*Disp(1,pid)+thN(2,pid)*Disp(2,pid))/b1/2.+
     !    ((2.+3.*wc(1,pid))*vN(1,pid)+vN(2,pid))/6.-b1/12.*
     !    (3.*g1(1,pid)+g1(2,pid)))
          F2=epsi*b1/6.*(2.*g0(1,pid)+g0(2,pid))
          F1=F1+F2

*       1.type BC
          if(kBotCh.eq.1) then
              D(1,pid)=1.
              E(1,pid)=0.
              F(1,pid)=cBot(pid,jS)
          end if

*       3. type  BC
          if(kBotCh.eq.-1) then
              if(vN(1,pid).gt.0..or.(lVapor.and.rBot(pid).eq.0.)) then
                  E(1,pid)=E1
                  D(1,pid)=D1
                  F(1,pid)=F(1,pid)+F2+epsi*cBot(pid,jS)*vN(1,pid)
              else
                  D(1,pid)=-1.
                  E(1,pid)=1.
                  F(1,pid)=0.
              end if
          end if

*       Free drainage
          if(kBotCh.eq.0) then
              D(1,pid)=-1.
              E(1,pid)=1.
              F(1,pid)=0.
          end if
      end if

      do 11 i=2,N-1
          a1=b1
          b1=x(i+1,pid)-x(i,pid)
          dx=(x(i+1,pid)-x(i-1,pid))/2.

          if(Level.eq.1) then

                               F(i,pid)=       Conc(jS,i-1,pid)*
     !       alf*((thO(i-1,pid)*Disp(i-1,pid)+thO(i,pid)*Disp(i,pid))/a1/2.+
     !            ((2.+3.*wc(i-1,pid))*vO(i-1,pid)+vO(i,pid))/6.+
     !            a1/12.*(g1(i-1,pid)+g1(i,pid)))+
     !                Conc(jS,i,pid)*
     !      (dx/dt(pid)*thO(i,pid)*Retard(i,pid)+
     !      alf*(-(thO(i-1,pid)*Disp(i-1,pid)+thO(i,pid)*Disp(i,pid))/a1/2.-
     !           (thO(i+1,pid)*Disp(i+1,pid)+thO(i,pid)*Disp(i,pid))/b1/2.-
     !           (vO(i+1,pid)+3.*(wc(i-1,pid)+wc(i,pid))*vO(i,pid)-
     !                   vO(i-1,pid))/6.+
     !           (a1*(g1(i-1,pid)+3.*g1(i,pid))+b1*(3.*g1(i,pid)+
     !                   g1(i+1,pid)))/12.))+
     !                Conc(jS,i+1,pid)*
     !      alf*((thO(i+1,pid)*Disp(i+1,pid)+thO(i,pid)*Disp(i,pid))/b1/2.-
     !           (vO(i,pid)+(2.-3.*wc(i,pid))*vO(i+1,pid))/6.+
     !           b1/12.*(g1(i,pid)+g1(i+1,pid)))+
     !              alf*(a1*(g0(i-1,pid)+2.*g0(i,pid))+b1*(2.*g0(i,pid)+
     !                   g0(i+1,pid)))/6.


          else
              B(i,pid)=epsi*(-(thN(i-1,pid)*Disp(i-1,pid)+thN(i,pid)*
     !        Disp(i,pid))/a1/2.-
     !        ((2.+3.*wc(i-1,pid))*vN(i-1,pid)+vN(i,pid))/6.-
     !        a1/12.*(g1(i-1,pid)+g1(i,pid)))
              D(i,pid)=dx/dt(pid)*thN(i,pid)*Retard(i,pid)+
     !        epsi*((thN(i-1,pid)*Disp(i-1,pid)+thN(i,pid)*Disp(i,pid))/a1/2.+
     !        (thN(i+1,pid)*Disp(i+1,pid)+thN(i,pid)*Disp(i,pid))/b1/2.+
     !        (vN(i+1,pid)+3.*(wc(i-1,pid)+wc(i,pid))*vN(i,pid)-vN(i-1,pid))/6.-
     !        (a1*(g1(i-1,pid)+3.*g1(i,pid))+b1*(3.*g1(i,pid)+g1(i+1,pid)))/12.)
              E(i,pid)=epsi*(-(thN(i+1,pid)*Disp(i+1,pid)+thN(i,pid)*
     !        Disp(i,pid))/b1/2.+(vN(i,pid)+(2.-3.*wc(i,pid))*vN(i+1,pid))/6.-
     !        b1/12.*(g1(i,pid)+g1(i+1,pid)))



              F(i,pid)=F(i,pid)+epsi*(a1*(g0(i-1,pid)+2.*g0(i,pid))+
     !        b1*(2.*g0(i,pid)+g0(i+1,pid)))/6.

          end if
          
            
          
          
11    continue
    
*    Upper boundary condition
      if(Level.eq.1) then
          FN=           Conc(jS,N-1,pid)*
     !    alf*((thO(N-1,pid)*Disp(N-1,pid)+thO(N,pid)*Disp(N,pid))/b1/2.+
     !    ((2.+3.*wc(N-1,pid))*vO(N-1,pid)+vO(N,pid))/6.+b1/12.*
     !    (g1(N-1,pid)+g1(N,pid)))+Conc(jS,N,pid)*
     !    (b1/2./dt(pid)*thO(N,pid)*Retard(N,pid)+
     !    alf*(-(thO(N-1,pid)*Disp(N-1,pid)+thO(N,pid)*Disp(N,pid))/b1/2.+
     !    (vO(N-1,pid)+(2.-3.*wc(N-1,pid))*vO(N,pid))/6.+
     !    b1/12.*(g1(N-1,pid)+3*g1(N,pid))))+
     !    alf*b1/6.*(g0(N-1,pid)+2.*g0(N,pid))

*       3. type BC
          if(kTopCh.le.0) then
              F(N,pid)=FN
              if(vO(N,pid).lt.0.) F(N,pid)=F(N,pid)-alf*vO(N,pid)*cTop(pid,js)
              if(kTopCh.eq.-2) then
                  M=MatNum(N,pid)
                  Tr=293.15
                  rr=8.314
                  jjj=(jS-1)*16
                  TT=(TempO(N,pid)+273.15-Tr)/rr/(TempO(N,pid)+273.15)/Tr
                  Dg=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
                  Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
                  F(N,pid)=F(N,pid)-alf*Dg/dSurf*Henry*Conc(jS,N,pid)+
     !            Dg/dSurf*cAtm
              end if
          end if
      else
          BN=epsi*(-(thN(N-1,pid)*Disp(N-1,pid)+thN(N,pid)*Disp(N,pid))/b1/2.-
     !    ((2.+3.*wc(N-1,pid))*vN(N-1,pid)+vN(N,pid))/6.-b1/12.*
     !    (g1(N-1,pid)+g1(N,pid)))
          DN=b1/2./dt(pid)*thN(N,pid)*Retard(N,pid)+
     !    epsi*((thN(N-1,pid)*Disp(N-1,pid)+thN(N,pid)*Disp(N,pid))/b1/2.-
     !    (vN(N-1,pid)+(2.-3.*wc(N-1,pid))*vN(N,pid))/6.-
     !    b1/12.*(g1(N-1,pid)+3.*g1(N,pid)))
          FE=epsi*b1/6.*(g0(N-1,pid)+2.*g0(N,pid))
          FN=FN+FE

*       1. type BC
          if(kTopCh.gt.0) then
              B(N,pid)=0.
              D(N,pid)=1.
              F(N,pid)=cTop(pid,jS)

*       3. type BC
          else
              B(N,pid)=BN
              D(N,pid)=DN
              F(N,pid)=F(N,pid)+FE
              if(vN(N,pid).lt.0.) F(N,pid)=F(N,pid)-epsi*vN(N,pid)*cTop(pid,jS)
              if(kTopCh.eq.-2) then
                  M=MatNum(N,pid)
                  Tr=293.15
                  rr=8.314
                  jjj=(jS-1)*16
                  TT=(TempN(N,pid)+273.15-Tr)/rr/(TempN(N,pid)+273.15)/Tr
                  Dg=ChPar(jjj+6,M)*exp(TDep(jjj+6)*TT)
                  Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
                  D(N,pid)=D(N,pid)+epsi*Dg/dSurf*Henry
              end if
          end if
      end if
      

      do 12 i=1,N
          M=MatNum(i,pid)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
              ThImob=ChPar(4,M)
              thN(i,pid)=thN(i,pid)+ThImob
              thO(i,pid)=thO(i,pid)+ThImob
          end if
          
     
12    continue




      return
      end


      real function rMD(NMatD,NSD,M,jS,jReact,DMoist,iReact,WDep,Theta,
     !iMoistDep)



*     Function expressing reaction rate dependence on the water content
*     ReacMin0  - relative minimum rate of reaction at low water contents
*     Theta0    - water content at which reaction rate start increasing
*     Theta1    - water content at which reaction rate stops increasing
*     Theta2    - water content at which reaction rate start decreasing  
*     Theta3    - water content at which reaction rate stops decreasing 
*     ReacMin1  - relative minimum rate of reaction at high water contents
*     If theta2=theta3=thetaS -> Anaerobic process
*     If theta0=theta1=0      -> Aerobic process
*     If theta2=0 -> no reduction

      double precision :: DMoist(NMatD,NSD,13,6),WDep(2+NMatD,NSD*9),theta
     !
      INTEGER imoistdep
      rMD=1.
      if(iMoistDep.eq.2) then
          if(jReact.eq.0) return
          ReacMin0=DMoist(M,jS,jReact,1)
          Theta0  =DMoist(M,jS,jReact,2)
          Theta1  =DMoist(M,jS,jReact,3)
          Theta2  =DMoist(M,jS,jReact,4)
          Theta3  =DMoist(M,jS,jReact,5)
          ReacMin1=DMoist(M,jS,jReact,6)
          if(abs(Theta2).lt.0.001) return
          if     (Theta.le.Theta0) then
              rMD=ReacMin0
          else if(Theta.le.Theta1) then
              rMD=ReacMin0+(Theta-Theta0)/(Theta1-Theta0)*(1.-ReacMin0)
          else if(Theta.le.Theta2) then
              rMD=1.
          else if(Theta.le.Theta3) then
              rMD=ReacMin1+(Theta-Theta3)/(Theta2-Theta3)*(1.-ReacMin1)
          else
              rMD=ReacMin1
          end if
      else if(iMoistDep.eq.1) then ! Walker's formula
          jjj=(jS-1)*9
          if(WDep(2+M,jjj+iReact).gt.Theta.and.WDep(2+M,jjj+iReact).gt.0.)
     !    rMD=(Theta/WDep(2+M,jjj+iReact))**WDep(1,jjj+iReact)
      end if

      return
      end


      subroutine Blocking(iPsi,SMax,psi,x,ss,Dc,SMax2)

      real Minf
      DOUBLE PRECISION X,SS,smax2,dc,ipsi,smax,psi
      psi=1.
      if     (iPsi.eq.1) then
          if(SMax.gt.0.) psi=1.-ss/SMax
      else if(iPsi.eq.2) then
          if(SMax.gt.0.) psi=max(ss**SMax,psi)
      else if(iPsi.eq.3) then
          Binf=1./SMax
          Sinf=.546
          Minf=Dc
          const=Sinf*Binf*ss
          if(ss.le.(0.8*SMax))
     !    psi=1.-(4.*const)+(3.08*const**2.)+(1.4069*const**3.)
          if(ss.gt.(0.8*SMax))
     !    psi=((1.-Binf*ss)**3.)/(2.*(Minf**2.)*(Binf**3.))
      else if(iPsi.eq.4) then
          if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
      else if(iPsi.eq.5) then
          if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
          if(SMax2.gt.0.) psi=psi*(1.-ss/SMax2)
      end if

      return
      end

*************************************************************************

*     Calculate the deposition coefficient for the bacteria transport,
*     All calculations within this subroutines are in meters and seconds
*     Conversions are needed

      subroutine Deposit(Ka1,Ka2,Dc1,Dp1,Alfa1,Alfa2,Theta,q,Temp,xConv,
     !tConv)

      real  mu,N_Pe,N_Lo,N_R,N_G
      double precision temp,theta,xconv,q,alfa1,alfa2,dc1,dp1,ka1,ka2,
     !tconv

*     Ka       - deposition coefficient (output) [1/T]
*     Dc       - diameter of the sand grains (m)
*     Dp       - diameter of the bacteria (0.95 microm) (m)
*     Alfa     - sticking efficiency (-)
*     Theta    - porosity (-)
*     q        - Darcy flux [L/T]
*     Temp     - Temperature in Celcius

      Dc=Dc1/xConv
      Dp=Dp1/xConv
      if(Dp.le.0.and.Dc.le.0.) then
          write(*,*) 'Both Dp and Dc are equal to zero !!!'
          write(*,*) 'Press Enter to continue'
          read(*,*)
          stop
      end if
      PI=3.1415                ! Ludolf's number
      mu=0.00093               ! fluid viscosity (Pa s)
      Bk=1.38048e-23           ! Boltzman constatnt (J/K)
      H=1.e-20                 ! Hamaker constant (J)
      g=9.81                   ! gravitational acceleration (m/s2)
      rop=1080.                ! bacterial density (kg/m3)
      rof=998.                 ! fluid density (kg/m3)
      Veloc=abs(q/xConv*tConv) ! absolute value of Darcy flux (converted to m/s)
      PVeloc=Veloc/Theta       ! pore velocity (converted to m/ )
      Dc=Dc1/xConv             ! conversion to m
      Dp=Dp1/xConv             ! conversion to m

      if(Veloc.gt.0.) then
          gamma=(1.-Theta)**(1./3.)
          As=2.*(1.-gamma**5)/(2.-3.*gamma+3.*gamma**5-2.*gamma**6) ! Correct.factor
          N_Pe=3.*PI*mu*Dp*Dc*Veloc/(Bk*(Temp+273.15)) ! Peclet number
          e_diff=4.*As**(1./3.)*N_Pe**(-2./3.)         ! removal by diffusion

          N_Lo=4.*H/(9.*PI*mu*Dp**2*Veloc)             ! London number
          N_R=Dp/Dc                                    ! Interception number
          e_inter=As*N_Lo**(1./8.)*N_R**(15./8.)       ! removal interception

          N_G=g*(rop-rof)*Dp**2/(18.*mu*Veloc)         ! Gravitation number
          e_grav=0.00338*As*N_G**1.2*N_R**(-0.4)       ! removal by gravitational
*                                                      sedimentation
      else
          e_diff =0.
          e_inter=0.
          e_grav =0.
      end if
      eta=e_diff+e_inter+e_grav                      ! single-collector efficiency

*     Original Filtration Theory
      Ka1=3.*(1.-Theta)/2./dc*eta*Alfa1*PVeloc
      Ka1=Ka1/tConv
      Ka2=3.*(1.-Theta)/2./dc*eta*Alfa2*PVeloc
      Ka2=Ka2/tConv

      return
      end

      !---------------------------------------------------------------------------
      subroutine Disper(pid,i,retard)

      USE GWFUNSFMODULE,ONLY: NMat,NSD,NumNP,npunsf,numnpd, 
     !i,dt,lTort,lArtD,lUpW,lMobIm,matnum,xKs,fExp,Frac,
     !iDualPor,Level,NLevel,ChPar,Disp,thSat,TauG,iTort,
     !ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,PeCr,lBact,
     !xNu,cMid,dSConc
      
      INTEGER,INTENT(IN) :: PID
      double precision retard(numnpd,npunsf)
      
      m=matnum(i,pid)
      if(lTort) then
          ThS=thSat(M)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact)
     !    ThS=max(thSat(M)-ThImob,0.001)
          if(              iDualPor.gt.0) ThS=    thSat(M)+ThImob
          if(iTort.eq.0) then
              TauW=ThW**(7./3.)/ThS**2
              TauG=ThG**(7./3.)/ThS**2
          else
              TauW=0.66*(ThW/ThS)**(8./3.)
              TauG=ThG**1.5/ThS
          end if
      else
          TauW=1.
          TauG=1.
      end if
      Disp(i,pid)=ChPar(2,M)*abs(v)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
      if(.not.lArtD.and..not.lUpW) then
          fi=0.
          if(cMid.gt.0.)
     !    fi=6.*ThW*ro*xKs*cMid**(fExp-1.)*
     !    (fExp/(1.+xNu*cMid**fExp)**2-1./(1.+xNu*cMid**fExp))
          DPom=amax1(dt(pid)/(6.*ThW*(ThW+ro*Frac*xKs*dSConc+ThG*Henry)+fi),0.)
          if(Level.ne.NLevel) then
              Disp(i,pid)=Disp(i,pid)+v*v*DPom
          else
              Disp(i,pid)=amax1(Disp(i,pid)-v*v*DPom,Disp(i,pid)/2.)
          end if
      end if
      if(lArtD) then
          DD=0.
          if(PeCr.ne.0.and.abs(v).gt.1.e-15) DD=v*v*dt(pid)/thW/thW/
     !    Retard(i,pid)/PeCr
          if(DD.gt.Disp(i,pid)) Disp(i,pid)=DD
      end if

      return
      end

*************************************************************************

      subroutine NEquil(pid,i,GamL,GamL1,GamLO,GamL1O,GamS,
     !GamS1,GamSO,GamS1O)

      use GWFUNSFMODULE, only: i,jS,NSD,NumNP,NMat,M,Conc,Sorb,Sorb2,
     !SorbN,SorbN2,SSorb,SSorb2,lMobIm,lLinear,lBact,GamS1Pi,xKsP,
     !Level,NLevel,dt,epsi,ro,xKs,xKsO,cc,SConc,ThImobO,iDualPor,
     !SConcO,SConcS,SConcOS,dSConcS,xMuL,xMuLO,xMuS,lNEquil,GamL1Pi,
     !xMuSO,dRetardS,Omega,OmegaO,rKa1,rKa1O,rKa2,psi1,psi1O,cNew,
     !rKa2O,rKd1,rKd1O,rKd2,rKd2O,ThW,ThWO,psi2,psi2O,DMobI,Frac,
     ! ThImob,SinkIm,FlMacro,SConcPS,lDualNEq,f_em,OmegaS,matnum

      INTEGER,INTENT(IN) :: PID
      double precision GamL,GamL1,GamLO,GamL1O,GamS,
     !GamS1,GamSO,GamS1O

      M=matnum(i,pid)
      FlMacro=0.
      if(iDualPor.gt.0) then            ! mobile-immobile model
          if(SinkIm(i,pid).gt.0) then
              FlMacro=SinkIm(i,pid)*Conc(jS,i,pid)
          else
              FlMacro=SinkIm(i,pid)*Sorb(jS,i,pid)
          end if
      end if
      SSorb =Sorb (jS,i,pid)
      SSorb2=Sorb2(jS,i,pid)
      if(Level.eq.NLevel) then

*       mobile-immobile model
      if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
          AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs*dSConcS
          dTheta=ThImob-ThImobO
          EMobI =ThImob*xMuL+(1.-Frac)*ro*xMuS-(1.-Frac)*ro*dRetardS
          EMobIO=ThImobO*xMuLO+(1.-Frac)*ro*xMuSO-(1.-Frac)*ro*dRetardS
          EMobI =EMobI +FlMacro
          EMobIO=EMobIO+FlMacro
          if(lNEquil) then
              cS=Sorb(jS-1,i,pid)*(ThImob*GamL1Pi+
     !        (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
              EMobI =EMobI +cS
              EMobIO=EMobIO+cS
          end if
          BMobI =ThImob *(GamL +GamL1 )+
     !    (1.-Frac)*ro*(GamS +GamS1 )*xKs *SConcS
          BMobIO=ThImobO*(GamLO+GamL1O)+
     !    (1.-Frac)*ro*(GamSO+GamS1O)*xKsO*SConcOS
          if(lLinear(jS)) then
              DMobI=  2.*AMobI+dt(pid)*(Omega +BMobI )+dTheta
              GMobI0=(2.*AMobI-dt(pid)*(OmegaO+BMobIO)-dTheta)/DMobI
              Sorb(jS,i,pid)=Sorb(jS,i,pid)*GMobI0+
     !        dt(pid)*(OmegaO*Conc(jS,i,pid)+EMobI+EMobIO)/DMobI
              SSorb=Sorb(jS,i,pid)
          else
              SorbN(i,pid)=Sorb(jS,i,pid)+dt(pid)/AMobI*( epsi *(Omega *
     !        (cNew(i,PID)   -SorbN(i,pid))  -BMobI *SorbN(i,PID)  +EMobI)+
     !        (1.-epsi)*(OmegaO*
     !        (Conc(jS,i,pid)-Sorb(jS,i,pid))-BMobIO*Sorb(jS,i,pid)+EMobIO))
              SSorb=SorbN(i,pid)
          end if
          if(lDualNEq) then
              cS=0.
              if(lNEquil) cS=GamS1Pi*Sorb2(jS-1,i,pid)*dt(pid)*2.
              if(lLinear(jS)) then
                  Sorb2(jS,i,pid)=((2.-(OmegaS+GamSO+GamS1O)*dt(pid))*
     !            Sorb2(jS,i,pid)+
     !            dt(pid)*Frac*(1.-f_em)*OmegaS*xKsO*Conc(jS,i,pid)+
     !            dt(pid)*(1.-f_em)*(xMuSO+xMuS)+cS)/
     !            (2.+dt(pid)*(OmegaS+GamS+GamS1))
                  SSorb2=Sorb2(jS,i,pid)
              else
                  SorbN2(i,pid)=Sorb2(jS,i,pid)+dt(pid)*
     !            (epsi* (OmegaS*(Frac*(1.-f_em)*SConc *xKs *cc-SorbN2(i,pid))-
     !            (GamS+GamS1)*SorbN2(i,pid)+(1.-f_em)*xMuS)+
     !            (1.-epsi)*(OmegaS*(Frac*(1.-f_em)*SConcO*xKsO*
     !            Conc(jS,i,pid)-SSorb2)-
     !            (GamSO+GamS1O)*SSorb2+(1.-f_em)*xMuSO))
                  SSorb2=SorbN2(i,pid)
              end if
          end if

*       two-site sorption model
      else if(.not.lBact) then
          cS=0.
          if(lNEquil) cS=GamS1Pi*Sorb(jS-1,i,pid)*dt(pid)*2.
          if(lLinear(jS)) then
              Sorb(jS,i,pid)=((2.-(OmegaO+GamSO+GamS1O)*dt(pid))*Sorb(jS,i,pid)+
     !        dt(pid)*(1.-Frac)*OmegaO*xKsO*Conc(jS,i,pid)+
     !        dt(pid)*(1.-Frac)*(xMuSO+xMuS)+cS)/
     !        (2.+dt(pid)*(Omega+GamS+GamS1))
              SSorb=Sorb(jS,i,pid)
          else
              SorbN(i,pid)=Sorb(jS,i,pid)+dt(pid)*
     !        (epsi* (Omega* ((1.-Frac)*SConc *xKs *cc-SorbN(i,pid))-
     !        (GamS+GamS1)*SorbN(i,pid)+(1.-Frac)*xMuS)+
     !        (1.-epsi)*(OmegaO*((1.-Frac)*SConcO*xKsO*Conc(jS,i,pid)-SSorb)-
     !        (GamSO+GamS1O)*SSorb+(1.-Frac)*xMuSO))
              SSorb=SorbN(i,pid)
          end if

*       filtration model
      else if(lBact) then
          if(lLinear(jS)) then
              Sorb(jS,i,pid)=((2.-dt(pid)*(rKd1O+GamSO+GamS1O))*Sorb(jS,i,pid)+
     !        dt(pid)*rKa1O*ThW*Conc(jS,i,pid)/ro)/
     !        (2.+dt(pid)*(rKd1+GamS+GamS1))
              SSorb=Sorb(jS,i,pid)
              Sorb2(jS,i,pid)=((2.-dt(pid)*(rKd2O+GamSO+GamS1O))*
     !        Sorb2(jS,i,pid)+
     !        dt(pid)*rKa2O*ThW*Conc(jS,i,pid)/ro)/
     !        (2.+dt(pid)*(rKd2+GamS+GamS1))
              SSorb2=Sorb2(jS,i,pid)
          else
              SorbN(i,pid)=Sorb(jS,i,pid)+dt(pid)*
     !        (epsi*    (rKa1*ThW/ro*psi1*cc-
     !        (rKd1+GamS+GamS1)*SorbN(i,pid))+
     !        (1.-epsi)*(rKa1O*ThWO/ro*psi1O*Conc(jS,i,pid)-
     !        (rKd1O+GamSO+GamS1O)*Sorb(jS,i,pid)))
              SSorb=SorbN(i,pid)
              SorbN2(i,pid)=Sorb2(jS,i,pid)+dt(pid)*
     !        (epsi*    (rKa2*ThW/ro*psi2*cc-
     !        (rKd2+GamS+GamS1)*SorbN2(i,pid))+
     !        (1.-epsi)*(rKa2O*ThWO/ro*psi2O*Conc(jS,i,pid)-
     !        (rKd2O+GamSO+GamS1O)*Sorb2(jS,i,pid)))
              SSorb2=SorbN2(i,pid)
          end if
      end if
      end if

      return
      end

*************************************************************************

*     Calculate the maximum local Peclet and Courant numbers

      subroutine PeCour(pid,i,j)

      use GWFUNSFMODULE,only: Level,NLevel,lUpW,lArtD,dt,x,v,wc,ThW,
     !vj,Thj,Disp,Retard,Peclet,Courant,CourMax,PeCr,Iter,epsi,dtMaxC,
     !numnp
      
      INTEGER,INTENT(IN) :: PID

      TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))

      if(i.ne.numnp(pid)) then
          dx=x(j,pid)-x(i,pid)
          vv=0.
          if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv=(abs(v)/ThW+abs(vj)/Thj)/2.
          vv1=0.
          if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv1=(v/ThW+vj/Thj)/2.
          DD=(Disp(i,pid)+Disp(j,pid))/2.
          if(Level.eq.NLevel) then
              Pec=99999.
              dtMax=1.e+30
c          vMax=amax1(abs(v)/ThW,abs(vj)/Thj)
              vMax=(abs(v)+abs(vj))/(ThW+Thj)
              RMin=amin1(Retard(i,pid),Retard(j,pid))
              if(DD.gt.0.) Pec=abs(vv)*dx/DD
              Cour=vMax*dt(pid)/dx/RMin
              Peclet=amax1(Peclet,Pec)
              Courant=amax1(Courant,Cour)
              Cour1=CourMax
              if(.not.lUpW.and..not.lArtD) then
                  if(Pec.ne.99999.) Cour1=amin1(1.,PeCr/amax1(0.5,Pec))
              end if
              if(epsi.lt.1..and.vMax.gt.1.e-20) dtMax=Cour1*dx*RMin/vMax
*         the von Neumann time step limit
c          RThE=(ThW+thj)/2.*RMin
c          if(abs(DD).gt.1.e-20)dtMax=amin1(dtMax,10.*RThE*dx*dx/2./DD)
              dtMaxC=amin1(dtMaxC,dtMax)

*       Calculate upstream weighting factors
          else if(lUpW.and.Iter(pid).eq.1) then
              Pe2=11.
              if(DD.gt.0.) Pe2=dx*vv1/DD/2.
              if(abs(vv).lt.1.e-30) then
                  wc(i,pid)=0.
              else if(abs(Pe2).gt.10.) then
                  if(vv1.gt.0.) wc(i,pid)=1.
                  if(vv1.lt.0.) wc(i,pid)=-1
              else
                  wc(i,pid)=1./TanH(Pe2)-1./Pe2
                  wc(i,pid)=amin1( 1.,wc(i,pid))
                  wc(i,pid)=amax1(-1.0,wc(i,pid))
              end if
          end if
      end if

      return
      end
*      *************************************************************************

*     Calculate mass-transfer fluxes at the end of the time interval

      subroutine MassTran(pid,jS,Temp,epsi,theta,Veloc)

      use GWFUNSFMODULE,only: NS,NSD,MatNum,lMobIm,lEquil,ChPar,
     !TDep,Sorb,Conc,NMat,x,cvCh0,cvCh1,cvChR,cvChIm,q0,q1,sSink,
     !lFiltr,iDualPor,SinkIm,xConv,tConv,lBact,Sorb2,smax1i,smaxij,
     !lDualNEq,STrans,lLinear,numnp,aa,alfa1,alfa2,dc,dp,numnpd,
     !iPsi1,ipsi2,psi1j,psi2i,psi2j,rka1,rka2,psi1i,smax1j,Smax2i,
     !smax2j,thetai,NPUNSF
      
      INTEGER,INTENT(IN) :: PID
      DOuble precision :: Temp(numnpd,NPUNSF),Veloc(numnpd),theta(NumNPd,npunsf)
      INTEGER:: N

      N=numnp(pid)
      Tr=293.15
      rr=8.314
      jjj=(jS-1)*16
      do 11 i=1,N
          Mi=matnum(i,pid)
          TTi=(Temp(i,pid)+273.15-Tr)/rr/(Temp(i,pid)+273.15)/Tr
          if(i.eq.N) goto 10
          j=i+1
          dx=x(j,pid)-x(i,pid)
          cvCh0(jS,pid)=cvCh0(jS,pid)+epsi*dx*(q0(i,pid)+q0(j,pid))/2.
          cvCh1(jS,pid)=cvCh1(jS,pid)+epsi*dx*(q1(i,pid)+q1(j,pid))/2.
          cvChR(jS,pid)=cvChR(jS,pid)+epsi*dx*(sSink(i,pid)+sSink(j,pid))/2.
          if(.not.lEquil) then
              Mj=MatNum(j,pid)
              TTj=(Temp(j,pid)+273.15-Tr)/rr/(Temp(j,pid)+273.15)/Tr
              Omegai=ChPar(jjj+20,Mi)*exp(TDep(jjj+20)*TTi)
              Omegaj=ChPar(jjj+20,Mj)*exp(TDep(jjj+20)*TTj)

*         mobile-immobile model
              if((lMobIm(Mi).or.iDualPor.gt.0).and..not.lBact) then
                  FlMacroi=0.
                  FlMacroj=0.
                  if(iDualPor.gt.0) then
                      if(SinkIm(i,pid).gt.0) then
                          FlMacroi=SinkIm(i,pid)*Conc(jS,i,pid)
                      else
                          FlMacroi=SinkIm(i,pid)*Sorb(jS,i,pid)
                      end if
                      if(SinkIm(j,pid).gt.0) then
                          FlMacroj=SinkIm(j,pid)*Conc(jS,j,pid)
                      else
                          FlMacroj=SinkIm(j,pid)*Sorb(jS,j,pid)
                      end if
                  end if
                  cvChIm(jS,pid)=cvChIm(jS,pid)+epsi*dx/2.*
     !            (Omegai*(Conc(jS,i,pid)-Sorb(jS,i,pid))+
     !            Omegaj*(Conc(jS,j,pid)-Sorb(jS,j,pid))+
     !            FlMacroi+FlMacroj)
                  if(jS.eq.1.and..not.lLinear(jS))
     !            STrans(i,pid)=epsi*(Omegai*(Conc(jS,i,pid)-Sorb(jS,i,pid))+
     !            FlMacroi)
                  if(lDualNEq) then
                      roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
                      Fraci =ChPar(3,     Mi)*exp(TDep(3)     *TTi)
                      xKsi  =ChPar(jjj+ 7,Mi)*exp(TDep(jjj+ 7)*TTi)
                      xNui  =ChPar(jjj+ 8,Mi)*exp(TDep(jjj+ 8)*TTi)
                      fExpi =ChPar(jjj+ 9,Mi) !*exp(TDep(jjj+ 9)*TTi)
                      roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
                      Fracj =ChPar(3,     Mj)*exp(TDep(3)     *TTj)
                      xKsj  =ChPar(jjj+ 7,Mj)*exp(TDep(jjj+ 7)*TTj)
                      xNuj  =ChPar(jjj+ 8,Mj)*exp(TDep(jjj+ 8)*TTj)
                      fExpj =ChPar(jjj+ 9,Mj) !*exp(TDep(jjj+ 9)*TTj)
                      f_emi =ChPar(jjj+13,Mi)*exp(TDep(jjj+13)*TTi)
                      f_emj =ChPar(jjj+13,Mj)*exp(TDep(jjj+13)*TTj)
                      OmegaSi=ChPar(jjj+16,Mi)*exp(TDep(jjj+16)*TTi)
                      OmegaSj=ChPar(jjj+16,Mj)*exp(TDep(jjj+16)*TTj)
                      cci   =Conc(jS,i,pid)
                      ccj   =Conc(jS,j,pid)
                      SorbEi=0.
                      SorbEj=0.
                      if(cci.gt.0.) SorbEi=
     !                Fraci*(1.-f_emi)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
                      if(ccj.gt.0.) SorbEj=
     !                Fracj*(1.-f_emj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
                      cvChIm(jS,pid)=cvChIm(jS,pid)+epsi*dx/2.*
     !                (roi*OmegaSi*(SorbEi-Sorb2(jS,i,pid))+
     !                roj*OmegaSj*(SorbEj-Sorb2(jS,j,pid)))
                      if(jS.eq.1.and..not.lLinear(jS))
     !                STrans(i,pid)=STrans(i,pid)+epsi*roi*OmegaSi*
     !                (SorbEi-Sorb2(jS,i,pid))
                  end if

*         two-site sorption model
              else if(.not.lBact) then
                  roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
                  Fraci =ChPar(3,     Mi)*exp(TDep(3)     *TTi)
                  xKsi  =ChPar(jjj+ 7,Mi)*exp(TDep(jjj+ 7)*TTi)
                  xNui  =ChPar(jjj+ 8,Mi)*exp(TDep(jjj+ 8)*TTi)
                  fExpi =ChPar(jjj+ 9,Mi) !*exp(TDep(jjj+ 9)*TTi)
                  roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
                  Fracj =ChPar(3,     Mj)*exp(TDep(3)     *TTj)
                  xKsj  =ChPar(jjj+ 7,Mj)*exp(TDep(jjj+ 7)*TTj)
                  xNuj  =ChPar(jjj+ 8,Mj)*exp(TDep(jjj+ 8)*TTj)
                  fExpj =ChPar(jjj+ 9,Mj) !*exp(TDep(jjj+ 9)*TTj)
                  cci   =Conc(jS,i,pid)
                  ccj   =Conc(jS,j,pid)
                  SorbEi=0.
                  SorbEj=0.
                  if(cci.gt.0.)
     !            SorbEi=(1.-Fraci)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
                  if(ccj.gt.0.)
     !            SorbEj=(1.-Fracj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
                  cvChIm(jS,pid)=cvChIm(jS,pid)+epsi*dx/2.*
     !            (roi*Omegai*(SorbEi-Sorb(jS,i,pid))+
     !            roj*Omegaj*(SorbEj-Sorb(jS,j,pid)))
                  if(jS.eq.1)
     !            STrans(i,pid)=epsi*roi*Omegai*(SorbEi-Sorb(jS,i,pid))

*         filtration model
              else if(lBact) then
                  roi   =ChPar(1,     Mi)*exp(TDep(1)     *TTi)
                  roj   =ChPar(1,     Mj)*exp(TDep(1)     *TTj)
                  SMax1i=ChPar(jjj+18,Mi)*exp(TDep(jjj+18)*TTi)
                  SMax1j=ChPar(jjj+18,Mj)*exp(TDep(jjj+18)*TTj)
                  rKa1i =ChPar(jjj+19,Mi)*exp(TDep(jjj+19)*TTi)
                  rKa1j =ChPar(jjj+19,Mj)*exp(TDep(jjj+19)*TTj)
                  rKd1i =ChPar(jjj+20,Mi)*exp(TDep(jjj+20)*TTi)
                  rKd1j =ChPar(jjj+20,Mj)*exp(TDep(jjj+20)*TTj)
                  SMax2i=ChPar(jjj+15,Mi)*exp(TDep(jjj+15)*TTi)
                  SMax2j=ChPar(jjj+15,Mj)*exp(TDep(jjj+15)*TTj)
                  rKa2i =ChPar(jjj+16,Mi)*exp(TDep(jjj+16)*TTi)
                  rKa2j =ChPar(jjj+16,Mj)*exp(TDep(jjj+16)*TTj)
                  rKd2i =ChPar(jjj+17,Mi)*exp(TDep(jjj+17)*TTi)
                  rKd2j =ChPar(jjj+17,Mj)*exp(TDep(jjj+17)*TTj)
                  ThImobi=ChPar(4,Mi)
                  ThImobj=ChPar(4,Mj)
                  Thetai=theta(i,pid)-ThImobi
                  Thetaj=theta(j,pid)-ThImobj
                  iPsi1=0
                  iPsi2=0
                  if(.not.lFiltr) iPsi2=int(ChPar(jjj+13,Mi))
                  if(.not.lFiltr) iPsi1=int(ChPar(jjj+14,Mj))
                  if(iPsi1.eq.0.and.SMax1i.gt.0.) iPsi1=1
                  if(iPsi2.eq.0.and.SMax2i.gt.0.) iPsi2=1
                  if(iPsi1.ge.3.or.iPsi2.ge.3)
     !            Dc=ChPar(jjj+6,Mi)*exp(TDep(jjj+6)*TTi)
                  if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(jjj+15,Mi)
                  psi1i=1.
                  psi1j=1.
                  psi2i=1.
                  psi2j=1.
                  if(iPsi1.gt.0) then
                      call Blocking(iPsi1,SMax1i,psi1i,x(i,pid),
     !                Sorb(jS,i,pid),Dc,aa)
                      call Blocking(iPsi1,SMax1j,psi1j,x(j,pid),
     !                Sorb(jS,j,pid),Dc,aa)
                  end if
                  if(iPsi2.gt.0) then
                      call Blocking(iPsi2,SMax2i,psi2i,x(i,pid),
     !                Sorb2(jS,i,pid),Dc,aa)
                      call Blocking(iPsi2,SMax2j,psi2j,x(j,pid),
     !                Sorb2(jS,j,pid),Dc,aa)
                  end if
                  if(lFiltr) then
                      Dc=ChPar(jjj+13,Mi)*exp(TDep(jjj+13)*TTi)
                      Dp=ChPar(jjj+14,Mi)*exp(TDep(jjj+14)*TTi)
                      Alfa1=rKa1i
                      Alfa2=rKa2i
                      call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,Thetai,
     !                Veloc(i),Temp(i,pid),xConv,tConv)
                      rKa1i=rKa1
                      rKa1j=rKa1
                      rKa2i=rKa2
                      rKa2j=rKa2
                  end if
                  cvChIm(jS,pid)=cvChIm(jS,pid)+epsi*dx/2.*
     !            (Conc(jS,i,pid)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)+
     !            Conc(jS,j,pid)*Thetaj*(psi1j*rKa1j+psi2j*rKa2j)-
     !            roi*(Sorb(jS,i,pid)*rKd1i+Sorb2(jS,i,pid)*rKd2i)-
     !            roj*(Sorb(jS,j,pid)*rKd1j+Sorb2(jS,j,pid)*rKd2j))
                  if(jS.eq.1) STrans(i,pid)=epsi*
     !            (Conc(jS,i,pid)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)-
     !            roi*(Sorb(jS,i,pid)*rKd1i+Sorb2(jS,i,pid)*rKd2i))
              end if
          end if
10        continue
11    continue

      return
      end

*************************************************************************
*     Solve matrix equation

      subroutine BanSol(pid,A,B,C,F)
      use GWFUNSFMODULE,only: numnp,numnpd,npunsf
      INTEGER,INTENT(IN) :: PID
      double precision A(Numnpd,npunsf),B(Numnpd,npunsf),
     !C(Numnpd,npunsf),F(Numnpd,npunsf)
      INTEGER i,j

      N=numnp(pid)
      do 11 i=2,N
          B(i,pid)=B(i,pid)-A(i,pid)*C(i-1,pid)/B(i-1,pid)
          F(i,pid)=F(i,pid)-A(i,pid)*F(i-1,pid)/B(i-1,pid)
11    continue
      F(N,pid)=F(N,pid)/B(N,pid)

      do 12 i=2,N
          j=N-i+1
          F(j,pid)=(F(j,pid)-C(j,pid)*F(j+1,pid))/B(j,pid)

12    continue


      return
      end

*************************************************************************

      subroutine SorbConc(pid,js,Temp,ThW, Veloc,ThIm)



      use GWFUNSFMODULE,only:NSD,MatNum,lMobIm,ChPar,TDep,Sorb,
     !Conc,dt,NMat,lBact,Sorb2,lFiltr,
     !iDualPor,ThOIm,SinkIm,STrans,iMoistDep,numnp,
     !NMatD,DMoist,WDep,xConv,tConv,lDualNEq,numnpd,npunsf,
     !alfa1, alfa2,dc,dp,rka1,rka2,theta,rkd1,rkd2,thimob

      INTEGER,INTENT(IN) :: PID
      INTEGER:: N
     !! logical lmobim(NMatd),lbact,lfiltr,ldualneq
     !! INTEGER  idualpor
      double precision thetas
      double precision  Temp(numnpd,NPUNSF),ThW(numnpd,NPUNSF),
     !Veloc(numnpd),ThIm(numnpd,NPUNSF)
      Tr=293.15
      rr=8.314
      jjj=(jS-1)*16
      N=numnp(pid)
      do 11 i=1,N
          M=MatNum(i,pid)
          TT=(Temp(i,pid)+273.15-Tr)/rr/(Temp(i,pid)+273.15)/Tr
          Frac  =ChPar(3,     M)*exp(TDep(3)     *TT)
          xKs   =ChPar(jjj+ 7,M)*exp(TDep(jjj+ 7)*TT)
          f1=1.
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThW(i,pid),iMoistDep)
          GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
          if(iMoistDep.gt.0)
     !    f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThW(i,pid),iMoistDep)
          GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
          Omega =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)

          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
              ro    =ChPar(1,     M)*exp(TDep(1 )     *TT)
              if(iDualPor.eq.0) then
                  ThImob=ChPar(4,M)*exp(TDep(4 )*TT)
                  ThImobO=ThImob
              else if(iDualPor.gt.0) then
                  ThImob =ThIm(i,pid)
                  ThImobO=ThOIm(i,pid)
              end if
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,10,DMoist,1,WDep,ThImob,iMoistDep)
              GamL  =ChPar(jjj+11,M)*exp(TDep(jjj+11)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,12,DMoist,4,WDep,ThImob,iMoistDep)
              GamL1 =ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,11,DMoist,2,WDep,ThImob,iMoistDep)
              GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)*f1
              if(iMoistDep.gt.0)
     !        f1=rMD(NMatD,NSD,M,jS,13,DMoist,5,WDep,ThImob,iMoistDep)
              GamS1 =ChPar(jjj+15,M)*exp(TDep(jjj+15)*TT)*f1
              dTheta=ThImob-ThImobO
              AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs
              BMobI=ThImob*(GamL+GamL1)+(GamS+GamS1)*ro*(1.-Frac)*xKs
              DMobI=2.*AMobI+dt(pid)*(Omega+BMobI)+dTheta
              Sorb(jS,i,pid)=Sorb(jS,i,pid)+dt(pid)*Omega*Conc(jS,i,pid)/DMobI
              FlMacro=0.
              if(iDualPor.gt.0) then
                  if(SinkIm(i,pid).gt.0) then
                      FlMacro=SinkIm(i,pid)*Conc(jS,i,pid)
                  else
                      FlMacro=SinkIm(i,pid)*Sorb(jS,i,pid)
                  end if
              end if
              if(jS.eq.1) STrans(i,pid)=Omega*(Conc(jS,i,pid)-
     !        Sorb(jS,i,pid))+FlMacro
              if(lDualNEq) then
                  f_em  =ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
                  OmegaS=ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
                  Sorb2(jS,i,pid)=
     !            Sorb2(jS,i,pid)+dt(pid)*OmegaS*Frac*(1.-f_em)*xKs*
     !            Conc(jS,i,pid)/
     !            (2.+dt(pid)*(OmegaS+GamS+GamS1))
              end if

          else if(.not.lBact) then                 ! two-site sorption model
              Sorb(jS,i,pid)=Sorb(jS,i,pid)+dt(pid)*Omega*(1.-Frac)*xKs*
     !        Conc(jS,i,PID)/(2.+dt(pid)*(Omega+GamS+GamS1))

          else if(lBact) then                      ! filtration model
              ro    =ChPar(1,     M)*exp(TDep(1)     *TT)
              ThImob=ChPar(4,M)
              Thetas=ThW(i,pid)-ThImob
              GamS  =ChPar(jjj+12,M)*exp(TDep(jjj+12)*TT)
              GamS1 =0.
              rKa1  =ChPar(jjj+19,M)*exp(TDep(jjj+19)*TT)
              rKd1  =ChPar(jjj+20,M)*exp(TDep(jjj+20)*TT)
              rKa2  =ChPar(jjj+16,M)*exp(TDep(jjj+16)*TT)
              rKd2  =ChPar(jjj+17,M)*exp(TDep(jjj+17)*TT)
              if(lFiltr) then
                  Dc=ChPar(jjj+13,M)*exp(TDep(jjj+13)*TT)
                  Dp=ChPar(jjj+14,M)*exp(TDep(jjj+14)*TT)
                  Alfa1=rKa1
                  Alfa2=rKa2
                  call Deposit(rKa1,rKa2,Dc,Dp,Alfa1,Alfa2,Thetas,Veloc(i),
     !            Temp(i,pid),xConv,tConv)
              end if
              Sorb(jS,i,pid) =Sorb(jS,i,pid) +dt(pid)*rKa1*Thetas*
     !        Conc(jS,i,pid)/ro/(2.+dt(pid)*(rKd1+GamS+GamS1))
              Sorb2(jS,i,pid)=Sorb2(jS,i,pid)+dt(pid)*rKa2*Thetas*
     !        Conc(jS,i,pid)/ro/(2.+dt(pid)*(rKd2+GamS+GamS1))
          end if
11    continue

      return
      end

*     Calculate flux concentration for the first solute

      subroutine FluxConc(pid,x,v,theta,Temp,ConcF,ThIm,jS)

      use GWFUNSFMODULE,only: TDep,thsat, conc,chpar,numnp,numnpd,matnum,
     ! npunsf,nmatd,nsd,nmat
      INTEGER,INTENT(IN) :: PID
      logical lTort,lMobIm(NMat)
      Double precision x(NumNPd,npunsf),v(NumNPd),theta(NumNPd,npunsf),
     !ThIm(NumNPd),Temp(NumNPd,npunsf),ConcF(NumNPd,npunsf)
      

      jjj=(jS-1)*16
      Tr=293.15
      Rr=8.314

      do 11 i=1,NumNP(pid)
          M=MatNum(i,pid)
          ThW=Theta(i,pid)
          ThG=amax1(0.,thSat(M)-ThW)
          if(lMobIm(M)) then
              if(iDualPor.eq.0) then
                  ThImob=ChPar(4,M)
                  ThW=max(ThW-ThImob,0.001)
              else if(iDualPor.gt.0) then
                  ThImob=ThIm(i)
              end if
          end if
          TT=(Temp(i,pid)+273.15-Tr)/RR/(Temp(i,pid)+273.15)/Tr
          Dw   =ChPar(jjj+ 5,M)*exp(TDep(jjj+ 5)*TT)
          Dg   =ChPar(jjj+ 6,M)*exp(TDep(jjj+ 6)*TT)
          Henry=ChPar(jjj+10,M)*exp(TDep(jjj+10)*TT)
          if(lTort) then
              ThS=thSat(M)
              if(lMobIm(M).and.iDualPor.eq.0) ThS=max(thSat(M)-ThImob,0.001)
              if(              iDualPor.gt.0) ThS=    thSat(M)+ThImob
              TauW=ThW**(7./3.)/ThS**2
              TauG=ThG**(7./3.)/ThS**2
          else
              TauW=1.
              TauG=1.
          end if
          qW=v(i)
          Disp=ChPar(2,M)*abs(qW)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
          cGrad=0.
          if(i.eq.1) then
              cGrad=(Conc(jS,i+1,pid)-Conc(jS,i,pid))/(x(i+1,pid)-x(i,pid))
          else if(i.eq.NumNP(pid)) then
              cGrad=(Conc(jS,i,pid)-Conc(jS,i-1,pid))/(x(i,pid)-x(i-1,pid))
          else
              cGrad=(Conc(jS,i+1,pid)-Conc(jS,i-1,pid))/(x(i+1,pid)-x(i-1,pid))
          end if
          ConcF(i,pid)=Conc(jS,i,pid)
          if(qW.ne.0.) ConcF(i,pid)=Conc(jS,i,pid)-Disp*ThW/qW*cGrad
11    continue

      return
      end

**********************************************************************

*     Subroutine calculating root solute uptake with and without compensation

      subroutine SetSSnk(pid,SinkS)

      use GWFUNSFMODULE,only:x,t,Beta,sink,Conc,numnp,numnpd,nmatd, nmat,nsd,
     1npunsf,ns,crootmax
      INTEGER,INTENT(IN) :: PID

      logical lActRSU,lLast
      double precision SinkS(numnpd,npunsf)
      N=numnp(pid)
*     Inputs:
*     SPot      - potential root solute uptake
*     OmegaS    - solute stress index
*     rKM       - Michaelis-Menten constant
*     lActRSU   - consider active root solute uptake
*     cRootMax  - maximum concentration for the passive solute uptake

*     From Water Flow
*     Sink(i)   - Root water uptake
*     OmegaW    - ratio of actual and potential transpiration

*     SPUptake  - passive root solute uptake (step 1)
*     SAUptakeP - potential active solute uptake (step 1)
*     SAUptakeA - uncompensated actual active solute uptake (step 2)
*     SAUptakeA - compensated actual active solute uptake (step 3)
*     SinkS(i)  - local active solute uptake

*     Initialization
      Compen=1.
      nStep=1
      if(lActRSU)                  nStep=2
      if(lActRSU.and.OmegaS.lt.1.) nStep=3
*     step 1: Passive uptake
*     step 2: Active uptake without compensation
*     step 3: Active uptake with compensation
      lLast=.false.                 ! Active uptake only for the last solute
      if(lLast.and.jS.lt.NS) nStep=1
      Omega=0.
      SPUptake=0.
      do 10 i=1,N
          SinkS(i,pid)=0.
10    continue

      do 12 iStep=1,nStep
          SAUptakeA=0.
          do 11 i=1,N
              if(Beta(i,pid).gt.0.) then
                  if(i.eq.N) then
                      dxM=(x(i,pid)-x(i-1,pid))/2.
                  else if(i.eq.1) then
                      dxM=(x(i,pid)-x(i+1,pid))/2.
                  else
                      dxM=(x(i+1,pid)-x(i-1,pid))/2.
                  end if
                  cc=amax1(Conc(jS,i,pid)-cMin,0.)
                  if(iStep.eq.1) then
                      SinkS(i,pid)=Sink(i,pid)*amax1(amin1(Conc(jS,i,pid)
     !                ,cRootMax(pid,js)),0.)
                      SPUptake=SPUptake+SinkS(i,pid)*dxM
*             This is needed only for the last node, but that node may not have beta
                      SAUptakeP=amax1(SPot*OmegaW-SPUptake,0.)
                  else if(iStep.eq.2) then
                      AUptakeA=cc/(rKM+cc)*Beta(i,pid)*SAUptakeP
                      Omega=Omega+AUptakeA*dxM
                      if(nStep.eq.2) SinkS(i,pid)=SinkS(i,pid)+AUptakeA
*             This is needed only for the last node, but that node may not have beta
                      SAUptakeA =Omega
                      SAUptakeAN=Omega
                      if(SAUptakeP.gt.0.) Omega1=Omega/SAUptakeP
                  else if(iStep.eq.3) then
*             This is needed only for the first node, but that node may not have beta
                      if(Omega1.lt.OmegaS.and.Omega1.gt.0.) Compen=OmegaS
                      if(Omega1.ge.OmegaS)                  Compen=Omega1
                      if(Compen.gt.0.) AUptakeA=cc/(rKM+cc)*
     !                Beta(i,pid)*SAUptakeP/Compen
                      SinkS(i,pid)=SinkS(i,pid)+AUptakeA
                      SAUptakeA=SAUptakeA+AUptakeA*dxM
                  end if
              else
                  SinkS(i,pid)=0.
              end if
11        continue
          if(iStep.eq.nStep.and.jS.eq.NS)
     !    write(78,100) t,SPUptake,SAUptakeP,SAUptakeA,SAUptakeAN ! the last is uncompensated
12    continue
      return

100   format(3x,e14.7,1x,4e12.4)
      end

************************************************************************