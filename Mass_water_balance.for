
      subroutine gwf2unsf1SubReg(pid,ThN,ThO)
      !I,hhNew,ThNew,ThOld,t,dt,PLevel, cPrevO
      USE GWFUNSFMODULE, ONLY:ths,sink,tlevel,ns,cmid,tpulse,nsd,llinear,
     !retard,wc,thsat,tdep,nsd, WCUMT,WCUMA,watin,wvoli,area, conc,ccuma,
     !ccumt,con,chpar, tdep, sorb,sorb2,laynum,cprevo, idualpor,cvoli,
     !solin,lLinear ,ns,ths,con,thnewim,tholdim,cprev,x,nmat,nmatd,
     !numnpd,npunsf,numnp,matnum,plevel,hhnew,dt,t,lmobim,lWat,lChem,
     !lEquil,lPrint,
     !lVapor, lBact, lDensity,lDualNEq,js,volume,xmasschange

      INTEGER,INTENT(IN) :: PID
      logical lcentrif,ltemp,lwtdep
      integer   i

      dimension cMean(11,10),  ConSub(11,10),
     !ConVolIm2(11), hMean(10),
     1TMean(10),SubVol(10),SubCha(10),ConVol(11),
     1cTot(11),SubT(10) ,ConVolIm(11),
     1ConSubIm(11,10),cMeanIm(11,10),cTotIm(11),ConSubIm2(11,10)
      double precision:: dxn, ThN(Numnpd,npunsf),WW,
     !ThO(Numnpd,npunsf) ,wbalt,temp(Numnpd),cnewi,vnewi,voldi,cEl,DELTW

      N=numnp(pid)

      fRE=1
      Grav=1
      ATot=0.
      Tr=293.15
      R=8.314
      cosalf=1
      lwat(pid)=1
      temp=20
       DeltW=0.
       VNewi=0.
       VOLDI=0.
       VOLUME=0.
       
      do 15 i=N-1,1,-1
          j=i+1
          cEl=0.
          Mi=MatNum(i,pid )
          Mj=MatNum(j,pid)
          Lay=LayNum(i,pid)
          dx=x(j,pid)-x(i,pid)
          Area(Lay)=Area(Lay)+dx
          ATot=ATot+dx
          TT=(Temp(i)+Temp(j))/2.+273.15
          hE=(hhNew(i,pid)+hhNew(j,pid))/2.
          VNewi=dx*(ThN(i,pid)+ThN(j,pid))/2.
          VOldi=dx*(ThO(i,pid)+ThO(j,pid))/2.
          Volume=Volume+VNewi
          Change=Change+(VNewi-VOldi)/dt(pid)
          SubCha(Lay)=SubCha(Lay)+(VNewi-VOldi)/dt(pid)
          SubVol(Lay)=SubVol(Lay)+VNewi
          hTot=hTot+hE*dx
          hMean(Lay)=hMean(Lay)+hE*dx
          DeltW=DeltW+abs(WatIn(i,pid) -vNewi)
         
15    continue

      if(ATot.gt.0.) then
          if(lWat(pid).or.PLevel(pid).eq.0) hTot=hTot/ATot
          if(lTemp)               TTot=TTot/ATot
      end if

      fRE=1
      Grav=1
      dx1=x(2,pid)-x(1,pid)
      v1=-(Con(1,pid)+Con(2,pid))/2.*((hhNew(2,pid)-hhNew(1,pid))/dx1+Grav*fRE)
      fRE=1
      Grav=1
      dxN=x(n,pid)-x(N-1,pid)
      vN=-(Con(N,pid)+Con(N-1,pid))/2*((hhNew(N,pid)-hhNew(N-1,pid))/
     !dxN+Grav*fRE)

      write(76,111,err=901) t(pid)
      write(76,120,err=901) (i,i=1,NLay)
      write(76,130,err=901)
      write(76,140,err=901)   ATot,  (Area(i),i=1,NLay)
      write(76,150,err=901) Volume,(SubVol(i),i=1,NLay)
      if(iDualPor.gt.0) write(76,151,err=901) VolumeIm
      write(76,160,err=901) Change,(SubCha(i),i=1,NLay)
      write(76,170,err=901) hTot,  ( hMean(i),i=1,NLay)
      write(76,220,err=901) vN,v1
      wBalT=Volume-wVolI(pid)-wCumT(pid)
      write(76,230,err=901) wBalT
      ww=amax1(DeltW,wCumA(pid))
      if(ww.gt.1.e-25) then
      wBalR=abs(wBalT)/ww*100.
      write(76,240,err=901) wBalR
      end if
      write(76,130,err=901)
      return

*     Error when writing into an output file 
901   ierr=1
      return

110   format(
     !    /'----------------------------------------------------------'/
     !        ' Time       [T]',f14.4/
     !     '----------------------------------------------------------')
111   format(
     !    /'----------------------------------------------------------'/
     !        ' Time       [T]',e15.8/
     !     '----------------------------------------------------------')
120   format( ' Sub-region num.               ',9(I7,6x))
130   format( 
     !      '----------------------------------------------------------')
140   format( ' Area     [L]      ',e13.5,9e13.5)
150   format( ' W-volume [L]      ',e13.5,9e13.5)
151   format( ' W-volumeI[L]      ',e13.5,9e13.5)
160   format( ' In-flow  [L/T]    ',e13.5,9e13.5)
170   format( ' h Mean   [L]      ',e13.5,9e13.5)
180   format( ' HeatVol  [M/T2]   ',e13.5,10e13.5)
190   format( ' tMean    [K]      ',f13.3,10f13.3)
200   format( ' ConcVol  [M/L2] ',i1,1x,e13.5,10e13.5)
201   format( ' ConcVolIm[M/L2] ',i1,1x,e13.5,10e13.5)
202   format( ' SorbVolIm[M/L2] ',i1,1x,e13.5,10e13.5)
203   format( ' SorbVolIm2[M/L2]',i1,1x,e13.5,10e13.5)
210   format( ' cMean    [M/L3] ',i1,1x,e13.5,10e13.5)
211   format( ' cMeanIm  [M/L3] ',i1,1x,e13.5,10e13.5)
212   format( ' sMeanIm  [-]    ',i1,1x,e13.5,10e13.5)
220   format( ' Top Flux [L/T]    ',e13.5/
     !             ' Bot Flux [L/T]    ',e13.5)
230   format( ' WatBalT  [L]      ',e13.5)
240   format( ' WatBalR  [%]      ',f13.3)
250   format( ' CncBalT  [M]    ',i1,1x,e13.5)
260   format( ' CncBalR  [%]    ',i1,1x,f13.3)
      end

      SUBROUTINE GWF2UNSF1WVOLI(PID,THN,THO)
      
      USE GWFUNSFMODULE, ONLY:ths,sink,tlevel,ns,cmid,tpulse,nsd,llinear,
     !retard,wc,thsat,tdep,nsd, WCUMT,WCUMA,watin,wvoli,area, conc,ccuma,
     !ccumt,con,chpar, tdep, sorb,sorb2,laynum,cprevo, idualpor,cvoli,
     ! solin,lLinear ,ns,ths,con,thnewim,tholdim,cprev,x,nmat,nmatd,
     ! numnpd,npunsf,numnp,matnum,plevel,hhnew,dt,t,lmobim,lWat,lChem, 
     ! lEquil,lPrint,  
     !lVapor, lBact, lDensity,lDualNEq,js,volume,xmasschange
      
      INTEGER,INTENT(IN) :: PID
       logical lcentrif,ltemp,lwtdep
      integer   i
      
      dimension cMean(11,10),  ConSub(11,10), 
     ! ConVolIm2(11), hMean(10),
     1          TMean(10),SubVol(10),SubCha(10),ConVol(11),
     1          cTot(11),SubT(10) ,ConVolIm(11),
     1          ConSubIm(11,10),cMeanIm(11,10),cTotIm(11),ConSubIm2(11,10)
      double precision:: dxn, ThN(Numnpd,npunsf),
     !  ThO(Numnpd,npunsf) ,wbalt,temp(Numnpd),cnewi,vnewi,voldi,cEl
      
        N=numnp(pid)
      
      fRE=1 
      Grav=1
      ATot=0.
      Tr=293.15
      R=8.314
      cosalf=1
      lwat(pid)=1
      temp=20
      Volume=0.
      VolumeIm=0.
      Change=0.
      hTot=0.
      DeltW=0.


      do 15 i=N-1,1,-1
        j=i+1
        cEl=0.
        Mi=MatNum(i,pid )
        Mj=MatNum(j,pid)
        Lay=LayNum(i,pid)
        dx=x(j,pid)-x(i,pid)
        Area(Lay)=Area(Lay)+dx
        ATot=ATot+dx
         TT=(Temp(i)+Temp(j))/2.+273.15
          hE=(hhNew(i,pid)+hhNew(j,pid))/2.
          VNewi=dx*(ThN(i,pid)+ThN(j,pid))/2.
          VOldi=dx*(ThO(i,pid)+ThO(j,pid))/2.
          Volume=Volume+VNewi
          Change=Change+(VNewi-VOldi)/dt(pid)
          SubCha(Lay)=SubCha(Lay)+(VNewi-VOldi)/dt(pid)
          SubVol(Lay)=SubVol(Lay)+VNewi
          hTot=hTot+hE*dx
          hMean(Lay)=hMean(Lay)+hE*dx
          WatIn(i,pid)=vNewi
 15   continue
      
      wVolI(pid)=Volume
      
      END SUBROUTINE 