**************************************************
*           USER SUBROUTINE KMAT                 *
**************************************************  
C
C ** HCP, BCC, FCC ALL IN ONE: LARGE DEFORMATION **
C ** SLIP DIRECTIONS AND NORMALS NOT UPDATED **
C
C
      SUBROUTINE kmat(dtime,nsvars,usvars,xI,jelem,kint,knsdv,time,F,
     + L,iphase,C,stressV,dstressinc,totstran,dtotstran,
     + TEMP,DTEMP,vms,pdot, pnewdt)
C
      INCLUDE 'ABA_PARAM.INC'
      
      INTEGER, parameter:: M=3,N=3,L0=18,L1=12,L2=12,L4=7,KM=6,KN=6
      REAL*8,parameter  :: zero=1.0e-8
C
      !scalars
      INTEGER,intent(in) :: nsvars,jelem,kint,knsdv,iphase
      REAL*8,intent(in) :: dtime
      
      !arrays
      REAL*8,intent(inout):: usvars(nsvars)
      REAL*8,intent(in) :: time(2),F(3,3),L(3,3),xI(3,3)
      REAL*8,intent(out) :: C(6,6)     


C     *** USER DEFINED ARRAYS ***
      REAL*8 :: stressM(3,3),plasStrainInc(3,3), totalStrainInc(6), 
     + prod(M),tempNorm(M), tempDir(M), plasStrainRate(3,3),
     + xRot(M,N),VRES1(6),VRES2(6),plasStrainInc2(6),Lp(3,3),
     + totplasstran(6),Le(3,3),
     + tSigma(6,6),tStran(6,6),tSigmainv(6,6),prod6(6,6),
     + xsnt(M,N),xsnv(KM),xnsv(KM),xsnnst(KM,KN),xIden6(KM,KN),
     + xnst(M,N),tmat(KM,KN),
     + trialstress(6),xjfai(KM,KN),
     + xjfaiinv(KM,KN),stressV(6),fai(6),trialstressM(3,3),
     + dstressinc(6),xu(3,3),xuinver(3,3),
     + xstressmdef(M,N),xstressdef(6),
     + xStiff(6,6),xStiffdef(6,6),elasspin(3,3),plasspin(3,3),
     + tSig(6,6),tStr(6,6),tSiginv(6,6), T(6,6),Tinv(6,6),
     + tStraninv(6,6),gmatinv(3,3),gmatinvold(3,3),devstress(3,3),tSigTranspose(6,6),
     + compliance(6,6), print1(3,3),print2(3,3),tmpvec(6),print3(3,3),
     + print4(3,3),xrote(3,3),dtotstran(6),totstran(6),
     + totStrainRate(3,3),spin(3,3), tempstrain(3,3),curlfp(3,3),
     + curlfe(3,3),fp(3,3),fe(3,3),feinv(3,3),fpinv(3,3),xfinv(3,3),fpold(3,3)
       
      REAL*8 ::  vresth(6), dstranth(6),thermat(3,3),expanse33(3,3),
     + sigthv(6),sigthm(3,3), gmatinvnew(3,3) 
        
      REAL*8,dimension(:,:),allocatable :: xNorm,xDir
      REAL*8,dimension(:),allocatable :: tau, gammadot, gndall, gndold,
     +  gnde, gnds, burgerv, tauc, tau2, signtau 
      INTEGER,dimension(:),allocatable :: ids
      
      
      REAL*8 gammast, rho, TAUCSSD !gamma star term for RhoSSD calcs
      REAL*8 LpFeinv(3,3), matrix(3,3), update(3,3), maxdS
!      REAL(selected_real_kind(8,70)) :: rhoc
      
      character(len=*),parameter :: fmt2 = "(24(' ',(I2,1X)))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt24 = "(24(' ',(ES11.3,1X)))"
      
     
      
      !common/therm/ytemp,ytemprate


!----------------Note for use---------------
C     Tempature must larger than 273K since the material properties are interploted 
C     for only those great than room temperature, i.e. 20 degree C
!--------------------------------------------------------------------------------
!   
C removed (kint-1)*knsdv+ from all the SDV allocation, SDV reduced from 712 to 89, with SDV(1) = sum over k=1,8 of SDV(1+89*(k-1))
! individual integration point values of Fp stored in common block   ET updated 10/04/2015

C     *** MATERIAL CONSTANTS, ETC ***       
      SELECT CASE(iphase)
      CASE(0) !hcp
      nSys = L0
      ns = 9 ! number of screw systems     
      burger1 = 2.28E-4                         
      caratio = 1.57 !Zr. For Ti=1.587

      E1 = 289.38E3
      E2 = E1
      E3 = 335.17E3
      G12 = 132.80E3
      G13 = 162.50E3
      G23 = G13
      v12 = 0.09 
      v13 = 0.04
      v23 = v13
      !rhossd = rhossd+(sum(gndall))
      !XTAUC1 = TAUCSSD + (1.0 * G12 * burger1 * sqrt(rhossd))
      XTAUC1 = 15.2 ! 3 <a> basal
      XTAUC2 = 67.7 ! 3 <a> prismatic
      XTAUC3 = 1E6 ! 6 <a> pyramidal
      XTAUC4 = 2000.0 ! 6 <c+a>  secondary pymramidal
      burger2 = 4.242E-4 ! sqrt(a^2 + c^2)
        
      alpha1 = 9.5D-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1
	gammast = 0
      
      allocate(xNorm(L0,M),xDir(L0,M),tau(L0),gammadot(L0),gndall(L0+ns),
     + gnde(L0),gnds(ns),gndold(L0+ns),burgerv(L0),tauc(L0),
     + ids(L0),tau2(L0),signtau(L0),STAT=ialloc)
     
      burgerv(1:12) = burger1; burgerv(13:18) = burger2
      tauc(1:3) = XTAUC1; tauc(4:6) = XTAUC2; tauc(7:12) = XTAUC3
      tauc(13:18) = XTAUC4
          
            
      case(1) !bcc
      nSys = L1
      ns = 4
      XTAUC = 300.
      burger = 2.48E-4  
      
      !bcc properties from Dave's paper
      E1 =  204E3; E2 = E1; E3 = E1
      v12 = 0.29; v13 = v12; v23 = v12
      G12 = 78E3; G13 = G12; G23 = G12
      
      !! makes code crash
      !E1 = 125.0E3; E2 = E1; E3 = E1
      !v12 = 0.35; v13 = v12; v23 = v12
      !G12 = 116.3E3; G13 = G12; G23 = G12
      
      
      alpha1 = 9.5e-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1
      
      allocate(xNorm(L1,M),xDir(L1,M),tau(L1),gammaDot(L1),gndall(L1+ns),
     + gnde(L1),gnds(ns),gndold(L1+ns),burgerv(L1),tauc(L1),
     + ids(L1),tau2(L1),signtau(L1),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
     
      case(2) !fcc
      nSys = L2
      XTAUC =  220.0 !230.0 ! (563.4/sqrt(6) = 230.0)
      burger = 3.5072e-4
      
   !   E1 = 207.0E+3; E3 = E1
   !   v12 = 0.28; v13 = v12
   !   G12 = E1/(2.0*(1.0+v12)); G13 = G12
      
       E1 = 123.3E3; E2 = E1; E3 = E1
      v12 = 0.38; v13 = v12; v23 = v12
      G12 = 126.2E3; G13 = G12; G23 = G12
      
      alpha1 = 13.0e-6; alpha2 = alpha1; alpha3 = alpha1
      
      gammast = 0.05
     
      
      allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),gndall(L2+ns),
     + gnde(L2),gnds(ns),gndold(L2+ns),burgerv(L2),tauc(L2),
     + ids(L2),tau2(L2),signtau(L2),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
      
      case(3) !carbide
      nSys = L2
      XTAUC = 230.0D10
      burger = 3.5072e-4
      
      E1 = 207.0E+4; E3 = E1
      v12 = 0.28; v13 = v12
      G12 = E1/(2.0*(1.0+v12)); G13 = G12; G23 = G12
      
      alpha1 = 4.5e-6; alpha2 = alpha1; alpha3 = alpha1
      
      allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),gndall(L2+ns),
     + gnde(L2),gnds(ns),gndold(L2+ns),burgerv(L2),tauc(L2),
     + ids(L2),tau2(L2),signtau(L2),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
      
      case(4) !olivine
      nSys = L4
      ns = 2
      XTAUC = 1720.0
      burger = 3.5072e-4
      ! see anisotropy.m
      E1 = 286.4E3; E2 = 165.1E3; E3 = 336.5E3
      v12 = 0.258; v13 = 0.223; v23 = 0.283
      G12 = 64.0E3; G13 =  77.0E3; G23 = 79.0E3
      
      alpha1 = 4.5e-6; alpha2 = alpha1; alpha3 = alpha1
      
      allocate(xNorm(L4,M),xDir(L4,M),tau(L4),gammaDot(L4),gndall(L4+ns),
     + gnde(L4),gnds(ns),gndold(L4+ns),burgerv(L4),tauc(L4),
     + ids(L4),tau2(L4),signtau(L4),STAT=ialloc)
     
      burgerv = burger
      tauc(1:L4) = xtauc*(/2.0,1.0,1.0,1.5,4.0,4.0,1.0/)
      
	  case(5) !bcc
      nSys = L1
      ns = 4
      XTAUC = 360.
      burger = 2.74E-4  
      
      !bcc properties from Dave's paper
      E1 =  421E3; E2 = E1; E3 = E1
      v12 = 0.28; v13 = v12; v23 = v12
      G12 = 164.4E3; G13 = G12; G23 = G12
      
      !! makes code crash
      !E1 = 125.0E3; E2 = E1; E3 = E1
      !v12 = 0.35; v13 = v12; v23 = v12
      !G12 = 116.3E3; G13 = G12; G23 = G12
      
      
      alpha1 = 9.5e-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1
      
      allocate(xNorm(L1,M),xDir(L1,M),tau(L1),gammaDot(L1),gndall(L1+ns),
     + gnde(L1),gnds(ns),gndold(L1+ns),burgerv(L1),tauc(L1),
     + ids(L1),tau2(L1),signtau(L1),STAT=ialloc)
     
      burgerv = burger
      tauc = xtauc
	  
      case default
      WRITE(6,*)
      WRITE(6,*)"Not sure what crystal type. Material constants."
      END SELECT

      
      h = 1E-5     
      kslip = 5
      !0 = original slip rule with no GND coupling.
      !5 = original slip rule with GND coupling, with a new parameter to make rhoc more physically reasonable, slip included.
      !6 = Slip rule with constant Activation Volume
      
      gndon =1 !1=on, 0=off



C     *** ZERO ARRAYS ***

      result = 0.0; totstran = 0.0; totplasstran = 0.0;devstress=0.
      plasStrainInc2=0.; plasStrainRate = 0.
      xStiff=0.0; xStiffdef=0.0; C=0.0; trialstressM=0.0
      tmat=0.0;xjfai=0.
      xjfaiinv=0.;plasStrainInc=0.;trialstress=0.; stressV=0.
      fai=0.;dstressinc=0.;xIden6=0.; xu=0.
      xuinver=0.;xRot=0.;prod=0.
      T=0; Tinv=0 !NEW ADDITIONS
      tSig=0.;tStr=0.;tSiginv=0.;tStraninv=0.
      tSigmainv=0.;gmatinv=0.;Lp = 0.;compliance=0.;Le = 0.
      print1=0.;print2=0.; tmpvec = 0.;print3 = 0.;print4 = 0.
      thermat =0.
      vresth=0.; dstranth=0.;expanse33=0.;sigthv=0.;sigthm=0.
      xrote=0.;tempstrain=0.;spin=0.;totStrainRate=0.
      curlfp=0.;curlfe=0.;fp=0.;fe=0.;feinv=0.;fpinv=0.;xfinv=0.
      gmatinvnew = 0.
      
      
      DO I=1,KM; xIden6(I,I)=1.; END DO      
      DO I=1,M; xRot(I,I) = 1.; END DO
      
      
      !Zero allocate arrays from above!
      xNorm=0.; xDir=0.; tau=0.; gammadot=0.; signtau=1.
      gndall=0.; gndold=0.; ids=0; tau2=0.
      gndcas=0.;gndcap=0.;gndapy=0.;gndapr=0.;gndab=0.;rhognd=0.
      
      !Ben's experimental data storage groups
      capyramedge=0.;capyramscrew=0.;apyramedge=0.;aprismedge=0.
      abasedge=0.;ascrew=0.  
      
      thermat(1,1) = alpha1; thermat(2,2)=alpha2; thermat(3,3) = alpha3
      
C     *** SET UP ELASTIC STIFFNESS MATRIX IN LATTICE SYSTEM ***   
      compliance(1,1:3) = (/1./E1,-v12/E1,-v13/E1/)
      compliance(2,2:3) =         (/1./E2,-v23/E2/)
      compliance(3,3:3) =                 (/1./E3/)
      compliance(4,4:4) =                       (/1./G12/)
      compliance(5,5:5) =                       (/1./G13/)
      compliance(6,6:6) =                       (/1./G23/)
C
      DO i=2,6
         DO j=1,i-1
            compliance(i,j)=compliance(j,i)
         END DO
      END DO
      
      CALL lapinverse(compliance,6,info,xStiff)
     
    
C     *** INITIALIZE USER ARRAYS ***

      DO i=1,3
        DO j=1,3
         gmatinv(i,j) = usvars(j+(i-1)*3)
        END DO
      END DO
C
      gmatinvold = gmatinv
      p = usvars(10)
C
      DO i=1,6
        totplasstran(i) = usvars(10+i)
      END DO
	  
	  
C
     
	  DO i=1,6
        totstran(i) = usvars(16+i)
      END DO
		
      DO i=1,6
        xstressdef(i) = usvars(47+i)
      END DO
C
      rhognd = usvars(26)  
C
      !r= usvars((kint-1)*knsdv+56)       
      rhossd = usvars(54)
      
      r = usvars(56)
                
      DO i=1,nSys
        gndold(i) = usvars(56+i)
      END DO
C
      DO i=1,3
        DO j=1,3
          fp(i,j) = usvars(80+j+((i-1)*3))
        END DO
      END DO              
C
      DO i=1,3
        DO j=1,3 
         curlfp(i,j) = usvars(37+j+(i-1)*3)
        END DO
      END DO
      
	  slip = usvars(47)
	   
	 
		rhodefect = usvars(35)
		soluteforce = usvars(36)
		!print*, 'reading in', rhodefect
	 
	
	 
	   

C     *** DIRECTIONS FROM LATTICE TO DEFORMED SYSTEM ***
        
      CALL kdirns(gmatinv,iphase,nSys,xDir,xNorm,caratio)
        

C     *** STIFFNESS FROM LATTICE TO DEFORMED SYSTEM ***

      !CALL rotord4sig(gmatinv,tSig)
      !CALL rotord4str(gmatinv,tStr)
      !CALL lapinverse(tSig,6,info2,tSiginv)
      !
      !prod6 = matmul(tSiginv,xStiff)      
      !xStiffdef = matmul(prod6,tStr)
      CALL rotord4sig(gmatinv,tSig)
  
      prod6 = matmul(tSig,xStiff)      
      tSigTranspose = transpose(tSig)
      
      xStiffdef = matmul(prod6,tSigtranspose)

      
!      expanse33 = matmul(matmul(gmatinv,thermat),transpose(gmatinv))
!C      expanse33 = expanse33*ytemprate*dtime !dstrain = alpha*dT
!      expanse33 = expanse33*DTEMP !dstrain = alpha*dT
      
      CALL kmatvec6(expanse33,dstranth)
      dstranth(4:6) = 2.0*dstranth(4:6)
            

    
C     *** DETERMINE INCREMENT IN TOTAL STRAIN (6X1 ***     

      tempstrain=(L+transpose(L))*0.5*dtime
      spin=(L-transpose(L))*0.5 

      CALL kmatvec6(tempstrain,dtotstran)
      dtotstran(4:6) = 2.0*dtotstran(4:6)


C     *** COMPUTE TRIAL STRESS ***

      stressV = xstressdef ! old stress
      trialstress = stressV+matmul(xStiffdef,dtotstran)-
     + matmul(xStiffdef,dstranth)            
      CALL kvecmat6(trialstress,trialstressM) 
            
      CALL kvecmat6(stressV,stressM) 
      trialstressM = trialstressM + (matmul(spin,stressM) - 
     + matmul(stressM,spin))*dtime 


 
C     *** CALCULATE RESOLVED SHEAR STRESS ON A SLIP SYSTEMS ***

      
      DO I=1,nSys
        tempNorm = xNorm(I,:); tempDir = xDir(I,:)
        prod = matmul(trialstressM,tempNorm)
        tau(I)= dot_product(prod,tempDir)
        signtau(I) = 1.d0      
        IF(tau(I) .LT. 0.0) THEN
          tau(I) = -1.E0*tau(I)
        !  DO K=1,3
        !    xDir(I,K)=-1.E0*xDir(I,K)
        !  END DO
            signtau(I) = -1.d0
        END IF
      END DO
        
      
		tauc = 300. + 0.1*19.344*sqrt(rhognd)!19.344 = G*burger
	 
	  
		
	  
	  
      xtau = maxval(tau/tauc)  

         
C     *** PLASTIC DEFORMATION ***

      IF (xtau >= 1.0 ) THEN
      
      debug = 0
      do while (debug == 1 .and. kint == 1 )           
          debugwait = 0
      end do
      

      
      faivalue=1.
      xacc=1.e-8
      iter=0
      fpold = fp


C     *** USE NEWTON METHOD TO DETERMINE STRESS INCREMENT ***

      DO WHILE (faivalue .gt. xacc)      
      
      iter=iter+1


      !============================================================================   
      !  Slip rule:
      !  Returns Lp and tmat required to define the material jacobian.  
      !============================================================================  

      IF (kslip == 0) THEN
      !Original slip rule with no GND coupling i.e. using alpha and beta
      CALL kslip0(xNorm,xDir,tau,tauc,caratio,dtime,nSys,r,iphase,Lp,
     + tmat) !need to pass in signtau and modify slip law
     
      ELSE IF (kslip == 5) THEN
      !Original slip rule with GND coupling            
      CALL kslip5ET(xNorm,xDir,tau,signtau,tauc,burgerv,rhossd,rhognd,caratio,
     + dtime,nSys,r,iphase,Lp,tmat,gammaDot)
          
      ELSE
      !Original slip rule with no GNDs eg no hardending 
      CALL kslip6(xNorm,xDir,tau,tauc,burgerv,caratio,      
     + dtime,nSys,r,iphase,Lp,tmat,TEMP) !need to pass in signtau and modify slip law
           
      END IF
           
      
     
      if(any(tmat /= tmat)) then
     
          call Mutexlock( 10 )   ! lock Mutex #5
         
          pnewdt = 0.5 ! if sinh( ) has probably blown up then try again with smaller dt
          write(*,*) "*** WARNING tmat  = NaN: jelem, kint, time: ", jelem, kint, time     
          call MutexUnlock( 10 )   ! unlock Mutex #5
               
          return
      end if   
   
      !============================================================================  


C     *** DETERMINE PLASTIC STRAIN INCREMENETS FOR UEL

      plasStrainInc = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainInc,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6)            



C     *** CALCULATE THE STRESS INCREMENT ***

      xjfai =  xIden6 + matmul(xStiffdef,tmat)
      CALL lapinverse(xjfai,6,info3,xjfaiinv)
!      IF(info3 /= 0) write(6,*) "inverse failure: xjfai in kmat"
      fai = trialstress - stressV - matmul(xStiffdef,plasStrainInc2)
      dstressinc = matmul(xjfaiinv,fai)
! for stability improvement keep relative values (and sign) of stress increment but scale dS to have max value of 10.0
      maxdS = maxval(abs(dstressinc))
      if (maxdS > 10.0) then
		dstressinc = dstressinc*10.0/maxdS
          !write(*,*) 'reducing dS by', 10.0/maxdS, jelem, kint, time(1)
      end if

      stressV = stressV + dstressinc
      CALL kvecmat6(stressV,stressM)      
      faivalue = sqrt(sum(fai*fai))        


C     *** UPDATE RESOLVED SHEAR STRESS ACCORDING TO NEW STRESS ***
      ids=0; tau2=0. 
      DO I=1,nSys
      
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)    
          prod = matmul(stressM,tempNorm)
          tau(I)= dot_product(prod,tempDir)
          signtau(I) = 1.d0   
          IF(tau(I) < 0.0) THEN
            tau(I) = -1.E0*tau(I)
          !  DO K=1,3
          !    xDir(I,K)=-1.E0*xDir(I,K)
          !  END DO
            signtau(I) = -1.d0
          END IF
          !=============================== ids
          IF(tau(i)/tauc(i) >= 1.0) THEN
            ids(i) = i; tau2(i) = tau(i)/tauc(i)
          END IF
          !===============================
          
      END DO
        
      xtau = maxval(tau/tauc) 
          
      
      IF (iter .gt. 100) THEN

          call Mutexlock( 11 )   ! unlock Mutex
          
          pnewdt = 0.5
          WRITE(*,*) "WARNING NEWTON LOOP NOT CONVERGED: jelem, kint, time:", jelem, kint, time(1)
          
          call MutexUnlock( 11)   ! unlock Mutex 
          return
          !CALL XIT 
      END IF
               
      !!*** THE END OF NEWTON ITERATION ***
      END DO

     
	 !DO I=1,nSys
	 
		!slip = slip + abs(gammaDot(I))*dtime
	 !END DO
	
      
	 
        !IF (tauc(1) .gt. 300.) THEN           
          !soluteforce = 1440. *(1.0-((slip-0.4)/(0.5-0.4)))
		 ! soluteforce =700. * exp(-slip/0.015)!((0.01-slip)/0.01)
        !END IF
	    
		! IF (soluteforce < 0.) THEN
			! soluteforce = 0.
		! END IF
	  

C     *** NOW CALCULATE THE JACOBIAN***

      C = matmul(xjfaiinv,xStiffdef)
       
       
C     *** ROTATE STRESS BACK TO GLOBAL SYSTEM *** 

      xstressmdef = stressM


C     *** UPDATE OUTPUT VARIABLES ***

      plasStrainrate=(Lp+transpose(Lp))*0.5       
      pdot=sqrt(2./3.*sum(plasStrainrate*plasStrainrate))
      
      dp = pdot*dtime 
      
                       
	p = p + pdot*dtime
      r = r + h*pdot*dtime   
      
      
C     *** UPDATE PLASTIC DEFORMATION GRADIENT
    
      print2 = 0.; print3 = 0.
      print2 = xI - Lp*dtime      
      CALL kdeter(print2,deter)      
      IF (deter /= 0.0) THEN
         CALL lapinverse(print2,3,info4,print3)
         fp = matmul(print3,fpold)
      ELSE
         fp = fp
      END IF  
      
      
    !=========================================================================
C     *** DETERMINE DENSITY OF GNDs
    ! Definitely needs to be inside plasticity loop
    ! And should use the directions that have been modified according to tau!
    !=========================================================================
    ! Edge-screw separated
      
      
      IF (gndon == 0) THEN !Switching GND evolution on and off!
       gndall = 0.
            
      ELSE IF (xtau >= 1.0 ) THEN  
      
      CALL kgndl2ET(curlfp,xNorm,xDir,tau2,ids,burgerv,iphase,nSys,jelem,
     +    kint,time, gndall,gnde,gnds,ns)
          !Extra storage has pure edge of pure screw groupings!
       
      ELSE
      
      gndall = 0.
      write(*,*) "warning:  no active slip systems ids = 0: jelem,kint,time",jelem,kint,time
      ! kgndl2 will crash as ids = 0 (no active slip systems) this happens if the Newton loop causes tau to drop below tauc
      END IF
      

      
      gndold = gndall
      rhognd = sum(gndall)
      
      if (iphase == 0) then          
         !edge
          gnde = sqrt(gnde*gnde)
          abasedge = sum(gnde(1:3))
          aprismedge = sum(gnde(4:6))
          apyramedge = sum(gnde(7:12))
          capyramedge = sum(gnde(13:18))
      !screw
          gnds = sqrt(gnds*gnds)
          ascrew = sum(gnds(1:3))
          capyramscrew = sum(gnds(4:9))
          
          gndab  = abasedge + ascrew !sum(gndall(1:3))
          gndapr = aprismedge !sum(gndall(4:6))
          gndapy = apyramedge !sum(gndall(7:12))
          gndcas = capyramedge+capyramscrew ! sum(gndall(13:18))
         end if              


    !=========================================================================
    ! SSD Evolution 
     
      rhossd = rhossd + (gammast*pdot*dtime)
      
    !=========================================================================
      
C     *** ELASTIC DEFORMATION ***     
      ELSE
          xstressmdef = trialstressM
          C = xStiffdef      
      END IF
      
!     CALL kvecmat6(xstressdef,stressM) 
!     xstressmdef = xstressmdef + (matmul(spin,stressM) - matmul(stressM,spin))*dtime
      
      CALL kmatvec6(xstressmdef,xstressdef) !output stress
      devstress = xstressmdef - 1./3.*trace(xstressmdef)*xI
      vms = sqrt(3./2.*(sum(devstress*devstress))) !von mises stress 

    !   call MutexLock( 1 )      ! lock Mutex #1   
      !write(*,*) gmatinv
      !write(*,*) kint, jelem 
    !  call MutexUnLock( 1 )      ! lock Mutex #1   

C     *** ORIENTATION UPDATE ***
      !Assuming that all rigid body rotatoin is lumped into Fe and that the elastic strians are small 
!     then the elastic spin is We = d(Fe)/dt inv(Fe)
      !L = We + Fe Lp inv(Fe) therefore 
      !We = L - Fe Lp inv(Fe)
      ! G(t+dt) = G(t) + We G(t)dt dt or an implicit update is G(t+dt)  = G(t)exp[We(t+dt)dt]  ~ inv[I - We(t+dt) dt] G(t) 
      
      ! We need Fe and inv(Fe) using F = Fe Fp gives Fe = F.inv(Fp)
      CALL kdeter(Fp,deter)      
      
      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info5,Fpinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"
         Fe = matmul(F,Fpinv)          
      ELSE
         write(*,*) "Error in orientation update: finding inv(Fp)",jelem,kint, kinc
         call XIT 
      
      END IF  
      
      
      CALL kdeter(Fe,deter)      
      
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)
!         IF(info5 /= 0) write(6,*) "inverse failure: print3 in kmat"         
      ELSE
          write(*,*) "Error in orientation update: finding inv(Fe)",jelem,kint, kinc
         call XIT 
      
      END IF        
            
      LpFeinv = 0.; 
      LpFeinv = matmul(Lp, Feinv)
      Le = L - matmul(Fe,LpFeinv)        
      elasspin=(Le-transpose(Le))*0.5
      matrix = xI - elasspin*dtime      
      CALL kdeter(matrix,deter)      
      
           
          
      
      
      IF (deter /= 0.) THEN
         update = 0.
         CALL lapinverse(matrix,3,info5,update)
         IF(info5 /= 0) write(*,*) "inverse failure: print3 in kmat"
         !print3 = 0.
         !print3 = gmatinv + dtime*matmul(elasspin,gmatinv)
         gmatinvnew = matmul(update,gmatinvold)                  
       
         !write(*,*) "gmatinv, print3", gmatinv, print3
      
      ELSE         
         gmatinvnew = gmatinvold
      write(*,*) "WARNING gmatinv not updated at jelem,kint, kinc:", jelem,kint, kinc
      END IF      

      gmatinv =  gmatinvnew            
      
       if (maxval(gmatinv) > 1) then
          write(*,*) "something very wrong with gmatinv"
          call XIT
       end if
       
  


C     ******************************************
C     Heat                     !Changed 19/06/2014
C      rho=7.8e-3
C      specHeat=460
C      beta=0.95
C      ytemprate=vms*pdot*beta/(rho*specHeat)
C     ******************************************      

C     *** UPDATE STATE VARIABLES *** ! Free: None

      DO i=1,3
        DO j=1,3
        usvars(j+(i-1)*3) = gmatinv(i,j)
        END DO
      END DO

      usvars(10) = p
 
      DO i=1,6
        usvars(10+i) = totplasstran(i) + 
     +                                  plasStrainInc2(i)
      END DO
C
      DO i=1,6
        usvars(16+i) = totstran(i) + dtotstran(i)
      END DO
      
         
      !23-25 are rotations stored in UEL
          
     
      usvars(26) = rhognd !all
      usvars(27) = gndab  !a basal
      usvars(28) = gndapr !a prismatic  
      usvars(29) = gndapy !a pyramidal
      usvars(30) = 0.0 !c+a secdonary pyramidal   
      usvars(31) = gndcas !c+a secdonary pyramidal 

      !38-46 are curlfp terms calulated in gradient routine (kcurl)
      
      DO i=1,6
       usvars(47+i) = xstressdef(i)
      END DO

      usvars(32) = maxval(plasStrainrate)
      usvars(33) = pdot
      usvars(34) = xtau !XTAUC1 
      usvars(54) = rhossd
      usvars(55) = vms
      usvars(56) = r

	  usvars(47) = slip
	  
      !GNDs on indiviual systems
      !max(nSys) is currently limited to 24. IF all 48 of bcc is needed, storage should be raised to match that!
      DO i=1,nSys+ns
       usvars(56+i) = gndold(i)
      END DO
!
      
      DO i=1,3
       DO j=1,3 
        usvars(80+j+(i-1)*3) = fp(i,j)
        usvars(89+j+(i-1)*3) = F(i,j)
       END DO
      END DO   
     
	  usvars(35) = rhodefect
	  usvars(36) = soluteforce
           


!C     *** WRITE RESULTS TO DUMMY UMAT!!
!
!      DO i=1,knsdv
!        ksdv(kint,i) = usvars(i)
!      END DO

           
      DEALLOCATE(xNorm,xDir,tau,gammadot,gndall,gndold,burgerv,
     + tauc,ids,tau2,signtau)

      RETURN

      END

