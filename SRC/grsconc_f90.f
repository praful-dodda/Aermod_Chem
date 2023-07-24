

      SUBROUTINE ZENANG(YR,JDAY,HR,OLAT,ZA)
! ----------------------------------------------------------------------------------- !
!     Written by Alejandro Valencia (AMV)
!
!     RLINE v1.6, October 2015
!
!     This code (ZENith ANGle) returns the zenith angle given latittude and time.
! ----------------------------------------------------------------------------------- !


        !-------------------------------------------------------------------------------
        ! argument list variables
        !-------------------------------------------------------------------------------

        use line_source_data, only: double, pi
        IMPLICIT NONE
        
        INTEGER,  INTENT(IN)   :: yr        ! Year
        INTEGER,  INTENT(IN)   :: jday      ! Day of Year
        INTEGER,  INTENT(IN)   :: hr        ! Hour 
        REAL(kind=double),     INTENT(IN)   :: olat      ! latitude (+ = north)
        REAL,     INTENT(OUT)  :: za        ! Zenith angle in deg at time
                    
        !-------------------------------------------------------------------------------
        ! internal variables
        !-------------------------------------------------------------------------------

        REAL,  PARAMETER  :: rdpdg  = 0.0174532925199433     ! radians per degree
        REAL    :: Dl,Eob,Njd,Lec,Lm,Gm,d,ha,ea,yr4

        ! Adjust for 2
            if (YR < 50)then
            YR4 = YR + 2000
            else
            YR4 = YR + 1900
            endif            

        ! Number of leap days since or before the year 2000 (Dl)
            if (YR4 >= 2001)then
            Dl=(YR4-2001)/4
            else
            Dl=(YR4-2000)/4-1
            endif

        ! Obliquity of the ecliptic (Eob)
            Njd=364.5+(YR4-2001)*365+Dl+JDAY
            Eob=(23.439-0.0000004*Njd)*rdpdg

        ! Ecliptic Longitude if the Sun
            Lm=(280.460+0.9856474*Njd)*rdpdg
            Gm=(357.528+0.9856003*Njd)*rdpdg
            Lec=Lm+1.915*SIN(Gm)+0.020*SIN(2*Gm)

        ! Sine of Declination Angle
            D = ASIN(SIN(Eob)*SIN(Lec))

        ! Local Hour Angle
            HA = ((2*180*rdpdg*(HR-12)/(24))/RDPDG)*rdpdg

        ! Zenith angle in degrees
            ZA=ACOS((SIN(OLAT*rdpdg)*SIN(D))+(COS(HA)*COS(OLAT*rdpdg)*COS(D)))
            ZA= ZA/rdpdg

        ! Elevation Angle    
        EA=90-ZA

      END SUBROUTINE zenang
! ----------------------------------------------------------------------------------- !
! The Research LINE source (R-LINE) model is in continuous development by various     !
! groups and is based on information from these groups: Federal Government employees, !
! contractors working within a United States Government contract, and non-Federal     !
! sources including research institutions.  These groups give the Government          !
! permission to use, prepare derivative works of, and distribute copies of their work !
! in the R-LINE model to the public and to permit others to do so.  The United States !
! Environmental Protection Agency therefore grants similar permission to use the      !
! R-LINE model, but users are requested to provide copies of derivative works or      !
! products designed to operate in the R-LINE model to the United States Government    !
! without restrictions as to use by others.  Software that is used with the R-LINE    !
! model but distributed under the GNU General Public License or the GNU Lesser        !
! General Public License is subject to their copyright restrictions.                  !
! ----------------------------------------------------------------------------------- !


      subroutine ierno2(zag,rk1)
      
! ----------------------------------------------------------------------------------- !
!     Code adapted by Alejandro M. Valencia (AMV)     
!    
!     RLINE v1.6, October 2015
!
!     This subroutine uses the zenith angle and calculates the surface short wave 
!     flux to compute the no2 photolysis rate coefficient.  The rate coefficient
!     is used in the integration of VOC for smog produced. This
!     subroutine is based on HYSPLIT's IERNO2 subroutine written at the
!     Air Resources Laboratory by Roland Draxler.
!
! ----------------------------------------------------------------------------------- !
        implicit none
      
        real, intent(in)  :: zag      ! zenith angle (deg)
        real, intent(out) :: rk1     ! rate coefficient for no2 photolysis
      
        real, parameter   :: rpdg = 0.01745329  ! radians per degree
        real, parameter :: solc     = 1104.0     ! cloud top solar const (w/m2)
        
        REAL  :: zar,sea,swf
      
      ! Sine of the elevation angle (rad)
        sea = sin((90-zag)*rpdg)
      
      ! zenith angle in radians from elevation angle
        zar=zag*rpdg
      
      
      ! Calculate SWF : short-wave flux (w/m2) 
      ! Solar Radiation at the Earth's Surface assuming clear skys
      
      ! SWF = SEA*SOLC*T
      ! where SEA is the sine of the elevation angle in radians, SOLC is the solar 
      ! constant incident at the top of the cloud layer and T is the fraction transmitted
      ! through the clouds (assumed as 1). Using the normal top-of-the-atmosphere solar constant to be about 1380 W
      ! m-2 and assuming that on average about 20% of the radiation is absorbed or reflected back into
      ! space (Sellers, 1972) by a clear atmosphere, then S = 1100 W m-2. 
      
      
        swf = max(0.0,sea*solc)
      
      ! rate coefficient for no2 photolysis
      
        if(zag.ge.0.0.and.zag.lt.47.0)then
           rk1=(4.23e-04+(1.09e-04/cos(zar)))*swf
        elseif(zag.ge.47.0.and.zag.lt.64.0)then
           rk1=5.82e-04*swf
        elseif(zag.ge.64.0.and.zag.le.90.0)then
           rk1=(-0.997e-04+1.2e-03*(1.0-cos(zar)))*swf
        else
           rk1=0.0
        end if
      
      end subroutine ierno2
! ----------------------------------------------------------------------------------- !
! The Research LINE source (R-LINE) model is in continuous development by various     !
! groups and is based on information from these groups: Federal Government employees, !
! contractors working within a United States Government contract, and non-Federal     !
! sources including research institutions.  These groups give the Government          !
! permission to use, prepare derivative works of, and distribute copies of their work !
! in the R-LINE model to the public and to permit others to do so.  The United States !
! Environmental Protection Agency therefore grants similar permission to use the      !
! R-LINE model, but users are requested to provide copies of derivative works or      !
! products designed to operate in the R-LINE model to the United States Government    !
! without restrictions as to use by others.  Software that is used with the R-LINE    !
! model but distributed under the GNU General Public License or the GNU Lesser        !
! General Public License is subject to their copyright restrictions.                  !
! ----------------------------------------------------------------------------------- !


        subroutine grates(zenith,temp,k,kp,activity)
 
! ----------------------------------------------------------------------------------- !
!       Code adapted by Alejandro M. Valencia (AMV)
!
!       This routine calculates all photolytic and molar reaction
!       rates for the GRS mechanism. It is based on the HYSPLIT's GRATES
!       subroutine.
! ----------------------------------------------------------------------------------- !

        real,       intent(out)     :: k(7)

        real kp
        parameter (t1=0.00316)

!       Compute rate constants
!       Set and adjust photolysis rates

        tempi=1.0/temp
        ! k(4)=9.24e+05*tempi*exp(-1450.0*tempi)
        k(4)=2643.0*exp(-1370.0*tempi)

        if(zenith.lt.90.0)then
            ! k(1)=kp*activity*exp(-1000.0*4.7*(tempi-t1))
            k(1)=kp*10000.0*exp(-4710.0*tempi)

            ! k(2)=(3.58e+06)*tempi
            k(2)=5482.0*exp(242.0*tempi)

            k(3)=kp
            k(5)=10000.
            k(6)=125.
            k(7)=k(6)

!       If nightime then set photolysis rates to zero
        else
            k(1)=0.
            k(2)=0.
            k(3)=0.
            k(5)=0.
            k(6)=0.
            k(7)=0.
        endif



        end subroutine GRATES
! ----------------------------------------------------------------------------------- !
! The Research LINE source (R-LINE) model is in continuous development by various     !
! groups and is based on information from these groups: Federal Government employees, !
! contractors working within a United States Government contract, and non-Federal     !
! sources including research institutions.  These groups give the Government          !
! permission to use, prepare derivative works of, and distribute copies of their work !
! in the R-LINE model to the public and to permit others to do so.  The United States !
! Environmental Protection Agency therefore grants similar permission to use the      !
! R-LINE model, but users are requested to provide copies of derivative works or      !
! products designed to operate in the R-LINE model to the United States Government    !
! without restrictions as to use by others.  Software that is used with the R-LINE    !
! model but distributed under the GNU General Public License or the GNU Lesser        !
! General Public License is subject to their copyright restrictions.                  !
! ----------------------------------------------------------------------------------- !

SUBROUTINE GCHEM(DIFTM,CON,CON1,K,active1)

!-----------------------------------------------------------------------------
! This file contains subroutines to compute ozone formation using the GRS
! chemistry model.  The routines were developed by Martin Cope, Victoria EPA,
! and later modified by G.D. Hess, BMRC, BoM, Melbourne, Australia.
! http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/lib/sorc/arl_hysplit/grseqn.f
! 
! Adapted to RLINE by Alejandro M. Valencia (AMV)
!
!-----------------------------------------------------------------------------
 
! number of species and reactions
  PARAMETER (NOSPEC=6,NOREAC=7)

! REACTION RATES AND TOTAL FORMATION/LOSS TERMS
  
  REAL A(NOREAC),AP(NOREAC),B(NOREAC),AOB(NOREAC),R(NOREAC)
  REAL ZERO(NOREAC)
  REAL LMBDA,active1
  
! CONCENTRATION ARRAYS. CON= INPUT ARRAY
!                       CON1=PREDICTED CONCENTRATION
  REAL CON(NOREAC),CON1(NOREAC),K(NOREAC)
  


  DATA ZERO /7*1.0E-15/

  SAVE AOB,A,AP,B,R,ZERO
  
  elpstm = 0.0
  R      = 0.0
  A      = 0.0
  AP     = 0.0
  AOB    = 0.0
  B      = 0.0

! ENTRY FOR INTEGRATION STAGE. (PECE MODE). CONTINUE UNTIL
! FULL TIME STEP COMPLETED
! ALSO ENTRY POINT FOR NEW PREDICTION

  DO WHILE (elpstm.LT.diftm)

     DT=diftm-elpstm
 

!    CALCUATE TIME RATE OF CHANGE OF ALL SPECIES.
!    CHECK D[NO]/DT FOR CONVERGENCE CRITERIA
     CALL ABCALC(NOSPEC,A,B,K,R,CON,active1)

     LMBDA=0.001
 
!    WORK OUT MAXIMUM TIME STEP. ONLY CONSIDER SPECIES
!    WITH 'APPRECIABLE' CONCENTRATIONS (NO+NO2+O3)
!    INTEGRATE FORWARD.
 
     CONSUM=CON(1)+CON(2)+CON(4)
     CONSUM=CONSUM*0.001

     DO I=1,NOSPEC
        IF(B(I).GT.0)THEN
           AOB(I)=A(I)/B(I)
           IF(CON(I).GT.CONSUM)THEN
               IF(con(i)-aob(i).ne.0.0) THEN
                 F=LMBDA*AMAX1(AOB(I),CON(I))/ABS(CON(I)-AOB(I))
               ENDIF
               IF(F.LT.1.) THEN
                 DT=AMIN1(DT,-ALOG(1.0-F)/B(I))
               ENDIF
          ENDIF
        ENDIF
     END DO           
     
     if(dt.le.0)THEN
        !write(*,*)'*ERROR*: GRS chems subroutine dt<0'
        dt=0.001
     END IF 

     CALL INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO)

!    CORRECTED SOLUTION
     CALL ABCALC(NOSPEC,AP,B,K,R,CON1,active1) 

     DO I=1,NOSPEC
        A(I)=(A(I)+AP(I))*0.5
        IF(B(I).GT.0)AOB(I)=A(I)/B(I)
     END DO

     CALL INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO) 
     DO I=1,NOSPEC
        CON(I)=CON1(I)
     END DO
     elpstm=elpstm+DT

  END DO
 
  END SUBROUTINE gchem

!=======================================================

SUBROUTINE ABCALC(NOSPEC,A,B,K,R,CON1,active1)
 
! ROUTINE CALCULATES FORMATION AND REACTION
! RATES FOR ALL SPECIES
 
  REAL NO,NO2
  REAL A(7),B(7),K(7),R(7)
  REAL CONSP(200),CON1(6)
 
! WORKING ARRAY EQUIVALENCED TO SPECIES MENOMICS FOR CONVIENCE
 
  EQUIVALENCE (consp(1),no2),(consp(2),no),(consp(3),roc),                  &
              (consp(4),o3),(consp(5),sgn),(consp(6),sngn)

! ASSIGN CURRENT ARRAY TO WORKING ARRAY

  DO L=1,NOSPEC
     CONSP(L)=CON1(L)
  END DO      


! EVALUATE RATE FUNCTIONS

  R(1)=K(1)*ROC+ k(1)/active1
  R(3)=K(3)*NO2
  R(4)=K(4)*NO*O3
 
! calculate concentrations of RSP
 
  rx2=k(2)*no
  rx6=k(6)*no2
  rx7=k(7)*no2
  if(k(3).gt.0)then
     bquad=rx2+rx6+rx7
     root2=bquad*bquad+4.*k(5)*r(1)
     if(root2.ge.0.0)then
        rsp=-bquad+sqrt(root2)
        rsp=rsp*0.5/k(5)
     else
        rsp=0.
     end if
  else
     rsp=0.
  endif

  R(2)=rx2*RSP
  R(5)=k(5)*RSP*RSP
  R(6)=Rx6*RSP
  R(7)=Rx7*RSP


! FORMATION (A) and LOSS TERMS (B)

! NO2
      a(1)=R(2)+R(4)
      b(1)=R(3)+R(6)+R(7)
! NO
      a(2)=R(3)
      b(2)=R(2)+R(4)
! VOC
      a(3)=0.
      b(3)=0.
! O3
      a(4)=R(3)
      b(4)=R(4)
! SGN
      a(5)=R(6)
      b(5)=0
! SNGN
      a(6)=R(7)
      b(6)=0.

  DO L=1,NOSPEC
     IF(CONSP(L).GT.0.0)B(L)=B(L)/CONSP(L)
  END DO      

END SUBROUTINE ABCALC

!========================================================

SUBROUTINE INTGRT(NOSPEC,A,B,AOB,DT,CON,CON1,ZERO)
 
! ROUTINE PERFORMS EXPONENTIAL PREDICTOR STEPCORRECTION STEP
! FOR CHEMISTRY.
 
  REAL A(7),B(7),AOB(7),CON(7),ZERO(7),CON1(7)
 
  DO I=1,NOSPEC
     T = B(I)*DT
 
!    CASE 1 B*DT<<1
     IF( T.LE.1.E-06)THEN
        CON1(I)=A(I)*DT+CON(I)
 
!    CASE 2 SMALL<B*DT<LARGE
     ELSEIF(T.LE.75.)THEN
        CON1(I)=AOB(I)+(CON(I)-AOB(I))*EXP(-T)
 
!    CASE 3 B*DT>>1
     ELSE
        CON1(I)=AOB(I)
     ENDIF

     CON1(I)=AMAX1(ZERO(I),CON1(I))

  END DO 

END SUBROUTINE INTGRT

!=======================================================

