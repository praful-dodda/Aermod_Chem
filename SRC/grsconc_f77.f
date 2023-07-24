      SUBROUTINE GRSConc(yr, jday, hr, lat, temp, dtchem, conc0, conc1,
     &                   PollutantName, recpt, GRSIndexMap)
      IMPLICIT NONE
      INTEGER yr, jday, hr, recpt, nospec, noreac
      REAL lat, temp, dtchem, rgc, za, rk1
      CHARACTER*40 PollutantName(6)
      REAL conc0(6), conc1(6)
      INTEGER GRSIndexMap(3)
      REAL ksr(7)
      INTEGER kk, kchem, kintr

C     Constants
      PARAMETER (nospec = 6)
      PARAMETER (noreac = 7)

C     Calculate Rgc = R*TEMP/PRES assume pressures 1 atm
      rgc = 8.3144 * temp / 101.325

C     Calculate Zenith Angle
      CALL zenang(yr, jday, hr, lat, za)

C     Calculate Rate Coefficients
      CALL ierno2(za, rk1)
      CALL grates(za, temp, ksr, rk1, 1.0)

C     Define species index numbers
      kk = 0
      IF (PollutantName(GRSIndexMap(1)) .NE. 'NO2  ') kk = -1
      IF (PollutantName(GRSIndexMap(2)) .NE. 'NO   ') kk = -1
      IF (PollutantName(GRSIndexMap(3)) .NE. 'VOC  ') kk = -1
      IF (kk .LT. 0) THEN
          WRITE(*, *) '*ERROR* grscon: inappropriate application'
          WRITE(*, *) ' Necessary species not correctly indexed'
          STOP
      END IF

C     Number of chemistry time steps per meteo time step
      kchem = 1

C     Run chemistry if any non-zero species exist
      IF (SUM(conc0) .GT. 0.0) THEN

C     Change units ug/m3 to ppm
          conc0(1) = conc0(1) * rgc / 46 / 1000
          conc0(2) = conc0(2) * rgc / 30 / 1000
          conc0(3) = conc0(3) * rgc / 36 / 1000
          conc0(4) = conc0(4) * rgc / 48 / 1000

C     Main chemistry routine
          DO kintr = 1, kchem
              CALL gchem(dtchem, conc0, conc1, ksr, 1.0)
          ENDDO

C     Change units back to ug/m3
          conc1(1) = conc1(1) * 1000 * 46 / rgc
          conc1(2) = conc1(2) * 1000 * 30 / rgc
          conc1(3) = conc1(3) * 1000 * 36 / rgc
          conc1(4) = conc1(4) * 1000 * 48 / rgc

      ENDIF

      END

      SUBROUTINE GRATES(zenith, temp, k, kp, activity)

C     Constants
      REAL kp, t1
      PARAMETER (t1 = 0.00316)

C     Variables
      REAL k(7), tempi
      REAL activity, zenith, temp

C     Compute rate constants
C     Set and adjust photolysis rates
      tempi = 1.0 / temp
C     k(4) = 9.24e+05 * tempi * EXP(-1450.0 * tempi)
      k(4) = 2643.0 * EXP(-1370.0 * tempi)

      IF (zenith .LT. 90.0) THEN
C         k(1) = kp * activity * EXP(-1000.0 * 4.7 * (tempi - t1))
          k(1) = kp * 10000.0 * EXP(-4710.0 * tempi)
C         k(2) = (3.58e+06) * tempi
          k(2) = 5482.0 * EXP(242.0 * tempi)

          k(3) = kp
          k(5) = 10000.0
          k(6) = 125.0
          k(7) = k(6)
      ELSE
C         If nighttime then set photolysis rates to zero
          k(1) = 0.0
          k(2) = 0.0
          k(3) = 0.0
          k(5) = 0.0
          k(6) = 0.0
          k(7) = 0.0
      ENDIF

      END SUBROUTINE GRATES

      SUBROUTINE ZENANG(YR, JDAY, HR, OLAT, ZA)

C     Constants
      REAL, PARAMETER :: rdpdg = 0.0174532925199433   ! radians per degree
      REAL, PARAMETER :: dtr = 1.0 / rdpdg

C     Arguments
      INTEGER YR, JDAY, HR
      REAL OLAT, ZA

C     Internal variables
      REAL Dl, Eob, Njd, Lec, Lm, Gm, d, ha, ea, yr4

C     Adjust for 2-digit year representation
      IF (YR < 50) THEN
          YR4 = YR + 2000
      ELSE
          YR4 = YR + 1900
      ENDIF

C     Number of leap days since or before the year 2000 (Dl)
      IF (YR4 >= 2001) THEN
          Dl = (YR4 - 2001) / 4
      ELSE
          Dl = (YR4 - 2000) / 4 - 1
      ENDIF

C     Obliquity of the ecliptic (Eob)
      Njd = 364.5 + (YR4 - 2001) * 365 + Dl + JDAY
      Eob = (23.439 - 0.0000004 * Njd) * rdpdg

C     Ecliptic Longitude if the Sun
      Lm = (280.460 + 0.9856474 * Njd) * rdpdg
      Gm = (357.528 + 0.9856003 * Njd) * rdpdg
      Lec = Lm + 1.915 * SIN(Gm) + 0.020 * SIN(2 * Gm)

C     Sine of Declination Angle
      d = ASIN(SIN(Eob) * SIN(Lec))

C     Local Hour Angle
      ha = (HR - 12) * 15.0 * dtr

C     Zenith angle in degrees
      ZA = ACOS(SIN(OLAT * rdpdg) * SIN(d) + COS(ha) * COS(OLAT * rdpdg) * COS(d)) / rdpdg

C     Elevation Angle    
      ea = 90.0 - ZA

      END SUBROUTINE ZENANG

      SUBROUTINE IERNO2(ZAG, RK1)

C     Constants
      REAL, PARAMETER :: rpdg = 0.01745329    ! radians per degree
      REAL, PARAMETER :: solc = 1104.0       ! cloud top solar constant (w/m2)
        
C     Arguments
      REAL ZAG, RK1

C     Internal variables
      REAL zar, sea, swf

C     Sine of the elevation angle (rad)
      sea = SIN((90.0 - ZAG) * rpdg)

C     Zenith angle in radians from elevation angle
      zar = ZAG * rpdg

C     Calculate SWF: short-wave flux (w/m2) 
C     Solar Radiation at the Earth's Surface assuming clear skies
C     SWF = SEA * SOLC * T
C     where SEA is the sine of the elevation angle in radians, SOLC is the solar 
C     constant incident at the top of the cloud layer and T is the fraction transmitted
C     through the clouds (assumed as 1). Using the normal top-of-the-atmosphere solar constant 
C     to be about 1380 W m-2 and assuming that on average about 20% of the radiation is absorbed 
C     or reflected back into space (Sellers, 1972) by a clear atmosphere, then S = 1100 W m-2. 
      swf = MAX(0.0, sea * solc)

C     Rate coefficient for NO2 photolysis
      IF (ZAG .GE. 0.0 .AND. ZAG .LT. 47.0) THEN
          RK1 = (4.23e-04 + (1.09e-04 / COS(zar))) * swf
      ELSEIF (ZAG .GE. 47.0 .AND. ZAG .LT. 64.0) THEN
          RK1 = 5.82e-04 * swf
      ELSEIF (ZAG .GE. 64.0 .AND. ZAG .LE. 90.0) THEN
          RK1 = (-0.997e-04 + 1.2e-03 * (1.0 - COS(zar))) * swf
      ELSE
          RK1 = 0.0
      ENDIF

      END SUBROUTINE IERNO2
