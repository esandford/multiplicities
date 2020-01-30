! =======================================================
 REAL(8) FUNCTION mass_forecast(radius)
 
 REAL(8), INTENT(IN) :: radius
 REAL(8), DIMENSION(4) :: C, S, T, Tr
 REAL(8) :: logR, logM, rando

! 10**c has units of Rearth. c1 is fit; c2-c4 are derived iteratively, 
! assuming segments meet at transition points.
 C(1) = 0.003460532109506489D0
 C(2) = -0.09252481979252211D0
 C(3) = 1.2489241453022766D0
 C(4) = -2.844669555228572D0

 ! power-law indices; unitless
 S(1) = 0.279
 S(2) = 0.589
 S(3) = -0.044
 S(4) = 0.881

 ! transition points, units of log(Mearth)
 ! T(1) = 0.3096301674258988 ! log10(2.04 Mearth)                        Terran -> Neptunian
 ! T(2) = 2.1191926778748793 ! log10(0.414 Mjup * 317.828 Mearth/Mjup)   Neptunian -> Jovian
 ! T(3) = 4.425506703276593  ! log10(0.0800 Msun * 333060.4 Mearth/Msun) Jovian -> Stellar
 ! T(4) = 0.0

 ! corresponding transition points in units of log(Rearth)
 ! computed from Jingjing's mr.Mstat2R() applied to transition points above
 Tr(1) = LOG10(1.25)  ! Terran -> Neptunian
 Tr(2) = LOG10(13.86) ! Neptunian -> Jovian
 Tr(3) = LOG10(11.59) ! Jovian -> Stellar
 
 logR = DLOG10(radius)
 ! if logR < Tr(1), unambiguously Terran   
 IF( logR .LT. Tr(1) ) THEN
   logM = ( logR - C(1) )/S(1)

 ! elif **logR < Tr(3)**, unambiguously Neptunian
 ELSE IF ( logR .LT. Tr(3) ) THEN
   logM = ( logR - C(2) )/S(2)

 ! elif **logR > Tr(2)**, unambiguously stellar
 ELSE IF ( logR .GT. Tr(2) ) THEN
   logM = ( logR - C(4) )/S(4)
 
 ELSE
 ! could be Neptunian, Jovian, or stellar. based on Monte Carlo
 ! tests from Jingjing's code, the breakdown of neptunians
 ! vs. jovians vs. stellars in this radius regime is 
 ! frac_neptunian = 0.08842
 ! frac_jovian = 0.839126
 ! frac_stellar = 0.072454
   
   call randomreal(.FALSE.,rando)

   IF ( rando .LE. 0.08842) THEN ! Neptunian
     logM = ( logR - C(2) )/S(2)
   ELSE IF ( rando .LE. 0.927546) THEN ! Jovian
     logM = ( logR - C(3) )/S(3)
   ELSE
     logM = ( logR - C(4) )/S(4)
   END IF

 END IF

 mass_forecast = 10.0D0**logM
 
 END FUNCTION mass_forecast
 ! =======================================================