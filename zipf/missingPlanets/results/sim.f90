PROGRAM sim
use librarymod
use modelmod

implicit none

 INTEGER :: i, j, p, v
 
 INTEGER, PARAMETER :: nchain = 1000
 INTEGER, PARAMETER :: vmax = 11 ! # of likelihood iterations
 INTEGER, PARAMETER :: vmed = 6 !1+FLOOR(vmax*0.5)
 REAL(8), DIMENSION(vmax) :: loglike_temp
 REAL(8), DIMENSION(nchain) :: loglike_c
 INTEGER :: proceed
 REAL(8) :: probaccept, randy
 
 ! Free parameters
 REAL(8), DIMENSION(nchain) :: asmall_c, abig_c, Rpcrit_c, sigmaR_c, azipf_c, sigmainc_c
 REAL(8), PARAMETER :: asmall_d = 0.051D0
 REAL(8), PARAMETER :: abig_d = 0.086D0
 REAL(8), PARAMETER :: Rpcrit_d = 0.042D0
 REAL(8), PARAMETER :: sigmaR_d = 0.01D0
 REAL(8), PARAMETER :: azipf_d = 0.038D0
 REAL(8), PARAMETER :: sigmainc_d = 0.12D0*0.017453292519943295D0
 REAL(8) :: asmall_t, abig_t, Rpcrit_t, sigmaR_t, azipf_t, sigmainc_t, loglike_t
 
 ! Simulation length
 ! number of kepler systems; runs each 11 times and takes the median (position 6)
 INTEGER, PARAMETER :: imax = 1966
 
 ! Observed Kepler counts
 INTEGER, PARAMETER :: nmax = 10
 INTEGER, PARAMETER :: lenkepcounts = 6 ! length of kepcounts_orig
 INTEGER, DIMENSION(lenkepcounts) :: kepcounts_orig
 INTEGER, DIMENSION(nmax) :: kepcounts, simcounts, gencounts, transitcounts
 INTEGER, DIMENSION(imax) :: ndetected_table
 
 ! Stellar properties
 INTEGER, PARAMETER :: nmathur = 108429
 REAL(8), DIMENSION(nmathur) :: Teff, logg, Rstar, Mstar, rhostar, CDPP, eCDPP, Kp
 LOGICAL, DIMENSION(18,nmathur) :: QQ
 
 ! Observed Kepler radii
 INTEGER, PARAMETER :: ikep = 1966
 REAL(8), DIMENSION(ikep) :: keplerradii
 
 ! Hardwired limits
 REAL(8), PARAMETER :: Rpmin = 0.50D0 ! hardwired
 REAL(8), PARAMETER :: Rpmax = 32.0D0 ! hardwired
 REAL(8), PARAMETER :: Pmin = 6.25D0  ! hardwired
 REAL(8), PARAMETER :: Pmax = 400.0D0 ! hardwired
 REAL(8) :: logPmin, logPmax, logRpmin, logRpmax
 INTEGER, PARAMETER :: ntransitsmin = 3
 
 ! Final likelihoods
 REAL(8) :: KSloglike, BINloglike, POIloglike
 
 ! File name parameters
 CHARACTER(LEN=250) :: part1, part1sim, part1gen, part1transits, part3
 CHARACTER(LEN=250) :: combinedfilename, combinedsimfilename
 CHARACTER(LEN=250) :: combinedgenfilename, combinedtransitfilename
 CHARACTER(LEN=3) :: filenamei

 ! Initiate
 call randomreal(.TRUE.,randy)
 
 ! log em
 logPmin = DLOG(Pmin)
 logPmax = DLOG(Pmax)
 logRpmin = DLOG(Rpmin)
 logRpmax = DLOG(Rpmax)
 
 ! # of 1, 2, 3,... planet systems discovered by Kepler
 kepcounts_orig = (/ 1225, 218, 76, 15, 1, 2 /)
 kepcounts(1:lenkepcounts) = kepcounts_orig(:)
 ! append on zero's as necessary
 IF( nmax .GT. lenkepcounts ) THEN
   DO i=lenkepcounts+1,nmax
     kepcounts(i) = 0
   END DO
 END IF
 
 ! get the stellar properties
 call getmathur(Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,&
                QQ(1,:),QQ(2,:),QQ(3,:),QQ(4,:),QQ(5,:),QQ(6,:),&
                QQ(7,:),QQ(8,:),QQ(9,:),QQ(10,:),QQ(11,:),QQ(12,:),&
                QQ(13,:),QQ(14,:),QQ(15,:),QQ(16,:),QQ(17,:),QQ(18,:))
 ! get the observed Kepler radii    
 call getkeplerradii(ikep,keplerradii)
 
 ! file name
 part1 = '/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/mcmc'
 part1sim = '/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/simcounts'
 part1gen = '/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/gencounts'
 part1transits = '/rigel/home/es3197/multiPlanetDist_fixedForecaster/zipf/missingPlanets/transitcounts'

 part3 = '.dat'
 CALL GETARG(1, filenamei)

 ! filenamei = 1
 combinedfilename = trim(part1) // trim(filenamei) // trim(part3)
 combinedsimfilename = trim(part1sim) // trim(filenamei) // trim(part3)
 combinedgenfilename = trim(part1gen) // trim(filenamei) // trim(part3)
 combinedtransitfilename = trim(part1transits) // trim(filenamei) // trim(part3)

 ! initialize the chain
 OPEN(unit=103,file=combinedfilename)
 OPEN(unit=104,file=combinedsimfilename)
 OPEN(unit=105,file=combinedgenfilename)
 OPEN(unit=106,file=combinedtransitfilename)

 asmall_c(1) = 0.63D0 ! BEST-FIT Pr(Rp) power-law index for small planets
 abig_c(1) = 4.69D0 !   BEST-FIT Pr(Rp) power-law index for big planets
 Rpcrit_c(1) = 2.62D0 ! BEST-FIT Pr(Rp) critical turn-over radius
 sigmaR_c(1) = 0.13D0 ! BEST-FIT radius spread in Earth radii
 azipf_c(1) = 0.80D0 !   BEST-FIT Zipfian slope in Pr(nplanets)
 sigmainc_c(1) = 2.28D0*0.017453292519943295D0 ! BEST-FIT inclination spread in rads
 
 ! make a "likelihood call"
 DO v=1,vmax
   call model(asmall_c(1),abig_c(1),Rpcrit_c(1),sigmaR_c(1),azipf_c(1),sigmainc_c(1),& ! free parameters
              ! === other stuff the model needs to know ===
              imax,nmax,kepcounts,&
              nmathur,Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,QQ,&
              ikep,keplerradii,&
              logPmin,logPmax,logRpmin,logRpmax,ntransitsmin,&
              ! final values
              KSloglike,BINloglike,POIloglike,simcounts,gencounts,&
              transitcounts)
   loglike_temp(v) = KSloglike + POIloglike
   IF( loglike_temp(v) .GT. HUGE(1.0) ) THEN
     loglike_temp(v) = HUGE(1.0)
   ELSE IF( loglike_temp(v) .LT. -HUGE(1.0) ) THEN
     loglike_temp(v) = -HUGE(1.0)
   END IF
 END DO
 call sorter(loglike_temp,vmax)
 loglike_c(1) = loglike_temp(vmed)

 j = 1 ! accepted chain indices
 !write(*,*) asmall_c(j),abig_c(j),Rpcrit_c(j),sigmaR_c(j),azipf_c(j),sigmainc_c(j)*57.29577951308232D0,loglike_c(j)
 !write(*,*) j
 write(104,*) simcounts
 write(105,*) gencounts
 write(106,*) transitcounts
 write(103,*) asmall_c(j),abig_c(j),Rpcrit_c(j),sigmaR_c(j),azipf_c(j),sigmainc_c(j)*57.29577951308232D0,loglike_c(j)

 p = 1 ! proposals to date
 DO WHILE( j .LT. nchain )
   
   ! == make a proposal ==
   
   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(asmall_c(j),asmall_d,asmall_t)
     IF( asmall_t .GT. -1.0D0 ) THEN
       proceed = 1
     END IF
   END DO
   
   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(abig_c(j),abig_d,abig_t)
     IF( abig_t .GT. -1.0D0 ) THEN
       proceed = 1
     END IF
   END DO
   
   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(Rpcrit_c(j),Rpcrit_d,Rpcrit_t)
     IF( Rpcrit_t .GT. Rpmin .AND. Rpcrit_t .LT. Rpmax ) THEN
       proceed = 1
     END IF
   END DO

   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(sigmaR_c(j),sigmaR_d,sigmaR_t)
     IF( sigmaR_t .GT. 0.0D0 ) THEN
       proceed = 1
     END IF
   END DO

   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(azipf_c(j),azipf_d,azipf_t)
     proceed = 1
   END DO
   
   proceed = 0
   DO WHILE( proceed .EQ. 0 )
     call random_norm(sigmainc_c(j),sigmainc_d,sigmainc_t)
     IF( sigmainc_t .GT. 0.0D0 ) THEN
       proceed = 1
     END IF
   END DO
   
   p = p + 1
   
   ! make a "likelihood call"
   DO v=1,vmax
     call model(asmall_t,abig_t,Rpcrit_t,sigmaR_t,azipf_t,sigmainc_t,& ! free parameters
                ! === other stuff the model needs to know ===
                imax,nmax,kepcounts,&
                nmathur,Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,QQ,&
                ikep,keplerradii,&
                logPmin,logPmax,logRpmin,logRpmax,ntransitsmin,&
                ! final values
                KSloglike,BINloglike,POIloglike,simcounts,gencounts,&
                transitcounts)
     loglike_temp(v) = KSloglike + POIloglike
     IF( loglike_temp(v) .GT. HUGE(1.0) ) THEN
       loglike_temp(v) = HUGE(1.0)
     ELSE IF( loglike_temp(v) .LT. -HUGE(1.0) ) THEN
       loglike_temp(v) = -HUGE(1.0)
     END IF
   END DO
   call sorter(loglike_temp,vmax)
   loglike_t = loglike_temp(vmed)
   
   probaccept = DEXP( loglike_t - loglike_c(j) )
   call randomreal(.FALSE.,randy)
   IF( randy .LT. probaccept ) THEN
     j = j + 1
     write(*,*) j,p
     asmall_c(j) = asmall_t
     abig_c(j) = abig_t
     Rpcrit_c(j) = Rpcrit_t
     sigmaR_c(j) = sigmaR_t
     azipf_c(j) = azipf_t
     sigmainc_c(j) = sigmainc_t
     loglike_c(j) = loglike_t
     !write(*,*) asmall_c(j),abig_c(j),Rpcrit_c(j),sigmaR_c(j),azipf_c(j),sigmainc_c(j)*57.29577951308232D0,loglike_c(j)
     write(*,*) j
     write(104,*) simcounts
     write(105,*) gencounts
     write(106,*) transitcounts
     write(103,*) asmall_c(j),abig_c(j),Rpcrit_c(j),sigmaR_c(j),azipf_c(j),sigmainc_c(j)*57.29577951308232D0,loglike_c(j)
   END IF

 END DO
 
 write(*,*) j,' accepted trials from ',p,' proposals (',100.0D0*DBLE(j)/DBLE(p),'%)'
 
 CLOSE(103)
 CLOSE(104)
 CLOSE(105)
 CLOSE(106)
 
! ======================================================
 CONTAINS
! ======================================================

! =======================================================
SUBROUTINE getmathur(Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,&
                      Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,&
                      Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17)

INTEGER :: i
INTEGER, PARAMETER :: nmathur = 108429
REAL(8), DIMENSION(nmathur), INTENT(OUT) :: Teff, logg, Rstar, Mstar, rhostar, CDPP, eCDPP, Kp
INTEGER, DIMENSION(nmathur) :: Q0i, Q1i, Q2i, Q3i, Q4i, Q5i, Q6i, Q7i, Q8i, Q9i
INTEGER, DIMENSION(nmathur) :: Q10i, Q11i, Q12i, Q13i, Q14i, Q15i, Q16i, Q17i
LOGICAL, DIMENSION(nmathur), INTENT(OUT) :: Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9
LOGICAL, DIMENSION(nmathur), INTENT(OUT) :: Q10, Q11, Q12, Q13, Q14, Q15, Q16, Q17

OPEN(unit=101,file='../../mathur.dat')

DO i = 1,nmathur
 
  READ(101,*) Teff(i),logg(i),Rstar(i),Mstar(i),rhostar(i),CDPP(i),eCDPP(i),Kp(i),&
              Q0i(i),Q1i(i),Q2i(i),Q3i(i),Q4i(i),Q5i(i),Q6i(i),Q7i(i),Q8i(i),Q9i(i),&
              Q10i(i),Q11i(i),Q12i(i),Q13i(i),Q14i(i),Q15i(i),Q16i(i),Q17i(i)

  ! quarter 0       
  IF( Q0i(i) .EQ. 1 ) THEN
    Q0(i) = .TRUE.
  ELSE
    Q0(i) = .FALSE.
  END IF
 
  ! quarter 1       
  IF( Q1i(i) .EQ. 1 ) THEN
    Q1(i) = .TRUE.
  ELSE
    Q1(i) = .FALSE.
  END IF

  ! quarter 2       
  IF( Q2i(i) .EQ. 1 ) THEN
    Q2(i) = .TRUE.
  ELSE
    Q2(i) = .FALSE.
  END IF
 
  ! quarter 3       
  IF( Q3i(i) .EQ. 1 ) THEN
    Q3(i) = .TRUE.
  ELSE
    Q3(i) = .FALSE.
  END IF
 
  ! quarter 4       
  IF( Q4i(i) .EQ. 1 ) THEN
    Q4(i) = .TRUE.
  ELSE
    Q4(i) = .FALSE.
  END IF
 
  ! quarter 5       
  IF( Q5i(i) .EQ. 1 ) THEN
    Q5(i) = .TRUE.
  ELSE
    Q5(i) = .FALSE.
  END IF
 
  ! quarter 6       
  IF( Q6i(i) .EQ. 1 ) THEN
    Q6(i) = .TRUE.
  ELSE
    Q6(i) = .FALSE.
  END IF
 
  ! quarter 7       
  IF( Q7i(i) .EQ. 1 ) THEN
    Q7(i) = .TRUE.
  ELSE
    Q7(i) = .FALSE.
  END IF
 
  ! quarter 8       
  IF( Q8i(i) .EQ. 1 ) THEN
    Q8(i) = .TRUE.
  ELSE
    Q8(i) = .FALSE.
  END IF
 
  ! quarter 9       
  IF( Q9i(i) .EQ. 1 ) THEN
    Q9(i) = .TRUE.
  ELSE
    Q9(i) = .FALSE.
  END IF
 
  ! quarter 10  
  IF( Q10i(i) .EQ. 1 ) THEN
    Q10(i) = .TRUE.
  ELSE
    Q10(i) = .FALSE.
  END IF
 
  ! quarter 11    
  IF( Q11i(i) .EQ. 1 ) THEN
    Q11(i) = .TRUE.
  ELSE
    Q11(i) = .FALSE.
  END IF
 
  ! quarter 12        
  IF( Q12i(i) .EQ. 1 ) THEN
    Q12(i) = .TRUE.
  ELSE
    Q12(i) = .FALSE.
  END IF
 
  ! quarter 13        
  IF( Q13i(i) .EQ. 1 ) THEN
    Q13(i) = .TRUE.
  ELSE
    Q13(i) = .FALSE.
  END IF
 
  ! quarter 14        
  IF( Q14i(i) .EQ. 1 ) THEN
    Q14(i) = .TRUE.
  ELSE
    Q14(i) = .FALSE.
  END IF
 
  ! quarter 15        
  IF( Q15i(i) .EQ. 1 ) THEN
    Q15(i) = .TRUE.
  ELSE
    Q15(i) = .FALSE.
  END IF
 
  ! quarter 16        
  IF( Q16i(i) .EQ. 1 ) THEN
    Q16(i) = .TRUE.
  ELSE
    Q16(i) = .FALSE.
  END IF
 
  ! quarter 17        
  IF( Q17i(i) .EQ. 1 ) THEN
    Q17(i) = .TRUE.
  ELSE
    Q17(i) = .FALSE.
  END IF
 
END DO

CLOSE(101)

END SUBROUTINE getmathur
! =======================================================

! =======================================================
SUBROUTINE getkeplerradii(nkep,keplerradii)

INTEGER :: i
INTEGER, INTENT(IN) :: nkep
REAL(8), DIMENSION(nkep), INTENT(OUT) :: keplerradii

OPEN(unit=104,file='../../keplerradii.dat')

DO i = 1,nkep
  READ(104,*) keplerradii(i)
END DO

CLOSE(104)

END SUBROUTINE getkeplerradii
! =======================================================

! =======================================================
CHARACTER(LEN=20) FUNCTION str(k)
 !   "Convert an integer to string."

 INTEGER, INTENT(IN) :: k
 write(str,*) k
 str = ADJUSTL(str)

END FUNCTION str
! =======================================================

END PROGRAM sim
