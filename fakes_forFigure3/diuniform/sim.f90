PROGRAM sim
use librarymod
use modelmod

implicit none

 INTEGER :: i, j, p
 
 INTEGER, PARAMETER :: nchain = 1
 REAL(8), DIMENSION(nchain) :: loglike_c
 INTEGER :: proceed
 REAL(8) :: probaccept, randy
 REAL(8) :: r = 1.0D0

 INTEGER :: chainindex_int, IOS, chaincounter
 CHARACTER(LEN=10) :: chainindex

 CHARACTER(LEN=200) :: chainfile, trimmedchainfile, radfile, trimmedradfile, detfile, trimmeddetfile, countfile, trimmedcountfile

 REAL(8), DIMENSION(8) :: lastlinearr, tmpline
 
 ! Free parameters
 REAL(8), DIMENSION(nchain) :: asmall_c, abig_c, Rpcrit_c, sigmaR_c, azipf_c, singfrac_c, sigmainc_c
 REAL(8), PARAMETER :: asmall_d = 0.05D0
 REAL(8), PARAMETER :: abig_d = 0.08D0
 REAL(8), PARAMETER :: Rpcrit_d = 0.04D0
 REAL(8), PARAMETER :: sigmaR_d = 0.01D0
 REAL(8), PARAMETER :: azipf_d = 0.035D0
 REAL(8), PARAMETER :: singfrac_d = 0.05D0
 REAL(8), PARAMETER :: sigmainc_d = 0.10D0*0.017453292519943295D0
 REAL(8) :: asmall_t, abig_t, Rpcrit_t, sigmaR_t, azipf_t, singfrac_t, sigmainc_t, loglike_t
 
 ! Simulation length
 INTEGER, PARAMETER :: imax = 1966
 
 ! Observed Kepler counts
 INTEGER, PARAMETER :: nmax = 10
 INTEGER, PARAMETER :: lenkepcounts = 6 ! length of kepcounts_orig
 INTEGER, DIMENSION(lenkepcounts) :: kepcounts_orig
 INTEGER, DIMENSION(nmax) :: kepcounts, simcounts
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

 ! Read in chain file name and random chain index from command line
 CALL GETARG(1, chainfile)
 trimmedchainfile = trim(chainfile)

 CALL GETARG(2, chainindex)
 READ(chainindex,*) chainindex_int

 CALL GETARG(3, radfile)
 trimmedradfile = trim(radfile)

 CALL GETARG(4, detfile)
 trimmeddetfile = trim(detfile)

 CALL GETARG(5, countfile)
 trimmedcountfile = trim(countfile)

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
 
 ! We don't just want the 50th percentile of the posteriors. We want to draw sample
 ! chainindex from the posteriors.

 ! asmall_c(1) = 0.31D0 ! Pr(Rp) power-law index for small planets
 ! abig_c(1) = 4.91D0 ! Pr(Rp) power-law index for big planets
 ! Rpcrit_c(1) = 2.47D0 ! Pr(Rp) critical turn-over radius
 ! azipf_c(1) = 0.74D0 ! Zipfian slope in Pr(nplanets)
 ! sigmainc_c(1) = 1.86D0*0.017453292519943295D0 ! inclination spread in rads
 
 OPEN(unit=105,file=trimmedchainfile,access='SEQUENTIAL')

 IOS = 0
 chaincounter = 0
 DO WHILE (chaincounter.LE.chainindex_int)
   READ(UNIT=105,FMT='(G26.16,G26.16,G26.16,G26.16,G26.16,G26.16,G26.16,G26.16)',IOSTAT=IOS) tmpline
   chaincounter = chaincounter + 1
   IF (chaincounter.EQ.chainindex_int) THEN
     lastlinearr = tmpline
   END IF
 END DO

 CLOSE(105)

 !PRINT * , chainindex_int
 !PRINT * , "lastlinearr: ",lastlinearr
 
 asmall_c(1) = lastlinearr(1)    ! 0.11D0 ! Pr(Rp) power-law index for small planets
 abig_c(1) = lastlinearr(2)      ! 4.63D0 ! Pr(Rp) power-law index for big planets
 Rpcrit_c(1) = lastlinearr(3)    ! 2.38D0 ! Pr(Rp) critical turn-over radius
 sigmaR_c(1) = lastlinearr(4)    ! 0.01D0 ! spread in radii
 azipf_c(1) = lastlinearr(5)     ! 1.0D0 ! Zipfian slope in Pr(nplanets)
 singfrac_c(1) = lastlinearr(6)  ! single fraction
 sigmainc_c(1) = lastlinearr(7) * 0.017453292519943295 ! 1.0D0*0.017453292519943295D0 ! inclination spread in rads
 
 !PRINT *, "sigmaR: ", sigmaR_c(1)
 ! call the model
 call model(asmall_c(1),abig_c(1),Rpcrit_c(1),sigmaR_c(1),azipf_c(1),singfrac_c(1),sigmainc_c(1),& ! free parameters
            ! === other stuff the model needs to know ===
            imax,nmax,kepcounts,&
            nmathur,Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,QQ,&
            ikep,keplerradii,&
            logPmin,logPmax,logRpmin,logRpmax,ntransitsmin,&
            ! final values
            KSloglike,BINloglike,POIloglike,simcounts,&
            trimmedradfile,trimmeddetfile,trimmedcountfile)
 loglike_c(1) = KSloglike + POIloglike
 
 
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

OPEN(unit=101,file='/Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/mathur.dat')

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

OPEN(unit=104,file='/Users/Emily/Documents/Columbia/planetary_linguistics/multiPlanetDist_fixedForecaster/keplerradii.dat')

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
