MODULE modelmod
use librarymod

implicit none

CONTAINS

 ! =======================================================
 SUBROUTINE model(a_small,a_big,Rpcrit,sigmaR,a_zipf,sigmainc,&
                  imax,nmax,kepcounts,&
                  nmathur,Teff,logg,Rstar,Mstar,rhostar,CDPP,eCDPP,Kp,QQ,&
                  ikep,keplerradii,&
                  logPmin,logPmax,logRpmin,logRpmax,ntransitsmin,&
                  KSloglike,BINloglike,POIloglike)

 INTEGER :: i, j, k, v, vmax
 INTEGER :: proceed
 
 ! Free parameters
 REAL(8), INTENT(IN) :: a_small, a_big, Rpcrit, sigmaR, a_zipf, sigmainc
 REAL(8) :: logRpcrit
 
 ! Hardwired limits
 REAL(8), INTENT(IN) :: logPmin, logPmax, logRpmin, logRpmax
 INTEGER, INTENT(IN) :: ntransitsmin
 
 ! Simulation length
 INTEGER, INTENT(IN) :: imax
 
 ! Observed Kepler counts
 INTEGER, INTENT(IN) :: nmax
 INTEGER, DIMENSION(nmax), INTENT(IN) :: kepcounts
 INTEGER, DIMENSION(nmax) :: simcounts
 INTEGER, DIMENSION(imax) :: ndetected_table
 
 ! Stellar properties
 INTEGER, INTENT(IN) :: nmathur
 REAL(8), DIMENSION(nmathur) :: Teff, logg, Rstar, Mstar, rhostar, CDPP, eCDPP, Kp
 LOGICAL, DIMENSION(18,nmathur) :: QQ
 
 ! Observed Kepler radii
 INTEGER, INTENT(IN) :: ikep
 REAL(8), DIMENSION(ikep), INTENT(IN) :: keplerradii
 
 ! Astrophysical constants
 REAL(8), PARAMETER :: REarth = 6371.1D3
 REAL(8), PARAMETER :: RSun = 6.955D8
 REAL(8), PARAMETER :: days = 86400.0D0
 
 ! Generartive variables
 REAL(8) :: random_re, lambda_temp, firstrad, rad_temp
 INTEGER :: row, maxepoch, totquarters, quarter_temp
 REAL(8), DIMENSION(:), ALLOCATABLE :: Pdays_tab, Rp_tab, aR_tab, cosi, b_tab, p_tab
 REAL(8), DIMENSION(:,:), ALLOCATABLE :: PdaysRp_tab
 REAL(8) :: p_threshold
 REAL(8) :: Pdays, Rp, stellardensity, stellarradius, stellarmass, tauzero
 REAL(8) :: p, aR, b, sigmaval, SNRtot, Pdetect
 REAL(8), DIMENSION(18) :: noise
 REAL(8), DIMENSION(:), ALLOCATABLE :: taus
 INTEGER, DIMENSION(:), ALLOCATABLE :: obsquarters
 LOGICAL, DIMENSION(18) :: quarters
 REAL(8), DIMENSION(imax) :: Pdays_accept, Rp_accept, b_accept, aR_accept, Mstar_accept
 INTEGER, DIMENSION(imax) :: row_accept
 
 ! How many attempts do we allow for stability checks?
 INTEGER :: attempts
 INTEGER, PARAMETER :: maxattempts = 1000
 
 ! Counting planets
 INTEGER :: nplanets, ntransiters, stable, allowed, ndetected
 INTEGER :: satisfied
 
 ! Final likelihoods
 REAL(8), INTENT(OUT) :: KSloglike, BINloglike, POIloglike
 
 ! Initiate
 call randomreal(.TRUE.,random_re)

 logRpcrit = DLOG(Rpcrit)
 p_threshold = smallplanet_prob(a_small,a_big,logRpmin,logRpcrit,logRpmax)
 !write(*,*) 'args = ',a_small,a_big,logRpmin,logRpcrit,logRpmax
 !write(*,*) 'p_threshold = ',p_threshold
 
 ALLOCATE(taus(1))
 ALLOCATE(obsquarters(1))
 ALLOCATE(Pdays_tab(1))
 ALLOCATE(Rp_tab(1))
 ALLOCATE(PdaysRp_tab(2,1))
 ALLOCATE(p_tab(1))
 ALLOCATE(aR_tab(1))
 ALLOCATE(b_tab(1))
 ALLOCATE(cosi(1))

 i = 0
 v = 0
 DO WHILE (i .LE. imax)
     
     allowed = 1
     ntransiters = 0
     
     DO WHILE( ntransiters .LT. 1 )
     
       satisfied = 0
     
       DO WHILE( satisfied .EQ. 0 )
            
         ! random star
         call randomreal(.FALSE.,random_re)
         row = NINT( random_re*nmathur )

         DO WHILE ( row .EQ. 0 )
            call randomreal(.FALSE.,random_re)
            row = NINT( random_re*nmathur )
         END DO
     
         ! random Mstar
         stellarmass = Mstar(row)
         
           ! random Rstar
           stellarradius = Rstar(row)
         
           ! random rhostar
           stellardensity = 1000.0D0*rhostar(row)
    
         ! random # of planets (verified this works!)
           call nplanet_generator(a_zipf,nmax,nplanets)
         
         ! generate planets
         DEALLOCATE(Pdays_tab)
         DEALLOCATE(Rp_tab)
         DEALLOCATE(PdaysRp_tab)
         ALLOCATE(Pdays_tab(nplanets))
         ALLOCATE(Rp_tab(nplanets))
       ALLOCATE(PdaysRp_tab(2,nplanets))
         stable = 0 ! assume planets are unstable until proven otherwise
       attempts = 0
         DO WHILE( stable .EQ. 0 .AND. attempts .LT. maxattempts )
           ! keep trying to generare systems until a stable one is made
        
             attempts = attempts + 1
             
             ! a) generate a trial system (which we will later test if stable or not)
           
           ! correlated radii: draw innermost Rp from the DSPL; perturb it by sigmaR to get Rp for the others.

           ! periods
           DO j=1,nplanets ! this system has nplanets
               ! random period
               call randomreal(.FALSE.,random_re)
               Pdays_tab(j) = DEXP( logPmin + random_re*(logPmax-logPmin) )
           END DO

           ! radii
           DO j=1,nplanets ! this system has nplanets
               IF( j .EQ. 1) THEN
                 call randomreal(.FALSE.,random_re)
                 IF( random_re .LT. p_threshold ) THEN
                   call randomreal(.FALSE.,random_re)
                   firstrad = random_smallR(a_small,logRpmin,logRpcrit,random_re)
                 ELSE
                   call randomreal(.FALSE.,random_re)
                   firstrad = random_bigR(a_big,logRpcrit,logRpmax,random_re)
                 END IF
                 Rp_tab(j) = DEXP(firstrad)

               ELSE
                 proceed = 0
                 DO WHILE( proceed .EQ. 0 )
                   ! random radius. make sure it's positive
                   call random_norm(Rp_tab(1),sigmaR,rad_temp)
                   IF( rad_temp .GT. 0.0D0 ) THEN
                     proceed = 1
                   END IF
                 END DO
                 Rp_tab(j) = rad_temp
               END IF

           END DO
                 
           ! b) sort them, this makes the stability checks easier later
           PdaysRp_tab(1,:) = Pdays_tab(:)
           PdaysRp_tab(2,:) = Rp_tab(:)
           call simplexsort(Pdays_tab,PdaysRp_tab,nplanets,2)
           Pdays_tab(:) = PdaysRp_tab(1,:)
           Rp_tab(:) = PdaysRp_tab(2,:)
                 
                 ! c) great, now we override the j>1 planetary radii using "peas-in-a-pod" logic
                 !IF( nplanets .GT. 1 ) THEN
                 !  write(*,*) 'orig Rp = ',Rp_tab
                 !  DO j=2,nplanets
                 !      call randomreal(.FALSE.,random_re)
                 !      Rp_tab(j) = Rp_tab(j-1)*lognormal(sigma_peas,random_re)
                 !  END DO
                 !  write(*,*) 'new Rp = ',Rp_tab
                 !END IF
             
           ! d) now check if stable...
           IF( nplanets .EQ. 1 ) THEN
               stable = 1 ! if there's only one planet, it's always stable
           ELSE IF( nplanets .EQ. 2 ) THEN
             IF( Delta_Hill(stellarmass,Rp_tab(1),Rp_tab(2),Pdays_tab(1),Pdays_tab(2)) &
                     .GT. 3.4641016151377544D0 ) THEN
                   stable = 1 ! from Eqn~(8) of 2014ApJ...790..146F
               END IF
           ELSE
                 stable = 1 ! assume stable until proven otherwise
                   ! test 1: go through pairs, from Eqn~(8) of 2014ApJ...790..146F
                   DO j=1,nplanets-1
                 IF( Delta_Hill(stellarmass,Rp_tab(j),Rp_tab(j+1),Pdays_tab(j),Pdays_tab(j+1)) &
                         .LE. 3.4641016151377544D0 ) THEN
                       stable = 0 ! from Eqn~(8) of 2014ApJ...790..146F
                   END IF
                   END DO
                   ! test 2: go through triples, from Eqn~(9) of 2014ApJ...790..146F
               DO j=1,nplanets-2
                   IF( (Delta_Hill(stellarmass,Rp_tab(j),Rp_tab(j+1),Pdays_tab(j),Pdays_tab(j+1)) + &
                      Delta_Hill(stellarmass,Rp_tab(j+1),Rp_tab(j+2),Pdays_tab(j+1),Pdays_tab(j+2))) &
                              .LE. 18.0D0 ) THEN
                              stable = 0
                   END IF
               END DO
           END IF   
            
         END DO ! done with stability check
             
           ! loop exited, which may be because of simply exeeding maxattempts
           ! this is not satisfactory...
           IF( stable .EQ. 0 ) THEN
               satisfied = 0
                 !write(*,*) 'smashed it'
           ELSE
               satisfied = 1
           END IF
         
         END DO
         ! by this point we have generated a stable system
         ! write(*,*) 'generated stable system'
         ! write(*,*) nplanets

         ! now we need to see how many of these stable planets transit
        
         ! assign p and aR
         DEALLOCATE(aR_tab)
         DEALLOCATE(p_tab)
         ALLOCATE(aR_tab(nplanets))
         ALLOCATE(p_tab(nplanets))
         DO j=1,nplanets
           p_tab(j) = ( Rp_tab(j)*REarth )/( stellarradius*RSun )
             aR_tab(j) = aR_fn(stellardensity,Pdays_tab(j))
         END DO
        
         ! impact parameters
       DEALLOCATE(b_tab)
         DEALLOCATE(cosi)
     ALLOCATE(b_tab(nplanets))
         ALLOCATE(cosi(nplanets))
         ntransiters = 0
         DO j=1,nplanets
             ! write(*,*) j
             IF( j .EQ. 1 ) THEN
                 call randomreal(.FALSE.,cosi(j)) ! returns cosi
                 b_tab(j) = aR_tab(j)*cosi(j)
             ELSE
                 call randomreal(.FALSE.,random_re)
                 
                 ! old version, circular normal
                 ! cosi(j) = circnorm(sigmainc,random_re) ! returns -pi<delta(inc)<pi, centered on 0 

                 ! new version, Rayleigh of Rayleighs
                 lambda_temp = rayleigh(sigmainc,random_re)
                 call randomreal(.FALSE.,random_re)
                 cosi(j) = rayleigh(lambda_temp,random_re)

                 ! cos(inc_1+inc_j) = cos(inc_1)*cos(inc_j) - sin(inc_1)*sin(inc_j)
                 cosi(j) = cosi(1)*DCOS(cosi(j)) - DSQRT(1.0D0-cosi(1)**2)*DSIN(cosi(j))
                 b_tab(j) = DABS(aR_tab(j)*cosi(j))
             END IF
             IF( b_tab(j) .LT. (1.0D0+p_tab(j)) .AND. &
                 aR_tab(j) .GT. (1.0D0+p_tab(j)) ) THEN
                 ntransiters = ntransiters + 1
             END IF
         END DO ! finished cycling through nplanets
     END DO ! successfully generated a system with at least 1 transiter
     ! write(*,*) 'nplanets = ',ntransiters,' / ',nplanets
     
     ! OK, worth generating noise now...
     ! random quarters
     quarters = QQ(:,row)
     ! random noise levels
     DO j=1,18
         call random_norm( CDPP(row), eCDPP(row), noise(j) )
         !write(*,*) CDPP(random_row), eCDPP(random_row), random_noise(j)
     END DO
     
     ! let's get SNRs...
     ndetected = 0
     DO k=1,nplanets
         IF( b_tab(k) .LT. (1.0D0+p_tab(k)) .AND. &
             aR_tab(k) .GT. (1.0D0+p_tab(k)) ) THEN
             
             ! assign this planet's parameters
             Pdays = Pdays_tab(k)
             Rp = Rp_tab(k)
             b = b_tab(k)
             p = p_tab(k)
             aR = aR_tab(k)
             
             ! random tau0
             call randomreal(.FALSE.,random_re)
             tauzero = 120.53881583872862D0 + random_re*Pdays
             
             ! random max epochs
             maxepoch = CEILING( 1.0 + (1591.001264764498D0 - tauzero)/Pdays )
             DEALLOCATE(taus)
             ALLOCATE(taus(maxepoch))
             DEALLOCATE(obsquarters)
             ALLOCATE(obsquarters(maxepoch))
             
             ! random taus
             totquarters = 0
             DO j=1,maxepoch
                 taus(j) = tauzero + (j-1)*Pdays
                 quarter_temp = quarterid_fn(taus(j))
                 IF( quarters(quarter_temp) ) THEN
                     totquarters = totquarters + 1
                     obsquarters(totquarters) = quarter_temp
                 END IF
             END DO
             ! Check at least 3 transits were caught
             IF( totquarters .LT. ntransitsmin ) THEN
                 allowed = 0
             ELSE
                 allowed = 1
             END IF
             
             ! random SNR sum
             SNRtot = 0.0D0
             IF( allowed .EQ. 1 ) THEN
                 DO j=1,totquarters
                     sigmaval = noise(obsquarters(j))
                     !write(*,*) 'sigma val = ',sigmaval,' p = ',p*p*1E6
                     SNRtot = SNRtot + (SNR_fn(p,aR,b,Pdays,sigmaval))**2
                 END DO
                 SNRtot = DSQRT(SNRtot)
             END IF
             
             ! random Pdetect
             IF( allowed .EQ. 1 ) THEN
                 Pdetect = Pdetect_fn(SNRtot)
             ELSE
                 Pdetect = 0.0D0
             END IF
             
             ! ok, so bernouilli it to evaluate if "detected"
             call randomreal(.FALSE.,random_re)
             IF( random_re .LT. Pdetect ) THEN
                 ! accept...
                 i = i + 1

                 IF( i .LE. imax) THEN
                     ndetected = ndetected + 1
                     Pdays_accept(i) = Pdays
                     Rp_accept(i) = Rp
                     b_accept(i) = b
                     aR_accept(i) = aR
                     Mstar_accept(i) = Mstar(row)
                     row_accept(i) = row
                 END IF
             END IF
         END IF
     END DO
     IF ( ndetected .GE. 1 ) THEN
         v = v + 1
         ndetected_table(v) = ndetected
     END IF
 END DO
 vmax = v
 
 ! Compute KS log-like
 call KStwo(imax,ikep,Rp_accept,keplerradii,KSloglike)
 ! spit out the simulated kepler radii
 OPEN(unit=204,file='fake_keplerradii_9.dat')
 DO i=1,imax
   write(204,*) Rp_accept(i)
 END DO
 CLOSE(204)
 
 ! Compute binomial loglike
 DO i=1,nmax
     ! how many one-planet systems did you find?
     simcounts(i) = 0
     DO j=1,vmax
         IF( ndetected_table(j) .EQ. i ) THEN
             simcounts(i) = simcounts(i) + 1
         END IF
     END DO
 END DO
 
 BINloglike = 0.0D0
 DO i=1,nmax
   ! now compute loglike
   BINloglike = BINloglike + binomloglike(SUM(kepcounts),SUM(simcounts),kepcounts(i),simcounts(i))
     !write(*,*) kepcounts(i),simcounts(i),BINloglike
 END DO
 
 POIloglike = 0.0D0
 DO i=1,nmax
   ! now compute loglike
   POIloglike = POIloglike + poissonloglike(kepcounts(i),simcounts(i))
     !write(*,*) kepcounts(i),simcounts(i),POIloglike
 END DO
 ! spit out fake counts
 OPEN(unit=205,file='fake_keplercounts_9.dat')
 write(205,*) simcounts
 CLOSE(205)
 
 END SUBROUTINE model
 ! =======================================================

END MODULE modelmod
