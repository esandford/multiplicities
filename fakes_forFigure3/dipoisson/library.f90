 MODULE librarymod
 use randommod

 implicit none
 
 REAL(8), PARAMETER :: pi = 3.14159265358979323846D0
 REAL(8), PARAMETER :: third = 0.333333333D0
 REAL(8), PARAMETER :: Grv = 6.674D-11
 
 CONTAINS
 
 ! =======================================================
 REAL(8) FUNCTION Pdetect_fn(x)

 REAL(8), INTENT(IN) :: x
 REAL(8), PARAMETER :: bl = 8.060D0
 REAL(8), PARAMETER :: cl = 8.110D0
 REAL(8), PARAMETER :: dl = 0.995D0

 IF( x .EQ. 0.0 ) THEN
	 Pdetect_fn = 0.0D0
 ELSE
	 Pdetect_fn = 1.0D0 + (x/cl)**bl
	 Pdetect_fn = dl*( 1.0 - (1.0D0/Pdetect_fn) )
 END IF

 END FUNCTION Pdetect_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION Ptran_fn(p,aR)

 REAL(8), INTENT(IN) :: p, aR

 Ptran_fn = (1.0D0+p)/aR

 END FUNCTION Ptran_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION aR_fn(rho,Pdays)

 REAL(8), INTENT(IN) :: rho, Pdays
 REAL(8), PARAMETER :: days = 86400.0

 aR_fn = ( rho*Grv*(Pdays*days)**2 )/( 3.0D0*pi )
 aR_fn = aR_fn**third

 END FUNCTION aR_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION T14_fn(p,aR,b,Pdays)

 REAL(8), INTENT(IN) :: p, aR, b, Pdays
 REAL(8) :: Psecs
 REAL(8), PARAMETER :: days = 86400.0

 Psecs = Pdays*days
 IF( b .LT. (1.0D0+p) ) THEN
	 T14_fn = (Psecs/pi)*DASIN( DSQRT( ( (1.0D0+p)**2 - b**2 )/( aR**2 - b**2 ) ) )
 ELSE
	 T14_fn = 0.0D0
 END IF

 END FUNCTION T14_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION T23_fn(p,aR,b,Pdays)
 
 REAL(8), INTENT(IN) :: p, aR, b, Pdays
 REAL(8) :: Psecs
 REAL(8), PARAMETER :: days = 86400.0

 Psecs = Pdays*days
 IF( b .LT. (1.0D0-p) ) THEN
	 T23_fn = (Psecs/pi)*DASIN( DSQRT( ( (1.0D0-p)**2 - b**2 )/( aR**2 - b**2 ) ) )
 ELSE
	 T23_fn = 0.0D0
 END IF

 END FUNCTION T23_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION Wdur_fn(p,aR,b,Pdays)
 
 REAL(8), INTENT(IN) :: p, aR, b, Pdays
 REAL(8), PARAMETER :: days = 86400.0

 Wdur_fn = 0.5D0*( T14_fn(p,aR,b,Pdays) + T23_fn(p,aR,b,Pdays) )

 END FUNCTION Wdur_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION SNR_fn(p,aR,b,Pdays,sigma6)
 
 REAL(8), INTENT(IN) :: p, aR, b, Pdays, sigma6

 SNR_fn = Wdur_fn(p,aR,b,Pdays)/21600.0D0
 SNR_fn = SNR_fn / DSQRT( (T14_fn(p,aR,b,Pdays)/21600.0D0) - p*p*SNR_fn )
 SNR_fn = 1.0D6*p*p*SNR_fn/sigma6

 END FUNCTION SNR_fn
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION smallplanet_prob(a_small,a_big,logRmin,logRcrit,logRmax)
 
 REAL(8), INTENT(IN) :: a_small, a_big, logRmin, logRcrit, logRmax
 
 smallplanet_prob = -logRmin*(1.0+a_big) + logRcrit*(a_big-a_small) + logRmax*(1.0+a_small)
 smallplanet_prob = ( (logRcrit-logRmin)*(1.0+a_big) )/smallplanet_prob
 
 END FUNCTION smallplanet_prob
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION random_smallR(a_small,logRmin,logRcrit,z)
 
 REAL(8), INTENT(IN) :: a_small, logRmin, logRcrit, z

 random_smallR = z*(logRcrit - logRmin)**(1.0D0+a_small)
 random_smallR = logRmin + random_smallR**( 1.0D0/(1.0D0+a_small) )
 
 END FUNCTION random_smallR
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION random_bigR(a_big,logRcrit,logRmax,z)
 
 REAL(8), INTENT(IN) :: a_big, logRcrit, logRmax, z

 random_bigR = (1.0D0-z)**( 1.0D0/(1.0D0+a_big) )
 random_bigR = logRmax - random_bigR*(logRmax - logRcrit)
 
 END FUNCTION random_bigR
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE Qmarks(Qstart,Qend)
	
 REAL(8), DIMENSION(18), INTENT(OUT) :: Qstart, Qend
 
 Qstart(1) = 120.53881583872862D0
 Qstart(2) = 131.51205364937778D0
 Qstart(3) = 169.7649547641049D0
 Qstart(4) = 260.24513508344535D0
 Qstart(5) = 352.3975477916392D0
 Qstart(6) = 443.5107354770662D0
 Qstart(7) = 539.4700248112713D0
 Qstart(8) = 630.1954667248137D0
 Qstart(9) = 735.384193657912D0
 Qstart(10) = 808.536385360414D0
 Qstart(11) = 906.8660123426234D0
 Qstart(12) = 1001.2287755899888D0
 Qstart(13) = 1099.4293814410048D0
 Qstart(14) = 1182.7573462175933D0
 Qstart(15) = 1274.159798465902D0
 Qstart(16) = 1373.5085059477424D0
 Qstart(17) = 1472.1177699850596D0
 Qstart(18) = 1559.2464908554612D0
 
 Qend(1) = 130.24512557836715D0
 Qend(2) = 164.9833790288394D0
 Qend(3) = 258.46743138637976D0
 Qend(4) = 349.49605601447547D0
 Qend(5) = 442.20297481342277D0
 Qend(6) = 538.1622656835898D0
 Qend(7) = 629.2963867858052D0
 Qend(8) = 719.5485955661061D0
 Qend(9) = 802.3449168577863D0
 Qend(10) = 905.9260796948656D0
 Qend(11) = 1000.2684258245572D0
 Qend(12) = 1098.3259974060275D0
 Qend(13) = 1182.0217741210654D0
 Qend(14) = 1273.056388007033D0
 Qend(15) = 1371.3221660849595D0
 Qend(16) = 1471.1369960530792D0
 Qend(17) = 1557.9591908484654D0
 Qend(18) = 1591.001264764498D0
 
 END SUBROUTINE Qmarks
 ! =======================================================
 
 ! =======================================================
 INTEGER FUNCTION quarterid_fn(t)

 REAL(8), INTENT(IN) :: t
 REAL(8), DIMENSION(18) :: Qstart, Qend
 
 call Qmarks(Qstart,Qend)

 IF( t .LT. Qend(1) ) THEN
	 quarterid_fn = 1
 ELSE IF( t .LT. Qend(2) ) THEN
   quarterid_fn = 2
 ELSE IF( t .LT. Qend(3) ) THEN
   quarterid_fn = 3
 ELSE IF( t .LT. Qend(4) ) THEN
   quarterid_fn = 4
 ELSE IF( t .LT. Qend(5) ) THEN
   quarterid_fn = 5
 ELSE IF( t .LT. Qend(6) ) THEN
   quarterid_fn = 6
 ELSE IF( t .LT. Qend(7) ) THEN
   quarterid_fn = 7
 ELSE IF( t .LT. Qend(8) ) THEN
   quarterid_fn = 8
 ELSE IF( t .LT. Qend(9) ) THEN
   quarterid_fn = 9
 ELSE IF( t .LT. Qend(10) ) THEN
   quarterid_fn = 10
 ELSE IF( t .LT. Qend(11) ) THEN
   quarterid_fn = 11
 ELSE IF( t .LT. Qend(12) ) THEN
   quarterid_fn = 12
 ELSE IF( t .LT. Qend(13) ) THEN
   quarterid_fn = 13
 ELSE IF( t .LT. Qend(14) ) THEN
   quarterid_fn = 14
 ELSE IF( t .LT. Qend(15) ) THEN
   quarterid_fn = 15
 ELSE IF( t .LT. Qend(16) ) THEN
   quarterid_fn = 16
 ELSE IF( t .LT. Qend(17) ) THEN
   quarterid_fn = 17
 ELSE !IF( t .LT.  ) THEN
   quarterid_fn = 18
 END IF
 
 END FUNCTION quarterid_fn
 ! =======================================================
 
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
 
 ! =======================================================
 REAL(8) FUNCTION Delta_Hill(Mstar,radius1,radius2,period1,period2)
 
 REAL(8), INTENT(IN) :: Mstar, radius1, radius2, period1, period2
 REAL(8), DIMENSION(4) :: C, S, T
 REAL(8) :: mass1, mass2, ratio, periodx_inner, periodx_outer, mass_star
 REAL(8), PARAMETER :: twothirds = 0.666666667D0
 REAL(8), PARAMETER :: onethird = 0.333333333D0
 
 mass1 = mass_forecast(radius1)
 mass2 = mass_forecast(radius2)
 mass_star = Mstar*332978.9015405224D0 ! convert Solar masses -> Earth masses
 ratio = 2.8844991406148166D0/( (mass1+mass2)/mass_star )**onethird
 periodx_inner = MIN(period1,period2)**twothirds
 periodx_outer = MAX(period1,period2)**twothirds
 Delta_Hill = ratio*( (periodx_outer-periodx_inner)/(periodx_outer+periodx_inner) )
 
 END FUNCTION Delta_Hill
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION rayleigh(lambda,z)
 
 REAL(8), INTENT(IN) :: lambda, z
 ! 1.1107207345395915 = pi/[2sqrt(2)]
 
 rayleigh = 1.0D0 - DEXP(-(1.1107207345395915D0/lambda)**2)
 rayleigh = lambda*DSQRT( -2.0D0*DLOG(1.0D0 - z*rayleigh) )
 
 END FUNCTION rayleigh
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE simplexsort(Xsort,Rsort,length,param)
 ! Bubble sorter for the simplex component of Jammi
 ! simplexsort is passed array R (parameters) and
 ! X (chi^2 values). These are combined into array
 ! Zsort=(X,R). simplexsort then sorts the list
 ! from smallest to largest.

 implicit none
 INTEGER, INTENT(in) :: length, param
 INTEGER :: j
 REAL(8), DIMENSION(length), INTENT(inout) :: Xsort
 REAL(8), DIMENSION(param,length), INTENT(inout) :: Rsort
 INTEGER :: swaps_made, counts
 REAL(8) :: tempX
 REAL(8), DIMENSION(param) :: tempR

  !write(*,*) 'Xsort = ',Xsort
  DO ! Repeat this loop until we break out
   swaps_made=0  ! Initially, we've made no swaps
   ! Make one pass of the bubble sort algorithm
   DO counts=1,(length-1)
    ! If item is greater than the one after it, then we initiate a swap
    IF( Xsort(counts) .GT. Xsort(counts+1) ) THEN
     ! Define initial temp
     tempX = Xsort(counts)
     DO j=1,param
      tempR(j) = Rsort(j,counts)
     END DO
     ! Displace by one
     Xsort(counts) = Xsort(counts+1)
     DO j=1,param
      Rsort(j,counts) = Rsort(j,counts+1)
     END DO
     ! Finalize swap
     Xsort(counts+1) = tempX
     DO j=1,param
      Rsort(j,counts+1) = tempR(j)
     END DO
     ! Increase the swaps_made count
     swaps_made = swaps_made+1
    END IF
   END DO
   ! If no swaps,55 break loop
   IF( swaps_made == 0 ) exit
  END DO

 END SUBROUTINE simplexsort
 ! =======================================================
 
 ! =======================================================
 SUBROUTINE sorter(asort,length)
 ! Bubble sorts vector asort(length) from smallest to largest

 implicit none
 INTEGER, INTENT(in) :: length
 REAL(8), DIMENSION(length), INTENT(inout) :: asort

 INTEGER :: swaps_made, counts
 REAL(8) :: temp

  DO ! Repeat this loop until we break out
   swaps_made = 0  ! Initially, we've made no swaps
   ! Make one pass of the bubble sort algorithm
   DO counts = 1,(length-1)
    ! If item is greater than the one after it,
    ! then we initiate a swap
    IF( asort(counts) > asort(counts+1) ) THEN
     temp = asort(counts)
     asort(counts) = asort(counts+1)
     asort(counts + 1) = temp
     swaps_made = swaps_made + 1
    END IF
   END DO
   ! If no swaps, break loop
   IF( swaps_made == 0 ) exit
  END DO

 END SUBROUTINE sorter
 ! =======================================================

 ! =======================================================
 SUBROUTINE KStwo(n1,n2,x1,x2,loglike)

  INTEGER :: i
  INTEGER, INTENT(IN) :: n1
  INTEGER, INTENT(IN) :: n2
  REAL(8), DIMENSION(n1), INTENT(IN) :: x1
  REAL(8), DIMENSION(n2), INTENT(IN) :: x2
  REAL(8), DIMENSION(n1+n2) :: sb, cdf1, cdf2, KSvals
  REAL(8), DIMENSION(n1) :: s1 ! sorted version of x1
  REAL(8), DIMENSION(n2) :: s2 ! sorted version of x2
  REAL(8) :: KScrit, KS_statistic
  REAL(8), INTENT(OUT) :: loglike

  ! combine x1 and x2
  DO i=1,n1+n2
    IF( i .LE. n1 ) THEN
      sb(i) = x1(i)
    ELSE
 	   sb(i) = x2(i-n1)
    END IF
  END DO
  call sorter(sb,n1+n2)

  ! sort x1 and x2
  s1 = x1
  s2 = x2
  call sorter(s1,n1)
  call sorter(s2,n2)

  ! compute cdf's at each unique x value
  DO i=1,n1+n2
    cdf1(i) = cdf(n1,s1,sb(i))
    cdf2(i) = cdf(n2,s2,sb(i))
    KSvals(i) = DABS( cdf1(i) - cdf2(i) )
  END DO
  KS_statistic = MAXVAL(KSvals)
 
  ! If KS_statistic < KScrit, loglikes go positive, so we'll truncate above this
  ! threshold to unity
  ! 0.5887050112577373 = SQRT( log(2)/2 )
  KScrit = 0.5887050112577373D0*DSQRT( DBLE(n1+n2)/DBLE(n1*n2) )
 
  IF( KS_statistic .GT. KScrit ) THEN
    ! 0.6931471805599453 = log(2)
    loglike = ( -2.0D0*KS_statistic**2*n1*n2 + 0.6931471805599453D0*(n1+n2) )/(n1+n2)
  ELSE
 	 loglike = 0.0D0
  END IF

 END SUBROUTINE KStwo
 ! =======================================================

 ! =======================================================
 REAL(8) FUNCTION cdf(n,x,xtest)

 ! x needs to be sorted before passing to this function

  INTEGER :: i
  INTEGER, INTENT(IN) :: n
  REAL(8), DIMENSION(n), INTENT(IN) :: x
  REAL(8), INTENT(IN) :: xtest

  IF( xtest .LT. x(1) ) THEN
    cdf = 0.0D0
  ELSE IF( xtest .GE. x(n) ) THEN
    cdf = 1.0D0
  ELSE
    DO i=1,n-1
 	   IF( xtest .LT. x(i+1) .AND. xtest .GE. x(i) ) THEN
 		   cdf = DBLE(i)/DBLE(n)
 	   END IF
    END DO
  END IF

 END FUNCTION cdf
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION binomloglike(ntotobs,ntotsim,nobs,nsim)

  INTEGER, INTENT(IN) :: ntotobs, ntotsim, nobs, nsim
	REAL(8) :: psim
  
	IF( nsim .EQ. 0 .AND. nobs .EQ. 0 ) THEN
		binomloglike = 0.0D0
	ELSE
    IF( nsim .EQ. 0 ) THEN
	    psim = 0.5D0/DBLE(ntotsim)
	  ELSE
      psim = DBLE(nsim)/DBLE(ntotsim)
	  END IF
	  binomloglike = (ntotobs-nobs)*DLOG(1.0D0-psim) + nobs*DLOG(psim)
	  binomloglike = binomloglike + logbinom(ntotobs,nobs)
	END IF

 END FUNCTION binomloglike
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION logbinom(n,k)

  INTEGER, INTENT(IN) :: n, k

  logbinom = LOG_GAMMA(REAL(n+1)) - LOG_GAMMA(REAL(n-k+1)) - LOG_GAMMA(REAL(k+1))

 END FUNCTION logbinom
 ! =======================================================
 
 ! =======================================================
 REAL(8) FUNCTION poissonloglike(nobs,nsim)

  INTEGER, INTENT(IN) :: nobs, nsim
	REAL(8) :: nsimx
	
	IF( nsim .EQ. 0 .AND. nobs .EQ. 0 ) THEN
		poissonloglike = 0.0D0
	ELSE
    IF( nsim .EQ. 0 ) THEN
		  nsimx = 0.5D0
	  ELSE
		  nsimx = nsim
	  END IF
    poissonloglike = nobs*DLOG(nsimx) - nsimx - LOG_GAMMA(REAL(nobs+1))
	END IF

 END FUNCTION poissonloglike
 ! =======================================================
 
 ! =======================================================
  SUBROUTINE random_norm(mu,stdev,r)
	 
  REAL(8), INTENT(IN) :: mu, stdev
  REAL(8) :: r
  REAL(8) :: random_real
  REAL(8), PARAMETER :: roottwo = 1.4142135623730951D0

  call randomreal(.FALSE.,random_real)
  r = mu + roottwo*stdev*inverf(-1.0D0+2.0D0*random_real)

  END SUBROUTINE random_norm
 ! =======================================================

 ! =======================================================
  SUBROUTINE randomreal(firsttime,r)

  INTEGER :: count
  LOGICAL :: firsttime
  REAL(8) :: r

  IF( firsttime ) THEN
    CALL SYSTEM_CLOCK(count)
    CALL srand(count)
    r = rand()
  ELSE
    r = rand()
  END IF

  END SUBROUTINE randomreal
 ! =======================================================

 ! =======================================================
  SUBROUTINE nplanet_generator(azipf,n_max,nplanets)

  INTEGER :: j
  INTEGER, INTENT(IN) :: n_max
  REAL(8), INTENT(IN) :: azipf
  INTEGER, INTENT(OUT) :: nplanets
  REAL(8) :: norm, randomre
  REAL(8), DIMENSION(n_max) :: pdf, cdf
	
	nplanets = 11
	DO WHILE( nplanets .GT. n_max)
		nplanets = random_Poisson(REAL(azipf), .FALSE.)
	END DO

  END SUBROUTINE nplanet_generator
 ! =======================================================

 ! ======================================================================
 FUNCTION inverf(x)

  implicit none

  REAL(8) :: x
  REAL(8), PARAMETER :: awil = 0.14001228868666646D0
  REAL(8), PARAMETER :: bwil = 4.546884979448289D0
  REAL(8) :: factor, xsq, inverf

  IF( x .LT. 0.0D0 ) THEN
   factor = -1.0D0
  ELSE
   factor = 1.0D0
  END IF

  xsq = 1.0D0 - x**2
  x = bwil + 0.5D0*DLOG(xsq)
  x = DSQRT( x**2 - (DLOG(xsq)/awil) ) - x
  inverf = factor*DSQRT(x)

 END FUNCTION
 ! ======================================================================

 ! =======================================================
 REAL(8) FUNCTION circnorm(sigma,z)

 REAL(8), INTENT(IN) :: sigma, z
 ! 2.22144146907918 = pi/[sqrt(2)]

 circnorm = (2.0D0*z - 1.0D0)*ERF(2.221441469079183D0/sigma)
 circnorm = 1.4142135623730951D0*sigma*inverf(circnorm)

 END FUNCTION circnorm
 ! =======================================================

 ! =======================================================
 REAL(8) FUNCTION lognormal(sigma,z)
 ! lognormal distribution with a mode of unity

 REAL(8), INTENT(IN) :: sigma, z

 lognormal = sigma**2 - 1.4142135623730951D0*sigma*inverf(1.0D0 - 2.0D0*z)
 lognormal = DEXP(lognormal)

 END FUNCTION lognormal
 ! =======================================================

 END MODULE librarymod