!***************************************************************************************
!*                              MODULE RandomNumberGenerator
!***************************************************************************************
!
!>  \brief     Random number generators
!>  \details   This class defines ...
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             15 January 2013
!>
!***************************************************************************************
!
!>  \pre              
!
!***************************************************************************************
!
!>   \remark     The function ran(...) is taken from  \n 
!>               (*) Press, Teukolsky, Vetterling, Flannery \n
!>                   Numerical Recipes, the art of scientific computing, \n
!>                   Pag. 1142 Vol. 2 (FORTRAN 90) \n
!>                   Cambridge University Press \n
!>                   http://apps.nrbook.com/fortran/index.html
!
!***************************************************************************************
MODULE RandomNumberGenerator
   USE MyConsts

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetSeed, UniformRandomNr, GaussianRandomNr
   PUBLIC :: TestGaussianDistribution, TestGaussianDistribution2

   ! Integer type for the random number generation
   INTEGER, PARAMETER :: K4B=selected_int_kind(9)

   ! The seed of the pseudo-random series of number is stored internally in the module
   INTEGER(K4B), SAVE :: StoredSeed

   ! Since the gaussian random nr are generated in couples, the following variable
   ! store the non used gaussian number for later calls of the subroutine
   REAL, SAVE    :: TempGaussian
   LOGICAL, SAVE :: GaussianAvail = .FALSE.

CONTAINS   

   SUBROUTINE TestGaussianDistribution( Average, Sigma, MaxN, Seed )
      IMPLICIT NONE
      REAL, INTENT(IN) :: Average, Sigma
      INTEGER, INTENT(IN) :: MaxN
!       INTEGER(K4B), INTENT(IN), OPTIONAL :: Seed
      INTEGER, INTENT(IN), OPTIONAL :: Seed

      INTEGER :: PrintStep, iN
      REAL :: Random1, Random2 
      REAL :: AverageEst1 = 0.0, SigmaEst1 = 0.0
      REAL :: AverageEst2 = 0.0, SigmaEst2 = 0.0

      ! In case initialize seed
      IF ( Present( Seed ) )   CALL SetSeed( Seed )

      ! Set the output steps
      PrintStep = MaxN / 200

      ! Print intestation of the table
      WRITE(123,*) "# Nr of random numbers, error in Average value, error in Standard deviation "
      WRITE(124,*) "# Nr of random numbers, error in Average value, error in Standard deviation "

      AverageEst1 = 0.0
      SigmaEst1 = 0.0
      AverageEst2 = 0.0
      SigmaEst2 = 0.0
      ! Cycle over nr of number to generate
      DO iN = 1, MaxN
            ! Generate gaussian rnd number
           Random1 = Average + GaussianRandomNr( Sigma )
           Random2 = Average + GaussianRandomNr( Sigma )
!             Random = Average + UniformRandomNr( -sqrt(3.)*Sigma, sqrt(3.)*Sigma )
            ! increment sum and squared sum
            AverageEst1 = AverageEst1 + Random1
            SigmaEst1 = SigmaEst1 + Random1**2
            AverageEst2 = AverageEst2 + Random2
            SigmaEst2 = SigmaEst2 + Random2**2
            ! IF it's a printing step, print average and st dev
            IF ( mod( iN, PrintStep ) == 0 ) THEN
               WRITE(123,*) iN, AverageEst1/iN-Average, sqrt(SigmaEst1/iN-(AverageEst1/iN)**2)-Sigma
               WRITE(124,*) iN, AverageEst2/iN-Average, sqrt(SigmaEst2/iN-(AverageEst2/iN)**2)-Sigma
            ENDIF
      END DO

      WRITE(123,*) " "
      WRITE(123,*) "# Exact   ", Average, Sigma
      WRITE(124,*) " "
      WRITE(124,*) "# Exact   ", Average, Sigma

   END SUBROUTINE TestGaussianDistribution


   SUBROUTINE TestGaussianDistribution2( Sigma, MaxN, Seed )
      IMPLICIT NONE
      REAL, INTENT(IN) :: Sigma
      INTEGER, INTENT(IN) :: MaxN
      INTEGER, INTENT(IN), OPTIONAL :: Seed

      INTEGER, DIMENSION(61) :: IntCounter
      REAL    :: Random
      INTEGER :: iN, RanIndex

      WRITE(122,*) "# Interval of x, Nr of random numbers "

      ! In case initialize seed
      IF ( Present( Seed ) )   CALL SetSeed( Seed )

      IntCounter(:) = 0

      DO iN = 1, MaxN 

         Random = GaussianRandomNr( Sigma )
         RanIndex = CEILING( 30.0 + 10.*Random/Sigma )
         IF ( RanIndex >= 1 .OR. RanIndex <= 61 ) THEN
            IntCounter(RanIndex) =  IntCounter(RanIndex) + 1
         END IF

      END DO

      DO iN = 1, 61
         WRITE(122,*)  (real(iN-31)/10.)*Sigma, IntCounter(iN)
      END DO

   END SUBROUTINE TestGaussianDistribution2

   !****************************************************************************

   SUBROUTINE SetSeed( Seed )
      IMPLICIT NONE
!       INTEGER(K4B), INTENT(IN) :: Seed
      INTEGER, INTENT(IN) :: Seed
      INTEGER, DIMENSION(10) :: SeedVec
      REAL :: Temp

      SeedVec = Seed
      CALL RANDOM_SEED( put= SeedVec )

   END SUBROUTINE SetSeed

   !****************************************************************************

   REAL FUNCTION UniformRandomNr( X0, X1 )  RESULT( RandNr )
      IMPLICIT NONE
      REAL, INTENT(IN), OPTIONAL  :: X0, X1

      CALL RANDOM_NUMBER( RandNr )  

      IF ( PRESENT(X1) .AND. PRESENT(X0) ) THEN    ! both input values are present: generate nr in (X0, X1)
            RandNr = X0 + RandNr * ( X1 - X0 )
      ELSE IF ( PRESENT(X0) ) THEN    ! only x0: generate nr in (0, X0)
            RandNr = RandNr * X0
      ELSE IF ( PRESENT(X1) ) THEN    ! only x1: generate nr in (0, X1)
            RandNr = RandNr * X1
      END IF
  
   END FUNCTION UniformRandomNr

   !****************************************************************************

   ! Generate random numbers distribuited according to a gaussian distribution
   ! with standard deviation sigma ( Box-Muller algorithm ) 
   REAL FUNCTION GaussianRandomNr( Sigma ) RESULT( RandNr )
      IMPLICIT NONE
      REAL, INTENT(IN)   ::  Sigma    ! standard deviation of the gaussian distrib
      REAL :: Theta, R, X, Y

      IF ( .NOT. GaussianAvail ) THEN

            ! Generate random number for 2D gaussian function
            Theta = UniformRandomNr( 0., 2.*MyConsts_PI )
            R     = SQRT( -2. * LOG( 1- UniformRandomNr() ) )

            ! TRansform in cartesian coordinate
            X = R * sin( Theta )
            Y = R * cos( Theta )

            ! Return one number and store the other
            RandNr = Sigma * X
            TempGaussian = Sigma * Y
            GaussianAvail = .TRUE.

      ELSE IF ( GaussianAvail ) THEN

            RandNr = TempGaussian
            GaussianAvail = .FALSE.

      END IF
  
   END FUNCTION GaussianRandomNr

!    !****************************************************************************
! 
!     ! taken from Numerical Recipes for Fortran (for details, see pag 1142)
!    FUNCTION ran(idum)
!       IMPLICIT NONE
!       INTEGER(K4B), INTENT(INOUT) :: idum
!       REAL :: ran
!       INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
!       REAL, SAVE :: am
!       INTEGER(K4B), SAVE :: ix=-1, iy=-1, k
! 
!       if (idum <= 0 .or. iy < 0) then 
!          am = nearest(1.0,-1.0)/IM
!          iy = ior(ieor(888889999,abs(idum)),1)
!          ix = ieor(777755555,abs(idum))
!          idum = abs(idum)+1
!       end if
!       ix = ieor(ix, ishft(ix,13))
!       ix = ieor(ix, ishft(ix,-17))
!       ix = ieor(ix, ishft(ix,5))
!       k=iy/IQ
!       iy=IA*(iy-k*IQ)-IR*k
!       if (iy < 0) iy=iy+IM
!       ran=am*ior(iand(IM,ieor(ix,iy)),1)
! 
!    END FUNCTION ran
  
END MODULE RandomNumberGenerator
