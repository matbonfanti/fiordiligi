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
!>   \remark     The function Random() is taken from  \n 
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

   ! The seed of the pseudo-random series of number is stored internally in the module
   INTEGER, SAVE :: StoredSeed = 99

   ! Since the gaussian random nr are generated in couples, the following variable
   ! store the non used gaussian number for later calls of the subroutine
   REAL, SAVE    :: TempGaussian
   LOGICAL, SAVE :: GaussianAvail = .FALSE.

CONTAINS   

   SUBROUTINE SetSeed( Seed )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Seed

      StoredSeed = Seed

   END SUBROUTINE SetSeed

   !****************************************************************************

   REAL FUNCTION UniformRandomNr( X0, X1 )  RESULT( RandNr )
      IMPLICIT NONE
      REAL, INTENT(IN), OPTIONAL  :: X0, X1

      ! Generate random number in ( 0, 1 )
      RandNr = Random( StoredSeed )

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
            R     = SQRT( -2. * Sigma**2 * LOG( 1- UniformRandomNr() ) )

            ! TRansform in cartesian coordinate
            X = R * sin( Theta )
            Y = R * cos( Theta )

            ! Return one number and store the other
            RandNr = X
            TempGaussian = Y
            GaussianAvail = .TRUE.

      ELSE IF ( GaussianAvail ) THEN

            RandNr = TempGaussian
            GaussianAvail = .FALSE.

      END IF
  
   END FUNCTION GaussianRandomNr

   !****************************************************************************

   REAL FUNCTION Random(idum)
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: idum
      INTEGER, PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
      REAL, SAVE :: am
      INTEGER, SAVE :: ix=-1, iy=-1, k

      if (idum <= 0 .or. iy < 0) then 
         am = nearest(1.0,-1.0)/IM
         iy = ior(ieor(888889999,abs(idum)),1)
         ix = ieor(777755555,abs(idum))
         idum = abs(idum)+1
      end if

      ix = ieor(ix, ishft(ix,13))
      ix = ieor(ix, ishft(ix,-17))
      ix = ieor(ix, ishft(ix,5))
      k=iy/IQ
      iy=IA*(iy-k*IQ)-IR*k
      if (iy < 0) iy=iy+IM
      Random = am*ior(iand(IM,ieor(ix,iy)),1)

   END FUNCTION Random
  
END MODULE RandomNumberGenerator
