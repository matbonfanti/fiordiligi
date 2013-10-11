!***************************************************************************************
!*                      MODULE IndependentOscillatorsModel
!***************************************************************************************
!
!>  \brief      Independent oscillator models
!>  \details    This class define the necessary subroutines to define and compute
!>              the potential and the derivatives of a discrete independent oscillator
!>              model, both in the standard and in the chain representation
!
!***************************************************************************************
!
!>  \author     Matteo Bonfanti
!>  \version    1.0
!>  \date       8 March 2013
!>
!***************************************************************************************
!
!>  \pre        To use the class the potential needs to be setup with the 
!>              couplings and the frequencies of the oscillators
!>              This data is read from input file by the setup subroutine      
!
!***************************************************************************************
!
!>  \remark     
!
!***************************************************************************************
!
!>  \todo    
!>            * Fix the units of input data for bath setup
!>            * put velocity in X and Y for the H atom
!
!***************************************************************************************

MODULE IndependentOscillatorsModel
   USE ErrorTrap
   USE MyConsts
   USE PotentialModule
   USE SplineInterpolator
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupIndepOscillatorsModel, PotentialIndepOscillatorsModel, DisposeIndepOscillatorsModel
   PUBLIC :: ZeroKelvinBathConditions, InitialBathConditions, GenericSystemAndBath
   PUBLIC :: CouplingEnergy, PotEnergyOfTheBath, GetDistorsionForce, FirstEffectiveMode
   PUBLIC :: HessianIndepOscillatorsModel

   INTEGER, PARAMETER :: STANDARD_BATH = 0        ! < normal bath, in which all the oscillators are coupled to the system
   INTEGER, PARAMETER :: CHAIN_BATH    = 1        ! < chain bath, in which all the oscillators are coupled in chain
    
   ! > Integer parameter to store the bath type (chain, standard ... )
   INTEGER :: BathType 

   ! > Number of oscillators of the bath
   INTEGER :: BathSize
   ! > Harmonic frequencies of the bath (stored in AU)
   REAL, DIMENSION(:), POINTER :: Frequencies
   ! > Coupling of the bath oscillators (stored in AU)
   REAL, DIMENSION(:), POINTER :: Couplings
   ! > Cutoff frequency of the normal bath (stored in AU)
   REAL :: CutOff
   ! > Frequency spacing of the normal bath (stored in AU)
   REAL :: DeltaOmega
   ! > Mass of the oscillators
   REAL :: OscillatorsMass
   ! > Force constant of the distorsion correction
   REAL :: DistorsionForce

   ! > Logical variable to set status of the class
   LOGICAL :: BathIsSetup = .FALSE.

   !> Max nr of iterations for potential optimization
   INTEGER, PARAMETER :: MaxIter = 10000
   !> Threshold for conjugate gradient convergence
   REAL, PARAMETER :: GradEps = 1.0E-4

CONTAINS   

! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
! SetupIndepOscillatorsModel
!*******************************************************************************
!>  This subroutine can be used to setup the parameters of the oscillator bath
!>  reading the couplings and the frequencies from input file
!> 
!> @param    FileName       Name of the file with the input parameters.
!*******************************************************************************
   SUBROUTINE SetupIndepOscillatorsModel( N, SetBathType, FileName, Mass, CutOffFreq )
      IMPLICIT NONE
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(IN)        :: SetBathType
      CHARACTER(*), INTENT(IN)   :: FileName
      REAL, INTENT(IN)           :: Mass
      REAL, INTENT(IN)           :: CutOffFreq

      INTEGER :: InpUnit, RdStatus
      INTEGER :: iBath, NData
      LOGICAL :: FileIsPresent
      REAL, DIMENSION(:), ALLOCATABLE :: RdFreq, RdSpectralDens
      TYPE(SplineType) :: SpectralDensitySpline
      REAL, DIMENSION(124) :: Coord
      REAL :: Value, D0

      ! If data already setup give a warning and deallocate memory
      CALL WARN( BathIsSetup, "IndependentOscillatorsModel.SetupIndepOscillatorsModel: overwriting bath data" )
      IF (BathIsSetup) CALL DisposeIndepOscillatorsModel( )

      ! Store number of bath degrees of freedom
      BathSize = N
      ! Set the type of bath
      BathType = SetBathType
      ! Set the mass of the oscillators
      OscillatorsMass = Mass

      ! Allocate memory
      ALLOCATE( Frequencies(BathSize), Couplings(BathSize) )

      ! Check if spectral density file exists
      INQUIRE( File = TRIM(ADJUSTL(FileName)), EXIST=FileIsPresent ) 
      CALL ERROR( .NOT. FileIsPresent,            &
                "IndependentOscillatorsModel.SetupIndepOscillatorsModel: spectral density file does not exists" )

      IF ( BathType == CHAIN_BATH ) THEN

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read force constant of the distorsion
            READ(InpUnit,*) DistorsionForce

            ! Read frequencies and couplings
            DO iBath = 1, BathSize
               READ(InpUnit,*, IOSTAT=RdStatus) Frequencies(iBath), Couplings(iBath)
               ! if data is missing, a rubin continuation of the chain is assumed 
               IF ( RdStatus /= 0 ) THEN
                  IF ( CutOffFreq > 0.0 ) THEN
                     Frequencies(iBath) = CutOffFreq / SQRT(2.)
                     Couplings(iBath) = 0.4999 * Frequencies(iBath)**2
                  ELSE 
                     Frequencies(iBath) = Frequencies(iBath-1)
                     Couplings(iBath) = Couplings(iBath-1)
                  ENDIF
               END IF
            END DO

            ! Scale with the masses
            Couplings(1) = Couplings(1) * sqrt( OscillatorsMass )
            DO iBath = 2, BathSize
               Couplings(iBath) = Couplings(iBath) * OscillatorsMass
            END DO
            D0 = Couplings(1)

            ! Close input file
            CLOSE( InpUnit )

      ELSE IF ( BathType == STANDARD_BATH ) THEN

            ! Count the number of data in the input file
            NData = CountLinesInFile( FileName )

            ! Allocate temporary arrays to store the input spectral density
            ALLOCATE( RdFreq(NData), RdSpectralDens(NData) )

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read frequencies and couplings
            DO iBath = 1, NData
               READ(InpUnit,*) RdFreq(iBath), RdSpectralDens(iBath)
               ! Transform frequencies from cm-1 to au
               RdFreq(iBath) = RdFreq(iBath) * MyConsts_cmmin1toAU
            END DO

            ! Close input file
            CLOSE( InpUnit )

            ! Set cutoff frequency
            IF ( CutOffFreq > 0.0 ) THEN
               CutOff = CutOffFreq               ! if is given from input
            ELSE 
               CutOff = RdSpectralDens( NData )  ! otherwise is the last freq in the input file
            ENDIF
            DeltaOmega = CutOff / real( BathSize)

            ! Spline interpolation of the input spectral density
            CALL SetupSpline( SpectralDensitySpline, RdFreq, RdSpectralDens)

            ! initialize distortion constant
            DistorsionForce = 0.0

            ! Compute frequencies and couplings 
            D0 = 0.0
            DO iBath = 1, BathSize
               Frequencies(iBath) = iBath * DeltaOmega
               Couplings(iBath) = SQRT( 2.0 * OscillatorsMass * Frequencies(iBath) * DeltaOmega *      & 
                    GetSpline( SpectralDensitySpline, Frequencies(iBath) ) / MyConsts_PI )
               ! Compute force constant of the distorsion
               DistorsionForce = DistorsionForce + Couplings(iBath)**2 / ( OscillatorsMass * Frequencies(iBath)**2 ) 
               D0 = D0 + Couplings(iBath)**2
            ENDDO
            D0 = sqrt(D0)

            ! Deallocate memory for spectral density interpolation
            CALL DisposeSpline ( SpectralDensitySpline )
            DEALLOCATE( RdFreq, RdSpectralDens )

      ENDIF

      DO iBath = 1, BathSize
         WRITE(888,*) Frequencies(iBath), Couplings(iBath) 
      END DO

      ! Module is setup
      BathIsSetup = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model potential has been setup"
      WRITE(*,*) " ... (details) ... "
      WRITE(*,*) " Distorsion frequency coefficient (atomic units): ", DistorsionForce
      WRITE(*,*) " D0 (atomic units): ", D0
      WRITE(*,*) " Oscillators mass (atomic units): ", OscillatorsMass
#endif

   END SUBROUTINE SetupIndepOscillatorsModel


! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
!> System potential plus bath of independent oscillators (normal or in chain form)
!> @ref RIFERIMENTI PER SURFACE OSCILLATOR MODEL
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom, 1 for 
!>                     the nearest Cu and N for the bath
!> @param Forces       Output array with the derivatives of the potential (in au)
!> @param vv           output potential in atomic units
!*******************************************************************************     
      REAL FUNCTION PotentialIndepOscillatorsModel( Positions, Forces ) RESULT(V) 
         IMPLICIT NONE

         REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
         REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 

         INTEGER :: IBath
         REAL    :: Coupl, OutOfEqZc
         REAL, DIMENSION(:), POINTER :: RH, Qbath, Dn
         REAL, POINTER :: D0
!          REAL, DIMENSION(BathSize) :: Qbath

         ! Check if the bath is setup
         CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.SystemAndIndepedentOscillators: bath is not setup" )

         ! Check the size of the position and forces arrays
         CALL ERROR( size(Positions) /= size(Forces), &
                        "IndependentOscillatorsModel.SystemAndIndepedentOscillators: array dimension mismatch" )
         CALL ERROR( size(Positions) /= BathSize+4, &
                        "IndependentOscillatorsModel.SystemAndIndepedentOscillators: wrong bath size" )

         V = 0.0
         Forces(:) = 0.0

         ! 4D POTENTIAL OF THE SYSTEM
         V = VHFourDimensional( Positions(1:4), Forces(1:4) ) 

         IF ( BathType == CHAIN_BATH ) THEN

               RH    => Positions(1:3)
               OutOfEqZc = Positions(4) - C1Puckering
               Qbath => Positions(5:BathSize+4)

               D0 => Couplings(1)
               Dn => Couplings(2:BathSize)

               ! POTENTIAL OF THE BATH OSCILLATORS PLUS COUPLING

               V = V - D0 * OutOfEqZc * Qbath(1) + 0.5 * DistorsionForce * OutOfEqZc**2
               DO iBath = 1, BathSize - 1
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Qbath(iBath) )**2 - &
                          Dn(iBath) * Qbath(iBath) * Qbath(iBath+1)
               END DO
               V = V + 0.5 * OscillatorsMass *( Frequencies(BathSize) * Qbath(BathSize) )**2

               Forces(4) = Forces(4) + Couplings(1)*Qbath(1) - DistorsionForce * OutOfEqZc
               Forces(5) = - OscillatorsMass*Frequencies(1)**2 * Qbath(1) + Couplings(1)*OutOfEqZc + Couplings(2)*Qbath(2)
               DO iBath = 2, BathSize-1
                  Forces(4+iBath) = - OscillatorsMass*Frequencies(iBath)**2 * Qbath(iBath) + &
                                    Couplings(iBath)*Qbath(iBath-1) + Couplings(iBath+1)*Qbath(iBath+1)
               END DO
               Forces(4+BathSize) = - OscillatorsMass*Frequencies(BathSize)**2 * Qbath(BathSize)  + &
                                      Couplings(BathSize) * Qbath(BathSize-1)

         ELSE IF ( BathType == STANDARD_BATH ) THEN

               OutOfEqZc = Positions(4) - C1Puckering

               ! POTENTIAL OF THE BATH OSCILLATORS PLUS COUPLING
               Coupl = 0.0
               DO iBath = 1, BathSize
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Positions(4+iBath) )**2
                  Forces(4+iBath) = - OscillatorsMass * ( Frequencies(iBath) )**2 * Positions(4+iBath)  &
                                   + OutOfEqZc * Couplings(iBath)
                  Coupl = Coupl + Couplings(iBath) * Positions(4+iBath)
               END DO
               V = V - Coupl * OutOfEqZc + 0.5*DistorsionForce*OutOfEqZc**2
               Forces(4) = Forces(4) + Coupl - DistorsionForce*OutOfEqZc

         END IF


      END FUNCTION PotentialIndepOscillatorsModel

! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
!> Generic subroutine for computing potential and its derivative of a system
!> bilinearly coupled to a bath. The coupling coordinate should be the last degree
!> of freedom of the system potential.
!>
!> @param  Positions    Array with the cartesian coordinates for the system and the bath 
!> @param  Forces       Output array with the derivatives of the potential (in au)
!> @param  NSystem      Nr of coordinates of the system
!> @result V            output potential in atomic units
!*******************************************************************************     
      REAL FUNCTION GenericSystemAndBath( Positions, Forces, SystemV, NSystem ) RESULT(V) 
         IMPLICIT NONE

         REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
         REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 
         INTERFACE
            REAL FUNCTION SystemV( X, Force )
               REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
               REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
            END FUNCTION SystemV
         END INTERFACE
         INTEGER, INTENT(IN)                     :: NSystem

         INTEGER :: IBath
         REAL    :: Coupl, CouplCoord
         REAL, DIMENSION(:), POINTER :: Qbath, Dn
         REAL, POINTER :: D0

         ! Check if the bath is setup
         CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.GenericSystemAndBath: bath is not setup" )

         ! Check the size of the position and forces arrays
         CALL ERROR( size(Positions) /= size(Forces), &
                        "IndependentOscillatorsModel.GenericSystemAndBath: array dimension mismatch" )
         CALL ERROR( size(Positions) /= BathSize+NSystem, &
                        "IndependentOscillatorsModel.GenericSystemAndBath: wrong bath size" )

         V = 0.0
         Forces(:) = 0.0

         ! 4D POTENTIAL OF THE SYSTEM
         V = SystemV( Positions(1:NSystem), Forces(1:NSystem) ) 

         IF ( BathType == CHAIN_BATH ) THEN

               CouplCoord = Positions(NSystem)
               Qbath => Positions(NSystem+1:NSystem+BathSize)
               D0 => Couplings(1)
               Dn => Couplings(2:BathSize)

               ! POTENTIAL OF THE BATH OSCILLATORS PLUS COUPLING

               V = V - D0 * CouplCoord * Qbath(1) + 0.5 * DistorsionForce * CouplCoord**2
               DO iBath = 1, BathSize - 1
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Qbath(iBath) )**2 - &
                          Dn(iBath) * Qbath(iBath) * Qbath(iBath+1)
               END DO
               V = V + 0.5 * OscillatorsMass *( Frequencies(BathSize) * Qbath(BathSize) )**2

               Forces(NSystem) = Forces(NSystem) + Couplings(1)*Qbath(1) - DistorsionForce * CouplCoord
               Forces(NSystem+1) = - OscillatorsMass*Frequencies(1)**2 * Qbath(1) + Couplings(1)*CouplCoord + Couplings(2)*Qbath(2)
               DO iBath = 2, BathSize-1
                  Forces(NSystem+iBath) = - OscillatorsMass*Frequencies(iBath)**2 * Qbath(iBath) + &
                                    Couplings(iBath)*Qbath(iBath-1) + Couplings(iBath+1)*Qbath(iBath+1)
               END DO
               Forces(NSystem+BathSize) = - OscillatorsMass*Frequencies(BathSize)**2 * Qbath(BathSize)  + &
                                      Couplings(BathSize) * Qbath(BathSize-1)

         ELSE IF ( BathType == STANDARD_BATH ) THEN

               ! POTENTIAL OF THE BATH OSCILLATORS PLUS COUPLING
               Coupl = 0.0
               DO iBath = 1, BathSize
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Positions(NSystem+iBath) )**2
                  Forces(NSystem+iBath) = - OscillatorsMass * ( Frequencies(iBath) )**2 * Positions(NSystem+iBath)  &
                                   + Positions(NSystem) * Couplings(iBath)
                  Coupl = Coupl + Couplings(iBath) * Positions(NSystem+iBath)
               END DO
               V = V - Coupl * Positions(NSystem) + 0.5*DistorsionForce*Positions(NSystem)**2
               Forces(NSystem) = Forces(NSystem) + Coupl - DistorsionForce*Positions(NSystem)

         END IF

      END FUNCTION GenericSystemAndBath

! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
!> Initial conditions for an independent oscillators model in normal or chain form
!> DETAILS ABOUT THE INITIAL CONDITIONS 
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom, 1 for 
!>                     the nearest Cu and N for the bath
!> @param Forces       Output array with the derivatives of the potential (in au)
!> @param vv           output potential in atomic units
!*******************************************************************************     
      SUBROUTINE InitialBathConditions( Positions, Velocities, Temperature )
         IMPLICIT NONE

         REAL, DIMENSION(:), INTENT(INOUT) :: Positions
         REAL, DIMENSION(:), INTENT(INOUT) :: Velocities 
         REAL, INTENT(IN) :: Temperature

         REAL :: SigmaQ, SigmaV, SigmaVH
         INTEGER :: iBath
         
         ! Check if the bath is setup
         CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.InitialBathConditions: bath not setup" )

         ! Check the size of the position and forces arrays
         CALL ERROR( size(Positions) /= size(Velocities), &
                        "IndependentOscillatorsModel.InitialBathConditions: array dimension mismatch" )
         CALL ERROR( size(Positions) /= BathSize+4, &
                        "IndependentOscillatorsModel.InitialBathConditions: wrong bath size" )

         IF ( BathType == CHAIN_BATH ) THEN

               ! Equilibrium position of H atom
               Positions(1) = 0.0000
               Positions(2) = 0.0000
               Positions(3) = 1.483 / MyConsts_Bohr2Ang
               ! Equilibrium position of C1 atom
               Positions(4) = C1Puckering
               ! Zero momentum of H and C1
               SigmaVH = sqrt( Temperature / MyConsts_mH )
               Velocities(1) = GaussianRandomNr( SigmaVH )
               Velocities(2) = GaussianRandomNr( SigmaVH )
               Velocities(3:4) = 0.0
      
               ! THE OSCILLATORS IN A CORRECT CANONICAL DISTRIBUTION FOR ZERO COUPLING
               SigmaQ = sqrt( Temperature / OscillatorsMass )
               DO iBath = 1, BathSize
                  SigmaV = SigmaQ / Frequencies(iBath)
                  Positions(4+iBath) = GaussianRandomNr( SigmaQ )
                  Velocities(4+iBath) = GaussianRandomNr( SigmaV )
               END DO

         ELSE IF ( BathType == STANDARD_BATH ) THEN

               ! Equilibrium position of H atom
               Positions(1) = 0.0000
               Positions(2) = 0.0000
               Positions(3) = 1.483 / MyConsts_Bohr2Ang
               ! Equilibrium position of C1 atom
               Positions(4) = C1Puckering
               ! Zero momentum of H and C1
               Velocities(1:4) = 0.0
      
               ! THE OSCILLATORS IN A CORRECT CANONICAL DISTRIBUTION FOR ZERO COUPLING
               SigmaQ = sqrt( Temperature / OscillatorsMass )
               DO iBath = 1, BathSize
                  SigmaV = SigmaQ / Frequencies(iBath)
                  Positions(4+iBath) = GaussianRandomNr( SigmaQ )
                  Velocities(4+iBath) = GaussianRandomNr( SigmaV )
               END DO
         
         END IF

      END SUBROUTINE InitialBathConditions


      ! Setup initial conditions for the system plus bath for a simulation of vibrational relaxation
      ! The bath is fixed in the equilibrium position with no momentum ( classical 0K )
      ! The initial position and momenta of C and H are randomly chosen among a set of 
      ! conditions which are given as input
      ! data are initialized in ATOMIC UNITS
      SUBROUTINE ZeroKelvinBathConditions( Positions, Velocities, CHInitialConditions )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
         REAL, DIMENSION(:,:), INTENT(IN) :: CHInitialConditions
         LOGICAL, DIMENSION( size(Positions ) ) :: MinimizeCheck
         INTEGER :: NRandom, iCoord, NInit
         REAL :: Value

         ! Check if the bath is setup
         CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.ZeroKelvinBathConditions: bath not setup" )

         ! Check the size of the position and forces arrays
         CALL ERROR( size(Positions) /= size(Velocities), &
                        "IndependentOscillatorsModel.ZeroKelvinBathConditions: array dimension mismatch" )
         CALL ERROR( size(Positions) /= BathSize+4, &
                        "IndependentOscillatorsModel.ZeroKelvinBathConditions: wrong bath size" )

         ! Check the nr of starting conditions given ( there should be 8 coordinates: 4 positions and 4 momenta )
         NRandom = size( CHInitialConditions, 1 )
         CALL ERROR( size( CHInitialConditions, 2 ) /= 8, &
                                          "IndependentOscillatorsModel.ZeroKelvinBathConditions: wrong number of coords " )
      
         ! Set the velocities to zero
         Velocities(:) = 0.0

         ! Set the coordinates of the bath
         IF ( BathType == CHAIN_BATH ) THEN
               Positions(1:2) = 0.0
               Positions(3) = HZEquilibrium
               Positions(4) = C1Puckering
               Positions(5:BathSize+4) = 0.0
               MinimizeCheck = (/ (.FALSE., iCoord=1,4) ,(.TRUE., iCoord=5,BathSize+4)  /)
               Value = MinimizeBathCoords( Positions, MinimizeCheck  )
         ELSE IF ( BathType == STANDARD_BATH ) THEN
               Positions(4) = C1Puckering
               DO iCoord = 1, BathSize
                  Positions(4+iCoord) = 0.0
                  Velocities(4+iCoord) = 0.0
               END DO
         END IF

         ! Choose a random initial set of coordinates
         NInit = CEILING( UniformRandomNr(0.0, real(NRandom))  )

         ! Accordingly set position and velocity
         Positions(1:4) = CHInitialConditions( NInit, 1:4 )
         Velocities(1:4) = CHInitialConditions( NInit, 5:8 )

      END SUBROUTINE ZeroKelvinBathConditions


   REAL FUNCTION MinimizeBathCoords( Coords, Mask ) RESULT( Pot )
      IMPLICIT NONE
      REAL, INTENT(INOUT), DIMENSION(:)            :: Coords
      LOGICAL, INTENT(IN), DIMENSION(size(Coords)) :: Mask

      INTEGER :: NrDimension, NrOptimization
      INTEGER :: iIter, iCoord
      REAL, DIMENSION(size(Coords)) :: Gradient
      REAL :: Norm

      ! Set dimension number
      NrDimension = size(Coords)
      ! Set optimization coordinates nr
      NrOptimization = count( Mask )
      ! Check if the nr of dimension is compatible with the slab maximum size
      CALL ERROR( NrDimension /= BathSize + 4, "IndependentOscillatorsModel.MinimizeBathCoords: wrong number of DoFs" )

      ! Cycle over steepest descent iterations
      DO iIter = 1, MaxIter

         ! compute negative of the gradient
         Pot = PotentialIndepOscillatorsModel( Coords, Gradient )

         ! compute norm of the gradient
         Norm = 0.0
         DO iCoord = 1, NrDimension
            IF ( Mask( iCoord ) ) THEN
               Norm = Norm + Gradient(iCoord)**2
            END IF
         END DO
         Norm = SQRT( Norm / NrOptimization )

         ! check convergence
         IF (Norm < GradEps) EXIT
   
         ! move geometry along gradient
         DO iCoord = 1, NrDimension
            IF ( Mask( iCoord ) ) THEN
               Coords(iCoord) = Coords(iCoord) + Gradient(iCoord)
            END IF
         END DO

      END DO

      IF ( iIter == MaxIter ) PRINT*, " NOT CONVERGED !!!!!"

   END FUNCTION MinimizeBathCoords

! ------------------------------------------------------------------------------------------------------------

   REAL FUNCTION PotEnergyOfTheBath( Coords ) RESULT( V )
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:)            :: Coords
      INTEGER :: iBath

      ! Check if the bath is setup
      CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.PotEnergyOfTheBath: bath is not setup" )
      ! Check the size of the position array
      CALL ERROR( size(Coords) /= BathSize+4, "IndependentOscillatorsModel.PotEnergyOfTheBath: wrong bath size" )

      V = 0.0

      ! COMPUTE POTENTIAL OF THE BATH OSCILLATORS
      IF ( BathType == CHAIN_BATH ) THEN

               DO iBath = 1, BathSize - 1
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Coords(4+iBath) )**2 - &
                          Couplings(iBath+1) * Coords(4+iBath) * Coords(5+iBath)
               END DO
               V = V + 0.5 * OscillatorsMass *( Frequencies(BathSize) * Coords(4+BathSize) )**2

         ELSE IF ( BathType == STANDARD_BATH ) THEN

            ! POTENTIAL OF THE BATH OSCILLATORS
            DO iBath = 1, BathSize
               V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Coords(4+iBath) )**2
            END DO

         END IF

   END FUNCTION PotEnergyOfTheBath


   REAL FUNCTION CouplingEnergy( Coords, CouplingCoordinate ) RESULT( V )
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:)  :: Coords
      REAL, INTENT(OUT), OPTIONAL     :: CouplingCoordinate
      INTEGER :: iBath
      REAL    :: Coupl

      ! Check if the bath is setup
      CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.CouplingEnergy: bath is not setup" )
      ! Check the size of the position array
      CALL ERROR( size(Coords) /= BathSize+4, "IndependentOscillatorsModel.CouplingEnergy: wrong bath size" )

      ! SYSTEM - BATH COUPLING
      IF ( BathType == CHAIN_BATH ) THEN

         V = - Couplings(1) * ( Coords(4) - C1Puckering ) * Coords(5)
         IF( PRESENT( CouplingCoordinate ) ) CouplingCoordinate = Couplings(1)*Coords(5)

      ELSE IF ( BathType == STANDARD_BATH ) THEN

         Coupl = 0.0
         DO iBath = 1, BathSize
            Coupl = Coupl + Couplings(iBath) * Coords(4+iBath)
         END DO
         V = - Coupl * ( Coords(4) - C1Puckering ) 
         IF( PRESENT( CouplingCoordinate ) ) CouplingCoordinate = Coupl

      END IF

   END FUNCTION CouplingEnergy


   REAL FUNCTION GetDistorsionForce() RESULT( Dist)
      IMPLICIT NONE
      Dist = DistorsionForce
   END FUNCTION GetDistorsionForce

   REAL FUNCTION FirstEffectiveMode( Q ) RESULT( Mode )
      IMPLICIT NONE
      REAL, DIMENSION( BathSize ), INTENT(IN) :: Q
      INTEGER :: iBath

      ! Check if the bath is setup
      CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.FirstEffectiveMode: bath is not setup" )

      ! BATH FIRST EFFECTIVE MODE
      IF ( BathType == CHAIN_BATH ) THEN
         Mode = Couplings(1)*Q(1)
      ELSE IF ( BathType == STANDARD_BATH ) THEN
         DO iBath = 1, BathSize
            Mode = Mode + Couplings(iBath) * Q(iBath)
         END DO
      END IF
   END FUNCTION

! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
!> Compute the Hessian of the potential of the bath
!>
!> @param Hessian    output matrix of size NxN, where N is the size of the bath
!*******************************************************************************     
   SUBROUTINE HessianIndepOscillatorsModel( Hessian, MassZH ) 
      IMPLICIT NONE
      REAL, DIMENSION(BathSize+4, BathSize+4), INTENT(INOUT) :: Hessian
      REAL, INTENT(IN) :: MassZH
      INTEGER :: IBath

      ! Check if the bath is setup
      CALL ERROR( .NOT. BathIsSetup, " IndependentOscillatorsModel.HessianIndepOscillatorsModel: bath is not setup" )

      ! Initialize
      Hessian(5:BathSize+4,5:BathSize+4) = 0.0
      Hessian(1:4,5:BathSize+4) = 0.0
      Hessian(5:BathSize+4,1:4) = 0.0

      IF ( BathType == CHAIN_BATH ) THEN

         ! Distorsion correction of the Zc coordinate
         Hessian(4,4) = Hessian(4,4) + DistorsionForce / MassZH
         ! Bath system coupling
         Hessian(4,5) = - Couplings(1) /  ( SQRT(MassZH) * SQRT(OscillatorsMass) )
         Hessian(5,4) = - Couplings(1) /  ( SQRT(MassZH) * SQRT(OscillatorsMass) )
         ! Diagonal elements: quadratic terms of the potential
         DO IBath = 1, BathSize
            Hessian(4+IBath,4+IBath) = Frequencies(IBath)**2
         END DO
         ! off-diagonal elements: couplings
         DO IBath = 1, BathSize-1
            Hessian(4+IBath,4+IBath+1) = -Couplings(IBath+1) / OscillatorsMass
            Hessian(4+IBath+1,4+IBath) = -Couplings(IBath+1) / OscillatorsMass
         END DO

      ELSE IF ( BathType == STANDARD_BATH ) THEN
         CALL ERROR( 1 == 1, " IndependentOscillatorsModel.HessianIndepOscillatorsModel: standard bath not yet implemented" )

      END IF

   END SUBROUTINE HessianIndepOscillatorsModel



! ------------------------------------------------------------------------------------------------------------

!*******************************************************************************
! DisposeIndepOscillatorsModel
!*******************************************************************************
!>  Deallocate memory used by this module
!*******************************************************************************
   SUBROUTINE DisposeIndepOscillatorsModel( )
      IMPLICIT NONE

      ! Do nothing if bath is not yet setup
      IF (.NOT. BathIsSetup)  RETURN

      ! deallocate memory
      DEALLOCATE( Frequencies, Couplings )

      ! Module is no longer setup
      BathIsSetup = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model data has been disposed"
#endif

   END SUBROUTINE DisposeIndepOscillatorsModel

! ------------------------------------------------------------------------------------------------------------


END MODULE IndependentOscillatorsModel
