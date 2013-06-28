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
   PUBLIC :: SetupIndepOscillatorsModel, PotentialIndepOscillatorsModel, InitialBathConditions, DisposeIndepOscillatorsModel

   INTEGER, PARAMETER, PUBLIC :: STANDARD_BATH = 0        ! < normal bath, in which all the oscillators are coupled to the system
   INTEGER, PARAMETER, PUBLIC :: CHAIN_BATH    = 1        ! < chain bath, in which all the oscillators are coupled in chain
    
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
   ! > Coordinate of the slab in the minimum
   REAL, DIMENSION(120), SAVE :: MinSlab

   ! > Logical variable to set status of the class
   LOGICAL :: BathIsSetup = .FALSE.


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
      REAL, OPTIONAL, INTENT(IN) :: CutOffFreq

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
                  IF ( PRESENT( CutOffFreq )) THEN
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
            IF ( PRESENT( CutOffFreq )) THEN
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

      ! Minimize slab potential
      Coord(1:124) = 0.0
      Coord(3) = C1Puckering + 2.0
      Coord(4) = C1Puckering
      Value =  MinimizePotential( Coord, (/ (.TRUE., iBath=1,124)  /) )      
      Value = (Coord(5)+Coord(6)+Coord(7))/3.0

      ! Translate to bring C3,C4,C5 in the Z=0 plane
      DO iBath= 3,124
          Coord(iBath) = Coord(iBath) - Value
      END DO
      
      ! Store the coordinate of the slab
      MinSlab(:) = Coord(5:124)

!       PRINT*, " "
!       DO iBath = 1,8
!          WRITE(*,500) iBath+1, MinSlab(iBath)
!       END DO
!       DO iBath = 9,98
!          WRITE(*,501) iBath+1, MinSlab(iBath)
!       END DO
!       DO iBath = 99,120
!          WRITE(*,502) iBath+1, MinSlab(iBath)
!       END DO
!       500 FORMAT( " z(",I1,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       501 FORMAT( " z(",I2,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       502 FORMAT( " z(",I3,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       PRINT*, " "
!       STOP

      ! Module is setup
      BathIsSetup = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model potential has been setup"
      WRITE(*,*) " ... (details) ... "
      WRITE(*,*) " Distorsion frequency coefficient (atomic units): ", DistorsionForce
      WRITE(*,*) " D0 (atomic units): ", D0
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
         REAL, DIMENSION(124) :: Dummy
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
         V = VHSticking( (/ Positions(1:4), MinSlab(:) /), Dummy(:) ) 
         Forces(1:4) = Dummy(1:4)

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
               DO iBath = 1, BathSize-1
                  SigmaV = SigmaQ / Frequencies(iBath)
                  Positions(4+iBath) = GaussianRandomNr( SigmaQ )
                  Velocities(4+iBath) = GaussianRandomNr( SigmaV )
               END DO
               Positions(4+iBath) = GaussianRandomNr( SigmaQ )
               Velocities(4+iBath) = GaussianRandomNr( SigmaV )

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
