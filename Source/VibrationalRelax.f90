!***************************************************************************************
!*                              MODULE VibrationalRelax
!***************************************************************************************
!
!>  \brief     Subroutines for the simulations of vibrational relaxation
!>  \details   This module contains subroutines to compute averages  \n
!>             in a vibrational relaxation simulation                \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             4 October 2013
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 23 October 2013: double chain bath has been implemented
!
!>  \todo         Implement energy flux along the chain
!>  \todo         The cluster bath has some problems... needs debugging
!>  |todo         Initial conditions for non collinear calculations should be defined
!>                 
!***************************************************************************************
MODULE VibrationalRelax
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: VibrationalRelax_ReadInput, VibrationalRelax_Initialize, VibrationalRelax_Run, VibrationalRelax_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

   !> Variables to set a model 2D harmonic potential for the system
   LOGICAL :: ModelHarmonic                             !< Model potential is enabled
   REAL    :: OmegaC, OmegaH, OmegaBend, CouplingCH     !< Parameters of the model potential
   REAL, DIMENSION(4,4) :: HessianMatrix                !< Hessian matrix of the potential in the minimum

   ! Initial conditions
   REAL    :: InitEnergy                !< Initial energy of the system
   INTEGER :: NInitNormalMode           !< Initial vibrational normal mode
   REAL    :: TimeBetweenSnaps          !< Nr of initial snapshots to randomize initial conditions
   INTEGER :: NrOfInitSnapshots         !< Time between each initial snapshot

   ! Initial equilibration
   REAL    :: Temperature               !< Temperature of the simulation
   REAL    :: EquilTStep                !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma                !< Friction parameter of the Langevin equation
   INTEGER :: NrEquilibSteps            !< Nr of time step of the equilibration
   INTEGER :: EquilibrationStepInterval !< Nr of time steps between each printing of equilibration info ( = NrEquilibSteps / NrOfPrintSteps )

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of time steps between each printing of propagation info ( = NrOfSteps / NrOfPrintSteps )

   ! Time evolution dataset
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results
   TYPE(Evolution),SAVE :: InitialConditions     !< Propagate in microcanonical ensamble to generate initial conditions
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in macrocanonical ensamble at given T to generate init conditions

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageEBath       !< Average energy of the bath vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageECoup       !< Average coupling energy vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageESys        !< Average energy of the system vs time
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord       !< Average i-th coordinate vs time
   REAL, DIMENSION(:,:), ALLOCATABLE    :: OscillCorrelations !< 1-n oscillator correlations in the chain

   TYPE(RNGInternalState), SAVE :: RandomNr

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific to the vibrational
!> relaxation problem.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE VibrationalRelax_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! Set model harmonic potential
      IF ( BathType /=  SLAB_POTENTIAL ) THEN
         CALL SetFieldFromInput( InputData, "ModelHarmonic", ModelHarmonic, .FALSE. )
      END IF
      IF ( ModelHarmonic ) THEN
         CALL SetFieldFromInput( InputData, "OmegaC",     OmegaC     )
         OmegaC = OmegaC * FreqConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "OmegaH",     OmegaH     )
         OmegaH = OmegaH * FreqConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "OmegaBend",  OmegaBend  )
         OmegaBend = OmegaBend * FreqConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "CouplingCH", CouplingCH )
         CouplingCH = CouplingCH * (FreqConversion(InputUnits, InternalUnits)**2)
      END IF

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Normal mode to be initially excited
      CALL SetFieldFromInput( InputData, "NInitNormalMode", NInitNormalMode )

      ! Initial energy of the system (input in eV and transform to AU)
      CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
      InitEnergy = InitEnergy * EnergyConversion(InputUnits, InternalUnits)

      ! Nr of initial snapshots
      CALL SetFieldFromInput( InputData, "NrOfInitSnapshots", NrOfInitSnapshots )

      ! Time between each snapshots (input in fs and transform to AU)
      CALL SetFieldFromInput( InputData, "TimeBetweenSnaps", TimeBetweenSnaps)
      TimeBetweenSnaps = TimeBetweenSnaps * TimeConversion(InputUnits, InternalUnits) 

      ! READ THE VARIABLES FOR THE RELAXATION SIMULATION 

      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )

      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      TimeStep = TimeStep * TimeConversion(InputUnits, InternalUnits) 

      ! Nr of steps of the propagation
      CALL SetFieldFromInput( InputData, "NrOfSteps", NrOfSteps )

      ! Nr of steps in which the output is written 
      CALL SetFieldFromInput( InputData, "NrOfPrintSteps", NrOfPrintSteps )
      ! Accordingly set the interval between each printing step
      PrintStepInterval = NrOfSteps / NrOfPrintSteps

      IF ( BathType ==  SLAB_POTENTIAL ) THEN

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

         ! Temperature of the equilibrium simulation
         CALL SetFieldFromInput( InputData, "Temperature", Temperature )
         Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)

         ! Set gamma of the equilibration Langevin dynamics
         CALL SetFieldFromInput( InputData, "EquilRelaxTime", EquilGamma)
         EquilGamma = 1. / ( EquilGamma * TimeConversion(InputUnits, InternalUnits) )

         ! Set the time step of the equilibration
         CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
         EquilTStep = EquilTStep * TimeConversion(InputUnits, InternalUnits)

         ! Set nr of steps of the equilibration
         CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(10.0*(1.0/EquilGamma)/EquilTStep) )
         EquilibrationStepInterval = CEILING( real(NrEquilibSteps) / real(NrOfPrintSteps) )

      END IF

      WRITE(*, 902) InitEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),                 &
                    (InitEnergy-MinimumEnergy)*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                    NrOfInitSnapshots, TimeBetweenSnaps*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)

      WRITE(*, 903) NrTrajs, TimeStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), NrOfSteps, NrOfPrintSteps

      IF ( BathType ==  SLAB_POTENTIAL ) &
         WRITE(*, 904) 1./EquilGamma*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrEquilibSteps*EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrEquilibSteps

   902 FORMAT(" * Initial conditions of the atom-surface system ", /,&
              " * Absolute initial energy:                     ",F10.4,1X,A,/,&  
              " *  - w.r.t. the bottom of the well:            ",F10.4,1X,A,/,&  
              " * Nr of initial system snapshots:              ",I10,       /,& 
              " * Time between snapshots:                      ",F10.4,1X,A,/ )

   903 FORMAT(" * Dynamical simulation variables               ",           /,&
              " * Nr of trajectories:                          ",I10,       /,& 
              " * Propagation time step:                       ",F10.4,1X,A,/,&  
              " * Nr of time steps of each trajectory:         ",I10,       /,& 
              " * Nr of print steps of each trajectory:        ",I10,       / )

   904 FORMAT(" * Bath equilibration variables                 ",           /,&
              " * Relaxation time of the Langevin dynamics:    ",F10.4,1X,A,/,& 
              " * Equilibration time step:                     ",F10.4,1X,A,/,& 
              " * Equilibration total time:                    ",F10.4,1X,A,/,&
              " * Nr of equilibration steps:                   ",I10,       / )


   END SUBROUTINE VibrationalRelax_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         ALLOCATE( X(3+NCarbon), V(3+NCarbon), A(3+NCarbon), MassVector(3+NCarbon), LangevinSwitchOn(3+NCarbon) )
         MassVector = (/ (MassH, iCoord=1,3), (MassC, iCoord=1,NCarbon) /)

      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ALLOCATE( X(4+NBath), V(4+NBath), A(4+NBath), MassVector(4+NBath), LangevinSwitchOn(4+NBath) )
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,NBath) /)

      ELSE IF ( BathType ==  DOUBLE_CHAIN ) THEN
         ALLOCATE( X(4+2*NBath), V(4+2*NBath), A(4+2*NBath), MassVector(4+2*NBath), LangevinSwitchOn(4+2*NBath) )
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,2*NBath) /)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ALLOCATE( X(4), V(4), A(4), MassVector(4), LangevinSwitchOn(4) )
         MassVector = (/ (MassH, iCoord=1,3), MassC /)
      END IF

      ! store the nr of dimensions
      NDim = size(X)

      ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector, TimeStep )

      ! Set canonical dynamics at the borders of the carbon slab
      IF ( BathType == SLAB_POTENTIAL .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .TRUE.
         LangevinSwitchOn( 1: MIN( 73, NCarbon )+3 ) = .FALSE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Set canonical dynamics at the end of the oscillator chain or chains
      IF ( BathType == CHAIN_BATH .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( NDim ) = .TRUE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF
      IF ( BathType == DOUBLE_CHAIN .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( 4+NBath ) = .TRUE.  ! end of the first chain
         LangevinSwitchOn( NDim )    = .TRUE.  ! end of the second chain
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Switch on Langevin relaxation when needed
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration of the 4D system only in the microcanonical ensamble 
      CALL EvolutionSetup( InitialConditions, 4, MassVector(1:4), TimeStep )
      
      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
         CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
         LangevinSwitchOn = (/ (.FALSE., iCoord=1,3), (.TRUE., iCoord=4,NDim ) /)
         CALL SetupThermostat( Equilibration, EquilGamma, Temperature, LangevinSwitchOn )
      END IF

      DEALLOCATE( LangevinSwitchOn )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the average energies over time only
      ! DEBUG - print also the average coordinates, the power spectrum of the autocorrelation function, 
      !         the oscillators correlations for a linear bath
      ! FULL - print detailed information about the initial conditions and about the trajectories

      ! Allocate and initialize the variables for the trajectory averages
      ALLOCATE( AverageESys(0:NrOfPrintSteps), AverageEBath(0:NrOfPrintSteps), AverageECoup(0:NrOfPrintSteps) )
      AverageESys(:)            = 0.0
      AverageEBath(:)           = 0.0
      AverageECoup(:)           = 0.0

      ! Average coordinates over time and power spectrum of the auto-correlation function
      IF ( PrintType >= FULL ) THEN
         ALLOCATE( AverageCoord(NDim,0:NrOfPrintSteps) )
         AverageCoord(1:NDim,0:NrOfPrintSteps) = 0.0
      END IF

      ! Correlations in the chain bath
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) THEN
         ALLOCATE( OscillCorrelations(NDim-4, 0:NrOfPrintSteps) )
         OscillCorrelations(:,:) = 0.0
      END IF

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

      ! Initialize hessian matrix for the model potential
      IF ( ModelHarmonic ) THEN
         HessianMatrix(:,:) = 0.0
         HessianMatrix(1,1) = OmegaBend**2
         HessianMatrix(2,2) = OmegaBend**2
         HessianMatrix(3,3) = OmegaH**2
         HessianMatrix(4,4) = OmegaC**2
         HessianMatrix(3,4) = CouplingCH
         HessianMatrix(4,3) = CouplingCH
      END IF

   END SUBROUTINE VibrationalRelax_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Run()
      IMPLICIT NONE
      INTEGER  ::  AvEnergyOutputUnit, AvCoordOutputUnit       ! UNITs FOR OUTPUT AND DEBUG
      INTEGER  ::  AvBathCoordUnit, OscillCorrUnit, InitialCondUnit, NormalModesUnit
      INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel, DebugUnitNormalModes
      REAL     ::  VSys, KSys, Ecoup, VBath, KBath               ! ISTANTANEOUS ENERGY VALUES
      REAL     ::  TotEnergy, PotEnergy, KinEnergy, IstTemperature
      REAL     ::  TempAverage, TempVariance
      INTEGER  ::  NTimeStepEachSnap, iCoord, iTraj, iStep,  iOmega, kStep, NInit
      CHARACTER(100) :: OutFileName
      REAL, DIMENSION(NrOfInitSnapshots,8) :: CHInitConditions
      REAL, DIMENSION(:), ALLOCATABLE      :: AverCoordOverTime
      REAL, DIMENSION(4) :: NormalModes

      AvEnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageEnergy.dat", UNIT=AvEnergyOutputUnit )
      WRITE(AvEnergyOutputUnit, "(A,I6,A,/)") "# average energy vs time (fs | eV) - ", NrTrajs, " trajectories "

      IF ( PrintType >= FULL ) THEN
         AvCoordOutputUnit = LookForFreeUnit()
         OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
         WRITE(AvCoordOutputUnit, "(A,I6,A,/)") "# average coordinate vs time (fs | Ang) - ", NrTrajs, " trajectories "

         AvBathCoordUnit = LookForFreeUnit()
         OPEN( FILE="AverageBathCoords.dat", UNIT=AvBathCoordUnit )
         WRITE(AvBathCoordUnit, "(A,I6,A,/)") "# average Q coordinate vs time (Ang | fs) - ", NrTrajs, " trajectories "
      ENDIF

      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) THEN
         OscillCorrUnit = LookForFreeUnit()
         OPEN( FILE="OscillCorrelations.dat", UNIT=OscillCorrUnit )
         WRITE(OscillCorrUnit, "(A,I6,A,/)") "# Correlations of oscillators - ", NrTrajs, " trajectories (fs | au)"
      END IF

      IF ( PrintType == DEBUG ) THEN
         InitialCondUnit = LookForFreeUnit()
         OPEN( FILE="InitialConditionsCH.dat", UNIT=InitialCondUnit )
         WRITE(InitialCondUnit, "(A,I6,A,/)") "# Initial Conditions of the system: X_H, Y_H, Z_H, Z_C and velocities (au)"

         NormalModesUnit = LookForFreeUnit()
         OPEN( FILE="NormalModesSystem.dat", UNIT=NormalModesUnit )
         WRITE(NormalModesUnit, "(A,I6,A,/)") "# normal modes of the system vs time: Q1, Q2, Q3 and Q4 (au)"
      ENDIF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      ! Define the set of initial conditions for CH
      ! atoms in the equilibrium geometry, energy in the given normal mode

      WRITE(*,700) SQRT(NormalModes4D_Freq(NInitNormalMode))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   NormalModes4D_Vecs(1,NInitNormalMode)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                   NormalModes4D_Vecs(2,NInitNormalMode)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                   NormalModes4D_Vecs(3,NInitNormalMode)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                   NormalModes4D_Vecs(4,NInitNormalMode)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)

      X(1:2) = 0.0
      X(3) = HZEquilibrium 
      X(4) = C1Puckering  
      V(1) = sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassH ) * NormalModes4D_Vecs(1,NInitNormalMode)
      V(2) = sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassH ) * NormalModes4D_Vecs(2,NInitNormalMode)
      V(3) = sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassH ) * NormalModes4D_Vecs(3,NInitNormalMode)
      V(4) = sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassC ) * NormalModes4D_Vecs(4,NInitNormalMode)

      ! Define when to store the dynamics snapshot
      NTimeStepEachSnap = INT( TimeBetweenSnaps / TimeStep )

      ! Compute starting potential and forces
      A(:) = 0.0
      PotEnergy = VHFourDimensional( X(1:4), A(1:4) )
      DO iCoord = 1,4
         A(iCoord) = A(iCoord) / MassVector(iCoord)
      END DO

      ! Cycle over the nr of snapshot to store
      DO iTraj = 1, NrOfInitSnapshots

         ! Propagate the 4D traj in the microcanonical ensamble
         DO iStep = 1, NTimeStepEachSnap
            ! Propagate for one timestep with Velocity-Verlet
            CALL EOM_LangevinSecondOrder( InitialConditions, X(1:4), V(1:4), A(1:4), VHFourDimensional, PotEnergy, RandomNr )
            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy(InitialConditions, V(1:4) )
            TotEnergy = PotEnergy + KinEnergy

         END DO

         ! Store snapshot
         CHInitConditions( iTraj, 1:4 ) = X(1:4)
         CHInitConditions( iTraj, 5:8 ) = V(1:4)
      END DO

      ! Print snapshots if there debug printing is required
      IF ( PrintType == DEBUG ) THEN
         DO iTraj = 1, NrOfInitSnapshots
            WRITE( InitialCondUnit, "(8F15.5)" ) CHInitConditions( iTraj, : )
         END DO
      ENDIF

      !run NrTrajs number of trajectories
      DO iTraj = 1,NrTrajs
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions
         IF ( BathType == SLAB_POTENTIAL ) THEN 
               CALL ZeroKelvinSlabConditions( X, V, CHInitConditions, RandomNr ) 
               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
         ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
               NInit = CEILING( UniformRandomNr(RandomNr)*real(NrOfInitSnapshots) )
               X(1:4) = CHInitConditions( NInit, 1:4 )
               V(1:4) = CHInitConditions( NInit, 5:8 )
               CALL ZeroKelvinBathConditions( Bath, X(5:), V(5:), ZPECorrection, RandomNr )
         ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
               NInit = CEILING( UniformRandomNr(RandomNr)*real(NrOfInitSnapshots) )
               X(1:4) = CHInitConditions( NInit, 1:4 )
               V(1:4) = CHInitConditions( NInit, 5:8 )
               CALL ZeroKelvinBathConditions( DblBath(1), X(5:NBath+4), V(5:NBath+4), ZPECorrection, RandomNr )
               CALL ZeroKelvinBathConditions( DblBath(2), X(NBath+5:2*NBath+4), V(NBath+5:2*NBath+4), ZPECorrection, RandomNr )
         ELSE IF ( BathType == LANGEVIN_DYN ) THEN
               NInit = CEILING( UniformRandomNr(RandomNr)*real(NrOfInitSnapshots) )
               X(1:4) = CHInitConditions( NInit, 1:4 )
               V(1:4) = CHInitConditions( NInit, 5:8 )
         END IF

         IF ( BathType == SLAB_POTENTIAL ) THEN 

            PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", Temperature / MyConsts_K2AU

            ! Initialize temperature average and variance
            TempAverage = 0.0
            TempVariance = 0.0

            ! Compute starting potential and forces
            A(:) = 0.0
            PotEnergy = VibrRelaxPotential( X, A )
            A(:) = A(:) / MassVector(:)

            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy(Equilibration, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

            ! Do an equilibration run
            EquilibrationCycle: DO iStep = 1, NrEquilibSteps

               ! PROPAGATION for ONE TIME STEP
               CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, VibrRelaxPotential, PotEnergy, RandomNr )

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

               ! store temperature average and variance
               TempAverage = TempAverage + IstTemperature
               TempVariance = TempVariance + IstTemperature**2

            END DO EquilibrationCycle

            PRINT "(A)", " Equilibration completed! "

            X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
            NInit = CEILING( UniformRandomNr(RandomNr)*real(NrOfInitSnapshots) )
            X(1:4) = CHInitConditions( NInit, 1:4 )
            V(1:4) = CHInitConditions( NInit, 5:8 )

         END IF 

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Energy of the system, Coupling energy and energy of the bath
         CALL IstantaneousEnergies( PotEnergy, KSys, VSys, Ecoup, VBath, KBath, iTraj, iStep )
         IF ( NDim > 4 ) THEN
            IstTemperature = 2.0*KBath/(MyConsts_K2AU*(NDim-4))
         ELSE
            IstTemperature = 0.0
         ENDIF

         ! Store starting values of the averages
         AverageESys(0)         = AverageESys(0)  + KSys + VSys
         AverageEBath(0)        = AverageEBath(0) + KBath + VBath
         AverageECoup(0)        = AverageECoup(0) + Ecoup
         IF ( PrintType >= FULL )  AverageCoord(1:NDim,0) = AverageCoord(1:NDim,0) + X(:)

         IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) THEN
            DO iCoord = 1, NDim-4
               OscillCorrelations(iCoord, 0) = OscillCorrelations(iCoord, 0) + X(4) * X(4+iCoord) 
            END DO
         END IF

         ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
         WRITE(*,600)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature

        ! Open unit for massive output, with detailed info on trajectories
         IF ( PrintType == DEBUG ) THEN
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Energy.dat"
            DebugUnitEn = LookForFreeUnit()
            OPEN( Unit=DebugUnitEn, File=OutFileName )

            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Coord.dat"
            DebugUnitCoord = LookForFreeUnit()
            OPEN( Unit=DebugUnitCoord, File=OutFileName )

            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Vel.dat"
            DebugUnitVel = LookForFreeUnit()
            OPEN( Unit=DebugUnitVel, File=OutFileName )

            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_NormalModes.dat"
            DebugUnitNormalModes = LookForFreeUnit()
            OPEN( Unit=DebugUnitNormalModes, File=OutFileName )

            ! Write initial values
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | ESys, ECoup, EBath, KSys, VSys, KBath, VBath / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, X(1), X(2), X(3), X(4), (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)

            WRITE( DebugUnitNormalModes, "(/,A)" ) "# NORMAL MODES time / fs | E1, E2, E3, E4 / Eh "
            WRITE(DebugUnitNormalModes,800) 0.0, SystemNormalModeEnergies( X(1:4), V(1:4) )

          ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kStep = 0

         PRINT "(/,A)", " Propagating the H-Graphene system in time... "
         
         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = VibrRelaxPotential( X, A )
         DO iCoord = 1,NDim
            A(iCoord) = A(iCoord) / MassVector(iCoord)
         END DO

         WRITE(883,*)

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, VibrRelaxPotential, PotEnergy, RandomNr )

            ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
            IF ( BathType == SLAB_POTENTIAL ) THEN 
               X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
            ENDIF

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 

               ! Energy of the system
               CALL IstantaneousEnergies( PotEnergy, KSys, VSys, Ecoup, VBath, KBath, iTraj, iStep )

               ! Update averages over time
               AverageESys(kStep)            = AverageESys(kStep)  + KSys + VSys
               AverageEBath(kStep)           = AverageEBath(kStep) + KBath + VBath
               AverageECoup(kStep)           = AverageECoup(kStep) + Ecoup

               ! store the autocorrelation function for power spectrum computation
               IF ( PrintType >= FULL )  THEN
                  AverageCoord(1:NDim,kStep) = AverageCoord(1:NDim,kStep) + X(:)
               END IF
               IF ( BathType == CHAIN_BATH .AND. PrintType >= FULL ) THEN
                  DO iCoord = 1, NDim-4
                     OscillCorrelations(iCoord, kStep) = OscillCorrelations(iCoord, kStep) + X(4) * X(4+iCoord) 
                  END DO
               END IF

               ! If massive level of output, print information on normal modes
               IF ( PrintType == DEBUG ) THEN
                  IF ( kStep == 1 ) WRITE(NormalModesUnit,*) " "
                  NormalModes(1:4) = X(1:4)
                  NormalModes(3) = NormalModes(3) - HZEquilibrium
                  NormalModes(4) = NormalModes(4) - C1Puckering
                  NormalModes(:) = NormalModes(:) * SQRT(MassVector(1:4))
                  NormalModes = TheOneWithMatrixVectorProduct( NormalModes4D_Vecs, NormalModes )
                  WRITE(NormalModesUnit,"(10F15.6)") TimeStep*real(iStep)/MyConsts_fs2AU, NormalModes(1:4)
               END IF 

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                        KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                          (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
                  WRITE(DebugUnitNormalModes,800) TimeStep*real(iStep)/MyConsts_fs2AU,  &
                                             SystemNormalModeEnergies( X(1:4), V(1:4) )
               END IF

            END IF 

         END DO

         PRINT "(A)", " Time propagation completed! "

         ! Print log information about the final condition of the trajectory
         WRITE(*,601)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature

         IF ( PrintType == DEBUG ) THEN
               CLOSE( Unit=DebugUnitEn )
               CLOSE( Unit=DebugUnitCoord )
               CLOSE( Unit=DebugUnitVel )
               CLOSE( Unit=DebugUnitNormalModes )
         ENDIF

      END DO
      
      PRINT "(A)"," Done! "

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize averages 
      AverageESys(:)            = AverageESys(:)    / real(NrTrajs)
      AverageEBath(:)           = AverageEBath(:)   / real(NrTrajs)
      AverageECoup(:)           = AverageECoup(:)   / real(NrTrajs)
      IF ( PrintType >= FULL )    AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH )    &
               OscillCorrelations(1:NDim-4,:) = OscillCorrelations(1:NDim-4,:)  / real(NrTrajs)

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(AvEnergyOutputUnit,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
          AverageESys(iStep)*MyConsts_Hartree2eV, AverageECoup(iStep)*MyConsts_Hartree2eV,  AverageEBath(iStep)*MyConsts_Hartree2eV
      END DO

      IF ( PrintType >= FULL ) THEN

         ! PRINT average coordinates
         DO iStep = 0, NrOfPrintSteps
            WRITE(AvCoordOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                       AverageCoord(1:4,iStep)*MyConsts_Bohr2Ang 
         END DO

         ! PRINT average bath coordinates
         DO iCoord = 1, NDim-4
            WRITE(AvBathCoordUnit,"(/,A,I5,/)") "#  Bath Coord # ", iCoord
            DO iStep = 0, NrOfPrintSteps
               WRITE(AvBathCoordUnit,"(F20.8,3F20.8)") real(iCoord)+AverageCoord(4+iCoord,iStep)*10, &
                        TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU
            END DO
         END DO

         IF ( BathType == CHAIN_BATH ) THEN
            ALLOCATE( AverCoordOverTime( NDim ) )
            AverCoordOverTime(:) = SUM ( AverageCoord, 2 ) / real( NrOfPrintSteps+1 )
            DO iCoord = 1, NDim-4
               OscillCorrelations(iCoord,:) = OscillCorrelations(iCoord,:) - &
                                    AverCoordOverTime(3+1)*AverCoordOverTime(3+iCoord+1) 
               WRITE(OscillCorrUnit, *) iCoord, SUM( OscillCorrelations(iCoord, :) )/real(NrOfPrintSteps+1)
            END DO
            DEALLOCATE( AverCoordOverTime )
         END IF

      ENDIF

      ! Close output files
      CLOSE( AvEnergyOutputUnit )
      IF ( PrintType >= FULL ) THEN
         CLOSE( AvCoordOutputUnit )
         CLOSE( AvBathCoordUnit )
      END IF
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) THEN
         CLOSE( OscillCorrUnit ) 
      END IF
      IF ( PrintType == DEBUG ) THEN
         CLOSE( InitialCondUnit )
      END IF

      800 FORMAT(F12.5,1000F15.8)
      600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                     " * Energy of the system (eV)        ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)  ",1F10.4,/    &
                     " * Istantaneous temperature (K)     ",1F10.4,/ ) 
      601 FORMAT (/, " Final condition of the MD trajectory ",/   &
                     " * Energy of the system (eV)        ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)  ",1F10.4,/    &
                     " * Istantaneous temperature (K)     ",1F10.4,/ ) 
      700 FORMAT (/, " Info on initially excited normal mode ",             /, &
                     " * Frequency of the normal mode        ",1F15.6,1X,A, /, &
                     " * 1st component of mass-scaled coords ",1F15.6,1X,A, /, &
                     " * 2nd component of mass-scaled coords ",1F15.6,1X,A, /, &
                     " * 3rd component of mass-scaled coords ",1F15.6,1X,A, /, &
                     " * 4th component of mass-scaled coords ",1F15.6,1X,A, /   ) 

   END SUBROUTINE VibrationalRelax_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( X, V, A, MassVector )
      DEALLOCATE( AverageESys, AverageEBath, AverageECoup )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord  )
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) DEALLOCATE( OscillCorrelations )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( InitialConditions )

   END SUBROUTINE VibrationalRelax_Dispose

!*************************************************************************************************

   REAL FUNCTION VibrRelaxPotential( Positions, Forces )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i       

      ! Check the number degrees of freedom
      NrDOF = size( Positions )
      CALL ERROR( size(Forces) /= NrDOF, "VibrationalRelax.VibrRelaxPotential: array dimension mismatch" )
      CALL ERROR( NrDOF /= NDim, "VibrationalRelax.VibrRelaxPotential: wrong number of DoFs" )

      ! Initialize forces and potential
      VibrRelaxPotential = 0.0
      Forces(:)          = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         VibrRelaxPotential = VHSticking( Positions, Forces )

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         IF ( ModelHarmonic ) THEN
            VibrRelaxPotential = ModelFourDimensionalPotential( Positions(1:4), Forces(1:4) )
         ELSE 
            VibrRelaxPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         END IF
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, Positions(4)-C1Puckering, Positions(5:), VibrRelaxPotential, &
                                                                               Forces(4), Forces(5:) ) 
         CALL BathPotentialAndForces( Bath, Positions(4)-C1Puckering, Positions(5:), VibrRelaxPotential, &
                                                                               Forces(4), Forces(5:) ) 
      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Compute potential and forces of the system
         IF ( ModelHarmonic ) THEN
            VibrRelaxPotential = ModelFourDimensionalPotential( Positions(1:4), Forces(1:4) )
         ELSE 
            VibrRelaxPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         END IF
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), VibrRelaxPotential, &
                                                                               Forces(4), Forces(5:NBath+4) ) 
         CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), VibrRelaxPotential, &
                                                                               Forces(4), Forces(NBath+5:2*NBath+4) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential and forces of the system
         IF ( ModelHarmonic ) THEN
            VibrRelaxPotential = ModelFourDimensionalPotential( Positions(1:4), Forces(1:4) )
         ELSE 
            VibrRelaxPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         END IF

      END IF

   END FUNCTION VibrRelaxPotential

!*************************************************************************************************

   SUBROUTINE IstantaneousEnergies( PotEnergy, KSys, VSys, Ecoup, VBath, KBath, iTraj, iStep )
      IMPLICIT NONE
      REAL, INTENT(IN)  :: PotEnergy
      REAL, INTENT(OUT) :: KSys, VSys, Ecoup, VBath, KBath
      REAL, DIMENSION(4) :: Dummy, TraslCoord
      INTEGER :: iTraj, iStep
      REAL :: KinEnergy, VBath1, VBath2, ECoup1, ECoup2

      ! Total kinetic energy
      KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )

      ! Energy of the system
      KSys = EOM_KineticEnergy( InitialConditions, V(1:4) )
      IF ( BathType == SLAB_POTENTIAL ) THEN 
         TraslCoord(1:2) = X(1:2)
         TraslCoord(3:4) = X(3:4) - (X(5)+X(6)+X(7))/3.0
         VSys = VHFourDimensional( TraslCoord(1:4), Dummy )
      ELSE
         VSys = VHFourDimensional( X(1:4), Dummy )
      END IF

      ! Coupling energy and energy of the bath
      IF ( BathType == SLAB_POTENTIAL ) THEN 
         Ecoup = 0.0
         VBath = PotEnergy - VSys
         KBath = KinEnergy - KSys
      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         CALL EnergyOfTheBath( Bath, X(4)-C1Puckering, X(5:), Ecoup, VBath )
         KBath = KinEnergy - KSys
      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         CALL EnergyOfTheBath( DblBath(1), X(4)-C1Puckering, X(5:NBath+4), Ecoup1, VBath1 )
         CALL EnergyOfTheBath( DblBath(2), X(4)-C1Puckering, X(NBath+5:2*NBath+4), Ecoup2, VBath2 )
!          WRITE(300+iTraj,*) TimeStep*real(iStep)/MyConsts_fs2AU, VBath1, VBath2
         Ecoup = Ecoup1 + Ecoup2
         VBath = VBath1 + VBath2
         KBath = KinEnergy - KSys
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         Ecoup = 0.0
         VBath = 0.0
         KBath = 0.0
      END IF

   END SUBROUTINE IstantaneousEnergies

!*************************************************************************************************

   FUNCTION SystemNormalModeEnergies( Pos, Vel ) RESULT( Energy )
      IMPLICIT NONE
      REAL, DIMENSION(4) :: Energy
      REAL, DIMENSION(4), INTENT(IN) :: Pos, Vel
      REAL, DIMENSION(4) :: Q, V

      ! Remove average value of the coordinate
      Q(:) = Pos(:) - (/ 0.0, 0.0, HZEquilibrium, C1Puckering /)
      V(:) = Vel(:) 

      ! Multiply by square root of mass
      Q(:) = SQRT( MassVector(1:4) ) * Q(:)
      V(:) = SQRT( MassVector(1:4) ) * V(:)

      ! Compute normal modes velocity and coordinates
      Q = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix( NormalModes4D_Vecs ), Q )
      V = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix( NormalModes4D_Vecs ), V )
!       Q = TheOneWithMatrixVectorProduct(  NormalModes4D_Vecs , Q )
!       V = TheOneWithMatrixVectorProduct(  NormalModes4D_Vecs , V )

      ! Compute potential energy associated to the normal modes
      Energy(:) = 0.5 * NormalModes4D_Freq(:) * Q(:)**2

      ! Compute potential energy associated to the normal modes
      Energy(:) = Energy(:) + 0.5 * V(:)**2

   END FUNCTION SystemNormalModeEnergies

!*************************************************************************************************

   REAL FUNCTION ModelFourDimensionalPotential( Coordinates, Forces ) RESULT(V)
      IMPLICIT NONE
      REAL, DIMENSION(4), INTENT(IN)  :: Coordinates
      REAL, DIMENSION(4), INTENT(OUT) :: Forces
      REAL, DIMENSION(4)  :: MassScaledX

      ! Shift coordinates with respect to the minimum
      MassScaledX(1:2) = Coordinates(1:2)
      MassScaledX(3)   = Coordinates(3) - HZEquilibrium
      MassScaledX(4)   = Coordinates(4) - C1Puckering
      ! Scale coordinates with masses
      MassScaledX(:) = MassScaledX(:) * SQRT( MassVector(1:4) )
      ! Compute potential and forces
      V = 0.5 * TheOneWithVectorDotVector( MassScaledX , TheOneWithMatrixVectorProduct( HessianMatrix, MassScaledX ) )
      Forces(:) = - SQRT( MassVector(1:4) ) * TheOneWithMatrixVectorProduct( HessianMatrix, MassScaledX )

   END FUNCTION ModelFourDimensionalPotential


!*************************************************************************************************

END MODULE VibrationalRelax
