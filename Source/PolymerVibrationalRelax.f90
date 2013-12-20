!***************************************************************************************
!*                              MODULE VibrationalRelax
!***************************************************************************************
!
!>  \brief     Subroutines for the simulations of vibrational relaxation
!>             using ring polymer dynamics
!>  \details   This module contains subroutines to compute averages  \n
!>             in a vibrational relaxation simulation                \n
!>             with N replicas of the system connected via harmonic springs
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             28 October 2013
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 
!
!>  \todo          * include constrained ring polymer relaxation for initial system coords
!>                 * clean the code and use unit conversion
!>                 * check OMP variables declaration
!>                 
!***************************************************************************************

MODULE PolymerVibrationalRelax
#include "preprocessoptions.cpp"
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator
   USE FFTWrapper

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: PolymerVibrationalRelax_ReadInput, PolymerVibrationalRelax_Initialize, &
             PolymerVibrationalRelax_Run, PolymerVibrationalRelax_Dispose

   REAL, PARAMETER    :: TimeBetweenSnaps  = 10*MyConsts_fs2AU      !< Nr of initial snapshots to randomize initial conditions
   INTEGER, PARAMETER :: NrOfInitSnapshots = 100                   !< Time between each initial snapshot

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim
   !> Nr of dimension of the system
   INTEGER :: NSystem

   !> Morse potential for the system
   LOGICAL :: MorsePotential
   REAL    :: MorseAlpha
   REAL    :: MorseDe
   REAL    :: HarmonicFreq

   !> variables of the ring polymer dynamics 
   INTEGER :: NBeads                    !< Nr of replicas of the system + bath 
   REAL    :: BeadsForceConst           !< Force constant of the harmonic spring between the beads
   REAL    :: BeadsFrequency            !< Harmonic frequency of the spring between the beads

   ! Initial conditions of the system
   REAL    :: InitEnergy                !< Initial energy of the system

   ! Initial conditions of the bath
   REAL    :: Temperature             !< Temperature of the simulation
   REAL    :: InitTemp

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of time steps between each printing of propagation info ( = NrOfSteps / NrOfPrintSteps )

   ! Time evolution dataset
   TYPE(Evolution), DIMENSION(:), POINTER :: MolecularDynamics   !< Propagate in micro/canonical ensamble to extract results
   TYPE(Evolution), DIMENSION(:), POINTER :: InitialConditions   !< Propagate in microcanonical ensamble to generate initial conditions
   TYPE(Evolution), SAVE :: MicrocanonicalInit                   !< Propagate in microcanonical ensamble to generate initial conditions

   ! Transform to ring normal modes 
   TYPE(FFTHalfComplexType), DIMENSION(:), POINTER  :: RingNormalModes   !< Transform ring coordinates to normal modes (and viceversa)

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: PositionCorrelation       !< Position correlation function
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord              !< Average i-th coordinate vs time
   INTEGER, PARAMETER                   :: NrOfEnergyAverages = 5    !< Nr of average values computed by the function EnergyAverages
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageE                  !< Traj averages of the energy in time

   ! Internal state of the random number generator
   TYPE(RNGInternalState), SAVE :: RandomNr                          !< Internal state of the random number generator
   
   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific to the vibrational
!> relaxation problem.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! Nr of replicas of the system for the ring polymer dynamics
      CALL SetFieldFromInput( InputData, "NBeads", NBeads )

      ! Morse toy model
      CALL SetFieldFromInput( InputData, "MorsePotential", MorsePotential, .FALSE. )

      IF ( MorsePotential ) THEN                ! Morse Parameters 
         CALL SetFieldFromInput( InputData, "MorseAlpha", MorseAlpha )
         MorseAlpha = MorseAlpha / LengthConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "MorseDe", MorseDe )
         MorseDe = MorseDe * EnergyConversion(InputUnits, InternalUnits)
         HarmonicFreq = SQRT( 2.0 * MorseAlpha**2 * MorseDe / MassH )
         NSystem = 1
      ELSE
         NSystem = 4
      END IF

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Initial energy of the system (input in eV and transform to AU)
      CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
      InitEnergy = InitEnergy * EnergyConversion(InputUnits, InternalUnits)

      ! READ THE VARIABLES TO SET THE EQUILIBRATION CONDITIONS

      ! Temperature of the bath
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 1.0 )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)

      ! Initial temperature of the system
      CALL SetFieldFromInput( InputData, "InitTemp", InitTemp )
      InitTemp = InitTemp * TemperatureConversion(InputUnits, InternalUnits)

      ! Frequency and force constant of the harmonic force between the beads
      BeadsFrequency = NBeads * Temperature
      BeadsForceConst = ( NBeads * Temperature )**2

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

      WRITE(*, 902) NBeads, BeadsFrequency/MyConsts_cmmin1toAU

      IF ( .NOT. MorsePotential ) THEN
           WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV
      ELSE 
           WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy+MorseDe)*MyConsts_Hartree2eV
      END IF

      WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps

      WRITE(*, 905) Temperature/MyConsts_K2AU


   902 FORMAT(" * Nr of replicas in the ring polymer dynamics: ", I10,/, &
              " * Frequency of the force between the beads     ", F10.4,/ )

   903 FORMAT(" * Initial conditions of the atom-surface system ", /,&
              " * Absolute initial energy (eV):                ",F10.4,/,& 
              " *  - w.r.t. the bottom of the well (eV):       ",F10.4,/ )

   904 FORMAT(" * Dynamical simulation variables               ", /,&
              " * Nr of trajectories:                          ",I10,  /,& 
              " * Propagation time step (fs):                  ",F10.4,/,& 
              " * Nr of time steps of each trajectory:         ",I10,  /,& 
              " * Nr of print steps of each trajectory:        ",I10,  / )

   905 FORMAT(" * Bath equilibration variables                 ", /,&
              " * Temperature of the surface:                  ",F10.4 )

   END SUBROUTINE PolymerVibrationalRelax_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord, iBead, NrOfThreads, iThread

      CALL ERROR( BathType ==  SLAB_POTENTIAL, " PolymerVibrationalRelax_Initialize: slab potential not implemented" )
      CALL ERROR( BathType ==  DOUBLE_CHAIN, " PolymerVibrationalRelax_Initialize: double chain not implemented" )
      CALL ERROR( BathType ==  LANGEVIN_DYN, " PolymerVibrationalRelax_Initialize: langevin dynamics not implemented" )

      ! Define dimension of the system+bath and mass vector
      NDim = NSystem + NBath
      ALLOCATE( MassVector(NDim) )
      IF ( MorsePotential ) THEN
         MassVector = (/ MassH, (MassBath, iCoord=1,NBath) /)
      ELSE 
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,NBath) /)
      END IF
      
      ! Allocate memory and initialize vectors for trajectory, acceleration and masses
      ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads) )

      !$OMP PARALLEL 
      ! Set total nr of OPENMP threads
      !$OMP MASTER
      NrOfThreads = __OMP_TotalNrOfThreads
      !$OMP END MASTER
      !$OMP END PARALLEL 
      
      ! Propagation data is replicated, so that in the parallel region of the code
      ! each threads has its own set of variables
      ALLOCATE( MolecularDynamics(NrOfThreads), RingNormalModes(NrOfThreads),  InitialConditions(NrOfThreads) )
      DO iThread = 1, NrOfThreads
         ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
         CALL EvolutionSetup( MolecularDynamics(iThread), NDim, MassVector, TimeStep )
         ! Set ring polymer molecular dynamics parameter
         CALL SetupRingPolymer( MolecularDynamics(iThread), NBeads, BeadsFrequency ) 
         
         ! Set transform from ring coordinates to normal modes
         CALL SetupFFT( RingNormalModes(iThread), NBeads ) 

         ! Set variables for EOM integration of the system only in the microcanonical ensamble 
         ! this will be done in a serial way, so no replicated data
         CALL EvolutionSetup( InitialConditions(iThread), NSystem, MassVector(1:NSystem), TimeStep )
         ! Set ring polymer molecular dynamics parameter
         CALL SetupRingPolymer( InitialConditions(iThread), NBeads, BeadsFrequency ) 
      END DO

      CALL EvolutionSetup( MicrocanonicalInit, NSystem, MassVector(1:NSystem), TimeStep )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the average energies over time only
      ! DEBUG - print also the average coordinates
      ! FULL - print detailed information about the initial conditions and about the trajectories

      ! Allocate and initialize the variable for the trajectory average of the position correlation
      ALLOCATE( PositionCorrelation(0:NrOfPrintSteps) )
      PositionCorrelation(:)  = 0.0
      
      ! Allocate and initialize the variable for the trajectory averages of the energy
      ALLOCATE( AverageE(NrOfEnergyAverages,0:NrOfPrintSteps) ) 
      AverageE(:,:) = 0.0
      
      ! Average coordinates over time 
      IF ( PrintType >= FULL ) THEN
         ALLOCATE( AverageCoord(NDim,0:NrOfPrintSteps) )
         AverageCoord(1:NDim,0:NrOfPrintSteps) = 0.0
      END IF

   END SUBROUTINE PolymerVibrationalRelax_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Run()
      IMPLICIT NONE
      INTEGER :: CurrentThread, NrOfThreads
      !> Output units: minimal and standard output
      INTEGER  ::  AvEnergyOutputUnit, AvCoordOutputUnit, AvBathCoordUnit, TLogUnit, PosCorrelationUnit
      !> Output units: debug output
      INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel
      !> Initial energies of the bath
      REAL     ::  InitKinAverage, InitKinVariance, InitPotAverage, InitPotVariance
      !> Energy averages
      REAL     ::  TotEnergy, PotEnergy, KinEnergy
      !> integer Counters
      INTEGER  ::  iTraj, iBead, iStep, kStep, iCoord
      !> Filename for output files
      CHARACTER(100) :: OutFileName
      !> Centroid of positions and velocities
      REAL, DIMENSION(NDim*NBeads) :: X0
      REAL, DIMENSION(NSystem*NBeads) :: InitX, InitV
      REAL, DIMENSION(NrOfInitSnapshots,NSystem*2) :: SystemInitConditions

      REAL, DIMENSION(NBath,NBeads) :: InitQBath, InitVBath

      PosCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="PositionCorrelation.dat", UNIT=PosCorrelationUnit )
      WRITE(PosCorrelationUnit, "(3A,I6,A,/)") "# <X(0)X(t)> vs time (",TimeUnit(InputUnits)," | au) - ", NrTrajs, " trajectories "
                                                                                  
      AvEnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="EnergyAverages.dat", UNIT=AvEnergyOutputUnit )
      WRITE(AvEnergyOutputUnit, "(5A,I6,A,/)") "# E vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(EnergyUnit(InputUnits)), ") - ", NrTrajs, " trajectories "
      
      IF ( PrintType >= FULL ) THEN
         AvCoordOutputUnit = LookForFreeUnit()
         OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
         WRITE(AvCoordOutputUnit, "(5A,I6,A,/)") "# System coordinate vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(LengthUnit(InputUnits)), ") - ", NrTrajs, " trajectories "

         AvBathCoordUnit = LookForFreeUnit()
         OPEN( FILE="AverageBathCoords.dat", UNIT=AvBathCoordUnit )
         WRITE(AvBathCoordUnit, "(5A,I6,A,/)") "# Bath coordinate vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(LengthUnit(InputUnits)), ") - ", NrTrajs, " trajectories "
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      ! Define a set of system initial conditions with microcanonical propagation of the system
      CALL MicrocanonicalSamplingOfTheSystem( SystemInitConditions )

      ! Initialize averages of the bath initial state
      InitKinAverage  = 0.0
      InitKinVariance = 0.0
      InitPotAverage  = 0.0
      InitPotVariance = 0.0

      !$OMP PARALLEL PRIVATE( CurrentThread, RandomNr, iBead, iStep, kStep, InitQBath, InitVBath,            &
      !$OMP&                  X, X0, V, A, PotEnergy, KinEnergy,  InitX, InitV,                   & 
      !$OMP&                        OutFileName, DebugUnitEn, DebugUnitCoord, DebugUnitVel )

      ! Set total nr of OPENMP threads and nr of the current thread
      !$OMP MASTER
      NrOfThreads = __OMP_TotalNrOfThreads
      !$OMP END MASTER
      CurrentThread = __OMP_CurrentThreadNum

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1-CurrentThread+1 )
      
      !run NrTrajs number of trajectories
      !$OMP DO REDUCTION( + : AverageE, PositionCorrelation, AverageCoord,    &
      !$OMP&                  InitKinAverage, InitKinVariance, InitPotAverage, InitPotVariance )
      DO iTraj = 1,NrTrajs
      
         __OMP_OnlyMasterBEGIN PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" __OMP_OnlyMasterEND

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! *************  Initial conditions of the system *****************

         ! Choose a random integer and choose a x,p in microcanonical distribution

         kStep = CEILING( UniformRandomNr(RandomNr)*real(NrOfInitSnapshots) )
         InitX = SystemInitConditions( kStep, 1 )
         InitV = SystemInitConditions( kStep, 2 )

         IF ( InitTemp > 0.0 ) THEN
             CALL RingPolymerInitialEquilibration(  InitX, InitV, CurrentThread )
         END IF

         DO iBead = 1, NBeads
            X((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = InitX((iBead-1)*NSystem+1: iBead*NSystem)
            V((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = InitV((iBead-1)*NSystem+1: iBead*NSystem)
         END DO
! 
!          DO iBead = 1, NBeads
!             X((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = SystemInitConditions( kStep, 1:NSystem )
!             V((iBead-1)*NDim+1:(iBead-1)*NDim+NSystem) = SystemInitConditions( kStep, NSystem+1:NSystem*2 )
!          END DO

         ! *************  Initial conditions of the bath *****************

         ! Use harmonic oscillator sampling for the normal modes of the ring polymer bath
         CALL BathOfRingsThermalConditions( Bath, NBeads, BeadsFrequency, Temperature*NBeads, InitQBath, InitVBath, & 
                     PotEnergy, KinEnergy, RandomNr, RingNormalModes(CurrentThread) )

         DO iBead = 1, NBeads
            X( (iBead-1)*NDim+NSystem+1 : iBead*NDim ) = InitQBath(:,iBead)
            V( (iBead-1)*NDim+NSystem+1 : iBead*NDim ) = InitVBath(:,iBead)
         END DO

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************
         
         ! Increment averages of the initial conditions of the bath
         InitKinAverage  = InitKinAverage + KinEnergy/real(NBath*NBeads) 
         InitKinVariance = InitKinVariance + (KinEnergy/real(NBath*NBeads))**2
         InitPotAverage  = InitPotAverage + PotEnergy/real(NBath*NBeads)
         InitPotVariance = InitPotVariance + (PotEnergy/real(NBath*NBeads))**2
         
         ! PRINT INFORMATION ON INITIAL CONDITIONS of THE BATH and THE SYSTEM
!         __OMP_OnlyMasterBEGIN
!          WRITE(*,600)  (KSys+VSys)*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, 2.0*KinEnergy/MyConsts_K2AU, &
!                                                         PotEnergy*MyConsts_Hartree2eV, 2.0*PotEnergy/MyConsts_K2AU
!         __OMP_OnlyMasterEND

         ! Store initial coordinates
         X0 = X
         ! Compute position correlation function
         PositionCorrelation(0) = PositionCorrelation(0) + CorrelationFunction( X, X0 )
         
         ! Various expectation values of the energy
         AverageE(:,0) = AverageE(:,0) + EnergyAverages( X, V, MassVector ) 

         ! Average values of the centroid coordinates
         IF ( PrintType >= FULL ) AverageCoord(1:NDim,0) = AverageCoord(1:NDim,0) + CentroidCoord( X )

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

            ! Write initial values
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | Ekin, Epot, Etot / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KinEnergy, PotEnergy, KinEnergy+PotEnergy

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, X(:)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)

          END IF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kStep = 0

         __OMP_OnlyMasterBEGIN 
         PRINT "(/,A)", " Propagating system and bath in the microcanonical ensamble... " __OMP_OnlyMasterEND
         
         ! Compute starting potential and forces
         CALL EOM_RPMSymplectic( MolecularDynamics(CurrentThread), X, V, A,  VibrRelaxPotential, PotEnergy, RandomNr, 1 )

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_RPMSymplectic( MolecularDynamics(CurrentThread), X, V, A, VibrRelaxPotential, PotEnergy, RandomNr )

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 

               ! Various expectation values of the energy
               AverageE(:,kStep) = AverageE(:,kStep) + EnergyAverages( X, V, MassVector ) 

               ! Compute position correlation function
               PositionCorrelation(kStep) = PositionCorrelation(kStep) + CorrelationFunction( X, X0 )
         
               ! Average values of the centroid coordinates
               IF ( PrintType >= FULL )    AverageCoord(1:NDim,kStep) = AverageCoord(1:NDim,kStep) + CentroidCoord( X )
               
               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, KinEnergy, PotEnergy, KinEnergy+PotEnergy
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(:)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
               END IF

            END IF 

         END DO

         __OMP_OnlyMasterBEGIN PRINT "(A)", " Time propagation completed! " __OMP_OnlyMasterEND

         IF ( PrintType == DEBUG ) THEN
            CLOSE( Unit=DebugUnitEn )
            CLOSE( Unit=DebugUnitCoord )
            CLOSE( Unit=DebugUnitVel )
         END IF

!          __OMP_OnlyMasterBEGIN
!          ! Print log information about the final condition of the trajectory
!          WRITE(*,601)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature
!          __OMP_OnlyMasterEND

      END DO
      !$OMP END DO
      !$OMP END PARALLEL
      
      PRINT "(A)"," Done! "

      InitKinAverage = InitKinAverage / real(NrTrajs) 
      InitKinVariance = SQRT( InitKinVariance/real(NrTrajs) - InitKinAverage**2 )
      InitPotAverage = InitPotAverage / real(NrTrajs) 
      InitPotVariance = SQRT( InitPotVariance/real(NrTrajs) - InitPotAverage**2 ) 
      
      TLogUnit = LookForFreeUnit()
      OPEN( UNIT = TLogUnit, FILE = "FinalAverTemperature.log" )
      WRITE(TLogUnit,500)  NrTrajs, 2.*InitKinAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits),  &
                                    2.*InitKinVariance*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
                                    2.*InitPotAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits),  &
                                    2.*InitPotVariance*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)
      CLOSE( UNIT = TLogUnit )

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize averages 
      PositionCorrelation(:)  =  PositionCorrelation(:)  / real(NrTrajs)
      AverageE(:,:)           = AverageE(:,:)            / real(NrTrajs)
      IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(AvEnergyOutputUnit,"(F14.8,100F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                                     AverageE(:,iStep)*MyConsts_Hartree2eV
      END DO

      DO iStep = 0, NrOfPrintSteps
         WRITE(PosCorrelationUnit,"(F14.8,100F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                                     PositionCorrelation(iStep)*MyConsts_Hartree2eV
      END DO

      IF ( PrintType >= FULL ) THEN

         ! PRINT average coordinates
         DO iStep = 0, NrOfPrintSteps
            WRITE(AvCoordOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                       AverageCoord(1:NSystem,iStep)*MyConsts_Bohr2Ang 
         END DO

         ! PRINT average bath coordinates
         DO iCoord = 1, NBath
            WRITE(AvBathCoordUnit,"(/,A,I5,/)") "#  Bath Coord # ", iCoord
            DO iStep = 0, NrOfPrintSteps
               WRITE(AvBathCoordUnit,"(F20.8,3F20.8)") real(iCoord)+AverageCoord(NSystem+iCoord,iStep)*10, &
                        TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU
            END DO
         END DO

      END IF

      ! Close output files
      CLOSE(PosCorrelationUnit)
      CLOSE( AvEnergyOutputUnit )
      IF ( PrintType >= FULL ) THEN
         CLOSE( AvCoordOutputUnit )
         CLOSE( AvBathCoordUnit )
      END IF

      800 FORMAT(F12.5,1000F15.8)
      600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                     " * Energy of the system (eV)          ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)    ",1F10.4,/    &
                     " * Bath translational temperature (K) ",1F10.4,/    &
                     " * Potential Energy of the bath (eV)  ",1F10.4,/    &
                     " * Bath vibrational temperature (K)   ",1F10.4,/ ) 
      601 FORMAT (/, " Final condition of the MD trajectory ",/   &
                     " * Energy of the system (eV)         ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)   ",1F10.4,/    &
                     " * Istantaneous temperature (K)      ",1F10.4,/ ) 
      500 FORMAT (/, " Nr of trajectories :                               ",1I10,       2/, &
                     " Average Translational Temperature :                ",1F15.6,1X,A, /, &
                     " Standard Deviation of Translational Temperature :  ",1F15.6,1X,A,2/, &
                     " Average Vibrational Temperature :                  ",1F15.6,1X,A, /, &
                     " Standard Deviation of Vibrational Temperature :    ",1F15.6,1X,A, / ) 

   END SUBROUTINE PolymerVibrationalRelax_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Dispose()
      IMPLICIT NONE
      INTEGER :: i 

      ! Deallocate memory 
      DEALLOCATE( X, V, A, MassVector )
      DEALLOCATE( AverageE )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord )

      ! Unset propagators 
      DO i = 1, size(MolecularDynamics)
         CALL DisposeEvolutionData( MolecularDynamics(i) )
         CALL DisposeFFT( RingNormalModes(i) )
         CALL DisposeEvolutionData( InitialConditions(i) )
      END DO

      DEALLOCATE( MolecularDynamics, RingNormalModes, InitialConditions )

   END SUBROUTINE PolymerVibrationalRelax_Dispose

!*************************************************************************************************

   REAL FUNCTION VibrRelaxPotential( Positions, Forces ) RESULT(V)
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      REAL, DIMENSION(NBath)   :: ForceQCoord
      REAL, DIMENSION(NSystem) :: ForceSystem
      REAL                     :: ForceCoupCoord, CouplingFunc

      ! Check the number degrees of freedom
      CALL ERROR( size(Forces) /= size( Positions ), "PolymerVibrationalRelax.VibrRelaxPotential: array dimension mismatch" )
      CALL ERROR( size( Positions ) /= NDim,         "PolymerVibrationalRelax.VibrRelaxPotential: wrong number of DoFs"     )
     
      V = 0.0
      Forces(:) = 0.0

      ! Compute potential and forces of the system replica
      IF ( MorsePotential ) THEN
         V = V + MorseV( Positions(1:1), ForceSystem )
      ELSE
         V = V + VHFourDimensional( Positions(1:4), ForceSystem )
      END IF
      Forces( 1:NSystem ) = Forces( 1:NSystem ) + ForceSystem(:)

      ! Compute potential and forces of the bath replica
      ForceCoupCoord = 0.0
      ForceQCoord(:) = 0.0
      IF ( MorsePotential ) THEN
         CouplingFunc = (exp( - MorseAlpha*Positions(1) ) - 1.0) / MorseAlpha
      ELSE 
         CouplingFunc = Positions(4) - C1Puckering
      END IF
      CALL BathPotentialAndForces( Bath, CouplingFunc, Positions( NSystem+1:NDim ), V, ForceCoupCoord, ForceQCoord(:) ) 
      IF ( MorsePotential ) THEN
         Forces(1) = Forces(1) + ForceCoupCoord * ( - exp( - MorseAlpha*Positions(1) ) )
      ELSE
         Forces(4) = Forces(4) + ForceCoupCoord
      END IF
      Forces( NSystem+1:NDim ) = Forces( NSystem+1:NDim ) + ForceQCoord(:)

      ! when collinear calculation, make sure that the forces on HX and HY are zero
      IF ( Collinear .AND. (.NOT. MorsePotential ))   Forces(1:2) = 0.0
      
   END FUNCTION VibrRelaxPotential

   REAL FUNCTION MorseV( Positions, Forces ) RESULT(V) 
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 

      V = MorseDe * ( exp(-2.0*MorseAlpha*Positions(1)) - 2.0 * exp(-MorseAlpha*Positions(1)) )  
      Forces(1) = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*Positions(1)) - exp(-MorseAlpha*Positions(1)) )  

   END FUNCTION MorseV
   
   SUBROUTINE PartitionedEnergyAndForces( Positions, Forces, SysV, CouplingV, BathV )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      REAL, INTENT(OUT)                       :: SysV, CouplingV, BathV
      REAL, DIMENSION(NBath)     :: ForceQCoord
      REAL, DIMENSION(NSystem)   :: ForceSystem
      REAL                       :: ForceCoupCoord, CouplingFunc

      ! Check the number degrees of freedom
      CALL ERROR( size(Forces) /= size(Positions), "PolymerVibrationalRelax.PartitionedEnergyAndForces: array dimension mismatch" )
      CALL ERROR( size(Positions) /= NDim,         "PolymerVibrationalRelax.PartitionedEnergyAndForces: wrong number of DoFs"     )
     
      CouplingV= 0.0; BathV = 0.0
      Forces(:) = 0.0

      ! Compute potential and forces of the system replica
      IF ( MorsePotential ) THEN
         SysV = MorseV( Positions(1:1), ForceSystem )
      ELSE
         SysV = VHFourDimensional( Positions(1:4), ForceSystem )
      END IF
      Forces( 1:NSystem ) = Forces( 1:NSystem ) + ForceSystem(:)

      ! Compute potential and forces of the bath replica
      ForceCoupCoord = 0.0
      ForceQCoord(:) = 0.0
      IF ( MorsePotential ) THEN
         CouplingFunc = (exp( - MorseAlpha*Positions(1) ) - 1.0) / MorseAlpha
      ELSE 
         CouplingFunc = Positions(4) - C1Puckering
      END IF
      CALL BathPotentialAndForces( Bath, CouplingFunc, Positions(NSystem+1:NDim), BathV, ForceCoupCoord, ForceQCoord(:), CouplingV ) 
      IF ( MorsePotential ) THEN
         Forces(1) = Forces(1) + ForceCoupCoord * ( - exp( - MorseAlpha*Positions(1) ) )
      ELSE
         Forces(4) = Forces(4) + ForceCoupCoord
      END IF
      Forces( NSystem+1:NDim ) = Forces( NSystem+1:NDim ) + ForceQCoord(:)

      ! when collinear calculation, make sure that the forces on HX and HY are zero
      IF ( Collinear .AND. (.NOT. MorsePotential ))   Forces(1:2) = 0.0
      
   END SUBROUTINE PartitionedEnergyAndForces

!*************************************************************************************************

   FUNCTION CentroidCoord( X )
      IMPLICIT NONE
      REAL, DIMENSION( NDim )   :: CentroidCoord
      REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X
      INTEGER :: i, iCoord

      CentroidCoord( 1:NDim ) = 0.0
      DO i = 1, NBeads
         DO iCoord = 1, NDim
            CentroidCoord( iCoord ) = CentroidCoord( iCoord ) + X( (i-1)*NDim + iCoord )
         END DO
      END DO
      CentroidCoord( 1:NDim ) = CentroidCoord( 1:NDim ) / NBeads

   END FUNCTION CentroidCoord
   
!*************************************************************************************************

   REAL FUNCTION CorrelationFunction( X, X0 )
      IMPLICIT NONE
      REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X, X0
      REAL, DIMENSION( NDim )                    :: XCentroid, X0Centroid

      XCentroid = CentroidCoord( X )
      X0Centroid = CentroidCoord( X0 )
      CorrelationFunction =  TheOneWithVectorDotVector( XCentroid(1:NSystem), X0Centroid(1:NSystem) )

   END FUNCTION CorrelationFunction
   
!*************************************************************************************************
   
   FUNCTION EnergyAverages( X, V, Mass ) 
      IMPLICIT NONE
      REAL, DIMENSION(5) :: EnergyAverages
      REAL, DIMENSION(NDim*NBeads), INTENT(IN) :: X, V
      REAL, DIMENSION(NDim), INTENT(IN)        :: Mass
      REAL                          :: SysV, CouplingV, BathV
      REAL, DIMENSION(NDim*NBeads)  :: SngForce
      REAL, DIMENSION(NDim)         :: CentroidV
      INTEGER                       :: iCoord, iBead

      EnergyAverages = 0.0
      
      ! (2,3,5) POTENTIAL ENERGY    *****************************************************
      ! System/bath potential and coupling averaged over the beads, and force acting on each bead
      SysV = 0.0; CouplingV = 0.0; BathV = 0.0
      DO iBead = 1, NBeads 
         CALL PartitionedEnergyAndForces( X((iBead-1)*NDim+1:iBead*NDim), SngForce((iBead-1)*NDim+1:iBead*NDim), &
                    SysV, CouplingV, BathV )
         EnergyAverages(3) = EnergyAverages(3) + CouplingV
         EnergyAverages(2) = EnergyAverages(2) + SysV 
         EnergyAverages(5) = EnergyAverages(5) + BathV 
      END DO
      EnergyAverages(3) = EnergyAverages(3) / real(NBeads)
      EnergyAverages(2) = EnergyAverages(2) / real(NBeads)
      EnergyAverages(5) = EnergyAverages(5) / real(NBeads)

      ! (1) VIRIAL TOTAL ENERGY    **************************************************
      ! centroid of the virial product (virial as average)
      EnergyAverages(1) = 0.0
      EnergyAverages(4) = 0.0
      DO iBead = 1, NBeads
         DO iCoord = 1, NSystem
            EnergyAverages(1) = EnergyAverages(1) - X((iBead-1)*NDim+iCoord) * SngForce((iBead-1)*NDim+iCoord)
         END DO
         DO iCoord = NSystem + 1, NDim
            EnergyAverages(4) = EnergyAverages(4) - X((iBead-1)*NDim+iCoord) * SngForce((iBead-1)*NDim+iCoord)
         END DO
      END DO
      EnergyAverages(1) = 0.5 * EnergyAverages(1) / real(NBeads)
      EnergyAverages(4) = 0.5 * EnergyAverages(4) / real(NBeads)

! 
!       ! 2) THERMODYNAMICS TOTAL ENERGY
!       ! sum up energy of the ring oscillators and subtract it to average energy
!       EnergyAverages(2) = 0.0 
!       DO iCoord = 1, NDim
!          IF ( NBeads == 2 ) THEN
!             EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(NDim+iCoord) - X(iCoord) )**2 
!          ELSE IF ( NBeads > 2 ) THEN
!             DO iBead = 1, NBeads-1
!                EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iBead*NDim+iCoord) - X((iBead-1)*NDim+iCoord) )**2 
!             END DO
!             EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iCoord) - X((NBeads-1)*NDim+iCoord) )**2 
!          END IF
!       END DO
!       EnergyAverages(2) = 0.5 * ( NBeads * NDim * Temperature - BeadsForceConst * EnergyAverages(2) )
!       ! add potential energy to get full energy
!       EnergyAverages(2) = EnergyAverages(2) + EnergyAverages(4)
! 
!       ! 3) POTENTIAL PLUS KINETIC ENERGY OF THE CENTROID
!       ! Compute centroid of the velocity
!       CentroidV = CentroidCoord( V )
!       EnergyAverages(3) = 0.0
!       DO iCoord = 1, NDim
!          EnergyAverages(3) = EnergyAverages(3) + Mass(iCoord) * CentroidV(iCoord)**2
!       END DO
!       EnergyAverages(3) =  0.5 * EnergyAverages(3)
!       ! add potential energy to get full energy
!       EnergyAverages(3) = EnergyAverages(3) + EnergyAverages(4)

   END FUNCTION EnergyAverages

!*************************************************************************************************

   SUBROUTINE MicrocanonicalSamplingOfTheSystem( InitConditions )
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(OUT) :: InitConditions
      REAL, DIMENSION(NSystem) :: Accel, Pos, Vel
      INTEGER :: NTimeStepEachSnap, iTraj, iStep
      REAL :: PotEnergy
      REAL, DIMENSION(1) :: Dummy

      CALL ERROR( SIZE(InitConditions,2) /= NSystem*2, "MicrocanonicalSamplingOfTheSystem: wrong nr of dimensions"  ) 

      IF ( MorsePotential ) THEN
         Pos(1) = 0.0
         PotEnergy = MorseV( Pos(1:1), Accel(1:1) )
         Vel(1) = sqrt( 2.0 * (InitEnergy-PotEnergy) / MassVector(1) )
         Accel(1:1) = Accel(1:1) / MassVector(1)
      ELSE
         Pos(1:2) = 0.0
         Pos(3) = HZEquilibrium 
         Pos(4) = C1Puckering  
         PotEnergy = VHFourDimensional( Pos(1:4), Accel(1:4) )
         Vel(1) = 0.0
         Vel(2) = 0.0
         Vel(3) = + sqrt( 2.0 * (InitEnergy-PotEnergy) / MassVector(3) ) * 0.958234410548192
         Vel(4) = - sqrt( 2.0 * (InitEnergy-PotEnergy) / MassVector(4) ) * 0.285983940880181 
         Accel(1:4) = Accel(1:4) / MassVector(1:4)
      END IF

      ! Define when to store the dynamics snapshot
      NTimeStepEachSnap = INT( TimeBetweenSnaps / TimeStep )

      ! Cycle over the nr of snapshot to store
      DO iTraj = 1, SIZE(InitConditions,1)

         ! Propagate the 4D traj in the microcanonical ensamble
         DO iStep = 1, NTimeStepEachSnap
            ! Propagate for one timestep with Velocity-Verlet
            IF ( MorsePotential ) THEN
               CALL EOM_LangevinSecondOrder( MicrocanonicalInit, Pos, Vel, Accel, MorseV, PotEnergy, RandomNr )
            ELSE
               CALL EOM_LangevinSecondOrder( MicrocanonicalInit, Pos, Vel, Accel, VHFourDimensional, PotEnergy, RandomNr )
            END IF
         END DO

         ! Store snapshot
         InitConditions( iTraj, 1:NSystem )           = Pos(1:NSystem)
         InitConditions( iTraj, NSystem+1:2*NSystem ) = Vel(1:NSystem)
      END DO

   END SUBROUTINE MicrocanonicalSamplingOfTheSystem

   SUBROUTINE RingPolymerInitialEquilibration( Pos, Vel, CurrentThread )
      IMPLICIT NONE
      REAL, DIMENSION(NSystem*NBeads), INTENT(INOUT) :: Pos, Vel
      INTEGER, INTENT(IN) :: CurrentThread
      REAL, DIMENSION(NSystem*NBeads) :: Accel
      INTEGER :: iStep, iBead, kStep
      REAL :: PotEnergy, KinEnergy, AverKinEnergy, VCentroid, XCentroid
      REAL, DIMENSION(1) :: Dummy

      CALL SetupThermostat( InitialConditions(CurrentThread), 0.01, InitTemp )
      CALL EOM_RPMSymplectic( InitialConditions(CurrentThread), Pos, Vel, Accel, MorseV, PotEnergy, RandomNr, 1 )

      DO iStep = 1, 3000

         ! Propagate for one timestep with Velocity-Verlet
         IF ( MorsePotential ) THEN
            CALL EOM_RPMSymplectic( InitialConditions(CurrentThread), Pos, Vel, Accel, MorseV, PotEnergy, RandomNr, 2 )
         ELSE
            CALL EOM_RPMSymplectic( InitialConditions(CurrentThread), Pos, Vel, Accel, VHFourDimensional, PotEnergy, RandomNr, 2 )
         END IF

      END DO

      CALL DisposeThermostat( InitialConditions(CurrentThread) )

      AverKinEnergy = 0.0
      kStep = 0
      DO iStep = 1, 50000

         IF ( MOD(iStep-1, 500) == 0.0 ) THEN

            kStep = kStep + 1

            XCentroid = 0.0
            PotEnergy = 0.0
            KinEnergy = 0.0
            DO iBead = 1, NBeads
               XCentroid = XCentroid + Pos(iBead)
               PotEnergy = PotEnergy + MorseV( Pos(iBead:iBead), Dummy )
               KinEnergy = KinEnergy - Dummy(1)*Pos(iBead)
            END DO
            XCentroid = XCentroid / NBeads
            PotEnergy = PotEnergy / NBeads
            KinEnergy = 0.5 * KinEnergy / NBeads

!            CALL ExecuteFFT( RingNormalModes(CurrentThread), Vel, DIRECT_FFT )

!            KinEnergy = 0.0
!            DO iBead = 1, NBeads
!               KinEnergy = KinEnergy + Vel(iBead)**2
!            END DO
!            KinEnergy = KinEnergy * 0.5 * MassVector(1) / (NBeads-1)
!            AverKinEnergy = AverKinEnergy + KinEnergy

!            VCentroid = Vel(1)

!            CALL ExecuteFFT( RingNormalModes(CurrentThread), Vel, INVERSE_FFT )
!            WRITE(745,"(I10,10F20.8)") kStep, 2.0*AverKinEnergy/real(kStep)*TemperatureConversion(InternalUnits,InputUnits), &
!                                   XCentroid, VCentroid
            WRITE(746,"(I10,10F20.8)") kStep, (PotEnergy+KinEnergy)*EnergyConversion(InternalUnits,InputUnits), &
                    PotEnergy*EnergyConversion(InternalUnits,InputUnits), KinEnergy*EnergyConversion(InternalUnits,InputUnits)
         END IF

         ! Propagate for one timestep with Velocity-Verlet
         IF ( MorsePotential ) THEN
            CALL EOM_RPMSymplectic( InitialConditions(CurrentThread), Pos, Vel, Accel, MorseV, PotEnergy, RandomNr )
         ELSE
            CALL EOM_RPMSymplectic( InitialConditions(CurrentThread), Pos, Vel, Accel, VHFourDimensional, PotEnergy, RandomNr )
         END IF

      END DO
      WRITE(745,"(/,A,/)") "# --------------------------------- "


   END SUBROUTINE RingPolymerInitialEquilibration

!*************************************************************************************************

END MODULE PolymerVibrationalRelax
