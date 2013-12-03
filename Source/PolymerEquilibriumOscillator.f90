!***************************************************************************************
!*                      MODULE PolymerEquilibriumOscillator
!***************************************************************************************
!
!>  \brief     
!>  \details   
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
!>  \todo          
!>                 
!***************************************************************************************
MODULE PolymerEquilibriumOscillator
#include "preprocessoptions.cpp"
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: PolymerEquilibriumOscillator_ReadInput, PolymerEquilibriumOscillator_Initialize, &
             PolymerEquilibriumOscillator_Run, PolymerEquilibriumOscillator_Dispose

   ! Nr of dimensions of the system
   INTEGER :: NDim

   !> Morse potential for the system
   LOGICAL :: MorsePotential
   REAL    :: MorseAlpha
   REAL    :: MorseDe

   REAL    :: HarmonicFreq

   !> Parameters for the system potential, if it's not a morse
   REAL ::  Quadratic, Cubic, Quartic

   !> System mass
   REAL :: MassSystem

   !> variables of the ring polymer dynamics 
   INTEGER :: NBeads                    !< Nr of replicas of the system + bath 
   REAL    :: BeadsForceConst           !< Force constant of the harmonic spring between the beads
   REAL    :: BeadsFrequency            !< Harmonic frequency of the spring between the beads

   ! Parameters for initial equilibration 
   REAL    :: Temperature               !< Temperature of the simulation
   REAL    :: EquilTStep                !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma                !< Friction parameter of the Langevin equation
   REAL    :: EquilibrationTime         !< Time of the equilibration
   INTEGER :: NrEquilibrationSteps      !< Nr of timesteps for equilibrating the initial conditions

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of time steps between each printing of propagation info ( = NrOfSteps / NrOfPrintSteps )

   ! Time evolution dataset
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/canonical ensamble to extract results
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in canonical ensamble to generate initial conditions of the bath

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: PositionCorrelation    !< Position correlation function
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord           !< Average i-th coordinate vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: KineticAverage         !< Average of the kinetic energy of the system
   REAL, DIMENSION(:), ALLOCATABLE      :: PotentialAverage       !< Average of the potential energy of the system
   REAL, DIMENSION(:), ALLOCATABLE      :: EquilibAveTvsTime      !< Average temperature at time T during equilibration

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific to the equilibrium
!> simulation cproblem.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE PolymerEquilibriumOscillator_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! PARAMETERS FOR THE SYSTEM DEFINITION 

      ! Nr of replicas of the system for the ring polymer dynamics
      CALL SetFieldFromInput( InputData, "NBeads", NBeads )

      ! Mass of the system ( in amu )
      CALL SetFieldFromInput( InputData, "MassSystem", MassSystem )
      MassSystem = MassSystem * MassConversion(InputUnits, InternalUnits)

      ! Potential is a morse function?
      CALL SetFieldFromInput( InputData, "MorsePotential", MorsePotential )

      IF ( MorsePotential ) THEN
         CALL SetFieldFromInput( InputData, "MorseAlpha", MorseAlpha )
         MorseAlpha = MorseAlpha / LengthConversion(InputUnits, InternalUnits)
         CALL SetFieldFromInput( InputData, "MorseDe", MorseDe )
         MorseDe = MorseDe * EnergyConversion(InputUnits, InternalUnits)
         HarmonicFreq = SQRT( 2.0 * MorseAlpha**2 * MorseDe / MassSystem )
      ELSE
         ! Parameters for the system potential ( in a.u. )
         CALL SetFieldFromInput( InputData, "Quadratic", Quadratic )
         CALL SetFieldFromInput( InputData, "Cubic", Cubic )
         CALL SetFieldFromInput( InputData, "Quartic", Quartic )
         IF ( Quadratic /= 0.0 ) THEN 
            HarmonicFreq = SQRT( 2.0 * Quadratic / MassSystem )
         ELSE 
            HarmonicFreq = 1.0
         END IF
      END IF

      ! READ THE VARIABLES FOR THE PROPAGATION 

      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )

      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      TimeStep = TimeStep * TimeConversion(InputUnits, InternalUnits) 

      ! Nr of steps of the propagation
      CALL SetFieldFromInput( InputData, "NrOfSteps", NrOfSteps )

      ! Nr of steps in which the output is written 
      CALL SetFieldFromInput( InputData, "PrintStepInterval", PrintStepInterval )
      ! Accordingly set the interval between each printing step
      NrOfPrintSteps = CEILING( real(NrOfSteps) / real(PrintStepInterval) )

      ! READ THE VARIABLES TO SET THE EQUILIBRATION CONDITIONS

      ! Temperature of the bath
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 1.0 )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)

      ! Frequency and force constant of the harmonic force between the beads
      BeadsFrequency = NBeads * Temperature
      BeadsForceConst = ( NBeads * Temperature )**2

      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilRelaxTime", EquilGamma)
      EquilGamma = 1. / ( EquilGamma * TimeConversion(InputUnits, InternalUnits) )

      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep", EquilTStep, TimeStep*TimeConversion(InternalUnits,InputUnits) )
      EquilTStep = EquilTStep * TimeConversion(InputUnits, InternalUnits)

      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "EquilibrationTime", EquilibrationTime, &
                                                                  10.0/EquilGamma*TimeConversion(InternalUnits,InputUnits) )
      EquilibrationTime = EquilibrationTime * TimeConversion(InputUnits, InternalUnits)
      NrEquilibrationSteps = CEILING( EquilibrationTime / EquilTStep )

!       WRITE(*, 902) NBeads, BeadsFrequency/MyConsts_cmmin1toAU

!       WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV, &
!                     NrOfInitSnapshots, TimeBetweenSnaps/MyConsts_fs2AU

!       WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps

!       WRITE(*, 905) Temperature/MyConsts_K2AU, NrOfSnapshots, EquilTStep/MyConsts_fs2AU,  &
!                     EquilGamma*MyConsts_fs2AU, EquilibrTimeBetweenSnap


   902 FORMAT(" * Nr of replicas in the ring polymer dynamics: ", I10,/,   &
              " * Frequency of the force between the beads     ", F10.4 )

   903 FORMAT(" * Initial conditions of the atom-surface system ", /,&
              " * Absolute initial energy (eV):                ",F10.4,/,& 
              " *  - w.r.t. the bottom of the well (eV):       ",F10.4,/,& 
              " * Nr of initial system snapshots:              ",I10,  /,& 
              " * Time between snapshots (fs):                 ",F10.4,/ )

   904 FORMAT(" * Dynamical simulation variables               ", /,&
              " * Nr of trajectories:                          ",I10,  /,& 
              " * Propagation time step (fs):                  ",F10.4,/,& 
              " * Nr of time steps of each trajectory:         ",I10,  /,& 
              " * Nr of print steps of each trajectory:        ",I10,  / )

   905 FORMAT(" * Bath equilibration variables                 ", /,&
              " * Temperature of the surface:                  ",F10.4,/,&
              " * Nr of initial bath snapshots:                ",I10,  /,& 
              " * Equilibration time step (fs):                ",F10.4,/,& 
              " * Gamma of the Langevin force (fs^-1):         ",F10.4,/,& 
              " * Equilibration time between each snapshot:    ",F10.4,/ )

   END SUBROUTINE PolymerEquilibriumOscillator_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE PolymerEquilibriumOscillator_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord, iBead

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses

      NDim = 1
      ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads), APre(NDim*NBeads), MassVector(NDim*NBeads) )
      MassVector( : ) = MassSystem

      ! Set variables for EOM integration in the microcanonical ensamble
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector(1:NDim), TimeStep )
      ! Set ring polymer molecular dynamics parameter
      CALL SetupRingPolymer( MolecularDynamics, NBeads, BeadsFrequency ) 

      ! Set variables for EOM integration for the Bath in the canonical ensamble
      CALL EvolutionSetup( Equilibration, NDim*NBeads, MassVector, EquilTStep )
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature*NBeads )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the average energies over time only
      ! DEBUG - print also the average coordinates
      ! FULL - print detailed information about the initial conditions and about the trajectories

      ! Allocate and initialize the variables for the trajectory averages
      ALLOCATE( PositionCorrelation(0:NrOfPrintSteps) )
      PositionCorrelation(:)  = 0.0
      ! Allocate memory to compute the average energies (computed as centroid of classical quantities)
      ALLOCATE( KineticAverage(0:NrOfPrintSteps),  PotentialAverage(0:NrOfPrintSteps) )
      KineticAverage(:)  = 0.0
      PotentialAverage(:)  = 0.0

      ! Average coordinates over time and power spectrum of the auto-correlation function
      IF ( PrintType >= FULL ) THEN
         ALLOCATE( AverageCoord(NDim,0:NrOfPrintSteps) )
         AverageCoord(1:NDim,0:NrOfPrintSteps) = 0.0
      END IF

      IF ( PrintType >= FULL ) THEN
         ALLOCATE(  EquilibAveTvsTime( NrEquilibrationSteps/PrintStepInterval ) )
         EquilibAveTvsTime(:) = 0.0
      END IF

      ! Initialize random number seed
      CALL SetSeed( 10 )

   END SUBROUTINE PolymerEquilibriumOscillator_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerEquilibriumOscillator_Run()
      IMPLICIT NONE
      !> Output unit numbers
      INTEGER :: PosCorrelationUnit, AvCoordOutputUnit, TEquilUnit, EnergyOutputUnit, InitTUnit
      !> Debug unit numbers
      INTEGER :: DebugUnitEn, DebugUnitCoord, DebugUnitVel
      !> Temperature averages during equilibration
      REAL    :: TempAverage, TempVariance
      !> Energy averages
      REAL    :: TotEnergy, PotEnergy, KinEnergy, IstTemperature
      !> Initial energies of the bath
      REAL     ::  TrajKinAverage, TrajKinVariance, TrajPotAverage, TrajPotVariance
      REAL     ::  TotalKinAverage, TotalKinVariance, TotalPotAverage, TotalPotVariance
      !> Initial coordinates for the autocorrelation function
      REAL, DIMENSION(NDim*NBeads) :: X0
      !> Counters
      INTEGER :: iTraj, iStep, kStep, NInit
      !> Output file name
      CHARACTER(100) :: OutFileName

      REAL :: K1, K2, K3, K4
      REAL, DIMENSION(0:NrOfPrintSteps,4)    :: AverageKin
      AverageKin(:,:) = 0.0

      PosCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="PositionCorrelation.dat", UNIT=PosCorrelationUnit )
      WRITE(PosCorrelationUnit, "(3A,I6,A,/)") "# <X(0)X(t)> vs time (", TimeUnit(InputUnits), &
                                                                                  " | au) - ", NrTrajs, " trajectories "
      EnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="EnergyAverages.dat", UNIT=EnergyOutputUnit )
      WRITE(EnergyOutputUnit, "(5A,I6,A,/)") "# Ekin, Epot, ETot vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(EnergyUnit(InputUnits)), ") - ", NrTrajs, " trajectories "
      IF ( PrintType >= FULL ) THEN
         AvCoordOutputUnit = LookForFreeUnit()
         OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
         WRITE(AvCoordOutputUnit, "(5A,I6,A,/)") "# Coordinate vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(LengthUnit(InputUnits)), ") - ", NrTrajs, " trajectories "
         InitTUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationFinalT.dat", UNIT=InitTUnit )
         WRITE(InitTUnit, "(5A,I6,A,/)") "# <Temperature> vs time (", trim(TimeUnit(InputUnits)), ",", &
                                                             trim(TemperUnit(InputUnits)), ") - ", NrTrajs, " trajectories "
      ENDIF

      IF ( PrintType == DEBUG ) THEN
         TEquilUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationTemp.dat", UNIT=TEquilUnit )
         WRITE(TEquilUnit, "(A,I6,5A,/)") "# ", NrTrajs, " Temperature for each equilibration (", &
                                                            trim(TimeUnit(InputUnits)), ",", trim(TemperUnit(InputUnits)), ")"
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      ! Initialize averages of the bath initial state
      TotalKinAverage  = 0.0
      TotalKinVariance = 0.0
      TotalPotAverage  = 0.0
      TotalPotVariance = 0.0

      !run NrTrajs number of trajectories
      DO iTraj = 1,NrTrajs
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! set reasonable initial conditions for the equilibration via a normal mode calculations of the ring polymer
         CALL InitialConditionsForNormalModes( X, V ) 

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at Tp = ", &
                             NBeads*Temperature*TemperatureConversion(InternalUnits,InputUnits)

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = Potential( X, A )
         A(:) = A(:) / MassVector(:)

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         kStep = 0

         ! Long equilibration step before starting the initial condition sampling
         DO iStep = 1, NrEquilibrationSteps 

            ! PROPAGATION for ONE TIME STEP
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, Potential, PotEnergy )

            IF ( PrintType >= FULL ) THEN
               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(NDim*NBeads) 

               ! store temperature average and variance
               TempAverage = TempAverage + IstTemperature
               TempVariance = TempVariance + IstTemperature**2
            END IF

            ! every PrintStepInterval steps, write output
            IF ( mod(iStep-1,PrintStepInterval) == 0 ) THEN
               kStep = kStep + 1                                   ! Increment counter of print interval
                
               IF ( PrintType >= FULL ) THEN               ! compute average temperature
                  EquilibAveTvsTime(kStep) = EquilibAveTvsTime(kStep) + TempAverage/iStep
               END IF
               
               IF ( PrintType == DEBUG ) THEN              ! write debug output
                  ! Temperature profile during equilibration
                  IF ( iStep == 1 )  WRITE(TEquilUnit,*) " " 
                  WRITE(TEquilUnit,851)  real(iStep)*EquilTStep * TimeConversion(InternalUnits,InputUnits), &
                        TempAverage/iStep*TemperatureConversion(InternalUnits,InputUnits), &
                        sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)*TemperatureConversion(InternalUnits,InputUnits)
               END IF
               851 FORMAT( F20.5, 2F20.6 )
            END IF

         END DO

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! compute initial kinetic energy for this traj
         KinEnergy = EOM_KineticEnergy(Equilibration, V )
         PotEnergy = Potential( X, A )
         ! Shift potential energy with respect to the bottom of the well
         IF ( MorsePotential )   PotEnergy = PotEnergy + MorseDe*NBeads

         ! Increment averages of Kin and Pot for the single traj
         TrajKinAverage  = KinEnergy/(NDim*NBeads)
         TrajKinVariance = (KinEnergy/(NDim*NBeads))**2
         TrajPotAverage  = PotEnergy/(NDim*NBeads)
         TrajPotVariance = (PotEnergy/(NDim*NBeads))**2

         ! PRINT INITIAL CONDITIONS of THE BATH and THE SYSTEM
         WRITE(*,600)  KinEnergy/(NDim*NBeads)*EnergyConversion(InternalUnits,InputUnits), &
                       2.0*KinEnergy/(NDim*NBeads)*TemperatureConversion(InternalUnits,InputUnits), &
                       PotEnergy/(NDim*NBeads)*EnergyConversion(InternalUnits,InputUnits), &
                       2.0*PotEnergy/(NDim*NBeads)*TemperatureConversion(InternalUnits,InputUnits)

         ! Energy of the system, Coupling energy and energy of the bath
         CALL IstantaneousEnergies( KinEnergy, PotEnergy )
         CALL KineticEnergies( K1, K2, K3, K4 )

         ! Store starting values of the averages
         KineticAverage(0)    = KineticAverage(0)  + KinEnergy
         PotentialAverage(0)  = PotentialAverage(0) + PotEnergy

         ! Store initial coordinates
         X0 = X
         ! Compute position correlation function
         PositionCorrelation(0) = PositionCorrelation(0) + CorrelationFunction( X, X0 )

         AverageKin(0,1) = AverageKin(0,1) + K1 + PotEnergy
         AverageKin(0,2) = AverageKin(0,2) + K2 + PotEnergy
         AverageKin(0,3) = AverageKin(0,3) + K3 + PotEnergy
         AverageKin(0,4) = AverageKin(0,4) + K4 + PotEnergy

         ! Compute average coordinates
         IF ( PrintType >= FULL )  AverageCoord(1,0) = AverageCoord(1,0) + CentroidCoord( X )

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
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / au | Ekin, Epot, Etot / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KinEnergy, PotEnergy, KinEnergy+PotEnergy

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / au | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, X(:)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / au | V(1) V(2) ... V(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)

          ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kStep = 0

         PRINT "(/,A)", " Propagating the system in time... "
         
         ! Compute starting potential and forces
         CALL EOM_RPMSymplectic( MolecularDynamics, X, V, A,  SystemPotential, PotEnergy, .TRUE. )

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_RPMSymplectic( MolecularDynamics, X, V, A,  SystemPotential, PotEnergy )

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 
 
               ! compute kinetic energy for this traj
               KinEnergy = EOM_KineticEnergy(MolecularDynamics, V )
               ! Shift potential energy with respect to the bottom of the well
               IF ( MorsePotential )   PotEnergy = PotEnergy + MorseDe*NBeads
 
               ! Increment averages of Kin and Pot for the single traj
               TrajKinAverage  = TrajKinAverage + KinEnergy/(NDim*NBeads)
               TrajKinVariance = TrajKinVariance + (KinEnergy/(NDim*NBeads))**2
               TrajPotAverage  = TrajPotAverage + PotEnergy/(NDim*NBeads)
               TrajPotVariance = TrajPotVariance + (PotEnergy/(NDim*NBeads))**2

               ! Energy of the system
               CALL IstantaneousEnergies( KinEnergy, PotEnergy )
               CALL KineticEnergies( K1, K2, K3, K4 )

               ! Update averages over time
               KineticAverage(kStep)    = KineticAverage(kStep)  + KinEnergy
               PotentialAverage(kStep)  = PotentialAverage(kStep) + PotEnergy

               ! Compute position correlation function
               PositionCorrelation(kStep) = PositionCorrelation(kStep) + CorrelationFunction( X, X0 )

               AverageKin(kStep,1) = AverageKin(kStep,1) + K1 + PotEnergy
               AverageKin(kStep,2) = AverageKin(kStep,2) + K2 + PotEnergy
               AverageKin(kStep,3) = AverageKin(kStep,3) + K3 + PotEnergy
               AverageKin(kStep,4) = AverageKin(kStep,4) + K4 + PotEnergy

               IF ( PrintType >= FULL )   AverageCoord(1,kStep) = AverageCoord(1,kStep) + CentroidCoord( X )

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep), KinEnergy, PotEnergy, KinEnergy+PotEnergy
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep), X(:) 
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep), V(:)
               END IF

            END IF 

         END DO

         PRINT "(A)", " Time propagation completed! "

         IF ( PrintType == DEBUG ) THEN
               CLOSE( Unit=DebugUnitEn )
               CLOSE( Unit=DebugUnitCoord )
               CLOSE( Unit=DebugUnitVel )
         ENDIF

         TrajKinAverage  = TrajKinAverage / real(NrOfPrintSteps+1)
         TrajKinVariance = TrajKinVariance / real(NrOfPrintSteps+1) - TrajKinAverage**2
         TrajPotAverage  = TrajPotAverage / real(NrOfPrintSteps+1)
         TrajPotVariance = TrajPotVariance / real(NrOfPrintSteps+1) - TrajPotAverage**2

         ! PRINT AVERAGES FOR THE TRAJ
         WRITE(*,700)  TrajKinAverage*EnergyConversion(InternalUnits,InputUnits), &
                       2.0*TrajKinAverage*TemperatureConversion(InternalUnits,InputUnits), &
                       TrajPotAverage*EnergyConversion(InternalUnits,InputUnits), &
                       2.0*TrajPotAverage*TemperatureConversion(InternalUnits,InputUnits)

         ! Increment averages over the whole set of trajs
         TotalKinAverage  = TotalKinAverage  + TrajKinAverage
         TotalKinVariance = TotalKinVariance + TrajKinVariance
         TotalPotAverage  = TotalPotAverage  + TrajPotAverage
         TotalPotVariance = TotalPotVariance + TrajPotVariance

      END DO
      
      PRINT "(A)"," Done! "

      TotalKinAverage = TotalKinAverage / real(NrTrajs) 
      TotalKinVariance = SQRT( TotalKinVariance/real(NrTrajs) )
      TotalPotAverage = TotalPotAverage / real(NrTrajs) 
      TotalPotVariance = SQRT( TotalPotVariance/real(NrTrajs) )

      OPEN( Unit=999, FILE="FinalAverTemperature.log")
      WRITE(999,"(/,A,I10,/)") " Nr of trajectories ", NrTrajs
      WRITE(999,"(A)") " --------------------------------------------------- " 
      WRITE(999,*) " Average T : ",2.*TotalKinAverage/MyConsts_K2AU, &
                          " St dev : ", 2.*TotalKinVariance/MyConsts_K2AU
      WRITE(999,"(A,/)") " --------------------------------------------------- " 
      WRITE(999,"(A)") " --------------------------------------------------- " 
      WRITE(999,*) " Average V : ",2.*TotalPotAverage/MyConsts_K2AU, &
                         " St dev : ", 2.*TotalPotVariance/MyConsts_K2AU
      WRITE(999,"(A,/)") " --------------------------------------------------- " 
      CLOSE( Unit = 999 )

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize averages 
      PositionCorrelation(:) = PositionCorrelation(:)   /  real(NrTrajs)
      KineticAverage(:)      = KineticAverage(:)        / real(NrTrajs)
      PotentialAverage(:)    = PotentialAverage(:)      / real(NrTrajs)

      IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL )   EquilibAveTvsTime(:)   = EquilibAveTvsTime(:) / real(NrTrajs) 

      AverageKin(0:NrOfPrintSteps,1:4) = AverageKin(0:NrOfPrintSteps,1:4) / real(NrTrajs)

      DO iStep = 0, NrOfPrintSteps
         WRITE(900, "(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
               AverageKin(iStep,1:4)*MyConsts_Hartree2eV
      END DO

      ! Print average temperature over time during equilibration
      DO iStep = 0, size(EquilibAveTvsTime)-1
          WRITE(InitTUnit,"(2F20.8)") EquilTStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                  EquilibAveTvsTime(iStep+1)*TemperatureConversion(InternalUnits,InputUnits)
      END DO

      ! Print the position autocorrelation function
      DO iStep = 0, NrOfPrintSteps
          WRITE(PosCorrelationUnit,"(2F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, PositionCorrelation(iStep)
      END DO

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(EnergyOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                        PotentialAverage(iStep)*MyConsts_Hartree2eV, KineticAverage(iStep)*MyConsts_Hartree2eV
      END DO

      IF ( PrintType >= FULL ) THEN

         ! PRINT average coordinates
         DO iStep = 0, NrOfPrintSteps
            WRITE(AvCoordOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU,  &
                                                  AverageCoord(1,iStep)*MyConsts_Bohr2Ang
         END DO

      ENDIF

      ! Close output files
      CLOSE(PosCorrelationUnit)
      CLOSE(EnergyOutputUnit)
      IF ( PrintType >= FULL ) THEN
         CLOSE( AvCoordOutputUnit )
      END IF
 
      800 FORMAT(F12.5,1000F15.8)
      600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                     " * Kinetic Energy (incl. beads) (eV)    ",1F10.4,/    &
                     " * Translational temperature (K)        ",1F10.4,/    &
                     " * Potential Energy (incl. beads) (eV)  ",1F10.4,/    &
                     " * Vibrational temperature (K)          ",1F10.4,/ ) 
      700 FORMAT (/, " Avarages for the single MD trajectory ",/   &
                     " * Kinetic Energy (incl. beads) (eV)    ",1F10.4,/    &
                     " * Translational temperature (K)        ",1F10.4,/    &
                     " * Potential Energy (incl. beads) (eV)  ",1F10.4,/    &
                     " * Vibrational temperature (K)          ",1F10.4,/ ) 

   END SUBROUTINE PolymerEquilibriumOscillator_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerEquilibriumOscillator_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( X, V, A, APre, MassVector )
      DEALLOCATE( PositionCorrelation )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE PolymerEquilibriumOscillator_Dispose

!*************************************************************************************************

   FUNCTION CentroidCoord( X )
      IMPLICIT NONE
      REAL   :: CentroidCoord
      REAL, DIMENSION( NBeads ), INTENT(IN) :: X
      INTEGER :: i

      CentroidCoord = 0.0
      DO i = 1, NBeads
        CentroidCoord = CentroidCoord + X( i )
      END DO
      CentroidCoord = CentroidCoord / NBeads

   END FUNCTION CentroidCoord

!*************************************************************************************************

   REAL FUNCTION CorrelationFunction( X, X0 )
      IMPLICIT NONE
      REAL, DIMENSION( NBeads ), INTENT(IN) :: X, X0

      CorrelationFunction = CentroidCoord( X ) * CentroidCoord( X0 )
   END FUNCTION CorrelationFunction

!*************************************************************************************************
   
   SUBROUTINE KineticEnergies( K1, K2, K3, K4 )
      IMPLICIT NONE
      REAL, INTENT(OUT)             :: K1, K2, K3, K4
      REAL, DIMENSION(NDim*NBeads)  :: SngForce
      REAL, DIMENSION(NDim)         :: CentroidX, CentroidF, CentroidV, SingleF
      INTEGER                       :: iCoord, iBead
      REAL :: Pot, ForceSystem

      SngForce(:) = 0.0
      ! System potential for each bead
      DO iBead = 1, NBeads 
         ! Compute potential and forces of the system replicas
         IF ( MorsePotential ) THEN
            Pot = MorseV( X(iBead), ForceSystem )
         ELSE
            Pot = OscillatorV( X(iBead), ForceSystem )
         END IF
         SngForce(iBead) = SngForce(iBead) + ForceSystem
      END DO

      ! Compute centroids
      CentroidX = CentroidCoord( X )
      CentroidF = CentroidCoord( SngForce )
      CentroidV = CentroidCoord( V )

      ! First definition: kinetic energy of the centroid
      K1 = 0.0
      DO iCoord = 1, NDim
         K1 = K1 + 0.5 * MassVector(iCoord) * CentroidV(iCoord)**2
      END DO

      ! Second definition: centroid of kinetic energies
      K2 = 0.0
      DO iBead = 1, NBeads
         DO iCoord = 1, NDim
            K2  = K2  + 0.5 * MassVector(iCoord) * V((iBead-1)*NDim+iCoord)**2
         END DO
      END DO
      K2 = K2 / real(NBeads)

      ! Third definition: virial of the centroids (virial as correlation product)
      K3 = 0.0
      DO iCoord = 1, NDim
         K3  = K3  -0.5 * CentroidX(iCoord) * CentroidF(iCoord)
      END DO

      ! Fourth definition: centroid of the virial product (virial as average)
      K4 = 0.0
      DO iBead = 1, NBeads
         DO iCoord = 1, NDim
            K4  = K4  -0.5 * X((iBead-1)*NDim+iCoord) * SngForce((iBead-1)*NDim+iCoord)
         END DO
      END DO
      K4 = K4 / real(NBeads)

   END SUBROUTINE KineticEnergies

!*************************************************************************************************

   REAL FUNCTION Potential( Coordinates, Forces )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Coordinates
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces

      INTEGER                  :: iBead
      REAL                     :: ForceSystem
      REAL, DIMENSION(NBeads)  :: RingForces

      ! Check the number degrees of freedom
      CALL ERROR( size( Forces ) /= size( Coordinates ), "PolymerVibrationalRelax.Potential: array dimension mismatch" )
      CALL ERROR( size( Coordinates ) /= NDim*NBeads, "PolymerVibrationalRelax.Potential: wrong number of DoFs" )

      ! Initialize forces and potential
      Potential          = 0.0
      Forces(:)          = 0.0

      ! System potential for each bead
      DO iBead = 1, NBeads 
         ! Compute potential and forces of the system replicas
         IF ( MorsePotential ) THEN
            Potential = Potential+ MorseV( Coordinates(iBead), ForceSystem )
         ELSE
            Potential = Potential+ OscillatorV( Coordinates(iBead), ForceSystem )
         END IF
         Forces(iBead) = Forces(iBead) + ForceSystem
      END DO

      IF ( NBeads > 1 ) THEN
         ! Harmonic forces between beads
         CALL PolymerRingPotential( MassSystem, Coordinates, Potential, RingForces )
         Forces(:) = Forces(:) + RingForces(:)
      END IF
      
   END FUNCTION Potential

! ********************************************************************************************************************

   REAL FUNCTION SystemPotential( X, Force )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
      
      CALL ERROR( size( Force ) /= size( X ), "PolymerVibrationalRelax.SystemPotential: array dimension mismatch" )
      CALL ERROR( size( X ) /= NDim, "PolymerVibrationalRelax.SystemPotential: wrong number of DoFs" )

      ! Compute potential and forces of the system replicas
      IF ( MorsePotential ) THEN
         SystemPotential = MorseV( X(1), Force(1) )
      ELSE
         SystemPotential = OscillatorV( X(1), Force(1) )
      END IF
   END FUNCTION SystemPotential

   REAL FUNCTION OscillatorV( X, dVdX )
      IMPLICIT NONE
      REAL, INTENT(IN)  :: X
      REAL, INTENT(OUT), OPTIONAL :: dVdX 

      OscillatorV = Quadratic * X**2 + Cubic * X**3 + Quartic * X**4
      IF (PRESENT(dVdX) )  dVdX = - 2.0 * Quadratic * X - 3.0 * Cubic * X**2 - 4.0 * Quartic * X**3
   END FUNCTION OscillatorV


   REAL FUNCTION MorseV( Positions, Forces ) RESULT(V) 
      IMPLICIT NONE
      REAL, INTENT(IN)  :: Positions
      REAL, INTENT(OUT), OPTIONAL  :: Forces 

      V = MorseDe * ( exp(-2.0*MorseAlpha*Positions) - 2.0 * exp(-MorseAlpha*Positions) )  
      IF (PRESENT(Forces) ) Forces = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*Positions) - exp(-MorseAlpha*Positions) )  

   END FUNCTION MorseV

! ********************************************************************************************************************

   SUBROUTINE PolymerRingPotential( Mass, X, V, ForceOnX )
      IMPLICIT NONE
      REAL                               :: Mass
      REAL, DIMENSION(:), INTENT(IN)     :: X
      REAL, DIMENSION(:), INTENT(OUT)    :: ForceOnX
      REAL, INTENT(INOUT)                :: V
      INTEGER :: i, NRing

      NRing = size(X)

      IF ( NRing == 1 ) THEN
         V = V
         ForceOnX = 0.0
      ELSE IF ( NRing == 2 ) THEN
         V = V + 0.5 * Mass * BeadsForceConst * ( X(2) - X(1) )**2 
         ForceOnX(1) = Mass * BeadsForceConst * ( X(2) - X(1) )  
         ForceOnX(2) = Mass * BeadsForceConst * ( - X(2) + X(1) )  
      ELSE
         V = V + 0.5 * Mass * BeadsForceConst * ( X(2) - X(1) )**2 
         ForceOnX(1) = - Mass * BeadsForceConst * ( - X(NRing) + X(1) * 2 - X(2) )  
         DO i = 2, NRing-1
            V = V + 0.5 * Mass * BeadsForceConst * ( X(i+1) - X(i) )**2 
            ForceOnX(i) = - Mass * BeadsForceConst * ( - X(i-1) + X(i) * 2 - X(i+1) )  
         END DO
         V = V + 0.5 * Mass * BeadsForceConst * ( X(1) - X(NRing) )**2 
         ForceOnX(NRing) = - Mass * BeadsForceConst * ( - X(NRing-1) + X(NRing) * 2 - X(1) )  
      END IF

   END SUBROUTINE PolymerRingPotential

!*************************************************************************************************

   SUBROUTINE IstantaneousEnergies( KEn, VEn )
      IMPLICIT NONE
      REAL, INTENT(OUT) :: KEn, VEn
      REAL  :: CentroidV
      INTEGER ::  iBead

      ! Centroid of the velocities
      CentroidV = CentroidCoord( V )

      ! Compute kinetic energy as the kinetic energy of the centroid
      KEn =  0.5 * MassSystem * CentroidV**2

!       ! Compute kinetic energy as average of the kinetic energy over the replicas
!       KEn = 0.0
!       DO iBead = 1, NBeads
!             KEn = KEn + 0.5 * MassSystem * V( iBead )**2
!       END DO
!       KEn = KEn / NBeads 

      ! Energy of the system - average of the potential energies of the system replicas
      VEn = 0.0
      DO iBead = 1, NBeads
         IF ( MorsePotential ) THEN
            VEn = VEn + MorseV( X(iBead) ) + MorseDe
         ELSE
            VEn = VEn + OscillatorV( X(iBead) )
         END IF
      END DO
      VEn = VEn / NBeads 
      
 
   END SUBROUTINE IstantaneousEnergies


!*************************************************************************************************

   SUBROUTINE InitialConditionsForNormalModes( Pos, Vel )
      IMPLICIT NONE
      REAL, DIMENSION(NDim*NBeads), INTENT(OUT) :: Pos, Vel 

      REAL, DIMENSION(NBeads,NBeads) :: RingPotentialMatrix, RingNormalModes
      REAL, DIMENSION(NBeads)        :: RingEigenvalues
      REAL, DIMENSION(NDim,NDim)     :: SysPotentialMatrix, SysNormalModes
      REAL, DIMENSION(NDim)          :: SysEigenvalues

      REAL, DIMENSION(NDim,NBeads) :: NormalQ, NormalV 
      REAL, DIMENSION(NDim,NBeads) :: FinalQ,  FinalV 
      INTEGER :: iSys, iBead, jSys, jBead
      REAL    :: SigmaV

      ! Setup system potential matrix
      SysPotentialMatrix(1,1) = HarmonicFreq**2
      SysNormalModes(1,1) = 1.0
      SysEigenvalues(1) = SysPotentialMatrix(1,1)

      ! Setup ring potential matrix
      RingPotentialMatrix(:,:) = 0.0
      DO iBead = 1, NBeads-1
         RingPotentialMatrix(iBead,iBead)   = 2.0 * BeadsFrequency**2
         RingPotentialMatrix(iBead+1,iBead) = - BeadsFrequency**2
         RingPotentialMatrix(iBead,iBead+1) = - BeadsFrequency**2
      END DO
      RingPotentialMatrix(NBeads,NBeads) = 2.0 * BeadsFrequency**2
      RingPotentialMatrix(1,NBeads) = - BeadsFrequency**2
      RingPotentialMatrix(NBeads,1) = - BeadsFrequency**2

      ! Diagonalize ring potential
      CALL TheOneWithDiagonalization(RingPotentialMatrix, RingNormalModes, RingEigenvalues)

      ! Set sigma of the maxwell boltzmann distribution
      SigmaV = sqrt( NBeads*Temperature / MassSystem )

      ! Define random coordinates and velocities in the normal modes representation
      DO iBead = 1, NBeads
         DO iSys = 1, NDim
            NormalQ(iSys,iBead) = GaussianRandomNr( 1.0 ) * SigmaV / SQRT( SysEigenvalues(iSys) + RingEigenvalues(iBead) )
            NormalV(iSys,iBead) = GaussianRandomNr( 1.0 ) * SigmaV            
         END DO
      END DO

      ! Transform coordinates to original representation
      DO iSys = 1, NDim
         DO iBead = 1, NBeads
            FinalQ(iSys,iBead) = 0.0
            FinalV(iSys,iBead) = 0.0
            DO jSys = 1, NDim
               DO jBead = 1, NBeads
                  FinalQ(iSys,iBead) = FinalQ(iSys,iBead) + &
                                          SysNormalModes(iSys,jSys) * RingNormalModes(iBead,jBead) * NormalQ(jSys,jBead)
                  FinalV(iSys,iBead) = FinalV(iSys,iBead) + &
                                          SysNormalModes(iSys,jSys) * RingNormalModes(iBead,jBead) * NormalV(jSys,jBead)
               END DO
            END DO
         END DO
      END DO

      ! Give back the initial conditions
      DO iBead = 1, NBeads
         Pos( (iBead-1)*NDim+1 : iBead*NDim ) = FinalQ(:,iBead)
         Vel( (iBead-1)*NDim+1 : iBead*NDim ) = FinalV(:,iBead)
      END DO

   END SUBROUTINE InitialConditionsForNormalModes

!*************************************************************************************************

END MODULE PolymerEquilibriumOscillator