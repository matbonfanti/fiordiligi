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

      IF ( MorsePotential ) THEN
         WRITE(*,800) MorseDe*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),            &
                      MorseAlpha/LengthConversion(InternalUnits,InputUnits), "1/"//LengthUnit(InputUnits),   &
                      HarmonicFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),           &
                      MassSystem*MassConversion(InternalUnits,InputUnits), MassUnit(InputUnits)
      ELSE
         WRITE(*,800) MorseDe*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),             &
                      Quadratic, "au", Cubic, "au", Quartic, "au",                                            &
                      MassSystem*MassConversion(InternalUnits,InputUnits), MassUnit(InputUnits)
      END IF

      WRITE(*, 901) Temperature*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits),     &
                    CEILING( HarmonicFreq / Temperature )

      WRITE(*, 902) NBeads, BeadsFrequency*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                    2.*BeadsFrequency*SIN(MyConsts_PI/real(NBeads))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)

      WRITE(*, 903) NrTrajs, TimeStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrOfSteps, NrOfPrintSteps

      WRITE(*, 904) 1./EquilGamma*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    EquilibrationTime*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), NrEquilibrationSteps

   800 FORMAT(" * Morse potential simulation                   ", /,   &
              " * Depth of the well:                           ",F10.4,1X,A,/,&  
              " * Length exponential parameter:                ",F10.4,1X,A,/,&
              " * Harmonic frequency at the minimum:           ",F10.4,1X,A,/,&
              " * Mass of the oscillator:                      ",F10.4,1X,A,/ )
   801 FORMAT(" * Polynomial potential simulation              ", /,   &
              " * Quadratic term:                              ",F10.4,1X,A,/,&
              " * Cubic term:                                  ",F10.4,1X,A,/,&
              " * Quartic term:                                ",F10.4,1X,A,/,&
              " * Mass of the oscillator:                      ",F10.4,1X,A,/ )

   901 FORMAT(" * Temperature of the simulation:               ",F10.4,1X,A,/,&
              " * Minimum nr of beads at this temperature:     ",I10,/ )

   902 FORMAT(" * Nr of replicas in the ring polymer dynamics: ", I10,       /,&
              " * Frequency of the force between the beads:    ", F10.4,1X,A,/,&
              " * Frequency of 2nd normal mode of the ring:    ", F10.4,1X,A,/ )

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
      ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads), MassVector(NDim) )
      MassVector( : ) = MassSystem

      ! Set variables for EOM integration in the microcanonical ensamble
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector, TimeStep )
      ! Set ring polymer molecular dynamics parameter
      CALL SetupRingPolymer( MolecularDynamics, NBeads, BeadsFrequency ) 

      ! Set variables for EOM integration for the Bath in the canonical ensamble
      CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
      ! Set ring polymer molecular dynamics parameter
      CALL SetupRingPolymer( Equilibration, NBeads, BeadsFrequency ) 
      ! Set thermostatting with PILE
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the average energies over time only
      ! DEBUG - print also the average coordinates
      ! FULL - print detailed information about the initial conditions and about the trajectories

      ! Allocate and initialize the variables for the trajectory averages
      ALLOCATE( PositionCorrelation(0:NrOfPrintSteps) )
      PositionCorrelation(:)  = 0.0

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

      INTEGER, PARAMETER   :: NrOfEnergyAverages = 4
      REAL, DIMENSION(NrOfEnergyAverages,0:NrOfPrintSteps)    :: AverageE
      AverageE(:,:) = 0.0

      PosCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="PositionCorrelation.dat", UNIT=PosCorrelationUnit )
      WRITE(PosCorrelationUnit, "(3A,I6,A,/)") "# <X(0)X(t)> vs time (", TimeUnit(InputUnits), &
                                                                                  " | au) - ", NrTrajs, " trajectories "
      EnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="EnergyAverages.dat", UNIT=EnergyOutputUnit )
      WRITE(EnergyOutputUnit, "(5A,I6,A,/)") "# Evir, Etherm, Ecentr, V vs time (", trim(TimeUnit(InputUnits)), ",", &
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

         PRINT "(/,A,F12.1)"," Equilibrating the initial conditions at Tp = ", &
                             NBeads*Temperature*TemperatureConversion(InternalUnits,InputUnits)

         ! Compute starting potential and forces
         CALL EOM_RPMSymplectic( Equilibration, X, V, A,  SystemPotential, PotEnergy, .TRUE. )

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         kStep = 0

         ! Long equilibration step before starting the initial condition sampling
         DO iStep = 1, NrEquilibrationSteps 

            ! PROPAGATION for ONE TIME STEP
            CALL EOM_RPMSymplectic( Equilibration, X, V, A,  SystemPotential, PotEnergy )

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
         CALL EOM_RPMSymplectic( Equilibration, X, V, A,  SystemPotential, PotEnergy, .TRUE. )
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

         ! Compute and store expectation values of the energy
         AverageE(:,0) = AverageE(:,0) + EnergyAverages( X, V, MassVector )

         ! Store initial coordinates
         X0 = X
         ! Compute position correlation function
         PositionCorrelation(0) = PositionCorrelation(0) + CorrelationFunction( X, X0 )

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

               ! Compute and store expectation values of the energy
               AverageE(:,kStep) = AverageE(:,kStep) + EnergyAverages( X, V, MassVector ) 

               ! Compute position correlation function
               PositionCorrelation(kStep) = PositionCorrelation(kStep) + CorrelationFunction( X, X0 )

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
      AverageE(:,0:NrOfPrintSteps) = AverageE(:,0:NrOfPrintSteps)  / real(NrTrajs)
      IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL )   EquilibAveTvsTime(:)   = EquilibAveTvsTime(:) / real(NrTrajs) 

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
         WRITE(EnergyOutputUnit,"(F14.8,100F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU,  &
                              AverageE(:,iStep)*MyConsts_Hartree2eV
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
      DEALLOCATE( X, V, A, MassVector )
      DEALLOCATE( PositionCorrelation )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE PolymerEquilibriumOscillator_Dispose

! ********************************************************************************************************************

   REAL FUNCTION SystemPotential( X, Force )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
      
      CALL ERROR( size( Force ) /= size( X ), "PolymerVibrationalRelax.SystemPotential: array dimension mismatch" )
      CALL ERROR( size( X ) /= NDim, "PolymerVibrationalRelax.SystemPotential: wrong number of DoFs" )

      ! Compute potential and forces of the system replicas
      IF ( MorsePotential ) THEN
         SystemPotential = MorseDe * ( 1.0 - exp(-MorseAlpha*X(1)) )**2  
         Force(1) = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*X(1)) - exp(-MorseAlpha*X(1)) )  
      ELSE
         SystemPotential = Quadratic * X(1)**2 + Cubic * X(1)**3 + Quartic * X(1)**4
         Force(1) = - 2.0 * Quadratic * X(1) - 3.0 * Cubic * X(1)**2 - 4.0 * Quartic * X(1)**3
      END IF

   END FUNCTION SystemPotential

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
   
   FUNCTION EnergyAverages( X, V, Mass ) 
      IMPLICIT NONE
      REAL, DIMENSION(4) :: EnergyAverages
      REAL, DIMENSION(NDim*NBeads), INTENT(IN) :: X, V
      REAL, DIMENSION(NDim), INTENT(IN)        :: Mass
      REAL                          :: PotEnergy
      REAL, DIMENSION(NDim*NBeads)  :: SngForce
      REAL, DIMENSION(NDim)         :: CentroidV
      INTEGER                       :: iCoord, iBead

      ! (4) POTENTIAL ENERGY    *****************************************************
      ! System potential averaged over the beads, and force acting on each bead
      PotEnergy = 0.0
      DO iBead = 1, NBeads 
         PotEnergy = PotEnergy + SystemPotential( X((iBead-1)*NDim+1:iBead*NDim), SngForce((iBead-1)*NDim+1:iBead*NDim) )
      END DO
      EnergyAverages(4) = PotEnergy / real(NBeads)

      ! (1) VIRIAL TOTAL ENERGY    **************************************************
      ! centroid of the virial product (virial as average)
      EnergyAverages(1) = 0.0
      DO iBead = 1, NBeads
         DO iCoord = 1, NDim
            EnergyAverages(1) = EnergyAverages(1) - X((iBead-1)*NDim+iCoord) * SngForce((iBead-1)*NDim+iCoord)
         END DO
      END DO
      EnergyAverages(1) = 0.5 * EnergyAverages(1) / real(NBeads)
      ! add potential energy to get full energy
      EnergyAverages(1) = EnergyAverages(1) + EnergyAverages(4)

      ! 2) THERMODYNAMICS TOTAL ENERGY
      ! sum up energy of the ring oscillators and subtract it to average energy
      EnergyAverages(2) = 0.0 
      DO iCoord = 1, NDim
         IF ( NBeads == 2 ) THEN
            EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(NDim+iCoord) - X(iCoord) )**2 
         ELSE IF ( NBeads > 2 ) THEN
            DO iBead = 1, NBeads-1
               EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iBead*NDim+iCoord) - X((iBead-1)*NDim+iCoord) )**2 
            END DO
            EnergyAverages(2) = EnergyAverages(2) + Mass(iCoord) * ( X(iCoord) - X((NBeads-1)*NDim+iCoord) )**2 
         END IF
      END DO
      EnergyAverages(2) = 0.5 * ( NBeads**2 * NDim * Temperature - BeadsForceConst * EnergyAverages(2) )
      ! add potential energy to get full energy
      EnergyAverages(2) = EnergyAverages(2) + EnergyAverages(4)

      ! 3) POTENTIAL PLUS KINETIC ENERGY OF THE CENTROID
      ! Compute centroid of the velocity
      CentroidV = CentroidCoord( V )
      EnergyAverages(3) = 0.0
      DO iCoord = 1, NDim
         EnergyAverages(3) = EnergyAverages(3) + MassVector(iCoord) * CentroidV(iCoord)**2
      END DO
      EnergyAverages(3) =  0.5 * EnergyAverages(3)
      ! add potential energy to get full energy
      EnergyAverages(3) = EnergyAverages(3) + EnergyAverages(4)


   END FUNCTION EnergyAverages

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