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
   USE MyConsts
   USE ErrorTrap
   USE SharedData
   USE InputField
   USE ClassicalEqMotion
   USE RandomNumberGenerator


   IMPLICIT NONE

   PRIVATE

   PUBLIC :: PolymerEquilibriumOscillator_ReadInput, PolymerEquilibriumOscillator_Initialize, &
             PolymerEquilibriumOscillator_Run, PolymerEquilibriumOscillator_Dispose

   ! Input is in atomic units
   LOGICAL :: AtomicUnits

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
   TYPE(Evolution) :: MolecularDynamics     !< Propagate in micro/canonical ensamble to extract results
   TYPE(Evolution) :: Equilibration         !< Propagate in canonical ensamble to generate initial conditions of the bath

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: PositionCorrelation    !< Position correlation function
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord           !< Average i-th coordinate vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: KineticAverage         !< Average of the kinetic energy of the system
   REAL, DIMENSION(:), ALLOCATABLE      :: PotentialAverage       !< Average of the potential energy of the system
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageVirial          !< Average kinetic energy via virial theorem
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

      ! Nr of replicas of the system for the ring polymer dynamics
      CALL SetFieldFromInput( InputData, "AtomicUnits", AtomicUnits, .FALSE. )

      ! PARAMETERS FOR THE SYSTEM DEFINITION 

      ! Nr of replicas of the system for the ring polymer dynamics
      CALL SetFieldFromInput( InputData, "NBeads", NBeads )

      ! Mass of the system ( in amu )
      CALL SetFieldFromInput( InputData, "MassSystem", MassSystem )
      IF ( .NOT. AtomicUnits ) MassSystem = MassSystem * MyConsts_Uma2Au

      ! Potential is a morse function?
      CALL SetFieldFromInput( InputData, "MorsePotential", MorsePotential )

      IF ( MorsePotential ) THEN
         CALL SetFieldFromInput( InputData, "MorseAlpha", MorseAlpha )
         CALL SetFieldFromInput( InputData, "MorseDe", MorseDe )
         IF ( .NOT. AtomicUnits ) MorseDe = MorseDe / MyConsts_Hartree2eV
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

      ! READ THE VARIABLES TO SET THE EQUILIBRATION CONDITIONS

      ! Temperature of the bath
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 1.0 )
      IF ( .NOT. AtomicUnits ) Temperature = Temperature * MyConsts_K2AU

      ! Frequency and force constant of the harmonic force between the beads
      BeadsFrequency = NBeads * Temperature
      BeadsForceConst = ( NBeads * Temperature )**2

      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilGamma", EquilGamma)
      IF ( .NOT. AtomicUnits ) EquilGamma = EquilGamma / MyConsts_fs2AU

      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
      IF ( .NOT. AtomicUnits ) EquilTStep = EquilTStep * MyConsts_fs2AU

      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "EquilibrationTime", EquilibrationTime, 10.0*(1.0/EquilGamma)/MyConsts_fs2AU )
      IF ( .NOT. AtomicUnits ) EquilibrationTime = EquilibrationTime * MyConsts_fs2AU
      NrEquilibrationSteps = CEILING( EquilibrationTime / EquilTStep )

      ! READ THE VARIABLES FOR THE PROPAGATION 

      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )

      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      IF ( .NOT. AtomicUnits ) TimeStep = TimeStep * MyConsts_fs2AU

      ! Nr of steps of the propagation
      CALL SetFieldFromInput( InputData, "NrOfSteps", NrOfSteps )

      ! Nr of steps in which the output is written 
      CALL SetFieldFromInput( InputData, "PrintStepInterval", PrintStepInterval )
      ! Accordingly set the interval between each printing step
      NrOfPrintSteps = CEILING( real(NrOfSteps) / real(PrintStepInterval) )


      WRITE(*, 902) NBeads, BeadsFrequency/MyConsts_cmmin1toAU

!       WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV, &
!                     NrOfInitSnapshots, TimeBetweenSnaps/MyConsts_fs2AU

      WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps

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

      ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
      CALL EvolutionSetup( MolecularDynamics, NDim*NBeads, MassVector, TimeStep )

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
      ! Allocate memory to compute the average kinetic energy as centroid of the virial of the forces
      ALLOCATE( AverageVirial(0:NrOfPrintSteps) )
      AverageVirial(:) = 0.0

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
      !> Variance of position and velocity for initial guess of the coordinate
      REAL    :: SigmaX, SigmaV
      !> Initial coordinates for the autocorrelation function
      REAL, DIMENSION(NDim*NBeads) :: X0
      !> Counters
      INTEGER :: iTraj, iStep, kStep, NInit
      !> Output file name
      CHARACTER(100) :: OutFileName

      PosCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="PositionCorrelation.dat", UNIT=PosCorrelationUnit )
      WRITE(PosCorrelationUnit, "(A,I6,A,/)") "# X(0)X(t) correlation vs time (fs | au) - ", NrTrajs, " trajectories "

      EnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="EnergyAverages.dat", UNIT=EnergyOutputUnit )
      WRITE(EnergyOutputUnit, "(A,I6,A,/)") "# Energy  vs time (au) | EKIN - EV - ETOT - ", NrTrajs, " trajectories "

      IF ( PrintType >= FULL ) THEN
         AvCoordOutputUnit = LookForFreeUnit()
         OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
         WRITE(AvCoordOutputUnit, "(A,I6,A,/)") "# average coordinate vs time (fs | Ang) - ", NrTrajs, " trajectories "
         InitTUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationFinalT.dat", UNIT=InitTUnit )
         WRITE(InitTUnit, "(A,I6,A,/)") "# average T over all the trajectories vs Time (K | fs )", NrTrajs, " trajectories "
      ENDIF

      IF ( PrintType == DEBUG ) THEN
         TEquilUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationTemp.dat", UNIT=TEquilUnit )
         WRITE(TEquilUnit, "(A,I6,A,/)") "# ", NrTrajs, " temperature average and variance for each equilibration (fs | K)"
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      SigmaX = sqrt( Temperature / MassSystem ) / HarmonicFreq
      SigmaV = sqrt( Temperature / MassSystem )

      !run NrTrajs number of trajectories
      DO iTraj = 1,NrTrajs
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! set reasonable initial conditions for the equilibration
         ! The position and velocity variance is estimated via a harmonic approximation of the potential
         V(:) = GaussianRandomNr( 1.0 ) * SigmaV
         X(:) = GaussianRandomNr( 1.0 ) * SigmaX

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at Tp = ", NBeads*Temperature/MyConsts_K2AU

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

            ! PROPAGATION for ONE TIME STEP - if first step, previous acceleration are not available, use VV
!             IF ( iTraj == 1 ) THEN
!                ! Store initial accelerations
!                APre(:) = A(:)
!                ! Propagate for one timestep with Velocity-Verlet and langevin thermostat
!                CALL EOM_VelocityVerlet( Equilibration, X, V, A, Potential, PotEnergy )
!             ELSE
!                ! Propagate for one timestep with Beeman's method and langevin thermostat
!                CALL EOM_Beeman( Equilibration, X, V, A, APre, Potential, PotEnergy )
!             END IF
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, Potential, PotEnergy )

            IF ( PrintType >= FULL ) THEN
               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IF ( AtomicUnits ) THEN 
                  IstTemperature = 2.0*KinEnergy/(NDim*NBeads)
               ELSE
                  IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*NDim*NBeads)
               END IF

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
                  IF ( AtomicUnits ) THEN
                     WRITE(TEquilUnit,851)  real(iStep)*EquilTStep, &
                                                               TempAverage/iStep, sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
                  ELSE 
                     WRITE(TEquilUnit,851)  real(iStep)*EquilTStep/MyConsts_fs2AU, &
                                                               TempAverage/iStep, sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
                  END IF
               END IF
               851 FORMAT( F20.5, 2F20.6 )
            END IF

         END DO

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Energy of the system, Coupling energy and energy of the bath
         CALL IstantaneousEnergies( KinEnergy, PotEnergy )

         ! Store starting values of the averages
         KineticAverage(0)    = KineticAverage(0)  + KinEnergy
         PotentialAverage(0)  = PotentialAverage(0) + PotEnergy

         ! Store initial coordinates
         X0 = X
         ! Compute position correlation function
         PositionCorrelation(0) = PositionCorrelation(0) + CorrelationFunction( X, X0 )
         ! Compute virial of the force 
         AverageVirial(0) = AverageVirial(0) + Virial(X)

         ! Compute average coordinates
         IF ( PrintType >= FULL )  AverageCoord(1,0) = AverageCoord(1,0) + CentroidCoord( X )

!          ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
!          WRITE(*,600)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature

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

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | V(1) V(2) ... V(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)

          ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kStep = 0

         PRINT "(/,A)", " Propagating the system in time... "
         
         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = Potential( X, A )
         A(:) = A(:) / MassVector(:)

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, Potential, PotEnergy )

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 
 
               ! Energy of the system
               CALL IstantaneousEnergies( KinEnergy, PotEnergy )

               ! Update averages over time
               KineticAverage(kStep)    = KineticAverage(kStep)  + KinEnergy
               PotentialAverage(kStep)  = PotentialAverage(kStep) + PotEnergy

               ! Compute position correlation function
               PositionCorrelation(kStep) = PositionCorrelation(kStep) + CorrelationFunction( X, X0 )
               ! Compute virial of the force 
               AverageVirial(kStep) = AverageVirial(kStep) + Virial(X)

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

!          ! Print log information about the final condition of the trajectory
!          WRITE(*,601)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature

         IF ( PrintType == DEBUG ) THEN
               CLOSE( Unit=DebugUnitEn )
               CLOSE( Unit=DebugUnitCoord )
               CLOSE( Unit=DebugUnitVel )
         ENDIF

      END DO
      
      PRINT "(A)"," Done! "

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize averages 
      PositionCorrelation(:) = PositionCorrelation(:)   /  real(NrTrajs)
      KineticAverage(:)      = KineticAverage(:)        / real(NrTrajs)
      PotentialAverage(:)    = PotentialAverage(:)      / real(NrTrajs)
      AverageVirial(:)       = - AverageVirial(:) * 0.5 / real(NrTrajs)

      IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL )   EquilibAveTvsTime(:)   = EquilibAveTvsTime(:) / real(NrTrajs) 

      ! Print average temperature over time during equilibration
      DO iStep = 0, size(EquilibAveTvsTime)-1
          WRITE(InitTUnit,"(2F20.8)") EquilTStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, EquilibAveTvsTime(iStep+1)
      END DO

      ! Print the position autocorrelation function
      DO iStep = 0, NrOfPrintSteps
          WRITE(PosCorrelationUnit,"(2F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, PositionCorrelation(iStep)
      END DO

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(EnergyOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                       AverageVirial(iStep)*MyConsts_Hartree2eV,  PotentialAverage(iStep)*MyConsts_Hartree2eV, &
                                                                                       KineticAverage(iStep)*MyConsts_Hartree2eV
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
                     " * Energy of the system (eV)        ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)  ",1F10.4,/    &
                     " * Istantaneous temperature (K)     ",1F10.4,/ ) 
      601 FORMAT (/, " Final condition of the MD trajectory ",/   &
                     " * Energy of the system (eV)        ",1F10.4,/    &
                     " * Kinetic Energy of the bath (eV)  ",1F10.4,/    &
                     " * Istantaneous temperature (K)     ",1F10.4,/ ) 

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

   REAL FUNCTION Virial( X )
      IMPLICIT NONE
      REAL, DIMENSION( NBeads ), INTENT(IN) :: X
      REAL, DIMENSION( NBeads ) :: Forces
      REAL :: V, ForceSystem
      INTEGER :: iBead
   
      ! Initialize forces and potential
      V          = 0.0
      Forces(:)  = 0.0

      ! System potential for each bead
      DO iBead = 1, NBeads 
         ! Compute potential and forces of the system replicas
         IF ( MorsePotential ) THEN
            V = V + MorseV( X(iBead), ForceSystem )
         ELSE
            V = V + OscillatorV( X(iBead), ForceSystem )
         END IF
         Forces(iBead) = Forces(iBead) + ForceSystem
      END DO

!       IF ( NBeads > 1 ) THEN
!          ! Harmonic forces between beads
!          CALL PolymerRingPotential( MassSystem, Coordinates, Potential, RingForces )
!          Forces(:) = Forces(:) + RingForces(:)
!       END IF
      
      Virial = 0.0
      DO iBead = 1, NBeads
         Virial = Virial + X(iBead) * Forces(iBead)
      END DO
      Virial = Virial / NBeads

!       Virial = CentroidCoord( X ) * CentroidCoord( Forces )

   END FUNCTION Virial

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
            VEn = VEn + MorseV( X(iBead) )
         ELSE
            VEn = VEn + OscillatorV( X(iBead) )
         END IF
      END DO
      VEn = VEn / NBeads 
      
 
   END SUBROUTINE IstantaneousEnergies

!*************************************************************************************************

END MODULE PolymerEquilibriumOscillator