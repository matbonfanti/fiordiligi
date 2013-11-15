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
!>  \todo          implement vibrational relaxation of the morse oscillator
!>                 
!***************************************************************************************
MODULE PolymerVibrationalRelax
   USE MyConsts
   USE ErrorTrap
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: PolymerVibrationalRelax_ReadInput, PolymerVibrationalRelax_Initialize, &
             PolymerVibrationalRelax_Run, PolymerVibrationalRelax_Dispose

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
   REAL    :: EquilTStep              !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma              !< Friction parameter of the Langevin equation
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
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageEBath       !< Average energy of the bath vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageECoup       !< Average coupling energy vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageESys        !< Average energy of the system vs time
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord       !< Average i-th coordinate vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: EquilibAveTvsTime  !< Average temperature at time T during equilibration

   REAL, DIMENSION(:), ALLOCATABLE      :: DissipatedPower    !< 


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
         CALL SetFieldFromInput( InputData, "MorseDe", MorseDe )
         MorseDe = MorseDe / MyConsts_Hartree2eV
         HarmonicFreq = SQRT( 2.0 * MorseAlpha**2 * MorseDe / MassH )
      END IF

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Initial energy of the system (input in eV and transform to AU)
      CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
      InitEnergy = InitEnergy / MyConsts_Hartree2eV

      ! READ THE VARIABLES TO SET THE EQUILIBRATION CONDITIONS

      ! Temperature of the bath
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 1.0 )
      Temperature = Temperature * MyConsts_K2AU

      ! Frequency and force constant of the harmonic force between the beads
      BeadsFrequency = NBeads * Temperature
      BeadsForceConst = ( NBeads * Temperature )**2

      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilGamma", EquilGamma)
      EquilGamma = EquilGamma / MyConsts_fs2AU

      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
      EquilTStep = EquilTStep * MyConsts_fs2AU

      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "EquilibrationTime", EquilibrationTime, 10.0*(1.0/EquilGamma)/MyConsts_fs2AU )
      EquilibrationTime = EquilibrationTime * MyConsts_fs2AU
      NrEquilibrationSteps = CEILING( EquilibrationTime / EquilTStep )

      ! READ THE VARIABLES FOR THE RELAXATION SIMULATION 

      ! Nr of the trajectories of the simulation
      CALL SetFieldFromInput( InputData, "NrTrajs", NrTrajs )

      ! Timestep of the propagation
      CALL SetFieldFromInput( InputData, "TimeStep", TimeStep )
      TimeStep = TimeStep * MyConsts_fs2AU

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

      WRITE(*, 905) Temperature/MyConsts_K2AU, EquilibrationTime/MyConsts_fs2AU, &
                    EquilTStep/MyConsts_fs2AU, EquilGamma*MyConsts_fs2AU


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
              " * Temperature of the surface:                  ",F10.4,/,&
              " * Equilibration time (fs):                     ",F10.4,/,& 
              " * Equilibration time step (fs):                ",F10.4,/,& 
              " * Gamma of the Langevin force (fs^-1):         ",F10.4,/ )

   END SUBROUTINE PolymerVibrationalRelax_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord, iBead
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn
      REAL, DIMENSION(:), ALLOCATABLE    :: SystemMasses

      CALL ERROR( BathType ==  SLAB_POTENTIAL, " PolymerVibrationalRelax_Initialize: slab potential not implemented" )
      CALL ERROR( BathType ==  DOUBLE_CHAIN, " PolymerVibrationalRelax_Initialize: double chain not implemented" )

      ! Define dimension of the system
      IF ( MorsePotential ) THEN
         NSystem = 1
         ALLOCATE( SystemMasses(NSystem) )
         SystemMasses = (/ MassH /)
      ELSE 
         NSystem = 4
         ALLOCATE( SystemMasses(NSystem) )
         SystemMasses = (/ (MassH, iCoord=1,3), MassC /)
      END IF

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses

      IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         NDim = (NSystem+NBath)
         ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads), APre(NDim*NBeads), MassVector(NDim*NBeads) )
         DO iBead = 1, NBeads
            MassVector((iBead-1)*NDim+1:iBead*NDim ) = (/ (SystemMasses(iCoord), iCoord=1,NSystem), (MassBath, iCoord=1,NBath) /)
         END DO

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! store the nr of dimensions
         NDim = NSystem
         ALLOCATE( X(NDim*NBeads), V(NDim*NBeads), A(NDim*NBeads), APre(NDim*NBeads), MassVector(NDim*NBeads) )
         DO iBead = 1, NBeads
            MassVector( (iBead-1)*NDim + 1 : iBead*NDim ) = SystemMasses(:)
         END DO
      END IF

      ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
      CALL EvolutionSetup( MolecularDynamics, NDim*NBeads, MassVector, TimeStep )

      ALLOCATE( LangevinSwitchOn( NDim*NBeads ) )

      ! Switch on Langevin relaxation when needed
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         DO iBead = 1, NBeads
            LangevinSwitchOn( (iBead-1)*NDim + 1 : iBead*NDim ) = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)   
         END DO
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration for the Bath in the canonical ensamble
      CALL EvolutionSetup( Equilibration, NDim*NBeads, MassVector, EquilTStep )
      LangevinSwitchOn( 1 : NBeads*NDim ) = .TRUE.
      IF ( .NOT. MorsePotential ) THEN
         DO iBead = 1, NBeads
            IF ( Collinear )  LangevinSwitchOn((iBead-1)*NDim+1:iBead*NDim) = (/ .FALSE. , .FALSE. , (.TRUE., iCoord=1, NDim-2) /)   
         END DO
      END IF
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature*NBeads, LangevinSwitchOn )

      DEALLOCATE( LangevinSwitchOn, SystemMasses )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the average energies over time only
      ! DEBUG - print also the average coordinates
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

         ALLOCATE(  EquilibAveTvsTime( NrEquilibrationSteps/PrintStepInterval + 1 ) )
         EquilibAveTvsTime(:) = 0.0
      END IF

      ! Initialize random number seed
      CALL SetSeed( 1 )


!       ALLOCATE( DissipatedPower(0:NrOfPrintSteps) )

   END SUBROUTINE PolymerVibrationalRelax_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Run()
      IMPLICIT NONE
      !> Output units: minimal output
      INTEGER  ::  AvEnergyOutputUnit
      !> Output units: standard output
      INTEGER  ::  AvCoordOutputUnit, AvBathCoordUnit, InitTUnit
      !> Output units: debug output
      INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel, TEquilUnit
      ! Standard deviations for initial coordinate and velocity
      REAL     ::  SigmaX, SigmaV
      !> Averages during equilibration
      REAL     ::  IstTemperature, TempAverage, TempVariance, AvailKin
      !> Energy averages
      REAL     ::  VSys, KSys, Ecoup, VBath, KBath
      REAL     ::  TotEnergy, PotEnergy, KinEnergy
      !> integer Counters
      INTEGER  ::  iTraj, iBead, iStep, kStep, iCoord
      !> Pointers to section of coordinate and velocity arrays
      REAL, DIMENSION(:), POINTER  :: SingleX, SingleV
      !> Filename for output files
      CHARACTER(100) :: OutFileName
      !> Check variable
      LOGICAL  :: ExitCheck
      !> Dummy array
      REAL, DIMENSION(4) :: Dummy
      !> Centroid of positions and velocities
      REAL, DIMENSION(NDim) :: CentroidX, CentroidV

      REAL :: KSys2, KBath2
      REAL, DIMENSION(0:NrOfPrintSteps) :: VirialAverage, VirialAverage2

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

      IF ( MorsePotential ) THEN
         SigmaX = sqrt( Temperature / MassVector(1) ) / HarmonicFreq
         SigmaV = sqrt( Temperature / MassVector(1) )
      END IF

      !run NrTrajs number of trajectories
      DO iTraj = 1,NrTrajs
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! set reasonable initial conditions for the equilibration
         DO iBead = 1, NBeads

            ! Point to the current bead of the polymer ring
            SingleX => X( (iBead-1)*NDim+1 : iBead*NDim )
            SingleV => V( (iBead-1)*NDim+1 : iBead*NDim )

            IF ( MorsePotential ) THEN
               SingleX(1) = GaussianRandomNr( 1.0 ) * SigmaX
               SingleV(1) = GaussianRandomNr( 1.0 ) * SigmaV
            ELSE
               ! Equilibrium geometry of C and H, no velocity
               SingleX(1:2) = 0.0
               SingleX(3) = HZEquilibrium 
               SingleX(4) = C1Puckering
               AvailKin = InitEnergy - MinimumEnergy
               SingleV(3) = + sqrt( 2.0 * AvailKin / MassH ) * 0.958234410548192
               SingleV(4) = - sqrt( 2.0 * AvailKin / MassC ) * 0.285983940880181 
            END IF

            ! Set bath initial conditions
            IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
               CALL ThermalEquilibriumBathConditions( Bath, SingleX(NSystem+1:NDim), SingleV(NSystem+1:NDim), NBeads*Temperature )
            END IF

         END DO

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at Tp = ", NBeads*Temperature / MyConsts_K2AU

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = VibrRelaxPotential( X, A )
         A(:) = A(:) / MassVector(:)

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         kStep = 0
         iStep = 0

         ! Long equilibration step before starting the dynamics
         DO
            iStep = iStep + 1

            ! PROPAGATION for ONE TIME STEP 
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, VibrRelaxPotential, PotEnergy )

            IF ( PrintType >= FULL ) THEN
               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*NDim*NBeads)

               ! store temperature average and variance
               TempAverage = TempAverage + IstTemperature
               TempVariance = TempVariance + IstTemperature**2
            END IF

            ! every PrintStepInterval steps, write output
            IF ( mod(iStep-1,PrintStepInterval) == 0 ) THEN
               kStep = kStep + 1                                   ! Increment counter of print interval
                
               IF ( PrintType >= FULL .AND. kStep <= size(EquilibAveTvsTime) ) THEN     ! compute average temperature
                  EquilibAveTvsTime(kStep) = EquilibAveTvsTime(kStep) + TempAverage/iStep
               END IF
               
               IF ( PrintType == DEBUG ) THEN              ! write debug output
                  ! Temperature profile during equilibration
                  IF ( iStep == 1 )  WRITE(TEquilUnit,*) " " 
                  WRITE(TEquilUnit,851)  real(iStep)*EquilTStep/MyConsts_fs2AU, &
                                                      TempAverage/iStep, sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
               END IF
               851 FORMAT( F20.5, 2F20.6 )
            END IF

            ! Check if the potential for all the beads of the system is lower than the available energy
            IF ( iStep >= NrEquilibrationSteps ) THEN
               ExitCheck = .TRUE.
               DO iBead = 1, NBeads
                  IF ( MorsePotential ) THEN
                     VSys = MorseV( X((iBead-1)*NDim+1), Dummy(1) )
                  ELSE 
                     VSys = VHFourDimensional( X((iBead-1)*NDim+1:(iBead-1)*NDim+4), Dummy(1:4) )
                  END IF
                  ExitCheck = ExitCheck .AND. ( VSys <= InitEnergy ) 
               END DO
               IF ( ExitCheck ) EXIT
            END IF

         END DO

         IF ( Collinear .AND. ( .NOT. MorsePotential ) ) THEN
            DO iBead = 1, NBeads
               X((iBead-1)*NDim+1:(iBead-1)*NDim+2) = 0.0
            END DO
         END IF

         ! change system velocities to fix the initial energy of the system
         DO iBead = 1, NBeads
            IF ( MorsePotential ) THEN
               VSys = MorseV( X((iBead-1)*NDim+1), Dummy(1) )
            ELSE 
               VSys = VHFourDimensional( X((iBead-1)*NDim+1:(iBead-1)*NDim+4), Dummy(1:4) )
            END IF
            AvailKin = InitEnergy-VSys
            IF ( MorsePotential ) THEN
               V((iBead-1)*NDim+1) = sqrt( 2.0 * AvailKin / MassH )
            ELSE
               V((iBead-1)*NDim+1) = 0.0
               V((iBead-1)*NDim+2) = 0.0
               V((iBead-1)*NDim+3) = + sqrt( 2.0 * AvailKin / MassH ) * 0.958234410548192
               V((iBead-1)*NDim+4) = - sqrt( 2.0 * AvailKin / MassC ) * 0.285983940880181 
            END IF
         END DO

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Energy of the system, Coupling energy and energy of the bath
         CALL IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
         CALL Virial( KSys2, KBath2 )

!          PRINT*, KSys*MyConsts_Hartree2eV, VSys*MyConsts_Hartree2eV, Ecoup*MyConsts_Hartree2eV, &
!                       VBath*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV
!          STOP

         IF ( NDim > NSystem ) THEN
            IstTemperature = 2.0*KBath/(MyConsts_K2AU*(NDim-NSystem))
         ELSE
            IstTemperature = 0.0
         ENDIF

         ! Store starting values of the averages
         AverageESys(0)         = AverageESys(0)  + KSys + VSys
         AverageEBath(0)        = AverageEBath(0) + KBath + VBath
         AverageECoup(0)        = AverageECoup(0) + Ecoup
         CentroidX = CentroidCoord( X )
         CentroidV = CentroidCoord( V )
         IF ( PrintType >= FULL )  AverageCoord(1:NDim,0) = AverageCoord(1:NDim,0) + CentroidX(:)

!          DissipatedPower(:) = DissipatedPower(:) - DissipativeTermIntegral( X )

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

            ! Write initial values
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | ESys, ECoup, EBath, KSys, VSys, KBath, VBath / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, (/ (CentroidX(iCoord), iCoord = 1,NSystem) /), &
                                                            (/ (CentroidX(iCoord+NSystem)+0.05*iCoord, iCoord = 1, NBath) /)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, CentroidV(:)

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
         A(:) = A(:) / MassVector(:)

         VirialAverage(:) = 0.0
         VirialAverage2(:) = 0.0

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VibrRelaxPotential, PotEnergy )

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 

!                ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
!                IF ( BathType == SLAB_POTENTIAL ) THEN 
!                   X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
!                ENDIF

               ! Energy of the system
               CALL IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
               CALL Virial( KSys2, KBath2 )


               ! Update averages over time
               AverageESys(kStep)            = AverageESys(kStep)  + KSys !+ VSys
               AverageEBath(kStep)           = AverageEBath(kStep) + KBath !+ VBath
               AverageECoup(kStep)           = AverageECoup(kStep) + Ecoup

!                WRITE(777,"(1000F20.5)") kStep*1.0, KSys, KSys2
               VirialAverage(kStep) = VirialAverage(kStep) + KSys2 !+ VSys
               VirialAverage2(kStep) = VirialAverage2(kStep) + KBath2 !+ VSys

               CentroidX = CentroidCoord( X )
               CentroidV = CentroidCoord( V )
               IF ( PrintType >= FULL )   AverageCoord(1:NDim,kStep) = AverageCoord(1:NDim,kStep) + CentroidX(:)


!                DissipatedPower(kStep)  = DissipatedPower(kStep) + DissipativeTermIntegral( X )
!                DissipatedPower(kStep:NrOfPrintSteps)  = DissipatedPower(kStep:NrOfPrintSteps) + DissipativeTermDifferential( X, V )

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                        KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, (/ (CentroidX(iCoord), iCoord = 1,NSystem) /), &
                                          (/ (CentroidX(iCoord+NSystem)+0.05*iCoord, iCoord = 1, NBath) /)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, CentroidV(:)
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
      IF ( PrintType >= FULL )   AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL )   EquilibAveTvsTime(:)   = EquilibAveTvsTime(:) / real(NrTrajs) 
      VirialAverage = VirialAverage  / real(NrTrajs)
      VirialAverage2 = VirialAverage2  / real(NrTrajs)

!       DissipatedPower(:) =  DissipatedPower(:) / real(NrTrajs)
!       DO iStep = 0, NrOfPrintSteps
!          WRITE(999,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, DissipatedPower(iStep )
!       END DO

      DO iStep = 1, NrOfPrintSteps
         WRITE(999,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                 VirialAverage(iStep)*MyConsts_Hartree2eV,   AverageESys(iStep )*MyConsts_Hartree2eV
      END DO
      DO iStep = 1, NrOfPrintSteps
         WRITE(998,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                 VirialAverage2(iStep)*MyConsts_Hartree2eV,   AverageEBath(iStep )*MyConsts_Hartree2eV
      END DO

      ! Print average temperature over time during equilibration
      DO iStep = 0, size(EquilibAveTvsTime)-1
          WRITE(InitTUnit,"(2F20.8)") EquilTStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, EquilibAveTvsTime(iStep+1)
      END DO

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(AvEnergyOutputUnit,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
         AverageESys(iStep)*MyConsts_Hartree2eV, AverageECoup(iStep)*MyConsts_Hartree2eV,  AverageEBath(iStep)*MyConsts_Hartree2eV
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

      ENDIF

      ! Close output files
      CLOSE( AvEnergyOutputUnit )
      IF ( PrintType >= FULL ) THEN
         CLOSE( AvCoordOutputUnit )
         CLOSE( AvBathCoordUnit )
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

   END SUBROUTINE PolymerVibrationalRelax_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE PolymerVibrationalRelax_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( X, V, A, APre, MassVector )
      DEALLOCATE( AverageESys, AverageEBath, AverageECoup )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE PolymerVibrationalRelax_Dispose

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

!    REAL FUNCTION DissipativeTermIntegral( X )  RESULT( DissipIntegral )
!       IMPLICIT NONE
!       REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X
!       INTEGER :: i
! 
!       DissipIntegral = 0.0
!       DO i = 1, NBeads
!          DissipIntegral = DissipIntegral + GetDissipatedPower( Bath, X((i-1)*NDim+NSystem+1 : i*NDim), X((i-1)*NDim+NSystem)  ) 
!       END DO
!       DissipIntegral = DissipIntegral / NBeads
! 
!    END FUNCTION DissipativeTermIntegral
! 
!    REAL FUNCTION DissipativeTermDifferential( X, V )  RESULT( DissipDiffer )
!       IMPLICIT NONE
!       REAL, DIMENSION( NDim*NBeads ), INTENT(IN) :: X, V
!       REAL, DIMENSION( NDim )   :: CentroidV
!       INTEGER :: i, iCoord
! 
!       DissipDiffer = 0.0
!       CentroidV = CentroidCoord( V )
!       DO i = 1, NBeads
!          DissipDiffer = DissipDiffer + GetDissipatedPower( Bath, CentroidV(NSystem+1:NDim), X((i-1)*NDim+NSystem)  ) 
!       END DO
!       DissipDiffer = DissipDiffer / NBeads

!    END FUNCTION DissipativeTermDifferential

!*************************************************************************************************



   REAL FUNCTION VibrRelaxPotential( Positions, Forces )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces

      INTEGER :: iBead, iCoord
      REAL, DIMENSION(NBeads)  :: RingCoord, RingForces 
      REAL, DIMENSION(NDim)    :: SingleF
      REAL :: SingleV

      ! Check the number degrees of freedom
      CALL ERROR( size(Forces) /= size( Positions ), "PolymerVibrationalRelax.VibrRelaxPotential: array dimension mismatch" )
      CALL ERROR( size( Positions ) /= NDim*NBeads, "PolymerVibrationalRelax.VibrRelaxPotential: wrong number of DoFs" )

      ! Initialize forces and potential
      VibrRelaxPotential = 0.0
      Forces(:)          = 0.0

      ! System-bath potential for each bead
      DO iBead = 1, NBeads 
         ! Compute potential and force for a single bead
         CALL  SingleBeadPotential( Positions((iBead-1)*NDim+1:iBead*NDim), SingleF, SingleV )
         ! Sum up to total potential and total force
         VibrRelaxPotential = VibrRelaxPotential + SingleV
         Forces((iBead-1)*NDim+1:iBead*NDim) = Forces((iBead-1)*NDim+1:iBead*NDim) + SingleF(:)
      END DO

      IF ( NBeads > 1 ) THEN
         ! Harmonic forces between beads
         DO iCoord = 1, NDim
            DO iBead = 1, NBeads
               RingCoord(iBead) = X( iCoord + (iBead-1)*NDim )
            END DO
            CALL PolymerRingPotential( MassVector(iCoord), RingCoord, VibrRelaxPotential, RingForces )
            DO iBead = 1, NBeads
               Forces( iCoord + (iBead-1)*NDim ) =  Forces( iCoord + (iBead-1)*NDim ) + RingForces(iBead)
            END DO
         END DO
      END IF

      ! when collinear calculation, make sure that the forces on HX and HY are zero
      IF ( Collinear .AND. (.NOT. MorsePotential )) THEN
         DO iBead = 1, NBeads 
            Forces((iBead-1)*NDim+1 : (iBead-1)*NDim+2) = 0.0
         END DO
      END IF

   END FUNCTION VibrRelaxPotential

   REAL FUNCTION MorseV( Positions, Forces ) RESULT(V) 
      IMPLICIT NONE
      REAL, DIMENSION(1), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(1), TARGET, INTENT(OUT) :: Forces 

      V = MorseDe * ( exp(-2.0*MorseAlpha*Positions(1)) - 2.0 * exp(-MorseAlpha*Positions(1)) )  
      Forces(1) = 2.0 * MorseAlpha * MorseDe * (  exp(-2.0*MorseAlpha*Positions(1)) - exp(-MorseAlpha*Positions(1)) )  

   END FUNCTION MorseV

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

   SUBROUTINE SingleBeadPotential( Positions, Forces, V )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      REAL, INTENT(OUT)                       :: V

      INTEGER                  :: iCoord
      REAL, DIMENSION(NBath)   :: ForceQCoord
      REAL, DIMENSION(NSystem) :: ForceSystem
      REAL                     :: ForceCoupCoord, CouplingFunc

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
      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ForceCoupCoord = 0.0
         ForceQCoord(:) = 0.0
         IF ( MorsePotential ) THEN
            CouplingFunc = Positions(1)
         ELSE 
            CouplingFunc = Positions(4) - C1Puckering
         END IF
         CALL BathPotentialAndForces( Bath, CouplingFunc, Positions( NSystem+1:NDim ), V, ForceCoupCoord, ForceQCoord(:) ) 
         IF ( MorsePotential ) THEN
            Forces(1) = Forces(1) + ForceCoupCoord
         ELSE
            Forces(4) = Forces(4) + ForceCoupCoord
         END IF
         Forces( NSystem+1:NDim ) = Forces( NSystem+1:NDim ) + ForceQCoord(:)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! DO NOTHING
      END IF

   END SUBROUTINE SingleBeadPotential


!*************************************************************************************************
   
   SUBROUTINE Virial( KSys, KBath )
      IMPLICIT NONE
      REAL, INTENT(OUT)      :: KSys, KBath
      REAL, DIMENSION(NDim)  :: SngForce
      INTEGER :: iBead, iCoord
      REAL :: V

      REAL, DIMENSION(NBath)   :: ForceQCoord
      REAL                     :: ForceCoupCoord, CouplingFunc

      KSys = 0.0
      KBath = 0.0

      ! System-bath potential for each bead
      DO iBead = 1, NBeads 
         ! Compute potential and force for a single bead
         CALL  SingleBeadPotential( X((iBead-1)*NDim+1:iBead*NDim), SngForce, V )

!          ForceCoupCoord = 0.0
!          ForceQCoord(:) = 0.0
!          IF ( MorsePotential ) THEN
!             CouplingFunc = X((iBead-1)*NDim+1) 
!          ELSE 
!             CouplingFunc = X((iBead-1)*NDim+4) - C1Puckering
!          END IF
!       CALL BathPotentialAndForces( Bath, CouplingFunc, X((iBead-1)*NDim+NSystem+1:iBead*NDim), V, ForceCoupCoord, ForceQCoord(:)) 
!             KSys  = KSys  -0.5 *  CouplingFunc * ForceCoupCoord

         ! Compute virial for the current bead, divided between the system and the bath
         KSys  = KSys  -0.5 * ( X((iBead-1)*NDim+3) - HZEquilibrium ) * SngForce(3)
         KSys  = KSys  -0.5 * ( X((iBead-1)*NDim+4) - C1Puckering ) * SngForce(4)
         DO iCoord = NSystem+1, NDim
            KBath  = KBath  -0.5 *  X((iBead-1)*NDim+iCoord) * SngForce(iCoord)
         END DO
      END DO

      ! Normalize for nr of beads
      KSys = KSys / NBeads
      KBath = KBath / NBeads

   END SUBROUTINE Virial

   SUBROUTINE IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
      IMPLICIT NONE
      REAL, INTENT(OUT) :: KSys, VSys, Ecoup, VBath, KBath
      REAL, DIMENSION(NDim) :: CentroidV, Dummy
      REAL    :: VBath1, Ecoup1, CouplingFunc
      INTEGER :: iCoord, iBead

      ! Centroid of the velocities
      CentroidV = CentroidCoord( V )

!       ! Compute kinetic energy as average of the kinetic energy over the replicas
!       KSys = 0.0
!       DO iBead = 1, NBeads
!          DO iCoord = 1, 4
!             KSys = KSys + 0.5 * MassVector(iCoord) * V( (iBead-1)*NDim+iCoord )**2
!          END DO
!       END DO
!       KSys = KSys / NBeads 

      ! Compute kinetic energy as the kinetic energy of the centroid
      KSys = 0.0
      DO iCoord = 1, NSystem
         KSys = KSys + 0.5 * MassVector(iCoord) * CentroidV(iCoord)**2
      END DO

      ! Energy of the system - average of the potential energies of the system replicas
      VSys = 0.0
      DO iBead = 1, NBeads
         IF ( MorsePotential ) THEN
            VSys = VSys + MorseV( X((iBead-1)*NDim+1), Dummy(1) )
         ELSE 
            VSys = VSys + VHFourDimensional( X((iBead-1)*NDim+1:(iBead-1)*NDim+4), Dummy(1:4) )
         END IF
      END DO
      VSys = VSys / NBeads 
      
      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN

         ! Kinetic energy - kinetic energy of the velocity centroid for the bath
         KBath = 0.0
         DO iCoord = 5, NDim
            KBath = KBath + 0.5 * MassVector(iCoord) * CentroidV(iCoord)**2
         END DO

!          ! Compute kinetic energy as average of the kinetic energy over the replicas
!          KBath = 0.0
!          DO iBead = 1, NBeads
!             DO iCoord = 5, NDim
!                KBath = KBath + 0.5 * MassVector(iCoord) * V( (iBead-1)*NDim+iCoord )**2
!             END DO
!          END DO
!          KBath = KBath / NBeads 

         ! Coupling energy and energy of the bath - average of the potential energies of the system replicas
         Ecoup = 0.0
         VBath = 0.0
         DO iBead = 1, NBeads
            IF ( MorsePotential ) THEN
               CouplingFunc = X((iBead-1)*NDim+1)
            ELSE 
               CouplingFunc = X((iBead-1)*NDim+4) - C1Puckering
            END IF
            CALL EnergyOfTheBath( Bath, CouplingFunc, X( (iBead-1)*NDim+NSystem+1 : iBead*NDim ), Ecoup1, VBath1 )
            Ecoup = Ecoup + Ecoup1
            VBath = VBath + VBath1
         END DO
         Ecoup = Ecoup / NBeads
         VBath = VBath / NBeads

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN

         ! The energies of the bath and the coupling cannot be defined
         Ecoup = 0.0
         VBath = 0.0
         KBath = 0.0

      END IF

   END SUBROUTINE IstantaneousEnergies

!*************************************************************************************************

END MODULE PolymerVibrationalRelax