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
!>
!>  \todo         Implement energy flux along the chain
!>  \todo         The cluster bath has some problems... needs debugging
!>  |todo         Initial conditions for non collinear calculations should
!>                be defined
!>                 
!***************************************************************************************
MODULE VibrationalRelax
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

   PUBLIC :: VibrationalRelax_ReadInput, VibrationalRelax_Initialize, VibrationalRelax_Run, VibrationalRelax_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

   ! Initial conditions
   REAL    :: InitEnergy                !< Initial energy of the system
   REAL    :: TimeBetweenSnaps          !< Nr of initial snapshots to randomize initial conditions
   INTEGER :: NrOfInitSnapshots         !< Time between each initial snapshot

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of time steps between each printing of propagation info ( = NrOfSteps / NrOfPrintSteps )

   ! Time evolution dataset
   TYPE(Evolution) :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results
   TYPE(Evolution) :: InitialConditions     !< Propagate in microcanonical ensamble to generate initial conditions

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageEBath       !< Average energy of the bath vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageECoup       !< Average coupling energy vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageESys        !< Average energy of the system vs time
   REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord       !< Average i-th coordinate vs time
   COMPLEX, DIMENSION(:), ALLOCATABLE   :: PowerSpec          !< Autocorrelation function and its power spectrum
   REAL, DIMENSION(:), ALLOCATABLE      :: XatT0              !< Initial position of the dynamics
   REAL, DIMENSION(:,:), ALLOCATABLE    :: OscillCorrelations !< 1-n oscillator correlations in the chain

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

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Initial energy of the system (input in eV and transform to AU)
      CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
      InitEnergy = InitEnergy / MyConsts_Hartree2eV

      ! Nr of initial snapshots
      CALL SetFieldFromInput( InputData, "NrOfInitSnapshots", NrOfInitSnapshots )

      ! Time between each snapshots (input in fs and transform to AU)
      CALL SetFieldFromInput( InputData, "TimeBetweenSnaps", TimeBetweenSnaps)
      TimeBetweenSnaps = TimeBetweenSnaps * MyConsts_fs2AU

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

      WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV, &
                    NrOfInitSnapshots, TimeBetweenSnaps/MyConsts_fs2AU

      WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps

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

   END SUBROUTINE VibrationalRelax_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord, NDim
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         ALLOCATE( X(3+NCarbon), V(3+NCarbon), A(3+NCarbon), APre(3+NCarbon), MassVector(3+NCarbon), LangevinSwitchOn(3+NCarbon) )
         MassVector = (/ (MassH, iCoord=1,3), (MassC, iCoord=1,NCarbon) /)

      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ALLOCATE( X(4+NBath), V(4+NBath), A(4+NBath), APre(4+NBath), MassVector(4+NBath), LangevinSwitchOn(4+NBath) )
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,NBath) /)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ALLOCATE( X(4), V(4), A(4), APre(4), MassVector(4), LangevinSwitchOn(4) )
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

      ! Set canonical dynamics at the end of the oscillator chain
      IF ( BathType == CHAIN_BATH .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( NDim ) = .TRUE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Switch on Langevin relaxation when needed
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration of the 4D system only in the microcanonical ensamble 
      CALL EvolutionSetup( InitialConditions, 4, MassVector(1:4), TimeStep )

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
         ALLOCATE( AverageCoord(NDim,0:NrOfPrintSteps), PowerSpec(0:NrOfPrintSteps), XatT0(NDim) )
         AverageCoord(1:NDim,0:NrOfPrintSteps) = 0.0
         PowerSpec(:) = cmplx( 0.0, 0.0 )
      END IF

      ! Correlations in the chain bath
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH ) THEN
         ALLOCATE( OscillCorrelations(NDim-4, 0:NrOfPrintSteps) )
         OscillCorrelations(:,:) = 0.0
      END IF

      ! Initialize random number seed
      CALL SetSeed( 1 )

   END SUBROUTINE VibrationalRelax_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Run()
      IMPLICIT NONE
      INTEGER  ::  AvEnergyOutputUnit, AvCoordOutputUnit       ! UNITs FOR OUTPUT AND DEBUG
      INTEGER  ::  AvBathCoordUnit, OscillCorrUnit, InitialCondUnit, PowerSpectrumUnit
      INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel
      REAL     ::  VSys, KSys, Ecoup, VBath, KBath               ! ISTANTANEOUS ENERGY VALUES
      REAL     ::  TotEnergy, PotEnergy, KinEnergy, dOmega, IstTemperature
      INTEGER  ::  NTimeStepEachSnap, iCoord, iTraj, iStep,  iOmega, kStep, NInit
      CHARACTER(100) :: OutFileName
      REAL, DIMENSION(NrOfInitSnapshots,8) :: CHInitConditions
      REAL, DIMENSION(:), ALLOCATABLE      :: AverCoordOverTime

      AvEnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageEnergy.dat", UNIT=AvEnergyOutputUnit )
      WRITE(AvEnergyOutputUnit, "(A,I6,A,/)") "# average energy vs time (fs | eV) - ", NrTrajs, " trajectories "

      AvCoordOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
      WRITE(AvCoordOutputUnit, "(A,I6,A,/)") "# average coordinate vs time (fs | Ang) - ", NrTrajs, " trajectories "

      IF ( PrintType >= FULL ) THEN
         AvBathCoordUnit = LookForFreeUnit()
         OPEN( FILE="AverageBathCoords.dat", UNIT=AvBathCoordUnit )
         WRITE(AvBathCoordUnit, "(A,I6,A,/)") "# average Q coordinate vs time (Ang | fs) - ", NrTrajs, " trajectories "

         PowerSpectrumUnit = LookForFreeUnit()
         OPEN( FILE="PowerSpectrum.dat", UNIT=PowerSpectrumUnit )
         WRITE(PowerSpectrumUnit, "(A,I6,A,/)") "# Power spectrum of the trajs - ", NrTrajs, " trajectories (fs | au)"
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
      ENDIF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      ! Define the set of initial conditions for CH
      ! atoms in the equilibrium geometry, energy in the stretching normal mode

      X(1:2) = 0.0
      V(1:2) = 0.0
      X(3) = HZEquilibrium 
      X(4) = C1Puckering  
      V(3) = + sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassH ) * 0.958234410548192
      V(4) = - sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / MassC ) * 0.285983940880181 

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
            CALL EOM_VelocityVerlet( InitialConditions, X(1:4), V(1:4), A(1:4), VHFourDimensional, PotEnergy )
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
               CALL ZeroKelvinSlabConditions( X, V, CHInitConditions ) 
               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
         ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
               CALL ZeroKelvinBathConditions( X, V, CHInitConditions )
         ELSE IF ( BathType == LANGEVIN_DYN ) THEN
               NInit = CEILING( UniformRandomNr(0.0, real(NrOfInitSnapshots) ) )
               X(1:4) = CHInitConditions( NInit, 1:4 )
               V(1:4) = CHInitConditions( NInit, 5:8 )
         END IF

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Energy of the system, Coupling energy and energy of the bath
         CALL IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
         IF ( NDim > 4 ) THEN
            IstTemperature = 2.0*KBath/(MyConsts_K2AU*(NDim-4))
         ELSE
            IstTemperature = 0.0
         ENDIF

         ! Store starting values of the averages
         AverageESys(0)         = AverageESys(0)  + KSys + VSys
         AverageEBath(0)        = AverageEBath(0) + KBath + VBath
         AverageECoup(0)        = AverageECoup(0) + Ecoup
         AverageCoord(1:NDim,0) = AverageCoord(1:NDim,0) + X(:)

         IF ( PrintType >= FULL ) THEN
            XatT0(:) = X(:) 
            PowerSpec(0) = PowerSpec(0) + TheOneWithVectorDotVector( XatT0(1:NDim), XatT0(1:NDim) )
         ENDIF
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

            ! Write initial values
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | ESys, ECoup, EBath, KSys, VSys, KBath, VBath / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, X(1), X(2), X(3), X(4), (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)

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

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,NrOfSteps

            ! Propagate for one timestep
            CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VibrRelaxPotential, PotEnergy )

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 

               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               IF ( BathType == SLAB_POTENTIAL ) THEN 
                  X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
               ENDIF

               ! Energy of the system
               CALL IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
               
               ! Update averages over time
               AverageESys(kStep)            = AverageESys(kStep)  + KSys + VSys
               AverageEBath(kStep)           = AverageEBath(kStep) + KBath + VBath
               AverageECoup(kStep)           = AverageECoup(kStep) + Ecoup

               ! store the autocorrelation function for power spectrum computation
               IF ( PrintType >= FULL )  THEN
                  AverageCoord(1:NDim,kStep) = AverageCoord(1:NDim,kStep) + X(:)
                  PowerSpec(kStep) = PowerSpec(kStep) + TheOneWithVectorDotVector( XatT0(1:NDim), X(1:NDim) )
               END IF
               IF ( BathType == CHAIN_BATH .AND. PrintType >= FULL ) THEN
                  DO iCoord = 1, NDim-4
                     OscillCorrelations(iCoord, kStep) = OscillCorrelations(iCoord, kStep) + X(4) * X(4+iCoord) 
                  END DO
               END IF

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                        KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                          (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
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
      IF ( PrintType >= FULL )    AverageCoord(1:NDim,:) = AverageCoord(1:NDim,:) / real(NrTrajs) 
      IF ( PrintType >= FULL .AND. BathType == CHAIN_BATH )    &
               OscillCorrelations(1:NDim-4,:) = OscillCorrelations(1:NDim-4,:)  / real(NrTrajs)

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, NrOfPrintSteps
         WRITE(AvEnergyOutputUnit,"(F14.8,3F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
          AverageESys(iStep)*MyConsts_Hartree2eV, AverageECoup(iStep)*MyConsts_Hartree2eV,  AverageEBath(iStep)*MyConsts_Hartree2eV
      END DO

      ! PRINT average coordinates
      DO iStep = 0, NrOfPrintSteps
         WRITE(AvCoordOutputUnit,"(F14.8,4F14.8)") TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU, &
                                    AverageCoord(1:4,iStep)*MyConsts_Bohr2Ang 
      END DO

      IF ( PrintType >= FULL ) THEN

         ! PRINT average bath coordinates
         DO iCoord = 1, NDim-4
            WRITE(AvBathCoordUnit,"(/,A,I5,/)") "#  Bath Coord # ", iCoord
            DO iStep = 0, NrOfPrintSteps
               WRITE(AvBathCoordUnit,"(F20.8,3F20.8)") real(iCoord)+AverageCoord(4+iCoord,iStep)*10, &
                        TimeStep*real(PrintStepInterval*iStep)/MyConsts_fs2AU
            END DO
         END DO

         ! normalize autocorrelation function, compute fourier transfrom and print
         ALLOCATE( AverCoordOverTime( NDim ) )
         AverCoordOverTime(:) = SUM ( AverageCoord, 2 ) / real( NrOfPrintSteps+1 )
         PowerSpec(:) = ( PowerSpec(:) /  real( NrTrajs ) ) - DOT_PRODUCT(AverCoordOverTime, AverCoordOverTime )
         dOmega =  2.*MyConsts_PI/( real(PrintStepInterval*NrOfPrintSteps) * TimeStep )
         CALL DiscreteFourier( PowerSpec )
         DO iOmega = 0, NrOfPrintSteps
            WRITE( PowerSpectrumUnit,"(F20.8,3F20.8)" )  iOmega*dOmega/MyConsts_cmmin1toAU, abs(PowerSpec(iOmega)), & 
                                                  real(PowerSpec(iOmega)), aimag( PowerSpec(iOmega))
         END DO
         IF ( BathType == CHAIN_BATH ) THEN
            DO iCoord = 1, NDim-4
               OscillCorrelations(iCoord,:) = OscillCorrelations(iCoord,:) - &
                                    AverCoordOverTime(3+1)*AverCoordOverTime(3+iCoord+1) 
               WRITE(OscillCorrUnit, *) iCoord, SUM( OscillCorrelations(iCoord, :) )/real(NrOfPrintSteps+1)
            END DO
         END IF
         DEALLOCATE( AverCoordOverTime )

      ENDIF

      ! Close output files
      CLOSE( AvEnergyOutputUnit )
      CLOSE( AvCoordOutputUnit )
      IF ( PrintType >= FULL ) THEN
         CLOSE( AvBathCoordUnit )
         CLOSE( PowerSpectrumUnit )
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

   END SUBROUTINE VibrationalRelax_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE VibrationalRelax_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( X, V, A, APre, MassVector )
      DEALLOCATE( AverageESys, AverageEBath, AverageECoup )
      IF ( PrintType >= FULL ) DEALLOCATE( AverageCoord, PowerSpec, XatT0 )
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

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         VibrRelaxPotential = VHSticking( Positions, Forces )

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces with independentoscillatormodel subroutine
         VibrRelaxPotential = GenericSystemAndBath( Positions, Forces, VHFourDimensional, 4 ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential using only the 4D subroutine
         VibrRelaxPotential = VHFourDimensional( Positions, Forces )

      END IF

   END FUNCTION VibrRelaxPotential

!*************************************************************************************************

   SUBROUTINE IstantaneousEnergies( KSys, VSys, Ecoup, VBath, KBath )
      IMPLICIT NONE
      REAL, INTENT(OUT) :: KSys, VSys, Ecoup, VBath, KBath
      REAL, DIMENSION(4) :: Dummy
      REAL :: KinEnergy

      ! Total kinetic energy
      KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )

      ! Energy of the system
      KSys = EOM_KineticEnergy( InitialConditions, V(1:4) )
      VSys = VHFourDimensional( X(1:4), Dummy )

      ! Coupling energy and energy of the bath
      IF ( BathType == SLAB_POTENTIAL ) THEN 
         Ecoup = 0.0
         VBath = PotEnergy - VSys
         KBath = KinEnergy - KSys
      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         Ecoup = CouplingEnergy( X )
         VBath = PotEnergyOfTheBath( X )
         KBath = KinEnergy - KSys
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         Ecoup = 0.0
         VBath = 0.0
         KBath = 0.0
      END IF

   END SUBROUTINE IstantaneousEnergies

!*************************************************************************************************

   SUBROUTINE DiscreteFourier( Vector )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:), INTENT(INOUT) :: Vector
      COMPLEX, DIMENSION(size(Vector)) :: Transform
      INTEGER :: N, i, j
      COMPLEX :: Factor

      N = size( Vector )

      Transform(:) = 0.0

      DO i = 0, N-1
         Factor = 2.0 * MyConsts_I * MyConsts_PI * real(i) / real(N)
         DO j = 0, N-1
            Transform(i+1) = Transform(i+1) + Vector(j+1) * exp( Factor * real(j) )
         END DO
      END DO

      Vector(:) = Transform(:)

   END SUBROUTINE DiscreteFourier

   SUBROUTINE DiscreteInverseFourier( Vector )
      IMPLICIT NONE
      COMPLEX, DIMENSION(:), INTENT(INOUT) :: Vector
      COMPLEX, DIMENSION(size(Vector)) :: Transform
      INTEGER :: N, i, j
      COMPLEX :: Factor

      N = size( Vector )

      Transform(:) = 0.0

      DO i = 0, N-1
         Factor = - MyConsts_I * 2.0 * MyConsts_PI * real(i) / real(N)
         DO j = 0, N-1
            Transform(i+1) = Transform(i+1) + Vector(j+1) * exp( Factor * real(j) )
         END DO
      END DO

      Vector(:) = Transform(:) / ( N-1 )

   END SUBROUTINE DiscreteInverseFourier

!*************************************************************************************************

END MODULE VibrationalRelax