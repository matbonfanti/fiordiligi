!***************************************************************************************
!*                              MODULE Scattering
!***************************************************************************************
!
!>  \brief     Subroutines for a scattering simulation of a atom surface collision
!>  \details   This module contains subroutines to compute averages  \n
!>             in an scattering simulation of a atom-surface system  \n
!>             at a given temperature of the surface and a given incident  \n
!>             energy of the projectile atom.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             10 October 2013
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg 
!                
!>  \todo             print XYZ file of scattering trajectory
!
!***************************************************************************************
MODULE ScatteringSimulation
#include "preprocessoptions.cpp"
!    USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Scattering_ReadInput, Scattering_Initialize, Scattering_Run, Scattering_Dispose

   ! Initial conditions of the H atom
   REAL    :: EKinZH                    !< Initial kinetic energy of the scattering H atom
   REAL    :: ZHInit                    !< Initial Z position of the scattering H atom
   INTEGER :: NRhoMax                   !< Nr of rho values to sample
   REAL    :: DeltaRho                  !< Grid spacing in rho
  
   ! Initial equilibration of the slab
   REAL    :: Temperature               !< Temperature of the simulation
   REAL    :: EquilTStep                !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma                !< Friction parameter of the Langevin equation
   INTEGER :: NrEquilibSteps            !< Nr of time step of the equilibration
   INTEGER :: EquilibrationStepInterval !< Nr of tsteps between each print of equil info (=NrEquilibSteps/NrOfPrintSteps )

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories per impact parameter value
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of tsteps between each print of propag info (=NrOfSteps/NrOfPrintSteps)

   ! Time evolution dataset
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in macrocanonical ensamble at given T to generate init conditions

   ! Averages computed during propagation
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageZHydro      !< Average position of zH vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageVHydro      !< Average velocity of zH vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageZCarbon     !< Average position of zC vs time
   REAL, DIMENSION(:), ALLOCATABLE      :: AverageVCarbon     !< Average velocity of zC vs time
   REAL, DIMENSION(:,:), ALLOCATABLE    :: TrappingProb       !< trapping prob vs time and impact parameter

   ! Random number generator internal status
   TYPE(RNGInternalState), SAVE :: RandomNr    !< spectial type dat to store internal status of random nr generator

   INTEGER :: NDim

   CONTAINS


!*******************************************************************************
!> Read from input unit the variable which are specific for a scattering
!> simulation at given T and Einc.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE Scattering_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! PARAMETERS FOR THE INITIAL CONDITIONS OF H 

      ! Initial kinetic energy of the scattering H atom
      CALL SetFieldFromInput( InputData, "EKinZH", EKinZH )
      EKinZH = EKinZH * EnergyConversion(InputUnits, InternalUnits)
      ! Initial Z position of the scattering H atoms
      CALL SetFieldFromInput( InputData, "ZHInit", ZHInit )
      ZHInit = ZHInit * LengthConversion(InputUnits, InternalUnits)

      ! Nr of rho values to sample
      CALL SetFieldFromInput( InputData, "NRhoMax", NRhoMax )
      ! Grid spacing in rho
      CALL SetFieldFromInput( InputData, "DeltaRho", DeltaRho )
      DeltaRho = DeltaRho * LengthConversion(InputUnits, InternalUnits)

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

      WRITE(*, 902) ZHInit*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                    EKinZH*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                    NRhoMax, DeltaRho*NRhoMax*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)

      WRITE(*, 903) NrTrajs, TimeStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrOfSteps, NrOfPrintSteps

      IF ( BathType ==  SLAB_POTENTIAL ) &
         WRITE(*, 904) 1./EquilGamma*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    EquilTStep*NrEquilibSteps*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrEquilibSteps

   902 FORMAT(" * Initial conditions of the atom-surface system ", /,&
              " * Initial height of H atom:                    ",F10.4,1X,A,/,&  
              " * Initial energy of H atom:                    ",F10.4,1X,A,/,&  
              " * Nr of impact parameter values:               ",I10,       /,& 
              " * Max value of impact parameter:               ",F10.4,1X,A,/ )

   903 FORMAT(" * Dynamical simulation variables               ",           /,&
              " * Nr of trajectories per impact parameter:     ",I10,       /,& 
              " * Propagation time step:                       ",F10.4,1X,A,/,& 
              " * Nr of time steps of each trajectory:         ",I10,       /,& 
              " * Nr of print steps of each trajectory:        ",I10,       / )

   904 FORMAT(" * Bath equilibration variables                 ",           /,&
              " * Relaxation time of the Langevin dynamics:    ",F10.4,1X,A,/,& 
              " * Equilibration time step:                     ",F10.4,1X,A,/,& 
              " * Equilibration total time:                    ",F10.4,1X,A,/,&
              " * Nr of equilibration steps:                   ",I10,       / )

   END SUBROUTINE Scattering_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the scattering simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Initialize()
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn
      INTEGER :: iCoord

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

      ! In case of Langevin relaxation, switch on gamma at the carbon atom
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, 0.0, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
         LangevinSwitchOn = (/ (.FALSE., iCoord=1,3), (.TRUE., iCoord=4,NDim ) /)
         CALL SetupThermostat( Equilibration, EquilGamma, Temperature, LangevinSwitchOn )
      END IF

      DEALLOCATE( LangevinSwitchOn )

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! Allocate and initialize the variables for the trajectory averages
      ALLOCATE( AverageVHydro(0:NrOfPrintSteps), AverageZHydro(0:NrOfPrintSteps) )
      ALLOCATE( AverageVCarbon(0:NrOfPrintSteps), AverageZCarbon(0:NrOfPrintSteps) )

      ! Allocate and initialize trapping probability
      ALLOCATE( TrappingProb(0:NrOfPrintSteps, 0:NRhoMax) )
      TrappingProb = 0.0

!          ! if XYZ files of the trajectories are required, allocate memory to store the traj
!          IF ( PrintType >= FULL  .AND. ( RunType == SCATTERING .OR. RunType == EQUILIBRIUM ) ) THEN
!                ALLOCATE( Trajectory( 16, ntime ) )
!          END IF
! 

   END SUBROUTINE Scattering_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the scattering simulation and compute trapping probability over time.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Run()
      IMPLICIT NONE
      INTEGER  ::  AvHydroOutputUnit, AvCarbonOutputUnit       ! UNITs FOR OUTPUT AND DEBUG
      INTEGER  ::  jRho, iTraj, iStep, kStep
      REAL     ::  ImpactPar, Time
      REAL, DIMENSION(1,8)  :: AsymptoticCH
      REAL     ::  TotEnergy, PotEnergy, KinEnergy, IstTemperature
      REAL     ::  TempAverage, TempVariance

! 
!       REAL, DIMENSION(121) :: vzsum, zsum
!       REAL, DIMENSION(:,:), ALLOCATABLE :: vz2av, zcav
!       REAL :: vav, vsqav, zav, zsqav
! 
!       ! Initialize average values for the INITIAL CONDITION of the lattice Z coordinates
!       zav    = 0.0         ! average position
!       zsqav  = 0.0         ! mean squared position
!       vav    = 0.0         ! average velocity
!       vsqav  = 0.0         ! mean squared velocity
! 
!       ALLOCATE( zcav(10,ntime), vz2av(10,ntime) )
!       zcav(:,:)  = 0.0
!       vz2av(:,:) = 0.0
! 
!       ! allocate and initialize trapping probability variable
!       ALLOCATE( ptrap( irho+1, ntime ) )
!       ptrap(:,:) = 0.0
! 

! 
      AvHydroOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageHydro.dat", UNIT=AvHydroOutputUnit )
      WRITE(AvHydroOutputUnit, "(A,I6,A,/)") "# average zH coord and vel vs time (fs | ang) - ",NrTrajs, " trajs"

      AvCarbonOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageCarbon.dat", UNIT=AvCarbonOutputUnit )
      WRITE(AvCarbonOutputUnit, "(A,I6,A,/)") "# average zC coord and vel vs time (fs | ang) - ",NrTrajs," trajs"


      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "         H-GRAPHITE SCATTERING SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A,I5,A)"," Running ", NrTrajs, " trajectories per ",NRhoMax," impact parameters ... "

      ! Initialize variables to print average results at a given impact parameter value
      AverageVCarbon(:)         = 0.0
      AverageZCarbon(:)         = 0.0
      AverageVHydro(:)          = 0.0
      AverageZHydro(:)          = 0.0

      ! scan over the impact parameter... 
      ImpactParameter: DO jRho = 0,NRhoMax

         ! Set impact parameter and print message
         ImpactPar = float(jRho)*DeltaRho

         PRINT "(2/,A)",    "***************************************************"
         PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar*LengthConversion(InternalUnits,InputUnits)
         PRINT "(A,/)" ,    "***************************************************"

         PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

         ! run NTrajs number of trajectories at the current impact parameter
         Trajectory: DO iTraj = 1, NrTrajs

            !*************************************************************
            ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
            !*************************************************************

            ! Set initial coordinates of C and H atoms
            AsymptoticCH(1,1) = ImpactPar
            AsymptoticCH(1,2) = 0.0
            AsymptoticCH(1,3) = ZHInit
            AsymptoticCH(1,4) = 0.0

            ! Set initial velocities of C and H atoms
            AsymptoticCH(1,5) = 0.0
            AsymptoticCH(1,6) = 0.0
            AsymptoticCH(1,7) = sqrt( 2.0 * EKinZH / MassH )
            AsymptoticCH(1,8) = 0.0

            ! Set initial conditions
            IF ( BathType == SLAB_POTENTIAL ) THEN 
               CALL ZeroKelvinSlabConditions( X, V, AsymptoticCH, RandomNr ) 
            ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
               CALL ZeroKelvinBathConditions( Bath, X(5:), V(5:), ZPECorrection, RandomNr )
            ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
               CALL ZeroKelvinBathConditions( DblBath(1), X(5:NBath+4), V(5:NBath+4), ZPECorrection, RandomNr )
               CALL ZeroKelvinBathConditions( DblBath(2), X(NBath+5:2*NBath+4), V(NBath+5:2*NBath+4), ZPECorrection, RandomNr )
            ELSE IF ( BathType == LANGEVIN_DYN ) THEN
               ! nothing to do
            END IF

            X(1:4) = AsymptoticCH( 1, 1:4 )
            V(1:4) = AsymptoticCH( 1, 5:8 )

            IF ( BathType /= LANGEVIN_DYN ) THEN 

               PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ",  &
                                            Temperature*TemperatureConversion(InternalUnits,InputUnits)

               ! Initialize temperature average and variance
               TempAverage = 0.0
               TempVariance = 0.0

               ! Compute starting potential and forces
               A(:) = 0.0
               PotEnergy = ScatteringPotential( X, A )
               A(:) = A(:) / MassVector(:)

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/size(X)

               ! Do an equilibration run
               EquilibrationCycle: DO iStep = 1, NrEquilibSteps

                  ! PROPAGATION for ONE TIME STEP
                  CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, ScatteringPotential, PotEnergy, RandomNr )

                  ! compute kinetic energy and total energy
                  KinEnergy = EOM_KineticEnergy(Equilibration, V )
                  TotEnergy = PotEnergy + KinEnergy
                  IstTemperature = 2.0*KinEnergy/size(X)

                  ! store temperature average and variance
                  TempAverage = TempAverage + IstTemperature
                  TempVariance = TempVariance + IstTemperature**2

               END DO EquilibrationCycle

               PRINT "(A)", " Equilibration completed! "

               X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
               X(1:4) = AsymptoticCH( 1, 1:4 )
               V(1:4) = AsymptoticCH( 1, 5:8 )

            END IF 

            ! Compute average and standard deviation
            TempAverage = TempAverage / NrEquilibSteps 
            TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
            ! output message with average values
            PRINT "(/,A,1F10.4,A,/,A,1F10.4,A,/)",  " * Average temperature: ", &
                        TempAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
                                                    " * Standard deviation: ", &
                    sqrt(TempVariance)*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)


!           !*************************************************************
!           ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
!           !*************************************************************

            ! Compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            TotEnergy = PotEnergy + KinEnergy

            ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
            WRITE(*,600)  PotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          KinEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

            ! Increment averages at starting conditions
            AverageVCarbon(0) = AverageVCarbon(0) + V(4)
            AverageZCarbon(0) = AverageZCarbon(0) + X(4)
            AverageVHydro(0)  = AverageVHydro(0)  + V(3)
            AverageZHydro(0)  = AverageZHydro(0)  + X(3)


            !*************************************************************
            !         TIME EVOLUTION OF THE TRAJECTORY
            !*************************************************************
    
            ! initialize counter for printing steps
            kStep = 0

            PRINT "(/,A)", " Propagating the H-Graphene system in time... "

            ! cycle over nstep velocity verlet iterations
            Propagation: DO iStep = 1, NrOfSteps

               ! Propagate for one timestep
               CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, ScatteringPotential, PotEnergy, RandomNr )

               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               IF ( BathType == SLAB_POTENTIAL ) THEN 
                  X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
               ENDIF

               ! Compute kin energy and temperature
               KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

               ! output to write every nprint steps 
               IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

                  ! increment counter for printing steps
                  kStep = kStep+1
                  IF ( kStep > NrOfPrintSteps ) CYCLE 

                  ! increment ch coords averages
                  AverageVCarbon(kStep) = AverageVCarbon(kStep) + V(4)
                  AverageZCarbon(kStep) = AverageZCarbon(kStep) + X(4)
                  AverageVHydro(kStep)  = AverageVHydro(kStep)  + V(3)
                  AverageZHydro(kStep)  = AverageZHydro(kStep)  + X(3)
                  
                  ! check if H is in the trapping region
                  IF ( X(3) <= AdsorpLimit ) TrappingProb(jRho,kStep) = TrappingProb(jRho,kStep)+1.0

   !                ! Store the trajectory for XYZ printing
   !                IF ( PrintType >= FULL  .AND. RunType == EQUILIBRIUM ) THEN
   !                      Trajectory( :, kstep ) = 0.0
   !                      Trajectory( 1:min(16,3+nevo) , kstep ) = X( 1:min(16,3+nevo) ) 
   !                      NrOfTrajSteps = kstep
   !                END IF

               END IF 
    
            END DO Propagation

!             ! At the end of the propagation, write the xyz file of the trajectory, if requested
!             IF ( PrintType >= FULL ) THEN
! 
!                ! print the full trajectory as output file
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",jRho,"_",iTraj,".xyz"
!                CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
!                                              GraphiteLatticeConstant()*MyConsts_Bohr2Ang )
! 

            PRINT "(A)", " Time propagation completed! "

            ! print the final values of this trajectory 
            WRITE(*,700) X(3)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InternalUnits), &
                         X(4)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InternalUnits), &
                         TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InternalUnits)

         END DO Trajectory

         PRINT "(A)"," Done! "                                    

         ! Normalize coordinates averages 
         AverageVCarbon(:) = AverageVCarbon(:)/NrTrajs
         AverageZCarbon(:) = AverageZCarbon(:)/NrTrajs 
         AverageVHydro(:)  = AverageVHydro(:)/NrTrajs
         AverageZHydro(:)  = AverageZHydro(:)/NrTrajs

         ! write impact parameter header in output files
         WRITE(AvHydroOutputUnit,"(/,A,F10.5)")  &
                    "# Impact parameter : ",ImpactPar*LengthConversion(InternalUnits,InputUnits)
         WRITE(AvCarbonOutputUnit,"(/,A,F10.5)")  &
                    "# Impact parameter : ",ImpactPar*LengthConversion(InternalUnits,InputUnits)

         ! Print time resolved data to output files
         DO iStep = 1, NrOfPrintSteps
            Time = FLOAT(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits)
            WRITE(AvHydroOutputUnit,800) Time, AverageZHydro(iStep)*LengthConversion(InternalUnits,InputUnits),  &
                                             AverageVHydro(iStep)*VelocityConversion(InternalUnits,InputUnits)
            WRITE(AvCarbonOutputUnit,800) Time, AverageZCarbon(iStep)*LengthConversion(InternalUnits,InputUnits), &     
                                             AverageVCarbon(iStep)*VelocityConversion(InternalUnits,InputUnits)
         END DO

      END DO ImpactParameter

      ! Normalize trapping probability 
      TrappingProb(jRho,kStep) = TrappingProb(jRho,kStep) / NrTrajs

      ! ****************************************
      ! PRINT TO FILE THE TRAPPING PROB !!!!!!
      ! ****************************************


!       ! write cross section data 
!       write(6,*)'cross section vs time(ps):'
! 
!       ! Cycle over analysis steps
!       DO iStep=1,ntime
! 
!          ! integrate over rxn parameter
!          crscn(iStep)=0.0
!          DO jRho=1,irho+1
!             ImpactPar = 0.000001 + delrho * real(jRho-1) * MyConsts_Bohr2Ang
!             crscn(iStep)=crscn(iStep)+(ptrap(jRho,iStep)/float(inum))* ImpactPar
!          END DO
!          crscn(iStep) = 2.0* MyConsts_PI * delrho*MyConsts_Bohr2Ang * crscn(iStep)
! 
!          ! print trapping p to file
!          write(6,605) dt*real(iStep*nprint)/MyConsts_fs2AU,  crscn(iStep)
!       END DO


   600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                  " * Potential Energy (eV)        ",1F10.4,/    &
                  " * Kinetic Energy (eV)          ",1F10.4,/    &
                  " * Total Energy (eV)            ",1F10.4,/ ) 

   700 FORMAT (/, " * Final zH coord (Ang)         ",1F10.4,/    &
                  " * Final energy (eV)          ",1F10.4,/    &
                  " * Final zC coord (Ang)         ",1F10.4,/ ) 

   800 FORMAT(F12.5,1000F15.8)

   END SUBROUTINE Scattering_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the simulation.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Dispose()
      IMPLICIT NONE


   END SUBROUTINE Scattering_Dispose

!*************************************************************************************************

   REAL FUNCTION ScatteringPotential( Positions, Forces )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i       

      ! Check the number degrees of freedom
      NrDOF = size( Positions )
      CALL ERROR( size(Forces) /= NrDOF, "ScatteringSimulation.ScatteringPotential: array dimension mismatch" )
      CALL ERROR( NrDOF /= NDim, "ScatteringSimulation.ScatteringPotential: wrong number of DoFs" )

      ! Initialize forces and potential
      ScatteringPotential = 0.0
      Forces(:)          = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         ScatteringPotential = VHSticking( Positions, Forces )

      ELSE 
         ! Compute potential and forces of the system
          ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )

         ! Add potential and forces of the bath and the coupling
         IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
            CALL BathPotentialAndForces( Bath, Positions(4)-C1Puckering, Positions(5:), &
                          ScatteringPotential, Forces(4), Forces(5:) ) 
         ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
            CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), &
                         ScatteringPotential, Forces(4), Forces(5:NBath+4) ) 
            CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), &
                         ScatteringPotential, Forces(4), Forces(NBath+5:2*NBath+4) ) 
         END IF

      END IF

   END FUNCTION ScatteringPotential

!*************************************************************************************************

END MODULE ScatteringSimulation
