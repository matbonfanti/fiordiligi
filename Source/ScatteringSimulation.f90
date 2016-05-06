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
!>  \version          2.0
!>  \date             10 October 2013
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg 28 January 2015 : the module has been completely revise according to the new
!>                          structure of the code (version 2.0)
!>  \arg 28 January 2015 : * collinear / non collinear distinction
!>                         * dispose subroutine completed
!>  \arg 3 February 2015 : * debug printing is now available
!>                         * zC-f(ZH) coupling function has been added 
!>                         * thermal initial conditions before equilibration
!>  \arg 25 February 2015: * quasiclassical simulation in normal bath case
!>  \arg 6 May 2016      : * restart run is now implemented
!                
!>  \todo         print XYZ file of scattering trajectory
!>  \todo         print coupling function and its derivative
!
!***************************************************************************************
!
!>  \attention   The random number generator is restarted with the status
!>               at the end of the previous run. When running a non collinear
!>               calculation (NRhoMax > 0) implies that a single run is different
!>               from two runs with the same total number of trajectories, 
!>               because the same sequence of random number is used for 
!>               trajectories of different rho in a different order.
!
!***************************************************************************************
MODULE ScatteringSimulation
#include "preprocessoptions.cpp"
   USE SharedData
   USE InputField
   USE UnitConversion
   USE ClassicalEqMotion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator
   USE PrintTools
   USE SplineInterpolator

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

   ! Variables for restart capabilities
   LOGICAL         :: Restart           !< Logical flag to switch restart on
   CHARACTER(100)  :: ProbRestartFile   !< Name of the file with probabilities of previous run
   CHARACTER(100)  :: RandRestartFile   !< Name of the file with random number generator internal state
   INTEGER         :: NrTrajRestart     !< Nr of trajectories of previous run

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

   ! VTK file format
   TYPE(VTKInfo), SAVE :: TrappingFile         !< derived datatype to print trapping probability

   ! Number of dimensions of the system+bath hamiltonian
   INTEGER :: NDim                             !< integer nr of dofs of the system+bath hamiltonian

#if defined(__ZH_DEPENDENT_COUPLING)
   ! ZCeq as a function of ZH
   CHARACTER(100)   :: ZCofZHFile              !< file storing the data ZC vs ZH (in atomic units)
   TYPE(SplineType), SAVE :: ZCofZHSpline      !< spline function defining ZC as a function of ZH
#endif

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
      IF ( NRhoMax > 0 ) THEN
         CALL SetFieldFromInput( InputData, "DeltaRho", DeltaRho )
         DeltaRho = DeltaRho * LengthConversion(InputUnits, InternalUnits)
      ELSE  ! for collinear calculation, no deltarho is needed
         DeltaRho = 0.0
      ENDIF

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

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Temperature of the equilibrium simulation
      CALL SetFieldFromInput( InputData, "Temperature", Temperature, 0.0 )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)
      ! Set temperature to zero in case of a quasiclassical simulation
      IF ( (BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH) .AND. ZPECorrection ) Temperature = 0.0

      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilRelaxTime", EquilGamma)
      EquilGamma = 1. / ( EquilGamma * TimeConversion(InputUnits, InternalUnits) )

      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
      EquilTStep = EquilTStep * TimeConversion(InputUnits, InternalUnits)

      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(10.0*(1.0/EquilGamma)/EquilTStep) )
      EquilibrationStepInterval = CEILING( real(NrEquilibSteps) / real(NrOfPrintSteps) )

#if defined(__ZH_DEPENDENT_COUPLING)
      ! EQUILIBRIUM POSITION OF C ATOM AS FUNCTION OF H
      IF ( BathType /= LANGEVIN_DYN .AND. BathType /= SLAB_POTENTIAL ) &
                       CALL SetFieldFromInput( InputData, "ZCofZHFile", ZCofZHFile )
#endif

      ! RESTART VARIABLES 

      ! flat to set the restart on
      CALL SetFieldFromInput( InputData, "Restart", Restart, .FALSE. )
      IF ( Restart ) THEN
         ! nr of trajectories of previous run
         CALL SetFieldFromInput( InputData, "NrTrajRestart", NrTrajRestart )
         IF ( NRhoMax > 0 ) THEN             ! depending of collinear/non collinear, default file name is different
            ! file with trapping probability as a function of impact parameter
            CALL SetFieldFromInput( InputData, "ProbRestartFile", ProbRestartFile, "Trapping.vtr" )
         ELSE
            ! trapping probability only for collinear approach
            CALL SetFieldFromInput( InputData, "ProbRestartFile", ProbRestartFile, "CollinearTrapping.dat" )
         END IF
         ! Random number generator restart file
         CALL SetFieldFromInput( InputData, "RandRestartFile", RandRestartFile, "RandomNr.dat" )
      ELSE
         ! with no restart, set NrTrajRestart to zero
         NrTrajRestart = 0
      END IF

      WRITE(*, 902) ZHInit*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                    EKinZH*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                    NRhoMax, DeltaRho*NRhoMax*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits)

      WRITE(*, 903) NrTrajs, TimeStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                    NrOfSteps, NrOfPrintSteps

      WRITE(*, 904) 1./EquilGamma*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 EquilTStep*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 EquilTStep*NrEquilibSteps*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits), &
                 NrEquilibSteps

#if defined(__ZH_DEPENDENT_COUPLING)
      IF ( BathType /= LANGEVIN_DYN .AND. BathType /= SLAB_POTENTIAL ) WRITE(*, 905) TRIM(ADJUSTL(ZCofZHFile))
#endif

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

   905 FORMAT(" * File with zC_eq as a function of zH:         ",A,/)

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
      INTEGER :: ZCEqUnit
      INTEGER :: iCoord, iData, NData
      REAL, DIMENSION(:), ALLOCATABLE :: ZHGrid, ZCValue

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
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set canonical dynamics at the end of the oscillator chain or chains
      IF ( BathType == CHAIN_BATH .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( NDim ) = .TRUE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF
      IF ( BathType == DOUBLE_CHAIN .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( 4+NBath ) = .TRUE.  ! end of the first chain
         LangevinSwitchOn( NDim )    = .TRUE.  ! end of the second chain
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! In case of Langevin relaxation, switch on gamma at the carbon atom
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
      CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
      LangevinSwitchOn = (/ (.FALSE., iCoord=1,3), (.TRUE., iCoord=4,NDim ) /)
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature, LangevinSwitchOn )

      DEALLOCATE( LangevinSwitchOn )

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! Allocate and initialize the variables for the trajectory averages
      ALLOCATE( AverageVHydro(0:NrOfPrintSteps), AverageZHydro(0:NrOfPrintSteps) )
      ALLOCATE( AverageVCarbon(0:NrOfPrintSteps), AverageZCarbon(0:NrOfPrintSteps) )

      ! Allocate and initialize trapping probability
      ALLOCATE( TrappingProb(0:NRhoMax,NrOfPrintSteps) )
      TrappingProb = 0.0

      ! Read from restart file the trapping probability and remove normalization
      IF ( Restart ) THEN
         CALL ERROR( NrTrajRestart > NrTrajs, "Scattering_Initialize: in a restart calculation NrTrajs > NrTrajRestart" )
         CALL ReadForRestart( TrappingProb, RandomNr ) 
         TrappingProb(:,:) = TrappingProb(:,:) * NrTrajRestart
      END IF

      ! Check that a non collinear V is considered when running a non collinear simulation
      IF (NRhoMax > 0) THEN
         CALL ERROR( Collinear, " Scattering_Initialize: a non collinear potential is needed!" )
      END IF 

#if defined(__ZH_DEPENDENT_COUPLING)
      ! Define spline function with ZCeq as a function of ZH
      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN

         ! Count available data and allocate memory
         NData = CountLinesInFile( TRIM(ADJUSTL(ZCofZHFile)) )
         ALLOCATE( ZHGrid(NData), ZCValue(NData) )

         ! open file and read content
         OPEN( FILE=TRIM(ADJUSTL(ZCofZHFile)) , UNIT= ZCEqUnit )
         DO iData = 1, NData
            READ(ZCEqUnit,*) ZHGrid(iData), ZCValue(iData)
         END DO 
         CLOSE( ZCEqUnit )

         ! Define spline interpolating function and deallocate memory
         CALL SetupSpline( ZCofZHSpline, ZHGrid, ZCValue )
         DEALLOCATE( ZHGrid, ZCValue )

      END IF
#endif

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
      INTEGER  ::  CollinearTrapUnit, CrossSectionUnit, OpacityUnit
      INTEGER  ::  jRho, iTraj, iStep, kStep, iCoord
      INTEGER  ::  DebugUnitEn, DebugUnitCoord, DebugUnitVel
      REAL     ::  ImpactPar, Time, CrossSection
      REAL     ::  RndValue, CarbonFreq
      REAL, DIMENSION(1,8)  :: AsymptoticCH
      REAL     ::  TotEnergy, PotEnergy, KinEnergy, KinHydro, KinCarbon
      REAL     ::  TempAverage, TempVariance, IstTemperature
      REAL, DIMENSION(NRhoMax+1) :: ImpactParameterGrid
      REAL, DIMENSION(NrOfPrintSteps) :: TimeGrid
      CHARACTER(100) :: OutFileName

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "         H-GRAPHITE SCATTERING SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A,I5,A)"," Running ", NrTrajs, " trajectories per ",NRhoMax+1," impact parameters ... "
      IF ( Restart ) PRINT "(A,I5,A,I5,A)"," Restarting run with ", NrTrajRestart, " and computing additional ",  &
                                            NrTrajs-NrTrajRestart," trajectories per impact parameter"

      ! Initialize variables to print average results at a given impact parameter value
      AverageVCarbon(:)         = 0.0
      AverageZCarbon(:)         = 0.0
      AverageVHydro(:)          = 0.0
      AverageZHydro(:)          = 0.0

      ! scan over the impact parameter... 
      ImpactParameter: DO jRho = 0,NRhoMax

         ! Set impact parameter and print message
         ImpactPar = float(jRho)*DeltaRho

         IF ( NRhoMax == 0 ) THEN 
            PRINT "(/,A,I5,A)"," Collinear scattering simulation ... "
         ELSE
            PRINT "(2/,A)",    "***************************************************"
            PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar*LengthConversion(InternalUnits,InputUnits)
            PRINT "(A,/)" ,    "***************************************************"
         END IF

         PRINT "(A,I5,A)"," Running ", NrTrajs-NrTrajRestart, " trajectories ... "

         ! run NTrajs number of trajectories at the current impact parameter
         Trajectory: DO iTraj = NrTrajRestart+1, NrTrajs

            PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

            !*************************************************************
            ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
            !*************************************************************

            ! Put zero point energy in each oscillator of the bath and in the carbon
            IF ( (BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH) .AND. ZPECorrection ) THEN

               ! Set initial coordinates of H atom
               X(1) = ImpactPar
               X(2) = 0.0
               X(3) = ZHInit

               ! Set initial velocities of H atom
               V(1) = 0.0
               V(2) = 0.0
               V(3) = -sqrt( 2.0 * EKinZH / MassH )
 
               ! Set initial coord and vel of C atom
               RndValue = UniformRandomNr( RandomNr )*2.0*MyConsts_PI
               CarbonFreq = SQRT( CarbonForceConstant() / MassC )
               X(4) = COS(RndValue)/SQRT(MassC*CarbonFreq)
               V(4) = SIN(RndValue)*SQRT(CarbonFreq/MassC)

               ! Set initial conditions of the bath
               CALL ZeroKelvinBathConditions( Bath, X(5:), V(5:), ZPECorrection, RandomNr )

               ! Compute starting potential and forces
               A(:) = 0.0
               PotEnergy = ScatteringPotential( X, A )
               A(:) = A(:) / MassVector(:)

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/size(X)

            ELSE  ! Thermalize bath + carbon

               ! Set initial coordinates of C and H atoms
               AsymptoticCH(1,1) = ImpactPar
               AsymptoticCH(1,2) = 0.0
               AsymptoticCH(1,3) = ZHInit
               AsymptoticCH(1,4) = 0.0

               ! Set initial velocities of C and H atoms
               AsymptoticCH(1,5) = 0.0
               AsymptoticCH(1,6) = 0.0
               AsymptoticCH(1,7) = -sqrt( 2.0 * EKinZH / MassH )
               AsymptoticCH(1,8) = 0.0

               ! Set initial conditions
               IF ( BathType == SLAB_POTENTIAL ) THEN 
                  CALL ZeroKelvinSlabConditions( X, V, AsymptoticCH, RandomNr ) 
               ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
                  CALL ThermalEquilibriumBathConditions( Bath, X(5:), V(5:), Temperature, RandomNr )
               ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
                  CALL ThermalEquilibriumBathConditions( DblBath(1), X(5:NBath+4), V(5:NBath+4), Temperature, RandomNr )
                  CALL ThermalEquilibriumBathConditions( DblBath(2), X(NBath+5:2*NBath+4), V(NBath+5:2*NBath+4), &
                                                                              Temperature, RandomNr )
               ELSE IF ( BathType == LANGEVIN_DYN ) THEN
                  ! nothing to do
               END IF

               ! During equilibration, fix H in asymptotic position with null velocity
               ! and initialize carbon atom in the equilibrium position with null velocity
               X(1:4) = AsymptoticCH( 1, 1:4 )
               V(1:4) = 0.0

               IF ( BathType /= LANGEVIN_DYN ) THEN 

                  PRINT "(/,A,F6.1,1X,A)",   " Equilibrating the initial conditions at T = ", &
                      Temperature*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)

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

                  ! Compute average and standard deviation
                  TempAverage = TempAverage / NrEquilibSteps 
                  TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
                  ! output message with average values
                  WRITE(*,500)  TempAverage*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits), &
                                sqrt(TempVariance)*TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)

               END IF 

               ! for the slab potential translate C2-C3-C4 plane at z=0
               IF ( BathType == SLAB_POTENTIAL )   X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
               ! give right initial conditions to the H atom
               X(1:3) = AsymptoticCH( 1, 1:3 )
               V(1:3) = AsymptoticCH( 1, 5:7 )

            END IF


!           !*************************************************************
!           ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
!           !*************************************************************

            ! Compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            KinHydro  = EOM_KineticEnergy( MolecularDynamics, V, 3 )
            KinCarbon = KinEnergy - KinHydro
            TotEnergy = PotEnergy + KinEnergy

            ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
            WRITE(*,600)  PotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          KinHydro*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),  &
                          KinCarbon*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          KinEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits), &
                          TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

            ! Increment averages at starting conditions
            AverageVCarbon(0) = AverageVCarbon(0) + V(4)
            AverageZCarbon(0) = AverageZCarbon(0) + X(4)
            AverageVHydro(0)  = AverageVHydro(0)  + V(3)
            AverageZHydro(0)  = AverageZHydro(0)  + X(3)

            ! Open unit for massive output, with detailed info on trajectories
            IF ( PrintType == DEBUG ) THEN
               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Energy.dat"
               DebugUnitEn = LookForFreeUnit()
               OPEN( Unit=DebugUnitEn, File=OutFileName )

               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Coord.dat"
               DebugUnitCoord = LookForFreeUnit()
               OPEN( Unit=DebugUnitCoord, File=OutFileName )

               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Rho_",jRho,"_Traj_",iTraj,"_Vel.dat"
               DebugUnitVel = LookForFreeUnit()
               OPEN( Unit=DebugUnitVel, File=OutFileName )

               ! Write initial values
               WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | V, KH, KC, Kin, E / Eh "
               WRITE(DebugUnitEn,800) 0.0,  PotEnergy, KinHydro, KinCarbon, KinEnergy, TotEnergy

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

            ! cycle over nstep velocity verlet iterations
            Propagation: DO iStep = 1, NrOfSteps

               ! Propagate for one timestep
               CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, ScatteringPotential, PotEnergy, RandomNr )

               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               IF ( BathType == SLAB_POTENTIAL ) THEN 
                  X(3:NDim) = X(3:NDim)  - (X(5)+X(6)+X(7))/3.0
               ENDIF

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
                  
                  ! Compute kinetic energy and total energy
                  KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
                  KinHydro  = EOM_KineticEnergy( MolecularDynamics, V, 3 )
                  KinCarbon = KinEnergy - KinHydro
                  TotEnergy = PotEnergy + KinEnergy

                  ! check if H is in the trapping region
                  IF ( X(3) <= AdsorpLimit ) TrappingProb(jRho,kStep) = TrappingProb(jRho,kStep)+1.0

   !                ! Store the trajectory for XYZ printing
   !                IF ( PrintType >= FULL  .AND. RunType == EQUILIBRIUM ) THEN
   !                      Trajectory( :, kstep ) = 0.0
   !                      Trajectory( 1:min(16,3+nevo) , kstep ) = X( 1:min(16,3+nevo) ) 
   !                      NrOfTrajSteps = kstep
   !                END IF

                  ! If massive level of output, print traj information to std out
                  IF ( PrintType == DEBUG ) THEN
                     WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                           PotEnergy, KinHydro, KinCarbon, KinEnergy, TotEnergy
                     WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                             (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, NDim-4) /)
                     WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
                  END IF

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
            WRITE(*,700) X(3)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                         X(4)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits), &
                         TotEnergy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits)

            IF ( PrintType == DEBUG ) THEN
                  CLOSE( Unit=DebugUnitEn )
                  CLOSE( Unit=DebugUnitCoord )
                  CLOSE( Unit=DebugUnitVel )
            ENDIF

         END DO Trajectory

         PRINT "(A)"," Done! "                                    

         ! Normalize coordinates averages 
         AverageVCarbon(:) = AverageVCarbon(:)/(NrTrajs-NrTrajRestart)
         AverageZCarbon(:) = AverageZCarbon(:)/(NrTrajs-NrTrajRestart)
         AverageVHydro(:)  = AverageVHydro(:)/(NrTrajs-NrTrajRestart)
         AverageZHydro(:)  = AverageZHydro(:)/(NrTrajs-NrTrajRestart)

         IF ( PrintType >= FULL ) THEN

            AvHydroOutputUnit = LookForFreeUnit()
            OPEN( FILE="AverageHydro.dat", UNIT=AvHydroOutputUnit )
            WRITE(AvHydroOutputUnit, "(A,I6,A,/)") "# average zH coord and vel vs time (fs | ang) - ",NrTrajs, " trajs"

            AvCarbonOutputUnit = LookForFreeUnit()
            OPEN( FILE="AverageCarbon.dat", UNIT=AvCarbonOutputUnit )
            WRITE(AvCarbonOutputUnit, "(A,I6,A,/)") "# average zC coord and vel vs time (fs | ang) - ",NrTrajs," trajs"

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

            CALL WARN( Restart, " Scattering_Run: average coord.s and vel.s computed only for the trajs of the restart run " )  

         END IF

      END DO ImpactParameter

      ! Normalize trapping probability 
      TrappingProb(:,:) = TrappingProb(:,:) / NrTrajs

      ! Prepare time grid and impact parameter grid
      TimeGrid = (/( float(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits), &
                                      iStep= 1,NrOfPrintSteps  )/)
      ImpactParameterGrid = (/( float(jRho)*DeltaRho*LengthConversion(InternalUnits,InputUnits), jRho = 0,NRhoMax  )/)

      IF ( NRhoMax > 0 ) THEN
         ! Print trapping to vtk file
         CALL VTK_NewRectilinearSnapshot( TrappingFile, X=ImpactParameterGrid, Y=TimeGrid, FileName="Trapping" )
         CALL VTK_AddScalarField( TrappingFile, "trapping_p", RESHAPE( TrappingProb(:,:), (/NrOfPrintSteps*(NRhoMax+1) /))  )
      END IF

      ! Print to file collinear trapping
      CollinearTrapUnit = LookForFreeUnit()
      OPEN( FILE="CollinearTrapping.dat", UNIT=CollinearTrapUnit )
      WRITE(CollinearTrapUnit, "(A,I6,A,/)") "# collinear trapping probability (" // TRIM(TimeUnit(InputUnits)) // &
          " | adim) -",NrTrajs, " trajs"
      DO iStep= 1,NrOfPrintSteps
         WRITE(CollinearTrapUnit,800) TimeGrid(iStep), TrappingProb(0,iStep)
      END DO

      ! Close open units
      CLOSE( CollinearTrapUnit )
      IF ( PrintType >= FULL ) CLOSE( AvCarbonOutputUnit )
      IF ( PrintType >= FULL ) CLOSE( AvHydroOutputUnit )

      ! write cross section data 
      IF ( NRhoMax > 0 ) THEN

         ! open file and write header
         CrossSectionUnit = LookForFreeUnit()
         OPEN( FILE="CrossSection.dat", UNIT=CrossSectionUnit )
         WRITE(CrossSectionUnit, "(A,I6,A,/)") "# cross section vs time (" // TRIM(TimeUnit(InputUnits)) // &
             " | " // TRIM(LengthUnit(InputUnits)) //"^2) -",NrTrajs, " trajs"

         ! loop over time steps, compute and write trapping cross section
         DO iStep= 1,NrOfPrintSteps

            ! integrate over impact parameter
            CrossSection=0.0
            DO jRho = 0, NRhoMax
               CrossSection = CrossSection + TrappingProb(jRho,iStep) * float(jRho)*DeltaRho
            END DO
            CrossSection = 2.0*MyConsts_PI * DeltaRho * CrossSection * LengthConversion(InternalUnits,InputUnits)**2

            ! print trapping cross section to file
            WRITE(CrossSectionUnit,800) TimeGrid(iStep), CrossSection

         END DO

         ! Print to file final opacity function
         OpacityUnit = LookForFreeUnit()
         OPEN( FILE="OpacityFunction.dat", UNIT=OpacityUnit )
         WRITE(OpacityUnit, "(A,F8.2,A,/)") "# opacity function @ time ", &
            float(PrintStepInterval*iStep)*TimeStep*TimeConversion(InternalUnits,InputUnits), " "//TRIM(TimeUnit(InputUnits))
         DO jRho = 0,NRhoMax
           WRITE(OpacityUnit,800) ImpactParameterGrid(jRho+1), TrappingProb(jRho,NrOfPrintSteps) 
         END DO
         CLOSE(OpacityUnit)

      END IF

      ! dump the random number generator internal status to file
      CALL DumpInternalState( RandomNr, "RandomNr.dat" )


   500 FORMAT (/, " Equilibration averages                 ",/     &
                  " * Average temperature          ",1F10.4,1X,A,/ &
                  "   Standard deviation           ",1F10.4,1X,A,/ ) 

   600 FORMAT (/, " Initial condition of the MD trajectory ",/     &
                  " * Potential Energy             ",1F10.4,1X,A,/ &
                  " * Kinetic Energy of H atom     ",1F10.4,1X,A,/ &
                  " * Kinetic Energy of C atoms    ",1F10.4,1X,A,/ &
                  " * Total Kinetic Energy         ",1F10.4,1X,A,/ &
                  " * Total Energy                 ",1F10.4,1X,A,/ ) 

   700 FORMAT (/, " * Final zH coord               ",1F10.4,1X,A,/ &
                  " * Final zC coord               ",1F10.4,1X,A,/ &
                  " * Final energy                 ",1F10.4,1X,A,/ ) 

   800 FORMAT(F12.5,1000F15.8)

   END SUBROUTINE Scattering_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the simulation.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Dispose()
      IMPLICIT NONE

      ! Deallocate memory
      DEALLOCATE( X, V, A, MassVector )
      DEALLOCATE( AverageVHydro, AverageZHydro, AverageVCarbon, AverageZCarbon, TrappingProb )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

#if defined(__ZH_DEPENDENT_COUPLING)
      ! Dispose spline function
      CALL DisposeSpline( ZCofZHSpline )
#endif

   END SUBROUTINE Scattering_Dispose

!*************************************************************************************************

   REAL FUNCTION ScatteringPotential( Positions, Forces )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i
      REAL :: CouplingFs, DTimesX

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

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Coupling function of the system
#if defined(__ZH_DEPENDENT_COUPLING)
         CouplingFs = Positions(4)-GetSpline1D(ZCofZHSpline, Positions(3))
#else
         CouplingFs = Positions(4)
#endif
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, CouplingFs, Positions(5:), &
                                                ScatteringPotential, Forces(4), Forces(5:), DTimesEffMode=DTimesX ) 
#if defined(__ZH_DEPENDENT_COUPLING)
         ! Add forces on ZH coming from ZH-dependent distortion correction
         Forces(3) = Forces(3) + GetDistorsionForce(Bath)*CouplingFs*DerivSpline(ZCofZHSpline, Positions(3))
         Forces(3) = Forces(3) - DTimesX * DerivSpline(ZCofZHSpline, Positions(3))
#endif

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Compute potential and forces of the system
         ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), ScatteringPotential, &
                                                                               Forces(4), Forces(5:NBath+4) ) 
         CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), ScatteringPotential, &
                                                                               Forces(4), Forces(NBath+5:2*NBath+4) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential and forces of the system
         ScatteringPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )

      END IF

   END FUNCTION ScatteringPotential

!*************************************************************************************************

!*******************************************************************************
!> Restart subroutine: check that the restart is allowed and read the 
!> probability of previous trajectories from file.
!>
!*******************************************************************************
   SUBROUTINE ReadForRestart( TrappingProb, RandomNr  ) 
      IMPLICIT NONE
      TYPE( RNGInternalState), INTENT(INOUT) :: RandomNr
      REAL, DIMENSION(:,:), INTENT(INOUT)    :: TrappingProb

      INTEGER :: PReadUnit        ! input unit
      CHARACTER(150) :: Dummy     ! dummy string
      INTEGER :: ReadDataNr       ! nr of data values in the input file
      INTEGER :: Stat             ! integer i/o status
      REAL, DIMENSION(6) :: Numbers  ! dummy real variable
      INTEGER :: i                ! integer counter

      ! Read probability resolved over time and impact parameters from vtr file 
      IF ( NRhoMax > 0 ) THEN    

         ! Open file and look for the starting line of the "<PointData>" section
         OPEN(Unit=PReadUnit, FILE=ProbRestartFile)
         DO i = 1, 100000
            READ(PReadUnit,*) Dummy
            IF ( TRIM(ADJUSTL(Dummy)) == "<PointData>" ) EXIT
         END DO
         ! skip another line
         READ(PReadUnit,*) Dummy

         ! COUNT NR OF DATA IN THIS SECTION
         ReadDataNr = 0
         DO i = 1, 100000
            ! read the entire line and try to convert it to 6 real numbers
            READ(PReadUnit,"(A150)",IOSTAT=Stat) Dummy
            READ(Dummy,*,IOSTAT=Stat) Numbers
            ! when the conversion fails, exit from the cycle
            IF (Stat/=0) EXIT
            ! increment the number of input data
            ReadDataNr = ReadDataNr + 6
         END DO
         ! the last line (which gave the failure) might contain a smaller number of values
         DO i = 6,1,-1
            READ(Dummy,*,IOSTAT=Stat) Numbers(1:i)
            IF ( Stat == 0 ) THEN
                ReadDataNr = ReadDataNr + i 
                EXIT
            END IF
         END DO

         ! CHECK THAT THE NR OF DATA CORRESPONDS TO WHAT EXPECTED
         CALL ERROR( ReadDataNr /= size(TrappingProb,1)*size(TrappingProb,2), &
                " ScatteringSimulation.ReadForRestart: wrong number of values in file "//TRIM(ADJUSTL(ProbRestartFile)) )

         ! Again go to the beginning of the "<PointData>" section
         REWIND(PReadUnit)
         DO i = 1, 100000
            READ(PReadUnit,*) Dummy
            IF ( TRIM(ADJUSTL(Dummy)) == "<PointData>" ) EXIT
         END DO
         READ(PReadUnit,*) Dummy

         ! Read the probability
         READ(PReadUnit,"(1(10X,6(E15.8,X)))") TrappingProb(:,:)

         ! close file
         CLOSE(Unit=PReadUnit)

      ! Read probability resolved over time for collinear approach from dat file 
      ELSE IF ( NRhoMax == 0 ) THEN

         ! Counts the number of lines with data in the input file
         ReadDataNr = CountLinesInFile( ProbRestartFile ) - 1

         PRINT*, ReadDataNr, size(TrappingProb,1)*size(TrappingProb,2)

         ! CHECK THAT THE NR OF DATA CORRESPONDS TO WHAT EXPECTED
         CALL ERROR( ReadDataNr /= size(TrappingProb,1)*size(TrappingProb,2), &
                " ScatteringSimulation.ReadForRestart: wrong number of values in file "//TRIM(ADJUSTL(ProbRestartFile)) )

         ! Open file and skip two lines
         OPEN(Unit=PReadUnit, FILE=ProbRestartFile)
         READ(PReadUnit,*) 
         READ(PReadUnit,*) 

         ! Read the collinear sticking probability
         DO i = 1, ReadDataNr
            READ(PReadUnit,*) Numbers(1:2)
            TrappingProb(i,1) = Numbers(2)
         END DO

      END IF

      ! read internal status of the randon number generator from file 
      CALL LoadInternalState( RandomNr, RandRestartFile )

   END SUBROUTINE ReadForRestart 

END MODULE ScatteringSimulation
