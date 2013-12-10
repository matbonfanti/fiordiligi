!***************************************************************************************
!*                              MODULE Harmonic1DModel
!***************************************************************************************
!
!>  \brief     Subroutines for an equilibrium simulation of a model 1D oscillator
!>  \details   This module contains subroutines to compute averages  \n
!>             in an equilibrium simulation of a brownian 1D harmonic oscillator \n
!>             at a given temperature   \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             9 October 2013
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 23 October 2013: double chain bath has been implemented
!
!>  \todo            fix the log of the input variables
!>                 
!***************************************************************************************
MODULE Harmonic1DModel
   USE MyConsts
   USE ErrorTrap
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE ClassicalEqMotion
   USE IndependentOscillatorsModel
   USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Harmonic1DModel_ReadInput, Harmonic1DModel_Initialize, Harmonic1DModel_Run, Harmonic1DModel_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

   ! Variables of the system
   REAL    :: QFreq                     !< Frequency of the harmonic oscillator
   REAL    :: ForceConstant             !< Force constants

   ! Variables of the propagation
   INTEGER :: NrTrajs                   !< Nr of trajectories
   INTEGER :: NrOfSteps                 !< Total nr of time step per trajectory
   REAL    :: TimeStep                  !< Time step for the integration of the classical EOM
   INTEGER :: NrOfPrintSteps            !< Nr of analysis steps
   INTEGER :: PrintStepInterval         !< Nr of time steps between each printing of propagation info ( = NrOfSteps / NrOfPrintSteps )

   ! Initial conditions
   REAL    :: Temperature          !< Temperature of the simulation
   REAL    :: EquilTStep           !< Time step for integrating equilibration dynamics
   REAL    :: EquilGamma           !< Friction parameter of the Langevin equation
   INTEGER :: NrEquilibSteps       !< Nr of time step of the equilibration

   ! Time evolution dataset
   TYPE(Evolution),SAVE :: Equilibration         !< Propagate in macrocanonical ensamble at given T to generate init conditions
   TYPE(Evolution),SAVE :: MolecularDynamics     !< Propagate in micro/macrocanonical ensamble to extract results

   ! Averages and other output data
   REAL     :: dOmega                                      !< Frequency spacing of fourier transformed data
   REAL     :: GlobalTemperature                           !< Global temperature average
   REAL, DIMENSION(:), ALLOCATABLE  ::  XatT0              !< Initial coordinates to compute autocorrelation functions
   REAL, DIMENSION(:), ALLOCATABLE  ::  QAutoCorr          !< Autocorrelation function of the oscillator: < Q(0) Q(t) >
   REAL, DIMENSION(:), ALLOCATABLE  ::  AverageCoord       !< Average coordinates of each trajectory
   REAL, DIMENSION(:), ALLOCATABLE  ::  AverCoordOverTrajs !< Average coordinates of all the trajectories
   REAL, DIMENSION(:), ALLOCATABLE  ::  StDevCoord         !< Standard deviation of the coordinates of each trajectory
   COMPLEX, DIMENSION(:), ALLOCATABLE  :: PowerSpectrum    !< Complex array to compute the power spectrum of autocorr functions

   TYPE(RNGInternalState), SAVE :: RandomNr

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific for an equilibrium
!> simulation at given T.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE Harmonic1DModel_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! READ THE VARIABLES TO DEFINE THE SYSTEM PROPERTIES

      ! Frequency of the system
      CALL SetFieldFromInput( InputData, "QFreq", QFreq )
      QFreq = QFreq * MyConsts_cmmin1toAU
      ! Define accordingly the force constant
      ForceConstant = MassH * QFreq**2 

      ! READ THE VARIABLES FOR THE EQUILIBRIUM SIMULATION 

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

      ! READ THE VARIABLES TO SET THE INITIAL CONDITIONS OF THE SYSTEM

      ! Temperature of the equilibrium simulation
      CALL SetFieldFromInput( InputData, "Temperature", Temperature )
      Temperature = Temperature * MyConsts_K2AU

      ! Set gamma of the equilibration Langevin dynamics
      CALL SetFieldFromInput( InputData, "EquilGamma", EquilGamma)
      EquilGamma = EquilGamma / MyConsts_fs2AU

      ! Set the time step of the equilibration
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, TimeStep/MyConsts_fs2AU )
      EquilTStep = EquilTStep * MyConsts_fs2AU

      ! Set nr of steps of the equilibration
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(10.0*(1.0/EquilGamma)/EquilTStep) )

! 
!       WRITE(*, 903) InitEnergy*MyConsts_Hartree2eV, (InitEnergy-MinimumEnergy)*MyConsts_Hartree2eV, &
!                     NrOfInitSnapshots, TimeBetweenSnaps/MyConsts_fs2AU
! 
!       WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps
! 
!    903 FORMAT(" * Initial conditions of the atom-surface system ", /,&
!               " * Absolute initial energy (eV):                ",F10.4,/,& 
!               " *  - w.r.t. the bottom of the well (eV):       ",F10.4,/,& 
!               " * Nr of initial system snapshots:              ",I10,  /,& 
!               " * Time between snapshots (fs):                 ",F10.4,/ )
! 
!    904 FORMAT(" * Dynamical simulation variables               ", /,&
!               " * Nr of trajectories:                          ",I10,  /,& 
!               " * Propagation time step (fs):                  ",F10.4,/,& 
!               " * Nr of time steps of each trajectory:         ",I10,  /,& 
!               " * Nr of print steps of each trajectory:        ",I10,  / )

   END SUBROUTINE Harmonic1DModel_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE Harmonic1DModel_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord
      LOGICAL, DIMENSION(:), ALLOCATABLE :: LangevinSwitchOn

      ! Allocate memory and initialize vectors for trajectory, acceleration and masses

      CALL ERROR( BathType ==  SLAB_POTENTIAL, "Harmonic1DModel_Initialize: slab potential is not allowed for this run type " )

      IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ALLOCATE( X(1+NBath), V(1+NBath), A(1+NBath), MassVector(1+NBath), LangevinSwitchOn(1+NBath) )
         MassVector = (/ MassH, (MassBath, iCoord=1,NBath) /)

      ELSE IF ( BathType ==  DOUBLE_CHAIN ) THEN
         ALLOCATE( X(1+2*NBath), V(1+2*NBath), A(1+2*NBath), MassVector(1+2*NBath), LangevinSwitchOn(1+2*NBath) )
         MassVector = (/ MassH, (MassBath, iCoord=1,2*NBath) /)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ALLOCATE( X(1), V(1), A(1), MassVector(1), LangevinSwitchOn(1) )
         MassVector = (/ MassH /)
      END IF

      ! temporarily store the nr of dimensions
      NDim = size(X)

      ! Set variables for EOM integration in the microcanonical ensamble (system + bath)
      CALL EvolutionSetup( MolecularDynamics, NDim, MassVector, TimeStep )

      ! Set canonical dynamics at the end of the oscillator chain or chains
      IF ( BathType == CHAIN_BATH .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( NDim ) = .TRUE.
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF
      IF ( BathType == DOUBLE_CHAIN .AND. DynamicsGamma /= 0.0 ) THEN 
         LangevinSwitchOn = .FALSE.
         LangevinSwitchOn( 1+NBath ) = .TRUE.  ! end of the first chain
         LangevinSwitchOn( NDim )    = .TRUE.  ! end of the second chain
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Switch on Langevin relaxation when needed
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
      CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
      LangevinSwitchOn = (/ (.TRUE., iCoord=1,NDim ) /)
      CALL SetupThermostat( Equilibration, EquilGamma, Temperature, LangevinSwitchOn )

      DEALLOCATE( LangevinSwitchOn )

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! PRINTTYPE
      ! MINIMAL - print the ZH realizations only
      ! DEBUG - print also the average coordinates and some interesting autocorrelation functions
      ! FULL - print detailed information about the initial conditions and about the trajectories

      ! frequency spacing of the fourier transform
      dOmega =  2.*MyConsts_PI/( real(NrOfSteps) * TimeStep )

      ! Allocate memory for the averages and st dev of the coordinates
      ALLOCATE( AverageCoord(NDim), StDevCoord(NDim), AverCoordOverTrajs(NDim) )
      AverCoordOverTrajs(:) = 0.0

      IF ( PrintType >= FULL ) THEN
         ! Allocate memory for autocorrelation functions
         ALLOCATE( XatT0(NDim),QAutoCorr(0:NrOfPrintSteps), PowerSpectrum(0:NrOfPrintSteps) )
         ! Initialize variables
         QAutoCorr(:) = 0.0
         PowerSpectrum(:) = cmplx( 0.0, 0.0 )
      END IF

      ! Global temperature average
      GlobalTemperature = 0.0

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

   END SUBROUTINE Harmonic1DModel_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE Harmonic1DModel_Run()
      IMPLICIT NONE
      INTEGER :: QRealizationsUnit, TEquilUnit, QSpectralDensUnit, QSpectralAnalytic
      INTEGER :: DebugUnitEn, DebugUnitCoord, DebugUnitVel, DebugUnitEquil
      INTEGER :: iTraj, iCoord, iStep, iOmega, kStep
      REAL    :: TempAverage, TempVariance, TotEnergy, PotEnergy, KinEnergy, IstTemperature
      REAL    :: NoiseVariance, DistortedFreq, MotionVariance
      CHARACTER(100) :: OutFileName

      ! Open output file to print the brownian realizations of ZH vs time
      QRealizationsUnit = LookForFreeUnit()
      OPEN( FILE="Trajectories_Q.dat", UNIT=QRealizationsUnit )
      WRITE(QRealizationsUnit, "(A,I6,A,/)") "# ", NrTrajs, " realizations of brownian HO equilibrium dynamics (fs | Angstrom)"

      IF ( PrintType >= FULL ) THEN
         ! Open output file to print the numerical spectral density of the Q autocorrelation function
         QSpectralDensUnit = LookForFreeUnit()
         OPEN( FILE="QSpectralDensity.dat", UNIT=QSpectralDensUnit )
         WRITE(QSpectralDensUnit, "(A,I6,A,/)") "# SpectralD of brownian HO autocorrelation - ",NrTrajs," trajectories (fs | au)"

         IF ( BathType == LANGEVIN_DYN ) THEN
            ! Open output file to print the analytic spectral density of the Q autocorrelation function
            QSpectralAnalytic = LookForFreeUnit()
            OPEN( FILE="QAnalyticSD.dat", UNIT=QSpectralAnalytic )
            WRITE(QSpectralAnalytic, "(A,I6,A,/)") "# Analytic SD of brownian HO autocorrelation - ",NrTrajs," trajs (fs | au)"
         ENDIF
      END IF

      IF ( PrintType == DEBUG ) THEN
         ! Open output file to print the brownian realizations of ZH vs time
         TEquilUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationTemp.dat", UNIT=TEquilUnit )
         WRITE(TEquilUnit, "(A,I6,A,/)") "# ", NrTrajs, " temperature average and variance for each equilibration (fs | K)"
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               HARMONIC BROWNIAN TEST"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      Trajectory: DO iTraj = 1,NrTrajs        !run NrTrajs number of trajectories
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions of the system
         X(1) = 0.0
         V(1) = 0.0

         ! Set initial conditions of the bath
         IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
            CALL ThermalEquilibriumBathConditions( Bath, X(5:), V(5:), Temperature, RandomNr )
         ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
            CALL ThermalEquilibriumBathConditions( DblBath(1), X(5:NBath+4), V(5:NBath+4), Temperature, RandomNr )
            CALL ThermalEquilibriumBathConditions( DblBath(2), X(NBath+5:2*NBath+4), V(NBath+5:2*NBath+4), Temperature, RandomNr )
         ELSE IF ( BathType == LANGEVIN_DYN ) THEN
               ! nothing to do
         END IF

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", Temperature / MyConsts_K2AU

         IF ( PrintType == DEBUG ) THEN
            ! Open file to store debug info about equilibration
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Equil.dat"
            DebugUnitEquil = LookForFreeUnit()
            OPEN( Unit=DebugUnitEquil, File=OutFileName )
            WRITE( DebugUnitEquil, * ) "# EQUILIBRATION: time, potential E, kinetic E, Q"
         ENDIF

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = VHarmonic( X, A )
         DO iCoord = 1,size(X)
            A(iCoord) = A(iCoord) / MassVector(iCoord)
         END DO

         ! Do an equilibration run
         EquilibrationCycle: DO iStep = 1, NrEquilibSteps

            ! PROPAGATION for ONE TIME STEP
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, VHarmonic, PotEnergy, RandomNr )

            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy(Equilibration, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

            ! store temperature average and variance
            TempAverage = TempAverage + IstTemperature
            TempVariance = TempVariance + IstTemperature**2

            ! every PrintStepInterval steps, write debug output
            IF ( PrintType == DEBUG .AND. mod(iStep,PrintStepInterval) == 0 ) THEN
               ! Equilibration debug
               IF ( mod(iStep,PrintStepInterval) == 0 )  WRITE(DebugUnitEquil,850) real(iStep)*EquilTStep/MyConsts_fs2AU, &
                        PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, X(1)
               ! Temperature profile during equilibration
               IF (iStep == 1) WRITE(TEquilUnit,"(/,/,A,I5,/)" ) "# Trajectory nr. ", iTraj
               IF ( mod(iStep,PrintStepInterval) == 0 ) WRITE(TEquilUnit,850)  real(iStep)*EquilTStep/MyConsts_fs2AU, &
                     TempAverage/iStep, sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
            END IF
            850 FORMAT( F12.5, 2F13.6 )

         END DO EquilibrationCycle

         ! add white line to the file with temperature over time of the single traj
         IF ( PrintType == DEBUG )  WRITE(TEquilUnit,"(/)" )

         ! Close file to store equilibration detailed info
         IF ( PrintType == DEBUG ) CLOSE( Unit=DebugUnitEquil )

         PRINT "(A)", " Equilibration completed! "

         ! Compute average and standard deviation
         TempAverage = TempAverage / NrEquilibSteps 
         TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
         ! output message with average values
         PRINT "(/,A,1F10.4,/,A,1F10.4,/)",  " * Average temperature (K): ", TempAverage, &
                                             " * Standard deviation (K):  ", sqrt(TempVariance)

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Compute kinetic energy and total energy
         KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
         TotEnergy = PotEnergy + KinEnergy
         IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

         ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
         WRITE(*,600)  PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, &
                        TotEnergy*MyConsts_Hartree2eV, IstTemperature

         ! Initialize coordinate average
         AverageCoord(:) = 0.0
         StDevCoord(:) = 0.0
         ! Temperature average
         TempAverage = 0.0
         TempVariance = 0.0

         ! Write the Q coordinate in output file
         WRITE(QRealizationsUnit,"(F14.8,F14.8)") 0.0, X(1)*MyConsts_Bohr2Ang

         IF ( PrintType >= FULL ) THEN
            ! store initial coordinate of the trajectory 
            XatT0(:) = X(:)
            ! First value of the correlation functions
            QAutoCorr(0)   = QAutoCorr(0)   + XatT0(1) * X(1)
         ENDIF

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
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | EKin, EPot, ETot /Eh, IstTemp / K "
            WRITE(DebugUnitEn,800) 0.0, KinEnergy, PotEnergy, TotEnergy, IstTemperature

            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, (/ (X(iCoord)+0.05*(iCoord-1), iCoord = 1, size(X)) /)

            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)
          ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kStep = 0

         PRINT "(/,A)", " Propagating the equilibrium H-Graphene in time... "

         ! cycle over nstep velocity verlet iterations
         Propagation: DO iStep = 1, NrOfSteps

            ! Propagate for one timestep
            CALL EOM_LangevinSecondOrder( MolecularDynamics, X, V, A, VHarmonic, PotEnergy, RandomNr )

            ! Compute kin energy and temperature
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*size(X))

            ! Store variables for averages
            AverageCoord(:) = AverageCoord(:) + X(:)
            StDevCoord(:)   = StDevCoord(:)   + X(:)**2
            TempAverage  = TempAverage  + IstTemperature                     ! Temperature
            TempVariance = TempVariance + IstTemperature**2

            ! Increment global temperature average
            GlobalTemperature = GlobalTemperature + IstTemperature

            ! output to write every nprint steps 
            IF ( mod(iStep,PrintStepInterval) == 0 ) THEN

               ! increment counter for printing steps
               kStep = kStep+1
               IF ( kStep > NrOfPrintSteps ) CYCLE 

               ! Write the ZH coordinate in output file
               WRITE(QRealizationsUnit,"(F14.8,F14.8)") TimeStep*real(iStep)/MyConsts_fs2AU, X(1)*MyConsts_Bohr2Ang

               ! autocorrelation functions are computed only for a full output
               IF ( PrintType >= FULL ) THEN
                  QAutoCorr(kStep)   = QAutoCorr(kStep)   + XatT0(1) * X(1)
               ENDIF

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                        KinEnergy, PotEnergy, TotEnergy, IstTemperature
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                          (/ (X(iCoord)+0.05*(iCoord-1), iCoord = 1, size(X)) /)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
               END IF

            END IF 
 
         END DO Propagation

         ! Normalize averages of this trajectory
         AverageCoord(:) = AverageCoord(:) / NrOfSteps 
         StDevCoord(:) = SQRT( (StDevCoord(:) / NrOfSteps) -  AverageCoord(:)**2 )
         TempAverage = TempAverage / NrOfSteps 
         TempVariance = SQRT( (TempVariance / NrOfSteps) - TempAverage**2 )

         ! Accumulate average of the coordinates over the trajectories   
         AverCoordOverTrajs(:) = AverCoordOverTrajs(:) + AverageCoord(:)

         ! Add while line to the realizations output to separate each traj
         WRITE(QRealizationsUnit,*)  " "

         PRINT "(A)", " Time propagation completed! "

         ! print the average values of that trajectory 
         WRITE(*,700)  TempAverage, TempVariance, AverageCoord(1)* MyConsts_Bohr2Ang, StDevCoord(1)* MyConsts_Bohr2Ang

         IF ( PrintType == DEBUG ) THEN
               CLOSE( Unit=DebugUnitEn )
               CLOSE( Unit=DebugUnitCoord )
               CLOSE( Unit=DebugUnitVel )
         ENDIF

      END DO Trajectory

      PRINT "(A)"," Done! "                                    

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize global average of the temperature
      GlobalTemperature = GlobalTemperature / real( NrOfSteps*NrTrajs )
      ! Normalize average coordinates
      AverCoordOverTrajs(:) = AverCoordOverTrajs(:) / real( NrTrajs )

      ! Print temperature average
      WRITE(*,701)  GlobalTemperature

      ! If full print, compute and write spectral densities of the autocorrelation functions
      IF ( PrintType >= FULL ) THEN
     
         ! Q SPECTRAL DENSITY
         ! normalize autocorrelation function
         PowerSpectrum(:) = QAutoCorr(:) / real(NrTrajs ) - ( AverCoordOverTrajs(1) )**2
         ! compute fourier transform
         CALL DiscreteFourier( PowerSpectrum )
         ! Print spectral density of the stocastic process X(t)
         DO iOmega = 0, NrOfPrintSteps
               WRITE(QSpectralDensUnit,"(F20.8,3F20.8)")  iOmega*dOmega/MyConsts_cmmin1toAU,  abs(PowerSpectrum(iOmega)), & 
                                                  real(PowerSpectrum(iOmega)), aimag( PowerSpectrum(iOmega))
         END DO

         IF ( BathType == LANGEVIN_DYN ) THEN
            ! For comparison, compute the spectral density from the analytic results of the Langevin HO
            ! Compute the variables needed for the analytical autocorrelation function
            NoiseVariance = 2.0 * Temperature * MassH * DynamicsGamma
            DistortedFreq = SQRT( QFreq**2 - DynamicsGamma**2 / 4.0 )
            MotionVariance = NoiseVariance / ( 2.0 * MassH**2 * DynamicsGamma * QFreq**2 )

            ! Print spectral density of the stocastic process X(t)
            DO iOmega = 0, NrOfPrintSteps
                  WRITE(QSpectralAnalytic,"(1F12.6,3F24.18)")  iOmega*dOmega/MyConsts_cmmin1toAU, &
           NoiseVariance/(2.0*MassH**2*MyConsts_PI) / ( (QFreq**2 - (iOmega*dOmega)**2)**2 + (DynamicsGamma*iOmega*dOmega)**2  )
            END DO
         END IF

      END IF

      ! Close output files
      CLOSE( QRealizationsUnit )
      IF ( PrintType == DEBUG ) CLOSE( TEquilUnit )
      IF ( PrintType >= FULL ) CLOSE( QSpectralDensUnit )
      IF ( BathType == LANGEVIN_DYN ) CLOSE( QSpectralAnalytic )

   600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                  " * Potential Energy (eV)        ",1F10.4,/    &
                  " * Kinetic Energy (eV)          ",1F10.4,/    &
                  " * Total Energy (eV)            ",1F10.4,/    &
                  " * Istantaneous temperature (K) ",1F10.4,/ ) 

   700 FORMAT (/, " * Average temperature (K)        ",1F10.4,/    &
                  "   Standard deviation (K)         ",1F10.4,/    &
                  " * Average Q height (Ang)         ",1F10.4,/    &
                  "   Standard deviation (Ang)       ",1F10.4,/ ) 

   701 FORMAT (/, " Average temperature for all trajs and time steps ",/   &
                  " * Average temperature (K)        ",1F10.4,/ ) 

   800 FORMAT(F12.5,1000F15.8)

   END SUBROUTINE Harmonic1DModel_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE Harmonic1DModel_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( AverageCoord, StDevCoord )
      IF ( PrintType >= FULL )  DEALLOCATE( XatT0, QAutoCorr, AverCoordOverTrajs, PowerSpectrum )
      DEALLOCATE( X, V, A, MassVector )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE Harmonic1DModel_Dispose

!*************************************************************************************************

!*******************************************************************************
!> test harmonic potential (input and output in AU).
!>
!> @param Positions    Array with 1 cartesian coordinates for the brownian HO
!> @param Forces       Output array with the derivatives of the potential 
!> @param vv           output potential
!*******************************************************************************     
   REAL FUNCTION VHarmonic( Positions, Forces )
      IMPLICIT NONE
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 
      INTEGER :: NrDOF, i       
      
      ! Check the number of degree of freedom
      NrDOF = size( Positions )
      CALL ERROR( size(Forces) /= NrDOF, "PotentialModule.VHarmonic: array dimension mismatch" )
      CALL ERROR( (NrDOF /= NDim), "PotentialModule.VHarmonic: wrong number of DoFs" )
      
      ! Initialize forces and potential
      VHarmonic = 0.0
      Forces(:) = 0.0

      ! Compute potential and forces for the system potential
      VHarmonic = VHarmonic + 0.5 * ForceConstant * Positions(1)**2
      Forces(1) = -  ForceConstant * Positions(1)

      IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, Positions(1), Positions(2:), VHarmonic, Forces(1), Forces(2:) ) 

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), Positions(1), Positions(2:NBath+1), VHarmonic, Forces(1), Forces(2:NBath+1) ) 
         CALL BathPotentialAndForces( DblBath(2), Positions(1), Positions(NBath+2:2*NBath+1), VHarmonic, &
                                                                                          Forces(1), Forces(NBath+2:2*NBath+1) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! nothing to do
      END IF
      
   END FUNCTION VHarmonic

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

END MODULE Harmonic1DModel
