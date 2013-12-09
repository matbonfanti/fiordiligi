!***************************************************************************************
!*                              MODULE ThermalEquilibrium
!***************************************************************************************
!
!>  \brief     Subroutines for an equilibrium simulation at a given temperature
!>  \details   This module contains subroutines to compute averages  \n
!>             in an equilibrium simulation at a given temperature   \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             8 October 2013
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 23 October 2013: double chain bath has been implemented
!>  \arg 8 Novembre 2013: implemented the log of the input variables
!                
!>  \todo     Fix FT of the autocorrelation functions (they should be cosine transform)
!>  \todo     Implement an option to print the xyz trajectory to visualize the traj
!
!***************************************************************************************
MODULE ThermalEquilibrium
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

   PUBLIC :: ThermalEquilibrium_ReadInput, ThermalEquilibrium_Initialize, ThermalEquilibrium_Run, ThermalEquilibrium_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

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
   REAL, DIMENSION(:), ALLOCATABLE  ::  ZHAutoCorr         !< Autocorrelation function of ZH: < ZH(0) ZH(t) >
   REAL, DIMENSION(:), ALLOCATABLE  ::  CHAutoCorr         !< Autocorrelation function of CH: < (ZH-ZC)(0) (ZH-ZC)(t) >
   REAL, DIMENSION(:), ALLOCATABLE  ::  FullAutoCorr       !< Full autocorrelation function: < Xvec(0) dot Xvec(t) >
   REAL, DIMENSION(:), ALLOCATABLE  ::  AverageCoord       !< Average coordinates of each trajectory
   REAL, DIMENSION(:), ALLOCATABLE  ::  AverCoordOverTrajs !< Average coordinates of all the trajectories
   REAL, DIMENSION(:), ALLOCATABLE  ::  StDevCoord         !< Standard deviation of the coordinates of each trajectory
   COMPLEX, DIMENSION(:), ALLOCATABLE  :: PowerSpectrum    !< Complex array to compute the power spectrum of autocorr functions

   TYPE(RNGInternalState) :: RandomNr

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific for an equilibrium
!> simulation at given T.
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE ThermalEquilibrium_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

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


      WRITE(*, 905) Temperature/MyConsts_K2AU, EquilTStep/MyConsts_fs2AU,  &
                    EquilGamma*MyConsts_fs2AU, NrEquilibSteps

      WRITE(*, 904) NrTrajs, TimeStep/MyConsts_fs2AU, NrOfSteps, NrOfPrintSteps


   904 FORMAT(" * Dynamical simulation variables               ", /,&
              " * Nr of trajectories:                          ",I10,  /,& 
              " * Propagation time step (fs):                  ",F10.4,/,& 
              " * Nr of time steps of each trajectory:         ",I10,  /,& 
              " * Nr of print steps of each trajectory:        ",I10,  / )

   905 FORMAT(" * Equilibration variables                      ", /,&
              " * Temperature of the system+surface:           ",F10.4,/,&
              " * Equilibration time step (fs):                ",F10.4,/,& 
              " * Gamma of the Langevin force (fs^-1):         ",F10.4,/,& 
              " * Nr of equilibration steps:                   ",I10,  / )

   END SUBROUTINE ThermalEquilibrium_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the vibrational relaxation simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE ThermalEquilibrium_Initialize()
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

      ! temporarily store the nr of dimensions
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

      ! Switch on Langevin relaxation when needed
      IF ( BathType == LANGEVIN_DYN .AND. DynamicsGamma /= 0.0 ) THEN
         LangevinSwitchOn = (/ .FALSE. , .FALSE. , .FALSE. , .TRUE. /)
         CALL SetupThermostat( MolecularDynamics, DynamicsGamma, Temperature, LangevinSwitchOn )
      END IF

      ! Set variables for EOM integration with Langevin thermostat, during initial equilibration
      CALL EvolutionSetup( Equilibration, NDim, MassVector, EquilTStep )
      LangevinSwitchOn = (/ (.FALSE., iCoord=1,3), (.TRUE., iCoord=4,NDim ) /)
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
      ALLOCATE( AverageCoord(NDim), StDevCoord(NDim) )

      IF ( PrintType >= FULL ) THEN
         ! Allocate memory for autocorrelation functions
         ALLOCATE( XatT0(NDim),ZHAutoCorr(0:NrOfPrintSteps), CHAutoCorr(0:NrOfPrintSteps),   &
                    AverCoordOverTrajs(NDim), FullAutoCorr(0:NrOfPrintSteps), PowerSpectrum(0:NrOfPrintSteps) )
         ! Initialize variables
         ZHAutoCorr(:) = 0.0
         CHAutoCorr(:) = 0.0
         FullAutoCorr(:) = 0.0
         PowerSpectrum(:) = cmplx( 0.0, 0.0 )
         AverCoordOverTrajs(:) = 0.0
      END IF

      ! Global temperature average
      GlobalTemperature = 0.0

      ! Initialize random number seed
      CALL SetSeed( RandomNr, -1 )

   END SUBROUTINE ThermalEquilibrium_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE ThermalEquilibrium_Run()
      IMPLICIT NONE
      INTEGER :: ZHRealizationsUnit, TEquilUnit
      INTEGER :: ZHSpectralDensUnit, CHSpectralDensUnit, FullSpectralDensUnit
      INTEGER :: DebugUnitEn, DebugUnitCoord, DebugUnitVel, DebugUnitEquil
      INTEGER :: iTraj, iCoord, iStep, iOmega, kStep
      REAL    :: TempAverage, TempVariance, TotEnergy, PotEnergy, KinEnergy, IstTemperature
      CHARACTER(100) :: OutFileName

      ! Open output file to print the brownian realizations of ZH vs time
      ZHRealizationsUnit = LookForFreeUnit()
      OPEN( FILE="Trajectories_Zh.dat", UNIT=ZHRealizationsUnit )
      WRITE(ZHRealizationsUnit, "(A,I6,A,/)") "# ", NrTrajs, " realizations of zH equilibrium dynamics (fs | Angstrom)"

      IF ( PrintType >= FULL ) THEN
         ! Open output file to print the spectral density of the ZH autocorrelation function
         ZHSpectralDensUnit = LookForFreeUnit()
         OPEN( FILE="ZHSpectralDensity.dat", UNIT=ZHSpectralDensUnit )
         WRITE(ZHSpectralDensUnit, "(A,I6,A,/)") "# Spectral density of ZH autocorrelation - ", NrTrajs, " trajectories (fs | au)"

         ! Open output file to print the spectral density of the C-H autocorrelation function
         CHSpectralDensUnit = LookForFreeUnit()
         OPEN( FILE="CHSpectralDensity.dat", UNIT=CHSpectralDensUnit )
         WRITE(CHSpectralDensUnit, "(A,I6,A,/)") "# Spectral density of C-H autocorrelation - ", NrTrajs, " trajectories (fs | au)"

         ! Open output file to print the spectral density of the full autocorrelation function
         FullSpectralDensUnit = LookForFreeUnit()
         OPEN( FILE="FullSpectralDensity.dat", UNIT=FullSpectralDensUnit )
         WRITE(FullSpectralDensUnit, "(A,I6,A,/)") "# Spectral density of the full autocorr - ", NrTrajs, " trajectories (fs | au)"
      END IF

      IF ( PrintType == DEBUG ) THEN
         ! Open output file to print the brownian realizations of ZH vs time
         TEquilUnit = LookForFreeUnit()
         OPEN( FILE="EquilibrationTemp.dat", UNIT=TEquilUnit )
         WRITE(TEquilUnit, "(A,I6,A,/)") "# ", NrTrajs, " temperature average and variance for each equilibration (fs | K)"
      END IF

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               EQUILIBRIUM SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", NrTrajs, " trajectories ... "

      Trajectory: DO iTraj = 1,NrTrajs        !run NrTrajs number of trajectories
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Equilibrium position of the system H and C atom
         X(1:2) = 0.0000
         X(3) = 1.483 / MyConsts_Bohr2Ang
         X(4) = C1Puckering
         V(1:4) = 0.0

         ! Set initial conditions of the bath
         IF ( BathType == SLAB_POTENTIAL ) THEN 
            CALL ThermalEquilibriumConditions( X, V, Temperature, MassH, MassC, RandomNr )
         ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
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
            WRITE( DebugUnitEquil, * ) "# EQUILIBRATION: time, potential E, kinetic E, X_H, Y_H, Z_H, Z_C1"
         ENDIF

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = ThermalEquilibriumPotential( X, A )
         A(:) = A(:) / MassVector(:)

         ! Do an equilibration run
         EquilibrationCycle: DO iStep = 1, NrEquilibSteps

            ! PROPAGATION for ONE TIME STEP
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, ThermalEquilibriumPotential, PotEnergy, RandomNr )

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
                        PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, X(1:4)
               ! Temperature profile during equilibration
               IF (iStep == 1) WRITE(TEquilUnit,"(/,/,A,I5,/)" ) "# Trajectory nr. ", iTraj
               IF ( mod(iStep,PrintStepInterval) == 0 ) WRITE(TEquilUnit,851)  real(iStep)*EquilTStep/MyConsts_fs2AU, &
                     TempAverage/iStep, sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
            END IF
            850 FORMAT( F12.5, 6F13.6 )
            851 FORMAT( F12.5, 2F13.6 )

         END DO EquilibrationCycle

         ! Close file to store equilibration detailed info
         IF ( PrintType == DEBUG )  CLOSE( Unit=DebugUnitEquil )

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

         ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
         IF ( BathType == SLAB_POTENTIAL )    X(3:size(X)) = X(3:size(X))  - (X(5)+X(6)+X(7))/3.0

         ! Write the ZH coordinate in output file
         WRITE(ZHRealizationsUnit,"(F14.8,F14.8)") 0.0, X(3)*MyConsts_Bohr2Ang

         IF ( PrintType >= FULL ) THEN
            ! store initial coordinate of the trajectory 
            XatT0(:) = X(:)
            ! First value of the correlation functions
            FullAutoCorr(0) = FullAutoCorr(0) + TheOneWithVectorDotVector( XatT0(:), X(:) )
            ZHAutoCorr(0)   = ZHAutoCorr(0)   + XatT0(3) * X(3)
            CHAutoCorr(0)   = CHAutoCorr(0)   + (XatT0(3)-XatT0(4)) * (X(3)-X(4))
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
            WRITE(DebugUnitCoord,800) 0.0, X(1), X(2), X(3), X(4), (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, size(X)-4) /)

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
            CALL EOM_LangevinSecondOrder( Equilibration, X, V, A, ThermalEquilibriumPotential, PotEnergy, RandomNr )

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

               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               IF ( BathType == SLAB_POTENTIAL )    X(3:size(X)) = X(3:size(X))  - (X(5)+X(6)+X(7))/3.0

               ! Write the ZH coordinate in output file
               WRITE(ZHRealizationsUnit,"(F14.8,F14.8)") TimeStep*real(iStep)/MyConsts_fs2AU, X(3)*MyConsts_Bohr2Ang

               ! autocorrelation functions are computed only for a full output
               IF ( PrintType >= FULL ) THEN
                  FullAutoCorr(kStep) = FullAutoCorr(kStep) + TheOneWithVectorDotVector( XatT0(:), X(:) )
                  ZHAutoCorr(kStep)   = ZHAutoCorr(kStep)   + XatT0(3) * X(3)
                  CHAutoCorr(kStep)   = CHAutoCorr(kStep)   + (XatT0(3)-XatT0(4)) * (X(3)-X(4))
               ENDIF

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) TimeStep*real(iStep)/MyConsts_fs2AU, &
                                        KinEnergy, PotEnergy, TotEnergy, IstTemperature
                  WRITE(DebugUnitCoord,800) TimeStep*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                          (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, size(X)-4) /)
                  WRITE(DebugUnitVel,800) TimeStep*real(iStep)/MyConsts_fs2AU, V(:)
               END IF

!                ! Store the trajectory for XYZ printing
!                IF ( PrintType >= FULL  .AND. RunType == EQUILIBRIUM ) THEN
!                      Trajectory( :, kstep ) = 0.0
!                      Trajectory( 1:min(16,3+nevo) , kstep ) = X( 1:min(16,3+nevo) ) 
!                      NrOfTrajSteps = kstep
!                END IF

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
         WRITE(ZHRealizationsUnit,*)  " "

         PRINT "(A)", " Time propagation completed! "

         ! print the average values of that trajectory 
         WRITE(*,700)  TempAverage, TempVariance, AverageCoord(1:3)* MyConsts_Bohr2Ang, StDevCoord(1:3)* MyConsts_Bohr2Ang, &
                     AverageCoord(4)* MyConsts_Bohr2Ang, StDevCoord(4)* MyConsts_Bohr2Ang

!             IF  ( RunType == EQUILIBRIUM ) THEN
!                ! write the xyz file of the trajectory, if requested
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",iTraj,".xyz"
!                CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
!                                              GraphiteLatticeConstant()*MyConsts_Bohr2Ang )
!             ENDIF
! 
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
         
         ! ZH SPECTRAL DENSITY
         ! normalize autocorrelation function
         PowerSpectrum(:) = ZHAutoCorr(:) / real(NrTrajs ) - ( AverCoordOverTrajs(3) )**2
         ! compute fourier transform
         CALL DiscreteFourier( PowerSpectrum )
         ! Print spectral density of the stocastic process X(t)
         DO iOmega = 0, NrOfPrintSteps
               WRITE(ZHSpectralDensUnit,"(F20.8,3F20.8)")  iOmega*dOmega/MyConsts_cmmin1toAU,  abs(PowerSpectrum(iOmega)), & 
                                                  real(PowerSpectrum(iOmega)), aimag( PowerSpectrum(iOmega))
         END DO

         ! CH SPECTRAL DENSITY
         ! normalize autocorrelation function
         PowerSpectrum(:) = CHAutoCorr(:) / real(NrTrajs ) - ( AverCoordOverTrajs(3)-AverCoordOverTrajs(4) )**2
         ! compute fourier transform
         CALL DiscreteFourier( PowerSpectrum )
         ! Print spectral density of the stocastic process X(t)
         DO iOmega = 0, NrOfPrintSteps
               WRITE(CHSpectralDensUnit,"(F20.8,3F20.8)")  iOmega*dOmega/MyConsts_cmmin1toAU,  abs(PowerSpectrum(iOmega)), & 
                                                  real(PowerSpectrum(iOmega)), aimag( PowerSpectrum(iOmega))
         END DO

         ! FULL SPECTRAL DENSITY
         ! normalize autocorrelation function
         PowerSpectrum(:) = FullAutoCorr(:) / real(NrTrajs )  - DOT_PRODUCT(AverCoordOverTrajs, AverCoordOverTrajs )
         ! compute fourier transform
         CALL DiscreteFourier( PowerSpectrum )
         ! Print spectral density of the stocastic process X(t)
         DO iOmega = 0, NrOfPrintSteps
               WRITE(FullSpectralDensUnit,"(F20.8,3F20.8)")  iOmega*dOmega/MyConsts_cmmin1toAU,  abs(PowerSpectrum(iOmega)), & 
                                                  real(PowerSpectrum(iOmega)), aimag( PowerSpectrum(iOmega))
         END DO
      END IF

      ! Close output files
      CLOSE( ZHRealizationsUnit )
      IF ( PrintType == DEBUG ) CLOSE( TEquilUnit )
      IF ( PrintType >= FULL ) THEN 
         CLOSE( ZHSpectralDensUnit )
         CLOSE( CHSpectralDensUnit )
         CLOSE( FullSpectralDensUnit )
      END IF

   600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                  " * Potential Energy (eV)        ",1F10.4,/    &
                  " * Kinetic Energy (eV)          ",1F10.4,/    &
                  " * Total Energy (eV)            ",1F10.4,/    &
                  " * Istantaneous temperature (K) ",1F10.4,/ ) 

   700 FORMAT (/, " * Average temperature (K)        ",1F10.4,/    &
                  "   Standard deviation (K)         ",1F10.4,/    &
                  " * Average H coordinates (Ang)    ",3F10.4,/    &
                  "   Standard deviation (Ang)       ",3F10.4,/    &
                  " * Average Z height (Ang)         ",1F10.4,/    &
                  "   Standard deviation (Ang)       ",1F10.4,/ ) 

   701 FORMAT (/, " Average temperature for all trajs and time steps ",/   &
                  " * Average temperature (K)        ",1F10.4,/ ) 

   800 FORMAT(F12.5,1000F15.8)


   END SUBROUTINE ThermalEquilibrium_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the vibrational relaxation simulation.
!>
!*******************************************************************************
   SUBROUTINE ThermalEquilibrium_Dispose()
      IMPLICIT NONE

      ! Deallocate memory 
      DEALLOCATE( AverageCoord, StDevCoord )
      IF ( PrintType >= FULL )  DEALLOCATE( XatT0, ZHAutoCorr, CHAutoCorr, AverCoordOverTrajs, FullAutoCorr, PowerSpectrum )
      DEALLOCATE( X, V, A, MassVector )

      ! Unset propagators 
      CALL DisposeEvolutionData( MolecularDynamics )
      CALL DisposeEvolutionData( Equilibration )

   END SUBROUTINE ThermalEquilibrium_Dispose

!*************************************************************************************************

   REAL FUNCTION ThermalEquilibriumPotential( Positions, Forces )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i       

      ! Check the number degrees of freedom
      NrDOF = size( Positions )
      CALL ERROR( size(Forces) /= NrDOF, "ThermalEquilibrium.ThermalEquilibriumPotential: array dimension mismatch" )
      CALL ERROR( NrDOF /= NDim, "ThermalEquilibrium.ThermalEquilibriumPotential: wrong number of DoFs" )

      ! Initialize forces and potential
      ThermalEquilibriumPotential = 0.0
      Forces(:)                   = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         ThermalEquilibriumPotential = VHSticking( Positions, Forces )

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         ThermalEquilibriumPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, Positions(4)-C1Puckering, Positions(5:), ThermalEquilibriumPotential, &
                                                                               Forces(4), Forces(5:) ) 
      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Compute potential and forces of the system
         ThermalEquilibriumPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), ThermalEquilibriumPotential, &
                                                                               Forces(4), Forces(5:NBath+4) ) 
         CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), &
                                                            ThermalEquilibriumPotential, Forces(4), Forces(NBath+5:2*NBath+4) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential using only the 4D subroutine
         ThermalEquilibriumPotential = VHFourDimensional( Positions, Forces )

      END IF

   END FUNCTION ThermalEquilibriumPotential

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

END MODULE ThermalEquilibrium
