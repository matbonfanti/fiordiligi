!***************************************************************************************
!*                              PROGRAM JK6
!***************************************************************************************
!>  \mainpage      Program JK6 - version 2
!>
!>  Classical simulations of H + graphene slab of 121 carbon atoms  \n
!>  * Model: 3D for H + 1D Z coordinate for 121 carbon atom         \n
!>  * Propagation: Velocity-Verlet in the microcanonical ensamble
!>
!>  \author        Bret Jackson, Jay Kerwin, Matteo Bonfanti
!>  \version       2.0
!>  \date          10 January 2013
!>
!>  \todo          Fix temperature boundaries in the istogram of temperature
!>                 sampling (in Langevin HO test)
!>  \todo          Fix normalization of DFT and inverse DFT 
!>  \todo          Introduce more C atoms in the slab when printing the traj
!>  \todo          PERFORMANCES OF FT !!!! \n
!>                 at least, since it's a matrix vector product with fixed
!>                 matrix, on-the-fly computation of exp coeff can be avoided
!>  \todo          deallocate arrays allocated in the main program!  \n
!>
!***************************************************************************************
PROGRAM JK6
   USE CommonData
   USE InputField
   USE PotentialModule
   USE MyConsts
   USE ErrorTrap
   USE RandomNumberGenerator
   USE ClassicalEqMotion
   USE PrintTools
   USE IndependentOscillatorsModel

   IMPLICIT NONE

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.
   
   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName
   ! Derived type to handle input data
   TYPE(InputFile) :: InputData

   ! Data type with evolution information
   TYPE(Evolution) :: Equilibration, MolecularDynamics

   ! Output variables
   CHARACTER(100) :: OutFileName         ! output file name
   INTEGER  :: NWriteUnit                ! output unit number
   INTEGER  :: TrajOutputUnit            ! unit nr for the trajectory output
   INTEGER  :: TimeCorrelationUnit       ! unit nr for the correlation function in time
   INTEGER  :: SpectrDensUnit            ! unit nr for the spectral density
   INTEGER  :: InitialTDistrib           ! unit nr to print the initial temperatures distribution
   INTEGER  :: Traj2OutputUnit

   ! Distribution of the initial istantaneous temperatures
   INTEGER, DIMENSION( 200:400 ) :: InitTempDistribution

   ! Real parameters to compute the analytical autocorrelation function
   REAL :: NoiseVariance
   REAL :: OsciFreq
   REAL :: DistortedFreq
   REAL :: MotionVariance
   REAL :: Time

   ! Autocorrelation computation
   REAL  :: ZHinTime
   ! Frequency spacing of the fourier transform
   REAL :: dOmega
   ! Store in time and in frequency the single realization
   COMPLEX, DIMENSION(:), ALLOCATABLE :: DeltaZ
   REAL, DIMENSION(:), ALLOCATABLE :: AverageDeltaZ
   COMPLEX, DIMENSION(:), ALLOCATABLE :: Analytic

   REAL, DIMENSION(501)      :: crscn

   INTEGER :: iCoord, iTraj, iStep, jRho, iCarbon, iOmega
   INTEGER :: kk, kstep

   REAL :: DynamicsGamma
   LOGICAL, DIMENSION(121) :: LangevinSwitchOn

   ! Print potential cuts in VTK format
   TYPE(VTKInfo) :: PotentialCH, PotentialCH_Model

   ! Parameters to define the grid for potential plots
   REAL :: ZHmin, ZHmax
   REAL :: ZCmin, ZCmax
   INTEGER :: NpointZH
   INTEGER :: NpointZC

   ! Parameters to define CH model potential
   REAL :: ZHFrequency, ZCFrequency, BilinearCoupl
   INTEGER :: FreezeGraphene

   ! Parameters for the model thermal bath
   REAL :: BathCutOffFreq                ! cutoff frequency of the bath
   CHARACTER(100) :: SpectralDensityFile         ! spectral density file name
   
   REAL, DIMENSION(:), ALLOCATABLE :: MassVector

   REAL, ALLOCATABLE :: XTEMP(:)
   INTEGER :: i

   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                                JK6_v2            ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                    Authors: Bret Jackson, Jay Kerwin  ')"
   PRINT "(       '         Translation in FORTRAN 90/95 by Matteo Bonfanti (January 2013) ',/)"


   !*************************************************************
   !         COMMAND LINE ARGUMENT
   !*************************************************************

   ! Check and read from command line the input file name
   NArgs = COMMAND_ARGUMENT_COUNT()
   IF (NArgs<1) THEN
      Help = .TRUE.
   ELSE
      CALL GET_COMMAND_ARGUMENT( 1, InputFileName )
      IF ( trim(InputFileName) == "help" ) Help = .TRUE.
   ENDIF
   IF (Help) THEN ! Call help
      PRINT*, ' Launch this program as:'
      PRINT*, ' % JK6 "InputFileName" '
      STOP
   ENDIF


   !*************************************************************
   !                 INPUT FILE 
   !*************************************************************

   !       ***** NOTE: units of input data ******
   !            energy      - electronVolt
   !            time        - femtosecond
   !            distance    - Angstrom
   !            temperature - Kelvin
   !            frequency   - cm-1
   !

   ! Open and read from input file the input parameters of the calculation
   CALL OpenFile( InputData, InputFileName )

   ! Check if a scattering simulation is required
   CALL SetFieldFromInput( InputData, "RunType", RunType, 1 )
   ! Set the level of output
   CALL SetFieldFromInput( InputData, "PrintType", PrintType )

   ! Read parameters common to all calculations and transform them in atomic units
   CALL SetFieldFromInput( InputData, "TimeStep", dt )
   dt = dt * MyConsts_fs2AU
   CALL SetFieldFromInput( InputData, "temp", temp )
   temp = temp * MyConsts_K2AU
   CALL SetFieldFromInput( InputData, "nstep", nstep )
   CALL SetFieldFromInput( InputData, "nprint", nprint )
   CALL SetFieldFromInput( InputData, "inum", inum )
   CALL SetFieldFromInput( InputData, "nevo", nevo )
   CALL SetFieldFromInput( InputData, "rmh", rmh )
   rmh = rmh * MyConsts_Uma2Au
   CALL SetFieldFromInput( InputData, "rmc", rmc )
   rmc = rmc * MyConsts_Uma2Au

   ! If scattering simulation, read all the relevant parameters
   IF (RunType == SCATTERING) THEN
      CALL SetFieldFromInput( InputData, "ezh", ezh )
      ezh = ezh / MyConsts_Hartree2eV
      CALL SetFieldFromInput( InputData, "zhi", zhi )
      zhi = zhi / MyConsts_Bohr2Ang
      CALL SetFieldFromInput( InputData, "irho", irho )
      CALL SetFieldFromInput( InputData, "delrho", delrho )
      delrho = delrho / MyConsts_Bohr2Ang

   ELSE IF (RunType == EQUILIBRIUM) THEN
      CALL SetFieldFromInput( InputData, "Gamma", Gamma)
      Gamma = Gamma / MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(5.0*(1.0/Gamma)/dt) )
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, dt/MyConsts_fs2AU )
      EquilTStep = EquilTStep * MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "DynamicsGamma",  DynamicsGamma, 0.0 ) 
      DynamicsGamma = DynamicsGamma / MyConsts_fs2AU

   ELSE IF ( ( RunType == OSCIBATH_EQUIL ) .OR. ( RunType == CHAINBATH_EQUIL ) ) THEN
      ! I AM NOT SURE THAT EQUILIBRATION IS STILL NEEDED!
      CALL SetFieldFromInput( InputData, "Gamma", Gamma)
      Gamma = Gamma / MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(5.0*(1.0/Gamma)/dt) )
      CALL SetFieldFromInput( InputData, "EquilTStep",  EquilTStep, dt/MyConsts_fs2AU )
      EquilTStep = EquilTStep * MyConsts_fs2AU
      ! Read cutoff frequency of the bath, only if normal oscillators bath
      BathCutOffFreq = 0.0
      IF ( RunType == OSCIBATH_EQUIL ) CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq )
      BathCutOffFreq = BathCutOffFreq * MyConsts_cmmin1toAU
      ! Read file with spectral density / normal modes freq and couplings
      CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
      ! No langevin oscillators in the thermal bath
      DynamicsGamma = 0.0
      
   ELSE IF (RunType == HARMONICMODEL) THEN
      CALL SetFieldFromInput( InputData, "Gamma", Gamma)
      Gamma = Gamma / MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(10.0*(1.0/Gamma)/dt) )

   ELSE IF (RunType == POTENTIALPRINT) THEN
      CALL SetFieldFromInput( InputData, "ZHmin", ZHmin, 0.8 )
      ZHmin = ZHmin / MyConsts_Bohr2Ang
      CALL SetFieldFromInput( InputData, "ZHmax", ZHmax, 4.0 )
      ZHmax = ZHmax / MyConsts_Bohr2Ang
      CALL SetFieldFromInput( InputData, "ZCmin", ZCmin, -0.5 )
      ZCmin = ZCmin / MyConsts_Bohr2Ang
      CALL SetFieldFromInput( InputData, "ZCmax", ZCmax, 1.0 )
      ZCmax = ZCmax / MyConsts_Bohr2Ang
      CALL SetFieldFromInput( InputData, "NpointZH", NpointZH, 400 )
      CALL SetFieldFromInput( InputData, "NpointZC", NpointZC, 100 )

      CALL SetFieldFromInput( InputData, "ZHFrequency", ZHFrequency )
      ZHFrequency = ZHFrequency * MyConsts_cmmin1toAU
      CALL SetFieldFromInput( InputData, "ZCFrequency", ZCFrequency )
      ZCFrequency = ZCFrequency * MyConsts_cmmin1toAU
      CALL SetFieldFromInput( InputData, "BilinearCoupl", BilinearCoupl )
      CALL SetFieldFromInput( InputData, "FreezeGraphene", FreezeGraphene, 1 )

   END IF

   ! close input file
   CALL CloseFile( InputData )


   !*************************************************************
   !       ALLOCATION, GENERAL INITIALIZATIONS
   !*************************************************************

   ! number of analysis step
   ntime=int(nstep/nprint)

   ! In case of a thermal bath, nevo is the number of bath dofs, so we have to
   ! add 1 to include also the C1 degree of freedom  
   IF (( RunType == OSCIBATH_EQUIL ) .OR. ( RunType == CHAINBATH_EQUIL ))  &
            nevo =  nevo + 1 
   
   IF ( (RunType == SCATTERING) .OR. (RunType == EQUILIBRIUM) .OR.   &
                               ( RunType == OSCIBATH_EQUIL ) .OR. ( RunType == CHAINBATH_EQUIL ) ) THEN


         ! Allocation of position, velocity, acceleration arrays
         ALLOCATE( X(nevo+3), V(nevo+3), A(nevo+3), APre(nevo+3) )
         
         ! if XYZ files of the trajectories are required, allocate memory to store the traj
         IF ( PrintType >= FULL  .AND. ( RunType == SCATTERING .OR. RunType == EQUILIBRIUM ) ) THEN
               ALLOCATE( Trajectory( 16, ntime ) )
         END IF

         ! Allocate and define masses
         ALLOCATE( MassVector( nevo + 3 ) )
         IF ( RunType == CHAINBATH_EQUIL ) THEN
               MassVector = (/ (rmh, iCoord=1,3), rmc, (1.0, iCoord=1,nevo-1) /)
         ELSE
               MassVector = (/ (rmh, iCoord=1,3), (rmc, iCoord=1,nevo) /)
         ENDIF

         ! Set variables for EOM integration in the microcanonical ensamble
         CALL EvolutionSetup( MolecularDynamics, nevo+3, MassVector, dt )

         IF ( DynamicsGamma /= 0.0 .AND. RunType /= OSCIBATH_EQUIL .AND.  RunType /= CHAINBATH_EQUIL  ) THEN 
            PRINT "(/,A,/)", " Setting langevin atoms at the border of the slab "
            ! Set canonical dynamics at the borders of the carbon slab
            LangevinSwitchOn = .TRUE.
            LangevinSwitchOn( 1: MIN( 73, nevo )+3 ) = .FALSE.
            CALL SetupThermostat( MolecularDynamics, DynamicsGamma, temp, LangevinSwitchOn(1:nevo+3) )
         END IF

         ! Set variables for EOM integration with Langevin thermostat
         CALL EvolutionSetup( Equilibration, nevo+3, MassVector, EquilTStep )
         CALL SetupThermostat( Equilibration, Gamma, temp, (/ (.FALSE., iCoord=1,4), (.TRUE., iCoord=1,nevo-1) /) )

   ELSE IF (RunType == HARMONICMODEL) THEN 

         ! Allocation of position, velocity, acceleration arrays
         ALLOCATE( X(nevo), V(nevo), A(nevo), APre(nevo) )

         ! Set variables for EOM integration with Langevin thermostat
         CALL EvolutionSetup( Equilibration, nevo, (/ (rmh, iCoord=1,nevo) /), dt )
         CALL SetupThermostat( Equilibration, Gamma, temp, (/ (.TRUE., iCoord=1,nevo) /) )

   ELSE IF (RunType == POTENTIALPRINT) THEN 

         ! Allocation of position, velocity, acceleration arrays
         ALLOCATE( X(nevo+3), A(nevo+3) )

   END IF

   !*************************************************************
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   CALL SetupPotential( .FALSE. )
   
   ! If needed setup bath frequencies and coupling for IO Model in normal form
   IF (  RunType == OSCIBATH_EQUIL ) THEN
      CALL SetupIndepOscillatorsModel( nevo-1, STANDARD_BATH, SpectralDensityFile, rmc, BathCutOffFreq )
   END IF
   IF (  RunType == CHAINBATH_EQUIL ) THEN
      CALL SetupIndepOscillatorsModel( nevo-1, CHAIN_BATH, SpectralDensityFile, rmc )
   END IF
   
   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! THIS LOG PRINTING PART SHOULD BE EXTENDED AND ORDERED !!!!!!!!!!!!!1
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! write to standard output the parameters of the calculation

   WRITE(*,898) rmh / MyConsts_Uma2Au, rmc / MyConsts_Uma2Au, temp / MyConsts_K2AU
   WRITE(*,899) dt / MyConsts_fs2AU, nstep, ntime
   IF (RunType == SCATTERING) WRITE(*,900) ezh * MyConsts_Hartree2eV, zhi * MyConsts_Bohr2Ang
   WRITE(*,901) nevo, inum
   IF (RunType == SCATTERING) WRITE(*,902) irho, delrho * MyConsts_Bohr2Ang
   IF (RunType == EQUILIBRIUM) WRITE(*,903) Gamma / MyConsts_Uma2Au * MyConsts_fs2AU,  EquilTStep/MyConsts_fs2AU,   &
                 EquilTStep*real(NrEquilibSteps)/MyConsts_fs2AU
   WRITE(*,"(/)")

   898 FORMAT(" * Mass of the H atom (UMA):                    ",F10.4,/, &
              " * Mass of the C atoms (UMA):                   ",F10.4,/, &
              " * Temperature of the system (Kelvin):          ",F10.4      )
   899 FORMAT(" * Time step of the evolution (fs):             ",F10.4,/, &
              " * Total nr of time step per trajectory:        ",I8   ,/, &
              " * Nr of analysis steps:                        ",I5         )
   900 FORMAT(" * Initial kinetic energy of the H atom (eV):   ",F10.4,/, &
              " * Initial Z position of the H atom (Ang):      ",F10.4      )
   901 FORMAT(" * Nr of evolving C atoms:                      ",I6,/, &
              " * Nr of trajectories (per rho value):          ",I4         )
   902 FORMAT(" * Nr of rho values to sample:                  ",I6,/, &
              " * Grid spacing in rho (Ang):                   ",F10.4      )
   903 FORMAT(" * Langevin friction constant (1/fs)            ",F10.4,/, &
              " * Equilibration time step (fs)                 ",F10.4,/, &
              " * Equilibration total time (fs)                ",F10.4      )

   !****************************************************************************
   IF (RunType == SCATTERING) THEN
   !****************************************************************************

         CALL ScatteringSimulation()

   !****************************************************************************
   ELSE IF (RunType == EQUILIBRIUM .OR. RunType == OSCIBATH_EQUIL .OR. RunType  == CHAINBATH_EQUIL ) THEN   
   !****************************************************************************

         CALL HGraphiteEquilibrium()
      
   !****************************************************************************
   ELSE IF (RunType == HARMONICMODEL) THEN   
   !****************************************************************************

         CALL BrownianHarmonicOscillator()

   !****************************************************************************
   ELSE IF (RunType == POTENTIALPRINT) THEN   
   !****************************************************************************

         CALL PotentialCuts()

   END IF


   ! Dispose memory for evolution data
   CALL DisposeEvolutionData( MolecularDynamics )
   CALL DisposeEvolutionData( Equilibration )


      CONTAINS

!***************************************************************************************************
!                  H-GRAPHITE SCATTERING SIMULATION
!***************************************************************************************************

   SUBROUTINE ScatteringSimulation()
      IMPLICIT NONE
      REAL, DIMENSION(121) :: vzsum, zsum
      REAL, DIMENSION(:,:), ALLOCATABLE :: vz2av, zcav
      REAL :: vav, vsqav, zav, zsqav

      !*************************************************************
      !       AVERAGES VARIABLES ALLOCATION AND INITIALIZATION
      !*************************************************************

      ! Initialize average values for the INITIAL CONDITION of the lattice Z coordinates
      zav    = 0.0         ! average position
      zsqav  = 0.0         ! mean squared position
      vav    = 0.0         ! average velocity
      vsqav  = 0.0         ! mean squared velocity

      ALLOCATE( zcav(10,ntime), vz2av(10,ntime) )
      zcav(:,:)  = 0.0
      vz2av(:,:) = 0.0

      ! allocate and initialize trapping probability variable
      ALLOCATE( ptrap( irho+1, ntime ) )
      ptrap(:,:) = 0.0

      ! Initialize random generation of initial conditions for scattering simulation
      call random_seed()

      ! scan over the impact parameter... 
      DO jRho = 1,irho+1

         ! Set impact parameter and print message
         ImpactPar = float(jRho-1)*delrho+0.000001

         PRINT "(2/,A)",    "***************************************************"
         PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar
         PRINT "(A,/)" ,    "***************************************************"

         PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "

         !run inum number of trajectories at the current rxn parameter
         DO iTraj = 1,inum
         
            PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

            !*************************************************************
            ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
            !*************************************************************

            ! Set initial conditions
            CALL ScatteringConditions( X, V, ImpactPar, zhi, ezh, temp, rmh, rmc )

            !*************************************************************
            ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
            !*************************************************************

            ! increment average values for the initial lattice condition
            DO iCarbon = 4,nevo+3
               zav   = zav   + X(iCarbon)
               zsqav = zsqav + X(iCarbon)**2
               vav   = vav   + V(iCarbon)
               vsqav = vsqav + V(iCarbon)**2
            END DO

            DO iCarbon = 1,121
               zsum(iCarbon)  = 0.0
               vzsum(iCarbon) = 0.0
            END DO

            ! Compute starting potential and forces
            PotEnergy = VHSticking( X, A )
            A(1:3) = A(1:3) / rmh
            A(4:nevo+3) = A(4:nevo+3) / rmc
            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

            !*************************************************************
            !         TIME EVOLUTION OF THE TRAJECTORY
            !*************************************************************
 
            ! initialize counter for printing steps
            kstep=0

            ! cycle over nstep velocity verlet iterations
            DO iStep = 1,nstep

               ! Propagate for one timestep only if in the interaction region
               IF ( X(3) <= zhi )  CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, PotEnergy )

               ! Compute kin energy and temperature
               KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

               ! compute average values of the carbon positions
               DO iCarbon = 1,nevo
                  zsum(iCarbon)  = zsum(iCarbon)  + X(iCarbon+3)
                  vzsum(iCarbon) = vzsum(iCarbon) + V(iCarbon+3)**2
               END DO

               ! every nprint steps, compute trapping and write output
               IF ( mod(iStep,nprint) == 0 ) THEN

                  ! increment counter for printing steps
                  kstep=kstep+1

                  ! check if H is in the trapping region
                  IF ( X(3) <= AdsorpLimit ) ptrap(jRho,kstep) = ptrap(jRho,kstep)+1.0

                  ! average lattice velocity and positions
                  CALL lattice( kstep, zsum, vzsum, zcav, vz2av  )

                  ! Store the trajectory for XYZ printing
                  IF ( PrintType >= FULL ) THEN
                        Trajectory( :, kstep ) = X(1:16)
                        NrOfTrajSteps = kstep
                  END IF

               END IF 

            END DO

            ! At the end of the propagation, write the xyz file of the trajectory, if requested
            IF ( PrintType >= FULL ) THEN

               ! print the full trajectory as output file
               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",jRho,"_",iTraj,".xyz"
               CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
                                             GraphiteLatticeConstant()*MyConsts_Bohr2Ang )

            END IF

         END DO
      
         PRINT "(A)"," Done! "

      END DO

      ! pritn info about the average for the initial lattice condition
      zav=zav/float(nevo*inum*(irho+1))
      zsqav=zsqav/float(nevo*inum*(irho+1))
      vav=vav/float(nevo*inum*(irho+1))
      vsqav=vsqav/float(nevo*inum*(irho+1))

      WRITE(*,"(/,A,/)") " * Average values for the initial lattice coordinates and velocities "
      write(*,*) "average init zc=",zav * MyConsts_Bohr2Ang, " Ang"
      write(*,*) "average init Vc=",0.5*zsqav*CarbonForceConstant() / MyConsts_K2AU,' K'
      write(*,*) "average init vc=",vav, " au of velocity"
      write(*,*) "average init Kc=",0.5*vsqav*rmc / MyConsts_K2AU,' K'
      WRITE(*,"(/,/)")

      DO kk=1,10
         DO iStep=1,ntime
            zcav(kk,iStep)=zcav(kk,iStep)/float(nprint*inum*(irho+1))
            vz2av(kk,iStep)=vz2av(kk,iStep)/float(nprint*inum*(irho+1))
            vz2av(kk,iStep)=0.5*rmc*vz2av(kk,iStep) / MyConsts_K2AU
         END DO
      END DO

      PRINT*, " "
      DO iStep=1,ntime
         write(6,435) dt*real(iStep*nprint)/MyConsts_fs2AU, ( zcav(kk,iStep),  kk=1,10 )
      END DO

      PRINT*, " "
      DO iStep=1,ntime
         write(6,436) dt*real(iStep*nprint)/MyConsts_fs2AU, ( vz2av(kk,iStep), kk=1,10 )
      END DO

      ! write trapping data to mathematica file
      open(unit=15,file='ptrap.nb',status='replace')
      write(15,'(A6)') 'data={'
      DO jRho=1,irho+1
         write(15,'(A1)') '{'
         DO iStep=1,ntime
            IF(iStep.eq.ntime)THEN
               write(15,'(F9.6)') ptrap(jRho,iStep)/inum
            END IF
            IF(iStep.ne.ntime)THEN 
               write(15,'(F9.6,A1)') ptrap(jRho,iStep)/inum,','
            END IF
         END DO
         IF(jRho.ne.irho+1)THEN
            write(15,'(A2)')'},'
         END IF
         IF(jRho.eq.irho+1)THEN
            write(15,'(A2)') '};'
         END IF
      END DO
      write(15,505)'ListPlot3D[data,ColorFunction->Hue,AxesLabel->', &
            '{"Time(ps)","Reaction Parameter(b)","Ptrap(t,b)"},Boxed->False]'
      CLOSE(15)

      ! write cross section data 
      write(6,*)'cross section vs time(ps):'

      ! Cycle over analysis steps
      DO iStep=1,ntime

         ! integrate over rxn parameter
         crscn(iStep)=0.0
         DO jRho=1,irho+1
            ImpactPar = 0.000001 + delrho * real(jRho-1) * MyConsts_Bohr2Ang
            crscn(iStep)=crscn(iStep)+(ptrap(jRho,iStep)/float(inum))* ImpactPar
         END DO
         crscn(iStep) = 2.0* MyConsts_PI * delrho*MyConsts_Bohr2Ang * crscn(iStep)

         ! print trapping p to file
         write(6,605) dt*real(iStep*nprint)/MyConsts_fs2AU,  crscn(iStep)
      END DO

      435 FORMAT(1f6.3,10f7.3)
      436 FORMAT(1f6.3,10f7.2)
      505 FORMAT(A46,A63)
      605 FORMAT(f8.4,3x,f10.5)

   END SUBROUTINE ScatteringSimulation



!***************************************************************************************************
!                       H-GRAPHITE EQUILIBRIUM SIMULATION
!***************************************************************************************************

   SUBROUTINE HGraphiteEquilibrium()
      IMPLICIT NONE

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! allocate arrays for deltaZ, its average and and its fourier transform
      ALLOCATE( DeltaZ(0:ntime), AverageDeltaZ(0:ntime) )
      AverageDeltaZ(:) = cmplx(0.0, 0.0)
      ! frequency spacing of the fourier transform
      dOmega =  2.*MyConsts_PI/( real(nstep) * dt )

      ! Initialize global temperature average
      GlobalTemperature = 0.0

      ! Open output file to print the brownian realizations vs time
      TrajOutputUnit = LookForFreeUnit()
      OPEN( FILE="Trajectories_CH.dat", UNIT=TrajOutputUnit )
      WRITE(TrajOutputUnit, "(A,I6,A,/)") "# ", inum, " realizations of H-graphite equilibrium dynamics (fs | Angstrom)"
      ! Open output file to print the second brownian realizations vs time
      Traj2OutputUnit = LookForFreeUnit()
      OPEN( FILE="Trajectories_Zh.dat", UNIT=Traj2OutputUnit )
      WRITE(Traj2OutputUnit, "(A,I6,A,/)") "# ", inum, " realizations of zH equilibrium dynamics (fs | Angstrom)"
      ! Open output file to print the spectral density of the brownian motion
      SpectrDensUnit = LookForFreeUnit()
      OPEN( FILE="SpectralDensity.dat", UNIT=SpectrDensUnit )
      WRITE(SpectrDensUnit, "(A,I6,A,/)") "# Spectral density of H-graphite equil dynamics - ", inum, " trajectories (fs | au)"
      ! Open output file to print the  autocorrelation function of the brownian motion
      TimeCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="Autocorrelation.dat", UNIT=TimeCorrelationUnit )
      WRITE(TimeCorrelationUnit, "(A,I6,A,/)") "# Autocorrelation of H-graphite equil dyn - ", inum, " trajectories (fs | Ang^2)"

      ! Initialize random number seed
      CALL SetSeed( 1 )

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               EQUILIBRIUM SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "

      IF ( RunType == OSCIBATH_EQUIL )  PRINT "(A)"," Bath represented as indipendent oscillators in normal form. "
      IF ( RunType == CHAINBATH_EQUIL ) PRINT "(A)"," Bath represented as indipendent oscillators in chain form. "
      IF ( RunType == EQUILIBRIUM )  PRINT "(A)"," Bath represented with atomistic force field. "
      
      !run inum number of trajectories at the current rxn parameter
      DO iTraj = 1,inum
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions
         IF ( RunType == EQUILIBRIUM ) THEN 
               CALL ThermalEquilibriumConditions( X, V, temp, rmh, rmc )
         ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL ) THEN
               CALL InitialBathConditions( X, V, temp )
         END IF
               
         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", temp / MyConsts_K2AU


         IF ( PrintType == DEBUG ) THEN
            ! Open file to store equilibration
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,".dat"
            NWriteUnit = LookForFreeUnit()
            OPEN( Unit=NWriteUnit, File=OutFileName )
            WRITE( NWriteUnit, * ) "# EQUILIBRATION: time, potential E, kinetic E, X_H, Y_H, Z_H, Z_C1"
         ENDIF

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         ! Compute starting potential and forces
         A(:) = 0.0
         IF ( RunType == EQUILIBRIUM ) THEN 
               PotEnergy = VHSticking( X, A )
         ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL ) THEN
               PotEnergy = PotentialIndepOscillatorsModel( X, A )
         END IF

         A(1:3) = A(1:3) / rmh
         IF (  RunType == CHAINBATH_EQUIL ) THEN
             A(4) = A(4) / rmc   
! !              A(5:nevo+3) = A(5:nevo+3) / rmc   
         ELSE
             A(4:nevo+3) = A(4:nevo+3) / rmc   
         ENDIF

         ! Do an equilibration run
         DO iStep = 1, NrEquilibSteps

               IF ( iStep == 1 ) THEN           ! If first step, previous acceleration are not available
                  ! Store initial accelerations
                  APre(:) = A(:)
                  ! Propagate for one timestep with Velocity-Verlet and langevin thermostat
                  IF ( RunType == EQUILIBRIUM ) THEN 
                        CALL EOM_VelocityVerlet( Equilibration, X, V, A, VHSticking, PotEnergy )
                  ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL ) THEN
                        CALL EOM_VelocityVerlet( Equilibration, X, V, A, PotentialIndepOscillatorsModel, PotEnergy )
                  END IF
               ELSE
                  ! Propagate for one timestep with Beeman's method and langevin thermostat
                  IF ( RunType == EQUILIBRIUM ) THEN 
                        CALL EOM_Beeman( Equilibration, X, V, A, APre, VHSticking, PotEnergy )
                  ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL ) THEN
                        CALL EOM_Beeman( Equilibration, X, V, A, APre, PotentialIndepOscillatorsModel, PotEnergy )
                  END IF
               END IF

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

               ! every nprint steps, compute trapping and write debug output
               IF ( mod(iStep,nprint) == 0 ) THEN
                  IF ( PrintType == DEBUG ) THEN
                        WRITE(NWriteUnit,850) real(iStep)*EquilTStep/MyConsts_fs2AU, PotEnergy*MyConsts_Hartree2eV, &
                                              KinEnergy*MyConsts_Hartree2eV,  X(1:4)
                  END IF
               END IF
               850 FORMAT( F12.5, 6F13.6 )

               ! If the step is after some time, store temperature averages
               TempAverage = TempAverage + IstTemperature
               TempVariance = TempVariance + IstTemperature**2
               IF ( PrintType == EQUILIBRDBG ) THEN
                  IF (iStep == 1) WRITE(567,* ) " "
                  IF ( mod(iStep,nprint) == 0 ) WRITE(567,* )  real(iStep)*EquilTStep/MyConsts_fs2AU, TempAverage/iStep, &
                                 sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
               ENDIF

         END DO

         IF ( PrintType == DEBUG ) THEN
            ! Close file to store equilibration
            CLOSE( Unit=NWriteUnit )
         ENDIF

         PRINT "(A)", " Equilibration completed! "
        
         ALLOCATE( XTEMP( size(X) ) )
         DO i = 1, size(X)
            PotEnergy = PotentialIndepOscillatorsModel( X, A )
            TotEnergy = -A(i)
!             print*, " analytical ", -A(i)
            PotEnergy = 0.0
            XTEMP = X + (/ (0.0, iStep=1,i-1 ), 0.00001, (0.0, iStep=i+1,size(X) ) /)
            PotEnergy = PotEnergy + PotentialIndepOscillatorsModel( XTEMP, A ) / ( 2 * 0.00001 )
            XTEMP = X - (/ (0.0, iStep=1,i-1 ), 0.00001, (0.0, iStep=i+1,size(X) ) /)
            PotEnergy = PotEnergy - PotentialIndepOscillatorsModel( XTEMP, A ) / ( 2 * 0.00001 )
!             print*, "numerical ", PotEnergy
            PRINT*, " deriv error ", TotEnergy - PotEnergy
         END DO
         STOP

         IF ( PrintType >= FULL ) THEN
            ! Compute average and standard deviation
            TempAverage = TempAverage / NrEquilibSteps 
            TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2
            
            ! output message with average values
            PRINT "(/,A,1F10.4,/,A,1F10.4,/)",  " * Average temperature (K): ", TempAverage, &
                                                " * Standard deviation (K):  ", sqrt(TempVariance)
         ENDIF

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Compute kinetic energy and total energy
         KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
         TotEnergy = PotEnergy + KinEnergy
         IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))
         
         ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
         IF ( PrintType >= FULL ) THEN
            WRITE(*,600)  PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV, &
                           TotEnergy*MyConsts_Hartree2eV, IstTemperature
            600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                           " * Potential Energy (eV)        ",1F10.4,/    &
                           " * Kinetic Energy (eV)          ",1F10.4,/    &
                           " * Total Energy (eV)            ",1F10.4,/    &
                           " * Istantaneous temperature (K) ",1F10.4,/ ) 
         END IF

         ! Initialize the variables for the trajectory averages
         IF ( PrintType >= FULL ) THEN
            TempAverage = 0.0      ! Temperature
            TempVariance = 0.0
            HPosAverage(:) = 0.0      ! H position
            HPosVariance(:) = 0.0
            C1Average = 0.0        ! C1 position
            C1Variance = 0.0
         ENDIF

         IF ( PrintType == DEBUG ) THEN
            ! Open file to store trajectory information
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,".dat"
            NWriteUnit = LookForFreeUnit()
            OPEN( Unit=NWriteUnit, File=OutFileName, Position="APPEND" )
            WRITE( NWriteUnit, "(/,A)" ) "# TRAJECTORY "
         ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kstep=0

         ! Initialize array with C-H distance
         DeltaZ(:) = 0.0

         ! Compute ZH at the beginning of the trajectory
         IF ( RunType == EQUILIBRIUM ) THEN
            DeltaZ(0) = X(3)  - (X(5)+X(6)+X(7))/3 
         ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL) THEN
            DeltaZ(0) = X(3)  
         END IF
        
         ! Write the CH bond distance in output file
         WRITE(TrajOutputUnit,"(F14.8,2F14.8)") 0.0, (X(3)- X(4))*MyConsts_Bohr2Ang

         PRINT "(/,A)", " Propagating the equilibrium H-Graphene in time... "
         
         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,nstep

            ! Propagate for one timestep
            IF ( RunType == EQUILIBRIUM ) THEN 
               CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, PotEnergy )
            ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL ) THEN
               CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, PotentialIndepOscillatorsModel, PotEnergy )
            END IF

!             DO i = 1, size(X)
!                IF ( MOD(iStep,1000) == 0 ) WRITE(1000+i,"(2F20.6)") real(iStep)*dt, X(i)
!             END DO

            ! Compute kin energy and temperature
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

            ! Store ZH for the autocorrelation function
            IF ( RunType == EQUILIBRIUM ) THEN
               ZHinTime = X(3)  - (X(5)+X(6)+X(7))/3
            ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL) THEN
               ZHinTime = X(3)
            END IF

            ! Necessary averages that are always computed
            HPosAverage(3)  = HPosAverage(3)  +   ZHinTime   ! H Z Coordinate
            HPosVariance(3) = HPosVariance(3) + ( ZHinTime )**2
            
            ! These averages  are computed only for a full output
            IF ( PrintType >= FULL ) THEN
               TempAverage = TempAverage + IstTemperature                     ! Temperature
               TempVariance = TempVariance + IstTemperature**2
               DO iCoord = 1, 2                                               ! H X and Y Coordinates
                  HPosAverage(iCoord)  = HPosAverage(iCoord)  + X(iCoord)
                  HPosVariance(iCoord) = HPosVariance(iCoord) + X(iCoord)**2
               END DO
               IF ( RunType == EQUILIBRIUM ) THEN
                  C1Average  = C1Average  + ( X(4) - (X(5)+X(6)+X(7))/3.0 )      ! C1 Z Coordinate
                  C1Variance = C1Variance + ( X(4) - (X(5)+X(6)+X(7))/3.0 )**2
               ELSE IF ( RunType == OSCIBATH_EQUIL .OR. RunType == CHAINBATH_EQUIL) THEN
                  C1Average  = C1Average  +    X(4)       ! C1 Z Coordinate
                  C1Variance = C1Variance +  ( X(4) )**2
               END IF
            ENDIF

            ! Increment global temperature average
            GlobalTemperature = GlobalTemperature + IstTemperature

            ! output to write every nprint steps 
            IF ( mod(iStep,nprint) == 0 ) THEN

               ! increment counter for printing steps
               kstep=kstep+1

               ! Store DeltaZ only in the print steps
               DeltaZ(kstep) = ZHinTime

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(NWriteUnit,800) dt*real(iStep)/MyConsts_fs2AU, PotEnergy*MyConsts_Hartree2eV,  &
                        KinEnergy*MyConsts_Hartree2eV, TotEnergy*MyConsts_Hartree2eV, IstTemperature                    
                  WRITE(NWriteUnit,801)  dt*real(iStep)/MyConsts_fs2AU, X(1:5)*MyConsts_Bohr2Ang, X(8)*MyConsts_Bohr2Ang, &
                                             X(14)*MyConsts_Bohr2Ang, X(17)*MyConsts_Bohr2Ang
                  800 FORMAT("T= ",F12.5," V= ",1F15.8," K= ",1F15.8," E= ",1F15.8," Temp= ", 1F15.8)
                  801 FORMAT("T= ",F12.5," COORDS= ", 8F8.4)
               END IF

               ! Store the trajectory for XYZ printing
               IF ( PrintType >= FULL  .AND. RunType == EQUILIBRIUM ) THEN
                     Trajectory( :, kstep ) = 0.0
                     Trajectory( 1:min(16,3+nevo) , kstep ) = X( 1:min(16,3+nevo) ) 
                     NrOfTrajSteps = kstep
               END IF

               WRITE(TrajOutputUnit,"(F14.8,2F14.8)")  dt*real(iStep)/MyConsts_fs2AU, (X(3)- X(4))*MyConsts_Bohr2Ang


            END IF 

         END DO

         ! Compute ZH averages and 
         HPosAverage(3)  = HPosAverage(3) / nstep
         HPosVariance(3) = (HPosVariance(3) / nstep) - ( HPosAverage(3) )**2

         ! correct the array with deltaZ vs time
         DeltaZ(:) = DeltaZ(:) -  HPosAverage(3)

         ! PRINT deltaZ vs TIME to output file
         DO iStep = 0, ntime
               WRITE(Traj2OutputUnit,"(F14.8,F14.8)")  dt*real(nprint*iStep)/MyConsts_fs2AU, real(DeltaZ(iStep))*MyConsts_Bohr2Ang
         END DO
         WRITE(Traj2OutputUnit,*)  " "
         WRITE(TrajOutputUnit,*)  " "
         
         ! Fourier transform
         CALL DiscreteInverseFourier( DeltaZ )

         ! Store average of DeltaZ fourier components
         AverageDeltaZ(:) = AverageDeltaZ(:) + abs(DeltaZ(:))**2

         PRINT "(A)", " Time propagation completed! "

         IF ( PrintType >= FULL ) THEN

            IF  ( RunType == EQUILIBRIUM ) THEN
               ! write the xyz file of the trajectory, if requested
               WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",iTraj,".xyz"
               CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
                                             GraphiteLatticeConstant()*MyConsts_Bohr2Ang )
            ENDIF

            ! Divide averages for number of time steps and compute variances
            TempAverage  = TempAverage  / nstep
            TempVariance = TempVariance / nstep - TempAverage**2
            HPosAverage(1) = HPosAverage(1) / nstep
            HPosAverage(2) = HPosAverage(2) / nstep
            HPosVariance(1) = HPosVariance(1) / nstep - HPosAverage(1)**2
            HPosVariance(2) = HPosVariance(2) / nstep - HPosAverage(2)**2
            C1Average  = C1Average  / nstep
            C1Variance = C1Variance / nstep - C1Average**2        
                                            
            ! And print the average values of that trajectory 
            WRITE(*,700)  TempAverage, sqrt(TempVariance), HPosAverage(:)* MyConsts_Bohr2Ang, &
                       sqrt(HPosVariance(:))* MyConsts_Bohr2Ang, C1Average* MyConsts_Bohr2Ang, sqrt(C1Variance)* MyConsts_Bohr2Ang
            700 FORMAT (/, " * Average temperature (K)        ",1F10.4,/    &
                           "   Standard deviation (K)         ",1F10.4,/    &
                           " * Average H coordinates (Ang)    ",3F10.4,/    &
                           "   Standard deviation (Ang)       ",3F10.4,/    &
                           " * Average Z height (Ang)         ",1F10.4,/    &
                           "   Standard deviation (Ang)       ",1F10.4,/ ) 
                                                
         END IF

         IF ( PrintType == DEBUG ) THEN
               ! Close file to store equilibration
               CLOSE( Unit=NWriteUnit )
         ENDIF

      END DO
      
      PRINT "(A)"," Done! "

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************
      
      ! Normalize global average of the temperature
      GlobalTemperature = GlobalTemperature / real( inum*nstep )

      WRITE(*,701)  GlobalTemperature
      701 FORMAT (/, " Average temperature for all trajs and time steps ",/   &
                     " * Average temperature (K)        ",1F10.4,/ ) 

      ! Normalize average of deltaZ fourier components
      AverageDeltaZ(:) = AverageDeltaZ(:) / real(inum )

      ! Print spectral density of the stocastic process X(t)
      DO iOmega = 0, ntime
            WRITE(SpectrDensUnit,*)  iOmega*dOmega*MyConsts_fs2AU,  AverageDeltaZ(iOmega)         
      END DO

      DeltaZ(:) = AverageDeltaZ(:)
      CALL DiscreteInverseFourier(DeltaZ)

      ! Print the autocorrelation function
      DO iStep = (ntime/2)+1, ntime
         WRITE(TimeCorrelationUnit,*)  dt*real(nprint*(iStep-ntime-1))/MyConsts_fs2AU,  real(DeltaZ(iStep))*MyConsts_Bohr2Ang**2
      END DO
      DO iStep = 0, ntime/2
         WRITE(TimeCorrelationUnit,*)  dt*real(nprint*iStep)/MyConsts_fs2AU,  real(DeltaZ(iStep))*MyConsts_Bohr2Ang**2
      END DO

      ! Close output files
      CLOSE( TrajOutputUnit )
      CLOSE( Traj2OutputUnit )
      CLOSE( TimeCorrelationUnit )
      CLOSE( SpectrDensUnit )


   END SUBROUTINE HGraphiteEquilibrium



!***************************************************************************************************
!                       BROWNIAN HARMONIC OSCILLATOR SIMULATION
!***************************************************************************************************

   SUBROUTINE BrownianHarmonicOscillator()
      IMPLICIT NONE

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ! allocate arrays for deltaZ, its average and and its fourier transform
      ALLOCATE( DeltaZ(0:ntime), AverageDeltaZ(0:ntime) )
      AverageDeltaZ(:) = cmplx(0.0, 0.0)
      ! frequency spacing of the fourier transform
      dOmega = 2.*MyConsts_PI/( real(nstep) * dt )

      ! Initialize global temperature average
      GlobalTemperature = 0.0

      ! Initialize temperatures distribution
      InitTempDistribution(:) = 0

      ! Open output file to print the distribution of initial temperatures after equilibration
      InitialTDistrib = LookForFreeUnit()
      OPEN( FILE="InitialTemp.dat", UNIT=InitialTDistrib )
      WRITE(InitialTDistrib, "(A,I6,A,/)") "# Distribution of initial temp after equilibration - ", inum, " trajectories (K)"
      ! Open output file to print the brownian realizations vs time
      TrajOutputUnit = LookForFreeUnit()
      OPEN( FILE="Trajectories.dat", UNIT=TrajOutputUnit )
      WRITE(TrajOutputUnit, "(A,I6,A,/)") "# ", inum, " realizations of harmonic brownian motion (fs | Angstrom)"
      ! Open output file to print the spectral density of the brownian motion
      SpectrDensUnit = LookForFreeUnit()
      OPEN( FILE="SpectralDensity.dat", UNIT=SpectrDensUnit )
      WRITE(SpectrDensUnit, "(A,I6,A,/)") "# Spectral density of harmonic brownian motion - ", inum, " trajectories (fs | au)"
      ! Open output file to print the  autocorrelation function of the brownian motion
      TimeCorrelationUnit = LookForFreeUnit()
      OPEN( FILE="Autocorrelation.dat", UNIT=TimeCorrelationUnit )
      WRITE(TimeCorrelationUnit, "(A,I6,A,/)") "# Autocorrelation func of harmonic brownian motion - ", inum, " trajs (fs | Ang^2)"


      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               HARMONIC BROWNIAN TEST"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "
      
      !run inum number of trajectories at the current rxn parameter
      DO iTraj = 1,inum
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions
         CALL HarmonicConditions( X, V )

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", temp / MyConsts_K2AU

         ! Initialize temperature average and variance
         TempAverage = 0.0
         TempVariance = 0.0

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = VHarmonic( X, A )
         A(1:nevo) = A(1:nevo) / rmh

         ! Do an equilibration run
         DO iStep = 1, NrEquilibSteps

               IF ( iStep == 1 ) THEN           ! If first step, previous acceleration are not available
                  ! Store initial accelerations
                  APre(:) = A(:)
                  ! Propagate for one timestep with Velocity-Verlet and langevin thermostat
                  CALL EOM_VelocityVerlet( Equilibration, X, V, A, VHarmonic, PotEnergy )
               ELSE
                  ! Propagate for one timestep with Beeman's method and langevin thermostat
                  CALL EOM_Beeman( Equilibration, X, V, A, APre, VHarmonic, PotEnergy )
               END IF

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy( Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*nevo)

               TempAverage = TempAverage + IstTemperature
               TempVariance = TempVariance + IstTemperature**2

               IF ( PrintType == EQUILIBRDBG ) THEN
                  IF (iStep == 1) WRITE(567,* ) " "
                  IF ( mod(iStep,nprint) == 0 ) WRITE(567,* )  real(iStep)*dt/MyConsts_fs2AU, TempAverage/iStep, &
                                 sqrt((TempVariance/iStep)-(TempAverage/iStep)**2)
               END IF

         END DO

         PRINT "(A)", " Equilibration completed! "

         IF ( PrintType >= FULL ) THEN
            
            ! Compute average and standard deviation
            TempAverage = TempAverage / NrEquilibSteps 
            TempVariance = (TempVariance/NrEquilibSteps) - TempAverage**2

            ! output message with average values
            PRINT "(/,A,1F10.4,/,A,1F10.4,/)",  " * Average temperature (K): ", TempAverage, &
                                                " * Standard deviation (K):  ", sqrt(TempVariance)
         ENDIF

         ! Store the initial temperature
         IF ( (INT( TempAverage ) >= 200) .AND. (INT( TempAverage ) <= 400) ) &
               InitTempDistribution( INT( TempAverage ) ) = InitTempDistribution( INT( TempAverage ) ) + 1 

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Compute kinetic energy and total energy
         KinEnergy = EOM_KineticEnergy( Equilibration, V )
         TotEnergy = PotEnergy + KinEnergy
         IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*nevo)

         
         ! Initialize the variables for the trajectory averages
         TempAverage  = 0.0      ! Temperature
         TempVariance = 0.0
         HPosAverage(1)  = 0.0
         HPosVariance(1) = 0.0

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kstep=0

         ! Initialize array with C-H distance
         DeltaZ(:) = 0.0
         DeltaZ(0) = X(1) 

         PRINT "(/,A)", " Propagating the harmonic brownian motion in time... "

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,nstep

            ! Propagate for one timestep
            CALL EOM_Beeman( Equilibration, X, V, A, APre, VHarmonic, PotEnergy )

            ! Compute kin energy and temperature
            KinEnergy = EOM_KineticEnergy( Equilibration, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*nevo)

            ! Store first dof for the autocorrelation function
            ZHinTime = X(1)

            ! Necessary averages that are always computed
            HPosAverage(1)  = HPosAverage(1)  +   ZHinTime   ! H Z Coordinate
            HPosVariance(1) = HPosVariance(1) + ( ZHinTime )**2
            TempAverage = TempAverage + IstTemperature                     ! Temperature
            TempVariance = TempVariance + IstTemperature**2

            ! Increment global temperature average
            GlobalTemperature = GlobalTemperature + IstTemperature

            ! output to write every nprint steps 
            IF ( mod(iStep,nprint) == 0 ) THEN
               ! increment counter for printing steps
               kstep=kstep+1
               ! Store DeltaZ only in the print steps
               DeltaZ(kstep) = ZHinTime
            END IF 

         END DO

         ! Compute ZH averages and 
         HPosAverage(1)  = HPosAverage(1) / nstep
         HPosVariance(1) = (HPosVariance(1) / nstep) - ( HPosAverage(1) )**2

         ! correct the array with deltaZ vs time
         DeltaZ(:) = DeltaZ(:) -  HPosAverage(1)

         ! PRINT deltaZ vs TIME to output file
         DO iStep = 0, ntime
               WRITE(TrajOutputUnit,"(F14.8,F14.8)")  dt*real(nprint*iStep)/MyConsts_fs2AU,  real(DeltaZ(iStep))*MyConsts_Bohr2Ang
         END DO
         WRITE(TrajOutputUnit,*)  " "

         ! Fourier transform
         CALL DiscreteInverseFourier( DeltaZ )

         ! Store average of DeltaZ fourier components
         AverageDeltaZ(:) = AverageDeltaZ(:) + abs(DeltaZ(:))**2

         ! Divide averages for number of time steps and compute variances
         TempAverage  = TempAverage  / nstep
         TempVariance = TempVariance / nstep - TempAverage**2
                                            
         ! And print the average values of that trajectory 
         PRINT "(A)", " Time propagation completed! "
         IF ( PrintType >= FULL ) THEN
            PRINT 710, TempAverage, sqrt(TempVariance), HPosAverage(1)* MyConsts_Bohr2Ang, sqrt(HPosVariance(1))* MyConsts_Bohr2Ang
            710 FORMAT (/,    " * Average temperature (K)        ",1F10.4,/    &
                              "   Standard deviation (K)         ",1F10.4,/    &
                              " * Average coordinates (Ang)      ",3F10.4,/    &
                              "   Standard deviation (Ang)       ",3F10.4,/ ) 
         END IF                                       

      END DO
      
      PRINT "(A)"," Done! "

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! PRint the distribution of initial istant temperatures
      DO iStep = 200,400
         WRITE(InitialTDistrib,*)  iStep, InitTempDistribution( iStep )
      END DO

      ! Normalize global average of the temperature
      GlobalTemperature = GlobalTemperature / real( inum*nstep )

      WRITE(*,702)  GlobalTemperature
      702 FORMAT (/, " Average temperature for all trajs and time steps ",/   &
                     " * Average temperature (K)        ",1F10.4,/ ) 

      ! Normalize average of deltaZ fourier components
      AverageDeltaZ(:) = AverageDeltaZ(:)  / real(inum)

      ! For comparison, compute the spectral density from the analytic results of the Langevin HO

      ! Compute the variables needed for the analytical autocorrelation function
      NoiseVariance = 2.0*temp*rmh*Gamma
!      OsciFreq = (MyConsts_PI*2./(100.*MyConsts_fs2AU))
      OsciFreq =  ( 0.484969 / MyConsts_fs2AU )
      DistortedFreq = sqrt( OsciFreq**2 - Gamma**2/4.0 )
      MotionVariance = NoiseVariance / ( 2 * rmh**2 * Gamma * OsciFreq**2 )

      ! Store the autocorrelation function
      ALLOCATE( Analytic( 0: ntime ) )
      DO iStep = (ntime/2)+1, ntime
          Time = dt*real(nprint*(iStep-ntime-1)) 
          Analytic( iStep ) = MotionVariance * ( cos(DistortedFreq*Time) + Gamma/(2.*DistortedFreq) & 
                                * sin( - DistortedFreq*Time) ) * exp( Gamma * Time /2.)
      END DO   
      DO iStep = 0, ntime/2
          Time = dt*real(nprint*iStep)
          Analytic( iStep ) = MotionVariance * ( cos(DistortedFreq*Time) + Gamma/(2.*DistortedFreq) & 
                                * sin( DistortedFreq*Time) ) * exp( - Gamma * Time /2.)
      ENDDO
      
      ! Fourier transform to get the spectral density
      CALL DiscreteInverseFourier( Analytic )

      ! Print spectral density of the stocastic process X(t)
      DO iOmega = 0, ntime
            WRITE(SpectrDensUnit,"(1F12.6,3F24.18)")  iOmega*dOmega*MyConsts_fs2AU,  AverageDeltaZ(iOmega)/dOmega,  &
                    real( Analytic(iOmega) )/dOmega, NoiseVariance/(2*rmh**2*MyConsts_PI) / &
                            ( (OsciFreq**2 - (iOmega*dOmega)**2)**2 + (Gamma*iOmega*dOmega)**2  )
      END DO

      ! fourier transform the spectral density to obtain the autocorrelation function
      DeltaZ(:) = AverageDeltaZ(:)
      CALL DiscreteFourier(DeltaZ)

      ! Print the autocorrelation function
      DO iStep = (ntime/2)+1, ntime
         WRITE(TimeCorrelationUnit,*)  dt*real(nprint*(iStep-ntime-1))/MyConsts_fs2AU,  real(DeltaZ(iStep))*MyConsts_Bohr2Ang**2
      END DO
      DO iStep = 0, ntime/2
         WRITE(TimeCorrelationUnit,*)  dt*real(nprint*iStep)/MyConsts_fs2AU,  real(DeltaZ(iStep))*MyConsts_Bohr2Ang**2
      END DO
      WRITE(TimeCorrelationUnit,"(/,A,/)") "# analytical result"

      PRINT*, " "
      PRINT*, " * Expected standard deviation of the brownian oscillator: ", sqrt(MotionVariance)*MyConsts_Bohr2Ang,  " Ang "
      PRINT*, " * Computed standard deviation of the brownian oscillator: ", sqrt(real(DeltaZ(0)))*MyConsts_Bohr2Ang, " Ang "
      PRINT*, " "

      ! print in the autocorrelation file the analytic prediction
      DO iStep = 0, ntime/2
         Time = dt*real(nprint*iStep)
         WRITE(TimeCorrelationUnit,*)  Time/MyConsts_fs2AU,  &
            MotionVariance * ( cos(DistortedFreq*Time) + Gamma/(2.*DistortedFreq) * sin(DistortedFreq*Time) )        &
                                                               * exp( - Gamma * Time /2.)*MyConsts_Bohr2Ang**2
      END DO

      ! Close output files
      CLOSE( TrajOutputUnit )
      CLOSE( TimeCorrelationUnit )
      CLOSE( SpectrDensUnit )
      CLOSE( InitialTDistrib )

   END SUBROUTINE BrownianHarmonicOscillator


!***************************************************************************************************
!                       PRINT POTENTIAL CUTS
!***************************************************************************************************

   SUBROUTINE PotentialCuts()
      IMPLICIT NONE

      REAL, PARAMETER :: dd = 0.001
      
      INTEGER :: iZH, iZC, nPoint
      INTEGER :: HCurveOutputUnit

      REAL, DIMENSION(:), ALLOCATABLE :: PotentialArray
      REAL, DIMENSION(:), ALLOCATABLE :: ZCArray, ZHArray

      REAL :: DeltaZH, DeltaZC
      REAL :: MinimumZH,          MinimumZC,          MinimumE
      REAL :: ZHFreq, ZCFreq, Coupling
      REAL :: Tmp1, Tmp2
      
      LOGICAL, DIMENSION(nevo+3) :: OptimizationMask

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               PES ANALYSIS "
      PRINT "(A,/)" ,    "***************************************************"

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointZH * NpointZC ), ZCArray( NpointZC ), ZHArray( NpointZH ) )
   
      ! Define grid spacings
      DeltaZH = (ZHmax-ZHmin)/(NpointZH-1)
      DeltaZC = (ZCmax-ZCmin)/(NpointZC-1)

      ! Define coordinate grids
      ZCArray = (/ ( ZCmin + DeltaZC*(iZC-1), iZC=1,NpointZC) /)
      ZHArray = (/ ( ZHmin + DeltaZH*(iZH-1), iZH=1,NpointZH) /)

      !*************************************************************
      ! EQUILIBRIUM POSITION OF THE H-GRAPHITE POTENTIAL
      !*************************************************************

      ! minimum of the CH 4D potential

      ! set optimization mask
      IF ( FreezeGraphene == 1 ) THEN
         OptimizationMask = (/ (.FALSE., iCoord=1,nevo+3) /)
         OptimizationMask(1:4) = .TRUE.
      ELSE 
         OptimizationMask = (/ (.TRUE., iCoord=1,nevo+3) /)
      ENDIF
         
      ! Set guess geometry: collinear H and other Cs in ideal geometry
      X(:) = 0.0
      X(4) = C1Puckering
      X(3) = C1Puckering + 2.2

      ! Compute potential
      PotEnergy = MinimizePotential( X,  OptimizationMask )

      ! Store values, with respect to the plane defined by C2,C3 and C4
      Tmp1 = X(3)
      Tmp2 = X(4)
      MinimumZH = X(3) - (X(5)+X(6)+X(7))/3.0
      MinimumZC = X(4) - (X(5)+X(6)+X(7))/3.0
      MinimumE = PotEnergy

      ! ZH frequency at the minimum of the potential
      X(3) = Tmp1 + dd
      PotEnergy = VHSticking( X, A )
      ZHFreq = - A(3) / ( 2.0 * dd )
      X(3) = Tmp1 - dd
      PotEnergy =  VHSticking( X, A )
      ZHFreq = ZHFreq + A(3) / ( 2.0 * dd )
      ZHFreq = sqrt( ZHFreq /  rmh )

      ! ZC frequency at the minimum of the potential
      X(3) = Tmp1
      X(4) = Tmp2 + dd
      PotEnergy = VHSticking( X, A )
      ZCFreq = - A(4) / ( 2.0 * dd )
      X(4) = Tmp2 - dd
      PotEnergy =  VHSticking( X, A )
      ZCFreq = ZCFreq + A(4) / ( 2.0 * dd )
      ZCFreq = sqrt( ZCFreq /  rmc )
      
     
      
      WRITE(*,801) MinimumE*MyConsts_Hartree2eV, MinimumZH*MyConsts_Bohr2Ang, MinimumZC*MyConsts_Bohr2Ang, &
                   (MinimumZH-MinimumZC)*MyConsts_Bohr2Ang, ZHFreq/MyConsts_cmmin1toAU, ZCFreq/MyConsts_cmmin1toAU
      WRITE(*,802) (MinimumE+5*MyConsts_K2AU)*MyConsts_Hartree2eV, (MinimumE+50*MyConsts_K2AU)*MyConsts_Hartree2eV,      &
                   (MinimumE+100*MyConsts_K2AU)*MyConsts_Hartree2eV, (MinimumE+300*MyConsts_K2AU)*MyConsts_Hartree2eV,   &
                   (MinimumE+500*MyConsts_K2AU)*MyConsts_Hartree2eV

      801 FORMAT (/,    " * Optimization with ideal graphite surface ",/  &
                        "   (Z=0 defined by C2,C3 and C4 plane)",         /  &
                        "   Energy at the minimum (eV)     ",1F10.4,/     &
                        "   Z coordinate of H atom (Ang)   ",1F10.4,/     &
                        "   Z coordinate of C1 atom (Ang)  ",1F10.4,/     &
                        "   H-C1 distance (Ang)            ",1F10.4,/     &
                        "   Frequency along ZH (cm-1)      ",1F10.1,/     &
                        "   Frequency along ZC (cm-1)      ",1F10.1           ) 

      802 FORMAT (/,    " * Average vibrational energy at given T (eV) ",/,  &
                        "   at T = 5K   :   ",1F10.4,/, & 
                        "   at T = 50K  :   ",1F10.4,/, &
                        "   at T = 100K :   ",1F10.4,/, &
                        "   at T = 300K :   ",1F10.4,/, & 
                        "   at T = 500K :   ",1F10.4,/ ) 


      !*************************************************************
      ! PRINT H VIBRATIONAL CURVE FOR C FIXED IN PUCKERING GEOMETRY
      !*************************************************************

      ! Open output file to print the H vibrational curve with fixed carbons
      HCurveOutputUnit = LookForFreeUnit()
      OPEN( FILE="GraphiteHBindingCurve.dat", UNIT=HCurveOutputUnit )
      WRITE(HCurveOutputUnit, "(A,I6,A,/)") "# H-Graphite binding at C fixed puckered geometry - ", NpointZH, " points (Ang | eV)"

      ! Set collinear H, C1 puckered and other Cs in ideal geometry
      X(:) = 0.0
      X(4) = MinimumZC

      ! Cycle over the ZH coordinate values
      DO iZH = 1, NpointZH
         ! Define the H Z coordinate
         X(3) = ZHArray(iZH)
         ! Compute potential
         PotEnergy = VHSticking( X, A )
         ! Print energy to dat file
         WRITE(HCurveOutputUnit,*) X(3)*MyConsts_Bohr2Ang, PotEnergy*MyConsts_Hartree2eV
      END DO
      CLOSE( HCurveOutputUnit )

      WRITE(*,"(/,A)") " * CH binding curve written to file GraphiteHBindingCurve.dat"



      !*************************************************************
      ! PRINT 2D C-H MODEL V from FREQUENCIES AND COUPLING
      !*************************************************************

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO iZC = 1, NpointZC
         X(4) = ZCArray(iZC)
         ! Cycle over the ZH coordinate values
         DO iZH = 1, NpointZH
            X(3) = ZHArray(iZH)
            nPoint = nPoint + 1
            ! Compute potential
            PotentialArray(nPoint) = MinimumE + 0.5*rmh*(ZHFrequency*(X(3)-MinimumZH))**2 + &
                  0.5*rmc*(ZCFrequency*(X(4)-MinimumZC))**2 + BilinearCoupl*(X(3)-MinimumZH)*(X(4)-MinimumZC)
         END DO
      END DO

      ! Print the potential to vtk file
      CALL VTK_NewRectilinearSnapshot ( PotentialCH_Model, X=ZHArray*MyConsts_Bohr2Ang, Y=ZCArray*MyConsts_Bohr2Ang, &
                            FileName="ModelPotential" )
      CALL VTK_AddScalarField (PotentialCH_Model, Name="CHPotential", Field=PotentialArray*MyConsts_Hartree2eV )

      WRITE(*,"(/,A,A)") " * Model PES from vibrational frequency and bilinear coupling written in VTR ", & 
                           "format to file ModelPotential.vtr"


      !*************************************************************
      ! PRINT 2D C-H V WITH OTHERS CARBONS IN OPTIMIZED GEOMETRY
      !*************************************************************

      ! set optimization mask
      IF ( FreezeGraphene == 1 ) THEN
         OptimizationMask = (/ (.TRUE., iCoord=1,2), (.FALSE., iCoord=1,nevo+1) /)
      ELSE 
         OptimizationMask = (/ (.TRUE., iCoord=1,2), (.FALSE., iCoord=1,2), (.TRUE., iCoord=1,nevo-4) /)
      ENDIF

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO iZC = 1, NpointZC
         ! Cycle over the ZH coordinate values
         DO iZH = 1, NpointZH
            nPoint = nPoint + 1

            ! Set collinear H and other Cs in ideal geometry
            X(:) = 0.0
            X(4) = ZCArray(iZC)
            X(3) = ZHArray(iZH)

            ! Compute potential at optimized geometry
            PotentialArray(nPoint) = MinimizePotential( X,  OptimizationMask )
         END DO
      END DO
      
      ! Print the potential to vtk file
      CALL VTK_NewRectilinearSnapshot ( PotentialCH, X=ZHArray*MyConsts_Bohr2Ang,  & 
                                Y=ZCArray*MyConsts_Bohr2Ang, FileName="GraphiteHSticking" )
      CALL VTK_AddScalarField (PotentialCH, Name="CHPotential", Field=PotentialArray*MyConsts_Hartree2eV )

      WRITE(*,"(/,A)") " * PES as a func of ZC and ZH with optim puckered graphite written as VTR to file GraphiteHSticking_.vtr"

      ! Deallocate memory
      DEALLOCATE( PotentialArray, ZHArray, ZCArray )

      ! OTHER FEATURE TO IMPLEMENT IN THIS SECTION:
      ! * calculation of the harmonic frequencies at the minimum by diagonalization
      !   of the hessian, in the 2D collinear model and the 4D non collinear model
      ! * by subtraction of the harmonic frequency, computation of the non harmonic coupling 
      !   between the normal modes
      ! * with respect to the minimum, plot of the energies corresponding to average
      !   thermic energy


   END SUBROUTINE PotentialCuts




!*********************************************************************************************************

   SUBROUTINE lattice( kstep, zsum, vzsum, zcav, vz2av )
      IMPLICIT NONE

      REAL, DIMENSION(121), INTENT(INOUT) :: vzsum, zsum
      REAL, DIMENSION(10,ntime), INTENT(INOUT) :: zcav, vz2av 
      INTEGER, INTENT(IN) :: kstep

      zcav(1,kstep)=zcav(1,kstep)+zsum(1)
      zcav(2,kstep)=zcav(2,kstep)+(zsum(2)+zsum(3)+zsum(4))/3.0
      zcav(3,kstep)=zcav(3,kstep)+(zsum(5)+zsum(6)+zsum(7)+zsum(8)+zsum(9)+zsum(10))/6.0
      zcav(4,kstep)=zcav(4,kstep)+(zsum(11)+zsum(12)+zsum(13))/3.0
      zcav(5,kstep)=zcav(5,kstep)+(zsum(14)+zsum(15)+zsum(16)+zsum(17)+zsum(18)+zsum(19))/6.0
      zcav(6,kstep)=zcav(6,kstep)+(zsum(26)+zsum(27)+zsum(28)+zsum(29)+zsum(30)+zsum(31))/6.0
      zcav(7,kstep)=zcav(7,kstep)+(zsum(41)+zsum(42)+zsum(43)+zsum(44)+zsum(45)+zsum(46))/6.0
      zcav(8,kstep)=zcav(8,kstep)+(zsum(62)+zsum(63)+zsum(64)+zsum(65)+zsum(66)+zsum(67))/6.0
      zcav(9,kstep)=zcav(9,kstep)+(zsum(86)+zsum(87)+zsum(88)+zsum(89)+zsum(90)+zsum(91))/6.0
      zcav(10,kstep)=zcav(10,kstep)+(zsum(116)+zsum(117)+zsum(118)+zsum(119)+zsum(120)+zsum(121))/6.0

      vz2av(1,kstep)=vz2av(1,kstep)+vzsum(1)
      vz2av(2,kstep)=vz2av(2,kstep)+(vzsum(2)+vzsum(3)+vzsum(4))/3.0
      vz2av(3,kstep)=vz2av(3,kstep)+(vzsum(5)+vzsum(6)+vzsum(7)+vzsum(8)+vzsum(9)+vzsum(10))/6.0
      vz2av(4,kstep)=vz2av(4,kstep)+(vzsum(11)+vzsum(12)+vzsum(13))/3.0
      vz2av(5,kstep)=vz2av(5,kstep)+(vzsum(14)+vzsum(15)+vzsum(16)+vzsum(17)+vzsum(18)+vzsum(19))/6.0
      vz2av(6,kstep)=vz2av(6,kstep)+(vzsum(26)+vzsum(27)+vzsum(28)+vzsum(29)+vzsum(30)+vzsum(31))/6.0
      vz2av(7,kstep)=vz2av(7,kstep)+(vzsum(41)+vzsum(42)+vzsum(43)+vzsum(44)+vzsum(45)+vzsum(46))/6.0
      vz2av(8,kstep)=vz2av(8,kstep)+(vzsum(62)+vzsum(63)+vzsum(64)+vzsum(65)+vzsum(66)+vzsum(67))/6.0
      vz2av(9,kstep)=vz2av(9,kstep)+(vzsum(86)+vzsum(87)+vzsum(88)+vzsum(89)+vzsum(90)+vzsum(91))/6.0
      vz2av(10,kstep)=vz2av(10,kstep)+(vzsum(116)+vzsum(117)+vzsum(118)+vzsum(119)+vzsum(120)+vzsum(121))/6.0

      zsum(:)  = 0.0
      vzsum(:) = 0.0

   END SUBROUTINE lattice


!*******************************************************************************
! WriteTrajectoryXYZ
!*******************************************************************************
!>  Given an array containing one trajectory over time, print it
!>  to file. 
!> 
!> @param   CoordinatesVsTime  Array with trajectory (size: 7 x N time steps).
!> @param   FileName           Output file name.
!> @param   LatticeConst       Lattice constant of the graphite slab. 
!*******************************************************************************
   SUBROUTINE WriteTrajectoryXYZ( CoordinatesVsTime, FileName, LatticeConst )
      IMPLICIT NONE

         REAL, DIMENSION( :,: ), INTENT(IN) :: CoordinatesVsTime
         CHARACTER(*), INTENT(IN)           :: FileName
         REAL, INTENT(IN)                   :: LatticeConst

         INTEGER :: NrTimeStep, UnitNr, i
         REAL    :: A1_x, A1_y, A2_x, A2_y

         ! check and store array dimensions
         NrTimeStep = size( CoordinatesVsTime, 2 )
         CALL ERROR( size( CoordinatesVsTime, 1 ) /= min(16,3+nevo) , " WriteTrajectoryXYZ: Invalid array dimension "  )

         ! Open input file in a free output unit
         UnitNr =  LookForFreeUnit()
         OPEN( FILE = trim(adjustl(FileName)), UNIT=UnitNr )
   
         A1_x = LatticeConst
         A1_y = 0.0 
         A2_x = LatticeConst/2.0
         A2_y = sqrt(3.)*LatticeConst/2.0

         ! write snapshot to output file
         DO i = 1, NrTimeStep
            WRITE( UnitNr, * ) " 14 "
            WRITE( UnitNr, * ) " Step nr ",i," of the trajectory "
            WRITE( UnitNr, "(A,3F15.6)" ) " H ", CoordinatesVsTime(1,i), CoordinatesVsTime(2,i)     , CoordinatesVsTime(3,i)
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , 0.0                        , CoordinatesVsTime(4,i) !C1
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , LatticeConst/sqrt(3.)      , CoordinatesVsTime(5,i) !C2
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", - LatticeConst/2.0    , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(6,i) !C3
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", LatticeConst/2.0      , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(7,i) !C4
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x+A2_x            , -A1_y+A2_y                 , CoordinatesVsTime(8,i) !C5
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A2_x                 , +A2_y                      , CoordinatesVsTime(9,i) !C6
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x                 , -A1_y                      , CoordinatesVsTime(10,i) !C7
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x                 , +A1_y                      , CoordinatesVsTime(11,i) !C8
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A2_x                 , -A2_y                      , CoordinatesVsTime(12,i) !C9
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x-A2_x            , +A1_y-A2_y                 , CoordinatesVsTime(13,i) !C10
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x                 , LatticeConst/sqrt(3.)-A1_y , CoordinatesVsTime(14,i) !C11
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x                 , LatticeConst/sqrt(3.)+A1_y , CoordinatesVsTime(15,i) !C12
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0+A1_x-2*A2_x , LatticeConst/sqrt(3.)+A1_y-2*A2_y, CoordinatesVsTime(16,i) !C13


         END DO

         ! close input file
         CLOSE( UnitNr )
   
   END SUBROUTINE WriteTrajectoryXYZ

!************************************************************************************************************

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

!************************************************************************************************************

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


   REAL FUNCTION ZCentorOfMass( X )
      IMPLICIT NONE
      REAL, DIMENSION(nevo+3), INTENT(IN) :: X
      INTEGER :: iAtom

      ZCentorOfMass = X(3)*rmh
      DO iAtom = 1, nevo
         ZCentorOfMass = ZCentorOfMass + X(3+iAtom)*rmc
      END DO
      ZCentorOfMass = ZCentorOfMass / ( nevo*rmc + 1.0*rmh )

   END FUNCTION ZCentorOfMass

!************************************************************************************************************

END PROGRAM JK6
