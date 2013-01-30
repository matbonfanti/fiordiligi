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
!***************************************************************************************
PROGRAM JK6
   USE CommonData
   USE InputField
   USE PotentialModule
   USE MyConsts
   USE ErrorTrap
   USE RandomNumberGenerator
   USE ClassicalEqMotion

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

!    REAL, DIMENSION(6) :: HCoordAndVel    ! temporary array to store initial conditions of H atom

   ! Autocorrelation computation
!    REAL, DIMENSION(:,:), ALLOCATABLE  :: ZHinTime
   REAL  :: ZHinTime
   ! Store in time and in frequency the single realization
   COMPLEX, DIMENSION(:), ALLOCATABLE :: DeltaZ
   COMPLEX, DIMENSION(:), ALLOCATABLE :: FreqGrid
   REAL :: dOmega
   REAL, DIMENSION(:), ALLOCATABLE :: AverageDeltaZ

!> \name latt
!> ???
!> @{
   REAL, DIMENSION(121) :: vzsum, zsum
   REAL, DIMENSION(:,:), ALLOCATABLE :: vz2av, zcav
   REAL :: vav, vsqav, zav, zsqav
!> @}

   REAL, DIMENSION(501)      :: crscn

   INTEGER :: iCoord, iTraj, iStep, jRho, iCarbon, iOmega
   INTEGER :: kk, kstep


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
      ! set irho to remove the cycle on rho
      irho = 0
      CALL SetFieldFromInput( InputData, "Gamma", Gamma)
      Gamma = Gamma * MyConsts_Uma2Au / MyConsts_fs2AU
      CALL SetFieldFromInput( InputData, "NrEquilibrSteps", NrEquilibSteps, int(5.0*(1.0/Gamma)/dt) )
      CALL SetFieldFromInput( InputData, "EquilibrInitAverage", EquilibrInitAverage, int(2.0*(1.0/Gamma)/dt) )
   END IF

   ! close input file
   CALL CloseFile( InputData )


   !*************************************************************
   !       ALLOCATION, GENERAL INITIALIZATIONS
   !*************************************************************

   ! number of analysis step
   ntime=int(nstep/nprint)

   ! Allocation of position, velocity, acceleration arrays
   ALLOCATE( X(nevo+3), V(nevo+3), A(nevo+3), APre(nevo+3) )
   
   ! if XYZ files of the trajectories are required, allocate memory to store the traj
   IF ( PrintType >= FULL ) THEN
         ALLOCATE( Trajectory( 7, ntime ) )
   END IF

   ! Set variables for EOM integration in the microcanonical ensamble
   CALL EvolutionSetup( MolecularDynamics, nevo+3, (/ (rmh, iCoord=1,3), (rmc, iCoord=1,nevo) /), dt )
   
   ! Set variables for EOM integration with Langevin thermostat
   CALL EvolutionSetup( Equilibration, nevo+3, (/ (rmh, iCoord=1,3), (rmc, iCoord=1,nevo) /), dt )
   CALL SetupThermostat( Equilibration, Gamma, temp )

   !*************************************************************
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   CALL SetupPotential( )
   
   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

   ! write to standard output the parameters of the calculation

   WRITE(*,898) rmh / MyConsts_Uma2Au, rmc / MyConsts_Uma2Au, temp / MyConsts_K2AU
   WRITE(*,899) dt / MyConsts_fs2AU, nstep, ntime
   IF (RunType == SCATTERING) WRITE(*,900) ezh * MyConsts_Hartree2eV, zhi * MyConsts_Bohr2Ang
   WRITE(*,901) nevo, inum
   IF (RunType == SCATTERING) WRITE(*,902) irho, delrho * MyConsts_Bohr2Ang
   IF (RunType == EQUILIBRIUM) WRITE(*,903) Gamma / MyConsts_Uma2Au * MyConsts_fs2AU,      &
                 dt*real(NrEquilibSteps)/MyConsts_fs2AU, real(EquilibrInitAverage)*dt/MyConsts_fs2AU
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
   903 FORMAT(" * Langevin friction constant (UMA/fs)          ",F10.4,/, &
              " * Equilibration time (fs)                      ",F10.4,/, &
              " * Equilibr averages are computed from (fs)     " F10.4      )

   ! Here the program have a fork: scattering and equilibrium simulations

   !****************************************************************************
   IF (RunType == SCATTERING) THEN
   !****************************************************************************

   !*************************************************************
   !       AVERAGES VARIABLES ALLOCATION AND INITIALIZATION
   !*************************************************************

      ! Initialize average values for the INITIAL CONDITION of the lattice Z coordinates
      zav    = 0.0         ! average position
      zsqav  = 0.0         ! mean squared position
      vav    = 0.0         ! average velocity
      vsqav  = 0.0         ! mean squared velocity

      ! ????????????
      ALLOCATE( zcav(10,ntime), vz2av(10,ntime) )
      zcav(:,:)  = 0.0
      vz2av(:,:) = 0.0

      ! allocate and initialize trapping probability variable
      ALLOCATE( ptrap( irho+1, ntime ) )
      ptrap(:,:) = 0.0

      ! Initialize random generation of initial conditions for scattering simulation
      call random_seed()

      ! scan over the impact parameter... 
      ! this cycle is a trivial DO j=1,1 in the case of an equilibrium simulation
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
                  CALL lattice( kstep )

                  ! Store the trajectory for XYZ printing
                  IF ( PrintType >= FULL ) THEN
                        Trajectory( :, kstep ) = X(1:7)
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



   !****************************************************************************
   ELSE IF (RunType == EQUILIBRIUM) THEN   
   !****************************************************************************

      ! ALLOCATE AND INITIALIZE THOSE DATA THAT WILL CONTAIN THE AVERAGES OVER 
      ! THE SET OF TRAJECTORIES

!       ALLOCATE( ZHinTime(nstep, inum) )
      ! allocate arrays
      ALLOCATE( DeltaZ(0:ntime), AverageDeltaZ(0:ntime) )
      ALLOCATE( FreqGrid(0:ntime) )
      ! frequency spacing
      dOmega =  MyConsts_PI/( real(nstep) * dt )
      ! initialize fourier coeffs
      DO iOmega = 0, ntime
         FreqGrid(iOmega) = MyConsts_I * (2.0*MyConsts_PI*real(iOmega)/real(ntime))
      END DO


      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               EQUILIBRIUM SIMULATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "
      
      !run inum number of trajectories at the current rxn parameter
      DO iTraj = 1,inum
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions
         CALL ThermalEquilibriumConditions( X, V, temp, rmh, rmc )

         PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", temp / MyConsts_K2AU

         IF ( PrintType == DEBUG ) THEN
            ! Open file to store equilibration
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,".dat"
            NWriteUnit = LookForFreeUnit()
            OPEN( Unit=NWriteUnit, File=OutFileName )
            WRITE( NWriteUnit, * ) "# EQUILIBRATION: time, temperature, potential, kinetic energy, total energy"
         ENDIF

         ! Nr of steps to average
         NrOfStepsAver = 0
         IF ( PrintType >= FULL ) THEN
            ! Initialize temperature average and variance
            TempAverage = 0.0
            TempVariance = 0.0
         ENDIF

         ! Compute starting potential and forces
         A(:) = 0.0
         PotEnergy = VHSticking( X, A )
         A(1:3) = A(1:3) / rmh
         A(4:nevo+3) = A(4:nevo+3) / rmc

         ! Do an equilibration run
         DO iStep = 1, NrEquilibSteps

               IF ( iStep == 1 ) THEN           ! If first step, previous acceleration are not available
                  ! Store initial accelerations
                  APre(:) = A(:)
                  ! Propagate for one timestep with Velocity-Verlet and langevin thermostat
                  CALL EOM_VelocityVerlet( Equilibration, X, V, A, VHSticking, PotEnergy )
!                   CALL EOM_VelocityVerlet( Equilibration, X, V, A, VHarmonic, PotEnergy )
               ELSE
                  ! Propagate for one timestep with Beeman's method and langevin thermostat
                  CALL EOM_Beeman( Equilibration, X, V, A, APre, VHSticking, PotEnergy )
!                   CALL EOM_Beeman( Equilibration, X, V, A, APre, VHarmonic, PotEnergy )
               END IF

               ! compute kinetic energy and total energy
               KinEnergy = EOM_KineticEnergy(Equilibration, V )
               TotEnergy = PotEnergy + KinEnergy
               IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

               ! every nprint steps, compute trapping and write debug output
               IF ( mod(iStep,nprint) == 0 ) THEN
                  IF ( PrintType == DEBUG ) THEN
                        WRITE(NWriteUnit,850) real(iStep)*dt/MyConsts_fs2AU, IstTemperature, &
                        PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV,  TotEnergy*MyConsts_Hartree2eV
                  END IF
               END IF
               850 FORMAT( F12.5, 5F15.8 )

               ! If the step is after 2.0/gamma time, store temperature averages
               IF ( iStep >= EquilibrInitAverage ) THEN
                     NrOfStepsAver = NrOfStepsAver + 1
                     IF ( PrintType >= FULL ) THEN
                        TempAverage = TempAverage + IstTemperature
                        TempVariance = TempVariance + IstTemperature**2
                     ENDIF
                     IF ( PrintType == EQUILIBRDBG ) THEN                     
                        IF ( mod(iStep,10000) == 0 ) WRITE(567,* )  real(iStep)*dt/MyConsts_fs2AU, TempAverage/NrOfStepsAver, &
                                       sqrt((TempVariance/NrOfStepsAver)-(TempAverage/NrOfStepsAver)**2)
                     ENDIF
               END IF

         END DO

         IF ( PrintType == DEBUG ) THEN
            ! Close file to store equilibration
            CLOSE( Unit=NWriteUnit )
         ENDIF

         IF ( PrintType >= FULL ) THEN
            ! Compute average and standard deviation
            TempAverage = TempAverage / NrOfStepsAver 
            TempVariance = (TempVariance/NrOfStepsAver) - TempAverage**2

            ! output message with average values
            PRINT "(/,A)",                    " Equilibration completed! "
            PRINT "(A,1F10.4,/,A,1F10.4,/)",  " Average temperature (K): ", TempAverage, &
                                              " Standard deviation (K):  ", sqrt(TempVariance)
         ENDIF


         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Compute starting potential and forces
         PotEnergy = VHSticking( X, A )
         A(1:3) = A(1:3) / rmh
         A(4:nevo+3) = A(4:nevo+3) / rmc
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

         DeltaZ(:) = 0.0
         GraphitePlaneZ = AverageCluster(X(4:))
         DeltaZ(0) = X(3) - GraphitePlaneZ

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,nstep

            ! Propagate for one timestep
            CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, PotEnergy )
!             CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHarmonic, PotEnergy )

            ! Compute kin energy and temperature
            KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
            TotEnergy = PotEnergy + KinEnergy
            IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))
            
            ! Compute the reference position of the graphite plane
            GraphitePlaneZ = AverageCluster(X(4:))

            ! Store ZH-ZHeq for the autocorrelation function
!             ZHinTime(iStep, iTraj) = X(3) - GraphitePlaneZ
!            ZHinTime = X(3) - GraphitePlaneZ
            ZHinTime = X(3) - X(4)

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
               C1Average  = C1Average  + (X(4) - GraphitePlaneZ)              ! C1 Z Coordinate
               C1Variance = C1Variance + (X(4) - GraphitePlaneZ)**2
            ENDIF

            ! output to write every nprint steps 
            IF ( mod(iStep,nprint) == 0 ) THEN

               ! increment counter for printing steps
               kstep=kstep+1

               DeltaZ(kstep) = ZHinTime

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(NWriteUnit,800) dt*real(iStep)/MyConsts_fs2AU, PotEnergy*MyConsts_Hartree2eV,  &
                        KinEnergy*MyConsts_Hartree2eV, TotEnergy*MyConsts_Hartree2eV, IstTemperature                    
                  WRITE(NWriteUnit,801)  dt*real(iStep)/MyConsts_fs2AU, X(1:5)*MyConsts_Bohr2Ang,  &
                           X(8)*MyConsts_Bohr2Ang, X(14)*MyConsts_Bohr2Ang, X(17)*MyConsts_Bohr2Ang
                  800 FORMAT("T= ",F12.5," V= ",1F15.8," K= ",1F15.8," E= ",1F15.8," Temp= ", 1F15.8)
                  801 FORMAT("T= ",F12.5," COORDS= ", 8F8.4)
               END IF

               ! Store the trajectory for XYZ printing
               IF ( PrintType >= FULL ) THEN
                     Trajectory( :, kstep ) = X(1:7)
                     NrOfTrajSteps = kstep
               END IF

!                DO iOmega = 1, ntime
! !                   FourierTrasf(iOmega) = FourierTrasf(iOmega) + ZHinTime(iStep, iTraj) * exp( MyConsts_I * (iStep*dt) * FreqGrid(iOmega))
!                   FourierTrasf(iOmega) = FourierTrasf(iOmega) + ZHinTime * cos( (iStep*dt) * FreqGrid(iOmega))
!                END DO

            END IF 



         END DO

         ! Compute ZH averages and correct the array with deltaZ vs time
         HPosAverage(3)  = HPosAverage(3) / nstep
         HPosVariance(3) = (HPosVariance(3) / nstep) - ( HPosAverage(3) )**2
         DeltaZ(:) = DeltaZ(:) -  HPosAverage(3)

         
         DO iStep = 0, ntime
               WRITE(888,*)  dt*real(nprint*iStep)/MyConsts_fs2AU,  real(DeltaZ(iStep)), aimag(DeltaZ(iStep))         
         END DO
         WRITE(888,*)  " "

         ! Fourier transform
         CALL DiscreteFourier( DeltaZ )

         ! Print fourier
         DO iOmega = 0, ntime
               WRITE(889,*)  iOmega*dOmega*MyConsts_fs2AU,  abs(DeltaZ(iOmega))         
         END DO
         WRITE(889,*)  " "

         ! Store average
         AverageDeltaZ(:) = AverageDeltaZ(:) + abs(DeltaZ(:))**2

         IF ( PrintType >= FULL ) THEN

            ! write the xyz file of the trajectory, if requested
            WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",iTraj,".xyz"
            CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
                                            GraphiteLatticeConstant()*MyConsts_Bohr2Ang )

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
            WRITE(*,700)  TempAverage, sqrt(TempVariance), HPosAverage(:)* MyConsts_Bohr2Ang, sqrt(HPosVariance(:))* MyConsts_Bohr2Ang, &
                          C1Average* MyConsts_Bohr2Ang, sqrt(C1Variance)* MyConsts_Bohr2Ang
            700 FORMAT (/, " Time propagation completed! ",/   &
                           " * Average temperature (K)        ",1F10.4,/    &
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
      
      AverageDeltaZ(:) = AverageDeltaZ(:) / real(inum )

      ! Print fourier
      DO iOmega = 0, ntime
            WRITE(891,*)  iOmega*dOmega*MyConsts_fs2AU,  AverageDeltaZ(iOmega)         
      END DO
      WRITE(891,*)  " "

      DeltaZ(:) = AverageDeltaZ(:)
      CALL DiscreteInverseFourier(DeltaZ)

      DO iStep = (ntime/2)+1, ntime
         WRITE(890,*)  dt*real(nprint*(iStep-ntime-1))/MyConsts_fs2AU,  real(DeltaZ(iStep)), aimag(DeltaZ(iStep))         
      END DO
      DO iStep = 0, ntime/2
         WRITE(890,*)  dt*real(nprint*iStep)/MyConsts_fs2AU,  real(DeltaZ(iStep)), aimag(DeltaZ(iStep))         
      END DO
      WRITE(890,*)  " "

!       DO iTraj = 1, inum, 50
!          CALL PrintCorrelation( ZHinTime, iTraj )
!             
!    !          DO n = 1, size( ZHinTime, 1 )
!    !             WRITE(879,*) dt*real(n)/MyConsts_fs2AU, ZHinTime(n,i)
!    !          END DO
!       END DO
!       CALL PrintCorrelation( ZHinTime, inum )
         
   END IF


      CONTAINS

!****************************************************************************************************************


   SUBROUTINE lattice( kstep )
      IMPLICIT NONE

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

         ! check and store array dimensions
         NrTimeStep = size( CoordinatesVsTime, 2 )
         CALL ERROR( size( CoordinatesVsTime, 1 ) /= 7 , " WriteTrajectoryXYZ: Invalid array dimension "  )

         ! Open input file in a free output unit
         UnitNr =  LookForFreeUnit()
         OPEN( FILE = trim(adjustl(FileName)), UNIT=UnitNr )
   
         ! write snapshot to output file
         DO i = 1, NrTimeStep
            WRITE( UnitNr, * ) " 5 "
            WRITE( UnitNr, * ) " Step nr ",i," of the trajectory "
            WRITE( UnitNr, "(A,3F15.6)" ) " H ", CoordinatesVsTime(1,i), CoordinatesVsTime(2,i)     , CoordinatesVsTime(3,i)
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , 0.0                        , CoordinatesVsTime(4,i)
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , LatticeConst/sqrt(3.)      , CoordinatesVsTime(5,i)
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", - LatticeConst/2.0    , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(6,i)
            WRITE( UnitNr, "(A,3F15.6)" ) " C ", LatticeConst/2.0      , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(7,i)
         END DO

         ! close input file
         CLOSE( UnitNr )
   
   END SUBROUTINE WriteTrajectoryXYZ

!****************************************************************************************************************

   SUBROUTINE PrintCorrelation( ZArray, NrOfTraj )
      IMPLICIT NONE
      REAL, INTENT(IN), DIMENSION(:,:) :: ZArray
      INTEGER, INTENT(IN) :: NrOfTraj
      
      INTEGER :: NrOfTStep, iStep, iTraj, OutUnit
      REAL :: Correlation
      CHARACTER(100) :: FileName, WarnMessage
      LOGICAL, DIMENSION( NrOfTraj ) :: TrajIsOK
      
      NrOfTStep = size( ZArray, 1 )
      CALL ERROR( NrOfTraj > size( ZArray, 2 ), "PrintCorrelation: wrong trajectories number " )
      
      OutUnit = LookForFreeUnit()
      WRITE(FileName,"(A,I4.4,A)") "Correlation_",NrOfTraj,".dat"
      OPEN( Unit=OutUnit, File=FileName )

      DO iTraj = 1, NrOfTraj
         TrajIsOK( iTraj ) = .TRUE.
         TStepCycle: DO iStep = 1, NrOfTStep
            IF ( abs(ZArray(iStep, iTraj)) > 5.0 )  THEN
               TrajIsOK( iTraj ) = .FALSE.
               WRITE(WarnMessage, *) " Excluding trajectory ",iTraj," from autocorrelation average "
               CALL ShowWarning( WarnMessage )
               EXIT TStepCycle
            END IF
         END DO TStepCycle
      END DO
      
      DO iStep = 1, NrOfTStep
          Correlation = 0.0
          DO iTraj = 1, NrOfTraj
             IF ( TrajIsOK( iTraj ) ) Correlation = Correlation + ZArray(1, iTraj)*ZArray(iStep, iTraj)
          END DO
          Correlation = Correlation / COUNT( TrajIsOK )
          WRITE(OutUnit,*) dt*real(iStep-1)/MyConsts_fs2AU, Correlation
      END DO
      
      CLOSE(OutUnit)
      
   END SUBROUTINE PrintCorrelation

   REAL FUNCTION AverageCluster( Z )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: Z

      INTEGER :: k

      IF (nevo == 121) THEN
         AverageCluster = Z(82)+Z(83)+Z(84)
         DO k = 96,101
            AverageCluster =   AverageCluster + Z(k)
         END DO
         DO k = 104,121
            AverageCluster =   AverageCluster + Z(k)
         END DO
         AverageCluster = AverageCluster / real(27)
      ELSE 
         AverageCluster = 0.0
      ENDIF

   END FUNCTION AverageCluster

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

      Vector(:) = Transform(:) / (N-1)

   END SUBROUTINE

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

      Vector(:) = Transform(:)

   END SUBROUTINE

!****************************************************************************************************************


END PROGRAM JK6


      




