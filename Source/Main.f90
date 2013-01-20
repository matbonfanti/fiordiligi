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

   ! Variable to store trapping probability, resolved in time and in rho
   REAL, DIMENSION(:,:), ALLOCATABLE :: ptrap

   REAL :: ImpactPar                     ! impact parameter of H
   REAL, DIMENSION(6) :: HCoordAndVel    ! temporary array to store initial conditions of H atom
   CHARACTER(100) :: OutFileName         ! output file name
   INTEGER  :: NWriteUnit                ! output unit number

   ! Energy
   REAL  :: KinEnergy, PotEnergy, TotEnergy, IstTemperature

   ! Equilibration analysis data
   REAL  :: EquilibrationAverage, EquilibrationVariance, ZHEquilibrium
   INTEGER :: NrOfEquilibrStepsAver
   
   ! Propagation analysis data
   REAL  :: TempAverage, TempVariance
   REAL, DIMENSION(3)  :: HAverage, HVariance
   REAL  :: C1Average, C1Variance
   
   ! Autocorrelation computation
   REAL, DIMENSION(:,:), ALLOCATABLE  :: ZHinTime

   REAL, DIMENSION(501)      :: crscn

   INTEGER :: iCoord
   INTEGER :: i, ipt, j, jrho, kk, kstep, lpt, n, nn, nrho

   REAL :: entot, rho, vv
   REAL :: vav, vsqav, zav, zsqav

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
      Gamma = Gamma / MyConsts_fs2AU
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
   ALLOCATE( X(nevo+3), V(nevo+3), A(nevo+3) )
   
   ! if XYZ files of the trajectories are required, allocate memory to store the traj
   IF ( PrintType >= FULL ) THEN
         ALLOCATE( Trajectory( 7, ntime ) )
   END IF

   ! Set variables for EOM integration in the microcanonical ensamble
   CALL EvolutionSetup( MolecularDynamics, nevo+3, (/ (rmh, n=1,3), (rmc, n=1,nevo) /), dt )
   
   ! Set variables for EOM integration with Langevin thermostat
   CALL EvolutionSetup( Equilibration, nevo+3, (/ (rmh, n=1,3), (rmc, n=1,nevo) /), dt )
   CALL SetupThermostat( Equilibration, Gamma, temp )

   
   !*************************************************************
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   CALL SetupPotential( rmh, rmc )
   
   
   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

   ! write to standard output the parameters of the calculation

   WRITE(*,898) rmh / MyConsts_Uma2Au, rmc / MyConsts_Uma2Au, temp / MyConsts_K2AU
   WRITE(*,899) dt / MyConsts_fs2AU, nstep, ntime
   IF (RunType == SCATTERING) WRITE(*,900) ezh * MyConsts_Hartree2eV, zhi * MyConsts_Bohr2Ang
   WRITE(*,901) nevo, inum
   IF (RunType == SCATTERING) WRITE(*,902) irho, delrho * MyConsts_Bohr2Ang
   IF (RunType == EQUILIBRIUM) WRITE(*,903) (1/gamma)/MyConsts_fs2AU,      &
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
   903 FORMAT(" * Langevin relaxation time (fs)                ",F10.4,/, &
              " * Equilibration time (fs)                      ",F10.4,/, &
              " * Equilibr averages are computed from (fs)     " F10.4      )


   !*************************************************************
   !       AVERAGES VARIABLES ALLOCATION AND INITIALIZATION
   !*************************************************************

   IF (RunType == SCATTERING) THEN

      ! Initialize average values for the INITIAL CONDITION of the lattice Z coordinates
      zav=0.0           ! average position
      zsqav=0.0         ! mean squared position
      vav=0.0           ! average velocity
      vsqav=0.0         ! mean squared velocity

      ! ????????????
      DO kk = 1,10
         DO nn = 1,ntime
            zcav(kk,nn)=0.0
            vz2av(kk,nn)=0.0
         END DO
      END DO

      ! allocate and initialize trapping probability variable
      ALLOCATE( ptrap( irho+1, ntime ) )
      ptrap(:,:) = 0.0

      ! Initialize random generation of initial conditions for scattering simulation
      call random_seed()

   ELSE IF (RunType == EQUILIBRIUM) THEN   

      ! ALLOCATE AND INITIALIZE THOSE DATA THAT WILL CONTAIN THE AVERAGES OVER 
      ! THE SET OF TRAJECTORIES
      ALLOCATE( ZHinTime(nstep, inum) )

   END IF  

   ! scan over the impact parameter... 
   ! this cycle is a trivial DO j=1,1 in the case of an equilibrium simulation
   DO j = 1,irho+1

      IF (RunType == SCATTERING) THEN
         ! Set impact parameter and print message
         ImpactPar = float(j-1)*delrho+0.000001

         PRINT "(2/,A)",    "***************************************************"
         PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar
         PRINT "(A,/)" ,    "***************************************************"

      ELSE IF (RunType == EQUILIBRIUM) THEN
         PRINT "(2/,A)",    "***************************************************"
         PRINT "(A,F10.5)", "               EQUILIBRIUM SIMULATION"
         PRINT "(A,/)" ,    "***************************************************"

      END IF  

      PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "
      
      !run inum number of trajectories at the current rxn parameter
      DO i=1,inum
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", i," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************
       
         IF (RunType == SCATTERING) THEN         ! in case of scattering simulation

            ! Set initial conditions
            CALL HInitScattering( ImpactPar, zhi, ezh )
            CALL LatticeInitialCondition( temp )

         ELSE IF (RunType == EQUILIBRIUM) THEN    ! in case of thermal equilibrium simulation

            ! Set initial conditions
            CALL ThermalEquilibriumConditions( X, V, temp, rmh, rmc )

            PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", temp / MyConsts_K2AU

            IF ( PrintType >= DEBUG ) THEN
               ! Open file to store equilibration
               WRITE(OutFileName,"(A,I4.4,A)") "Traj_",i,".dat"
               NWriteUnit = LookForFreeUnit()
               OPEN( Unit=NWriteUnit, File=OutFileName )
               WRITE( NWriteUnit, * ) "# EQUILIBRATION: time, temperature, potential, kinetic energy, total energy"
            ENDIF

            ! Equilibrium position of Z
            ZHEquilibrium = 0.0
            ! Nr of steps to average
            NrOfEquilibrStepsAver = 0
            IF ( PrintType >= FULL ) THEN
               ! Initialize temperature average and variance
               EquilibrationAverage = 0.0
               EquilibrationVariance = 0.0
            ENDIF

            ! Compute starting potential and forces
            PotEnergy = VHSticking( X, A )
            A(1:3) = A(1:3) / rmh
            A(4:nevo+3) = A(4:nevo+3) / rmc
            
            ! Do an equilibration run
            DO n = 1, NrEquilibSteps

                  ! Propagate for one timestep with langevin thermostat
                  CALL EOM_VelocityVerlet( Equilibration, X, V, A, VHSticking, PotEnergy )

                  ! compute kinetic energy and total energy
                  KinEnergy = EOM_KineticEnergy(Equilibration, V )
                  TotEnergy = PotEnergy + KinEnergy
                  IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))

                  ! every nprint steps, compute trapping and write debug output
                  IF ( mod(n,nprint) == 0 ) THEN
                     IF ( PrintType >= DEBUG ) THEN
                         WRITE(NWriteUnit,850) real(n)*dt/MyConsts_fs2AU, IstTemperature, &
                         PotEnergy*MyConsts_Hartree2eV, KinEnergy*MyConsts_Hartree2eV,  TotEnergy*MyConsts_Hartree2eV
                     END IF
                  END IF
                  850 FORMAT( F12.5, 5F15.8 )

                  ! If the step is after 2.0/gamma time, store temperature averages
                  IF ( n >= EquilibrInitAverage ) THEN
                        ZHEquilibrium = ZHEquilibrium + ( X(3) - (X(5)+X(6)+X(7))/3.0 )
                        IF ( PrintType >= FULL ) THEN
                           EquilibrationAverage = EquilibrationAverage + IstTemperature
                           EquilibrationVariance = EquilibrationVariance + IstTemperature**2
                        END IF
                        NrOfEquilibrStepsAver = NrOfEquilibrStepsAver + 1
                  END IF

            END DO

            IF ( PrintType >= DEBUG ) THEN
               ! Close file to store equilibration
               CLOSE( Unit=NWriteUnit )
            ENDIF

            ZHEquilibrium = ZHEquilibrium / NrOfEquilibrStepsAver
            IF ( PrintType >= FULL ) THEN
               ! Compute average and standard deviation
               EquilibrationAverage = EquilibrationAverage / NrOfEquilibrStepsAver 
               EquilibrationVariance = (EquilibrationVariance/NrOfEquilibrStepsAver) - EquilibrationAverage**2

               ! output message with average values
               PRINT "(/,A)",                    " Equilibration completed! "
               PRINT "(A,1F10.4,/,A,1F10.4,/)",  " Average temperature (K): ", EquilibrationAverage, &
                                                 " Standard deviation (K):  ", sqrt(EquilibrationVariance), &
                                                 " H Equilibrium Z (Ang):   ", ZHEquilibrium*MyConsts_Bohr2Ang
            ENDIF

         END IF

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         IF (RunType == SCATTERING) THEN         ! in case of scattering simulation

            ! increment average values for the initial lattice condition
            DO nn=1,nevo
               zav=zav+z(nn)
               zsqav=zsqav+z(nn)**2
               vav=vav+vzc(nn)
               vsqav=vsqav+vzc(nn)**2
            END DO

            DO nn=1,121
               zsum(nn)=0.0
               vzsum(nn)=0.0
            END DO

            ! compute potential, and total energy
            CALL ComputePotential( xh, yh, zh, z, vv, axh, ayh, azh, azc )
            entot = vv + KineticEnergy() 

         ELSE IF (RunType == EQUILIBRIUM) THEN    ! in case of thermal equilibrium simulation

            ! Compute starting potential and forces
            PotEnergy = VHSticking( X, A )
            A(1:3) = A(1:3) / rmh
            A(4:nevo+3) = A(4:nevo+3) / rmc
!             CALL ComputePotential( xh, yh, zh, z, vv, axh, ayh, azh, azc )
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

            ! INITIALIZE THE RELEVANT VARIABLES FOR THE SINGLE TRAJECTORY AVERAGE
            IF ( PrintType >= FULL ) THEN
               TempAverage = 0.0      ! Temperature
               TempVariance = 0.0
               HAverage(:) = 0.0      ! H position
               HVariance(:) = 0.0
               C1Average = 0.0        ! C1 position
               C1Variance = 0.0
            ENDIF

            IF ( PrintType >= DEBUG ) THEN
               ! Open file to store trajectory information
               WRITE(OutFileName,"(A,I4.4,A)") "Traj_",i,".dat"
               NWriteUnit = LookForFreeUnit()
               OPEN( Unit=NWriteUnit, File=OutFileName, Position="APPEND" )
               WRITE( NWriteUnit, "(/,A)" ) "# TRAJECTORY "
            ENDIF

         END IF


         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kstep=0

         ! cycle over nstep velocity verlet iterations
         DO n = 1,nstep

            ! propagate the traj
            IF (RunType == SCATTERING) THEN

                  ! if scattering simulation, propagate only if in the interaction region
                  IF ( zh <= zhi )   CALL VelocityVerlet(  )

                  ! compute average values of the carbon positions
                  DO nn = 1,121
                     zsum(nn)=zsum(nn)+z(nn)
                     vzsum(nn)=vzsum(nn)+vzc(nn)**2
                  END DO

            ELSE IF (RunType == EQUILIBRIUM) THEN

                  ! Propagate for one timestep
                  CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, vv )
!                   CALL VelocityVerlet(  )

                  ! Compute kin energy and temperature
                  PotEnergy = vv
                  KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
                  TotEnergy = PotEnergy + KinEnergy
                  IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))
                  
                  ! Store ZH-ZHeq for the autocorrelation function
                  ZHinTime(n, i) = ( X(3) - (X(5)+X(6)+X(7))/3.0 ) - ZHEquilibrium
                  
                  ! Increment variables for averages
                  IF ( PrintType >= FULL ) THEN
                     TempAverage = TempAverage + IstTemperature
                     TempVariance = TempVariance + IstTemperature**2
                     DO iCoord = 1, 2 
                        HAverage(iCoord) = HAverage(iCoord) + X(iCoord)
                        HVariance(iCoord) = HVariance(iCoord) + X(iCoord)**2
                     END DO
                     HAverage(3) = HAverage(3) + X(3)-X(4)
                     HVariance(3) = HVariance(3) + ( X(3)-X(4) )**2
                     C1Average = C1Average + (X(4) - (X(5)+X(6)+X(7))/3.0)
                     C1Variance = C1Variance + (X(4) - (X(5)+X(6)+X(7))/3.0)**2
                  ENDIF
                  
            END IF

            ! every nprint steps, compute trapping and write output
            IF ( mod(n,nprint) == 0 ) THEN

               ! increment counter for printing steps
               kstep=kstep+1

               IF (RunType == SCATTERING) THEN

                  ! check if H is in the trapping region
                  IF ( zh <= AdsorpLimit ) ptrap(j,kstep) = ptrap(j,kstep)+1.0

                  ! average lattice velocity
                  CALL lattice( kstep )

               ELSE IF (RunType == EQUILIBRIUM) THEN

                  ! If massive level of output, print traj information to std out
                  IF ( PrintType >= DEBUG ) THEN
                     WRITE(NWriteUnit,800) dt*real(n)/MyConsts_fs2AU, PotEnergy*MyConsts_Hartree2eV,  &
                         KinEnergy*MyConsts_Hartree2eV, TotEnergy*MyConsts_Hartree2eV, IstTemperature                    
                     WRITE(NWriteUnit,801)  dt*real(n)/MyConsts_fs2AU, X(1:5)*MyConsts_Bohr2Ang,  &
                             X(8)*MyConsts_Bohr2Ang, X(14)*MyConsts_Bohr2Ang, X(17)*MyConsts_Bohr2Ang
                     800 FORMAT("T= ",F12.5," V= ",1F15.8," K= ",1F15.8," E= ",1F15.8," Temp= ", 1F15.8)
                     801 FORMAT("T= ",F12.5," COORDS= ", 8F8.4)
                  END IF

               END IF

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
            WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",j,"_",i,".xyz"
            CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
                                            GraphiteLatticeConstant()*MyConsts_Bohr2Ang )

            ! Divide averages for number of time steps and compute variances
            TempAverage = TempAverage / nstep
            TempVariance = TempVariance / nstep - TempAverage**2
            HAverage(:) = HAverage(:) / nstep
            HVariance(1) = HVariance(1) / nstep - HAverage(1)**2
            HVariance(2) = HVariance(2) / nstep - HAverage(2)**2
            HVariance(3) = HVariance(3) / nstep - HAverage(3)**2
            C1Average = C1Average / nstep
            C1Variance = C1Variance / nstep - C1Average**2        
                                            
            ! And print the average values of that trajectory 
            WRITE(*,700)  TempAverage, sqrt(TempVariance), HAverage(:)* MyConsts_Bohr2Ang, sqrt(HVariance(:))* MyConsts_Bohr2Ang, &
                          C1Average* MyConsts_Bohr2Ang, sqrt(C1Variance)* MyConsts_Bohr2Ang
            700 FORMAT (/, " Time propagation completed! ",/   &
                           " * Average temperature (K)        ",1F10.4,/    &
                           "   Standard deviation (K)         ",1F10.4,/    &
                           " * Average H coordinates (Ang)    ",3F10.4,/    &
                           "   Standard deviation (Ang)       ",3F10.4,/    &
                           " * Average Z height (Ang)         ",1F10.4,/    &
                           "   Standard deviation (Ang)       ",1F10.4,/ ) 
                                                
         END IF

         IF ( (RunType == EQUILIBRIUM) .AND. ( PrintType >= DEBUG ) ) THEN
               ! Close file to store equilibration
               CLOSE( Unit=NWriteUnit )
         ENDIF

      END DO
      
      PRINT "(A)"," Done! "

   END DO

   !*************************************************************
   !         OUTPUT OF THE RELEVANT AVERAGES 
   !*************************************************************
   
   IF (RunType == SCATTERING) THEN
            
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
            DO nn=1,ntime
               zcav(kk,nn)=zcav(kk,nn)/float(nprint*inum*(irho+1))
               vz2av(kk,nn)=vz2av(kk,nn)/float(nprint*inum*(irho+1))
               vz2av(kk,nn)=0.5*rmc*vz2av(kk,nn) / MyConsts_K2AU
            END DO
         END DO

         PRINT*, " "
         DO ipt=1,ntime
            write(6,435) dt*real(n)/MyConsts_fs2AU,(zcav(kk,ipt),kk=1,10)
         END DO

         PRINT*, " "
         DO ipt=1,ntime
            write(6,436) dt*real(n)/MyConsts_fs2AU,(vz2av(kk,ipt),kk=1,10)
         END DO

         ! write trapping data to mathematica file
         open(unit=15,file='ptrap.nb',status='replace')
         write(15,'(A6)') 'data={'
         DO nrho=1,irho+1
            write(15,'(A1)') '{'
            DO lpt=1,ntime
               IF(lpt.eq.ntime)THEN
                  write(15,'(F9.6)') ptrap(nrho,lpt)/inum
               END IF
               IF(lpt.ne.ntime)THEN 
                  write(15,'(F9.6,A1)') ptrap(nrho,lpt)/inum,','
               END IF
            END DO
            IF(nrho.ne.irho+1)THEN
               write(15,'(A2)')'},'
            END IF
            IF(nrho.eq.irho+1)THEN
               write(15,'(A2)') '};'
            END IF
         END DO
         write(15,505)'ListPlot3D[data,ColorFunction->Hue,AxesLabel->', &
               '{"Time(ps)","Reaction Parameter(b)","Ptrap(t,b)"},Boxed->False]'
         CLOSE(15)

         ! write cross section data 
         write(6,*)'cross section vs time(ps):'
         DO ipt=1,ntime
            crscn(ipt)=0.0
            ! integrate over rxn parameter
            DO jrho=1,irho+1
               rho=0.000001 + delrho * real(jrho-1) * MyConsts_Bohr2Ang
               crscn(ipt)=crscn(ipt)+(ptrap(jrho,ipt)/float(inum))*rho
            END DO
            crscn(ipt) = 2.0* MyConsts_PI * delrho*MyConsts_Bohr2Ang * crscn(ipt)
            write(6,605) dt*real(ipt*nprint)/MyConsts_fs2AU,  crscn(ipt)
         END DO

         435 FORMAT(1f6.3,10f7.3)
         436 FORMAT(1f6.3,10f7.2)
         505 FORMAT(A46,A63)
         605 FORMAT(f8.4,3x,f10.5)

   ELSE IF (RunType == EQUILIBRIUM) THEN

      DO i = 1, inum, 50
         CALL PrintCorrelation( ZHinTime, i )
      END DO
      CALL PrintCorrelation( ZHinTime, inum )
      
   END IF


      CONTAINS

!****************************************************************************************************************

   SUBROUTINE HInitScattering( ImpPar, Z, IncEnergy )
      IMPLICIT NONE
      REAL, INTENT(IN) :: ImpPar, Z, IncEnergy

      ! Parallel position
      xh = ImpPar
      yh = 0.000001
      ! Parallel velocity
      vxh = 0.0
      vyh = 0.0
      ! Perpendicular position and velocity
      zh = Z
      vzh = - sqrt( 2.0* IncEnergy / rmh )

!       HInitScatterin(1) = 0.000001+delrho*float(j-1)
!       HInitScatterin(3) = zhi
!       HInitScatterin(6) = -sqrt(2.0*ezh/(ev_tu*rmh))

   END SUBROUTINE HInitScattering

!****************************************************************************************************************

   SUBROUTINE LatticeInitialCondition( Temperature )
      IMPLICIT NONE

      REAL, INTENT(IN) :: Temperature
      ! NO NEED for T <-> E conversion factor: T is already in Hartree!

      INTEGER :: nn
      REAL :: stdevz, stdevp
      REAL :: gaus1, gaus2, gvar1, gvar2, gr1, gr2, gs1, gs2

      ! standard deviation of the position distribution
      stdevz = sqrt( Temperature / CarbonForceConstant() )
      ! standard deviation of the momentum distribution
      stdevp = sqrt( rmc * Temperature )

      DO nn = 1,121
         z(nn)   = 0.0
         vzc(nn) = 0.0
      END DO

      DO nn=1,nevo
         DO 
            call random_number(gaus1)
            call random_number(gaus2)
            gvar1=2.0*gaus1-1
            gvar2=2.0*gaus2-1
            gr1=gvar1**2+gvar2**2
            IF (gr1 < 1.0) EXIT 
         END DO

         gr2=sqrt(-2.0*alog(gr1)/gr1)
         gs1=gvar1*gr2
         gs2=gvar2*gr2
         z(nn)=stdevz*gs1
         vzc(nn)=stdevp*gs2/rmc
      END DO

   END SUBROUTINE LatticeInitialCondition

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

      DO nn = 1, 121
         zsum(nn)  = 0.0
         vzsum(nn) = 0.0
      END DO

   END SUBROUTINE lattice

! *****************************************************************************************
   ! IF the optional argument is given, a Langevin thermostat is added to the system
   SUBROUTINE VelocityVerlet( )
      IMPLICIT NONE
         REAL, DIMENSION(124) :: Position, Acceleration
         INTEGER :: jDoF

         ! (1) FULL TIME STEP FOR THE POSITIONS
         xh=xh+vxh*dt+0.5*axh*(dt**2)
         yh=yh+vyh*dt+0.5*ayh*(dt**2)
         zh=zh+vzh*dt+0.5*azh*(dt**2)
         DO nn=1,nevo
            z(nn)=z(nn)+vzc(nn)*dt+0.5*azc(nn)*(dt**2)
         END DO

         ! (2) HALF TIME STEP FOR THE VELOCITIES
         vxh=vxh+0.5*axh*dt
         vyh=vyh+0.5*ayh*dt
         vzh=vzh+0.5*azh*dt
         DO nn = 1,nevo
            vzc(nn)=vzc(nn)+0.5*azc(nn)*dt
         END DO

         ! (3) NEW FORCES AND ACCELERATIONS 
         Position(1) = xh
         Position(2) = yh
         Position(3) = zh
         Acceleration(1) = axh
         Acceleration(2) = ayh
         Acceleration(3) = azh
         DO jDoF = 1, 121
            Position(3+jDoF) = z(jDoF)
            Acceleration(3+jDoF) = azc(jDoF)
         END DO            
         vv = VHSticking( Position, Acceleration )  
         xh = Position(1)
         yh = Position(2)
         zh = Position(3)
         axh = Acceleration(1) / rmh
         ayh = Acceleration(2) / rmh
         azh = Acceleration(3) / rmh
         DO jDoF = 1, 121
            z(jDoF) = Position(3+jDoF)
            azc(jDoF) = Acceleration(3+jDoF) /  rmc
         END DO                 
!        CALL ComputePotential( xh, yh, zh, z, vv, axh, ayh, azh, azc )

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
         vxh=vxh+0.5*axh*dt
         vyh=vyh+0.5*ayh*dt
         vzh=vzh+0.5*azh*dt
         DO nn = 1,nevo
            vzc(nn)=vzc(nn)+0.5*azc(nn)*dt
         END DO

   END SUBROUTINE VelocityVerlet

!****************************************************************************************************************

   REAL FUNCTION KineticEnergy()
      IMPLICIT NONE

         KineticEnergy = 0.5 * rmh * (vzh**2+vyh**2+vxh**2)
         DO nn = 1,nevo
            KineticEnergy = KineticEnergy + 0.5 * rmc * vzc(nn)**2
         END DO

   END FUNCTION KineticEnergy

!****************************************************************************************************************

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
      CHARACTER(100) :: FileName
      
      NrOfTStep = size( ZArray, 1 )
      CALL ERROR( NrOfTraj > size( ZArray, 2 ), "PrintCorrelation: wrong trajectories number " )
      
      OutUnit = LookForFreeUnit()
      WRITE(FileName,"(A,I4.4,A)") "Correlation_",NrOfTraj,".dat"
      OPEN( Unit=OutUnit, File=FileName )
      
      DO iStep = 1, NrOfTStep
          Correlation = 0.0
          DO iTraj = 1, NrOfTraj
             Correlation = Correlation + ZArray(1, iTraj)*ZArray(iStep, iTraj)
          END DO
          Correlation = Correlation / NrOfTraj
          WRITE(OutUnit,*) dt*real(iStep-1)/MyConsts_fs2AU, Correlation
      END DO
      
      CLOSE(OutUnit)
      
   END SUBROUTINE PrintCorrelation

!****************************************************************************************************************


END PROGRAM JK6


      




