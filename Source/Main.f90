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

   IMPLICIT NONE

   ! Variable to define adsorption condition
   REAL, PARAMETER :: AdsorpLimit = 2.8

   ! Various variable to handle the input
   TYPE(InputFile) :: InputData
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.
   CHARACTER(120) :: InputFileName

   ! Variable to define which kind of calculation is required 
   INTEGER :: RunType
   INTEGER, PARAMETER :: SCATTERING = 1, &
                         EQUILIBRIUM = 2

   ! Variable to set print of the trajectory in XYZ format
   INTEGER :: PrintType
   INTEGER, PARAMETER :: DEBUG = 3,   &   ! fully detailed information about the trajs
                         FULL = 2,    &   ! files to plot the make animations 
                         MINIMAL = 1      ! minimal level of output

   ! GENERAL GUIDE LINES of STANDARD OUTPUT
   ! MINIMAL: print only the final average values per impact parameter and integrated along rho
   ! FULL: print also the relevant averages per single trajectory, print the trajectory to other dat files
   ! DEBUG: print average values along all the steps of the trajectory

   ! Real array to store the single trajectory
   REAL, DIMENSION(:,:), ALLOCATABLE :: Trajectory 
   INTEGER :: NrOfTrajSteps

   ! Variable to store trapping probability, resolved in time and in rho
   REAL, DIMENSION(:,:), ALLOCATABLE :: ptrap

   REAL :: ImpactPar                     ! impact parameter of H
   REAL, DIMENSION(6) :: HCoordAndVel    ! temporary array to store initial conditions of H atom
   CHARACTER(100) :: OutFileName         ! output file name
   INTEGER  :: NWriteUnit                ! output unit number

   ! Langevin equilibration parameters
   REAL  ::  Gamma,  LangIntegr,  GammaNoiseH, GammaNoiseC  
   INTEGER  ::   NrEquilibSteps, EquilibrInitAverage

   ! Energy
   REAL  :: KinEnergy, PotEnergy, TotEnergy, IstTemperature

   ! Equilibration analysis data
   REAL  :: EquilibrationAverage, EquilibrationVariance
   INTEGER :: NrOfEquilibrStepsAver

   ! Initial kinetic energy and position for scattering hydrogen
   REAL  :: ezh = 0.0, zhi = 0.0

   REAL, DIMENSION(501)      :: crscn

   INTEGER :: i, inum, ipt, irho, j, jrho, kk, kstep, lpt, n, nn, nprint, nrho
   INTEGER :: nstep, ntime

   REAL :: delrho, dt, entot, rho, temp, vv
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

   ! read parameters common to all calculations and transform them in atomic units
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

   ! output parameters
   CALL SetFieldFromInput( InputData, "PrintType", PrintType )

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

   ! if XYZ files of the trajectories are required, allocate memory to store the traj
   IF ( PrintType >= FULL ) THEN
         ALLOCATE( Trajectory( 7, ntime ) )
   END IF

   ! Set variables for Langevin dynamics integration
   LangIntegr  =  1.0 - 0.5 * Gamma * dt  
   GammaNoiseH = sqrt( 2.0 * temp * rmh * Gamma / dt )
   GammaNoiseC = sqrt( 2.0 * temp * rmc * Gamma / dt )

   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

   ! write to standard output the parameters of the calculation

   WRITE(*,898) rmh / MyConsts_Uma2Au, rmc / MyConsts_Uma2Au, temp / MyConsts_K2AU
   WRITE(*,899) dt / MyConsts_fs2AU, nstep, ntime
   IF (RunType == SCATTERING) WRITE(*,900) ezh * MyConsts_Hartree2eV, zhi * MyConsts_Bohr2Ang
   WRITE(*,901) nevo, inum
   WRITE(*,902) irho, delrho * MyConsts_Bohr2Ang
   WRITE(*,903) (1/gamma)/MyConsts_fs2AU, real(NrEquilibSteps)*dt/MyConsts_fs2AU, real(EquilibrInitAverage)*dt/MyConsts_fs2AU
   WRITE(*,"(/)")

   898 FORMAT(" * Mass of the H atom (UMA):                    ",F10.4,/, &
              " * Mass of the C atoms (UMA):                   ",F10.4,/, &
              " * Temperature of the system (Kelvin):          ",F10.4      )
   899 FORMAT(" * Time step of the evolution (fs):             ",F10.4,/, &
              " * Total nr of time step per trajectory:        ",I5   ,/, &
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
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   CALL SetupPotential( rmh, rmc )


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

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************
       
         IF (RunType == SCATTERING) THEN         ! in case of scattering simulation

            ! Set initial conditions
            CALL HInitScattering( ImpactPar, zhi, ezh )
            CALL LatticeInitialCondition( temp )

         ELSE IF (RunType == EQUILIBRIUM) THEN    ! in case of thermal equilibrium simulation

            ! Set initial conditions
            CALL ThermalEquilibriumConditions( temp )

            PRINT "(/,A,F6.1)"," Equilibrating the initial conditions at T = ", temp

            IF ( PrintType >= DEBUG ) THEN
               ! Open file to store equilibration
               WRITE(OutFileName,"(A,I4.4,A)") "Traj_",i,".dat"
               NWriteUnit = LookForFreeUnit()
               OPEN( Unit=NWriteUnit, File=OutFileName )
               WRITE( NWriteUnit, * ) "# EQUILIBRATION: time, temperature, potential, kinetic energy, total energy"
            ENDIF

            IF ( PrintType >= FULL ) THEN
               ! Initialize temperature average and variance
               EquilibrationAverage = 0.0
               EquilibrationVariance = 0.0
               NrOfEquilibrStepsAver = 0
            ENDIF

            ! Do an equilibration run
            DO n = 1, NrEquilibSteps

                  ! Propagate for one timestep with langevin thermostat
                  vv = 0.0
                  CALL VelocityVerlet( .TRUE. )
                  ! compute kinetic energy and total energy
                  PotEnergy = vv
                  KinEnergy = KineticEnergy()
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
                  IF (( n >= EquilibrInitAverage ) .AND. ( PrintType >= FULL )) THEN
                        EquilibrationAverage = EquilibrationAverage + IstTemperature
                        EquilibrationVariance = EquilibrationVariance + IstTemperature**2
                        NrOfEquilibrStepsAver = NrOfEquilibrStepsAver + 1
                  END IF

            END DO

            IF ( PrintType >= DEBUG ) THEN
               ! Close file to store equilibration
               CLOSE( Unit=NWriteUnit )
            ENDIF

            IF ( PrintType >= FULL ) THEN
               ! Compute average and standard deviation
               EquilibrationAverage = EquilibrationAverage / NrOfEquilibrStepsAver 
               EquilibrationVariance = (EquilibrationVariance/NrOfEquilibrStepsAver) - EquilibrationAverage**2

               ! output message with average values
               PRINT "(/,A)",                    " Equilibration completed! "
               PRINT "(A,1F10.4,/,A,1F10.4,/)",  " Average temperature (K): ", EquilibrationAverage, &
                                                " Standard deviation (K): ", sqrt(EquilibrationVariance)
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

            ! compute potential, and total energy
            CALL ComputePotential( xh, yh, zh, z, vv, axh, ayh, azh, azc )
            entot = vv + KineticEnergy()

            ! INITIALIZE THE RELEVANT VARIABLES FOR THE SINGLE TRAJECTORY AVERAGE

            ! PRINT INITIAL CONDITIONS
            IF ( PrintType >= FULL ) THEN
               ! RENDERE PIU' ESPLICATIVO QUESTO OUTPUT!!!!!!!!!!!!!1
               !   WRITE(*,"(/,A10,A25,A25)") "Trajectory", "Initial V", "Initial Total E"
                  WRITE(*,600) i, vv*MyConsts_Hartree2eV, entot*MyConsts_Hartree2eV
                  600 FORMAT( I10, F25.4, F25.4 )
            END IF

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
                  IF ( zh <= zhi )   CALL VelocityVerlet( .FALSE. )

                  ! compute average values of the carbon positions
                  DO nn = 1,121
                     zsum(nn)=zsum(nn)+z(nn)
                     vzsum(nn)=vzsum(nn)+vzc(nn)**2
                  END DO

            ELSE IF (RunType == EQUILIBRIUM) THEN

                  ! Propagate for one timestep
                  CALL VelocityVerlet( .FALSE. )

            END IF

            ! every nprint steps, compute trapping and write output
            IF ( mod(n,nprint) == 0 ) THEN

               ! increment counter for printing steps
               kstep=kstep+1

               IF (RunType == SCATTERING) THEN

                  ! check if H is in the trapping region
                  IF ( zh <= AdsorpLimit ) ptrap(j,kstep) = ptrap(j,kstep)+1.0

                  ! ???????????????????????????????
                  CALL lattice( kstep )

               ELSE IF (RunType == EQUILIBRIUM) THEN

                  ! If massive level of output, print traj information to std out
                  IF ( PrintType >= DEBUG ) THEN
                     entot = vv + KineticEnergy() 
                     WRITE(NWriteUnit,800) dt*real(n)/MyConsts_fs2AU, vv*MyConsts_Hartree2eV, (entot-vv)*MyConsts_Hartree2eV, &
                                  entot*MyConsts_Hartree2eV
                     WRITE(NWriteUnit,801)  dt*real(n)/MyConsts_fs2AU, xh*MyConsts_Bohr2Ang, yh*MyConsts_Bohr2Ang,  &
                                   zh*MyConsts_Bohr2Ang, z(1)*MyConsts_Bohr2Ang, z(2)*MyConsts_Bohr2Ang, z(5)*MyConsts_Bohr2Ang,   &
                                   z(11)*MyConsts_Bohr2Ang, z(14)*MyConsts_Bohr2Ang
                     800 FORMAT( "T= ",F12.5," ENERGY= ", 3F15.8 )
                     801 FORMAT( "T= ",F12.5," COORDS= ", 8F8.4  )
                  END IF

               END IF

               ! Store the trajectory for XYZ printing
               IF ( PrintType >= FULL ) THEN
                     Trajectory( :, kstep ) = (/ xh,yh,zh,z(1),z(2),z(3),z(4) /)
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

            ! And print the average values of that trajectory (?)
            !        -----

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
                  write(15,'(F8.6)') ptrap(nrho,lpt)/inum
               END IF
               IF(lpt.ne.ntime)THEN 
                  write(15,'(F8.6,A1)') ptrap(nrho,lpt)/inum,','
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

   SUBROUTINE ThermalEquilibriumConditions( Temperature )
      IMPLICIT NONE

      REAL, INTENT(IN) :: Temperature
      INTEGER          :: nCarbon
      REAL             :: SigmaCarbonVelocity, SigmaHydroVelocity

      ! All the atoms are initially at the equilibrium position for stable chemisorption 
      ! Value for the puckering are taken from J. Phys. Chem. B, 2006, 110, 18811-18817
      ! Equilibrium position of zH obtained instead from plot of the PES
      ! Velocities are sampled according to a Maxwell-Boltzmann distribution at temperature T

      ! Equilibrium position of H atom
      xh = 0.00001
      yh = 0.00001
      zh = 1.483 / MyConsts_Bohr2Ang

      ! Equilibrium position of C1 atom
      z(1) = 0.37 / MyConsts_Bohr2Ang

      ! Equilibrium position of the other carbon atoms 
      DO nCarbon = 2,121
         z(nCarbon)   = 0.0
      END DO

      ! Compute st deviation of Maxwell-Boltzmann distribution ( for the VELOCITY, not momenta!)
      SigmaCarbonVelocity = sqrt( Temperature / rmc )
      SigmaHydroVelocity  = sqrt( Temperature / rmh )

      ! Random velocities according to Maxwell-Boltzmann
      vxh = GaussianRandomNr( SigmaHydroVelocity ) 
      vyh = GaussianRandomNr( SigmaHydroVelocity ) 
      vzh = GaussianRandomNr( SigmaHydroVelocity ) 
      DO nCarbon = 1,121
         vzc(nCarbon) = GaussianRandomNr( SigmaCarbonVelocity ) 
      END DO

   END SUBROUTINE ThermalEquilibriumConditions

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
   SUBROUTINE VelocityVerlet( LangevinFriction )
      IMPLICIT NONE

         LOGICAL, INTENT(IN) :: LangevinFriction

         ! (1) FULL TIME STEP FOR THE POSITIONS

         xh=xh+vxh*dt+0.5*axh*(dt**2)
         yh=yh+vyh*dt+0.5*ayh*(dt**2)
         zh=zh+vzh*dt+0.5*azh*(dt**2)

         DO nn=1,nevo
            z(nn)=z(nn)+vzc(nn)*dt+0.5*azc(nn)*(dt**2)
         END DO

         ! (2) HALF TIME STEP FOR THE VELOCITIES

         IF ( .NOT. ( LangevinFriction ) ) THEN
            vxh=vxh+0.5*axh*dt
            vyh=vyh+0.5*ayh*dt
            vzh=vzh+0.5*azh*dt
            DO nn = 1,nevo
               vzc(nn)=vzc(nn)+0.5*azc(nn)*dt
            END DO
         ELSE IF ( LangevinFriction ) THEN
            vxh= 0.5 * axh * dt + LangIntegr * vxh
            vyh= 0.5 * ayh * dt + LangIntegr * vyh
            vzh= 0.5 * azh * dt + LangIntegr * vzh
            DO nn = 1,nevo
               vzc(nn)= 0.5 * azc(nn) * dt + LangIntegr * vzc(nn)
            END DO
         END IF

         ! (3) NEW FORCES AND ACCELERATIONS 

         CALL ComputePotential( xh, yh, zh, z, vv, axh, ayh, azh, azc )

         IF ( LangevinFriction ) THEN
!             axh = axh + GaussianRandomNr( GammaNoiseH ) / rmh
!             ayh = ayh + GaussianRandomNr( GammaNoiseH ) / rmh
!             azh = azh + GaussianRandomNr( GammaNoiseH ) / rmh
            axh = axh + UniformRandomNr( -sqrt(3.)*GammaNoiseH, sqrt(3.)*GammaNoiseH ) / rmh
            ayh = ayh + UniformRandomNr( -sqrt(3.)*GammaNoiseH, sqrt(3.)*GammaNoiseH ) / rmh
            azh = azh + UniformRandomNr( -sqrt(3.)*GammaNoiseH, sqrt(3.)*GammaNoiseH ) / rmh
            DO nn = 1,nevo
!                azc(nn) = azc(nn) + GaussianRandomNr( GammaNoiseC ) / rmc 
               azc(nn) = azc(nn) + UniformRandomNr( -sqrt(3.)*GammaNoiseC, sqrt(3.)*GammaNoiseC ) / rmc 
            END DO
         END IF

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES

         IF ( .NOT. LangevinFriction ) THEN
            vxh=vxh+0.5*axh*dt
            vyh=vyh+0.5*ayh*dt
            vzh=vzh+0.5*azh*dt
            DO nn = 1,nevo
               vzc(nn)=vzc(nn)+0.5*azc(nn)*dt
            END DO
         ELSE IF ( LangevinFriction ) THEN
            vxh= 0.5 * axh * dt + LangIntegr * vxh
            vyh= 0.5 * ayh * dt + LangIntegr * vyh
            vzh= 0.5 * azh * dt + LangIntegr * vzh
            DO nn = 1,nevo
               vzc(nn)= 0.5 * azc(nn) * dt + LangIntegr * vzc(nn)
            END DO
         END IF

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

END PROGRAM JK6


      




