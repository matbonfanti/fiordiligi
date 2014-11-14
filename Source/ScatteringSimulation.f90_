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
!>  \todo            everything
!>                 
!***************************************************************************************
MODULE ScatteringSimulation
   USE MyConsts
   USE ErrorTrap
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
!    USE ClassicalEqMotion
!    USE IndependentOscillatorsModel
!    USE RandomNumberGenerator


!    ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR A SCATTERING CALCULATION ONLY
!    
!    REAL    :: ezh    !< Initial kinetic energy of the scattering H atom
!    REAL    :: zhi    !< Initial Z position of the scattering H atom
!    INTEGER :: irho   !< Nr of rho values to sample
!    REAL    :: delrho !< Grid spacing in rho

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Scattering_ReadInput, Scattering_Initialize, Scattering_Run, Scattering_Dispose


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


!       CALL SetFieldFromInput( InputData, "ezh", ezh )
!       ezh = ezh / MyConsts_Hartree2eV
!       CALL SetFieldFromInput( InputData, "zhi", zhi )
!       zhi = zhi / MyConsts_Bohr2Ang
!       CALL SetFieldFromInput( InputData, "irho", irho )
!       CALL SetFieldFromInput( InputData, "delrho", delrho )
!       delrho = delrho / MyConsts_Bohr2Ang

   END SUBROUTINE Scattering_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the scattering simulation:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Initialize()
      IMPLICIT NONE

!          ! Allocation of position, velocity, acceleration arrays
!          ALLOCATE( X(nevo+3), V(nevo+3), A(nevo+3), APre(nevo+3) )
!          
!          ! if XYZ files of the trajectories are required, allocate memory to store the traj
!          IF ( PrintType >= FULL  .AND. ( RunType == SCATTERING .OR. RunType == EQUILIBRIUM ) ) THEN
!                ALLOCATE( Trajectory( 16, ntime ) )
!          END IF
! 
!          ! Allocate and define masses
!          ALLOCATE( MassVector( nevo + 3 ) )
!          MassVector = (/ (rmh, iCoord=1,3), (rmc, iCoord=1,nevo) /)
! 
!          ! Set variables for EOM integration in the microcanonical ensamble
!          CALL EvolutionSetup( MolecularDynamics, nevo+3, MassVector, dt )
! 
!          IF ( DynamicsGamma /= 0.0 .AND. RunType /= OSCIBATH_EQUIL .AND.  RunType /= CHAINBATH_EQUIL  ) THEN 
!             PRINT "(/,A,/)", " Setting langevin atoms at the border of the slab "
!             ! Set canonical dynamics at the borders of the carbon slab
!             LangevinSwitchOn = .TRUE.
!             LangevinSwitchOn( 1: MIN( 73, nevo )+3 ) = .FALSE.
!             CALL SetupThermostat( MolecularDynamics, DynamicsGamma, temp, LangevinSwitchOn(1:nevo+3) )
!          END IF
! 
!          ! Set variables for EOM integration with Langevin thermostat
!          CALL EvolutionSetup( Equilibration, nevo+3, MassVector, EquilTStep )
!          CALL SetupThermostat( Equilibration, Gamma, temp, (/ (.FALSE., iCoord=1,4), (.TRUE., iCoord=1,nevo-1) /) )


   END SUBROUTINE Scattering_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the scattering simulation and compute trapping probability over time.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Run()
      IMPLICIT NONE

! !***************************************************************************************************
! !                  H-GRAPHITE SCATTERING SIMULATION
! !***************************************************************************************************
! 
!    SUBROUTINE ScatteringSimulation()
!       IMPLICIT NONE
!       REAL, DIMENSION(121) :: vzsum, zsum
!       REAL, DIMENSION(:,:), ALLOCATABLE :: vz2av, zcav
!       REAL :: vav, vsqav, zav, zsqav
! 
!       !*************************************************************
!       !       AVERAGES VARIABLES ALLOCATION AND INITIALIZATION
!       !*************************************************************
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
!       ! Initialize random generation of initial conditions for scattering simulation
!       call random_seed()
! 
!       ! scan over the impact parameter... 
!       DO jRho = 1,irho+1
! 
!          ! Set impact parameter and print message
!          ImpactPar = float(jRho-1)*delrho+0.000001
! 
!          PRINT "(2/,A)",    "***************************************************"
!          PRINT "(A,F10.5)", "       IMPACT PARAMETER = ", ImpactPar
!          PRINT "(A,/)" ,    "***************************************************"
! 
!          PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "
! 
!          !run inum number of trajectories at the current rxn parameter
!          DO iTraj = 1,inum
!          
!             PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 
! 
!             !*************************************************************
!             ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
!             !*************************************************************
! 
!             ! Set initial conditions
!             CALL ScatteringConditions( X, V, ImpactPar, zhi, ezh, temp, rmh, rmc )
! 
!             !*************************************************************
!             ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
!             !*************************************************************
! 
!             ! increment average values for the initial lattice condition
!             DO iCarbon = 4,nevo+3
!                zav   = zav   + X(iCarbon)
!                zsqav = zsqav + X(iCarbon)**2
!                vav   = vav   + V(iCarbon)
!                vsqav = vsqav + V(iCarbon)**2
!             END DO
! 
!             DO iCarbon = 1,121
!                zsum(iCarbon)  = 0.0
!                vzsum(iCarbon) = 0.0
!             END DO
! 
!             ! Compute starting potential and forces
!             PotEnergy = VHSticking( X, A )
!             A(1:3) = A(1:3) / rmh
!             A(4:nevo+3) = A(4:nevo+3) / rmc
!             ! compute kinetic energy and total energy
!             KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
!             TotEnergy = PotEnergy + KinEnergy
!             IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))
! 
!             !*************************************************************
!             !         TIME EVOLUTION OF THE TRAJECTORY
!             !*************************************************************
!  
!             ! initialize counter for printing steps
!             kstep=0
! 
!             ! cycle over nstep velocity verlet iterations
!             DO iStep = 1,nstep
! 
!                ! Propagate for one timestep only if in the interaction region
!                IF ( X(3) <= zhi )  CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, PotEnergy )
! 
!                ! Compute kin energy and temperature
!                KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
!                TotEnergy = PotEnergy + KinEnergy
!                IstTemperature = 2.0*KinEnergy/(MyConsts_K2AU*(nevo+3))
! 
!                ! compute average values of the carbon positions
!                DO iCarbon = 1,nevo
!                   zsum(iCarbon)  = zsum(iCarbon)  + X(iCarbon+3)
!                   vzsum(iCarbon) = vzsum(iCarbon) + V(iCarbon+3)**2
!                END DO
! 
!                ! every nprint steps, compute trapping and write output
!                IF ( mod(iStep,nprint) == 0 ) THEN
! 
!                   ! increment counter for printing steps
!                   kstep=kstep+1
! 
!                   ! check if H is in the trapping region
!                   IF ( X(3) <= AdsorpLimit ) ptrap(jRho,kstep) = ptrap(jRho,kstep)+1.0
! 
!                   ! average lattice velocity and positions
!                   CALL lattice( kstep, zsum, vzsum, zcav, vz2av  )
! 
!                   ! Store the trajectory for XYZ printing
!                   IF ( PrintType >= FULL ) THEN
!                         Trajectory( :, kstep ) = X(1:16)
!                         NrOfTrajSteps = kstep
!                   END IF
! 
!                END IF 
! 
!             END DO
! 
!             ! At the end of the propagation, write the xyz file of the trajectory, if requested
!             IF ( PrintType >= FULL ) THEN
! 
!                ! print the full trajectory as output file
!                WRITE(OutFileName,"(A,I4.4,A,I4.4,A)") "Traj_",jRho,"_",iTraj,".xyz"
!                CALL WriteTrajectoryXYZ( Trajectory(:,1:NrOfTrajSteps)*MyConsts_Bohr2Ang, OutFileName, &
!                                              GraphiteLatticeConstant()*MyConsts_Bohr2Ang )
! 
!             END IF
! 
!          END DO
!       
!          PRINT "(A)"," Done! "
! 
!       END DO
! 
!       ! pritn info about the average for the initial lattice condition
!       zav=zav/float(nevo*inum*(irho+1))
!       zsqav=zsqav/float(nevo*inum*(irho+1))
!       vav=vav/float(nevo*inum*(irho+1))
!       vsqav=vsqav/float(nevo*inum*(irho+1))
! 
!       WRITE(*,"(/,A,/)") " * Average values for the initial lattice coordinates and velocities "
!       write(*,*) "average init zc=",zav * MyConsts_Bohr2Ang, " Ang"
!       write(*,*) "average init Vc=",0.5*zsqav*CarbonForceConstant() / MyConsts_K2AU,' K'
!       write(*,*) "average init vc=",vav, " au of velocity"
!       write(*,*) "average init Kc=",0.5*vsqav*rmc / MyConsts_K2AU,' K'
!       WRITE(*,"(/,/)")
! 
!       DO kk=1,10
!          DO iStep=1,ntime
!             zcav(kk,iStep)=zcav(kk,iStep)/float(nprint*inum*(irho+1))
!             vz2av(kk,iStep)=vz2av(kk,iStep)/float(nprint*inum*(irho+1))
!             vz2av(kk,iStep)=0.5*rmc*vz2av(kk,iStep) / MyConsts_K2AU
!          END DO
!       END DO
! 
!       PRINT*, " "
!       DO iStep=1,ntime
!          write(6,435) dt*real(iStep*nprint)/MyConsts_fs2AU, ( zcav(kk,iStep),  kk=1,10 )
!       END DO
! 
!       PRINT*, " "
!       DO iStep=1,ntime
!          write(6,436) dt*real(iStep*nprint)/MyConsts_fs2AU, ( vz2av(kk,iStep), kk=1,10 )
!       END DO
! 
!       ! write trapping data to mathematica file
!       open(unit=15,file='ptrap.nb',status='replace')
!       write(15,'(A6)') 'data={'
!       DO jRho=1,irho+1
!          write(15,'(A1)') '{'
!          DO iStep=1,ntime
!             IF(iStep.eq.ntime)THEN
!                write(15,'(F9.6)') ptrap(jRho,iStep)/inum
!             END IF
!             IF(iStep.ne.ntime)THEN 
!                write(15,'(F9.6,A1)') ptrap(jRho,iStep)/inum,','
!             END IF
!          END DO
!          IF(jRho.ne.irho+1)THEN
!             write(15,'(A2)')'},'
!          END IF
!          IF(jRho.eq.irho+1)THEN
!             write(15,'(A2)') '};'
!          END IF
!       END DO
!       write(15,505)'ListPlot3D[data,ColorFunction->Hue,AxesLabel->', &
!             '{"Time(ps)","Reaction Parameter(b)","Ptrap(t,b)"},Boxed->False]'
!       CLOSE(15)
! 
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
! 
!       435 FORMAT(1f6.3,10f7.3)
!       436 FORMAT(1f6.3,10f7.2)
!       505 FORMAT(A46,A63)
!       605 FORMAT(f8.4,3x,f10.5)
! 
!    END SUBROUTINE ScatteringSimulation


   END SUBROUTINE Scattering_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the simulation.
!>
!*******************************************************************************
   SUBROUTINE Scattering_Dispose()
      IMPLICIT NONE


   END SUBROUTINE Scattering_Dispose


END MODULE ScatteringSimulation
