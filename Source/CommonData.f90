!***************************************************************************************
!*                              MODULE CommonData
!***************************************************************************************
!
!>  \brief     Common data
!>  \details   This module include the data to be shared by the
!>             all the other modules of the code
!
!***************************************************************************************
MODULE CommonData
   IMPLICIT NONE

   PUBLIC          ! This module contain data that is supposed to be shared

   ! PARAMETERS 
   
   ! Variable to define adsorption condition
   REAL, PARAMETER :: AdsorpLimit = 2.8

   ! VARIABLES SET FROM INPUT
   
   !> Variable to define which kind of calculation is required 
   INTEGER :: RunType
   INTEGER, PARAMETER :: SCATTERING = 1, &  ! Scattering calculation with H coming from gas-phase
                         EQUILIBRIUM = 2    ! Equilibrium calculation with H already adsorbed

   !> Variable to set the print level of the calculation
   INTEGER :: PrintType
   INTEGER, PARAMETER :: EQUILIBRDBG = 4,   &   ! fully detailed information about the trajs
                         DEBUG = 3,   &   ! fully detailed information about the trajs
                         FULL = 2,    &   ! files to plot the make animations, averages for each traj
                         MINIMAL = 1      ! minimal level of output, only final averages

   ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR ALL THE KIND OF CALCULATIONS
   
   REAL    :: dt     !< Time step for the integration of the classical EOM
   REAL    :: temp   !< Temperature of the surface (and of H when adsorbed)
   INTEGER :: nstep  !< Total nr of time step per trajectory
   INTEGER :: nprint !< Nr of time steps between each printing of propagation info
   INTEGER :: inum   !< Nr of trajectories ( per rho value, in case of scattering simulation )  
   INTEGER :: nevo   !< Nr of evolving carbon atoms in the slab ( ranging from 4 to 121)
   REAL    :: rmh    !< Mass of the hydrogen atom
   REAL    :: rmc    !< Mass of the carbon atoms
   INTEGER :: ntime  !< Nr of analysis steps ( = nstep / nprint )
   
   ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR A SCATTERING CALCULATION ONLY
   
   REAL    :: ezh    !< Initial kinetic energy of the scattering H atom
   REAL    :: zhi    !< Initial Z position of the scattering H atom
   INTEGER :: irho   !< Nr of rho values to sample
   REAL    :: delrho !< Grid spacing in rho
      
   ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR AN EQUILIBRIUM CALCULATION ONLY
   
   REAL    :: Gamma                !< Friction parameter of the Langevin equation
   INTEGER :: NrEquilibSteps       !< Nr of time step of the equilibration
   INTEGER :: EquilibrInitAverage  !< Nr of time step to skip in computing averages during equilibration
   
   ! POSITION, VELOCITY, ACCELERATION 

   REAL, DIMENSION(:), SAVE, ALLOCATABLE :: X    ! Position at given timestep
   REAL, DIMENSION(:), SAVE, ALLOCATABLE :: V    ! Velocity at given timestep
   REAL, DIMENSION(:), SAVE, ALLOCATABLE :: A    ! Acceleration at ginve timestep
   REAL, DIMENSION(:), SAVE, ALLOCATABLE :: APre ! Acceleration from the previous time step

   ! ISTANTANEOUS PROPERTIES

   REAL  :: KinEnergy                   ! Kinetic energy
   REAL  :: PotEnergy                   ! Potential energy
   REAL  :: TotEnergy                   ! Total energy
   REAL  :: IstTemperature              ! Istantaneous temperature
   REAL  :: GraphitePlaneZ              ! Estimate on the reference Z for the slab

   ! TRAJECTORY AVERAGES

   INTEGER :: NrOfStepsAver             ! Number of steps of the average computation
   REAL  :: TempAverage                 ! To accumulate average temperature over time
   REAL  :: TempVariance                ! To accumulate squared temperature over time
   REAL, DIMENSION(3) :: HPosAverage    ! To accumulate average displacement of the H atom
   REAL, DIMENSION(3) :: HPosVariance   ! To accumulate squared displacement of the H atom
   REAL  :: C1Average                   ! To accumulate average displacement of the C1 atom
   REAL  :: C1Variance                  ! To accumulate squared displacement of the C1 atom

   ! OTHER DATA
     
   REAL, DIMENSION(:,:), ALLOCATABLE :: Trajectory   ! Real array to store a trajectory for XYZ printing
   INTEGER :: NrOfTrajSteps                          ! Nr of snapshot in the trajectory
   REAL :: ImpactPar                                 ! impact parameter of H in scattering calculations
   REAL, DIMENSION(:,:), ALLOCATABLE :: ptrap        ! Variable to store trapping probability, resolved in time and in rho

   


   !********************************************************************
   
   
END MODULE CommonData