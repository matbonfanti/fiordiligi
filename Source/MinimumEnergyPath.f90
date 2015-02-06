!***************************************************************************************
!*                              MODULE MinimumEnergyPath
!***************************************************************************************
!
!>  \brief     Subroutines for computing the MEP of the H-carbon sticking PES
!>  \details   This module contains subroutines to compute and plot the minimum \n
!>             energy path of the potential in different setup conditions
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             4 February 2015
!>
!***************************************************************************************
MODULE MinimumEnergyPath
#include "preprocessoptions.cpp"
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE UnitConversion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE PrintTools
!    USE SplineInterpolator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: MinimumEnergyPath_ReadInput, MinimumEnergyPath_Initialize, MinimumEnergyPath_Run, MinimumEnergyPath_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

   !> Restart file for full D PES transition state
   CHARACTER(100) :: FullPES_TSRestartFile
   !> TS Restart is enabled
   LOGICAL :: TSRestart = .FALSE.

   ! Parameters for optimizations
   INTEGER :: MaxOptSteps  !< Maximum number of steps for optimizations
   REAL :: OptThreshold    !< Convergence criterium for optimizations

   ! Parameters for Minimum Energy Path computation
   INTEGER :: MaxMEPNrSteps  !< Max number of MEP steps
   REAL    :: MEPStep        !< Size of each MEP step

   ! Arrays to store info about MEPS
   REAL, DIMENSION(:,:), ALLOCATABLE :: MEP_2D_Plus,   MEP_2D_Minus
   REAL, DIMENSION(:,:), ALLOCATABLE :: MEP_Full_Plus, MEP_Full_Minus

   !> Data for printing trajectory in VTK format
   TYPE(VTKInfo) :: MEPTraj_2D_Plus, MEPTraj_2D_Minus, MEPTraj_Full_Plus, MEPTraj_Full_Minus

!    !> Parameters to plot potential cuts
!    REAL :: GridSpacing          !< Spacing of the grid to plot the potential
!    INTEGER :: NGridPoint        !< Nr of points around the equilibrium value
! 
!    !> Parameters to plot 2D potential
!    REAL, DIMENSION(:), ALLOCATABLE :: PotentialArray                     !< array to store 2D potential
!    REAL, DIMENSION(:), ALLOCATABLE :: ZCArray, ZHArray, QArray, OptZC    !< array with coordinates grid
!    REAL :: ZHmin, ZHmax, ZCmin, ZCmax, Qmin, Qmax                        !< boundaries of the grid
!    INTEGER :: NpointZH, NpointZC, NpointQ                                !< nr of points of the grid
!    TYPE(VTKInfo) :: PotentialCH
!    TYPE(VTKInfo) :: CouplingV
!    TYPE(SplineType) :: SplineData                           !< Data for spline interpolation

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific for the
!> analysis of the potential
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE MinimumEnergyPath_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData

      ! Set variables for Transition State restart
      CALL SetFieldFromInput( InputData, "FullPES_TSRestartFile", FullPES_TSRestartFile, "NO_RESTART" )
      IF ( TRIM(FullPES_TSRestartFile) /= "NO_RESTART" ) TSRestart = .TRUE.

      ! Set variables for minima and TS optimization
      CALL SetFieldFromInput( InputData, "MaxOptSteps", MaxOptSteps, 10**6 )
      CALL SetFieldFromInput( InputData, "OptThreshold", OptThreshold, 1.E-6 )

      ! Set variables for Minimum Energy Path computation
      CALL SetFieldFromInput( InputData, "MaxMEPNrSteps", MaxMEPNrSteps, 1000 )
      CALL SetFieldFromInput( InputData, "MEPStep", MEPStep, 0.001 )
      MEPStep = MEPStep * LengthConversion(InputUnits, InternalUnits)

   END SUBROUTINE MinimumEnergyPath_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the analysis of the Potential:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE MinimumEnergyPath_Initialize()
      IMPLICIT NONE
      INTEGER :: iCoord

      ! Allocate memory and initialize vectors for positions, forces and masses

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         ALLOCATE( X(3+NCarbon), A(3+NCarbon), MassVector(3+NCarbon) )
         MassVector = (/ (MassH, iCoord=1,3), (MassC, iCoord=1,NCarbon) /)

      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ALLOCATE( X(4+NBath), A(4+NBath), MassVector(4+NBath) )
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,NBath) /)

      ELSE IF ( BathType ==  DOUBLE_CHAIN ) THEN
         ALLOCATE( X(4+2*NBath), A(4+2*NBath), MassVector(4+2*NBath) )
         MassVector = (/ (MassH, iCoord=1,3), MassC, (MassBath, iCoord=1,2*NBath) /)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ALLOCATE( X(4), A(4), MassVector(4) )
         MassVector = (/ (MassH, iCoord=1,3), MassC /)
      END IF

      ! store the nr of dimensions
      NDim = size(X)

   END SUBROUTINE MinimumEnergyPath_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the Potential Analysis.
!>
!*******************************************************************************
   SUBROUTINE MinimumEnergyPath_Run()
      IMPLICIT NONE
      INTEGER :: MEP2DUnit, MEPFullUnit, TS2DUnit, TSFullUnit

      REAL, DIMENSION(NDim)      :: XMin, XAsy, XTS1, XTS2, MEPX
      REAL :: Emin, Easy, ETS1, ETS2, EMEP
      LOGICAL, DIMENSION(NDim)   :: OptMask
      INTEGER :: iStep, Max2DPlus, Max2DMinus, MaxFullPlus, MaxFullMinus
      LOGICAL :: Check

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "            MINIMUM ENERGY PATH ANALYSIS  "
      PRINT "(A,/)" ,    "***************************************************"

      ! ===================================================================================================
      ! (1) identify the global minimum and the asymptotic geometry of the PES
      ! ===================================================================================================

      PRINT "(/,A)",    " **** Computing the global minumim of the PES ****"

      ! guess reasonable coordinates of the minimum of the PES
      X(1:2) = 0.0
      X(3) = HZEquilibrium
      X(4) = C1Puckering

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         X(5:NDim) = MinSlab(1:NDim-4)
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH .OR. BathType ==  DOUBLE_CHAIN ) THEN
         X(5:NDim) = 0.0
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! nothing to set
      END IF
 
      ! Find minimum by Newton's optimization
      XMin = SteepLocator( X, MaxOptSteps, OptThreshold )
      XMin(3:) = XMin(3:) - (XMin(5)+XMin(6)+XMin(7))/3.0
      ! Computing the energy at this geometry
      EMin = HStickPotential( XMin, A )

      PRINT "(2/,A)",    " **** Computing the asymptotic configuration of the PES ****"

      ! guess reasonable coordinates of the asymptote of the PES
      X(1:2) = 0.0
      X(3) = 100.0
      X(4) = 0.0

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         X(5:NDim) = 0.0
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH .OR. BathType ==  DOUBLE_CHAIN ) THEN
         X(5:NDim) = 0.0
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! nothing to set
      END IF
 
      ! Find minimum by Newton's optimization
      XAsy = SteepLocator( X, MaxOptSteps, OptThreshold )
      ! Computing the energy at this geometry
      Easy = HStickPotential( XAsy, A )

      ! Write info on geometry and energy as output
      WRITE(*,500) EAsy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),              &
                   EMin*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),              &
                   (EMin-EAsy)*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),       &
                   XMin(3)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits),           &
                   XMin(4)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits),           &
                   (XMin(3)-XMin(4))*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits) 


      ! ===================================================================================================
      ! (2) look for the transition state of the potential energy surface 
      ! ===================================================================================================


      PRINT "(2/,A)",    " **** Transition state of 2D (zH,zC) PES, bath in asymptotic geometry **** "

      ! guess reasonable coordinates of the TS of the PES
      X(:) = XAsy(:); X(3) = 3.4; X(4) = 0.190

      ! Set appropriate constraints and find TS by Newton's optimization
      OptMask = .FALSE.; OptMask(3:4) = .TRUE.
      XTS1 = NewtonLocator( X, MaxOptSteps, OptThreshold, OptThreshold, Mask=OptMask, TransitionState=.TRUE. )
      ETS1 = HStickPotential( XTS1, A )

      ! Print results to screen
      WRITE(*,501) XTS1(3)*LengthConversion(InternalUnits, InputUnits), LengthUnit(InputUnits), &
                   XTS1(4)*LengthConversion(InternalUnits, InputUnits), LengthUnit(InputUnits), &
                   ETS1*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),    &
                   (ETS1-EAsy)*EnergyConversion(InternalUnits, InputUnits),  EnergyUnit(InputUnits)

      ! Write to output file the coordinates of the TS
      TS2DUnit = LookForFreeUnit()
      OPEN( FILE="TS_Reduced2DPES.dat", UNIT=TS2DUnit )
      WRITE( TS2DUnit, * ) XTS1
      CLOSE( TS2DUnit )

      PRINT "(2/,A)",    " **** Transition state of full PES **** "

      ! Set appropriate constraints and find TS by Newton's optimization
      OptMask = .TRUE.
      OptMask(1:4) = (/ .FALSE., .FALSE., .TRUE., .TRUE. /)
      IF ( BathType ==  SLAB_POTENTIAL ) OptMask(5:7) = .FALSE.

      ! Set starting conditions
      IF ( TSRestart ) THEN
         TSFullUnit = LookForFreeUnit()
         OPEN( FILE=TRIM(FullPES_TSRestartFile), UNIT=TSFullUnit )
         READ( TSFullUnit, * ) XTS2
         XTS2 = XTS2*LengthConversion(InputUnits,InternalUnits)
         CLOSE( TSFullUnit )
      ELSE
         XTS2 = XTS1
      END IF

      ! Run optimization
      XTS2 = NewtonLocator( XTS2, MaxOptSteps, OptThreshold, OptThreshold, Mask=OptMask, TransitionState=.TRUE. )
      ETS2 = HStickPotential( XTS2, A )

      ! Print results to screen
      WRITE(*,501) XTS2(3)*LengthConversion(InternalUnits, InputUnits), LengthUnit(InputUnits), &
                   XTS2(4)*LengthConversion(InternalUnits, InputUnits), LengthUnit(InputUnits), &
                   ETS2*EnergyConversion(InternalUnits, InputUnits), EnergyUnit(InputUnits),    &
                   (ETS2-EAsy)*EnergyConversion(InternalUnits, InputUnits),  EnergyUnit(InputUnits)

      ! Write to output file the coordinates of the TS
      TSFullUnit = LookForFreeUnit()
      OPEN( FILE="TS_FullPES.dat", UNIT=TSFullUnit )
      WRITE( TSFullUnit, * ) XTS2*LengthConversion(InternalUnits, InputUnits)
      CLOSE( TSFullUnit )

      ! ===================================================================================================
      ! (3) compute the minimum energy path from the TS
      ! ===================================================================================================

      ! Allocate memory to store the MEPS
      ALLOCATE( MEP_2D_Plus(3,MaxMEPNrSteps), MEP_2D_Minus(3,MaxMEPNrSteps),  &
                MEP_Full_Plus(3,MaxMEPNrSteps), MEP_Full_Minus(3,MaxMEPNrSteps)    )

      PRINT "(/,A,/)",    " **** MEP from the 2D (zH,zC) PES TS **** "

      ! Define constraints mask
      OptMask = .FALSE.; OptMask(3:4) = .TRUE.

      ! Write header lines screen 
      WRITE(*,700) TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), TRIM(EnergyUnit(InputUnits))

      ! Write to output file and to screen the starting point
      WRITE(*,601) 0, XTS1(3)*LengthConversion(InternalUnits, InputUnits), &
         XTS1(4)*LengthConversion(InternalUnits, InputUnits), ETS1*EnergyConversion(InternalUnits, InputUnits)

      ! First step from the TS, with plus sign
      MEPX = XTS1
      Check = DoMEPStep( MEPX, +0.01, FirstStep=.TRUE., Mask=OptMask )
      MEP_2D_Plus(:,1) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)

      DO iStep = 2, MaxMEPNrSteps
         ! Following steps
         Check = DoMEPStep( MEPX, MEPStep, FirstStep=.FALSE., Mask=OptMask )
         MEP_2D_Plus(:,iStep) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)
         IF ( .NOT. Check )  EXIT
      END DO
      Max2DPlus = iStep - 1

      ! Write to screen final step
      WRITE(*,601) Max2DPlus, MEPX(3)*LengthConversion(InternalUnits, InputUnits), &
                   MEPX(4)*LengthConversion(InternalUnits, InputUnits), &
                   MEP_2D_Plus(3,Max2DPlus)*EnergyConversion(InternalUnits, InputUnits)

      ! First step from the TS, with minus sign
      MEPX = XTS1
      Check = DoMEPStep( MEPX, -0.01, FirstStep=.TRUE., Mask=OptMask )
      MEP_2D_Minus(:,1) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)

      DO iStep = 2, MaxMEPNrSteps
         ! Following steps
         Check = DoMEPStep( MEPX, MEPStep, FirstStep=.FALSE., Mask=OptMask )
         MEP_2D_Minus(:,iStep) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)
         IF ( .NOT. Check ) EXIT
      END DO
      Max2DMinus = iStep - 1

      ! Write to screen final step
      WRITE(*,601) -Max2DMinus, MEPX(3)*LengthConversion(InternalUnits, InputUnits), &
                   MEPX(4)*LengthConversion(InternalUnits, InputUnits), &
                   MEP_2D_Minus(3,Max2DMinus)*EnergyConversion(InternalUnits, InputUnits)


      ! Open output file and write path to file
      MEP2DUnit = LookForFreeUnit()
      OPEN( FILE="MEP_Reduced2DPES.dat", UNIT=MEP2DUnit )

      ! Header of the output file
      WRITE(MEP2DUnit,700) TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), TRIM(EnergyUnit(InputUnits))

      ! Write minus section of the MEP
      DO iStep = Max2DMinus, 1, -1
         WRITE(MEP2DUnit,601) -iStep, MEP_2D_Minus(1,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_2D_Minus(2,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_2D_Minus(3,iStep)*EnergyConversion(InternalUnits, InputUnits)
      END DO
      ! Write TS properties
      WRITE(MEP2DUnit,601) 0, XTS1(3)*LengthConversion(InternalUnits, InputUnits), &
         XTS1(4)*LengthConversion(InternalUnits, InputUnits), ETS1*EnergyConversion(InternalUnits, InputUnits)
      ! Write plus section of the MEP
      DO iStep = 1, Max2DPlus
         WRITE(MEP2DUnit,601) iStep, MEP_2D_Plus(1,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_2D_Plus(2,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_2D_Plus(3,iStep)*EnergyConversion(InternalUnits, InputUnits)
      END DO

      ! Close file
      CLOSE(MEP2DUnit)

      PRINT "(/,A,/)",    " **** MEP from the global PES TS **** "

      ! Define constraints mask
      OptMask = .TRUE.; OptMask(1:4) = (/ .FALSE., .FALSE., .TRUE., .TRUE. /)
      IF ( BathType ==  SLAB_POTENTIAL ) OptMask(5:7) = .FALSE.

      ! Write header lines screen 
      WRITE(*,700) TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), TRIM(EnergyUnit(InputUnits))

      ! Write to output file and to screen the starting point
      WRITE(*,601) 0, XTS2(3)*LengthConversion(InternalUnits, InputUnits), &
         XTS2(4)*LengthConversion(InternalUnits, InputUnits), ETS2*EnergyConversion(InternalUnits, InputUnits)

      ! First step from the TS, with plus sign
      MEPX = XTS2
      Check = DoMEPStep( MEPX, +0.01, FirstStep=.TRUE., Mask=OptMask )
      MEP_Full_Plus(:,1) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)

      DO iStep = 2, MaxMEPNrSteps
         ! Following steps
         Check = DoMEPStep( MEPX, MEPStep, FirstStep=.FALSE., Mask=OptMask )
         MEP_Full_Plus(:,iStep) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)
         IF ( .NOT. Check )  EXIT
      END DO
      MaxFullPlus = iStep - 1

      ! Write to screen final step
      WRITE(*,601) MaxFullPlus, MEPX(3)*LengthConversion(InternalUnits, InputUnits), &
                   MEPX(4)*LengthConversion(InternalUnits, InputUnits), &
                   MEP_Full_Plus(3,Max2DPlus)*EnergyConversion(InternalUnits, InputUnits)

      ! First step from the TS, with minus sign
      MEPX = XTS2
      Check = DoMEPStep( MEPX, -0.01, FirstStep=.TRUE., Mask=OptMask )
      MEP_Full_Minus(:,1) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)

      DO iStep = 2, MaxMEPNrSteps
         ! Following steps
         Check = DoMEPStep( MEPX, MEPStep, FirstStep=.FALSE., Mask=OptMask )
         MEP_Full_Minus(:,iStep) = (/ MEPX(3), MEPX(4), HStickPotential( MEPX, A ) /)
         IF ( .NOT. Check ) EXIT
      END DO
      MaxFullMinus = iStep - 1

      ! Write to screen final step
      WRITE(*,601) -MaxFullMinus, MEPX(3)*LengthConversion(InternalUnits, InputUnits), &
                   MEPX(4)*LengthConversion(InternalUnits, InputUnits), &
                   MEP_Full_Minus(3,Max2DMinus)*EnergyConversion(InternalUnits, InputUnits)


      ! Open output file and write path to file
      MEPFullUnit = LookForFreeUnit()
      OPEN( FILE="MEP_FullPES.dat", UNIT=MEPFullUnit )

      ! Header of the output file
      WRITE(MEPFullUnit,700) TRIM(LengthUnit(InputUnits)), TRIM(LengthUnit(InputUnits)), TRIM(EnergyUnit(InputUnits))

      ! Write minus section of the MEP
      DO iStep = MaxFullMinus, 1, -1
         WRITE(MEPFullUnit,601) -iStep, MEP_Full_Minus(1,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_Full_Minus(2,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_Full_Minus(3,iStep)*EnergyConversion(InternalUnits, InputUnits)
      END DO
      ! Write TS properties
      WRITE(MEPFullUnit,601) 0, XTS2(3)*LengthConversion(InternalUnits, InputUnits), &
         XTS2(4)*LengthConversion(InternalUnits, InputUnits), ETS2*EnergyConversion(InternalUnits, InputUnits)
      ! Write plus section of the MEP
      DO iStep = 1, MaxFullPlus
         WRITE(MEPFullUnit,601) iStep, MEP_Full_Plus(1,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_Full_Plus(2,iStep)*LengthConversion(InternalUnits, InputUnits), &
                              MEP_Full_Plus(3,iStep)*EnergyConversion(InternalUnits, InputUnits)
      END DO

      ! Close file
      CLOSE(MEPFullUnit)
 
      ! Write the MEP to VTV files
      CALL VTK_WriteTrajectory ( MEPTraj_2D_Plus, &
         MEP_2D_Plus(1:2,1:Max2DPlus)*LengthConversion(InternalUnits, InputUnits), "MEPTraj_2D_Plus" )
      CALL VTK_WriteTrajectory ( MEPTraj_2D_Minus, &
         MEP_2D_Minus(1:2,1:Max2DMinus)*LengthConversion(InternalUnits, InputUnits),     "MEPTraj_2D_Minus" )
      CALL VTK_WriteTrajectory ( MEPTraj_Full_Plus, &
         MEP_Full_Plus(1:2,1:MaxFullPlus)*LengthConversion(InternalUnits, InputUnits),   "MEPTraj_Full_Plus" )
      CALL VTK_WriteTrajectory ( MEPTraj_Full_Minus, &
         MEP_Full_Minus(1:2,1:MaxFullMinus)*LengthConversion(InternalUnits, InputUnits), "MEPTraj_Full_Minus" )



   500 FORMAT (/, " H-Graphite adsorption geometry and energy   ",/  &
                  " * Asymptotic energy:           ",1F15.6,1X,A, /, &
                  " * Energy at the minimum:       ",1F15.6,1X,A, /, &
                  " * Adsorption energy:           ",1F15.6,1X,A, /, &
                  " * ZH at the minimum:           ",1F15.6,1X,A, /, &
                  " * ZC at the minimum:           ",1F15.6,1X,A, /, &
                  " * C-H equilibrium distance:    ",1F15.6,1X,A     )

   501 FORMAT (/, " H-Graphite 2D cut Transition State          ",/  &
                  " * ZH at the TS:                ",1F15.6,1X,A, /, &
                  " * ZC at the TS:                ",1F15.6,1X,A, /, &
                  " * Energy at the TS:            ",1F15.6,1X,A, /, &
                  " * Barrier energy:              ",1F15.6,1X,A  )

   502 FORMAT (/, " H-Graphite global PES Transition State      ",/  &
                  " * ZH at the TS:                ",1F15.6,1X,A, /, &
                  " * ZC at the TS:                ",1F15.6,1X,A, /, &
                  " * Energy at the TS:            ",1F15.6,1X,A, /, &
                  " * Barrier energy:              ",1F15.6,1X,A   )

   600 FORMAT ( A12,A20,A20,A20 )
   601 FORMAT ( I12,F20.6,F20.6,F20.6 )

   700 FORMAT ( "#  N Step   ", "   zH Coord / ", A6, "   zC Coord / ", A6, "     Energy / ", A6, /,  &
                "#---------------------------------------------------------------------------------" )
   701 FORMAT ( "#---------------------------------------------------------------------------------" )

   END SUBROUTINE MinimumEnergyPath_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the simulation.
!>
!*******************************************************************************
   SUBROUTINE MinimumEnergyPath_Dispose()
      IMPLICIT NONE

      ! Deallocate memory
      DEALLOCATE( X, A, MassVector )

   END SUBROUTINE MinimumEnergyPath_Dispose


!*************************************************************************************************

!*******************************************************************************
!> Compute the potential in the different available forms.
!>
!> @param Positions   Input vector with the coordinates where to compute V
!> @param Forces      Output vectors with the forces (=-derivatives of V)
!> @returns           The real value of the potential in Positions
!*******************************************************************************
   REAL FUNCTION HStickPotential( Positions, Forces ) RESULT(V)
      REAL, DIMENSION(NDim), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(NDim), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i    
      REAL    :: CoupFs, DTimesX

      ! Initialize forces and potential
      V         = 0.0
      Forces(:) = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         V = VHSticking( Positions, Forces )

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         V = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Coupling function of the system
         CoupFs = Positions(4)-C1Puckering
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, CoupFs, Positions(5:), V, Forces(4), Forces(5:) ) 

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Compute potential and forces of the system
         V = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Coupling function of the system
         CoupFs = Positions(4)-C1Puckering
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), CoupFs, Positions(5:NBath+4), V, Forces(4), Forces(5:NBath+4) ) 
         CALL BathPotentialAndForces( DblBath(2), CoupFs, Positions(NBath+5:2*NBath+4), V, Forces(4), Forces(NBath+5:2*NBath+4) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential using only the 4D subroutine
         V = VHFourDimensional( Positions, Forces )

      END IF

   END FUNCTION HStickPotential

!*************************************************************************************************

!*******************************************************************************
!> Compute the Hessian of the potential in the different available forms.
!>
!> @param AtPoint   Input vector with the coordinates where to compute H
!> @returns         Hessian matrix of the potential in AtPoint
!*******************************************************************************
   FUNCTION HessianOfThePotential( AtPoint ) RESULT( Hessian )
      REAL, DIMENSION(NDim,NDim)        :: Hessian
      REAL, DIMENSION(NDim), INTENT(IN) :: AtPoint
      REAL, DIMENSION(NBath) :: CouplingsCoeffs
      REAL    :: DistortionCoeff
      INTEGER :: i

      IF ( BathType == SLAB_POTENTIAL ) THEN 

         ! Compute hessian of the full slab potential
         Hessian(1:NDim,1:NDim) = HessianOfTheFullPotential( AtPoint(1:NDim), MassH, MassC )

      ELSE IF ( BathType == NORMAL_BATH ) THEN

         ! Compute hessian of the system
         Hessian(1:4,1:4) = HessianOfTheSystem( AtPoint(1:4), MassH, MassC )

         ! add the part of the hessian of the independent oscillator model
         CALL  HessianOfTheBath( Bath, Hessian(5:NDim, 5:NDim) ) 

         ! add the contribution from SYSTEM-BATH coupling and from DISTORTION CORRECTION
         CALL  CouplingAndDistortionHessian( Bath, CouplingsCoeffs, DistortionCoeff )
         DO i = 1, NBath
            Hessian(4,4+i) = Hessian(4,4+i) + CouplingsCoeffs(i) / SQRT( MassVector(4) * MassVector(4+i) ) 
            Hessian(4+i,4) = Hessian(4+i,4) + CouplingsCoeffs(i) / SQRT( MassVector(4) * MassVector(4+i) ) 
         END DO
         Hessian(4,4) = Hessian(4,4) + DistortionCoeff / MassVector(4)

      ELSE IF ( BathType == CHAIN_BATH ) THEN

         ! Compute hessian of the system
         Hessian(1:4,1:4) = HessianOfTheSystem( AtPoint(1:4), MassH, MassC )

         ! add the part of the hessian of the independent oscillator model
         CALL  HessianOfTheBath( Bath, Hessian(5:NDim, 5:NDim) ) 

         ! add the contribution from SYSTEM-BATH coupling and from DISTORTION CORRECTION
         CALL  CouplingAndDistortionHessian( Bath, CouplingsCoeffs, DistortionCoeff )
         Hessian(4,5) = Hessian(4,5) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(5,4) = Hessian(5,4) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(4,4) = Hessian(4,4) + DistortionCoeff / MassVector(4)

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN

         ! Compute hessian of the system
         Hessian(1:4,1:4) = HessianOfTheSystem( AtPoint(1:4), MassH, MassC )

         ! add the part of the hessian of the independent oscillator model
         CALL HessianOfTheBath( DblBath(1), Hessian(5:4+NBath, 5:4+NBath) ) 
         CALL HessianOfTheBath( DblBath(2), Hessian(5+NBath:NDim, 5+NBath:NDim) ) 

         ! add the contribution from SYSTEM-BATH coupling and from DISTORTION CORRECTION
         CALL  CouplingAndDistortionHessian( DblBath(1), CouplingsCoeffs, DistortionCoeff )
         Hessian(4,5) = Hessian(4,5) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(5,4) = Hessian(4,5) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(4,4) = Hessian(4,4) + DistortionCoeff / MassVector(4)
         CALL  CouplingAndDistortionHessian( DblBath(2), CouplingsCoeffs, DistortionCoeff )
         Hessian(4,5+NBath) = Hessian(4,5+NBath) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(5+NBath,4) = Hessian(5+NBath,4) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(4,4) = Hessian(4,4) + DistortionCoeff / MassVector(4)

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN

         ! Compute hessian of the system
         Hessian = HessianOfTheSystem( AtPoint, MassH, MassC )

      END IF

   END FUNCTION HessianOfThePotential
 
!*************************************************************************************************

! J. Chem. Phys., Vol. 121, No. 20, 22 November 2004
   FUNCTION NewtonLocator( StartX, NMaxIter, GradThresh, DisplThresh, Mask, TransitionState ) RESULT( StationaryPoint )
      IMPLICIT NONE
      REAL, DIMENSION(NDim), INTENT(IN) :: StartX
      INTEGER, INTENT(IN)               :: NMaxIter
      REAL, INTENT(IN)                  :: GradThresh, DisplThresh
      LOGICAL, INTENT(IN), OPTIONAL     :: TransitionState
      LOGICAL, DIMENSION(NDim), INTENT(IN), OPTIONAL :: Mask
      REAL, DIMENSION(NDim)             :: StationaryPoint

      REAL, DIMENSION(NDim) :: CurrentX, Forces
      REAL, DIMENSION(NDim,NDim) :: Hessian

      REAL, DIMENSION(:), ALLOCATABLE :: EigenValues, Factors, WrkForces
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenVectors, WrkHessian

      REAL :: V, DisplNorm, GradNorm, Factor
      INTEGER :: NIter, i, n, NOpt
      LOGICAL :: TSCheck, SteepestDescent

      ! Decide whether to look for a TS
      IF (PRESENT( TransitionState )) THEN
         TSCheck = TransitionState
      ELSE
         TSCheck = .FALSE.
      ENDIF

      ! Print info to screen
      WRITE(*,"(/,A,/)") " Stationary point locator with Newton's method "
      IF (.NOT. TSCheck) THEN
         WRITE(*,"(A)") " Looking for a minimum of the potential"
      ELSE
         WRITE(*,"(A)") " Looking for a first order saddle point"
      END IF

      ! Number of non constrained variables
      IF ( PRESENT(Mask) ) THEN
         NOpt = COUNT(Mask)
         WRITE(*,"(A,I3,A,/)") " Optimization of ",NOpt," variables "
      ELSE
         NOpt = NDim
         WRITE(*,"(A,/)") " All the variables will be optimized "
      END IF

      ! Allocate memory
      ALLOCATE( EigenValues(NOpt), Factors(NOpt), WrkForces(NOpt) )
      ALLOCATE( EigenVectors(NOpt, NOpt), WrkHessian(NOpt, NOpt) )

      ! Start at initial position
      CurrentX = StartX
 
      WRITE(*,"(A12,A20,A20)") "N Iteration", "Displacement Norm", "Gradient Norm"
      WRITE(*,*)              "-------------------------------------------------------------"

      Iterations: DO NIter = 1, NMaxIter

         IF ( PRESENT(Mask) ) THEN
            ! Compute constrained forces at current position
            V = HStickPotential( CurrentX, Forces )
            WrkForces = ConstrainedVector( Forces, Mask )
            
            ! Compute constrained Hessian at current position
            Hessian = HessianOfThePotential( CurrentX )
            WrkHessian = ConstrainedHessian( Hessian, Mask )
         ELSE
            ! compute Hessian and forces at current position
            V = HStickPotential( CurrentX, WrkForces )
            WrkHessian = HessianOfThePotential( CurrentX )
         END IF

         ! Compute norm of the gradient
         GradNorm  = SQRT(TheOneWithVectorDotVector(WrkForces, WrkForces))

         ! When the gradient is large, switch off newton and use gradient only
         IF ( GradNorm < 0.25E-04 ) THEN
            SteepestDescent = .TRUE.
         ELSE
            SteepestDescent = .FALSE.
         END IF

         ! Diagonalize Hessian to transform coords to normal modes
         CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
         ! Transform forces to normal modes
         WrkForces = TheOneWithMatrixVectorProduct( TheOneWithTransposeMatrix(EigenVectors), WrkForces )

         ! Weigths forces with eigenvalues (hence obtain displacements)
         IF (.NOT.SteepestDescent ) THEN
            DO i = 1, NOpt
               IF ( TSCheck ) THEN 
                  Factors(i) = 0.5 * ( ABS(EigenValues(i)) + SQRT( EigenValues(i)**2 + 4.0 * WrkForces(i)**2  ) )
               ELSE
                  Factors(i) = EigenValues(i)
               END IF
            END DO
            WrkForces(:) =  WrkForces(:) / Factors(:)
         END IF

         ! In case of TS search, change the sign of the step in the direction of the eigenvector with lowest eigenvalue 
         IF ( TSCheck ) THEN
            i = MINLOC( EigenValues, 1 )
            WrkForces(i) = - WrkForces(i) 
         END IF

         ! Tranform displacements back to original coordinates
         WrkForces = TheOneWithMatrixVectorProduct( EigenVectors, WrkForces )

         ! Compute norm of the displacement
         IF (.NOT. SteepestDescent ) THEN
            WrkForces(:) = 0.0001 * WrkForces(:)
         ELSE
            WrkForces(:) = WrkForces(:)
         END IF
         DisplNorm = SQRT(TheOneWithVectorDotVector(WrkForces, WrkForces))

         ! Move to new coordinates
         IF ( PRESENT(Mask) ) THEN
            n = 0
            DO i = 1, size(CurrentX)
               IF ( Mask(i) ) THEN
                  n = n + 1 
                  CurrentX(i) = CurrentX(i) + WrkForces(n)
               END IF
            END DO
         ELSE
            CurrentX(:) = CurrentX(:) + WrkForces(:)
         END IF

         ! Print info to screen
!          IF ( MOD(NIter-1,NMaxIter/20) == 0 ) WRITE(*,"(I12,E20.6,E20.6)") NIter, DisplNorm, GradNorm
         IF ( MOD(NIter-1,5000) == 0 ) WRITE(*,"(I12,E20.6,E20.6)") NIter, DisplNorm, GradNorm

         ! Check convergence criteria
         IF ( GradNorm < GradThresh .AND. DisplNorm < DisplThresh ) THEN
            EXIT Iterations
         END IF

      END DO Iterations

      WRITE(*,"(I12,E20.6,E20.6)") NIter, DisplNorm, GradNorm
      ! Check max number of iterations
      CALL WARN( NIter == NMaxIter+1, " NewtonLocator: max number of iterations reached " )

      ! Store final point
      StationaryPoint = CurrentX
 
      ! Check the number of imaginary frequencies

      IF ( PRESENT(Mask) ) THEN
         ! Compute constrained Hessian at current position
         Hessian = HessianOfThePotential( CurrentX )
         WrkHessian = ConstrainedHessian( Hessian, Mask )
      ELSE
         ! compute Hessian at current position
         WrkHessian = HessianOfThePotential( CurrentX )
      END IF

      ! Diagonalize Hessian to transform coords to normal modes
      CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )
      WRITE(*,"(/,A,I3,A,/)") " Final stationary point has ",COUNT( EigenValues < 0.0 )," imaginary frequency/ies "

      DEALLOCATE( EigenValues, Factors, EigenVectors, WrkHessian, WrkForces )

   END FUNCTION NewtonLocator

!*************************************************************************************************

   FUNCTION SteepLocator( StartX, NMaxIter, GradThresh ) RESULT( StationaryPoint )
      IMPLICIT NONE
      REAL, DIMENSION(NDim), INTENT(IN) :: StartX
      INTEGER, INTENT(IN)               :: NMaxIter
      REAL, INTENT(IN)                  :: GradThresh
      REAL, DIMENSION(NDim)             :: StationaryPoint

      REAL, DIMENSION(NDim) :: CurrentX
      REAL, DIMENSION(NDim,NDim) :: EigenVectors, Hessian
      REAL, DIMENSION(NDim) :: Forces, EigenValues
      REAL :: V, GradNorm
      INTEGER :: NIter

      ! Print info to screen
      WRITE(*,"(/,A,/)") " Stationary point locator with steepest descent method "
      WRITE(*,"(A,/)") " Looking for a minimum of the potential"

      ! Start at initial position
      CurrentX = StartX
 
      WRITE(*,"(A12,A20)") "N Iteration",  "Gradient Norm"
      WRITE(*,*)              "---------------------------------------------"

      Iterations: DO NIter = 1, NMaxIter

         ! Compute Hessian and forces at current position
         V = HStickPotential( CurrentX, Forces )

         ! Compute norm of the gradient
         GradNorm  = SQRT(TheOneWithVectorDotVector(Forces, Forces))

         ! Move to new coordinates
         CurrentX(:) = CurrentX(:) + Forces(:)

         ! Print info to screen
         IF ( MOD(NIter-1,NMaxIter/20) == 0 ) WRITE(*,"(I12,E20.6,E20.6)") NIter, GradNorm

         ! Check convergence criteria
         IF ( GradNorm < GradThresh ) THEN
            EXIT Iterations
         END IF

      END DO Iterations

      WRITE(*,"(I12,E20.6,E20.6)") NIter, GradNorm

      ! Check max number of iterations
      CALL WARN( NIter == NMaxIter, " NewtonLocator: max number of iterations reached " )

      ! Store final point
      StationaryPoint = CurrentX

      ! Check the number of imaginary frequencies
      Hessian = HessianOfThePotential( StationaryPoint )
      CALL TheOneWithDiagonalization( Hessian, EigenVectors, EigenValues )

   END FUNCTION SteepLocator

!*************************************************************************************************

   LOGICAL FUNCTION DoMEPStep( X, StepLength, FirstStep, Mask )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT)  :: X
      REAL, INTENT(IN)                   :: StepLength
      LOGICAL, INTENT(IN)                :: FirstStep
      LOGICAL, DIMENSION(SIZE(X)), INTENT(IN), OPTIONAL :: Mask

      REAL, DIMENSION(:), ALLOCATABLE :: EigenValues
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenVectors, WrkHessian

      REAL, DIMENSION(SIZE(X)) :: Forces
      REAL, DIMENSION(SIZE(X),SIZE(X)) :: Hessian
      REAL :: V, FMax
      INTEGER :: i, n, NOpt, NegativeEig

      ! Number of non constrained variables
      IF ( PRESENT(Mask) ) THEN
         NOpt = COUNT(Mask)
      ELSE
         NOpt = NDim
      END IF

      IF ( FirstStep ) THEN

         ! Allocate memory
         ALLOCATE( EigenValues(NOpt), EigenVectors(NOpt, NOpt), WrkHessian(NOpt, NOpt) )

         ! Compute Hessian (constrained or not)
         IF ( PRESENT(Mask) ) THEN
            ! Compute constrained Hessian at current position
            Hessian = HessianOfThePotential( X )
            WrkHessian = ConstrainedHessian( Hessian, Mask )
         ELSE
            ! compute Hessian at current position
            WrkHessian = HessianOfThePotential( X )
         END IF

         ! Diagonalize Hessian to transform coords to normal modes
         CALL TheOneWithDiagonalization( WrkHessian, EigenVectors, EigenValues )

         ! Check that only one mode has a negative eigenvalue
         CALL ERROR( COUNT( EigenValues < 0.0 ) /= 1, " MEPStep: this is not a 1st order saddle point!" )

         ! Identify position of the negative eigenvalue
         NegativeEig = MINLOC( EigenValues, 1 )

         ! Move current position in the direction of the negative eigenvalue 
         IF ( PRESENT(Mask) ) THEN
            n = 0
            DO i = 1, size(X)
               IF ( Mask(i) ) THEN
                  n = n + 1 
                  X(i) = X(i) + StepLength * EigenVectors(n,NegativeEig) 
               END IF
            END DO
         ELSE
            X(:) = X(:) + StepLength * EigenVectors(:,NegativeEig)
         END IF

         ! Deallocate memory
         DEALLOCATE( EigenValues, EigenVectors, WrkHessian )

         ! Return TRUE value
         DoMEPStep = .TRUE.

      ELSE IF ( .NOT. FirstStep ) THEN

         ! compute forces at current position and the maximum component
         V = HStickPotential( X, Forces )
         IF ( PRESENT(Mask) ) THEN
            FMax = MAXVAL( ABS(Forces), Mask )
         ELSE 
            FMax = MAXVAL( ABS(Forces) )
         END IF

         ! Check if forces would result in a negligible step
         IF ( FMax < 1E-5 )  THEN
            ! Do not take any step
            DoMEPStep = .FALSE.
         ELSE
            ! move uncostrained coordinates in the direction of the forces
            IF ( PRESENT(Mask) ) THEN
               DO i = 1, size(X)
                  IF ( Mask(i) ) X(i) = X(i) + ABS(StepLength) * Forces(i)
               END DO
            ELSE
               X(:) = X(:) + ABS(StepLength) * Forces(:)
            END IF
            ! Return TRUE value
            DoMEPStep = .TRUE.
         END IF
      END IF

   END FUNCTION DoMEPStep

!*************************************************************************************************

   FUNCTION ConstrainedHessian( Hessian, Mask )
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(IN) :: Hessian
      LOGICAL, DIMENSION(SIZE(Hessian,1)), INTENT(IN) :: Mask
      REAL, DIMENSION(COUNT(Mask),COUNT(Mask)) :: ConstrainedHessian
      INTEGER :: i,j, n1, n2
      
      n1 = 0
      DO i = 1, SIZE(Hessian,2)
         IF ( .NOT. Mask(i) ) CYCLE
         n1 = n1+1
         n2 = 0
         DO j = 1, SIZE(Hessian,1)
            IF ( .NOT. Mask(j) ) CYCLE
            n2 = n2+1
            ConstrainedHessian(n2,n1) = Hessian(j,i)
         END DO
      END DO

   END FUNCTION ConstrainedHessian

!*************************************************************************************************

   FUNCTION ConstrainedVector( Vector, Mask )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: Vector
      LOGICAL, DIMENSION(SIZE(Vector)), INTENT(IN) :: Mask
      REAL, DIMENSION(COUNT(Mask)) :: ConstrainedVector
      INTEGER :: i, n

      n = 0
      DO i = 1, size(Vector)
         IF ( Mask(i) ) THEN
            n = n + 1
            ConstrainedVector(n) = Vector(i)
         END IF
      END DO
      
   END FUNCTION ConstrainedVector

!*************************************************************************************************

END MODULE MinimumEnergyPath
