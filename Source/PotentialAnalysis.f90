!***************************************************************************************
!*                              MODULE PotentialAnalysis
!***************************************************************************************
!
!>  \brief     Subroutines for an Potential Analysis of an atom-surface system
!>  \details   This module contains subroutines to analyze the potential (bath included)  \n
!>             for an atom-surface system \n
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
MODULE PotentialAnalysis
   USE MyConsts
   USE ErrorTrap
   USE MyLinearAlgebra
   USE SharedData
   USE InputField
   USE UnitConversion
   USE PotentialModule
   USE IndependentOscillatorsModel
   USE DIfferenceDerivatives

   IMPLICIT NONE

   REAL, PARAMETER :: SmallDelta = 1.E-04

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

!    ! THE FOLLOWING INPUT DATA ARE RELEVANT FOR A POTENTIAL PLOT RUN 
! 
!    INTEGER :: FreezeGraphene            !< Freeze graphene in planar geometry ( 1 = yes, 0 = no )
!    REAL    :: ZHmin, ZHmax              !< range of the ZH grid in the plots
!    REAL    :: ZCmin, ZCmax              !< range of the ZC grid in the plots
!    INTEGER :: NpointZH                  !< Nr of points in ZH
!    INTEGER :: NpointZC                  !< Nr of points in ZC

   CONTAINS

!*******************************************************************************
!> Read from input unit the variable which are specific for the
!> analysis of the potential
!>
!> @param InputData     Datatype with an already setup input unit
!*******************************************************************************
   SUBROUTINE PotentialAnalysis_ReadInput( InputData )
      IMPLICIT NONE
      TYPE(InputFile), INTENT(INOUT) :: InputData


!    ELSE IF (RunType == POTENTIALPRINT) THEN
!       CALL SetFieldFromInput( InputData, "ZHmin", ZHmin, 0.8 )
!       ZHmin = ZHmin / MyConsts_Bohr2Ang
!       CALL SetFieldFromInput( InputData, "ZHmax", ZHmax, 4.0 )
!       ZHmax = ZHmax / MyConsts_Bohr2Ang
!       CALL SetFieldFromInput( InputData, "ZCmin", ZCmin, -0.5 )
!       ZCmin = ZCmin / MyConsts_Bohr2Ang
!       CALL SetFieldFromInput( InputData, "ZCmax", ZCmax, 1.0 )
!       ZCmax = ZCmax / MyConsts_Bohr2Ang
!       CALL SetFieldFromInput( InputData, "NpointZH", NpointZH, 400 )
!       CALL SetFieldFromInput( InputData, "NpointZC", NpointZC, 100 )
!       CALL SetFieldFromInput( InputData, "FreezeGraphene", FreezeGraphene, 1 )
! 
!    ELSE IF (RunType == CHAINEIGEN) THEN 
! 
!       ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
!       BathCutOffFreq = 0.0
!       CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, BathCutOffFreq )
!       BathCutOffFreq = BathCutOffFreq * MyConsts_cmmin1toAU
!       ! Read file with spectral density / normal modes freq and couplings
!       CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
! 
!    END IF
! 

   END SUBROUTINE PotentialAnalysis_ReadInput

!-********************************************************************************************

!*******************************************************************************
!> Initialization of the data for the analysis of the Potential:
!> memory allocation, variable initialization and data type setup.
!>
!*******************************************************************************
   SUBROUTINE PotentialAnalysis_Initialize()
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

   END SUBROUTINE PotentialAnalysis_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the Potential Analysis.
!>
!*******************************************************************************
   SUBROUTINE PotentialAnalysis_Run()
      IMPLICIT NONE
 
      REAL, DIMENSION(NDim)      :: XMin, XAsy
      REAL, DIMENSION(NDim,NDim) :: HessianAtMinimum 
      REAL, DIMENSION(4,4)       :: HessianSystem

      REAL, DIMENSION(:), ALLOCATABLE   :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenModes

      REAL :: PotEnergy
      REAL :: ZHFreq, ZCFreq, EMin, EAsy

      INTEGER :: i, j

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               PES ANALYSIS "
      PRINT "(A,/)" ,    "***************************************************"

      ! ===================================================================================================
      ! (1) Coordinates and energies of the asymptotic and adsorption geometries
      ! ===================================================================================================

      PRINT "(A,/)"," Fixing the asymptotic geometry and the geometry of the adsorption minimum... "

      ! identify coordinates of the minimum of the PES
      XMin(1:2) = 0.0
      XMin(3) = HZEquilibrium
      XMin(4) = C1Puckering

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         XMin(5:NDim) = MinSlab(1:NBath)
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH .OR. BathType ==  DOUBLE_CHAIN ) THEN
         XMin(5:NDim) = 0.0
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! nothing to set
      END IF

      ! Computing the energy at this geometry
      EMin = VibrRelaxPotential( XMin, A )

      ! set asymptotic coordinates for the H atom
      XAsy(1:2) = 0.0
      XAsy(3) = 10.0

      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         XAsy(4) = 0.0
         XAsy(5:NDim) = 0.0
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH .OR. BathType ==  DOUBLE_CHAIN ) THEN
         XAsy(4) = C1Puckering
         XAsy(5:NDim) = 0.0
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         XAsy(4) = 0.0
      END IF

      ! Computing the energy at this geometry
      EAsy = VibrRelaxPotential( XAsy, A )

      PRINT "(A,/)"," Computing the hessian of the potential at the minimum... "

      ! compute hessian at the minimum of the PES
      HessianAtMinimum = HessianOfThePotential( XMin )

      ! Write info on geometry and energy as output
      WRITE(*,500) EAsy*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),              &
                   EMin*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),              &
                   (EMin-EAsy)*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),       &
                   XMin(3)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits),           &
                   XMin(4)*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits),           &
                   (XMin(3)-XMin(4))*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits) 

      ! ===================================================================================================
      ! (2) Frequencies along the ZC and ZH cut of the potential 
      ! ===================================================================================================

      ! Set starting geometry as the geometry of the minimum
      X = XMin 

      ! ZH frequency at the minimum of the potential
      X(3) = XMin(3) + SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZHFreq = - A(3) / ( 2.0 * SmallDelta )
      X(3) = XMin(3) - SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZHFreq = ZHFreq + A(3) / ( 2.0 * SmallDelta )
      ZHFreq = sqrt( ZHFreq /  MassVector(3) )
      X(3) = XMin(3)

      ! ZC frequency at the minimum of the potential
      X(4) = XMin(4) + SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZCFreq = - A(4) / ( 2.0 * SmallDelta )
      X(4) = XMin(4) - SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZCFreq = ZCFreq + A(4) / ( 2.0 * SmallDelta )
      ZCFreq = sqrt( ZCFreq /  MassVector(4) )
      X(4) = XMin(4)

      ! Write info on geometry and energy as output
      WRITE(*,501) ZHFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),              &
                   ZCFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)

      ! ===================================================================================================
      ! (3) Frequencies and normal modes of the 4D potential 
      ! ===================================================================================================

      ! Numerical mass scaled hessian of the 4D system potential
      ! (in general cannot be extracted from the full hessian for the distortion contribution in
      !  the independent oscillator model cases )
      CALL FivePointDifferenceHessian( FourDPotentialOnly, XMin(1:4), SmallDelta, HessianSystem )
      DO j = 1, 4
         DO i = 1, 4
            HessianSystem(i,j) = HessianSystem(i,j) / SQRT( MassVector(i)*MassVector(j) )
         END DO
      END DO

      ! Diagonalize the hessian
      ALLOCATE( EigenFreq(4), EigenModes(4,4) )
      CALL TheOneWithDiagonalization( HessianSystem, EigenModes, EigenFreq )

      WRITE(*,502) SQRT(EigenFreq(1))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:, 1), &
                   SQRT(EigenFreq(2))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:, 2), &
                   SQRT(EigenFreq(3))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:, 3), &
                   SQRT(EigenFreq(4))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:, 4)

      DEALLOCATE( EigenFreq, EigenModes )

      ! PHONON SPECTRUM OF THE BATH WITH CARBON ATOMS IN EQUILIBRIUM GEOMETRY FOR CH ADSORPTION
      ! AND H NEAR THE SURFACE AND FAR FROM THE SURFACE 

      500 FORMAT (/, " H-Graphite adsorption geometry and energy   ",/  &
                     " * Asymptotic energy:           ",1F15.6,1X,A, /, &
                     " * Energy at the minimum:       ",1F15.6,1X,A, /, &
                     " * Adsorption energy:           ",1F15.6,1X,A, /, &
                     " * ZH at the minimum:           ",1F15.6,1X,A, /, &
                     " * ZC at the minimum:           ",1F15.6,1X,A, /, &
                     " * C-H equilibrium distance:    ",1F15.6,1X,A     )

      501 FORMAT (/, " Frequencies along Zc and Zh cuts            ",/  &
                     " * Frequency along ZH:          ",1F15.2,1X,A, /, &
                     " * Frequency along ZC:          ",1F15.2,1X,A     ) 

      502 FORMAT (/, " 4D system potential normal modes             ",/, &
                     " 1) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass Scaled Coords / ",A," : ",4F12.6,      /, &
                     " 2) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass Scaled Coords / ",A," : ",4F12.6,      /, &
                     " 3) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass Scaled Coords / ",A," : ",4F12.6,      /, &
                     " 4) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass Scaled Coords / ",A," : ",4F12.6          )


! 
! !***************************************************************************************************
! !                       PRINT POTENTIAL CUTS
! !***************************************************************************************************
! 
!    SUBROUTINE PotentialCuts()
!       IMPLICIT NONE
! 
!       ! Allocate temporary array to store potential data and coord grids
!       ALLOCATE( PotentialArray( NpointZH * NpointZC ), ZCArray( NpointZC ), ZHArray( NpointZH ) )
!    
!       ! Define grid spacings
!       DeltaZH = (ZHmax-ZHmin)/(NpointZH-1)
!       DeltaZC = (ZCmax-ZCmin)/(NpointZC-1)
! 
!       ! Define coordinate grids
!       ZCArray = (/ ( ZCmin + DeltaZC*(iZC-1), iZC=1,NpointZC) /)
!       ZHArray = (/ ( ZHmin + DeltaZH*(iZH-1), iZH=1,NpointZH) /)
! 
!       !*************************************************************
!       ! PRINT H VIBRATIONAL CURVE FOR C FIXED IN PUCKERING GEOMETRY
!       !*************************************************************
! 
!       ! Open output file to print the H vibrational curve with fixed carbons
!       HCurveOutputUnit = LookForFreeUnit()
!       OPEN( FILE="GraphiteHBindingCurve.dat", UNIT=HCurveOutputUnit )
!       WRITE(HCurveOutputUnit, "(A,I6,A,/)") "# H-Graphite binding at C fixed puckered geometry - ", NpointZH, " points (Ang | eV)"
! 
!       ! Set collinear H, C1 puckered and other Cs in ideal geometry
!       X(:) = 0.0
!       X(4) = MinimumZC
! 
!       ! Cycle over the ZH coordinate values
!       DO iZH = 1, NpointZH
!          ! Define the H Z coordinate
!          X(3) = ZHArray(iZH)
!          ! Compute potential
!          PotEnergy = VHSticking( X, A )
!          ! Print energy to dat file
!          WRITE(HCurveOutputUnit,*) X(3)*MyConsts_Bohr2Ang, PotEnergy*MyConsts_Hartree2eV
!       END DO
!       CLOSE( HCurveOutputUnit )
! 
!       WRITE(*,"(/,A)") " * CH binding curve written to file GraphiteHBindingCurve.dat"
! 
!       !*************************************************************
!       ! PRINT 2D C-H V WITH OTHERS CARBONS IN OPTIMIZED GEOMETRY
!       !*************************************************************
! 
!       ! set optimization mask
!       IF ( FreezeGraphene == 1 ) THEN
!          OptimizationMask = (/ (.TRUE., iCoord=1,2), (.FALSE., iCoord=1,nevo+1) /)
!       ELSE 
!          OptimizationMask = (/ (.TRUE., iCoord=1,2), (.FALSE., iCoord=1,2), (.TRUE., iCoord=1,nevo-4) /)
!       ENDIF
! 
!       nPoint = 0
!       ! Cycle over the ZC coordinate values
!       DO iZC = 1, NpointZC
!          ! Cycle over the ZH coordinate values
!          DO iZH = 1, NpointZH
!             nPoint = nPoint + 1
! 
!             ! Set collinear H and other Cs in ideal geometry
!             X(:) = 0.0
!             X(4) = ZCArray(iZC)
!             X(3) = ZHArray(iZH)
! 
!             ! Compute potential at optimized geometry
!             PotentialArray(nPoint) = MinimizePotential( X,  OptimizationMask )
!          END DO
!       END DO
!       
!       ! Print the potential to vtk file
!       CALL VTK_NewRectilinearSnapshot ( PotentialCH, X=ZHArray*MyConsts_Bohr2Ang,  & 
!                                 Y=ZCArray*MyConsts_Bohr2Ang, FileName="GraphiteHSticking" )
!       CALL VTK_AddScalarField (PotentialCH, Name="CHPotential", Field=PotentialArray*MyConsts_Hartree2eV )
! 
!       WRITE(*,"(/,A)") " * PES as a func of ZC and ZH with optim puckered graphite written as VTR to file GraphiteHSticking_.vtr"
! 
!       ! Deallocate memory
!       DEALLOCATE( PotentialArray, ZHArray, ZCArray )
! 
!    END SUBROUTINE PotentialCuts

!       PRINT "(2/,A)",    "***************************************************"
!       PRINT "(A,F10.5)", "          CHAIN EIGENFREQUENCIES ANALYSIS "
!       PRINT "(A,/)" ,    "***************************************************"
! 
!       ALLOCATE( Hessian(nevo+3,nevo+3), EigenModes(nevo+3,nevo+3), EigenFreq(nevo+3) )
!       Hessian(:,:) = 0.0
! 
!       ! Open output file to print the spectrum
!       SpectrumOutUnit = LookForFreeUnit()
!       OPEN( FILE="EigenNormalModes.dat", UNIT=SpectrumOutUnit )
!       WRITE(SpectrumOutUnit, "(A,I6,A,/)") "# ", inum, " eigenfrequency spectrum of the normal modes (cm-1 | arbitrary)"
! 
!       ! Compute the part of the potential related to the system potential
!       CALL FivePointDifferenceHessian( NonCollinearOnlyV, (/ 0.0, 0.0, HZEquilibrium, C1Puckering /), dd, Hessian(1:4,1:4) )
!       ! mass scaling
!       Hessian(1:3,1:3) = Hessian(1:3,1:3) / rmh
!       Hessian(4,4) = Hessian(4,4) / rmc
!       Hessian(1:3,4) = Hessian(1:3,4) / SQRT(rmh*rmc)
!       Hessian(4,1:3) = Hessian(4,1:3) / SQRT(rmh*rmc)
! 
!       ! add the part of the hessian of the independent oscillator model
!       CALL  HessianIndepOscillatorsModel( Hessian, rmc )
! 
!       DO iFreq = 1, size(EigenFreq) 
!          WRITE(*,"(100F10.6)") Hessian(iFreq,:)
!       END DO
! 
!       ! Diagonalize hessian matrix
!       CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )
! 
!       PRINT*, " "
!       PRINT*, " frequencies (au) "
!       WRITE(*,"(5F20.12)") SQRT(EigenFreq)
!       PRINT*, " "
!       PRINT*, " frequencies (cm-1) "
!       WRITE(*,"(5F20.12)") SQRT(EigenFreq)/MyConsts_cmmin1toAU
!       PRINT*, " "
! 
!       MaxFreq = 5000.*MyConsts_cmmin1toAU
!       DeltaFreq = 1.*MyConsts_cmmin1toAU
!       SigmaPeak = 2.0*MyConsts_cmmin1toAU
! 
!       NFreq = INT( MaxFreq / DeltaFreq )
!       DO iFreq = 1, NFreq
! 
!          ! initialize frequency value and value of the signal
!          Frequency = iFreq*DeltaFreq
!          Spectrum = 0.0
! 
!          ! sum over eigenfreq centered gaussians
!          DO iEigen = 1, size(EigenFreq)
!             Spectrum = Spectrum + EXP(- (Frequency-SQRT(EigenFreq(iEigen)))**2 / 2.0 / SigmaPeak**2 )
!          END DO
! 
!          ! PRINT
!          WRITE(SpectrumOutUnit,*) Frequency/MyConsts_cmmin1toAU, Spectrum
! 
!       END DO
! 
!       DEALLOCATE( Hessian, EigenModes, EigenFreq )
!       
!       ! Close output files
!       CLOSE( SpectrumOutUnit )
! 
!    END SUBROUTINE ChainEigenfrequencies
! 


   END SUBROUTINE PotentialAnalysis_Run

!*************************************************************************************************

!*******************************************************************************
!> Free the memory which has been used for the simulation.
!>
!*******************************************************************************
   SUBROUTINE PotentialAnalysis_Dispose()
      IMPLICIT NONE

   END SUBROUTINE PotentialAnalysis_Dispose


!*************************************************************************************************

   REAL FUNCTION VibrRelaxPotential( Positions, Forces )
      REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
      REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces
      INTEGER :: NrDOF, i       

      ! Check the number degrees of freedom
      NrDOF = size( Positions )
      CALL ERROR( size(Forces) /= NrDOF, "PotentialAnalysis.VibrRelaxPotential: array dimension mismatch" )
      CALL ERROR( NrDOF /= NDim, "PotentialAnalysis.VibrRelaxPotential: wrong number of DoFs" )

      ! Initialize forces and potential
      VibrRelaxPotential = 0.0
      Forces(:)          = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 
         ! Compute potential using the potential subroutine
         VibrRelaxPotential = VHSticking( Positions, Forces )

      ELSE IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
         ! Compute potential and forces of the system
         VibrRelaxPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( Bath, Positions(4)-C1Puckering, Positions(5:), VibrRelaxPotential, &
                                                                               Forces(4), Forces(5:) ) 

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
         ! Compute potential and forces of the system
         VibrRelaxPotential = VHFourDimensional( Positions(1:4), Forces(1:4) )
         ! Add potential and forces of the bath and the coupling
         CALL BathPotentialAndForces( DblBath(1), Positions(4)-C1Puckering, Positions(5:NBath+4), VibrRelaxPotential, &
                                                                               Forces(4), Forces(5:NBath+4) ) 
         CALL BathPotentialAndForces( DblBath(2), Positions(4)-C1Puckering, Positions(NBath+5:2*NBath+4), VibrRelaxPotential, &
                                                                               Forces(4), Forces(NBath+5:2*NBath+4) ) 

      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! Compute potential using only the 4D subroutine
         VibrRelaxPotential = VHFourDimensional( Positions, Forces )

      END IF

   END FUNCTION VibrRelaxPotential

   ! ************************************************************************************************

   REAL FUNCTION FullPotentialOnly( AtCoords )
      REAL, DIMENSION(:) :: AtCoords
      REAL, DIMENSION(NDim) :: Dummy

      FullPotentialOnly = VHSticking( AtCoords, Dummy )
   END FUNCTION FullPotentialOnly

   REAL FUNCTION FourDPotentialOnly( AtCoords )
      REAL, DIMENSION(:) :: AtCoords
      REAL, DIMENSION(4) :: Dummy

      FourDPotentialOnly = VHFourDimensional( AtCoords(1:4), Dummy )
   END FUNCTION FourDPotentialOnly

   ! ************************************************************************************************

   FUNCTION HessianOfThePotential( AtPoint ) RESULT( Hessian )
      REAL, DIMENSION(NDim,NDim) :: Hessian
      REAL, DIMENSION(NDim), INTENT(IN) :: AtPoint
      REAL, DIMENSION(NBath) :: CouplingsCoeffs
      REAL :: DistortionCoeff
      INTEGER :: i, j

      Hessian(:,:) = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 

         ! Numerical mass scaled hessian of the full system+slab potential
         CALL FivePointDifferenceHessian( FullPotentialOnly, AtPoint, SmallDelta, Hessian )
         DO j = 1, NDim
            DO i = 1, NDim
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

      ELSE IF ( BathType == NORMAL_BATH ) THEN

         ! Numerical mass scaled hessian of the 4D system potential
         CALL FivePointDifferenceHessian( FourDPotentialOnly, AtPoint(1:4), SmallDelta, Hessian(1:4,1:4) )
         DO j = 1, 4
            DO i = 1, 4
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

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

         ! Numerical mass scaled hessian of the 4D system potential
         CALL FivePointDifferenceHessian( FourDPotentialOnly, AtPoint(1:4), SmallDelta, Hessian(1:4,1:4) )
         DO j = 1, 4
            DO i = 1, 4
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

         ! add the part of the hessian of the independent oscillator model
         CALL  HessianOfTheBath( Bath, Hessian(5:NDim, 5:NDim) ) 

         ! add the contribution from SYSTEM-BATH coupling and from DISTORTION CORRECTION
         CALL  CouplingAndDistortionHessian( Bath, CouplingsCoeffs, DistortionCoeff )
         Hessian(4,5) = Hessian(4,5) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(5,4) = Hessian(5,4) + CouplingsCoeffs(1) / SQRT( MassVector(4) * MassVector(5) ) 
         Hessian(4,4) = Hessian(4,4) + DistortionCoeff / MassVector(4)

      ELSE IF ( BathType == DOUBLE_CHAIN ) THEN

         ! Numerical mass scaled hessian of the 4D system potential
         CALL FivePointDifferenceHessian( FourDPotentialOnly, AtPoint(1:4), SmallDelta, Hessian(1:4,1:4) )
         DO j = 1, 4
            DO i = 1, 4
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

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

         ! Numerical hessian of the 4D system potential
         CALL FivePointDifferenceHessian( FourDPotentialOnly, AtPoint, SmallDelta, Hessian )
         DO j = 1, NDim
            DO i = 1, NDim
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

      END IF

   END FUNCTION HessianOfThePotential


!*************************************************************************************************

END MODULE PotentialAnalysis
