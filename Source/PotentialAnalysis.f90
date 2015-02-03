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
   USE PrintTools
   USE SplineInterpolator

   IMPLICIT NONE

   REAL, PARAMETER :: SmallDelta = 1.E-04

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

   !> Nr of dimension of the system + bath 
   INTEGER :: NDim

   !> Parameters to plot spectra in the frequency domain
   REAL :: MaxFreq              !< Maximum frequency of the spectrum
   REAL :: DeltaFreq            !< Delta frequency of the discretization
   REAL :: SigmaPeak            !< Width parameter of the peaks in the spectrum

   !> Parameters to plot potential cuts
   REAL :: GridSpacing          !< Spacing of the grid to plot the potential
   INTEGER :: NGridPoint        !< Nr of points around the equilibrium value

   !> Parameters to plot 2D potential
   REAL, DIMENSION(:), ALLOCATABLE :: PotentialArray                     !< array to store 2D potential
   REAL, DIMENSION(:), ALLOCATABLE :: ZCArray, ZHArray, QArray, OptZC    !< array with coordinates grid
   REAL :: ZHmin, ZHmax, ZCmin, ZCmax, Qmin, Qmax                        !< boundaries of the grid
   INTEGER :: NpointZH, NpointZC, NpointQ                                !< nr of points of the grid
   TYPE(VTKInfo) :: PotentialCH
   TYPE(VTKInfo) :: CouplingV
   TYPE(SplineType) :: SplineData                           !< Data for spline interpolation

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


      CALL SetFieldFromInput( InputData, "PlotSpectrum_MaxFreq", MaxFreq )
      MaxFreq = MaxFreq * FreqConversion(InputUnits, InternalUnits)

      CALL SetFieldFromInput( InputData, "PlotSpectrum_DeltaFreq", DeltaFreq )
      DeltaFreq = DeltaFreq * FreqConversion(InputUnits, InternalUnits)

      CALL SetFieldFromInput( InputData, "PlotSpectrum_SigmaPeak", SigmaPeak )
      SigmaPeak = SigmaPeak * FreqConversion(InputUnits, InternalUnits)

      CALL SetFieldFromInput( InputData, "PlotSpectrum_NGridPoint", NGridPoint )

      CALL SetFieldFromInput( InputData, "PlotSpectrum_GridSpacing", GridSpacing )
      GridSpacing = GridSpacing * LengthConversion(InputUnits, InternalUnits)

      CALL SetFieldFromInput( InputData, "PlotSpectrum_ZHmin", ZHmin )
      ZHmin = ZHmin * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "PlotSpectrum_ZHmax", ZHmax )
      ZHmax = ZHmax * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "PlotSpectrum_ZCmin", ZCmin )
      ZCmin = ZCmin * LengthConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "PlotSpectrum_ZCmax", ZCmax )
      ZCmax = ZCmax * LengthConversion(InputUnits, InternalUnits)

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

      INTEGER :: AsyPhononSpectrumUnit, MinPhononSpectrumUnit, SDCouplingUnit
      INTEGER :: NormalModeCutsUnit, CartesianCutsUnit, ZHZCPathUnit

      REAL, DIMENSION(:), ALLOCATABLE   :: OmegaK, CoeffK, IntensityK
      INTEGER :: NrOfModes

      REAL, DIMENSION(NDim)      :: XMin, XAsy, XMinWithoutH, SystemCoord, BathCoord
      REAL, DIMENSION(NDim,NDim) :: HessianSysAndSlab
      REAL, DIMENSION(4,4)       :: HessianSystem
      LOGICAL, DIMENSION(NDim)   :: Mask

      REAL, DIMENSION(:), ALLOCATABLE   :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenModes

      REAL, DIMENSION(4) :: Dummy, XSysMin
      REAL :: PotEnergy, Norm
      REAL :: ZHFreq, ZCFreq, ZCAsyFreq, EMin, EAsy, ZCeq, ZH, Rho
      REAL :: Frequency, Spectrum, D0

      INTEGER :: iFreq, iEigen
      INTEGER :: i, j, nPoint

      AsyPhononSpectrumUnit = LookForFreeUnit()
      OPEN( FILE="AsymptoticPhononSpectrum.dat", UNIT=AsyPhononSpectrumUnit )
      WRITE(AsyPhononSpectrumUnit, "(A,/)") "# Slab harmonic phonon spectrum, asymptotic geometry "
      MinPhononSpectrumUnit = LookForFreeUnit()
      OPEN( FILE="InteractionPhononSpectrum.dat", UNIT=MinPhononSpectrumUnit )
      WRITE(MinPhononSpectrumUnit, "(A,/)") "# Slab harmonic phonon spectrum, adsorption minimum geometry "
      SDCouplingUnit = LookForFreeUnit()
      OPEN( FILE="SpectralDensityCoupling.dat", UNIT=SDCouplingUnit )
      WRITE(SDCouplingUnit, "(A,/)") "# Spectral density of the coupling "

      NormalModeCutsUnit = LookForFreeUnit()
      OPEN( FILE="NormalModeCuts.dat", UNIT=NormalModeCutsUnit )
      WRITE(NormalModeCutsUnit, "(A,/)") "# Potential cuts along the normal modes "
      CartesianCutsUnit = LookForFreeUnit()
      OPEN( FILE="CartesianCuts.dat", UNIT=CartesianCutsUnit )
      WRITE(CartesianCutsUnit, "(A,/)") "# Potential cuts along the cartesian coordinates "

      ZHZCPathUnit = LookForFreeUnit()
      OPEN( FILE="ZHZCPathUnit.dat", UNIT=ZHZCPathUnit )
      WRITE(ZHZCPathUnit, "(A,/)") "# ZC of the minimum potential for fixed ZH "

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
      XAsy(3) = 30.0

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

      ! Set starting geometry as the geometry of the asymptote
      X = XAsy

      ! ZC frequency with H far from the surface
      X(4) = XAsy(4) + SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZCAsyFreq = - A(4) / ( 2.0 * SmallDelta )
      X(4) = XAsy(4) - SmallDelta
      PotEnergy = VibrRelaxPotential( X, A )
      ZCAsyFreq = ZCAsyFreq + A(4) / ( 2.0 * SmallDelta )
      ZCAsyFreq = sqrt( ZCAsyFreq /  MassVector(4) )
      X(4) = XAsy(4)

      ! Write info on geometry and energy as output
      WRITE(*,501) ZHFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),              &
                   ZCFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),              &
                   ZCAsyFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)

      ! ===================================================================================================
      ! (2b) Parameters of the 4D potential with C slab at the minimum 
      ! ===================================================================================================

      IF ( BathType ==  SLAB_POTENTIAL ) THEN

         ! Set starting geometry as the geometry of the minimum
         X = XMin 
         ! set asymptotic coordinates for the H atom
         X(1:2) = 0.0
         X(3) = 30.0

         ! Set mask for optimization
         Mask = .FALSE.
         Mask(4) = .TRUE.

         ! Computing the energy at this geometry
         EMin = MinimizePotential( X, Mask ) 
         ZCeq = X(4)

         ! Computing the C1 frequency at this geometry
         X(4) = ZCeq + SmallDelta
         PotEnergy = VibrRelaxPotential( X, A )
         ZCAsyFreq = - A(4) / ( 2.0 * SmallDelta )
         X(4) = ZCeq - SmallDelta
         PotEnergy = VibrRelaxPotential( X, A )
         ZCAsyFreq = ZCAsyFreq + A(4) / ( 2.0 * SmallDelta )
         ZCAsyFreq = sqrt( ZCAsyFreq /  MassVector(4) )
         X(4) = ZCeq

         ! Write info  as output
         WRITE(*,506)   ZCeq*LengthConversion(InternalUnits,InputUnits), LengthUnit(InputUnits),      &
                        EMin*EnergyConversion(InternalUnits,InputUnits), EnergyUnit(InputUnits),      &
                        ZCAsyFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),     &
                        MassVector(4)*ZCAsyFreq**2, "au"

      END IF

      ! ===================================================================================================
      ! (3) Frequencies and normal modes of the 4D potential 
      ! ===================================================================================================

      ! Numerical mass scaled hessian of the 4D system potential
      ! (in general cannot be extracted from the full hessian for the distortion contribution in
      !  the independent oscillator model cases )
      HessianSystem = HessianOfTheSystem( XMin, MassH, MassC )

      WRITE(*,505) trim(FreqUnit(InputUnits)), HessianSystem*(FreqConversion(InternalUnits,InputUnits)**2)

      ! Diagonalize the hessian
      ALLOCATE( EigenFreq(4), EigenModes(4,4) )
      CALL TheOneWithDiagonalization( HessianSystem, EigenModes, EigenFreq )

      WRITE(*,502) SQRT(EigenFreq(1))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:,1), &
                   SQRT(EigenFreq(2))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:,2), &
                   SQRT(EigenFreq(3))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:,3), &
                   SQRT(EigenFreq(4))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits), &
                   TRIM(LengthUnit(InternalUnits)), EigenModes(:,4)

      DEALLOCATE( EigenFreq, EigenModes )

      ! ===================================================================================================
      ! (4) Phonon spectrum of the bath - ASYMPTOTIC GEOMETRY
      ! ===================================================================================================

      PRINT "(/,A,/)"," * Computing and diagonalizing the hessian of the potential at the asymptotic geometry... "

      ! compute hessian at the asymptotic geometry of the PES
      HessianSysAndSlab = HessianOfThePotential( XAsy )

      ! Compute eigenvectors and eigevalues of the hessian
      ALLOCATE( EigenFreq(NDim), EigenModes(NDim,NDim) )
      CALL TheOneWithDiagonalization( HessianSysAndSlab, EigenModes, EigenFreq )

      ! list eigenvalues
      DO iEigen = 1, NDim
         IF ( EigenFreq(iEigen) < 0.0 ) THEN
            WRITE(*,503)  iEigen, SQRT(-EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)
         ELSE
            WRITE(742,*) 0.0, SQRT(EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits)
            WRITE(742,*) 1.0, SQRT(EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits)
            WRITE(742,*) " "   
         END IF
      END DO

      DO iFreq = 1, INT( MaxFreq / DeltaFreq )

         ! initialize frequency value and value of the signal
         Frequency = REAL(iFreq)*DeltaFreq
         Spectrum = 0.0

         ! sum over eigenfreq centered gaussians
         DO iEigen = 1, size(EigenFreq)
            IF ( .NOT. EigenFreq(iEigen) < 0.0 ) & 
               Spectrum = Spectrum + LorentzianFunction( Frequency, SQRT(EigenFreq(iEigen)), SigmaPeak )
         END DO

         ! PRINT
         WRITE(AsyPhononSpectrumUnit,*) Frequency*FreqConversion(InternalUnits,InputUnits), Spectrum

      END DO

      PRINT "(/,A,/)"," Phonon spectrum written to file AsymptoticPhononSpectrum.dat"

      DEALLOCATE( EigenFreq, EigenModes )

      ! ===================================================================================================
      ! (5) Phonon spectrum of the bath - INTERACTION GEOMETRY
      ! ===================================================================================================

      PRINT "(/,A,/)"," * Computing and diagonalizing the hessian of the potential at the minimum... "

      ! compute hessian at the asymptotic geometry of the PES
      HessianSysAndSlab = HessianOfThePotential( XMin )

      ! Compute eigenvectors and eigevalues of the hessian
      ALLOCATE( EigenFreq(NDim), EigenModes(NDim,NDim) )
      CALL TheOneWithDiagonalization( HessianSysAndSlab, EigenModes, EigenFreq )

      ! list eigenvalues
      DO iEigen = 1, NDim
         IF ( EigenFreq(iEigen) < 0.0 ) THEN
            WRITE(*,503)  iEigen, SQRT(-EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)
         ELSE
            WRITE(743,*) 1.0, SQRT(EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits)
            WRITE(743,*) 2.0, SQRT(EigenFreq(iEigen))*FreqConversion(InternalUnits,InputUnits)
            WRITE(743,*) " "   
         END IF
      END DO

      DO iFreq = 1, INT( MaxFreq / DeltaFreq )

         ! initialize frequency value and value of the signal
         Frequency = REAL(iFreq)*DeltaFreq
         Spectrum = 0.0

         ! sum over eigenfreq centered gaussians
         DO iEigen = 1, size(EigenFreq)
            IF ( EigenFreq(iEigen) > 0.0 ) & 
               Spectrum = Spectrum + + LorentzianFunction( Frequency, SQRT(EigenFreq(iEigen)), SigmaPeak )
         END DO

         ! PRINT
         WRITE(MinPhononSpectrumUnit,*) Frequency*FreqConversion(InternalUnits,InputUnits), Spectrum

      END DO

      PRINT "(/,A,/)"," Phonon spectrum written to file InteractionPhononSpectrum.dat"

      DEALLOCATE( EigenFreq, EigenModes )

      ! ===================================================================================================
      ! (6) SD of the coupling
      ! ===================================================================================================

      IF ( BathType == SLAB_POTENTIAL ) THEN 

         PRINT "(/,A,/)"," * Computing the spectral density of the coupling for the zC coordinate... "

         XMinWithoutH = XMin
         XMinWithoutH(3) = 30.0 

         ! compute hessian at the asymptotic geometry of the PES
         HessianSysAndSlab = HessianOfThePotential( XMinWithoutH )

         ! Compute eigenvectors and eigevalues of the hessian
         ALLOCATE( EigenFreq(NDim), EigenModes(NDim,NDim) )
         EigenFreq = 0.0
         EigenModes = 0.0
         CALL TheOneWithDiagonalization( HessianSysAndSlab(5:NDim,5:NDim), EigenModes(5:NDim,5:NDim), EigenFreq(5:NDim) )
!          CALL TheOneWithDiagonalization( HessianSysAndSlab, EigenModes, EigenFreq )

         ! Count positive frequencies
         NrOfModes = 0
         DO iEigen = 1, NDim
            IF ( EigenFreq(iEigen) < 1.E-10 ) THEN
               WRITE(*,504)  iEigen, EigenFreq(iEigen)*(FreqConversion(InternalUnits,InputUnits))**2, FreqUnit(InputUnits)//"^2"
            ELSE
               NrOfModes = NrOfModes + 1
            END IF
         END DO

         PRINT "(/,A,I4,A,/)"," * Building spectral density with ", NrOfModes, " normal modes"

         ALLOCATE( OmegaK(NrOfModes), CoeffK(NrOfModes), IntensityK(NrOfModes) )

         SystemCoord(:) = 0.0
         SystemCoord(4) = 1.0

         ! Trasform normal mode definition to non mass scaled coordinate
!          DO iEigen = 1, size(EigenFreq)
!             Norm = 0.0
!             DO i = 1, 3 
!                EigenModes(i,iEigen) = EigenModes(i,iEigen) * SQRT(MassH)
!                Norm = Norm + EigenModes(i,iEigen)**2
!             END DO
!             DO i = 4, size(EigenFreq)
!                EigenModes(i,iEigen) = EigenModes(i,iEigen) * SQRT(MassC)
!                Norm = Norm + EigenModes(i,iEigen)**2
!             END DO
!             IF ( Norm > 1.E-12 )  EigenModes(:,iEigen) = EigenModes(:,iEigen) / SQRT(Norm)
!          END DO
         PRINT "(A)", " Data is written for peaks with I > 1.E-12 "

         PRINT "(/,A,A,A)","  # mode   |  frequency (",TRIM(FreqUnit(InputUnits)),") |   intensity (au) "
         i = 0
         D0 = 0.0
         DO iEigen = 1, size(EigenFreq)
            IF ( EigenFreq(iEigen) > 1.E-10 ) THEN
               i = i+1 
               OmegaK(i) = SQRT( EigenFreq(iEigen) )
               CoeffK(i) = DirectionalSecondDerivative( XMin, SystemCoord, EigenModes(:,iEigen) )
               IntensityK(i) = MyConsts_PI / 2.0  * CoeffK(i)**2 / OmegaK(i) / MassC
!                IF ( IntensityK(i) > 1.E-12 )  &
                           WRITE(*,"(I6,1F18.2,1E25.8)") i, OmegaK(i)*FreqConversion(InternalUnits,InputUnits), IntensityK(i)
               D0 = D0 + IntensityK(i)
            END IF
         END DO

         WRITE(*,"(/,A,1E16.8,/)") " The sum of the peak intensities is ", D0

         DO iFreq = 1, INT( MaxFreq / DeltaFreq )

            ! initialize frequency value and value of the signal
            Frequency = REAL(iFreq)*DeltaFreq
            Spectrum = 0.0

            ! sum over eigenfreq centered gaussians
            DO iEigen = 1, NrOfModes
               IF ( EigenFreq(iEigen) > 1.E-10 .AND. OmegaK(iEigen) > 30.*FreqConversion(InputUnits,InternalUnits) ) & 
                  Spectrum = Spectrum + MyConsts_PI / 2.0  * CoeffK(iEigen)**2 / OmegaK(iEigen) / MassC * &
                                 LorentzianFunction( Frequency, OmegaK(iEigen), SigmaPeak )
            END DO

            ! PRINT
            WRITE(SDCouplingUnit,*) Frequency*FreqConversion(InternalUnits,InputUnits), Spectrum

         END DO

         PRINT "(/,A,/)"," Phonon spectrum written to file SpectralDensityCoupling.dat"

         DEALLOCATE( OmegaK, CoeffK, IntensityK )
         CLOSE(SDCouplingUnit)

      ! ===================================================================================================
      ! (7) Cuts of the 4D potential
      ! ===================================================================================================

         PRINT "(/,A,/)"," * Writing 4D potential energy cuts to output files ... "

         XSysMin(:) = 0.0
         XSysMin(3) = HZEquilibrium
         XSysMin(4) = C1Puckering

         WRITE(CartesianCutsUnit, "(/,A,I6,A)") "# Cut along zH at equilibrium slab geom - ", 2*NGridPoint+1, " pts (Ang | eV)"
         X(1:4) = XSysMin(1:4)
         DO i = -NGridPoint, NGridPoint
            X(3) = HZEquilibrium + REAL(i) * GridSpacing
            WRITE(CartesianCutsUnit,*) X(3)*LengthConversion(InternalUnits,InputUnits), &
                                       VHFourDimensional( X(1:4), Dummy )*EnergyConversion(InternalUnits,InputUnits)
         END DO
         PRINT "(A)"," Written cut along zH to file NormalModeCuts.dat "

         X(1:4) = XSysMin(1:4)
         WRITE(CartesianCutsUnit, "(/,A,I6,A)") "# Cut along zC at equilibrium slab geom - ", 2*NGridPoint+1, " pts (Ang | eV)"
         X(:) = XMin(:) 
         DO i = -NGridPoint, NGridPoint
            X(4) = C1Puckering + REAL(i) * GridSpacing
            WRITE(CartesianCutsUnit,*) X(4)*LengthConversion(InternalUnits,InputUnits), &
                                       VHFourDimensional( X(1:4), Dummy )*EnergyConversion(InternalUnits,InputUnits)
         END DO
         PRINT "(A)"," Written cut along zH to file CartesianCuts.dat "


         WRITE(NormalModeCutsUnit, "(/,A,I6,A)") "# Cut along normal mode 3 - ", 2*NGridPoint+1, " pts (Ang | eV)"
         X(1:4) = XSysMin(1:4)
         DO i = -NGridPoint, NGridPoint
            X(1:4) = XSysMin(1:4) + NormalModes4D_Vecs(:,3) * REAL(i) * GridSpacing / SQRT( MassVector(1:4 ))
            WRITE(NormalModeCutsUnit,*) REAL(i) * GridSpacing, & !*LengthConversion(InternalUnits,InputUnits), &
                                        VHFourDimensional( X(1:4), Dummy ), &!*EnergyConversion(InternalUnits,InputUnits), &
            ( MinimumEnergy + 0.5 * NormalModes4D_Freq(3) * (REAL(i) * GridSpacing)**2 )!*EnergyConversion(InternalUnits,InputUnits) 
         END DO
         PRINT "(A)"," Written cut along 3rd normal mode to file NormalModeCuts.dat "

         WRITE(NormalModeCutsUnit, "(/,A,I6,A)") "# Cut along normal mode 4 - ", 2*NGridPoint+1, " pts (Ang | eV)"
         X(1:4) = XSysMin(1:4)
         DO i = -NGridPoint, NGridPoint
            X(1:4) = XSysMin(1:4) + NormalModes4D_Vecs(:,4) * REAL(i) * GridSpacing / SQRT( MassVector(1:4 ))
            WRITE(NormalModeCutsUnit,*) REAL(i) * GridSpacing, & !*LengthConversion(InternalUnits,InputUnits), &
                                        VHFourDimensional( X(1:4), Dummy ), & !*EnergyConversion(InternalUnits,InputUnits), &
            ( MinimumEnergy + 0.5 * NormalModes4D_Freq(4) * (REAL(i) * GridSpacing)**2 )!*EnergyConversion(InternalUnits,InputUnits) 
         END DO
         PRINT "(A)"," Written cut along 4th normal mode to file NormalModeCuts.dat "

      END IF

      ! ===================================================================================================
      ! (8) 2D pes in VTK format
      ! ===================================================================================================

      ! Set grid dimensions
      NpointZH = INT((ZHmax-ZHmin)/GridSpacing) + 1
      NpointZC = INT((ZCmax-ZCmin)/GridSpacing) + 1

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointZH * NpointZC ), ZCArray( NpointZC ), ZHArray( NpointZH ) )
    
      ! Define coordinate grids
      ZCArray = (/ ( ZCmin + GridSpacing*(i-1), i=1,NpointZC) /)
      ZHArray = (/ ( ZHmin + GridSpacing*(j-1), j=1,NpointZH) /)

      ! fix other coordinates 
      X(1:2) = 0.0
      IF ( BathType ==  SLAB_POTENTIAL ) THEN
         X(5:NDim) = MinSlab(1:NBath)
      ELSE IF ( BathType ==  NORMAL_BATH .OR. BathType == CHAIN_BATH .OR. BathType ==  DOUBLE_CHAIN ) THEN
         X(5:NDim) = 0.0
      ELSE IF ( BathType == LANGEVIN_DYN ) THEN
         ! nothing to set
      END IF

      ! Open VTK file
      CALL VTK_NewRectilinearSnapshot ( PotentialCH, X=ZHArray*LengthConversion(InternalUnits,InputUnits),  & 
                                Y=ZCArray*LengthConversion(InternalUnits,InputUnits), FileName="GraphiteHSticking" )

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO i = 1, NpointZC
         ! Cycle over the ZH coordinate values
         DO j = 1, NpointZH
            nPoint = nPoint + 1
            ! Set collinear H and other Cs in ideal geometry
            X(4) = ZCArray(i)
            X(3) = ZHArray(j)
            ! Compute potential at optimized geometry
            PotentialArray(nPoint) = VibrRelaxPotential( X, A )
            ! Remove asymptotic contributions of the energy
            PotentialArray(nPoint) = PotentialArray(nPoint)
         END DO
      END DO
      ! Print the potential to vtk file
      CALL VTK_AddScalarField (PotentialCH, Name="CHPotential", &
                        Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits), LetFileOpen=.TRUE. )

      nPoint = 0
      ! Cycle over the ZC coordinate values
      DO i = 1, NpointZC
         ! Cycle over the ZH coordinate values
         DO j = 1, NpointZH
            nPoint = nPoint + 1
            ! Set collinear H and other Cs in ideal geometry
            X(4) = ZCArray(i)
            X(3) = ZHArray(j)
            ! Compute potential at optimized geometry
            PotentialArray(nPoint) = VibrRelaxPotential( X, A )

            ! Remove costant asymptotic energy
            PotentialArray(nPoint) = PotentialArray(nPoint) - EMin
            ! Remove asymptotic contributions of the carbon strain
            PotentialArray(nPoint) = PotentialArray(nPoint) - 0.5 * MassVector(4) * ZCAsyFreq**2 * (X(4)-ZCeq)**2
            ! Remove asymptotic CH interaction at fixed Cz
	    X(4) = ZCeq
            PotentialArray(nPoint) = PotentialArray(nPoint) - (VibrRelaxPotential( X, A ) - EMin)

         END DO
      END DO
      ! Print the potential to vtk file
      CALL VTK_AddScalarField (PotentialCH, Name="CouplingOnly", & 
                    Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits) )

      WRITE(*,"(/,A)") " * PES as a func of ZC and ZH with reference graphite geom written as VTR to file GraphiteHSticking.vtr"

      ! Deallocate memory
      DEALLOCATE( PotentialArray, ZHArray, ZCArray )

      ! ===================================================================================================
      ! (9) ZCarbon - ZHydrogen PATH
      ! ===================================================================================================

      Rho = 0.0

      ! Set grid dimensions
      NpointZH = INT((ZHmax-ZHmin)/GridSpacing) + 1
      NpointZC = INT((0.8-0.0)/0.0001) + 1

      ! Allocate temporary array to store potential data and coord grids
      ALLOCATE( PotentialArray( NpointZC ), ZCArray( NpointZC ), OptZC( NpointZH ), ZHArray( NpointZH ) )
    
      ! Define coordinate grids
      ZCArray = (/ ( 0.0 + 0.0001*(i-1), i=1,NpointZC) /)
      ZHArray = (/ ( ZHmin + GridSpacing*(j-1), j=1,NpointZH) /)

      DO i = 1, NpointZH

         DO j = 1, NpointZC
            ! Set current geometry
            X = Xmin
            X(1) = Rho
            X(2) = 0.0
            X(4) = ZCArray(j) 
            X(3) = SQRT(ZHArray(i)**2 - Rho**2) + ZCArray(j)
            ! Compute potential at current geometry
            PotentialArray(j) = VibrRelaxPotential( X, A )
         END DO

         j = MINLOC( PotentialArray, 1 )
         OptZC(i) = ZCArray(j)

          WRITE(ZHZCPathUnit, *) SQRT((ZHArray(i)+OptZC(i))**2 + Rho**2)*LengthConversion(InternalUnits,InputUnits), &
                                OptZC(i)*LengthConversion(InternalUnits,InputUnits), &
                                PotentialArray(j)*EnergyConversion(InternalUnits,InputUnits)

      END DO

!       ! At each fixed sqrt(rho^2 + zH^2), find minimum of the potential along ZC
! 
!       ! Set ZH grid dimensions
!       NpointZH = INT((ZHmax-MAX(ZHmin,1.5))/(0.1*GridSpacing)) + 1
!       ! Allocate temporary array to store the ZH grid and the ZC of the minimum of the potential
!       ALLOCATE( PotentialArray( NpointZC ), ZCArray( NpointZH ), ZHArray( NpointZH ), OptZC( NpointZH ) )
!       ! Define coordinate grids
!       ZHArray = (/ ( MAX(ZHmin,1.5) + (0.1*GridSpacing)*(j-1), j=1,NpointZH) /)
! 
!       DO i = 1, NpointZH
! 
!          ! Set starting geometry as the geometry of the minimum
!          X = XMin 
!          ! set coordinates for the H atom
!          X(1:2) = 0.0
!          X(3) = ZHArray(i)
!          ! Set mask for optimization
!          Mask = .FALSE.
!          Mask(4) = .TRUE.
! 
!          ! Computing the energy at this geometry
!          EMin = MinimizePotential( X, Mask ) 
!          OptZC(i) = X(4)
! 
!          WRITE(ZHZCPathUnit, *) ZHArray(i)*LengthConversion(InternalUnits,InputUnits), &
!                                 OptZC(i)*LengthConversion(InternalUnits,InputUnits), &
!                                 Emin*EnergyConversion(InternalUnits,InputUnits)
! 
!       END DO

      ! Set spline interpolant through the points
      CALL SetupSpline( SplineData, ZHArray(:)+OptZC(:), OptZC )
!       CALL SetupSpline( SplineData, ZHArray(:), OptZC )

      ! Print spline interpolant
      WRITE(ZHZCPathUnit, "(/,A)") "# Spline interpolant "
      DO i = -NpointZH, NpointZH*20
         ZH =  ZHmin + REAL(i)*0.05*GridSpacing
         
         WRITE(ZHZCPathUnit, *) SQRT(ZH**2+Rho**2)*LengthConversion(InternalUnits,InputUnits), &
                                GetSpline( SplineData, ZH)*LengthConversion(InternalUnits,InputUnits)
      END DO

      ! Deallocate memory
      DEALLOCATE( ZHArray, ZCArray, PotentialArray )

      ! ===================================================================================================
      ! (10) 2D coupling potential in VTK format
      ! ===================================================================================================

      IF ( BathType ==  SLAB_POTENTIAL ) THEN

         CALL VTK_NewCollection ( CouplingV, 2, "CouplingV" )

         Qmin = -2.0
         Qmax = 2.0

         ! Set grid dimensions
         NpointZH = INT((ZHmax-ZHmin)/GridSpacing) + 1
         NpointZC = INT((ZCmax-ZCmin)/GridSpacing) + 1
         NpointQ  = INT((Qmax-Qmin)  /GridSpacing) + 1

         ! Allocate temporary array to store coord grids
         ALLOCATE( ZCArray( NpointZC ), ZHArray( NpointZH ), QArray( NpointQ ) )
         ! Define coordinate grids
         ZCArray = (/ ( ZCmin + GridSpacing*(i-1), i=1,NpointZC) /)
         ZHArray = (/ ( ZHmin + GridSpacing*(j-1), j=1,NpointZH) /)
         QArray  = (/ ( Qmin  + GridSpacing*(j-1), j=1,NpointQ)  /)

         ! A) Coupling potential as function of Q and ZH

         ! Allocate temporary array to store potential data
         ALLOCATE( PotentialArray( NpointZH * NpointQ ) )

         ! fix coordinates 
         X(1:4) = (/ 0.0, 0.0, HZEquilibrium, C1Puckering /)
         X(5:NDim) = MinSlab(1:NBath)
         PotEnergy = VibrRelaxPotential( X, A )

         ! Open VTK file
         CALL VTK_NewRectilinearSnapshot ( CouplingV, X=ZHArray*LengthConversion(InternalUnits,InputUnits),  & 
                                   Y=QArray*LengthConversion(InternalUnits,InputUnits), FileName="GraphiteHCouplingV_ZH" )

         nPoint = 0
         ! Cycle over the ZC coordinate values
         DO i = 1, NpointQ
            ! Cycle over the ZH coordinate values
            DO j = 1, NpointZH
               nPoint = nPoint + 1
               ! Set collinear H and other Cs in ideal geometry
               X(5:7) = QArray(i)
               X(3) = ZHArray(j)
               ! Compute potential at current geometry
               PotentialArray(nPoint) = VibrRelaxPotential( X, A )
               ! Remove potential of the system
               X(5:NDim) = MinSlab(1:NBath)
               X(3) = ZHArray(j)
               PotentialArray(nPoint) = PotentialArray(nPoint) - VibrRelaxPotential( X, A )
               ! Remove potential of the slab
               X(1:4) = (/ 0.0, 0.0, HZEquilibrium, C1Puckering /)
               X(5:7) = QArray(i)
               PotentialArray(nPoint) = PotentialArray(nPoint) - VibrRelaxPotential( X, A )
               ! shoft the V of  the equilibrium geometry to zero
               PotentialArray(nPoint) = PotentialArray(nPoint) + PotEnergy
            END DO
         END DO

         ! Print the potential to vtk file
         CALL VTK_AddScalarField (CouplingV, Name="CHPotential", &
                           Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits) )

         DEALLOCATE( PotentialArray )

         ! A) Coupling potential as function of Q and ZH

         ! Allocate temporary array to store potential data
         ALLOCATE( PotentialArray( NpointZC * NpointQ ) )

         ! fix coordinates 
         X(1:4) = (/ 0.0, 0.0, HZEquilibrium, C1Puckering /)
         X(5:NDim) = MinSlab(1:NBath)

         ! Open VTK file
         CALL VTK_NewRectilinearSnapshot ( CouplingV, X=ZCArray*LengthConversion(InternalUnits,InputUnits),  & 
                                   Y=QArray*LengthConversion(InternalUnits,InputUnits), FileName="GraphiteHCouplingV_ZC" )

         nPoint = 0
         ! Cycle over the ZC coordinate values
         DO i = 1, NpointQ
            ! Cycle over the ZH coordinate values
            DO j = 1, NpointZC
               nPoint = nPoint + 1
               ! Set collinear H and other Cs in ideal geometry
               X(5:7) = QArray(i)
               X(4) = ZCArray(j)
               ! Compute potential at current geometry
               PotentialArray(nPoint) = VibrRelaxPotential( X, A )
               ! Remove potential of the system
               X(5:NDim) = MinSlab(1:NBath)
               X(4) = ZCArray(j)
               PotentialArray(nPoint) = PotentialArray(nPoint) - VibrRelaxPotential( X, A )
               ! Remove potential of the slab
               X(1:4) = (/ 0.0, 0.0, HZEquilibrium, C1Puckering /)
               X(5:7) = QArray(i)
               PotentialArray(nPoint) = PotentialArray(nPoint) - VibrRelaxPotential( X, A )
               ! shoft the V of  the equilibrium geometry to zero
               PotentialArray(nPoint) = PotentialArray(nPoint) + PotEnergy
            END DO
         END DO

         ! Print the potential to vtk file
         CALL VTK_AddScalarField (CouplingV, Name="CHPotential", &
                           Field=PotentialArray*EnergyConversion(InternalUnits,InputUnits) )

         DEALLOCATE( PotentialArray )
     
         WRITE(*,"(/,A)") " * Coupling V written as VTR to file GraphiteHCouplingV_ZH.vtr and GraphiteHCouplingV_ZC.vtr"

         ! Deallocate memory
         DEALLOCATE( QArray, ZHArray, ZCArray )

      ENDIF

      CLOSE(MinPhononSpectrumUnit)
      CLOSE(AsyPhononSpectrumUnit)
      CLOSE(CartesianCutsUnit)
      CLOSE(NormalModeCutsUnit)
      CLOSE(ZHZCPathUnit)

      500 FORMAT (/, " H-Graphite adsorption geometry and energy   ",/  &
                     " * Asymptotic energy:           ",1F15.6,1X,A, /, &
                     " * Energy at the minimum:       ",1F15.6,1X,A, /, &
                     " * Adsorption energy:           ",1F15.6,1X,A, /, &
                     " * ZH at the minimum:           ",1F15.6,1X,A, /, &
                     " * ZC at the minimum:           ",1F15.6,1X,A, /, &
                     " * C-H equilibrium distance:    ",1F15.6,1X,A     )

      501 FORMAT (/, " Frequencies in the minimum geometry         ",/  &
                     " * Frequency along ZH cut:       ",1F15.2,1X,A,/, &
                     " * Frequency along ZC cut:       ",1F15.2,1X,A,2/, &
                     " Frequencies in the asymptotic geometry      ",/ &
                     " * Frequency along ZC:           ",1F15.2,1X,A    ) 

      502 FORMAT (/, " 4D system potential normal modes             ",/, &
                     " 1) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
                     " 2) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
                     " 3) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /, &
                     " 4) Frequency:                   ",1F15.2,1X,A, /, &
                     "    Mass-scaled coords of the normal mode / ",A," : ",4F12.6, /)

      503 FORMAT (   " Excluding normal mode # ",I5," with imag freq = ",1F12.6,1X,A )
      504 FORMAT (   " Excluding normal mode # ",I5," with eigenvalue = ",1F12.6,1X,A )

      505 FORMAT (/, " Mass-scaled Hessian of the system potential in the minimum (",A," ^2)",/, &
                     4E20.8,/,4E20.8,/,4E20.8,/,4E20.8,/ )

      506 FORMAT (/, " H-Graphite 4D potential parameters   ",/  &
                     " * ZC equilibrium position      ",1F15.6,1X,A, /, &
                     " * Asymptotic energy:           ",1F15.6,1X,A, /, &
                     " * ZC asymptotic frequency      ",1F15.6,1X,A, /, &
                     " * ZC Force constant            ",1F15.6,1X,A     )

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

   !*************************************************************************************************

   FUNCTION HessianOfThePotential( AtPoint ) RESULT( Hessian )
      REAL, DIMENSION(NDim,NDim) :: Hessian
      REAL, DIMENSION(NDim), INTENT(IN) :: AtPoint
      REAL, DIMENSION(NDim) :: Coordinates, FirstDerivative
      REAL, DIMENSION(NBath) :: CouplingsCoeffs
      REAL :: DistortionCoeff, Potential
      INTEGER :: i, j, k

      REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
      REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 

      REAL, DIMENSION(3), PARAMETER :: ForwardDeltas = (/  0.0,   +1.0,  +2.0   /)
      REAL, DIMENSION(3), PARAMETER :: ForwardCoeffs = (/ -3./2., +2.0,  -1./2. /) 

      Hessian(:,:) = 0.0

      IF ( BathType == SLAB_POTENTIAL ) THEN 

         ! Compute the second derivatives for displacements of x and y
         ! IMPORTANT!!! since rho = 0 is a singular value of the function, 
         ! the derivative is computed slightly off the minimum, and is computed for x,y > 0
         DO i = 1, 2
            DO k = 1, size(ForwardDeltas)

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               IF ( Coordinates(i) < 0.0 ) THEN
                  Coordinates(i) = - Coordinates(i)
               END IF

               IF ( Coordinates(i) < 0.001 ) THEN
                  Coordinates(i) = Coordinates(i) + 0.001 + ForwardDeltas(k)*SmallDelta
               ELSE
                  Coordinates(i) = Coordinates(i) + ForwardDeltas(k)*SmallDelta
               END IF

               ! Compute potential and forces in the displaced coordinate
               Potential = VHSticking( Coordinates, FirstDerivative )
               FirstDerivative = - FirstDerivative

               ! Increment numerical derivative of the analytical derivative
               Hessian(i,:) = Hessian(i,:) + ForwardCoeffs(k)*FirstDerivative(:)/SmallDelta

            END DO
         END DO

         DO i = 3, NDim
            DO k = 1, size(Deltas)

               ! Define small displacement from the point where compute the derivative
               Coordinates(:) = AtPoint(:)
               Coordinates(i) = Coordinates(i) + Deltas(k)*SmallDelta

               ! Compute potential and forces in the displaced coordinate
               Potential = VHSticking( Coordinates, FirstDerivative )
               FirstDerivative = - FirstDerivative

               ! Increment numerical derivative of the analytical derivative
               Hessian(i,:) = Hessian(i,:) + Coeffs(k)*FirstDerivative(:)/SmallDelta

            END DO
         END DO

         ! Numerical mass scaled hessian of the full system+slab potential
         DO j = 1, NDim
            DO i = 1, NDim
               Hessian(i,j) = Hessian(i,j) / SQRT( MassVector(i)*MassVector(j) )
            END DO
         END DO

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

!       CALL TheOneWithMatrixPrintedLineAfterLine( Hessian )

   END FUNCTION HessianOfThePotential

!*************************************************************************************************

   REAL FUNCTION DirectionalSecondDerivative( AtPoint, Direction, NormalMode )
      IMPLICIT NONE
      REAL, DIMENSION(NDim), INTENT(IN) :: AtPoint, Direction, NormalMode

      REAL, DIMENSION(4), PARAMETER :: Deltas = (/ -2.0,    -1.0,    +1.0,    +2.0    /)
      REAL, DIMENSION(4), PARAMETER :: Coeffs = (/ +1./12., -8./12., +8./12., -1./12. /) 

      REAL, DIMENSION(NDim) :: Coordinates, FirstDerivative
      REAL :: Potential, Norm
      INTEGER :: k

      Norm = 0.0
      DO k = 1, size(Direction)
         Norm = Norm + Direction(k)**2
      END DO
      CALL ERROR( ABS(SQRT(Norm)-1.0) > 1E-12 , " DirectionalSecondDerivative: Direction is not normalized" )
      Norm = 0.0
      DO k = 1, size(NormalMode)
         Norm = Norm + NormalMode(k)**2
      END DO
      CALL ERROR( ABS(SQRT(Norm)-1.0) > 1E-12 , " DirectionalSecondDerivative: NormalMode is not normalized" )

      DirectionalSecondDerivative = 0.0

      DO k = 1, size(Deltas)
         Coordinates(:) = AtPoint(:) + SmallDelta*Deltas(k)*NormalMode(:)
         Potential = VHSticking( Coordinates, FirstDerivative )
         DirectionalSecondDerivative = DirectionalSecondDerivative - &
                          Coeffs(k)* TheOneWithVectorDotVector(Direction,FirstDerivative) /SmallDelta
      END DO

   END FUNCTION DirectionalSecondDerivative

!*************************************************************************************************

   SUBROUTINE PrintCouplingPotential
      IMPLICIT NONE

   END SUBROUTINE PrintCouplingPotential

!*************************************************************************************************

   REAL FUNCTION LorentzianFunction( x, x0, gamma )
      IMPLICIT NONE
      REAL, INTENT(IN) :: x, x0, gamma
      LorentzianFunction = ( gamma / ( (x-x0)**2 + gamma**2 ) ) / MyConsts_PI
   END FUNCTION LorentzianFunction

!*************************************************************************************************

END MODULE PotentialAnalysis
