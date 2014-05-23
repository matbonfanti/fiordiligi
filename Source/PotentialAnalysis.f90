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

   !> Parameters to plot spectra in the frequency domain
   REAL :: MaxFreq              !< Maximum frequency of the spectrum
   REAL :: DeltaFreq            !< Delta frequency of the discretization
   REAL :: SigmaPeak            !< Width parameter of the peaks in the spectrum

   !> Parameters to plot potential cuts
   REAL :: GridSpacing          !< Spacing of the grid to plot the potential
   INTEGER :: NGridPoint        !< Nr of points around the equilibrium value


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
      INTEGER :: NormalModeCutsUnit, CartesianCutsUnit

      REAL, DIMENSION(:), ALLOCATABLE   :: OmegaK, CoeffK, IntensityK
      INTEGER :: NrOfModes

      REAL, DIMENSION(NDim)      :: XMin, XAsy, XMinWithoutH, SystemCoord, BathCoord
      REAL, DIMENSION(NDim,NDim) :: HessianSysAndSlab
      REAL, DIMENSION(4,4)       :: HessianSystem

      REAL, DIMENSION(:), ALLOCATABLE   :: EigenFreq
      REAL, DIMENSION(:,:), ALLOCATABLE :: EigenModes

      REAL, DIMENSION(4) :: Dummy, XSysMin
      REAL :: PotEnergy, Norm
      REAL :: ZHFreq, ZCFreq, EMin, EAsy
      REAL :: Frequency, Spectrum, D0

      INTEGER :: iFreq, iEigen
      INTEGER :: i, j

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

      ! Write info on geometry and energy as output
      WRITE(*,501) ZHFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits),              &
                   ZCFreq*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)

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
         CALL TheOneWithDiagonalization( HessianSysAndSlab, EigenModes, EigenFreq )

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
         DO iEigen = 1, size(EigenFreq)
            Norm = 0.0
            DO i = 1, 3 
               EigenModes(i,iEigen) = EigenModes(i,iEigen) * SQRT(MassH)
               Norm = Norm + EigenModes(i,iEigen)**2
            END DO
            DO i = 4, size(EigenFreq)
               EigenModes(i,iEigen) = EigenModes(i,iEigen) * SQRT(MassC)
               Norm = Norm + EigenModes(i,iEigen)**2
            END DO
            EigenModes(:,iEigen) = EigenModes(:,iEigen) / SQRT(Norm)
         END DO
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
               IF ( IntensityK(i) > 1.E-12 )  &
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

      CLOSE(MinPhononSpectrumUnit)
      CLOSE(AsyPhononSpectrumUnit)
      CLOSE(CartesianCutsUnit)
      CLOSE(NormalModeCutsUnit)

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

   REAL FUNCTION LorentzianFunction( x, x0, gamma )
      IMPLICIT NONE
      REAL, INTENT(IN) :: x, x0, gamma
      LorentzianFunction = ( gamma / ( (x-x0)**2 + gamma**2 ) ) / MyConsts_PI
   END FUNCTION LorentzianFunction

!*************************************************************************************************

END MODULE PotentialAnalysis
