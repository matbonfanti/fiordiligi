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
   USE PotentialModule
!    USE ClassicalEqMotion
!    USE IndependentOscillatorsModel
!    USE RandomNumberGenerator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: PotentialAnalysis_ReadInput, PotentialAnalysis_Initialize, PotentialAnalysis_Run, PotentialAnalysis_Dispose

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

!    ELSE IF (RunType == POTENTIALPRINT) THEN 
! 
!          ! Allocation of position and acceleration arrays
!          ALLOCATE( X(nevo+3), A(nevo+3) )
! 
!    ELSE IF (RunType == CHAINEIGEN) THEN 
! 
!          ! Allocation of position and acceleration arrays
!          ALLOCATE( X(nevo+3), A(nevo+3) )
! 
!          ! Allocate and define masses
!          ALLOCATE( MassVector( nevo + 3 ) )
!          MassVector = (/ (rmh, iCoord=1,3), rmc, (1.0*MyConsts_Uma2Au, iCoord=1,nevo-1) /)
! 
!    END IF


   END SUBROUTINE PotentialAnalysis_Initialize

!-********************************************************************************************

!*******************************************************************************
!> Run the Potential Analysis.
!>
!*******************************************************************************
   SUBROUTINE PotentialAnalysis_Run()
      IMPLICIT NONE

! 
! !***************************************************************************************************
! !                       PRINT POTENTIAL CUTS
! !***************************************************************************************************
! 
!    SUBROUTINE PotentialCuts()
!       IMPLICIT NONE
! 
!       REAL, PARAMETER :: dd = 0.001
!       
!       INTEGER :: iZH, iZC, nPoint
!       INTEGER :: HCurveOutputUnit
! 
!       REAL, DIMENSION(:), ALLOCATABLE :: PotentialArray
!       REAL, DIMENSION(:), ALLOCATABLE :: ZCArray, ZHArray
! 
!       REAL :: DeltaZH, DeltaZC
!       REAL :: MinimumZH,          MinimumZC,          MinimumE
!       REAL :: ZHFreq, ZCFreq, Coupling
!       REAL :: Tmp1, Tmp2
!       
!       LOGICAL, DIMENSION(nevo+3) :: OptimizationMask
! 
!       REAL, DIMENSION(:,:), ALLOCATABLE :: Hessian, EigenModes
!       REAL, DIMENSION(:), ALLOCATABLE :: EigenFreq
! 
!       PRINT "(2/,A)",    "***************************************************"
!       PRINT "(A,F10.5)", "               PES ANALYSIS "
!       PRINT "(A,/)" ,    "***************************************************"
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
!       ! EQUILIBRIUM POSITION OF THE H-GRAPHITE POTENTIAL
!       !*************************************************************
! 
!       ! minimum of the CH 4D potential
! 
!       ! set optimization mask
!       IF ( FreezeGraphene == 1 ) THEN
!          OptimizationMask = (/ (.FALSE., iCoord=1,nevo+3) /)
!          OptimizationMask(1:4) = .TRUE.
!       ELSE 
!          OptimizationMask = (/ (.TRUE., iCoord=1,nevo+3) /)
!       ENDIF
!          
!       ! Set guess geometry: collinear H and other Cs in ideal geometry
!       X(:) = 0.0
!       X(4) = C1Puckering
!       X(3) = C1Puckering + 2.2
! 
!       ! Compute potential
!       PotEnergy = MinimizePotential( X,  OptimizationMask )
! 
!       ! Store values, with respect to the plane defined by C2,C3 and C4
!       Tmp1 = X(3)
!       Tmp2 = X(4)
!       MinimumZH = X(3) - (X(5)+X(6)+X(7))/3.0
!       MinimumZC = X(4) - (X(5)+X(6)+X(7))/3.0
!       MinimumE = PotEnergy
! 
!       ! ZH frequency at the minimum of the potential
!       X(3) = Tmp1 + dd
!       PotEnergy = VHSticking( X, A )
!       ZHFreq = - A(3) / ( 2.0 * dd )
!       X(3) = Tmp1 - dd
!       PotEnergy =  VHSticking( X, A )
!       ZHFreq = ZHFreq + A(3) / ( 2.0 * dd )
!       ZHFreq = sqrt( ZHFreq /  rmh )
! 
!       ! ZC frequency at the minimum of the potential
!       X(3) = Tmp1
!       X(4) = Tmp2 + dd
!       PotEnergy = VHSticking( X, A )
!       ZCFreq = - A(4) / ( 2.0 * dd )
!       X(4) = Tmp2 - dd
!       PotEnergy =  VHSticking( X, A )
!       ZCFreq = ZCFreq + A(4) / ( 2.0 * dd )
!       ZCFreq = sqrt( ZCFreq /  rmc )
!       
!       WRITE(*,801) MinimumE*MyConsts_Hartree2eV, MinimumZH*MyConsts_Bohr2Ang, MinimumZC*MyConsts_Bohr2Ang, &
!                    (MinimumZH-MinimumZC)*MyConsts_Bohr2Ang, ZHFreq/MyConsts_cmmin1toAU, ZCFreq/MyConsts_cmmin1toAU
!       WRITE(*,802) (MinimumE+5*MyConsts_K2AU)*MyConsts_Hartree2eV, (MinimumE+50*MyConsts_K2AU)*MyConsts_Hartree2eV,      &
!                    (MinimumE+100*MyConsts_K2AU)*MyConsts_Hartree2eV, (MinimumE+300*MyConsts_K2AU)*MyConsts_Hartree2eV,   &
!                    (MinimumE+500*MyConsts_K2AU)*MyConsts_Hartree2eV
! 
!       801 FORMAT (/,    " * Optimization with ideal graphite surface ",/  &
!                         "   (Z=0 defined by C2,C3 and C4 plane)",/,/      &
!                         "   Energy at the minimum (eV)     ",1F10.4,/     &
!                         "   Z coordinate of H atom (Ang)   ",1F10.4,/     &
!                         "   Z coordinate of C1 atom (Ang)  ",1F10.4,/     &
!                         "   H-C1 distance (Ang)            ",1F10.4,/     &
!                         "   Frequency along ZH (cm-1)      ",1F10.1,/     &
!                         "   Frequency along ZC (cm-1)      ",1F10.1           ) 
! 
!       802 FORMAT (/,    " * Average vibrational energy at given T (eV) ",/,  &
!                         "   at T = 5K   :   ",1F10.4,/, & 
!                         "   at T = 50K  :   ",1F10.4,/, &
!                         "   at T = 100K :   ",1F10.4,/, &
!                         "   at T = 300K :   ",1F10.4,/, & 
!                         "   at T = 500K :   ",1F10.4,/ ) 
! 
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
!       !*************************************************************
!       !  COMPUTE HARMONIC FREQUENCIES IN THE MINIMUM
!       !*************************************************************
! 
!       ALLOCATE( Hessian(2,2), EigenModes(2,2), EigenFreq(2) )
! 
!       CALL FivePointDifferenceHessian( CollinearOnlyV, (/ MinimumZH, MinimumZC /), 1.E-4, Hessian )
! 
!       Hessian(1,1) = Hessian(1,1) / rmh
!       Hessian(2,2) = Hessian(2,2) / rmc
!       Hessian(1,2) = Hessian(1,2) / SQRT(rmh*rmc)
!       Hessian(2,1) = Hessian(2,1) / SQRT(rmh*rmc)
! 
!       CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )
! 
!       PRINT*, " "
!       PRINT*, " frequencies "
!       PRINT*, SQRT(EigenFreq), " au"
!       PRINT*, SQRT(EigenFreq)/MyConsts_cmmin1toAU, " cm-1"
!       PRINT*, " "
!       PRINT*, " eigenvalues "
!       DO iZC = 1, 2
!          PRINT*, " eig ", iZC, " : ", EigenModes(:, iZC)
!       END DO
!       DEALLOCATE( Hessian, EigenModes, EigenFreq )
! 
!       ALLOCATE( Hessian(4,4), EigenModes(4,4), EigenFreq(4) )
! 
!       CALL FivePointDifferenceHessian( NonCollinearOnlyV, (/ 0.0, 0.0, MinimumZH, MinimumZC /), 1.E-4, Hessian )
! 
!       Hessian(1:3,1:3) = Hessian(1:3,1:3) / rmh
!       Hessian(4,4) = Hessian(4,4) / rmc
!       Hessian(1:3,4) = Hessian(1:3,4) / SQRT(rmh*rmc)
!       Hessian(4,1:3) = Hessian(4,1:3) / SQRT(rmh*rmc)
! 
!       CALL TheOneWithDiagonalization( Hessian, EigenModes, EigenFreq )
! 
!       PRINT*, " "
!       PRINT*, " frequencies "
!       WRITE(*,"(4F20.12, A)") SQRT(EigenFreq), " au"
!       WRITE(*,"(4F20.12, A)") SQRT(EigenFreq)/MyConsts_cmmin1toAU, " cm-1"
!       PRINT*, " "
!       PRINT*, " eigenvalues "
!       DO iZC = 1, 4
!          WRITE(*,"(A,I3,A,4F20.12)") " eig ", iZC, " : ", EigenModes(:, iZC)
!       END DO
!       DEALLOCATE( Hessian, EigenModes, EigenFreq )
! 
!       ! OTHER FEATURE TO IMPLEMENT IN THIS SECTION:
!       ! * calculation of the harmonic frequencies at the minimum by diagonalization
!       !   of the hessian, in the 2D collinear model and the 4D non collinear model
!       ! * by subtraction of the harmonic frequency, computation of the non harmonic coupling 
!       !   between the normal modes
!       ! * with respect to the minimum, plot of the energies corresponding to average
!       !   thermic energy
! 
! 
!    END SUBROUTINE PotentialCuts

! 
! !***************************************************************************************************
! !                       CHAIN EIGEN-FREQUENCIES
! !***************************************************************************************************
! 
!    SUBROUTINE ChainEigenfrequencies()
!       IMPLICIT NONE
! 
!       REAL, PARAMETER :: dd = 0.0001
! 
!       REAL, DIMENSION(:,:), ALLOCATABLE :: Hessian, EigenModes
!       REAL, DIMENSION(:), ALLOCATABLE :: EigenFreq
! 
!       INTEGER :: SpectrumOutUnit
! 
!       REAL :: MinimumZH, MinimumZC
!       REAL :: Frequency, Spectrum
!       REAL :: MaxFreq, DeltaFreq, SigmaPeak
!       INTEGER :: NFreq, iFreq
!       INTEGER :: iEigen
! 
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

!*******************************************************************************
!> 4D adiabatic H-Graphene potential, only potential computation 
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom and 
!>                     1 Z coordinates for the first carbon atoms (in au)
!> @return V           output potential in atomic units
!*******************************************************************************   
      REAL FUNCTION NonCollinearOnlyV( Positions ) RESULT(V) 
         IMPLICIT NONE
         REAL, DIMENSION(:) :: Positions
         REAL, DIMENSION(4) :: Dummy

         CALL ERROR( size(Positions) /= 4, " NonCollinearOnlyV: wrong nr of dofs " )
         ! use the full potential with the carbon atoms in the equilibrium geometry
         V = VHFourDimensional( Positions(:), Dummy(:) ) 
      END FUNCTION NonCollinearOnlyV

!*******************************************************************************
!> 2D adiabatic H-Graphene potential, only potential computation 
!>
!> @param Positions    Array with 1 Z coordinate for the H atom and 
!>                     1 Z coordinates for the first carbon atoms (in au)
!> @return V           output potential in atomic units
!*******************************************************************************   
      REAL FUNCTION CollinearOnlyV( Positions ) RESULT(V) 
         IMPLICIT NONE
         REAL, DIMENSION(:) :: Positions
         REAL, DIMENSION(4) :: Dummy

         CALL ERROR( size(Positions) /= 2, " CollinearOnlyV: wrong nr of dofs " )
         ! use the full potential with the carbon atoms in the equilibrium geometry
         V = VHFourDimensional( (/ 0.0, 0.0, Positions /), Dummy(:) ) 
      END FUNCTION CollinearOnlyV


END MODULE PotentialAnalysis
