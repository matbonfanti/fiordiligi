!***************************************************************************************
!*                              PROGRAM JK6_v3
!***************************************************************************************
!>  \mainpage      Program JK6 - version 3
!>
!>  Classical simulations of H + graphene system + dissipative bath       \n
!>  * Model: 3D for H + 1D Z coordinate for carbon atom                   \n
!>  * Bath:  121 C atoms slab, normal bath, chain bath, langevin dynamics \n
!>  * Propagation: Velocity-Verlet in the microcanonical ensamble         \n
!>                 Beeman's algorithm in the canonical ensamble
!>
!>  \author        Matteo Bonfanti (Bret Jackson and Jay Kerwin for the PES sub)
!>  \version       3.0
!>  \date          7 October 2013
!>
!>  \todo          Fix temperature boundaries in the istogram of temperature
!>                 sampling (in Langevin HO test)
!>  \todo          Fix normalization of DFT and inverse DFT 
!>  \todo          Introduce more C atoms in the slab when printing the traj
!>  \todo          PERFORMANCES OF FT !!!! \n
!>                 at least, since it's a matrix vector product with fixed
!>                 matrix, on-the-fly computation of exp coeff can be avoided
!>  \todo          deallocate arrays allocated in the main program!  \n
!>
!***************************************************************************************
PROGRAM JK6_v3
   USE MyConsts
   USE ErrorTrap
   USE SharedData
   USE InputField
   USE UnitConversion
   USE VibrationalRelax
   USE PolymerVibrationalRelax
   USE PolymerEquilibriumOscillator
   USE ThermalEquilibrium
   USE Harmonic1DModel
   USE IndependentOscillatorsModel
   USE PotentialModule

   IMPLICIT NONE

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.
   
   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName

   ! Derived type to handle input data
   TYPE(InputFile) :: InputData
   ! Units of input data, defined from the input file
   INTEGER     :: InputLength, InputEnergy, InputMass, InputTime, InputTemp, InputFreq


   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                                JK6_v3            ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                    Author: Matteo Bonfanti  ')"
   PRINT "(       '         ( Potential originally implemented by B.Jackson and J.Kerwin )',/)"


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
      PRINT*, ' % JK6_v3 "InputFileName" '
      STOP
   ENDIF

   !*************************************************************
   !         INPUT SECTION 
   !*************************************************************

   ! Open and read from input file the input parameters of the calculation
   CALL OpenFile( InputData, InputFileName )

   ! read input units ( or set them to default value )

   !      ***** DEFAULT VALUES ******
   !      distance    - Angstrom
   !      energy      - electronVolt
   !      mass        - AMU
   !      time        - femtosecond
   !      temperature - Kelvin
   !      frequency   - cm-1

   CALL SetFieldFromInput( InputData, "InputLength", InputLength,  1 )
   CALL SetFieldFromInput( InputData, "InputEnergy", InputEnergy,  3 )
   CALL SetFieldFromInput( InputData, "InputMass",   InputMass,    8 )
   CALL SetFieldFromInput( InputData, "InputTime",   InputTime,   13 )
   CALL SetFieldFromInput( InputData, "InputTemp",   InputTemp,   16 )
   CALL SetFieldFromInput( InputData, "InputFreq",   InputFreq,   18 )
   CALL Initialize_UnitConversion( InputUnits, InputLength, InputEnergy, InputMass, 11, InputTime, InputTemp, InputFreq )

   ! Define the kind of simulation to do
   CALL SetFieldFromInput( InputData, "RunType", RunType, 1 )
   CALL CheckRunType( RunType )
   ! Set the representation of the bath
   CALL SetFieldFromInput( InputData, "BathType", BathType )
   CALL CheckBathType( BathType )
   ! Set the level of output
   CALL SetFieldFromInput( InputData, "PrintType", PrintType )
   CALL CheckPrintType( PrintType )
   ! Define whether the calculation is collinear or not
   CALL SetFieldFromInput( InputData, "Collinear", Collinear, .TRUE. )

   ! Hydrogen and carbon masses
   CALL SetFieldFromInput( InputData, "MassH", MassH )
   MassH = MassH * MassConversion(InputUnits, InternalUnits)
   CALL SetFieldFromInput( InputData, "MassC", MassC )
   MassC = MassC * MassConversion(InputUnits, InternalUnits)

   ! SET THE BATH-REPRESENTATION DEPENDENT VARIABLES 

   IF ( BathType == SLAB_POTENTIAL ) THEN
      ! Langevin relaxation at the border of the slab
      CALL SetFieldFromInput( InputData, "RelaxAtBorders",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )
      ! Nr of carbon atoms in the slab
      CALL SetFieldFromInput( InputData, "NCarbon",  NCarbon, 121 ) 

   ELSE IF ( BathType == NORMAL_BATH ) THEN
      ! No langevin oscillators in the normal bath
      DynamicsGamma = 0.0
      ! Mass of the bath oscillator
      CALL SetFieldFromInput( InputData, "MassBath", MassBath )
      MassBath = MassBath * MassConversion(InputUnits, InternalUnits)
      ! Nr of bath degrees of freedom
      CALL SetFieldFromInput( InputData, "NBath",  NBath ) 
      ! Read ohmic spectral density, when zero read spectral density from file
      CALL SetFieldFromInput( InputData, "OhmicRelaxT", OhmicGamma, 0.0 )
      IF ( OhmicGamma /= 0 ) THEN
         OhmicGamma = 1. / ( OhmicGamma * TimeConversion(InputUnits, InternalUnits) )
         ! For an ohmic SD, cutoff freq is compulsory
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq )
      ELSE IF ( OhmicGamma == 0.0 ) THEN
         ! Read file with spectral density
         CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
         ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, 0.0 )
      END IF
      BathCutOffFreq = BathCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Quasi-classical correction of the initial conditions of the bath (ZPE), relevant only for 0K
      CALL SetFieldFromInput( InputData, "ZPECorrection", ZPECorrection, .FALSE. )

   ELSE IF ( BathType == CHAIN_BATH ) THEN
      ! Langevin relaxation at the end of the chain
      CALL SetFieldFromInput( InputData, "RelaxAtChainEnd",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )
      ! Mass of the bath oscillator
      CALL SetFieldFromInput( InputData, "MassBath", MassBath )
      MassBath = MassBath * MassConversion(InputUnits, InternalUnits)
      ! Nr of bath degrees of freedom
      CALL SetFieldFromInput( InputData, "NBath",  NBath )
      ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
      CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, 0.0 )
      BathCutOffFreq = BathCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Read file with normal modes freq and couplings
      CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
      ! Quasi-classical correction of the initial conditions of the bath (ZPE), relevant only for 0K
      CALL SetFieldFromInput( InputData, "ZPECorrection", ZPECorrection, .FALSE. )

   ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
      ! Langevin relaxation at the end of the chain
      CALL SetFieldFromInput( InputData, "RelaxAtChainEnd",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )
      ! Mass of the bath oscillator
      CALL SetFieldFromInput( InputData, "MassBath", MassBath )
      MassBath = MassBath * MassConversion(InputUnits, InternalUnits)
      ! Nr of bath degrees of freedom
      CALL SetFieldFromInput( InputData, "NBath",  NBath )
      ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
      CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, 0.0 )
      BathCutOffFreq = BathCutOffFreq * FreqConversion(InputUnits, InternalUnits) 
      ! Read file with normal modes freq and couplings
      CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
      CALL SetFieldFromInput( InputData, "SpectralDensityFile2", SpectralDensityFile2 )

   ELSE IF ( BathType == LANGEVIN_DYN ) THEN
      ! Langevin relaxation of the system (at the carbon atom)
      CALL SetFieldFromInput( InputData, "RelaxAtCarbon",  DynamicsGamma, 0.0 ) 
      IF ( DynamicsGamma /= 0 ) &
         DynamicsGamma = 1. / ( DynamicsGamma * TimeConversion(InputUnits, InternalUnits) )

   END IF


   !*************************************************************
   !       PRINT OF THE INPUT DATA TO STD OUT
   !*************************************************************

   ! Write info about the kind of calculation
   SELECT CASE( RunType )
      CASE( EQUILIBRIUM )
         WRITE(*,"(/,A)") " * Atom-surface equilibrium simulation "
      CASE( HARMONICMODEL ) 
         WRITE(*,"(/,A)") " * Test harmonic brownian dynamics "
      CASE( RELAXATION )
         WRITE(*,"(/,A)") " * Atom-surface vibrational relaxation simulation "
      CASE( SCATTERING )
         WRITE(*,"(/,A)") " * Atom-surface sticking simulation "
      CASE( RPMD_RELAXATION )
         WRITE(*,"(/,A)") " * Atom-surface vibrational relaxation simulation with Ring Polymer MD"
      CASE( RPMD_EQUILIBRIUM )
         WRITE(*,"(/,A)") " * Equilibrium simulation with Ring Polymer MD"
      CASE( POTENTIALPRINT )
         WRITE(*,"(/,A)") " * Analysis of the potential energy surfaces "
   END SELECT

   IF ( (RunType /= HARMONICMODEL) .AND. (RunType /= RPMD_EQUILIBRIUM) ) THEN
      IF ( Collinear )  WRITE(*,"(/,A)") " * The atom is fixed in the collinear geometry "
      WRITE(*,898) MassH*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                   MassC*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits)
   END IF

   ! Write info about the bath representation
   SELECT CASE( BathType )
      CASE( SLAB_POTENTIAL )
         WRITE(*,904) NCarbon, 1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits)
      CASE( NORMAL_BATH ) 
         IF ( OhmicGamma == 0.0 ) THEN
            WRITE(*,900) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                         BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),  &
                         trim(adjustl(SpectralDensityFile))
         ELSE
            WRITE(*,910) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits), &
                         BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),  &
                         1.0/OhmicGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits)
         END IF
      CASE( CHAIN_BATH )
         WRITE(*,901) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits),   &
                      BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),    &
                      1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits), &
                      trim(adjustl(SpectralDensityFile))
      CASE( DOUBLE_CHAIN )
         WRITE(*,903) NBath, MassBath*MassConversion(InternalUnits, InputUnits), MassUnit(InputUnits),   &
                      BathCutOffFreq*FreqConversion(InternalUnits, InputUnits), FreqUnit(InputUnits),    &
                      1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits), &
                      trim(adjustl(SpectralDensityFile)), trim(adjustl(SpectralDensityFile2))
      CASE( LANGEVIN_DYN )
         IF (DynamicsGamma /= 0. ) THEN
            WRITE(*,902) 1.0/DynamicsGamma*TimeConversion(InternalUnits, InputUnits), TimeUnit(InputUnits)
         ELSE
            WRITE(*,802) 
         END IF
   END SELECT

   ! Write info about the kind of output
   SELECT CASE( PrintType )
      CASE( MINIMAL )
         WRITE(*,"(A,/)") " * Minimal output will be written "
      CASE( FULL )
         WRITE(*,"(A,/)") " * All the averages will be written to output files "
      CASE( DEBUG )
         WRITE(*,"(A,/)") " * Detailed information on each trajectory will be printed "
   END SELECT

   898 FORMAT(" * Mass of the H atom:                          ",F10.4,1X,A,/,&
              " * Mass of the C atom:                          ",F10.4,1X,A,/ )

   904 FORMAT(" * Bath is a slab of C atoms (force field potential) ", /,&
              " * Nr of Carbon atoms:                          ",I10,  /,&
              " * Langevin relax time at the edges:            ",F10.4,1X,A,/)

   900 FORMAT(" * Bath is a set of independent HO coupled to the system ",/,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * File with the spectral density:  "            ,A22,/ )

   910 FORMAT(" * Bath is a set of independent HO coupled to the system, ohmic SD ",/,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Relaxation time of the ohmic SD              ",F10.4,1X,A,/ )

   901 FORMAT(" * Bath is is a linear chain of harmonic oscillators ", /,&
              " * Nr of bath oscillators:                      ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff frequency of the bath:                ",F10.1,1X,A,/,& 
              " * Langevin relax time at the end of the chain: ",F10.4,1X,A,/,&
              " * File with the spectral density:  "            ,A22,  / )

   902 FORMAT(" * Bath is effectively represented by Langevin dynamics ", /,&
              " * Relaxation time of Langevin dynamics:        ",F10.4,1X,A,/ )
   802 FORMAT(" * Bath is effectively represented by Langevin dynamics ", /,&
              " * Infinite relaxation time                             ", / )

   903 FORMAT(" * Bath is double linear chain of harmonic oscillators ", /,&
              " * Nr of bath oscillators per chain:            ",I10,  /,& 
              " * Mass of the bath oscillator:                 ",F10.4,1X,A,/,& 
              " * Cutoff freq of the low freq chain:           ",F10.1,1X,A,/,& 
              " * Langevin relax time at the end of the chain: ",F10.4,1X,A,/,&
              " * File with the first SD:  ",                    A30,  /,&
              " * File with the second SD: ",                    A30,  / )

   !*************************************************************
   !       POTENTIAL SETUP 
   !*************************************************************

   ! Setup potential energy surface
   CALL SetupPotential( Collinear )
   
   ! If needed setup bath frequencies and coupling for oscillator bath models
   IF (  BathType == NORMAL_BATH ) THEN
         IF ( OhmicGamma == 0.0 ) THEN
            CALL SetupIndepOscillatorsModel( Bath, NBath, 0, SpectralDensityFile, MassBath, BathCutOffFreq )
         ELSE
            CALL SetupOhmicIndepOscillatorsModel( Bath, NBath, 0, OhmicGamma, MassBath, BathCutOffFreq )
         END IF
   ELSE IF (  BathType == CHAIN_BATH ) THEN
         CALL SetupIndepOscillatorsModel( Bath, NBath, 1, SpectralDensityFile, MassBath, BathCutOffFreq )
   ELSE IF ( BathType == DOUBLE_CHAIN ) THEN 
         CALL SetupIndepOscillatorsModel( DblBath(1), NBath, 1, SpectralDensityFile, MassBath, BathCutOffFreq )
         CALL SetupIndepOscillatorsModel( DblBath(2), NBath, 1, SpectralDensityFile2, MassBath, 0.0 )
   END IF
   
   IF  ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
      PRINT "(/,A,F10.6,/)"," * Bath distorsion force constant:              ", GetDistorsionForce( Bath ), " au"
   ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
      PRINT "(/,A,F10.6)"," * Bath 1 distorsion force constant:              ", GetDistorsionForce( DblBath(1) ), " au"
      PRINT "(A,F10.6,/)"," * Bath 2 distorsion force constant:              ", GetDistorsionForce( DblBath(2) ), " au"
   END IF

   !*************************************************************
   !       SPECIFIC INPUT SECTION 
   !*************************************************************

   SELECT CASE( RunType )
      CASE( EQUILIBRIUM )
         CALL ThermalEquilibrium_ReadInput( InputData )
      CASE( HARMONICMODEL ) 
         CALL Harmonic1DModel_ReadInput( InputData )
      CASE( RELAXATION )
         CALL VibrationalRelax_ReadInput( InputData )
      CASE( RPMD_RELAXATION )
         CALL PolymerVibrationalRelax_ReadInput( InputData )
      CASE( RPMD_EQUILIBRIUM )
         CALL PolymerEquilibriumOscillator_ReadInput( InputData )
      CASE( SCATTERING )
         ! ...
      CASE( POTENTIALPRINT )
         ! ...
   END SELECT

   CALL CloseFile( InputData )

   !*************************************************************
   !       INITIALIZE AND RUN CALCULATION 
   !*************************************************************

   SELECT CASE( RunType )
      CASE( EQUILIBRIUM )
         CALL ThermalEquilibrium_Initialize( )
         CALL ThermalEquilibrium_Run()
      CASE( HARMONICMODEL ) 
         CALL Harmonic1DModel_Initialize( )
         CALL Harmonic1DModel_Run()
      CASE( RELAXATION )
         CALL VibrationalRelax_Initialize( )
         CALL VibrationalRelax_Run()
      CASE( RPMD_RELAXATION )
         CALL PolymerVibrationalRelax_Initialize( )
         CALL PolymerVibrationalRelax_Run()
      CASE( RPMD_EQUILIBRIUM )
         CALL PolymerEquilibriumOscillator_Initialize( )
         CALL PolymerEquilibriumOscillator_Run( )
      CASE( SCATTERING )
         ! ...
      CASE( POTENTIALPRINT )
         ! ...
   END SELECT

   !*************************************************************
   !       DISPOSE MEMORY AND TERMINATE EXECUTION
   !*************************************************************

   SELECT CASE( RunType )
      CASE( EQUILIBRIUM )
         CALL ThermalEquilibrium_Dispose()
      CASE( HARMONICMODEL ) 
         CALL Harmonic1DModel_Dispose()
      CASE( RELAXATION )
         CALL VibrationalRelax_Dispose()
      CASE( RPMD_RELAXATION )
         CALL PolymerVibrationalRelax_Dispose()
      CASE( RPMD_EQUILIBRIUM )
         CALL PolymerEquilibriumOscillator_Dispose()
      CASE( SCATTERING )
         ! ...
      CASE( POTENTIALPRINT )
         ! ...
   END SELECT

   IF ( BathType == NORMAL_BATH .OR. BathType == CHAIN_BATH ) THEN
      CALL DisposeIndepOscillatorsModel( Bath )
   ELSE IF ( BathType == DOUBLE_CHAIN ) THEN
      CALL DisposeIndepOscillatorsModel( DblBath(1) )
      CALL DisposeIndepOscillatorsModel( DblBath(2) )
   END IF

END PROGRAM JK6_v3


!          ! if XYZ files of the trajectories are required, allocate memory to store the traj
!          IF ( PrintType >= FULL  .AND. ( RunType == SCATTERING .OR. RunType == EQUILIBRIUM ) ) THEN
!                ALLOCATE( Trajectory( 16, ntime ) )
!          END IF

! !*******************************************************************************
! ! WriteTrajectoryXYZ
! !*******************************************************************************
! !>  Given an array containing one trajectory over time, print it
! !>  to file. 
! !> 
! !> @param   CoordinatesVsTime  Array with trajectory (size: 7 x N time steps).
! !> @param   FileName           Output file name.
! !> @param   LatticeConst       Lattice constant of the graphite slab. 
! !*******************************************************************************
!    SUBROUTINE WriteTrajectoryXYZ( CoordinatesVsTime, FileName, LatticeConst )
!       IMPLICIT NONE
! 
!          REAL, DIMENSION( :,: ), INTENT(IN) :: CoordinatesVsTime
!          CHARACTER(*), INTENT(IN)           :: FileName
!          REAL, INTENT(IN)                   :: LatticeConst
! 
!          INTEGER :: NrTimeStep, UnitNr, i
!          REAL    :: A1_x, A1_y, A2_x, A2_y
! 
!          ! check and store array dimensions
!          NrTimeStep = size( CoordinatesVsTime, 2 )
!          CALL ERROR( size( CoordinatesVsTime, 1 ) /= min(16,3+nevo) , " WriteTrajectoryXYZ: Invalid array dimension "  )
! 
!          ! Open input file in a free output unit
!          UnitNr =  LookForFreeUnit()
!          OPEN( FILE = trim(adjustl(FileName)), UNIT=UnitNr )
!    
!          A1_x = LatticeConst
!          A1_y = 0.0 
!          A2_x = LatticeConst/2.0
!          A2_y = sqrt(3.)*LatticeConst/2.0
! 
!          ! write snapshot to output file
!          DO i = 1, NrTimeStep
!             WRITE( UnitNr, * ) " 14 "
!             WRITE( UnitNr, * ) " Step nr ",i," of the trajectory "
!             WRITE( UnitNr, "(A,3F15.6)" ) " H ", CoordinatesVsTime(1,i), CoordinatesVsTime(2,i)     , CoordinatesVsTime(3,i)
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , 0.0                        , CoordinatesVsTime(4,i) !C1
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0                   , LatticeConst/sqrt(3.)      , CoordinatesVsTime(5,i) !C2
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", - LatticeConst/2.0    , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(6,i) !C3
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", LatticeConst/2.0      , -LatticeConst/(2*sqrt(3.)) , CoordinatesVsTime(7,i) !C4
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x+A2_x            , -A1_y+A2_y                 , CoordinatesVsTime(8,i) !C5
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A2_x                 , +A2_y                      , CoordinatesVsTime(9,i) !C6
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x                 , -A1_y                      , CoordinatesVsTime(10,i) !C7
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x                 , +A1_y                      , CoordinatesVsTime(11,i) !C8
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A2_x                 , -A2_y                      , CoordinatesVsTime(12,i) !C9
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x-A2_x            , +A1_y-A2_y                 , CoordinatesVsTime(13,i) !C10
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", -A1_x                 , LatticeConst/sqrt(3.)-A1_y , CoordinatesVsTime(14,i) !C11
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", +A1_x                 , LatticeConst/sqrt(3.)+A1_y , CoordinatesVsTime(15,i) !C12
!             WRITE( UnitNr, "(A,3F15.6)" ) " C ", 0.0+A1_x-2*A2_x , LatticeConst/sqrt(3.)+A1_y-2*A2_y, CoordinatesVsTime(16,i) !C13
! 
! 
!          END DO
! 
!          ! close input file
!          CLOSE( UnitNr )
!    
!    END SUBROUTINE WriteTrajectoryXYZ
