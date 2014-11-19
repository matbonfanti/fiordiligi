PROGRAM NewWiener
   USE MyConsts
   USE ErrorTrap
   USE InputField
   USE UnitConversion
   USE FFTWrapper

   IMPLICIT NONE

   ! Variable to handle the command line
   INTEGER :: NArgs
   LOGICAL :: Help = .FALSE.

   ! Input file name, set from command line arguments
   CHARACTER(120) :: InputFileName

   ! Discrete Fourier TRansform data type
   TYPE(FFTComplexType) :: DFTData

   ! Name of the file storing the process realizations
   CHARACTER(120) :: SetFile
   ! Nr of realizations to include in the average, counter on the realizations
   INTEGER        :: NSet, iSet
   ! Nr of steps of each realization
   INTEGER        :: NStep, iStep, NCutOff
   ! Unit to read input sets of data, Unit to write output results
   INTEGER        :: SetsUnit, OutUnit
   ! Time step of the input data sets, Frequency spacing, Omega cut off for hilbert transform
   REAL           :: DeltaT, DeltaOmega, OmegaCutOff, LowOmegaCutOff
   ! Time and frequency grid
   REAL, DIMENSION(:), ALLOCATABLE :: TimeGrid, OmegaGrid
   ! Array to store dataset and average
   REAL, DIMENSION(:), ALLOCATABLE :: Hilbert, JOmega, DataSet
   COMPLEX, DIMENSION(:), ALLOCATABLE :: SetTransf, SOmega, CorrT, Cauchy
   ! Data for low frequency filter
   REAL :: OmegaMax, OmegaRatio, DampFactor
   INTEGER :: NPow
   ! Data for damping in time
   REAL :: TauDamp, T0Damp
   ! Data for removing LorentzDistrib from Spectral Density
!    REAL :: LorentzI, LorentzOmega, LorentzGamma, LorentzF
   LOGICAL :: NormalMode
   REAL    :: NormalModeCoeff
   ! Various real data
   REAL :: Time, SetAver, SetSquared, Temperature
   LOGICAL :: VelocitySet
   CHARACTER(100) :: ErrorMsg
   REAL :: D0, Omega0, Dw0, w02, dom02
   REAL :: MassHydro, MassCarbon

   PRINT "(/,     '                    ==============================')"
   PRINT "(       '                                NewWiener       ')"
   PRINT "(       '                    ==============================',/)"
   PRINT "(       '                    Author: R.Martinazzo and M.Bonfanti  '/)"


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
      PRINT*, ' % NewWiener "InputFileName" '
      STOP
   ENDIF

   !*************************************************************
   !                        READ INPUT FILE
   !*************************************************************

   CALL ReadInput( InputFileName )

   !*************************************************************
   !                    PROCESS INPUT DATA
   !*************************************************************

   ! DEFINE TIME AND FREQUENCY GRID
   ALLOCATE( TimeGrid(NStep), OmegaGrid(NStep) )

   ! Define frequency spacing
   DeltaOmega = 2.0*MyConsts_PI/(REAL(NStep)*DeltaT)
   ! Number of data in the interval 0, OmegaCutoff
   NCutOff = INT( OmegaCutOff/DeltaOmega ) 

   ! Build grids in time and frequency
   DO iStep = 1, NStep
      TimeGrid(iStep) = REAL(iStep-1)*DeltaT
      OmegaGrid(iStep) = REAL(iStep-1)*DeltaOmega
      IF (OmegaGrid(iStep) > MyConsts_PI/DeltaT )      &
                       OmegaGrid(iStep) = OmegaGrid(iStep) - 2.0 * MyConsts_PI / DeltaT
   END DO
   
   ! ALLOCATE MEMORY
   ALLOCATE( DataSet(NStep), SetTransf(NStep), SOmega(NStep), CorrT(NStep) ) 
   ALLOCATE( Hilbert(NCutOff), Cauchy(NCutOff), JOmega(NCutOff) )
   SOmega( : ) = CMPLX( 0.0 , 0.0 )

   ! Setup fourier transform
   CALL SetupFFT( DFTData, NStep )

   !*************************************************************
   !                    PRINT INFO ON INPUT DATA
   !*************************************************************
   
   write(*,*)' '
   write(*,*)'   ************** New Wiener ************* '
   write(*,*)'    '
   write(*,*)'   From Position Autocorrelation Function  '
   write(*,*)'     to Spectral Density                   '
   write(*,*)' '
   write(*,*)'   *************************************** '
   write(*,*)' '
   write(*,*) ' Data read from file ', TRIM(ADJUSTL(SetFile))
   write(*,*) ' '
   write(*,*) ' Successful read input: '
   write(*,*) ' '
   write(*,*) ' Number of points          =', NStep
   write(*,*) ' Time step                 =', DeltaT * TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)
   write(*,*) ' Temperature               =', Temperature * TemperatureConversion(InternalUnits,InputUnits), TemperUnit(InputUnits)
   write(*,*) ' Damping relaxation time   =', TauDamp * TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)
   write(*,*) ' Starting relaxation time  =', T0Damp * TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)
   write(*,*) '  '
   write(*,*) ' Number of trajectories    =', NSet
   write(*,*) ' Start of frequency filter =', OmegaMax * FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)
   write(*,*) ' Filter power              =', NPow
   write(*,*) ' Frequency cutoff          =', OmegaCutOff * FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)

   write(*,*) '  '
!    if(chain) then
!       write(*,*) ' Compute first effective mode chain parameters and density '
!       write(*,*) ' Effective mode mass  (u.a.m.)  =', mass / xmass
!       write(*,*) ' Effective mode cutoff  (cm-1)  =', wcut*etocm
!       write(*,*) '  '
!    endif 
   write(*,*) ' '
   write(*,*) '             *** Time analysis ***                '
   write(*,*) ' '
   write(*,"(x,a,f12.2,1x,a,f10.2,1x,a)") ' T final    = ', REAL(NStep)*DeltaT, TimeUnit(InternalUnits), &
                         REAL(NStep)*DeltaT*TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)
   write(*,"(x,a,e12.5,1x,a,f10.2,1x,a)") ' Del-omega  = ', DeltaOmega, FreqUnit(InternalUnits), &
                        DeltaOmega*FreqConversion(InternalUnits,InputUnits), FreqUnit(InputUnits)
   write(*,*) ' Frequency cutoff   = ', OmegaCutOff*FreqConversion(InternalUnits,InputUnits), FreqUnit(InternalUnits)
   write(*,*) '          ....  t_c = ', 2.*MyConsts_PI/OmegaCutOff * TimeConversion(InternalUnits,InputUnits), TimeUnit(InputUnits)
   write(*,*) ' '

   !*************************************************************
   !    READ DATA SET AND STORE AVERAGE FOURIER TRANSFORM
   !*************************************************************

   SetsUnit = LookForFreeUnit()
   OPEN( UNIT=SetsUnit, FILE=TRIM(ADJUSTL(SetFile)) )
   ! Skip two lines
   READ( SetsUnit, * )
   READ( SetsUnit, * )

   ! Cycle over the input sets
   DO iSet = 1, NSet

      ! Read a single set
      DO iStep = 1, NStep
         ! Read data
         READ( SetsUnit, * ) Time, DataSet(iStep)
         ! Check time
         WRITE( ErrorMsg, * ) " Incorrect time grid for set ",iSet, " at step ",iStep
         CALL ERROR( Time*TimeConversion(InputUnits, InternalUnits) /= TimeGrid(iStep), ErrorMsg )
         ! Transform to internal coordinates
         IF ( VelocitySet ) THEN
               DataSet(iStep) = DataSet(iStep) * &
                                      LengthConversion(InputUnits, InternalUnits) / TimeConversion(InputUnits, InternalUnits)
         ELSE 
               DataSet(iStep) = DataSet(iStep) * LengthConversion(InputUnits, InternalUnits)
         END IF
      END DO
      READ( SetsUnit , * )

      ! Calculate average and shift
      SetAver = SUM( DataSet ) / REAL(NStep)
      SetSquared = DOT_PRODUCT( DataSet,DataSet )/REAL(NStep)  
      IF ( .NOT. VelocitySet ) THEN
	 DataSet(:) = DataSet(:) - SetAver
      ENDIF

      ! fourier transform of the noise, x(w)
      SetTransf(:) = CMPLX( DataSet(:)*DeltaT, 0.0 )

      CALL ExecuteFFT( DFTData, SetTransf, INVERSE_FFT  )

      ! Divide by w^2 if velocity data was given
!       IF ( VelocitySet ) THEN
! 	 DO iStep = 1, NStep 
! 	    WRITE(600+iSet,"(I5,10E20.6)") iStep, REAL(SetTransf( iStep ))**2, MAX(OmegaGrid(iStep)**2,1.E-6) &
!                    ,REAL(SetTransf( iStep ))**2 / MAX(OmegaGrid(iStep)**2,1.E-6)
! 	    SetTransf( iStep ) = SetTransf( iStep ) / MAX(OmegaGrid(iStep)**2, 1.0)
! 	 END DO
!       ENDIF 

      ! Spectral density of the noise, S(w)=C(w)/2.d0/pi   
      IF ( VelocitySet ) THEN
	 DO iStep = 1, NStep              
	    IF ( OmegaGrid(iStep) > 5.0 * FreqConversion(InputUnits,InternalUnits) ) THEN
	       SOmega( iStep ) = SOmega( iStep ) +  &
               ABS(SetTransf(iStep))**2 * REAL(NStep) / ( 2.0 * MyConsts_PI * DeltaT ) / OmegaGrid(iStep)**2
	    ELSE
	       SOmega( iStep ) = 0.0 
	    END IF
	 END DO
      ELSE
	 DO iStep = 1, NStep              
	    SOmega( iStep ) = SOmega( iStep ) + ABS(SetTransf(iStep))**2 * REAL(NStep) / ( 2.0 * MyConsts_PI * DeltaT )  
	 END DO
      END IF

      write(*,"(/,A,I5,A)") ' ******* Set number ', iSet, ' *********'
      write(*,*) ' Average value (a.u.) =', SetAver
      write(*,*) ' <X(0)^2>  from data         = ', SetSquared 
!       write(*,*) '           form WK           = ', real(corr(1))
   
   END DO
   CLOSE( SetsUnit )

   ! Normalize average
   SOmega( : ) = SOmega( : ) / REAL( NSet )

   ! If Normal mode, multiply by factor
   IF (NormalMode) THEN
      SOmega( : ) = SOmega( : ) * NormalModeCoeff
   END IF

   ! Apply low frequency filter
   DO iStep = 1, NStep
      IF ( ABS( OmegaGrid(iStep) ) >= ABS(OmegaMax) ) THEN
         OmegaRatio = OmegaMax / OmegaGrid(iStep)
         DampFactor = EXP(1.0)*EXP( - 1.0/ ( 1.0-OmegaRatio**NPow ))
         SOmega( iStep ) = (1.0-DampFactor) * SOmega( iStep ) + DampFactor * SOmega( iStep ) * OmegaRatio**NPow
      ENDIF
   END DO 

   CorrT( : ) = SOmega( : ) * DeltaOmega
   CALL ExecuteFFT( DFTData, CorrT, DIRECT_FFT  )

   ! Print correlation function in time
   OutUnit = LookForFreeUnit()
   OPEN( UNIT=OutUnit, FILE='acorr.dat' )
   DO iStep = 1, NStep
      Time = REAL(iStep-1)*DeltaT
      IF ( Time > REAL(NStep)*DeltaT/2.0 )   Time = Time - REAL(NStep)*DeltaT
      WRITE( OutUnit,"(10E20.6)" ) Time*TimeConversion(InternalUnits,InputUnits), & 
               REAL( CorrT(iStep) )!,& !*LengthConversion(InternalUnits,InputUnits)**2 , &
               !EXP( -(T0Damp-abs(Time))**2 / (2.0*TauDamp**2) )
   ENDDO
   CLOSE( OutUnit )

   ! Apply gaussian damping in time
   DO iStep = 1, NStep 
      Time = REAL(iStep-1)*DeltaT
      IF ( Time > REAL(NStep)*DeltaT/2.0 )   Time = Time - REAL(NStep)*DeltaT
      IF ( ABS(Time) >= T0Damp )   CorrT(iStep) = CorrT(iStep) * EXP( -(T0Damp-abs(Time))**2 / (2.0*TauDamp**2) )
!       IF ( ABS(Time) >= T0Damp )   CorrT(iStep) = CorrT(iStep) * EXP( (T0Damp-abs(Time)) / TauDamp )
   ENDDO

   ! Print correlation function in time
   OutUnit = LookForFreeUnit()
   OPEN( UNIT=OutUnit, FILE='acorr-damped.dat' )
   DO iStep = 1, NStep
      Time = REAL(iStep-1)*DeltaT
      IF ( Time > REAL(NStep)*DeltaT/2.0 )   Time = Time - REAL(NStep)*DeltaT
      WRITE( OutUnit,901 ) Time*TimeConversion(InternalUnits,InputUnits), &
                           REAL( CorrT(iStep)) !*LengthConversion(InternalUnits,InputUnits)**2
   ENDDO
   CLOSE( OutUnit )

    w02  = Temperature / MassHydro / real(CorrT(1))                   ! this gives the square of the HO frequency

   ! Now compute the spectral density of the correlation function C(w)
   SOmega( : ) = CorrT( : ) * DeltaT * REAL(NStep)
   CALL ExecuteFFT( DFTData, SOmega, INVERSE_FFT  )

   ! and write the spectral density to file
   OutUnit = LookForFreeUnit()
   OPEN( UNIT=OutUnit, FILE='acorrw.dat' )
   DO iStep = 1, NStep/2-1
      IF ( OmegaGrid(iStep) >= MyConsts_PI / DeltaT ) EXIT
      WRITE( OutUnit,902 ) OmegaGrid(iStep)*FreqConversion(InternalUnits,InputUnits), &
                           REAL(SOmega(iStep)) !* LengthConversion(InternalUnits,InputUnits)**2
   ENDDO
   CLOSE( OutUnit )

   ! Now compute the cauchy transform in the range 0, Wcutoff 
   DO iStep = 1, NCutOff
      IF (OmegaGrid(iStep) > LowOmegaCutOff) THEN
	 Hilbert(iStep) = REAL(SOmega(iStep)) * OmegaGrid(iStep) / 2.0   ! this is the generating function, w*C(w)/2
      ELSE
	 Hilbert(iStep) =  0.0
      END IF
   END DO
   CALL HilbertTransform( NCutOff, Hilbert )
   DO iStep = 1, NCutOff
      IF (OmegaGrid(iStep) > LowOmegaCutOff) THEN
	 Cauchy(iStep) = CMPLX( Hilbert(iStep), 0.5*REAL(SOmega(iStep))*OmegaGrid(iStep) ) 
	 JOmega(iStep) = Temperature * (0.5 * REAL(SOmega(iStep)) * OmegaGrid(iStep)) / ABS(Cauchy(iStep))**2
      ELSE
	 Cauchy(iStep) = 0.0
	 JOmega(iStep) = 0.0
      ENDIF
   END DO     

!    DO iStep = 1, NCutOff
!       LorentzF = LorentzI * LorentzGamma**2 / ( (OmegaGrid(iStep)-LorentzOmega)**2 + LorentzGamma**2 )
! !       IF (LorentzF < JOmega(iStep)) JOmega(iStep) =  JOmega(iStep) - LorentzF
!       JOmega(iStep) =  JOmega(iStep) - LorentzF
!    END DO

   ! write w / cm^-1 and j(w)
   OutUnit = LookForFreeUnit()
   OPEN( UNIT=OutUnit, FILE='spectral.dat' )
   DO iStep = 1, NCutOff
      WRITE( OutUnit,902 ) OmegaGrid(iStep)*FreqConversion(InternalUnits,InputUnits), JOmega(iStep)   !, &
!          JOmega(iStep)+LorentzI * LorentzGamma**2 / ( (OmegaGrid(iStep)-LorentzOmega)**2 + LorentzGamma**2 )
   ENDDO
   CLOSE( OutUnit )

   ! Extract spectral density for the first effective mode of mass mass

   ! Compute D0
   D0 = 0.0
   Omega0 = 0.0
   Dw0 = 0.0
   DO iStep = 1, NCutOff
      D0 = D0 +  2.0 / MyConsts_PI * JOmega(iStep) * OmegaGrid(iStep) * DeltaOmega                   ! effective mode coupling
      IF (iStep > 1) Dw0 = Dw0 + 2.0/ MyConsts_PI * JOmega(iStep) / OmegaGrid(iStep) * DeltaOmega    ! HO renormalization
      Omega0 = Omega0 + 2.0 / MyConsts_PI * JOmega(iStep) * OmegaGrid(iStep)**3 * DeltaOmega        ! effective mode frequency
   END DO
   D0 = D0 * MassCarbon 
   Dw0 = Dw0 / MassHydro                                          ! dw^2_H
   Omega0 = Omega0 * MassCarbon / D0                              !  W^2_C

   ! Compute cauchy transform of JOmega
   DO iStep = 1, NCutOff
      Hilbert(iStep) = JOmega(iStep)
   END DO
   CALL HilbertTransform( NCutOff, Hilbert )
   DO iStep = 1, NCutOff
      Cauchy(iStep) = CMPLX( Hilbert(iStep), JOmega(iStep) ) 
   END DO

   ! Compute the jomega on the carbon
   DO iStep = 1, NCutOff
      JOmega(iStep) = D0 * JOmega(iStep) / ABS( Cauchy(iStep) )**2
   END DO

   ! write w / cm^-1 and j(w)
   OutUnit = LookForFreeUnit()
   OPEN( UNIT=OutUnit, FILE='spectralC.dat' )
   DO iStep = 1, NCutOff
      WRITE( OutUnit,902 ) OmegaGrid(iStep)*FreqConversion(InternalUnits,InputUnits), JOmega(iStep)
   ENDDO
   CLOSE( OutUnit )

   dom02 = 0.0
   DO iStep = 1, NCutOff
      IF (iStep > 1) dom02 = dom02 + 2.0/ MyConsts_PI * JOmega(iStep) / OmegaGrid(iStep) * DeltaOmega    ! HO renormalization
   END DO
   dom02 = dom02 / MassCarbon                                          ! dW^2_C

   write(*,*) '  '
   write(*,*) '           ***  Ensemble average ***           '
   write(*,*) ' HO frequency                          (cm-1) = ', sqrt(w02)*FreqConversion(InternalUnits,InputUnits)
   write(*,*) ' HO corrected frequency                (cm-1) = ', sqrt(w02 + Dw0)*FreqConversion(InternalUnits,InputUnits)
   write(*,*) ' Effective mode coupling constant D0^2 (a.u.) = ', D0
   write(*,*) ' Effective mode frequency              (cm-1) = ', sqrt(Omega0)*FreqConversion(InternalUnits,InputUnits)
   write(*,*) ' Effective mode renormalization        (cm-1) = ', sqrt(dom02)*FreqConversion(InternalUnits,InputUnits)
   write(*,*) ' Effective mode final frequency        (cm-1) = ', sqrt(Omega0-dom02)*FreqConversion(InternalUnits,InputUnits)
   write(*,*) '  '


901 FORMAT(f10.4,15(1x,g11.4))
902 FORMAT(f23.16,15(1x,g23.16))

   DEALLOCATE( TimeGrid, OmegaGrid, DataSet, SetTransf, SOmega, CorrT, Hilbert, Cauchy, JOmega )
   CALL DisposeFFT( DFTData )


CONTAINS

   SUBROUTINE ReadInput( InputFileName )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: InputFileName

      ! Derived type to handle input data
      TYPE(InputFile) :: InputData
      ! Units of input data, defined from the input file
      INTEGER  :: InputLength, InputEnergy, InputMass, InputTime, InputTemp, InputFreq
      ! Integer unit to read traj file
      INTEGER  :: UnitTraj
      ! Reading status
      INTEGER  :: RdStat
      LOGICAL :: FileIsPresent
      ! Dummy variables
      REAL :: Dummy1, Dummy2

      ! Open and read from input file the input parameters of the calculation
      CALL OpenFile( InputData, TRIM(ADJUSTL(InputFileName)) )

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

      ! Read file with data sets
      CALL SetFieldFromInput( InputData, "SetFile", SetFile )

      ! Check if SetFile file exists
      INQUIRE( File = TRIM(ADJUSTL(SetFile)), EXIST=FileIsPresent ) 
      CALL ERROR( .NOT. FileIsPresent, "NewWiener.ReadInput: dataset file does not exists" )

      ! Set the number of steps of the realizations by reading the first one from file SetFile
      UnitTraj = LookForFreeUnit()
      OPEN( FILE=TRIM(ADJUSTL(SetFile)), UNIT=UnitTraj )
      ! Skip two lines
      READ( UnitTraj, * )
      READ( UnitTraj, * )
      NStep = 0
      ! Read until a blank line is found and count lines which have been read
      DO 
         READ( UnitTraj, *, IOSTAT=RdStat ) Dummy1, Dummy2
         ! if second step, store delta T (assuming t(0) = 0 )
         IF ( NStep == 1 ) DeltaT = Dummy1 * TimeConversion(InputUnits, InternalUnits)    
         IF ( RdStat /= 0 ) EXIT
         IF ( NStep > 1 .AND. Dummy1 == 0.0 ) EXIT
         NStep = NStep + 1 
      END DO
      CLOSE( UnitTraj )

      ! Read number of datasets to include in the average
      CALL SetFieldFromInput( InputData, "NSet", NSet )

      ! Read data for low frequency filter: omega max and power decay
      CALL SetFieldFromInput( InputData, "OmegaMax", OmegaMax )
      OmegaMax = OmegaMax * FreqConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "NPow", NPow )

      ! Read data for gaussian damping in time
      CALL SetFieldFromInput( InputData, "T0Damp", T0Damp )
      T0Damp = T0Damp * TimeConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "TauDamp", TauDamp ) 
      TauDamp = TauDamp * TimeConversion(InputUnits, InternalUnits)

      ! Read frequency cutoff
      CALL SetFieldFromInput( InputData, "OmegaCutOff", OmegaCutOff )
      OmegaCutOff = OmegaCutOff * FreqConversion(InputUnits, InternalUnits)
      CALL SetFieldFromInput( InputData, "LowOmegaCutOff", LowOmegaCutOff, 0.0 )
      LowOmegaCutOff = LowOmegaCutOff * FreqConversion(InputUnits, InternalUnits)

      ! Read temperature
      CALL SetFieldFromInput( InputData, "Temperature", Temperature )
      Temperature = Temperature * TemperatureConversion(InputUnits, InternalUnits)

      ! Input datasets are velocities or positions
      CALL SetFieldFromInput( InputData, "VelocitySet", VelocitySet )

      ! Input masses
      CALL SetFieldFromInput( InputData, "MassCarbon", MassCarbon )
      CALL SetFieldFromInput( InputData, "MassHydro", MassHydro )
      MassCarbon = MassCarbon * MassConversion(InputUnits, InternalUnits)
      MassHydro = MassHydro * MassConversion(InputUnits, InternalUnits)

!       ! Removal of lorentz distribution
!       CALL SetFieldFromInput( InputData, "LorentzGamma", LorentzGamma )
!       CALL SetFieldFromInput( InputData, "LorentzI",     LorentzI )
!       CALL SetFieldFromInput( InputData, "LorentzOmega", LorentzOmega )
!       LorentzOmega = LorentzOmega * FreqConversion(InputUnits, InternalUnits)
!       LorentzGamma = LorentzGamma * FreqConversion(InputUnits, InternalUnits)

      ! Normal modes
      CALL  SetFieldFromInput( InputData, "NormalMode", NormalMode, .FALSE. )
      IF ( NormalMode ) THEN
	 CALL  SetFieldFromInput( InputData, "NormalModeCoeff", NormalModeCoeff )
      ENDIF

      ! Close input file
      CALL CloseFile( InputData )

   END SUBROUTINE ReadInput

   SUBROUTINE HilbertTransform( N, Func )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N), INTENT(INOUT)  :: Func
      REAL, DIMENSION(N) :: Transform

      INTEGER :: i, j

      Transform = 0.0
      DO i = 1, N
         DO j = 1, N
            IF (i /= j)  Transform(i) = Transform(i) + Func(j) * real(j) / ( real(j)**2 - real(i)**2  )
         END DO
      END DO
      Func(:) = Transform(:) * 2.0 / MyConsts_PI
   END SUBROUTINE HilbertTransform

END PROGRAM NewWiener
