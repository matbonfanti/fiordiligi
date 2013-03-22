!***************************************************************************************
!*                      MODULE IndependentOscillatorsModel
!***************************************************************************************
!
!>  \brief      Independent oscillator models
!>  \details    This class define the necessary subroutines to define and compute
!>              the potential and the derivatives of a discrete independent oscillator
!>              model, both in the standard and in the chain representation
!
!***************************************************************************************
!
!>  \author     Matteo Bonfanti
!>  \version    1.0
!>  \date       8 March 2013
!>
!***************************************************************************************
!
!>  \pre        To use the class the potential needs to be setup with the 
!>              couplings and the frequencies of the oscillators
!>              This data is read from input file by the setup subroutine      
!
!***************************************************************************************
!
!>  \remark     
!
!***************************************************************************************
!
!>  \todo     * Implement the bath in chain mode
!>            * Fix the units of input data for bath setup
!
!***************************************************************************************

MODULE IndependentOscillatorsModel
   USE ErrorTrap
   USE MyConsts
   USE PotentialModule
   USE SplineInterpolator

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupIndepOscillatorsModel

   INTEGER, PARAMETER, PUBLIC :: STANDARD_BATH = 0        ! < normal bath, in which all the oscillators are coupled to the system
   INTEGER, PARAMETER, PUBLIC :: CHAIN_BATH    = 1        ! < chain bath, in which all the oscillators are coupled in chain
    
   ! > Integer parameter to store the bath type (chain, standard ... )
   INTEGER :: BathType 

   ! > Number of oscillators of the bath
   INTEGER :: BathSize
   ! > Harmonic frequencies of the bath (stored in AU)
   REAL, DIMENSION(:), ALLOCATABLE :: Frequencies
   ! > Coupling of the bath oscillators (stored in AU)
   REAL, DIMENSION(:), ALLOCATABLE :: Couplings
   ! > Cutoff frequency of the normal bath (stored in AU)
   REAL :: CutOff
   ! > Frequency spacing of the normal bath (stored in AU)
   REAL :: DeltaOmega
   ! > Mass of the oscillators
   REAL :: OscillatorsMass

   ! > Logical variable to set status of the class
   LOGICAL :: BathIsSetup = .FALSE.


CONTAINS   

! ----------------------------------------------------------------------------------------------------------------

!*******************************************************************************
! SetupIndepOscillatorsModel
!*******************************************************************************
!>  This subroutine can be used to setup the parameters of the oscillator bath
!>  reading the couplings and the frequencies from input file
!> 
!> @param    FileName       Name of the file with the input parameters.
!*******************************************************************************
   SUBROUTINE SetupIndepOscillatorsModel( N, SetBathType, FileName, Mass, CutOffFreq )
      IMPLICIT NONE
      INTEGER, INTENT(IN)        :: N
      INTEGER, INTENT(IN)        :: SetBathType
      CHARACTER(*), INTENT(IN)   :: FileName
      REAL, INTENT(IN)           :: Mass
      REAL, OPTIONAL, INTENT(IN) :: CutOffFreq

      INTEGER :: InpUnit, RdStatus
      INTEGER :: iBath
      INTEGER :: NData
      REAL, DIMENSION(:), ALLOCATABLE :: RdFreq, RdSpectralDens
      TYPE(SplineType) :: SpectralDensitySpline

      ! If data already setup give a warning and deallocate memory
      CALL WARN( BathIsSetup, "IndependentOscillatorsModel.SetupIndepOscillatorsModel: overwriting bath data" )
      IF (BathIsSetup) CALL DisposeIndepOscillatorsModel( )

      ! Store number of bath degrees of freedom
      BathSize = N
      ! Set the type of bath
      BathType = SetBathType
      ! Set the mass of the oscillators
      OscillatorsMass = Mass

      ! Allocate memory
      ALLOCATE( Frequencies(BathSize), Couplings(BathSize) )

      IF ( BathType == CHAIN_BATH ) THEN

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read frequencies and couplings
            DO iBath = 1, BathSize
               READ(InpUnit,*, IOSTAT=RdStatus) Frequencies(iBath), Couplings(iBath)
               ! if data is missing, a markovian continuation of the chain is assumed 
               IF ( RdStatus /= 0 ) THEN
                  Frequencies(iBath) = Frequencies(iBath-1)
                  Couplings(iBath) = Couplings(iBath-1)
               END IF
            END DO

            ! Close input file
            CLOSE( InpUnit )

      ELSE IF ( BathType == STANDARD_BATH ) THEN

            ! Count the number of data in the input file
            NData = CountLinesInFile( FileName )

            ! Allocate temporary arrays to store the input spectral density
            ALLOCATE( RdFreq(NData), RdSpectralDens(NData) )

            ! Open input file
            InpUnit = LookForFreeUnit()
            OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

            ! Read frequencies and couplings
            DO iBath = 1, NData
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !  FORSE LE UNITA' DI MISURA SONO DA AGGIUSTARE ????
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               READ(InpUnit,*) RdFreq(iBath), RdSpectralDens(iBath)
            END DO

            ! Close input file
            CLOSE( InpUnit )

            ! Set cutoff frequency
            IF ( PRESENT( CutOffFreq )) THEN
               CutOff = CutOffFreq               ! if is given from input
            ELSE 
               CutOff = RdSpectralDens( NData )  ! otherwise is the last freq in the input file
            ENDIF
            DeltaOmega = CutOff / real( BathSize)

            ! Spline interpolation of the input spectral density
            CALL SetupSpline( SpectralDensitySpline, RdFreq, RdSpectralDens)

            ! Compute frequencies and couplings 
            DO iBath = 1, BathSize
               Frequencies(iBath) = iBath * DeltaOmega
               Couplings(iBath) = SQRT( 2.0 * OscillatorsMass * Frequencies(iBath) * DeltaOmega *      & 
                    GetSpline( SpectralDensitySpline, Frequencies(iBath) ) / MyConsts_PI )
            ENDDO

            ! Deallocate memory for spectral density interpolation
            CALL DisposeSpline ( SpectralDensitySpline )
            DEALLOCATE( RdFreq, RdSpectralDens )

      ENDIF

      ! Module is setup
      BathIsSetup = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model potential has been setup"
      WRITE(*,*) " ... (details) ... "
#endif

   END SUBROUTINE SetupIndepOscillatorsModel


! ----------------------------------------------------------------------------------------------------------------

!*******************************************************************************
!> System potential plus bath of independent oscillators (normal or in chain form)
!> @ref RIFERIMENTI PER SURFACE OSCILLATOR MODEL
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom, 1 for 
!>                     the nearest Cu and N for the bath
!> @param Forces       Output array with the derivatives of the potential (in au)
!> @param vv           output potential in atomic units
!*******************************************************************************     
      REAL FUNCTION SystemAndIndepedentOscillators( Positions, Forces ) RESULT(V) 
         IMPLICIT NONE

         REAL, DIMENSION(:), INTENT(IN)  :: Positions
         REAL, DIMENSION(:), INTENT(OUT) :: Forces 

         INTEGER :: IBath
         REAL    :: Coupl

         ! Check if the bath is setup
         CALL ERROR( .NOT. BathIsSetup, "IndependentOscillatorsModel.SystemAndIndepedentOscillators: bath is not setup" )

         ! Check the size of the position and forces arrays
         CALL ERROR( size(Positions) /= size(Forces), &
                        "IndependentOscillatorsModel.SystemAndIndepedentOscillators: array dimension mismatch" )
         CALL ERROR( size(Positions) /= BathSize+4, &
                        "IndependentOscillatorsModel.SystemAndIndepedentOscillators: wrong bath size" )

         ! 4D POTENTIAL OF THE SYSTEM
         V = CH_4Dimensional( Positions(1:4), Forces(1:4) ) 

         IF ( BathType == CHAIN_BATH ) THEN

               CALL ERROR ( 1>0, " IndependentOscillatorsModel.SystemAndIndepedentOscillators: chain modo not implemented" )

         ELSE IF ( BathType == STANDARD_BATH ) THEN

               ! POTENTIAL OF THE BATH OSCILLATORS PLUS COUPLING
               Coupl = 0.0
               DO iBath = 1, BathSize
                  V = V + 0.5 * OscillatorsMass * ( Frequencies(iBath) * Positions(4+iBath) )**2
                  Forces(iBath) = - OscillatorsMass * ( Frequencies(iBath) )**2 * Positions(4+iBath) + Positions(4) * Couplings(iBath)
                  Coupl = Coupl + Couplings(iBath) * Positions(4+iBath)
               END DO
               V = V + Coupl * Positions(4)
               Forces(4) = Forces(4) + Coupl

         END IF


      END FUNCTION SystemAndIndepedentOscillators


!*******************************************************************************
! DisposeIndepOscillatorsModel
!*******************************************************************************
!>  Deallocate memory used by this module
!*******************************************************************************
   SUBROUTINE DisposeIndepOscillatorsModel( )
      IMPLICIT NONE

      ! Do nothing if bath is not yet setup
      IF (.NOT. BathIsSetup)  RETURN

      ! deallocate memory
      DEALLOCATE( Frequencies, Couplings )

      ! Module is no longer setup
      BathIsSetup = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model data has been disposed"
#endif

   END SUBROUTINE DisposeIndepOscillatorsModel

! ----------------------------------------------------------------------------------------------------------------





END MODULE IndependentOscillatorsModel
