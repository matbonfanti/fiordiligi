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
!>            * Smart action when not enough data is provided in normal bath
!
!***************************************************************************************

MODULE IndependentOscillatorsModel
   USE ErrorTrap
   USE MyConsts
   USE PotentialModule

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupIndepOscillatorsModel

   INTEGER, PARAMETER :: STANDARD_BATH = 0        ! < normal bath, in which all the oscillators are coupled to the system
   INTEGER, PARAMETER :: CHAIN_BATH    = 1        ! < chain bath, in which all the oscillators are coupled in chain
    
   ! > Integer parameter to store the bath type (chain, standard ... )
   INTEGER :: BathType 

   ! > Maximum number of oscillators of the bath (defined by the nr of freq.s and coupl.s read from input file)
   INTEGER :: BathMaxSize
   ! > Harmonic frequencies of the bath (stored in AU)
   REAL, DIMENSION(:), ALLOCATABLE :: Frequencies
   ! > Coupling of the bath oscillators (stored in AU)
   REAL, DIMENSION(:), ALLOCATABLE :: Couplings

   ! > Logical variable to set status of the class
   LOGICAL :: BathIsSetup = .FALSE.


   LOGICAL, PARAMETER :: CollinearPES = .TRUE.

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
   SUBROUTINE SetupIndepOscillatorsModel( FileName )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: FileName
      INTEGER :: InpUnit
      INTEGER :: iBath
      INTEGER :: RdStatus

      ! If data already setup give a warning
      CALL WARN( BathIsSetup, "IndependentOscillatorsModel.SetupIndepOscillatorsModel: overwriting bath data" )

      ! Open input file
      InpUnit = LookForFreeUnit()
      OPEN( File = TRIM(ADJUSTL(FileName)), Unit = InpUnit )

      ! Read intestation of the file
      READ( InpUnit, * ) BathMaxSize, BathType

      ! Allocate memory
      ALLOCATE( Frequencies(BathMaxSize), Couplings(BathMaxSize) )

      ! Read frequencies and couplings
      DO iBath = 1, BathMaxSize
         READ(InpUnit,*, IOSTAT=RdStatus) Frequencies(iBath), Couplings(iBath)
         ! if data is missing, a markovian continuation of the chain is assumed 
         IF ( RdStatus /= 0 ) THEN
            IF (BathType == CHAIN_BATH) THEN
               Frequencies(iBath) = Frequencies(iBath-1)
               Couplings(iBath) = Couplings(iBath-1)
            ELSE
               CALL ERROR( (1>0), "IndependentOscillatorsModel.SetupIndepOscillatorsModel: missing input data " )
            END IF
         END IF
      END DO

      ! Close input file
      CLOSE( InpUnit )

      ! Module is setup
      BathIsSetup = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,*) " Independent oscillator model potential has been setup"
      WRITE(*,*) " ... (details) ... "
#endif

   END SUBROUTINE SetupIndepOscillatorsModel


! ----------------------------------------------------------------------------------------------------------------


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
