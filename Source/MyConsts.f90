!***************************************************************************************
!*                              MODULE MyConsts
!***************************************************************************************
!
!>  \brief     Often used mathematical and physical constants
!>  \details   This class defines commonly used mathematical and physical constants,
!>             and the numerical precision of the machine, that can be set dynamically
!>             at runtime with the subroutine Calculate_Constant_EPS().  \n
!>             The module also contains the integer function LookForFreeUnit(), that
!>             gives back the smallest unit that is not open for i/o. 
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             June 2012
!>
!***************************************************************************************
!
!>  \pre              Call the subroutine Calculate_Constant_EPS to set the 
!>                    machine precision  at runtime 
!
!***************************************************************************************
!
!>   \remark         Sources of the constants: \n 
!>                   1) The NIST Reference on Constants, Units and
!>                      Uncertainty ( http://physics.nist.gov/cuu/index.html ) \n
!>                   2) The NIST Atomic Weights and Isotopic
!>                      Compositions ( http://www.nist.gov/pml/data/comp.cfm )  \par
!>   \remark         The module is based on Mark Somers' constants.f90
!
!***************************************************************************************
MODULE MyConsts
   USE ErrorTrap

   IMPLICIT NONE

!> \name NUMERICAL PRECISION
!> Numerical precision of the machine \n
!> It can be set dynamically through the subroutine Calculate_Constant_EPS() 
!> @{
REAL               :: MyConsts_EPS      = 1E-12
!> @}

!> \name MATHEMATICAL CONSTANTS
!> Values of frequently used mathematical constants
!> @{
REAL, PARAMETER    :: MyConsts_SQRT_2    = SQRT( 2.0 )                                   !< Square Root of 2
REAL, PARAMETER    :: MyConsts_SQRT_3    = SQRT( 3.0 )                                   !< Square Root of 3
REAL, PARAMETER    :: MyConsts_INVSQRT_3 = 1./MyConsts_SQRT_3                            !< Inverse of the square root of 3
REAL, PARAMETER    :: MyConsts_PI        = 3.1415926535897932384626433832795028841971    !< Greek Pi
COMPLEX, PARAMETER :: MyConsts_I         = CMPLX( 0.0, 1.0 )                             !< Imaginary Unit
!> @}

!> \name PHYSICAL CONSTANTS
!> Values of frequently used physical constants
!> @{
REAL, PARAMETER    :: MyConsts_kb = 8.617332478e-5              !<  Boltzmann's constant in eV/K
REAL, PARAMETER    :: MyConsts_uma = 1.660538921e-27            !<  Atomic Mass Unit (AMU) in kg (from NIST reference)
REAL, PARAMETER    :: MyConsts_mel = 9.10938291e-31             !<  Electron mass in kg (from NIST reference)
REAL, PARAMETER    :: MyConsts_NAvo = 6.02214129E23             !<  Avogadro's number (from NIST reference)
REAL, PARAMETER    :: MyConsts_ThermoCal = 4.184                !<  Thermochemical Calorie ( from NIST reference )
REAL, PARAMETER    :: MyConsts_AUofTime = 2.418884326502e-2     !<  Atomic unit of time in fs ( from NIST reference )
REAL, PARAMETER    :: MyConsts_SpeedOfLight = 299792458         !<  Speed of light in the vacuum, in m/s ( from NIST reference )

!> @}

!> \name CONVERSION FACTORS
!> Conversions factor between common units
!> @{
REAL, PARAMETER    :: MyConsts_Hartree2eV     = 27.21138505                 !< Conversion factor from Hartree to eV (from NIST reference)
REAL, PARAMETER    :: MyConsts_Hatree2joule   = 4.35974434e-18              !< Conversion factor from Hartree to Joule (from NIST reference)
REAL, PARAMETER    :: MyConsts_Hatree2kjmol   = MyConsts_Hatree2joule*MyConsts_NAvo/1000.   !< Conversion factor from Hartree to kJ/mol
REAL, PARAMETER    :: MyConsts_Hatree2kcalmol = MyConsts_Hatree2kjmol/ MyConsts_ThermoCal   !< Conversion factor from Hartree to kcal/mol
REAL, PARAMETER    :: MyConsts_Deg2Rad        = MyConsts_PI / 180.0         !< Conversion factor from Decimal degrees to radiants
REAL, PARAMETER    :: MyConsts_Bohr2Ang       = 0.52917721092               !< Conversion factor from Bohr to Ang (from NIST reference) 
REAL, PARAMETER    :: MyConsts_Uma2Au         = MyConsts_uma / MyConsts_mel !< Conversion factor from UMA to Atomic Units
REAL, PARAMETER    :: MyConsts_fs2AU          = 1.0 / MyConsts_AUofTime     !< Conversion factor from femtosecond to Atomic Units
REAL, PARAMETER    :: MyConsts_K2AU           = MyConsts_kb / MyConsts_Hartree2eV   !< Conversion factor from Kelvin to Atomic Units
REAL, PARAMETER    :: MyConsts_cmmin1tofsmin1 = 2.0*MyConsts_PI*MyConsts_SpeedOfLight*1e-13  !< Conversion from cm-1 to fs-1 (2PI for angular freq conversion)
REAL, PARAMETER    :: MyConsts_cmmin1toAU     = MyConsts_cmmin1tofsmin1/MyConsts_fs2AU !< Conversion from cm-1 to Atomic Units (freq=energy)
!> @}

!> \name ELEMENTS MASSES
!>  Masses of the elements (both pure and weighted for the isotopic composition) \n
!>  Values have been taken from NIST:\n http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
!> @{
REAL, PARAMETER    :: MyConsts_mH  =  1.00782503207 * MyConsts_Uma2Au    !< Hydrogen atomic mass (from NIST)
REAL, PARAMETER    :: MyConsts_mD  =  2.0141017778 * MyConsts_Uma2Au     !< Deuterium atomic mass (from NIST)
REAL, PARAMETER    :: MyConsts_mCMix =  12.01078 * MyConsts_Uma2Au       !< Isotopic mix Carbon atomic mass (from NIST)
REAL, PARAMETER    :: MyConsts_mCuMix = 63.5463 * MyConsts_Uma2Au        !< Isotopic mix Cupper atomic mass (from NIST)
!> @}

!> \name OTHER
!> Other useful parameters, like special characters or numeric standard types \n
!> @{
CHARACTER(len=1), PARAMETER  :: NewLine = achar(10)
! The following kind parameters define the standard numeric kind, REGARDLESS of compiler options
! and coherently with the kind definition of the compiler
! hence SINGLE_PRECISION_KIND will be the kind of a true single precision real data for the compiler used
! and analogously for the other parameters
INTEGER, PARAMETER           :: SINGLE_PRECISION_KIND = SELECTED_REAL_KIND(5,20)
INTEGER, PARAMETER           :: DOUBLE_PRECISION_KIND = SELECTED_REAL_KIND(10,40)
INTEGER, PARAMETER           :: LONG_INTEGER_KIND     = SELECTED_INT_KIND( 8 )
INTEGER, PARAMETER           :: SHORT_INTEGER_KIND    = SELECTED_INT_KIND( 16 )
!> @}

CONTAINS

!*******************************************************************************
! Calculate_Constant_EPS
!*******************************************************************************
!>  Calculates the machines precision at runtime and stores 
!>  it into the Constant_EPS variable so it can be used on any system
!>  or byte sizes of real types...
!*******************************************************************************
   SUBROUTINE Calculate_Constant_EPS( )
   IMPLICIT NONE
   REAL A

      MyConsts_EPS = 1.0
      DO 
         A = 1.00 - MyConsts_EPS
         IF ( A >= 1.00 ) EXIT
         MyConsts_EPS = MyConsts_EPS / 2.00
      END DO 

      MyConsts_EPS = ABS( MyConsts_EPS )
      IF( MyConsts_EPS <= 0.0000000 ) MyConsts_EPS = 1.0000E-12   ! added this to be sure ...

   END SUBROUTINE Calculate_Constant_EPS

!*******************************************************************************
! LookForFreeUnit
!*******************************************************************************
!>  Set the smallest integer equal or greater than 20 that is 
!>  available as unit i/o. 
!> 
!> @returns           The integer number of the smallest free unit.
!*******************************************************************************
   INTEGER FUNCTION LookForFreeUnit()
   IMPLICIT NONE
      LOGICAL               :: File_Opened
      INTEGER, PARAMETER    :: UnitMax = 300

      DO LookForFreeUnit = 20, UnitMax   ! Look for a non opened I/O unit
         INQUIRE( UNIT=LookForFreeUnit, OPENED=File_Opened )
         IF (.NOT. File_Opened) EXIT
      END DO
      ! If an available unit has been found, use as function results 
      CALL ERROR((LookForFreeUnit == UnitMax), &
                                      "PrintTools: No free I/O unit available")
   END FUNCTION LookForFreeUnit

!*******************************************************************************
! NumberToString
!*******************************************************************************
!>  Given an integer number within 0 and 999, five back a 3 character string
!>  of the number with zeros instead of null characters.
!> 
!> @param    N        Input integer number.
!> @returns           A string with the number.
!*******************************************************************************
   FUNCTION NumberToString( N )      RESULT( String )
   IMPLICIT NONE
      INTEGER, INTENT(IN)   :: N
      CHARACTER(3)          :: String

      SELECT CASE ( N )
         CASE ( 0:9 )
            WRITE(String,"(A2,I1)")  "00", N
         CASE ( 10:99 )
            WRITE(String,"(A1,I2)")  "0", N
         CASE ( 100:999 )
            WRITE(String,"(I3)")      N
         CASE DEFAULT
            CALL AbortwithError( "NumberToString: N too large " )
      END SELECT

   END FUNCTION NumberToString
   
!*******************************************************************************
! CountLinesInFile
!*******************************************************************************
!>  Count the number of lines in a file with a given filename.
!> 
!> @param    FileName     Input string with the file name.
!> @returns               The integer number of rows of file FileName.
!*******************************************************************************
   INTEGER FUNCTION CountLinesInFile( FileName )
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: FileName
      INTEGER             :: Stat, UnitNr
      CHARACTER(1)        :: String

      UnitNr=LookForFreeUnit()
      OPEN( FILE=TRIM(ADJUSTL(FileName)), UNIT=UnitNr )
      DO CountLinesInFile = 0,100000
         READ( UnitNr, *, IOSTAT=Stat ) String
         IF (Stat /= 0) EXIT
      END DO
      CLOSE( UnitNr )

   END FUNCTION CountLinesInFile



END MODULE MyConsts

!********************************************* END OF FILE *******************************
