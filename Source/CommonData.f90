!***************************************************************************************
!*                              MODULE CommonData
!***************************************************************************************
!
!>  \brief     Common data
!>  \details   This module include the data that is shared by the
!>             all the other modules of the code
!
!***************************************************************************************
MODULE CommonData
   IMPLICIT NONE

!    REAL, PARAMETER :: ev_tu = 1.0364314
!    REAL, PARAMETER :: pi = 3.1415926535897932384626433832795028841971 

!> \name MASSES
!> Masses of the atoms
!> @{
   REAL, SAVE :: rmh                  !< mass of the H atom
   REAL, SAVE :: rmc                  !< mass of the C atoms
!> @}

!> \name POSITIONS
!> Positions of the atoms
!> @{
   REAL, SAVE :: xh, yh, zh           !< X,Y,Z positions of the H atom
   REAL, DIMENSION(121), SAVE :: z    !< Z positions of the C atoms
!> @}

!> \name VELOCITIES
!> Positions of the atoms
!> @{
   REAL, SAVE :: vxh, vyh, vzh        !< X,Y,Z velocity of the H atom
   REAL, DIMENSION(121), SAVE :: vzc  !< Z velocity of the C atoms
!> @}

!> \name ACCELERATIONS
!> Positions of the atoms
!> @{
   REAL, SAVE :: axh, ayh, azh        !< X,Y,Z acceleration of the H atom
   REAL, DIMENSION(121), SAVE :: azc  !< Z acceleration of the C atoms
!> @}

!> \name latt
!> ???
!> @{
   REAL, DIMENSION(121), SAVE :: vzsum, zsum
   REAL, DIMENSION(10,501), SAVE :: vz2av, zcav
!> @}

   INTEGER :: nevo                    !<  ????

END MODULE CommonData