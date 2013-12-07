!***************************************************************************************
!*                              MODULE PrintTraj
!***************************************************************************************
!
!>  \brief         Plot of a nuclear trajectory in XYZ format
!>  \details      ...................
!>                \arg  .........
!
!***************************************************************************************
!
!>  \author        Matteo Bonfanti
!>  \version       1.0
!>  \date          March 2012
!>
!***************************************************************************************
!
!>  \par 
!>  \arg 
!
!***************************************************************************************
!
!>  \todo 
!
!***************************************************************************************
MODULE PrintTraj
   USE ErrorTrap
   USE MyConsts

   IMPLICIT NONE

   PRIVATE

   !> Trajectory data type
   TYPE :: TrajOutput
      PRIVATE
      LOGICAL           :: IsOpen = .FALSE.     !< Actual status of the trajectory output
      CHARACTER(len=30) :: FileName             !< Name of the output file 
      INTEGER           :: NrOfAtoms            !< Nr of atoms in the trajectory
      INTEGER           :: Counter = 0          !< Counter of the snapshots being written
      INTEGER           :: Unit                 !< The unit to which it is being written
      CHARACTER(len=3), DIMENSION(:), POINTER  :: AtomLabels  !< Labels of the atoms 
   END TYPE TrajOutput



   
!********************************************************************************************************
   CONTAINS
!********************************************************************************************************



!*******************************************************************************
!          Traj_NewTrajectory 
!*******************************************************************************
!> Initialize a trajectory file.
!> In detail, the subroutine initialize the data and open the output
!> unit to print the trajectory 
!>
!> @param TrajData           Trajectory data type
!> @param CollectionName     Name of the output file
!> @param CollectionSize     Number of atoms for each snapshot
!*******************************************************************************
   SUBROUTINE Traj_NewTrajectory ( TrajData, FileName, AtomsLabels )
      IMPLICIT NONE
      TYPE(TrajOutput), INTENT(inout)            :: TrajData
      CHARACTER(len=*), INTENT(in)               :: FileName
      CHARACTER(len=*), DIMENSION(:), INTENT(in) :: AtomsLabels
      INTEGER :: i

      ! Check if current data is already being written
      CALL ERROR( TrajData%IsOpen, "Traj_NewTrajectory: TrajData already begin written" )

      ! Store file name and atom number
      TrajData%FileName = TRIM(ADJUSTL( FileName ))
      TrajData%NrOfAtoms = size(AtomsLabels)

      ! Store atomic labels
      ALLOCATE( TrajData%AtomLabels(TrajData%NrOfAtoms) )
      DO i=1, TrajData%NrOfAtoms
         TrajData%AtomLabels(i) = TRIM(ADJUSTL(AtomsLabels(i)))
      END DO

      ! Open output unit
      TrajData%Unit = LookForFreeUnit()
      OPEN( FILE=TRIM(ADJUSTL(TrajData%FileName)), UNIT=TrajData%Unit )

      ! Initialize snapshot counter
      TrajData%Counter = 0

      ! Set actual status 
      TrajData%IsOpen = .TRUE.
      
   END SUBROUTINE Traj_NewTrajectory


!*******************************************************************************
!          Traj_EndTrajectory 
!*******************************************************************************
!> Finalize a trajectory file.
!>
!> @param TrajData           Trajectory data type
!*******************************************************************************
   SUBROUTINE Traj_EndTrajectory ( TrajData, FileName, NrOfAtoms )
      IMPLICIT NONE
      TYPE(TrajOutput), INTENT(inout)   :: TrajData
      CHARACTER(len=*), INTENT(in)      :: FileName
      INTEGER, INTENT(in)               :: NrOfAtoms
   
      ! Check if current data is already being written
      CALL WARN( .NOT. TrajData%IsOpen, "Traj_EndTrajectory: TrajData is not open" )

      ! Close output unit
      CLOSE( UNIT=TrajData%Unit )

      ! Deallocate memory
      DEALLOCATE( TrajData%AtomLabels )

      ! Set actual status 
      TrajData%IsOpen = .FALSE.
      
   END SUBROUTINE Traj_EndTrajectory


!*******************************************************************************
!          Traj_WriteSnapshot
!*******************************************************************************
!> Initialize a new VTR file and set the rectilinear mesh.
!>
!> @param VTKType            VTK collection data type
!*******************************************************************************
   SUBROUTINE Traj_WriteSnapshot( TrajData, Coordinates )
      IMPLICIT NONE
         TYPE(TrajOutput), INTENT(inout)   :: TrajData
         REAL, DIMENSION(:,:), INTENT(in)  :: Coordinates

         INTEGER :: i

         ! check and store array dimensions
         CALL ERROR( size( Coordinates, 1 ) /= TrajData%NrOfAtoms, "Traj_WriteSnapshot: invalid number of atoms ")
         CALL ERROR( size( Coordinates, 2 ) /= 3 , " Traj_WriteSnapshot: Invalid nr of dimensions "  )

         ! check that the trajdata is open
         CALL ERROR( .NOT. TrajData%IsOpen, "Traj_WriteSnapshot: TrajData is not open" )

         ! increment snapshot counter
         TrajData%Counter = TrajData%Counter + 1

         ! print snapshot
         WRITE( TrajData%Unit, "(I4)" ) TrajData%NrOfAtoms
         WRITE( TrajData%Unit, * ) " Step nr ",TrajData%Counter," of the trajectory "
         DO i = 1, TrajData%NrOfAtoms 
           WRITE( TrajData%Unit, "(A1,A3,A1,3F15.6)" ) " ",TrajData%AtomLabels(i)," ",Coordinates(1:3,i)
         END DO
   
   END SUBROUTINE Traj_WriteSnapshot


   
END MODULE PrintTraj
