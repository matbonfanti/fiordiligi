!***************************************************************************************
!*                              MODULE ClassicalEqMotion
!***************************************************************************************
!
!>  \brief     Subroutines for integration of classical equations of motion
!>  \details   This module contains subroutines to evolve a classical system \n
!>             for one timestep. Some parameters (like timestep, thermostat options ...) \n
!>             are setup in a derived data type and setup as initialization of the module. \n
!>             Variables changing at each  timestep (positions, velocities, accelerations) \n
!>             instead are given as input of the propagation subroutines. \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             20 January 2013
!>
!***************************************************************************************

MODULE ClassicalEqMotion
   USE RandomNumberGenerator
   USE ErrorTrap
   IMPLICIT NONE
   
      PRIVATE
      PUBLIC :: Evolution
      PUBLIC :: EvolutionSetup, SetupThermostat, DisposeThermostat
      PUBLIC :: EOM_VelocityVerlet, EOM_KineticEnergy

      LOGICAL :: GaussianNoise = .TRUE.
      
      !> Evolution datatype, storing all the general data required
      !> for integration of the EoM, with or without Langevin thermostat
      TYPE Evolution
         INTEGER   :: NDoF                            !< Nr of degrees of freedom
         REAL, DIMENSION(:), POINTER :: Mass          !< Vector with the masses of the system
         REAL  :: dt                                  !< Time step of integration
         REAL  :: FrictionCoeff                       !< Coefficient including Langevin friction
         REAL, DIMENSION(:), POINTER :: ThermalNoise  !< Vector of the thermal noise sigma
         LOGICAL :: HasThermostat = .FALSE.           !< Thermostat data has been setup
         LOGICAL :: IsSetup = .FALSE.                 !< Evolution data has been setup
      END TYPE Evolution

   CONTAINS
   
!*******************************************************************************
!          EvolutionSetup
!*******************************************************************************
!> Setup evolution general data, store them in the Evolution datatype.
!>
!> @param EvolutionData    Evolution data type to setup
!> @param NDoF             Integer var with the number of degree of freedom
!> @param MassVector       NDoF-dim real vector with the masses
!*******************************************************************************
   SUBROUTINE EvolutionSetup( EvolutionData, NDoF, MassVector, TimeStep )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)  :: EvolutionData
      INTEGER, INTENT(IN)               :: NDoF
      REAL, DIMENSION(NDoF), INTENT(IN) :: MassVector
      REAL, INTENT(IN)                  :: TimeStep

      ! warn user if overwriting previously setup data
      CALL WARN( EvolutionData%IsSetup, "ClassicalEqMotion.EvolutionSetup: overwriting evolution data" ) 

      ! Store number of degree of freedom
      EvolutionData%NDoF = NDoF
      
      ! If necessary, deallocate array and then store masses
      IF (EvolutionData%IsSetup)  DEALLOCATE( EvolutionData%Mass )
      ALLOCATE( EvolutionData%Mass( NDoF ) )
      EvolutionData%Mass = MassVector
      
      ! Store timestep
      EvolutionData%dt = TimeStep
      
      ! evolutiondata is now setup
      EvolutionData%IsSetup = .TRUE.
      
#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,I3,A,1F8.3)") "Evolution data type is setup: Nr DoF is ",NDoF," and TimeStep ", TimeStep
#endif
      
   END SUBROUTINE EvolutionSetup
   
!*******************************************************************************
!> Setup thermostat data, store them in the Evolution datatype.
!>
!> @param EvolData     Evolution data type to setup
!> @param Gamma        Langevin friction coefficient
!> @param Temperature  Temperature of thermostat
!*******************************************************************************
   SUBROUTINE SetupThermostat( EvolData, Gamma, Temperature )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      REAL, INTENT(IN)                  :: Gamma, Temperature

      INTEGER :: iDoF
      
      ! error if trying to setup thermostat of a non-setup evolution type
      CALL ERROR( .NOT. EvolData%IsSetup, "ClassicalEqMotion.SetupThermostat: evolution data not setup" )
      ! warn user if overwriting previously setup data
      CALL WARN( EvolData%HasThermostat, "ClassicalEqMotion.SetupThermostat: overwriting thermostat data" ) 

      ! Set coefficient for velocity integration with Langevin friction
      EvolData%FrictionCoeff = 1.0 - 0.5 * Gamma * EvolData%dt
      
      ! Set standard deviations of the thermal noise
      IF ( .NOT. EvolData%HasThermostat )  ALLOCATE( EvolData%ThermalNoise( EvolData%NDoF ) )
      DO iDoF = 1, EvolData%NDoF 
           EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*EvolData%Mass(iDoF)*Gamma/EvolData%dt )
      END DO
      
      ! Themostat data is now setup
      EvolData%HasThermostat = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,1F8.3,A,1F8.3)") "Themostat is setup with Gamma = ",Gamma," and Temperature = ", Temperature
#endif
      
   END SUBROUTINE SetupThermostat
   
!*******************************************************************************
!> Dispose thermostat data.
!>
!> @param EvolData     Evolution data type with the thermostat to dispose
!*******************************************************************************
   SUBROUTINE DisposeThermostat( EvolData, Gamma, Temperature )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      REAL, INTENT(IN)                  :: Gamma, Temperature

      INTEGER :: iDoF
      
      ! continue if trying to dispose thermostat that is not setup
      IF ((.NOT. EvolData%IsSetup) .OR. (.NOT. EvolData%HasThermostat))  RETURN

      ! Set coefficient with zero friction
      EvolData%FrictionCoeff = 1.0
      
      ! Deallocate standard deviations of the thermal noise
      DEALLOCATE( EvolData%ThermalNoise )
      
      ! Themostat data is now disposed
      EvolData%HasThermostat = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Themostat has been disposed"
#endif
   END SUBROUTINE DisposeThermostat

   
   
!*******************************************************************************
!> Compute kinetic energy corresponding to a given velocity vector.
!>
!> @param EvolData     Evolution data type
!> @param Velocity     Array containing the velocity at given time step
!*******************************************************************************
   REAL FUNCTION EOM_KineticEnergy( EvolData, Vel ) RESULT( KinEnergy )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Vel
      INTEGER :: iDoF
      
      KinEnergy = 0.0
      DO iDoF = 1, EvolData%NDoF
         KinEnergy = KinEnergy + 0.5 * EvolData%Mass(iDoF) * Vel(iDoF)**2
      END DO

   END FUNCTION EOM_KineticEnergy


!*******************************************************************************
!> Propagate trajectory with Velocity-Verlet algorith.
!> If the Langevin parameters are setup, propagation is done in the
!> canonical ensamble with a Langevin thermostat
!> @ref http://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
!>
!> @param EvolData     Evolution data type
!*******************************************************************************
   SUBROUTINE EOM_VelocityVerlet( EvolData, Pos, Vel, Acc, GetPotential, V )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                :: V

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE
      
      INTEGER :: iDoF
 
      ! (1) FULL TIME STEP FOR THE POSITIONS
       
      DO iDoF = 1, EvolData%NDoF
          Pos(iDoF) = Pos(iDoF) + Vel(iDoF)*EvolData%dt + 0.5*Acc(iDoF)*(EvolData%dt**2)
      END DO
      
      ! (2) HALF TIME STEP FOR THE VELOCITIES

      IF ( .NOT. ( EvolData%HasThermostat ) ) THEN        ! Integration without Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END DO
            
      ELSE IF ( EvolData%HasThermostat ) THEN             ! Integration with Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = EvolData%FrictionCoeff*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END DO
      END IF

      ! (3) NEW FORCES AND ACCELERATIONS 
      
      V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
       
      ! Divide by mass and if it's  Langevin dynamics, add random noise
      IF ( .NOT. ( EvolData%HasThermostat ) ) THEN        ! only potential forces
         DO iDoF = 1, EvolData%NDoF
            Acc(iDoF) = Acc(iDoF)/EvolData%Mass(iDoF)
         END DO
      ELSE IF (( EvolData%HasThermostat ) .AND. ( GaussianNoise )) THEN    ! add gaussian noise
         DO iDoF = 1, EvolData%NDoF
            Acc(iDoF) = ( Acc(iDoF)+GaussianRandomNr( EvolData%ThermalNoise(iDoF) ) ) / EvolData%Mass(iDoF)
         END DO
      ELSE IF (( EvolData%HasThermostat ) .AND. ( .NOT. GaussianNoise )) THEN    ! add uniform noise
         DO iDoF = 1, EvolData%NDoF
            Acc(iDoF) = ( Acc(iDoF) + UniformRandomNr( -sqrt(3.)*EvolData%ThermalNoise(iDoF),     &
                                                       sqrt(3.)*EvolData%ThermalNoise(iDoF) ) )  / EvolData%Mass(iDoF)
         END DO
      END IF
      
      ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES

      IF ( .NOT. ( EvolData%HasThermostat ) ) THEN        ! Integration without Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END DO
            
      ELSE IF ( EvolData%HasThermostat ) THEN             ! Integration with Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = EvolData%FrictionCoeff*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END DO
      END IF

   END SUBROUTINE EOM_VelocityVerlet   
   
! SUBROUTINE beeman(pos,vel,force,force1,cell,link,ncel,rcel,len,potential,ksi,trace)
! !integrate the EOM using Beeman's predictor-corrector algorithm with Nose-
! Hoover thermostat
! USE constants
! IMPLICIT NONE
! REAL*8 pos(dim,num)
! REAL*8 vel(dim,num),velp(dim,num) !velp:the predicted velocity
! !forces at the previous, current and the next time point
! REAL*8 force(dim,num),force1(dim,num),force2(dim,num)
! INTEGER cell(0:ncel-1,0:ncel-1,0:ncel-1) !cell list
! INTEGER link(num)!link list
! INTEGER ncel !number of cells per dimension
! REAL*8 rcel,len,potential(num)
! REAL*8 ksi !friction factor
! INTEGER i,j !loop control
! REAL*8 trace(dim,num) !trace of the atoms
! REAL*8 delta_pos !relative motion
! DO i = 1,num
! DO j = 1,dim
! !calculate relative motion
! delta_pos = vel(j,i)*time_step +
! time_step*time_step/6*(4*force1(j,i)-force(j,i))/mass
! pos(j,i) = pos(j,i) + delta_pos
! trace(j,i) = trace(j,i) + delta_pos
! !periodic boundary conditions
! IF (pos(j,i) .LT. 0) pos(j,i) = pos(j,i) - INT(pos(j,i)/len)*len + len
! IF (pos(j,i) .GT. len) pos(j,i) = pos(j,i) - INT(pos(j,i)/len)*len
! !predicted velocities
! velp(j,i) = vel(j,i) + time_step/2*(3*force1(j,i)-force(j,i))/mass
! END DO
! END DO
! !update cell list and calculate forces
! CALL cell_list(pos,cell,ncel,rcel,link)
! CALL force_calculation (pos,cell,link,ncel,rcel,force2,len,potential)
! DO i = 1,num
! DO j = 1,dim
! !calculate forces incorporating thermostat
! force2(j,i)=force2(j,i)-ksi*mass*velp(j,i)
! !corrected velocity
! vel(j,i) = vel(j,i) + time_step/6*(2*force2(j,i)+5*force1(j,i)-
! force(j,i))/mass
! END DO
! END DO
! RETURN
! END SUBROUTINE beeman

END MODULE ClassicalEqMotion