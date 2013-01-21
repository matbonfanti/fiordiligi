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
      PUBLIC :: EOM_VelocityVerlet, EOM_KineticEnergy, EOM_Beeman

      LOGICAL :: GaussianNoise = .TRUE.
      
      !> Evolution datatype, storing all the general data required
      !> for integration of the EoM, with or without Langevin thermostat
      TYPE Evolution
         INTEGER   :: NDoF                            !< Nr of degrees of freedom
         REAL, DIMENSION(:), POINTER :: Mass          !< Vector with the masses of the system
         REAL  :: dt                                  !< Time step of integration
         REAL  :: FrictionCoeff_FullDt                !< Coefficient including Langevin friction for one timestep
         REAL  :: FrictionCoeff_HalfDt                !< Coefficient including Langevin friction for half timestep
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
      EvolData%FrictionCoeff_HalfDt = 1.0 - 0.5 * Gamma * EvolData%dt
      EvolData%FrictionCoeff_FullDt = 1.0 - Gamma * EvolData%dt
      
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
      EvolData%FrictionCoeff_HalfDt = 1.0
      EvolData%FrictionCoeff_FullDt = 1.0
      
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
            Vel(iDoF) = EvolData%FrictionCoeff_HalfDt*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
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
            Vel(iDoF) = EvolData%FrictionCoeff_HalfDt*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END DO
      END IF

   END SUBROUTINE EOM_VelocityVerlet   
   





!*******************************************************************************
!> Propagate trajectory with Beeman's algorith.
!> If the Langevin parameters are setup, propagation is done in the
!> canonical ensamble with a Langevin thermostat
!> @ref 
!>
!> @param EvolData     Evolution data type
!*******************************************************************************
   SUBROUTINE EOM_Beeman( EvolData, Pos, Vel, Acc, PreAcc, GetPotential, V )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc, PreAcc
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
          Pos(iDoF) = Pos(iDoF) + Vel(iDoF)*EvolData%dt + (4.*Acc(iDoF)-PreAcc(iDof))*(EvolData%dt**2)/6.0
      END DO
      
      ! (2) VELOCITIES: TIME STEP INCLUDING ACCELERATION AT PRESENT AND PREVIOUS T

      IF ( .NOT. ( EvolData%HasThermostat ) ) THEN        ! Integration without Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = Vel(iDoF) +  ( 5.0*Acc(iDoF) - PreAcc(iDof) ) * EvolData%dt / 6.0
         END DO
            
      ELSE IF ( EvolData%HasThermostat ) THEN             ! Integration with Langevin thermostat
         DO iDoF = 1, EvolData%NDoF
            Vel(iDoF) = EvolData%FrictionCoeff_FullDt * Vel(iDoF) +  ( 5.0*Acc(iDoF) - PreAcc(iDof) ) * EvolData%dt / 6.0
         END DO
      END IF

      ! (3) NEW FORCES AND ACCELERATIONS
 
      PreAcc(:) = Acc(:)                 ! Store present accelerations for later use
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
      
      ! (4) VELOCITIES: TIME STEP INCLUDING ACCELERATION AT NEW T 

      DO iDoF = 1, EvolData%NDoF
         Vel(iDoF) = Vel(iDoF) + Acc(iDoF)*EvolData%dt / 3.0
      END DO  
    
   END SUBROUTINE EOM_Beeman   



END MODULE ClassicalEqMotion