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
      PUBLIC :: EOM_VelocityVerlet, EOM_KineticEnergy, EOM_Beeman, EOM_ImpulseIntegrator

      LOGICAL :: GaussianNoise = .TRUE.
      
      !> Evolution datatype, storing all the general data required
      !> for integration of the EoM, with or without Langevin thermostat
      TYPE Evolution
         INTEGER   :: NDoF                            !< Nr of degrees of freedom
         REAL, DIMENSION(:), POINTER :: Mass          !< Vector with the masses of the system
         REAL  :: dt                                  !< Time step of integration
         REAL  :: Gamma                            !< Integrated Langevin friction for one timestep
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

      EvolutionData%HasThermostat = .FALSE.
      
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
      EvolData%Gamma = Gamma
      
      ! Set standard deviations of the thermal noise
!       IF ( .NOT. EvolData%HasThermostat )  ALLOCATE( EvolData%ThermalNoise( EvolData%NDoF ) )
!       DO iDoF = 1, EvolData%NDoF 
!            EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*EvolData%Mass(iDoF)*Gamma/EvolData%dt )
!            EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*Gamma/EvolData%dt )
!       END DO
      IF ( .NOT. EvolData%HasThermostat )  ALLOCATE( EvolData%ThermalNoise( EvolData%NDoF ) )
      DO iDoF = 1, 4                               ! H atom and C1 atom
           EvolData%ThermalNoise(iDoF) = 0.0
      END DO                                       ! other carbon atoms
      DO iDoF = 5, EvolData%NDoF 
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
      EvolData%Gamma = 0.0
      
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
!> NOTE: THIS INTEGRATOR IS BETTER SUITED FOR MICROCANONICAL DYNAMICS 
!>       in case of Langevin MD, use Beeman's algorithm!
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
      Pos(:) = Pos(:) + Vel(:)*EvolData%dt + 0.5*Acc(:)*(EvolData%dt**2)
 
      IF ( .NOT. EvolData%HasThermostat ) THEN        ! Integration without Langevin thermostat
   
         ! (2) HALF TIME STEP FOR THE VELOCITIES
         Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt

         ! (3) NEW FORCES AND ACCELERATIONS 
         V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
         Acc(:) = Acc(:)/EvolData%Mass(:)   ! only potential forces

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
         Vel(:) = Vel(:) + 0.5*Acc(:)*EvolData%dt


      ELSE IF ( ( EvolData%HasThermostat ) ) THEN     ! Integration with Langevin thermostat
            
         ! (2) HALF TIME STEP FOR THE VELOCITIES
         Vel(:) = EvolData%FrictionCoeff_HalfDt*Vel(:) + 0.5*Acc(:)*EvolData%dt

         ! (3) NEW FORCES AND ACCELERATIONS 
         V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value

         IF ( GaussianNoise ) THEN    ! add gaussian noise
            DO iDoF = 1, EvolData%NDoF
               Acc(iDoF) = ( Acc(iDoF)+GaussianRandomNr( EvolData%ThermalNoise(iDoF) ) ) / EvolData%Mass(iDoF)
            END DO
         ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise
            DO iDoF = 1, EvolData%NDoF
               Acc(iDoF) = ( Acc(iDoF) + UniformRandomNr( -sqrt(3.)*EvolData%ThermalNoise(iDoF),     &
                                                         sqrt(3.)*EvolData%ThermalNoise(iDoF) ) )  / EvolData%Mass(iDoF)
            END DO
         END IF

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
         Vel(:) = EvolData%FrictionCoeff_HalfDt*Vel(:) + 0.5*Acc(:)*EvolData%dt

      END IF

   END SUBROUTINE EOM_VelocityVerlet   
   


!*******************************************************************************
!> Propagate trajectory with Beeman's algorith.
!> If the Langevin parameters are setup, propagation is done in the
!> canonical ensamble with a Langevin thermostat
!> NOTE: THIS INTEGRATOR IS BETTER SUITED FOR LANGEVIN DYNAMICS 
!>       in case of microcanonical MD, use Velocity-Verlet!
!> @ref http://en.wikipedia.org/wiki/Beeman%27s_algorithm
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
      ! Temporary array for predicted velocity and new accelerations
      REAL, DIMENSION( EvolData%NDoF ) :: NewPos, NewVel, NewAcc

      ! (1) PREDICTED POSITIONS
      NewPos(:) = Pos(:) + Vel(:)*EvolData%dt + (4.*Acc(:)-PreAcc(:))*(EvolData%dt**2)/6.0

      ! (2) PREDICTED VELOCITY
      NewVel(:) = Vel(:) + (3.*Acc(:)-PreAcc(:))*EvolData%dt/2.0
 
      IF ( .NOT. EvolData%HasThermostat ) THEN        ! Integration without Langevin thermostat

            ! (2) NEW ACCELERATION
            V = GetPotential( NewPos, NewAcc )         ! Compute new forces and store the potential value
            NewAcc(:) =  NewAcc(:) / EvolData%Mass(:)   ! Devide by the mass

      ELSE IF ( ( EvolData%HasThermostat ) ) THEN
            
            ! (3) NEW ACCELERATION
            V = GetPotential( NewPos, NewAcc )         ! Compute new forces and store the potential value

            DO iDoF = 1, 4     ! H and C1 atoms are not Langevin particles
               NewAcc(iDoF) = NewAcc(iDoF)  / EvolData%Mass(iDoF)
            END DO

            IF ( GaussianNoise ) THEN               ! Add gaussian noise and friction
               DO iDoF = 5, EvolData%NDoF
                  NewAcc(iDoF) = ( NewAcc(iDoF) + GaussianRandomNr(EvolData%ThermalNoise(iDoF)) ) / EvolData%Mass(iDoF) &
                                                                                   - EvolData%Gamma*NewVel(iDoF)
               END DO
            ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise and friction
               DO iDoF = 5, EvolData%NDoF
                  NewAcc(iDoF) = ( NewAcc(iDoF) + UniformRandomNr(-sqrt(3.)*EvolData%ThermalNoise(iDoF),sqrt(3.)*EvolData%ThermalNoise(iDoF)) ) &
                                        / EvolData%Mass(iDoF) - EvolData%Gamma*NewVel(iDoF)
               END DO
            END IF   

!             IF ( GaussianNoise ) THEN               ! Add gaussian noise and friction
!                DO iDoF = 1, EvolData%NDoF
!                   NewAcc(iDoF) = ( NewAcc(iDoF) + GaussianRandomNr(EvolData%ThermalNoise(iDoF)) ) / EvolData%Mass(iDoF) &
!                                                                                    - EvolData%Gamma*NewVel(iDoF)
! !                   NewAcc(iDoF) = ( NewAcc(iDoF) + GaussianRandomNr(EvolData%ThermalNoise(iDoF)) - EvolData%Gamma*NewVel(iDoF) ) &
! !                                                         / EvolData%Mass(iDoF)
!                END DO
!             ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise and friction
!                DO iDoF = 1, EvolData%NDoF
!                   NewAcc(iDoF) = ( NewAcc(iDoF) + UniformRandomNr(-sqrt(3.)*EvolData%ThermalNoise(iDoF),sqrt(3.)*EvolData%ThermalNoise(iDoF)) &
!                                               - EvolData%Gamma*NewVel(iDoF) ) / EvolData%Mass(iDoF)
!                END DO
!             END IF   


      END IF

      ! CORRECTED POSITIONS
      Pos(:) = Pos(:) + Vel(:)*EvolData%dt + (NewAcc(:)+2*Acc(:))*(EvolData%dt**2)/6.0
! 
      ! (4) CORRECTED VELOCITIES
      Vel(:) = Vel(:) + (Acc(:) + NewAcc(:))*EvolData%dt/2.0    

      ! Store new acceleration
      PreAcc(:) = Acc(:)
      Acc(:) = NewAcc(:)

   END SUBROUTINE EOM_Beeman   

   SUBROUTINE EOM_ImpulseIntegrator( EvolData, Pos, Vel, Acc, PrePos, GetPotential, V )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc, PrePos
      REAL, INTENT(OUT)                                :: V

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), INTENT(IN)  :: X
            REAL, DIMENSION(:), INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE
      
      INTEGER :: iDoF

      ! Temporary array for predicted velocity and new accelerations
      REAL, DIMENSION( EvolData%NDoF ) :: NewPos
      REAL  :: ExpGtau

      ExpGtau = exp(-EvolData%Gamma*EvolData%dt)

 
      IF ( .NOT. EvolData%HasThermostat ) THEN        ! Integration without Langevin thermostat
          STOP

      ELSE IF ( ( EvolData%HasThermostat ) ) THEN

            ! (3) NEW ACCELERATION
            V = GetPotential( Pos, Acc )         ! Compute new forces and store the potential value

            IF ( GaussianNoise ) THEN               ! Add gaussian noise and friction
               DO iDoF = 1, EvolData%NDoF
                  Acc(iDoF) = ( Acc(iDoF) + GaussianRandomNr(EvolData%ThermalNoise(iDoF)) ) / EvolData%Mass(iDoF) 
               END DO
            ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise and friction
               DO iDoF = 1, EvolData%NDoF
                  Acc(iDoF) = ( Acc(iDoF) + UniformRandomNr(-sqrt(3.)*EvolData%ThermalNoise(iDoF),sqrt(3.)*EvolData%ThermalNoise(iDoF)) ) &
                                    / EvolData%Mass(iDoF) 
               END DO
            END IF   

            Vel(:) = (EvolData%Gamma*ExpGtau/(1-ExpGtau))*(Pos(:)-PrePos(:)) + EvolData%dt* &
                     ( 1.0- (ExpGtau-1.+EvolData%Gamma*EvolData%dt) / (EvolData%Gamma*EvolData%dt*(1-ExpGtau)) ) * Acc(:) 
            NewPos(:) = (1+ExpGtau)*Pos(:) - ExpGtau*PrePos(:) + EvolData%dt/EvolData%Gamma*(1-ExpGtau)*Acc(:)

      END IF
    
      ! Store new acceleration
      PrePos(:) = Pos(:)
      Pos(:) = NewPos(:)

   END SUBROUTINE EOM_ImpulseIntegrator   


END MODULE ClassicalEqMotion