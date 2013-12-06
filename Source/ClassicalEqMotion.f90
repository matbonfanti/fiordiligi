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
!
!>  \par Updates
!>  \arg 8 Novembre 2013: thermoswitch is now an optional argument of the
!>                        SetupThermostat subroutine (default: TRUE for all entries)
!>  \arg 28 November 2013: ring polymer propagation implemented, with symplectic
!>                         integrator
!
!>  \todo   clean up the code: leave only 1 propagator for RPMD and 1 propagator for
!>          normal MD, with or without langevin friction ( in case of RPMD, Parrinello
!>          algorithm, in case of regular MD, stochastic integrated algorith - both have
!>          Velocity Verlet as limit for gamma = 0 ) then rename subroutines for propagation
!>          and clean the datatype with the only data that are necessary finally fix 
!>          the setup subroutines including all the necessary checks and setup options
!>                 
!***************************************************************************************

MODULE ClassicalEqMotion
#include "preprocessoptions.cpp"
   USE RandomNumberGenerator
   USE FFTWrapper
   IMPLICIT NONE
   
      PRIVATE
      PUBLIC :: Evolution
      PUBLIC :: EvolutionSetup, SetupThermostat, SetupRingPolymer
      PUBLIC :: DisposeEvolutionData, DisposeThermostat, DisposeRingPolymer
      PUBLIC :: EOM_KineticEnergy
      PUBLIC :: EOM_VelocityVerlet, EOM_Beeman, EOM_LangevinSecondOrder, EOM_RPMSymplectic

      REAL, PARAMETER :: Over2Sqrt3 = 1.0 / ( 2.0 * SQRT(3.0) )

      LOGICAL :: GaussianNoise = .TRUE.
      
      !> Evolution datatype, storing all the general data required
      !> for integration of the EoM, with or without Langevin thermostat and with or without Ring Polymer MD
      TYPE Evolution
         INTEGER   :: NDoF                               !< Nr of degrees of freedom
         REAL, DIMENSION(:), POINTER :: Mass             !< Vector with the masses of the system
         REAL  :: dt                                     !< Time step of integration
         REAL  :: Gamma                                  !< Integrated Langevin friction for one timestep
         REAL  :: FrictionCoeff_HalfDt                   !< Coefficient including Langevin friction for half timestep
         REAL, DIMENSION(:), POINTER :: ThermalNoise     !< Vector of the thermal noise sigma
         REAL, DIMENSION(:), POINTER :: ThermalNoise2     !< Vector of the thermal noise sigma
         LOGICAL, DIMENSION(:), POINTER :: ThermoSwitch  !< Vector to set the thermostat on or off for the dof

         INTEGER :: NBeads                                 !< Nr of system replicas 
         REAL  :: RPFreq                                   !< Frequency of the harmonic force between the beads
         REAL, DIMENSION(:), POINTER :: NormModesFreq      !< Vector with the frequencies of the normal modes 
         REAL, DIMENSION(:,:), POINTER :: NormModesPropag  !< Propagation coefficients of normal coordinates and velocities
         TYPE(FFTHalfComplexType) :: RingNormalModes       !< FFT transform to compute normal modes of the RP 
         REAL, DIMENSION(:), POINTER :: GammaLang          !< Friction coeffs of PILE
         REAL, DIMENSION(:), POINTER :: AlphaLang          !< Set of coeffs for PILE integration
         REAL, DIMENSION(:), POINTER :: BetaLang           !< Set of coeffs for PILE integration

         LOGICAL :: HasThermostat = .FALSE.              !< Thermostat data has been setup
         LOGICAL :: HasRingPolymer = .FALSE.             !< Ring polymer data has been setup
         LOGICAL :: IsSetup = .FALSE.                    !< Evolution data has been setup
      END TYPE Evolution

   CONTAINS
   

!================================================================================================================================
!                                  SETUP SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                        EvolutionSetup
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
      REAL, DIMENSION(:), INTENT(IN)    :: MassVector
      REAL, INTENT(IN)                  :: TimeStep

      ! Check dimension of mass vector
      CALL ERROR( size(MassVector) /= NDoF, "ClassicalEqMotion.EvolutionSetup: MassVector array dimension mismatch" )

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
      EvolutionData%HasRingPolymer = .FALSE.
      EvolutionData%HasThermostat = .FALSE.
      
#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,I5,A,1F8.3)") "Evolution data type is setup: Nr DoF is ",NDoF," and TimeStep ", TimeStep
#endif
      
   END SUBROUTINE EvolutionSetup

   
!*******************************************************************************
!                       SetupThermostat
!*******************************************************************************
!> Setup thermostat data, store them in the Evolution datatype.
!>
!> @param EvolData     Evolution data type to setup
!> @param Gamma        Langevin friction coefficient
!> @param Temperature  Temperature of thermostat
!*******************************************************************************
   SUBROUTINE SetupThermostat( EvolData, Gamma, Temperature, ThermoSwitch )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)            :: EvolData
      REAL, INTENT(IN)                            :: Gamma, Temperature
      LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: ThermoSwitch
      INTEGER :: iDoF
      
      ! error if trying to setup thermostat of a non-setup evolution type
      CALL ERROR( .NOT. EvolData%IsSetup, "ClassicalEqMotion.SetupThermostat: evolution data not setup" )
      ! warn user if overwriting previously setup data
      CALL WARN( EvolData%HasThermostat, "ClassicalEqMotion.SetupThermostat: overwriting thermostat data" ) 
      ! error if the thermostat switch has wrong dimension
      IF ( PRESENT(ThermoSwitch) )  CALL ERROR( size(ThermoSwitch) /= EvolData%NDoF , &
                                          "ClassicalEqMotion.SetupThermostat: thermostat switch mismatch" )

      ! Store gamma value
      EvolData%Gamma = Gamma
      ! Set coefficient for velocity integration (in Vel-Verlet) with Langevin friction
      EvolData%FrictionCoeff_HalfDt = 1.0 - 0.5 * Gamma * EvolData%dt
      
      ! Store the Thermostat switch array
      ALLOCATE( EvolData%ThermoSwitch(EvolData%NDoF) )
      IF ( PRESENT(ThermoSwitch) ) THEN
         EvolData%ThermoSwitch(:) = ThermoSwitch(:)
      ELSE
         EvolData%ThermoSwitch(:) = .TRUE.
      END IF

      ! Set standard deviations of the thermal noise
      IF ( .NOT. EvolData%HasThermostat ) &
                 ALLOCATE( EvolData%ThermalNoise( EvolData%NDoF ), EvolData%ThermalNoise2( EvolData%NDoF ) )
      EvolData%ThermalNoise(:) = 0.0
      EvolData%ThermalNoise2(:) = 0.0
      DO iDoF = 1, EvolData%NDoF
         IF ( EvolData%ThermoSwitch(iDoF) ) THEN
            EvolData%ThermalNoise(iDoF) = sqrt( 2.0*Temperature*EvolData%Mass(iDoF)*Gamma/EvolData%dt )
            EvolData%ThermalNoise2(iDoF) = sqrt( 2.0*Temperature*Gamma/EvolData%Mass(iDoF) )
         END IF
      END DO

      ! Set arrays for PILE integration when setting thermostat for a ring polymer propagator
      IF ( EvolData%HasRingPolymer ) THEN
         IF ( .NOT. EvolData%HasThermostat )     ALLOCATE( EvolData%GammaLang(EvolData%NBeads), &
                                                EvolData%AlphaLang(EvolData%NBeads), EvolData%BetaLang(EvolData%NBeads) )

         ! Set friction parameters for the ring normal modes
         EvolData%GammaLang(1) = Gamma
         DO iDoF = 2, EvolData%NBeads
            EvolData%GammaLang(iDoF) = 2.0 * EvolData%NormModesFreq(iDoF)
         END DO

         ! Set integration coefficients
         DO iDoF = 1, EvolData%NBeads
            EvolData%AlphaLang(iDoF) = EXP( - 0.5 * EvolData%dt * EvolData%GammaLang(iDoF) )
            EvolData%BetaLang(iDoF)  = SQRT( (1.0 - EvolData%AlphaLang(iDoF)**2) * Temperature * EvolData%NBeads )
         END DO
      END IF

      ! Themostat data is now setup
      EvolData%HasThermostat = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,1F8.3,A,1F8.3)") "Thermostat is setup with Gamma = ",Gamma," and Temperature = ", Temperature
      WRITE(*,*) " Langevin DoFs: ", EvolData%ThermoSwitch(:)
#endif
      
   END SUBROUTINE SetupThermostat
   

!*******************************************************************************
!                          SetupRingPolymer
!*******************************************************************************
!> Setup data for ring polymer progation, store them in the Evolution datatype.
!>
!> @param EvolData     Evolution data type to setup
!> @param NBeads       Nr of replicas of the system for the RPMD
!*******************************************************************************
   SUBROUTINE SetupRingPolymer( EvolData, NBeads, PolymerFreq )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)    :: EvolData
      INTEGER, INTENT(IN)                 :: NBeads
      REAL, INTENT(IN)                    :: PolymerFreq
      INTEGER :: i

      ! error if trying to setup thermostat of a non-setup evolution type
      CALL ERROR( .NOT. EvolData%IsSetup, "ClassicalEqMotion.SetupRingPolymer: evolution data not setup" )
      ! warn user if overwriting previously setup data
      CALL WARN( EvolData%HasRingPolymer, "ClassicalEqMotion.SetupRingPolymer: overwriting ring polymer data" ) 

      ! Store nr of beads and ring polymer harmonic frequency
      EvolData%NBeads = NBeads
      EvolData%RPFreq = PolymerFreq

      ! Setup FFT to compute normal modes of the free ring polymer
      CALL SetupFFT( EvolData%RingNormalModes, EvolData%NBeads ) 

      ! Set frequencies of the normal modes
      ALLOCATE( EvolData%NormModesFreq( EvolData%NBeads ) )
      DO i = 1, EvolData%NBeads
         EvolData%NormModesFreq(i) = 2.0 * EvolData%RPFreq * SIN( MyConsts_PI * real(i-1) / real(EvolData%NBeads) )
      END DO

      ! Set free propagator of the normal modes
      ALLOCATE( EvolData%NormModesPropag( 4,EvolData%NBeads ) )
      EvolData%NormModesPropag(:,1) = (/ 1.0, EvolData%dt, 0.0, 1.0 /)  ! Free evolution of the w=0 normal mode 
      DO i = 2, EvolData%NBeads                                         ! Harmonic oscillator evolution of other modes
         EvolData%NormModesPropag(1,i) = COS( EvolData%NormModesFreq(i)*EvolData%dt )
         EvolData%NormModesPropag(2,i) = SIN( EvolData%NormModesFreq(i)*EvolData%dt ) / EvolData%NormModesFreq(i)
         EvolData%NormModesPropag(3,i) = - SIN( EvolData%NormModesFreq(i)*EvolData%dt ) * EvolData%NormModesFreq(i)
         EvolData%NormModesPropag(4,i) = COS( EvolData%NormModesFreq(i)*EvolData%dt )
      END DO

      ! RPMD data is now setup
      EvolData%HasRingPolymer = .TRUE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A,1I8,A)") "RPMD is setup with ",EvolData%NBeads," replicas "
#endif
      
   END SUBROUTINE SetupRingPolymer



!================================================================================================================================
!                              DISPOSE SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                          DisposeEvolutionData
!*******************************************************************************
!> Dispose evolution data.
!>
!> @param EvolData     Evolution data type  to dispose
!*******************************************************************************
   SUBROUTINE DisposeEvolutionData( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF

      ! continue if trying to dispose data that is not setup
      IF (.NOT. EvolData%IsSetup)  RETURN

      ! dispose thermostat and ringpolymer if it is the case
      IF ( EvolData%HasThermostat )  CALL DisposeThermostat( EvolData )
      IF ( EvolData%HasRingPolymer ) CALL DisposeRingPolymer( EvolData )

      ! Deallocate standard deviations of the thermal noise
      DEALLOCATE( EvolData%Mass )

      ! Themostat data is now disposed
      EvolData%HasThermostat = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Themostat has been disposed"
#endif
   END SUBROUTINE DisposeEvolutionData


!*******************************************************************************
!                          DisposeThermostat
!*******************************************************************************
!> Dispose thermostat data.
!>
!> @param EvolData     Evolution data type with the thermostat to dispose
!*******************************************************************************
   SUBROUTINE DisposeThermostat( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF
      
      ! continue if trying to dispose thermostat that is not setup
      IF ((.NOT. EvolData%IsSetup) .OR. (.NOT. EvolData%HasThermostat))  RETURN

      ! Set coefficient with zero friction
      EvolData%FrictionCoeff_HalfDt = 1.0
      EvolData%Gamma = 0.0
      
      ! Deallocate standard deviations of the thermal noise
      DEALLOCATE( EvolData%ThermalNoise, EvolData%ThermoSwitch )
      
      ! Themostat data is now disposed
      EvolData%HasThermostat = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Themostat has been disposed"
#endif
   END SUBROUTINE DisposeThermostat


!*******************************************************************************
!                          DisposeRingPolymer
!*******************************************************************************
!> Dispose ring polymer data.
!>
!> @param EvolData     Evolution data type with the thermostat to dispose
!*******************************************************************************
   SUBROUTINE DisposeRingPolymer( EvolData )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)  :: EvolData
      INTEGER :: iDoF
      
      ! continue if trying to dispose thermostat that is not setup
      IF ((.NOT. EvolData%IsSetup) .OR. (.NOT. EvolData%HasRingPolymer))  RETURN

      ! Dispose data for FFT
      CALL DisposeFFT( EvolData%RingNormalModes )
   
      ! Ring polymer data is now disposed
      EvolData%HasRingPolymer = .FALSE.

#if defined(VERBOSE_OUTPUT)
      WRITE(*,"(/,A)") "Ring polymer has been disposed"
#endif
   END SUBROUTINE DisposeRingPolymer


!================================================================================================================================
!                         PROPAGATION SUBROUTINES
!================================================================================================================================


!*******************************************************************************
!                          EOM_VelocityVerlet
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
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
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
         DO iDoF = 1, EvolData%NDoF
            IF ( EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = EvolData%FrictionCoeff_HalfDt*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            END IF
         END DO

         ! (3) NEW FORCES AND ACCELERATIONS 
         V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
         IF ( GaussianNoise ) THEN    ! add gaussian noise
            DO iDoF = 1, EvolData%NDoF
               IF ( EvolData%ThermoSwitch(iDoF) ) THEN
                  Acc(iDoF) = ( Acc(iDoF)+GaussianRandomNr( EvolData%ThermalNoise(iDoF) ) ) / EvolData%Mass(iDoF)
               ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
                  Acc(iDoF) = Acc(iDoF)  / EvolData%Mass(iDoF)
               END IF
            END DO
         ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise
            DO iDoF = 1, EvolData%NDoF
               IF ( EvolData%ThermoSwitch(iDoF) ) THEN
                  Acc(iDoF) = ( Acc(iDoF) + UniformRandomNr( -sqrt(3.)*EvolData%ThermalNoise(iDoF),     &
                                                         sqrt(3.)*EvolData%ThermalNoise(iDoF) ) )  / EvolData%Mass(iDoF)
               ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
                  Acc(iDoF) = Acc(iDoF)  / EvolData%Mass(iDoF)
               END IF
            END DO
         END IF

         ! (4) HALF TIME STEP AGAIN FOR THE VELOCITIES
         DO iDoF = 1, EvolData%NDoF
            IF ( EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = EvolData%FrictionCoeff_HalfDt*Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
               Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
            END IF
         END DO

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
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
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

         IF ( GaussianNoise ) THEN    ! add gaussian noise
            DO iDoF = 1, EvolData%NDoF
               IF ( EvolData%ThermoSwitch(iDoF) ) THEN
                  NewAcc(iDoF) = ( NewAcc(iDoF) + GaussianRandomNr(EvolData%ThermalNoise(iDoF)) ) / EvolData%Mass(iDoF) &
                                                                                   - EvolData%Gamma*NewVel(iDoF)
               ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
                  NewAcc(iDoF) = NewAcc(iDoF)  / EvolData%Mass(iDoF)
               END IF
            END DO
         ELSE IF ( .NOT. GaussianNoise ) THEN    ! add uniform noise
            DO iDoF = 1, EvolData%NDoF
               IF ( EvolData%ThermoSwitch(iDoF) ) THEN
                  NewAcc(iDoF) = ( NewAcc(iDoF) + UniformRandomNr(-sqrt(3.)*EvolData%ThermalNoise(iDoF),   &
                             sqrt(3.)*EvolData%ThermalNoise(iDoF)) ) / EvolData%Mass(iDoF) - EvolData%Gamma*NewVel(iDoF)
               ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
                  NewAcc(iDoF) = NewAcc(iDoF)  / EvolData%Mass(iDoF)
               END IF
            END DO
         END IF

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


!*******************************************************************************
!> Propagate trajectory with Vanden-Eijnden and Ciccoti algorith.
!>
!> @param EvolData     Evolution data type
!*******************************************************************************
   SUBROUTINE EOM_LangevinSecondOrder( EvolData, Pos, Vel, Acc, GetPotential, V )
      IMPLICIT NONE

      TYPE( Evolution ), INTENT(INOUT)                 :: EvolData
      REAL, DIMENSION( EvolData%NDoF ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                :: V

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE

      INTEGER :: iDoF
      REAL, DIMENSION( EvolData%NDoF ) :: A
      REAL, DIMENSION( EvolData%NDoF ) :: Xi, Eta 

      ! (0) COMPUTE NECESSARY RANDOM VALUES
      DO iDoF = 1, EvolData%NDoF
         IF ( EvolData%ThermoSwitch(iDoF) ) THEN
            Xi(iDoF)  = GaussianRandomNr(1.0)
            Eta(iDoF) = GaussianRandomNr(1.0)
            A(iDoF) = 0.5 * EvolData%dt**2 * ( Acc(iDoF) - EvolData%Gamma*Vel(iDoF) ) + &
                     EvolData%ThermalNoise2(iDoF) * EvolData%dt**(1.5) * ( 0.5 * Xi(iDoF) + Over2Sqrt3 * Eta(iDoF)  )
         ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
            A(iDoF) = 0.5 * EvolData%dt**2 * Acc(iDoF)
         END IF
      END DO

      ! (1) PARTIAL UPDATE OF THE VELOCITIES
      DO iDoF = 1, EvolData%NDoF
         IF ( EvolData%ThermoSwitch(iDoF) ) THEN
            Vel(iDoF) = (1.0 - EvolData%Gamma*EvolData%dt) * Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
            Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END IF
      END DO

      ! (2) UPDATE POSITION
      Pos(:) = Pos(:) + Vel(:)*EvolData%dt + A(:)

      ! (3) NEW FORCES AND ACCELERATIONS 
      V = GetPotential( Pos, Acc )       ! Compute new forces and store the potential value
      Acc(:) = Acc(:)  / EvolData%Mass(:)

      ! (4) FINAL UPDATE OF THE VELOCITIES
      DO iDoF = 1, EvolData%NDoF
         IF ( EvolData%ThermoSwitch(iDoF) ) THEN
            Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt + SQRT(EvolData%dt) * EvolData%ThermalNoise2(iDoF) * Xi(iDoF) &
                                  - EvolData%Gamma * A(iDoF)
         ELSE IF ( .NOT. EvolData%ThermoSwitch(iDoF) ) THEN
            Vel(iDoF) = Vel(iDoF) + 0.5*Acc(iDoF)*EvolData%dt
         END IF
      END DO

   END SUBROUTINE EOM_LangevinSecondOrder   


!*******************************************************************************
!> Propagate trajectory with symplectic algorithm for Ring-Polymer MD. 
!>
!> @param EvolData     Evolution data type
!*******************************************************************************
   SUBROUTINE EOM_RPMSymplectic( EvolData, Pos, Vel, Acc, GetPotential, V, InitializeAcceleration )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)                                   :: EvolData
      REAL, DIMENSION( EvolData%NDoF * EvolData%NBeads ), INTENT(INOUT)  :: Pos, Vel, Acc
      REAL, INTENT(OUT)                                                  :: V
      LOGICAL, OPTIONAL                                                  :: InitializeAcceleration

      INTERFACE
         REAL FUNCTION GetPotential( X, Force )
            REAL, DIMENSION(:), TARGET, INTENT(IN)  :: X
            REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Force
         END FUNCTION GetPotential
      END INTERFACE

      REAL, DIMENSION(EvolData%NBeads)  :: BeadQAt0, BeadVAt0, BeadQAtT, BeadVAtT
      INTEGER :: iDoF, iBead, iStart, iEnd
      REAL    :: VBead
      LOGICAL :: DoPropagation

      CALL ERROR( .NOT. EvolData%HasRingPolymer, " EOM_RPMSymplectic: data for RPMD are needed "  ) 

      ! the subroutine can be used for initial definition of acceleration
      IF ( PRESENT(InitializeAcceleration) ) THEN
         DoPropagation = .NOT. InitializeAcceleration      ! in such case, no propagation is performed
      ELSE
         DoPropagation = .TRUE.
      END IF

      IF ( DoPropagation ) THEN
      ! (1)  HALF TIME STEP FOR THERMOSTATTING THE SYSTEM (only when langevin dynamics)
         IF ( EvolData%HasThermostat ) THEN
            DO iDoF = 1, EvolData%NDoF

               DO iBead = 1, EvolData%NBeads                        ! extract single bead positions and velocities
                  BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
               END DO
               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT ) ! transform to normal modes coords

               DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
                  BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                    GaussianRandomNr( 1.0 ) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
               END DO

               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT ) ! transform back to original coords
               DO iBead = 1, EvolData%NBeads               ! copy single bead positions and velocities to input arrays
                  Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
               END DO

            END DO
         END IF

      ! (2) HALF TIME STEP FOR THE VELOCITIES
         DO iBead = 1, EvolData%NBeads
            iStart = (iBead-1) * EvolData%NDoF + 1
            iEnd   = iBead * EvolData%NDoF
            Vel( iStart:iEnd ) = Vel( iStart:iEnd ) + 0.5*Acc( iStart:iEnd )*EvolData%dt
         END DO
      END IF

      ! (3) EXACT PROPAGATION WITH THE INTERBEADS POTENTIAL
      DO iDoF = 1, EvolData%NDoF

         DO iBead = 1, EvolData%NBeads      ! extract single bead positions and velocities
            BeadQAt0(iBead) = Pos( (iBead-1) * EvolData%NDoF + iDoF )
            BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
         END DO

         CALL ExecuteFFT( EvolData%RingNormalModes, BeadQAt0, DIRECT_FFT ) ! transform to normal modes coords
         CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT )

         IF ( DoPropagation ) THEN
            DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
               BeadQAtT(iBead) = EvolData%NormModesPropag(1,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(2,iBead) * BeadVAt0(iBead)
               BeadVAtT(iBead) = EvolData%NormModesPropag(3,iBead) * BeadQAt0(iBead) + &
                                                                          EvolData%NormModesPropag(4,iBead) * BeadVAt0(iBead)
            END DO
         ELSE
            BeadQAtT(:) = BeadQAt0(:)
            BeadVAtT(:) = BeadVAt0(:)
         END IF

         V = 0           ! Compute potential energy of the normal modes
         DO iBead = 2, EvolData%NBeads
            V = V + 0.5 * EvolData%Mass(iDoF) * EvolData%NormModesFreq(iBead)**2 * BeadQAtT(iBead)**2
         END DO

         CALL ExecuteFFT( EvolData%RingNormalModes, BeadQAtT, INVERSE_FFT ) ! transform back to original coords
         CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT )

         DO iBead = 1, EvolData%NBeads      ! copy single bead positions and velocities to input arrays
            Pos( (iBead-1) * EvolData%NDoF + iDoF ) = BeadQAtT(iBead)
            Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
         END DO

      END DO

      ! (4) UPDATE FORCES
      V = 0.0
      DO iBead = 1, EvolData%NBeads
         iStart = (iBead-1) * EvolData%NDoF + 1
         iEnd   = iBead * EvolData%NDoF
         VBead = GetPotential( Pos( iStart:iEnd ), Acc( iStart:iEnd ) )   ! Compute new forces and store the potential value
         Acc( iStart:iEnd ) = Acc( iStart:iEnd ) / EvolData%Mass(:)       ! only potential forces 
         V = V + VBead
      END DO

      IF ( DoPropagation ) THEN
      ! (5) HALF TIME STEP FOR THE VELOCITIES
         DO iBead = 1, EvolData%NBeads
            iStart = (iBead-1) * EvolData%NDoF + 1
            iEnd   = iBead * EvolData%NDoF
            Vel( iStart:iEnd ) = Vel( iStart:iEnd ) + 0.5*Acc( iStart:iEnd )*EvolData%dt
         END DO

      ! (6)  HALF TIME STEP FOR THERMOSTATTING THE SYSTEM (only when langevin dynamics)
         IF ( EvolData%HasThermostat ) THEN
            DO iDoF = 1, EvolData%NDoF

               DO iBead = 1, EvolData%NBeads                        ! extract single bead positions and velocities
                  BeadVAt0(iBead) = Vel( (iBead-1) * EvolData%NDoF + iDoF )
               END DO
               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAt0, DIRECT_FFT ) ! transform to normal modes coords

               DO iBead = 1, EvolData%NBeads                                   ! evolve normal modes
                  BeadVAtT(iBead) = BeadVAt0(iBead) * EvolData%AlphaLang(iBead) + &
                                    GaussianRandomNr( 1.0 ) * EvolData%BetaLang(iBead) / SQRT( EvolData%Mass(iDoF) )
               END DO

               CALL ExecuteFFT( EvolData%RingNormalModes, BeadVAtT, INVERSE_FFT ) ! transform back to original coords
               DO iBead = 1, EvolData%NBeads               ! copy single bead positions and velocities to input arrays
                  Vel( (iBead-1) * EvolData%NDoF + iDoF ) = BeadVAtT(iBead)
               END DO

            END DO
         END IF

      END IF

   END SUBROUTINE EOM_RPMSymplectic


!================================================================================================================================
!                              OTHER SUBROUTINES
!================================================================================================================================
   
!*******************************************************************************
!> Compute kinetic energy corresponding to a given velocity vector.
!>
!> @param EvolData     Evolution data type
!> @param Velocity     Array containing the velocity at given time step
!*******************************************************************************
   REAL FUNCTION EOM_KineticEnergy( EvolData, Vel, NMax ) RESULT( KinEnergy )
      IMPLICIT NONE
      TYPE( Evolution ), INTENT(INOUT)    :: EvolData
      REAL, DIMENSION(:), INTENT(INOUT)   :: Vel
      INTEGER, INTENT(IN), OPTIONAL       :: NMax
      INTEGER :: iDoF, iBead, N, NDoF
      
      IF (.NOT. PRESENT( NMax )) THEN
         NDoF = EvolData%NDoF
      ELSE 
         NDoF = MIN(NMax,EvolData%NDoF)
      END IF

      KinEnergy = 0.0
      IF ( EvolData%HasRingPolymer ) THEN
         N = 0
         DO iBead = 1, EvolData%NBeads
            DO iDoF = 1, NDoF
               N = N + 1
               KinEnergy = KinEnergy + 0.5 * EvolData%Mass(iDoF) * Vel(N)**2
            END DO
         END DO

      ELSE IF ( .NOT. EvolData%HasRingPolymer ) THEN
         DO iDoF = 1, NDoF
            KinEnergy = KinEnergy + 0.5 * EvolData%Mass(iDoF) * Vel(iDoF)**2
         END DO

      ENDIF

   END FUNCTION EOM_KineticEnergy


!================================================================================================================================
!                              END OF MODULE
!================================================================================================================================

END MODULE ClassicalEqMotion
