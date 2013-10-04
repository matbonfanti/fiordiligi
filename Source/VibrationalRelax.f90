!***************************************************************************************
!*                              MODULE VibrationalRelax
!***************************************************************************************
!
!>  \brief     Subroutines for the simulations of vibrational relaxation
!>  \details   This module contains subroutines ... \n
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             4 October 2013
!>
!***************************************************************************************

   


   CONTAINS

   SUBROUTINE VibrationalRelax_ReadInput()
      IMPLICIT NONE

      IF (RunType == RELAXATION ) THEN
         ! Initial energy of the system (input in eV and transform to AU)
         CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
         InitEnergy = InitEnergy / MyConsts_Hartree2eV
         ! Nr of initial snapshots
         CALL SetFieldFromInput( InputData, "NrOfInitSnapshots", NrOfInitSnapshots )
         ! Time between each snapshots (input in fs and transform to AU)
         CALL SetFieldFromInput( InputData, "TimeBetweenSnaps", TimeBetweenSnaps)
         TimeBetweenSnaps = TimeBetweenSnaps * MyConsts_fs2AU
         ! Langevin relaxation at the border of the slab
         CALL SetFieldFromInput( InputData, "DynamicsGamma",  DynamicsGamma, 0.0 ) 
         DynamicsGamma = DynamicsGamma / MyConsts_fs2AU

      ELSE IF ( RunType == OSCIBATH_RELAX ) THEN
         ! Initial energy of the system (input in eV and transform to AU)
         CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
         InitEnergy = InitEnergy / MyConsts_Hartree2eV
         ! Nr of initial snapshots
         CALL SetFieldFromInput( InputData, "NrOfInitSnapshots", NrOfInitSnapshots )
         ! Time between each snapshots (input in fs and transform to AU)
         CALL SetFieldFromInput( InputData, "TimeBetweenSnaps", TimeBetweenSnaps)
         TimeBetweenSnaps = TimeBetweenSnaps * MyConsts_fs2AU
         ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
         BathCutOffFreq = 0.0
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, BathCutOffFreq )
         BathCutOffFreq = BathCutOffFreq * MyConsts_cmmin1toAU
         ! Read file with spectral density / normal modes freq and couplings
         CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
         ! No langevin oscillators in the thermal bath
         DynamicsGamma = 0.0

      ELSE IF ( RunType == CHAINBATH_RELAX ) THEN
         ! Initial energy of the system (input in eV and transform to AU)
         CALL SetFieldFromInput( InputData, "InitEnergy", InitEnergy)
         InitEnergy = InitEnergy / MyConsts_Hartree2eV
         ! Nr of initial snapshots
         CALL SetFieldFromInput( InputData, "NrOfInitSnapshots", NrOfInitSnapshots )
         ! Time between each snapshots (input in fs and transform to AU)
         CALL SetFieldFromInput( InputData, "TimeBetweenSnaps", TimeBetweenSnaps)
         TimeBetweenSnaps = TimeBetweenSnaps * MyConsts_fs2AU
         ! Read cutoff frequency of the bath, if BathCutOffFreq is not present, it is set to zero
         BathCutOffFreq = 0.0
         CALL SetFieldFromInput( InputData, "BathCutOffFreq", BathCutOffFreq, BathCutOffFreq )
         BathCutOffFreq = BathCutOffFreq * MyConsts_cmmin1toAU
         ! Read file with spectral density / normal modes freq and couplings
         CALL SetFieldFromInput( InputData, "SpectralDensityFile", SpectralDensityFile )
         ! Langevin relaxation at the end of the chain
         CALL SetFieldFromInput( InputData, "DynamicsGamma",  DynamicsGamma, 0.0 ) 
         DynamicsGamma = DynamicsGamma / MyConsts_fs2AU

      END IF

   END SUBROUTINE VibrationalRelax_ReadInput

   SUBROUTINE VibrationalRelax_ReadInput()
      IMPLICIT NONE
   END SUBROUTINE 

   SUBROUTINE VibrationalRelax()
      IMPLICIT NONE
      REAL, DIMENSION(:), ALLOCATABLE      :: AverageEBath, AverageECoup, AverageESys
      REAL, DIMENSION(:,:), ALLOCATABLE    :: AverageCoord
      REAL, DIMENSION(:), ALLOCATABLE      :: AverCoordOverTime
      REAL, DIMENSION(:,:), ALLOCATABLE    :: OscillCorrelations
      INTEGER                              :: AvEnergyOutputUnit, AvCoordOutputUnit
      INTEGER                              :: AvBathCoordUnit, OscillCorrUnit, InitialCondUnit
      REAL, DIMENSION(NrOfInitSnapshots,8) :: CHInitConditions
      INTEGER                              :: NTimeStepEachSnap
      REAL, DIMENSION(4)                   :: Dummy
      INTEGER                              :: DebugUnitEn, DebugUnitCoord, DebugUnitVel

      REAL :: VSys, KSys, Ecoup, VBath, KBath

      ! ALLOCATE AND INITIALIZE DATA WITH THE AVERAGES OVER THE SET OF TRAJECTORIES

      ALLOCATE( AverageESys(0:ntime), AverageEBath(0:ntime), AverageECoup(0:ntime) )
      ALLOCATE( AverageCoord(size(X),0:ntime) )

      ! allocate arrays for the power spectrum computation
      IF ( PrintType >= FULL ) THEN
         ALLOCATE( PowerSpec(0:ntime), XatT0(size(X)) )
         ALLOCATE( AverCoordOverTime(size(X)) )
         PowerSpec(:) = cmplx( 0.0, 0.0 )
         IF ( RunType == CHAINBATH_RELAX ) THEN
            ALLOCATE( OscillCorrelations(size(X)-4, 0:ntime) )
            OscillCorrelations(:,:) = 0.0
         ENDIF
      ENDIF

      ! Initialize the variables for the trajectory averages
      AverageESys(0:ntime)            = 0.0
      AverageEBath(0:ntime)           = 0.0
      AverageECoup(0:ntime)           = 0.0
      AverageCoord(1:size(X),0:ntime) = 0.0

      ! Open output file to print the brownian realizations vs time
      AvEnergyOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageEnergy.dat", UNIT=AvEnergyOutputUnit )
      WRITE(AvEnergyOutputUnit, "(A,I6,A,/)") "# average energy vs time (fs | eV)"
      AvCoordOutputUnit = LookForFreeUnit()
      OPEN( FILE="AverageCoords.dat", UNIT=AvCoordOutputUnit )
      WRITE(AvCoordOutputUnit, "(A,I6,A,/)") "# average coordinate vs time (fs | Ang)"
      AvBathCoordUnit = LookForFreeUnit()
      OPEN( FILE="AverageBathCoords.dat", UNIT=AvBathCoordUnit )
      WRITE(AvBathCoordUnit, "(A,I6,A,/)") "# average Q coordinate vs time (Ang | fs)"

      IF ( PrintType >= FULL ) THEN
         PowerSpectrumUnit = LookForFreeUnit()
         OPEN( FILE="PowerSpectrum.dat", UNIT=PowerSpectrumUnit )
         WRITE(PowerSpectrumUnit, "(A,I6,A,/)") "# Power spectrum of the trajs - ", inum, " trajectories (fs | au)"
         IF ( RunType == CHAINBATH_RELAX ) THEN
            OscillCorrUnit = LookForFreeUnit()
            OPEN( FILE="OscillCorrelations.dat", UNIT=OscillCorrUnit )
            WRITE(OscillCorrUnit, "(A,I6,A,/)") "# Correlations of oscillators - ", inum, " trajectories (fs | au)"
         END IF
      ENDIF

      ! Initialize random number seed
      CALL SetSeed( 1 )

      PRINT "(2/,A)",    "***************************************************"
      PRINT "(A,F10.5)", "               VIBRATIONAL RELAXATION"
      PRINT "(A,/)" ,    "***************************************************"

      PRINT "(A,I5,A)"," Running ", inum, " trajectories ... "

      IF ( RunType == OSCIBATH_RELAX )  PRINT "(A)"," Bath represented as indipendent oscillators in normal form. "
      IF ( RunType == CHAINBATH_RELAX ) PRINT "(A)"," Bath represented as indipendent oscillators in chain form. "
      IF ( RunType == RELAXATION )  PRINT "(A)"," Bath represented with atomistic force field. "
      IF  ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) &
               PRINT "(A,F15.6)"," Distorsion force constant: ", GetDistorsionForce() 

      ! Define the set of initial conditions for CH:
      ! atoms in the equilibrium geometry, energy in the stretching normal mode

      X(1:2) = 0.0
      V(1:2) = 0.0
      X(3) = HZEquilibrium 
      X(4) = C1Puckering  
      V(3) = + sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / rmh ) * 0.958234410548192
      V(4) = - sqrt( 2.0 * ( InitEnergy-MinimumEnergy ) / rmc ) * 0.285983940880181 

      ! Define when to store the dynamics snapshot
      NTimeStepEachSnap = INT( TimeBetweenSnaps / dt )

      ! Compute starting potential and forces
      A(:) = 0.0
      PotEnergy = VHFourDimensional( X(1:4), A(1:4) )
      DO iCoord = 1,4
         A(iCoord) = A(iCoord) / MassVector(iCoord)
      END DO

      ! Cycle over the nr of snapshot to store
      DO iTraj = 1, NrOfInitSnapshots

         ! Propagate the 4D traj in the microcanonical ensamble
         DO iStep = 1, NTimeStepEachSnap
            ! Propagate for one timestep with Velocity-Verlet
            CALL EOM_VelocityVerlet( InitialConditions, X(1:4), V(1:4), A(1:4), VHFourDimensional, PotEnergy )
            ! compute kinetic energy and total energy
            KinEnergy = EOM_KineticEnergy(InitialConditions, V(1:4) )
            TotEnergy = PotEnergy + KinEnergy
         END DO

         ! Store snapshot
         CHInitConditions( iTraj, 1:4 ) = X(1:4)
         CHInitConditions( iTraj, 5:8 ) = V(1:4)
      END DO

      IF ( PrintType == EQUILIBRDBG ) THEN
         InitialCondUnit = LookForFreeUnit()
         OPEN( FILE="InitialConditionsCH.dat", UNIT=InitialCondUnit )
         WRITE(InitialCondUnit, "(A,I6,A,/)") "# Initial Conditions of the system: X_H, Y_H, Z_H, Z_C and velocities (au)"
         DO iTraj = 1, NrOfInitSnapshots
            WRITE( InitialCondUnit, "(8F15.5)" ) CHInitConditions( iTraj, : )
         END DO
      ENDIF

      !run inum number of trajectories
      DO iTraj = 1,inum
      
         PRINT "(/,A,I6,A)"," **** Trajectory Nr. ", iTraj," ****" 

         !*************************************************************
         ! INITIALIZATION OF THE COORDINATES AND MOMENTA OF THE SYSTEM
         !*************************************************************

         ! Set initial conditions
         IF ( RunType == RELAXATION ) THEN 
               CALL ZeroKelvinSlabConditions( X, V, CHInitConditions ) 
         ELSE IF ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) THEN
               CALL ZeroKelvinBathConditions( X, V, CHInitConditions )
         END IF

         ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
         IF ( RunType == RELAXATION ) THEN 
            X(3:size(X)) = X(3:size(X))  - (X(5)+X(6)+X(7))/3.0
         ENDIF

         !*************************************************************
         ! INFORMATION ON INITIAL CONDITIONS, INITIALIZATION, OTHER...
         !*************************************************************

         ! Compute kinetic energy of system and bath and istantaneous temperature
         KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
         ! Energy of the system
         KSys = EOM_KineticEnergy( InitialConditions, V(1:4) )
         VSys = VHFourDimensional( X(1:4), Dummy )
         ! Coupling energy and energy of the bath
         IF ( RunType == RELAXATION ) THEN 
            CALL WARN( 1>0,  " NOT IMPLEMENTED YET ! " )
         ELSE IF ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) THEN
            Ecoup = CouplingEnergy( X )
            VBath = PotEnergyOfTheBath( X )
            KBath = KinEnergy - KSys
         END IF
         IstTemperature = 2.0*KBath/(MyConsts_K2AU*(size(X)-4))

         ! Store starting values of the averages
         AverageESys(0)            = AverageESys(0)  + KSys + VSys
         AverageEBath(0)           = AverageEBath(0) + KBath + VBath
         AverageECoup(0)           = AverageECoup(0) + Ecoup
         AverageCoord(1:size(X),0) = AverageCoord(1:size(X),0) + X(:)
         IF ( PrintType >= FULL ) THEN
            XatT0(:) = X(:) 
            PowerSpec(0) = PowerSpec(0) + CrossCorrelation( XatT0(1:size(X)), XatT0(1:size(X)) )
            IF ( RunType == CHAINBATH_RELAX ) THEN
               DO iCoord = 1, size(X)-4
                  OscillCorrelations(iCoord, 0) = OscillCorrelations(iCoord, 0) + X(3+iCoord) * X(3+iCoord+1) 
               END DO
            END IF
         ENDIF

         ! PRINT INITIAL CONDITIONS of THE TRAJECTORY
         IF ( PrintType >= FULL ) THEN
            WRITE(*,600)  (KSys+VSys)*MyConsts_Hartree2eV, KBath*MyConsts_Hartree2eV, IstTemperature
            600 FORMAT (/, " Initial condition of the MD trajectory ",/   &
                           " * Energy of the system (eV)        ",1F10.4,/    &
                           " * Kinetic Energy of the bath (eV)  ",1F10.4,/    &
                           " * Istantaneous temperature (K)     ",1F10.4,/ ) 
         END IF

         IF ( PrintType == DEBUG ) THEN
            ! Open unit for massive output, with detailed info on trajectories
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Energy.dat"
            DebugUnitEn = LookForFreeUnit()
            OPEN( Unit=DebugUnitEn, File=OutFileName, Position="APPEND" )
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Coord.dat"
            DebugUnitCoord = LookForFreeUnit()
            OPEN( Unit=DebugUnitCoord, File=OutFileName, Position="APPEND" )
            WRITE(OutFileName,"(A,I4.4,A)") "Traj_",iTraj,"_Vel.dat"
            DebugUnitVel = LookForFreeUnit()
            OPEN( Unit=DebugUnitVel, File=OutFileName, Position="APPEND" )

            ! Write initial values
            WRITE( DebugUnitEn, "(/,A)" ) "# TRAJECTORY ENERGY: time / fs | ESys, ECoup, EBath, KSys, VSys, KBath, VBath / Eh "
            WRITE(DebugUnitEn,800) 0.0,  KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath
            WRITE( DebugUnitCoord, "(/,A)" ) "# TRAJECTORY COORD: time / fs | X(1) X(2) ... X(N) / bohr "
            WRITE(DebugUnitCoord,800) 0.0, X(1), X(2), X(3), X(4), (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, size(X)-4) /)
            WRITE( DebugUnitVel, "(/,A)" ) "# TRAJECTORY VELOCITIES: time / fs | X(1) X(2) ... X(N) / au "
            WRITE(DebugUnitVel,800) 0.0, V(:)
            800 FORMAT(F12.5,1000F15.8)
          ENDIF

         !*************************************************************
         !         TIME EVOLUTION OF THE TRAJECTORY
         !*************************************************************
 
         ! initialize counter for printing steps
         kstep=0

         PRINT "(/,A)", " Propagating the H-Graphene system in time... "
         
         ! Compute starting potential and forces
         A(:) = 0.0
         IF ( RunType == RELAXATION ) THEN 
               PotEnergy = VHSticking( X, A )
         ELSE IF ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) THEN
               PotEnergy = PotentialIndepOscillatorsModel( X, A )
         END IF
         DO iCoord = 1,size(X)
            A(iCoord) = A(iCoord) / MassVector(iCoord)
         END DO

         ! cycle over nstep velocity verlet iterations
         DO iStep = 1,nstep

            ! Propagate for one timestep
            IF ( RunType == RELAXATION ) THEN 
               CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, VHSticking, PotEnergy )
            ELSE IF ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) THEN
               CALL EOM_VelocityVerlet( MolecularDynamics, X, V, A, PotentialIndepOscillatorsModel, PotEnergy )
            END IF

            ! output to write every nprint steps 
            IF ( mod(iStep,nprint) == 0 ) THEN

               ! increment counter for printing steps
               kstep=kstep+1

               ! Total kinetic energy, total energy
               KinEnergy = EOM_KineticEnergy( MolecularDynamics, V )
               TotEnergy = PotEnergy + KinEnergy

               ! for atomistic model of the bath, move the slab so that C2-C3-C4 plane has z=0
               IF ( RunType == RELAXATION ) THEN 
                  X(3:size(X)) = X(3:size(X))  - (X(5)+X(6)+X(7))/3.0
               ENDIF

               ! Energy of the system
               KSys = EOM_KineticEnergy( InitialConditions, V(1:4) )
               VSys = VHFourDimensional( X(1:4), Dummy )
               ! Coupling energy and energy of the bath
               IF ( RunType == RELAXATION ) THEN 
                  CALL WARN( 1>0,  " NOT IMPLEMENTED YET ! " )
               ELSE IF ( RunType == OSCIBATH_RELAX .OR. RunType == CHAINBATH_RELAX ) THEN
                  Ecoup = CouplingEnergy( X )
                  VBath = PotEnergyOfTheBath( X )
                  KBath = KinEnergy - KSys
               END IF
               
               ! Update averages over time
               AverageESys(kstep)            = AverageESys(kstep)  + KSys + VSys
               AverageEBath(kstep)           = AverageEBath(kstep) + KBath + VBath
               AverageECoup(kstep)           = AverageECoup(kstep) + Ecoup
               AverageCoord(1:size(X),kstep) = AverageCoord(1:size(X),kstep) + X(:)
               ! store the autocorrelation function for power spectrum computation
               IF ( PrintType >= FULL )  THEN
                  PowerSpec(kstep) = PowerSpec(kstep) + CrossCorrelation( XatT0(1:size(X)), X(1:size(X)) )
                  IF ( RunType == CHAINBATH_RELAX ) THEN
                     DO iCoord = 1, size(X)-4
                        OscillCorrelations(iCoord, kstep) = OscillCorrelations(iCoord, kstep) + X(3+iCoord) * X(3+iCoord+1) 
                     END DO
                  END IF
               END IF

               ! If massive level of output, print traj information to std out
               IF ( PrintType == DEBUG ) THEN
                  WRITE(DebugUnitEn,800) dt*real(iStep)/MyConsts_fs2AU, KSys+VSys, Ecoup, KBath+VBath, KSys, VSys, KBath, VBath
                  WRITE(DebugUnitCoord,800) dt*real(iStep)/MyConsts_fs2AU, X(1:4), &
                                                        (/ (X(iCoord+4)+0.05*iCoord, iCoord = 1, size(X)-4) /)
                  WRITE(DebugUnitVel,800) dt*real(iStep)/MyConsts_fs2AU, V(:)
               END IF

            END IF 

         END DO

         PRINT "(A)", " Time propagation completed! "

         IF ( PrintType == DEBUG ) THEN
               CLOSE( Unit=DebugUnitEn )
               CLOSE( Unit=DebugUnitCoord )
               CLOSE( Unit=DebugUnitVel )
         ENDIF

      END DO
      
      PRINT "(A)"," Done! "

      !*************************************************************
      !         OUTPUT OF THE RELEVANT AVERAGES 
      !*************************************************************

      ! Normalize averages 
      AverageESys(:)            = AverageESys(:)            / real(inum)
      AverageEBath(:)           = AverageEBath(:)           / real(inum)
      AverageECoup(:)           = AverageECoup(:)           / real(inum)
      AverageCoord(1:size(X),:) = AverageCoord(1:size(X),:) / real(inum) 
      OscillCorrelations(1:size(X)-4,:) = OscillCorrelations(1:size(X)-4,:)  / real(inum)

      ! PRINT average energy of the system, of the coupling, of the bath
      DO iStep = 0, ntime
         WRITE(AvEnergyOutputUnit,"(F14.8,3F14.8)") dt*real(nprint*iStep)/MyConsts_fs2AU, AverageESys(iStep)*MyConsts_Hartree2eV, &
                 AverageECoup(iStep)*MyConsts_Hartree2eV,  AverageEBath(iStep)*MyConsts_Hartree2eV
      END DO

      ! PRINT average coordinates
      DO iStep = 0, ntime
         WRITE(AvCoordOutputUnit,"(F14.8,3F14.8)") dt*real(nprint*iStep)/MyConsts_fs2AU, AverageCoord(3,iStep)*MyConsts_Bohr2Ang, &
                    AverageCoord(4,iStep)*MyConsts_Bohr2Ang, FirstEffectiveMode(AverageCoord(5:size(X),iStep))*MyConsts_Bohr2Ang
      END DO

      ! PRINT average bath coordinates
      DO iCoord = 1, nevo-1
         WRITE(AvBathCoordUnit,"(/,A,I5,/)") "#  Bath Coord # ", iCoord
         DO iStep = 0, ntime
            WRITE(AvBathCoordUnit,"(F20.8,3F20.8)") real(iCoord)+AverageCoord(4+iCoord,iStep)*10, &
                      dt*real(nprint*iStep)/MyConsts_fs2AU
                    
         END DO
      END DO

      ! normalize autocorrelation function, compute fourier transfrom and print
      IF ( PrintType >= FULL ) THEN
            AverCoordOverTime(:) = SUM ( AverageCoord, 2 ) / real(ntime+1)

            PowerSpec(:) = PowerSpec(:) /  real(inum )
            PowerSpec(:) = PowerSpec(:) - DOT_PRODUCT(AverCoordOverTime(1:size(X)), AverCoordOverTime(1:size(X)) )
            dOmega =  2.*MyConsts_PI/( real(nstep) * dt )

            CALL DiscreteInverseFourier( PowerSpec )
            DO iOmega = 0, ntime
               WRITE( PowerSpectrumUnit,* )  iOmega*dOmega/MyConsts_cmmin1toAU,  real(PowerSpec(iOmega)), aimag( PowerSpec(iOmega))
            END DO
            CLOSE( PowerSpectrumUnit )

            IF ( RunType == CHAINBATH_RELAX ) THEN
               DO iCoord = 1, size(X)-4
                  OscillCorrelations(iCoord,:) = OscillCorrelations(iCoord,:) - &
                                          AverCoordOverTime(3+1)*AverCoordOverTime(3+iCoord+1) 
                  WRITE(OscillCorrUnit, *) iCoord, SUM( OscillCorrelations(iCoord, :) )/real(ntime+1)
               END DO
            END IF

            DEALLOCATE( PowerSpec, AverCoordOverTime )
      ENDIF

      ! Close output files
      CLOSE( AvEnergyOutputUnit )
      CLOSE( AvCoordOutputUnit )


   END SUBROUTINE HGraphiteVibrationalRelax
