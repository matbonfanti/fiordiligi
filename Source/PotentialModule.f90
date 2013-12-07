MODULE PotentialModule
   USE MyConsts
   USE ErrorTrap
   USE RandomNumberGenerator

   PRIVATE
   PUBLIC :: SetupPotential, VHSticking, VHFourDimensional, MinimizePotential
   PUBLIC :: ThermalEquilibriumConditions, ScatteringConditions, ZeroKelvinSlabConditions
   PUBLIC :: CarbonForceConstant, GraphiteLatticeConstant

   !> Setup variable for the potential
   LOGICAL :: PotentialModuleIsSetup = .FALSE.

   !> Fix the collinear geometry
   LOGICAL :: CollinearPES = .FALSE.

   !> C1 Puckering displacement of the H-graphite potential (in bohr)
   REAL, PUBLIC :: C1Puckering =  0.3673 / MyConsts_Bohr2Ang
   !> H Z equilibrium coordinate (in bohr)
   REAL, PUBLIC :: HZEquilibrium = 1.4811 / MyConsts_Bohr2Ang
   !> Energy of the minimum
   REAL, PUBLIC :: MinimumEnergy = 0.0000
   !> H Z vibration frequency at the minimum
   REAL, PUBLIC :: HZForceConst = 0.0000

   !> Max nr of iterations for potential optimization
   INTEGER, PARAMETER :: MaxIter = 10000
   !> Threshold for conjugate gradient convergence
   REAL, PARAMETER :: GradEps = 1.0E-4
   !> Parameter for finite difference computation
   REAL, PARAMETER :: Delta = 1.0E-4

!> \name STICKING POTENTIAL
!> Parameters of the 4D Potential for C-H
!> @{
   ! parameters for the asymptotic Morse form ( rho -> infty ) of the PES
   REAL, PARAMETER :: di = 0.00775           !< potential well (eV)
   REAL, PARAMETER :: alphi = 0.954          !< curvature (Ang-1)
   REAL, PARAMETER :: ai = 4.01              !< equilibrium position (Ang)
   ! parameters to describe the behaviour along rho at large rho
   REAL, PARAMETER :: alp = 1.4              !< a_f in the text (Ang-1) 
   REAL, PARAMETER :: alp2 = 1.65            !< b_f in the text (Ang-1)
   REAL, PARAMETER :: ba = 2.9               !< c_f in the text (Ang-1)
   REAL, PARAMETER :: rhoa = 1.0             !< rho_f in the text (Ang)
   ! switching function from small rho to large rho
   REAL, PARAMETER :: bs = 6.0               !< b_s in the text (Ang-1)
   REAL, PARAMETER :: rhos = 0.53            !< rho_s in the text (Ang)
!> @}

!> \name GRAPHITE POTENTIAL
!> Parameters of the force field for the graphite slab
!> @{
   ! parameters for the carbon lattice force field
   REAL, PARAMETER :: delt = 2.651           !< twisting potential parameter delta (eV)
   REAL, PARAMETER :: gam2 = 1.958           !< puckering paramter gamma_2 (eV)
   REAL, PARAMETER :: bndprm = 2.4612        !< lattice constant (in Ang)
   REAL :: rkc                               !< force constant for single carbon displacement, setup at runtime
!> @}

   ! > Coordinate of the slab in the minimum
   REAL, DIMENSION(120), SAVE :: MinSlab

   CONTAINS


! ************************************************************************************


      SUBROUTINE SetupPotential( Collinear )
         IMPLICIT NONE
         LOGICAL, OPTIONAL    :: Collinear
         REAL, DIMENSION(124) :: Positions
         INTEGER :: iCoord
         REAL    :: Value
         REAL, DIMENSION(size(Positions)) :: Gradient

         ! exit if module is setup
         IF ( PotentialModuleIsSetup ) RETURN

         ! setup force constant ( in eV per Ang^2 )
         rkc = (36.0*gam2+6.0*delt)/(bndprm**2)

         ! Setup if potential is collinear or not
         IF ( PRESENT( Collinear ) ) THEN
            CollinearPES =  Collinear
         ELSE
            CollinearPES = .FALSE.
         END IF

         PotentialModuleIsSetup = .TRUE.

         ! Define the reference 4D potential for the graphite in minimum E

         ! Set guess starting coordinate for the minimization of the slab potential
         Positions(1:124) = 0.0
         Positions(3) = C1Puckering + 2.0   ! reasonable guess for H Z coordinate
         Positions(4) = C1Puckering

         ! Minimize potential
         MinimumEnergy =  MinimizePotential( Positions, (/ (.TRUE., iCoord=1,124)  /) )      

         ! Translate to bring C3,C4,C5 in the Z=0 plane
         Value = (Positions(5)+Positions(6)+Positions(7))/3.0
         DO iCoord= 3,124
            Positions(iCoord) = Positions(iCoord) - Value
         END DO

         ! Store the coordinate of the slab
         MinSlab(:) = Positions(5:124)

         ! Store the carbon puckering and the H Z at equilibrium
         C1Puckering = Positions(4)
         HZEquilibrium = Positions(3)

         ! HZFrequency initialization
         HZForceConst = - 2. * MinimumEnergy
         ! add displacement in +delta
         Positions(3) = HZEquilibrium + Delta
         HZForceConst = HZForceConst + VHSticking( Positions, Gradient )
         ! add displacement in -delta
         Positions(3) = HZEquilibrium - Delta
         HZForceConst = HZForceConst + VHSticking( Positions, Gradient )
         ! normalize force constant
         HZForceConst = HZForceConst / Delta**2

!       PRINT*, " "
!       DO iBath = 1,8
!          WRITE(*,500) iBath+1, MinSlab(iBath)
!       END DO
!       DO iBath = 9,98
!          WRITE(*,501) iBath+1, MinSlab(iBath)
!       END DO
!       DO iBath = 99,120
!          WRITE(*,502) iBath+1, MinSlab(iBath)
!       END DO
!       500 FORMAT( " z(",I1,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       501 FORMAT( " z(",I2,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       502 FORMAT( " z(",I3,")=",F13.10, " * MyConsts_Bohr2Ang " ) 
!       PRINT*, " "
!       STOP

#if defined(VERBOSE_OUTPUT)
            WRITE(*,*) " Potential has been setup"
            WRITE(*,*) " "
#endif
      END SUBROUTINE


! ************************************************************************************

      ! Setup initial conditions for the H atom + C slab for 
      ! a simulation of vibrational relaxation
      ! The slab is fixed in the equilibrium position with no momentum ( classical 0K )
      ! The initial position and momenta of C and H are randomly chosen among a set of 
      ! conditions which are given as input
      ! data are initialized in ATOMIC UNITS
      SUBROUTINE ZeroKelvinSlabConditions( Positions, Velocities, CHInitialConditions, RandomNr )
         IMPLICIT NONE
         REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
         REAL, DIMENSION(:,:), INTENT(IN) :: CHInitialConditions
         TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
         INTEGER :: NDoF, iBath, NRandom, NInit
         REAL :: Value

         ! Check the number of non frozen degree of freedom
         NDoF = size( Positions )
         CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ZeroKelvinSlabConditions: array dimension mismatch" )

         ! Check if the nr of dimension is compatible with the slab maximum size
         CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ZeroKelvinSlabConditions: wrong number of DoFs" )

         ! Check the nr of starting conditions given ( there should be 8 coordinates: 4 positions and 4 momenta )
         NRandom = size( CHInitialConditions, 1 )
         CALL ERROR( size( CHInitialConditions, 2 ) /= 8, "PotentialModule.ZeroKelvinSlabConditions: wrong number of coords " )
      
         ! Set the velocities to zero
         Velocities(:) = 0.0

         ! Set the slab in the equilibrium geometry
         Positions(5:NDoF) = MinSlab(1:NDoF-4)

         ! Choose a random initial set of coordinates
         NInit = CEILING( UniformRandomNr(RandomNr)*real(NRandom)  )

         ! Accordingly set position and velocity
         Positions(1:4) = CHInitialConditions( NInit, 1:4 )
         Velocities(1:4) = CHInitialConditions( NInit, 5:8 )

      END SUBROUTINE ZeroKelvinSlabConditions


      ! Setup initial conditions for the H atom + C slab
      ! data are initialized in ATOMIC UNITS
      SUBROUTINE ThermalEquilibriumConditions( Positions, Velocities, Temperature, MassHydro, MassCarb, RandomNr )
         IMPLICIT NONE

         REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
         REAL, INTENT(IN)  :: Temperature, MassCarb, MassHydro
         TYPE(RNGInternalState), INTENT(INOUT) :: RandomNr
         INTEGER           :: nCarbon, NDoF
         REAL              :: SigmaCarbonVelocity, SigmaHydroVelocity

         ! All the atoms are initially at the equilibrium position for stable chemisorption 
         ! Value for the puckering are taken from J. Phys. Chem. B, 2006, 110, 18811-18817
         ! Equilibrium position of zH obtained instead from plot of the PES
         ! Velocities are sampled according to a Maxwell-Boltzmann distribution at temperature T

         ! Check the number of non frozen degree of freedom
         NDoF = size( Positions )
         CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ThermalEquilibriumConditions: array dimension mismatch" )

         ! Check if the nr of dimension is compatible with the slab maximum size
         CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ThermalEquilibriumConditions: wrong number of DoFs" )
            
         ! Equilibrium position of H atom
         Positions(1) = 0.0000
         Positions(2) = 0.0000
         Positions(3) = 1.483 / MyConsts_Bohr2Ang

         ! Equilibrium position of C1 atom
         Positions(4) = C1Puckering

         ! Equilibrium position of the other carbon atoms 
         DO nCarbon = 5,NDoF
            Positions(nCarbon)   = 0.0
         END DO

         ! Compute st deviation of Maxwell-Boltzmann distribution ( for the VELOCITY, not momenta!)
         SigmaCarbonVelocity = sqrt( Temperature / MassCarb )
         SigmaHydroVelocity  = sqrt( Temperature / MassHydro )

         ! Random velocities according to Maxwell-Boltzmann
         IF ( CollinearPES ) THEN
               Velocities(1) = 0.0
               Velocities(2) = 0.0
         ELSE 
               Velocities(1) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
               Velocities(2) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
         END IF
         Velocities(3) = GaussianRandomNr( RandomNr ) * SigmaHydroVelocity
         DO nCarbon = 4,NDoF
            Velocities(nCarbon) = GaussianRandomNr( RandomNr ) * SigmaCarbonVelocity 
         END DO
!         Velocities(4) = 0.0      ! TO FIX EVEN C1 ATOM

      END SUBROUTINE ThermalEquilibriumConditions


! ************************************************************************************

      ! Setup initial conditions for the scattering of H atom on a thermalized C slab
      ! data are initialized in ATOMIC UNITS
      SUBROUTINE ScatteringConditions( Positions, Velocities, ImpactParam, InitZ, IncEnergy, Temperature, MassHydro, MassCarb )
         IMPLICIT NONE

         REAL, DIMENSION(:), INTENT(OUT) :: Positions, Velocities
         REAL, INTENT(IN)  :: Temperature, MassCarb, MassHydro
         REAL, INTENT(IN)  :: ImpactParam, InitZ, IncEnergy
         INTEGER           :: nCarbon, NDoF
         REAL              :: SigmaMomentum, SigmaPosition
         REAL              :: gaus1, gaus2, gvar1, gvar2, gr1, gr2, gs1, gs2

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.CarbonForceConstant : Module not Setup" )

         ! Check the number of non frozen degree of freedom
         NDoF = size( Positions )
         CALL ERROR( size(Velocities) /= NDoF, "PotentialModule.ScatteringConditions: array dimension mismatch" )

         ! Check if the nr of dimension is compatible with the slab maximum size
         CALL ERROR( (NDoF > 124) .OR. (NDoF < 4), "PotentialModule.ScatteringConditions: wrong number of DoFs" )

         ! Scattering position of H atom
         Positions(1) = ImpactParam
         Positions(2) = 0.00001
         Positions(3) = InitZ

         ! Velocity of the H atom
         Velocities(1) = 0.0
         Velocities(2) = 0.0
         Velocities(3) = - sqrt( 2.0* IncEnergy / MassHydro )

         ! standard deviation of the position distribution (force constant needs to be in AU)
         SigmaPosition = sqrt( Temperature / (rkc*(MyConsts_Bohr2Ang)**2/MyConsts_Hartree2eV) )
         ! standard deviation of the momentum distribution
         SigmaMomentum = sqrt( MassCarb * Temperature )

         ! Cycle over carbon atoms in the slab
         DO nCarbon = 4,NDoF

            ! Initialization
            Positions(nCarbon)   = 0.0
            Velocities(nCarbon)   = 0.0

            ! Generate gaussian random numbers for position and velocity
            DO 
               call random_number(gaus1)
               call random_number(gaus2)
               gvar1=2.0*gaus1-1
               gvar2=2.0*gaus2-1
               gr1=gvar1**2+gvar2**2
               IF (gr1 < 1.0) EXIT
            END DO
            gr2=sqrt(-2.0*alog(gr1)/gr1)
            gs1=gvar1*gr2
            gs2=gvar2*gr2

            ! Set position and velocity of the carbon atom
            Positions(nCarbon)  = SigmaPosition*gs1
            Velocities(nCarbon) = SigmaMomentum*gs2/MassCarb

         END DO

      END SUBROUTINE ScatteringConditions


! ************************************************************************************

      REAL FUNCTION GraphiteLatticeConstant()
         IMPLICIT NONE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GraphiteLatticeConstant : Module not Setup" )
         ! Convert force constant from eV per Ang^2 to Hartree per bohr^2
         GraphiteLatticeConstant = bndprm / MyConsts_Bohr2Ang 

      END FUNCTION GraphiteLatticeConstant

! ******************************************************************************************      

      REAL FUNCTION CarbonForceConstant()
         IMPLICIT NONE

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.GraphiteLatticeConstant : Module not Setup" )
         ! Convert force constant from eV per Ang^2 to Hartree per bohr^2
         CarbonForceConstant = rkc * (MyConsts_Bohr2Ang)**2 / MyConsts_Hartree2eV

      END FUNCTION CarbonForceConstant

! ******************************************************************************************      

      REAL FUNCTION MinimizePotential( Coords, Mask ) RESULT( Pot )
         IMPLICIT NONE
         REAL, INTENT(INOUT), DIMENSION(:)               :: Coords
         LOGICAL, INTENT(IN), DIMENSION(size(Coords)) :: Mask

         INTEGER :: NrDimension, NrOptimization
         INTEGER :: iIter, iCoord
         REAL, DIMENSION(size(Coords)) :: Gradient
         REAL :: Norm

         ! Set dimension number
         NrDimension = size(Coords)
         ! Set optimization coordinates nr
         NrOptimization = count( Mask )
         ! Check if the nr of dimension is compatible with the slab maximum size
         CALL ERROR( (NrDimension > 124) .OR. (NrDimension < 4), "PotentialModule.MinimizePotential: wrong number of DoFs" )

         ! Cycle over steepest descent iterations
         DO iIter = 1, MaxIter

            ! compute negative of the gradient
            Pot = VHSticking( Coords, Gradient )

            ! compute norm of the gradient
            Norm = 0.0
            DO iCoord = 1, NrDimension
               IF ( Mask( iCoord ) ) THEN
                  Norm = Norm + Gradient(iCoord)**2
               END IF
            END DO
            Norm = SQRT( Norm / NrOptimization )

            ! check convergence
            IF (Norm < GradEps) EXIT
      
            ! move geometry along gradient
            DO iCoord = 1, NrDimension
               IF ( Mask( iCoord ) ) THEN
                  Coords(iCoord) = Coords(iCoord) + Gradient(iCoord)
               END IF
            END DO

         END DO

#if defined(VERBOSE_OUTPUT)
         WRITE(*,"(/,A,I6,A)") "Convergence in ", iIter, " steps"
#endif
         CALL WARN( iIter == MaxIter, "PotentialModule. MinimizePotential: convergence not reached" )

      END FUNCTION MinimizePotential


! ******************************************************************************************      

!*******************************************************************************
!> 4D adiabatic H-Graphene potential 
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom and 
!>                     1 Z coordinates for the first carbon atoms (in au)
!> @param Forces       Output array with the derivatives of the potential (in au)
!> @param vv           output potential in atomic units
!*******************************************************************************   
      REAL FUNCTION VHFourDimensional( Positions, Forces ) RESULT(V) 
         IMPLICIT NONE

         REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
         REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 
         REAL, DIMENSION(124) :: Dummy

         Dummy = 0.0
         ! use the full potential with the carbon atoms in the equilibrium geometry
         V = VHSticking( (/ Positions, MinSlab(:) /), Dummy(:) ) 
         Forces(1:4) = Dummy(1:4)
         IF ( CollinearPES ) THEN 
            Forces(1:2) = 0.0
         END IF

      END FUNCTION VHFourDimensional


!*******************************************************************************
!> H-Graphene potential by Jackson and coworkers.
!> @ref http://pubs.acs.org/doi/abs/10.1021/jp057136%2B [ J.Phys.Chem.B 2006,110,18811 ]
!> @ref other papers jackson
!>
!> @param Positions    Array with 3 cartesian coordinates for the H atom and 
!>                     121 Z coordinates for the carbon atoms (in au)
!> @param Forces       Output array with the derivatives of the potential (in au)
!> @param vv           output potential in atomic units
!*******************************************************************************     
      REAL FUNCTION VHSticking( Positions, Forces ) RESULT(vv) 
         IMPLICIT NONE

         REAL, DIMENSION(:), TARGET, INTENT(IN)  :: Positions
         REAL, DIMENSION(:), TARGET, INTENT(OUT) :: Forces 

         ! temporary variables to store positions in Ang units
         REAL :: xh, yh, zh    
         REAL, DIMENSION(121) :: z

         ! Number of non frozen degrees of freedom
         INTEGER :: NrNonFrozen        
         
         REAL, DIMENSION(121) :: ddz

         REAL :: a, b, c1, c2, d0

         REAL :: dbdzc, dc1dzc, dc2dzc, dd0dzc, df2dr, dkrdzc, dkrdzh
         REAL :: dswdr, dvdzc, dvdzh, dvidzc, dvidzh, dvqdr, dvqdzc
         REAL :: dvqdzh, dvtdrho, dvtds1, dvtds2
         REAL :: dzgdzc, dzmdzc

         REAL :: fexp, ff1, ff2, ff3
         REAL :: qqq, rho, rkrho, sub1, sub2, sw
         REAL :: v, vi, vlatt, vpuck, vq, vt, vtors
         REAL :: zg, zm

         INTEGER :: nn

         ! Error if module not have been setup yet
         CALL ERROR( .NOT. PotentialModuleIsSetup, "PotentialModule.VHSticking : Module not Setup" )

         ! Check the number of non frozen degree of freedom
         NrNonFrozen = size( Positions )
         CALL ERROR( size(Forces) /= NrNonFrozen, "PotentialModule.VHSticking: array dimension mismatch" )

         ! Check if the nr of dimension is compatible with the slab maximum size
         CALL ERROR( (NrNonFrozen > 124) .OR. (NrNonFrozen < 4), "PotentialModule.VHSticking: wrong number of DoFs" )
         
         ! Transform input coordinate from AU to Angstrom
         IF ( CollinearPES ) THEN
            xh = 0.0
            yh = 0.0
         ELSE
            xh = Positions(1) * MyConsts_Bohr2Ang
            yh = Positions(2) * MyConsts_Bohr2Ang
         ENDIF
         zh = Positions(3) * MyConsts_Bohr2Ang
         z(:) = 0.0 ! frozen coordinates are set to zero
         z(1:NrNonFrozen-3) = Positions(4:NrNonFrozen) * MyConsts_Bohr2Ang

         ! Compute some relevant coordinates for the PES

         ! Q: average Z of C2, C3 and C4 - system-lattice coupling coordinate
         qqq=(z(2)+z(3)+z(4))/3.0

         ! positions of H and C1 with respect to Q
         sub1=zh-qqq
         sub2=z(1)-qqq

         ! rho
         rho=sqrt(xh**2+yh**2)

   ! **************************************************************************************************
   !                           POTENTIAL FOR THE C-H SYSTEM
   ! **************************************************************************************************

         ! Compute the parameters for the morse + gaussian 
         ! functional form for the collinear potential

         ! D_0 is the morse depth (eV)
         d0=0.474801+0.9878257*sub2-1.3921499*sub2**2              &
         +0.028278*sub2**3-1.8879928*sub2**4                       &
         +0.11*exp(-8.0*(sub2-0.28)**2)
         dd0dzc=0.9878257-2.7842998*sub2+0.084834*sub2**2-         &
         7.5519712*sub2**3+(                                       &
         -1.76*(sub2-0.28)*exp(-8.0*(sub2-0.28)**2))

         ! A is the morse curvature (Ang)
         a=2.276211

         ! Z_M is the morse equilibrium distance (Ang)
         zm=0.87447*(sub2)+1.17425
         dzmdzc=0.87447

         ! C_1 is the asympotic harmonic potential for C1 vibrations (eV)
         c1=0.5*rkc*(sub2)**2-0.00326
         dc1dzc=rkc*(sub2)

         ! C_2 is the gaussian height (eV)
         c2= 0.3090344*exp(-2.741813*(sub2-0.2619756))+            &
         0.03113325*exp(3.1844857*(sub2-0.186741))+                &
         0.02*exp(-20.0*(sub2-0.1)**2) 
         dc2dzc=-0.8473145*exp(-2.741813*(sub2-0.2619756))+        &
         0.0991434*exp(3.1844857*(sub2-0.186741))+(                &
         -0.8*(sub2-0.1)*exp(-20.0*(sub2-0.1)**2))

         ! B is the gaussian curvature (Ang-2)
         b=4.00181*exp(1.25965*(sub2-0.58729)) 
         dbdzc=5.0408799*exp(1.25965*(sub2-0.58729))

         ! Z_G is the center of the gaussian (Ang)
         zg=1.99155*(sub2)+1.46095
         dzgdzc=1.99155

         ! Compute the potential and the derivatives of the 
         ! collinear potential V_0

         v=c1+(c1+d0)*(exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+       &
         c2*exp(-b*(sub1-zg)**2)

         dvdzh=-b*c2*(sub1-zg)*2.0*exp(-b*(sub1-zg)**2) +                  &
         (c1+d0)*2.0*(-a*exp(-2.0*a*(sub1-zm))+a*exp(-a*(sub1-zm)))

         dvdzc=dc1dzc+(dc1dzc+dd0dzc)*                             &
         (exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+            &
         (c1+d0)*(a*dzmdzc*2.0*exp(-2.0*a*(sub1-zm))+              &
         (-a*dzmdzc)*2.0*exp(-a*(sub1-zm)))+                       &
         dc2dzc*exp(-b*(sub1-zg)**2)+                              &
         (-c2*dbdzc*(sub1-zg)**2+c2*b*2.0*(sub1-zg)*dzgdzc)*       &
         exp(-b*(sub1-zg)**2)

         ! Compute the force constant (and derivatives) for small rho 
         ! potential rkrho(zh-q,zc-q)

         rkrho=3.866259*exp(-17.038588*(sub2-0.537621)**2+         &
         0.312355*(sub2-0.537621)*(sub1-2.003753)-                 &
         4.479864*(sub1-2.003753)**2)+                             &
         4.317415*exp(-11.931770*(sub2-0.286858)**2+               &
         18.540974*(sub2-0.286858)*(sub1-1.540947)-                &
         14.537321*(sub1-1.540947)**2)

         dkrdzc=(-34.077176*(sub2-0.537621)+0.312355*              &
         (sub1-2.003753))                                          &
         *3.866259*exp(-17.038588*(sub2-0.537621)**2+              &
         0.312355*(sub2-0.537621)*(sub1-2.003753)-                 &
         4.479864*(sub1-2.003753)**2)+                             &
         (-23.86354*(sub2-0.286858)+18.540974*(sub1-1.540947))     &
         *4.317415*exp(-11.931770*(sub2-0.286858)**2+              &
         18.540974*(sub2-0.286858)*(sub1-1.540947)-                &
         14.537321*(sub1-1.540947)**2)

         dkrdzh=(0.312355*(sub2-0.537621)-8.959728*(sub1-2.003753))        &
         *3.866259*exp(-17.038588*(sub2-0.537621)**2+                      &
         0.312355*(sub2-0.537621)*(sub1-2.003753)-                         &
         4.479864*(sub1-2.003753)**2)+                                     &
         (18.540974*(sub2-0.286858)-29.074642*(sub1-1.540947))             &
         *4.317415*exp(-11.931770*(sub2-0.286858)**2+                      &
         18.540974*(sub2-0.286858)*(sub1-1.540947)-                        &
         14.537321*(sub1-1.540947)**2)

         ! Compute the small rho potential/derivatives

         vq=v+0.5*rkrho*rho**2
         dvqdzh=dvdzh+0.5*dkrdzh*rho**2
         dvqdzc=dvdzc+0.5*dkrdzc*rho**2
         dvqdr=rkrho*rho

         ! Compute the  "infinite" rho potential/derivatives

         vi=0.5*rkc*sub2**2-0.00326+                                       &
         di*(exp(-2.0*alphi*(sub1-ai))-2.0*exp(-alphi*(sub1-ai)))
         dvidzh=di*(-2.0*alphi*exp(-2.0*alphi*(sub1-ai))+                  &
               2.0*alphi*exp(-alphi*(sub1-ai)))
         dvidzc=rkc*sub2

         ! Switching function and associated functions

         fexp=exp(-ba*(rho-rhoa))
         ff1=1.0+fexp
         ff2=exp(-2.0*alp*rho)-2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)
         ff3=(vi-v)*ff2+vi
         sw=1.0/(1.0+exp(-bs*(rho-rhos)))

         ! Total H,C1 potential/derivatives

         vt=vq*(1.0-sw)+ff3*sw
      
         df2dr=-2.0*alp*exp(-2.0*alp*rho)-                &
            2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)*       &
            (-alp-(alp2/ff1)-alp2*rho*ba*fexp/(ff1**2))
         dswdr=(bs*exp(-bs*(rho-rhos)))/((1.0+exp(-bs*(rho-rhos)))**2)

         dvtds1=dvqdzh*(1.0-sw)+sw*((dvidzh-dvdzh)*ff2+dvidzh)
         dvtds2=dvqdzc*(1.0-sw)+sw*((dvidzc-dvdzc)*ff2+dvidzc)
         dvtdrho=dvqdr*(1.0-sw)+vq*(-dswdr)+sw*(vi-v)*df2dr+ff3*dswdr

   ! **************************************************************************************************
   !                           POTENTIAL FOR THE GRAPHITE LATTICE
   ! **************************************************************************************************

         DO nn = 1, 121
            ddz(nn) = 0.0
         END DO

         vpuck=0.0
         vtors=0.0

         vtors=vtors+(z(1)-z(5)-z(12)+z(15))**2
         ddz(1)=ddz(1)+(delt/bndprm**2)*(z(1)-z(5)-z(12)+z(15))
         ddz(5)=ddz(5)-(delt/bndprm**2)*(z(1)-z(5)-z(12)+z(15))
         ddz(12)=ddz(12)-(delt/bndprm**2)*(z(1)-z(5)-z(12)+z(15))
         ddz(15)=ddz(15)+(delt/bndprm**2)*(z(1)-z(5)-z(12)+z(15))

         vtors=vtors+(z(1)-z(7)-z(13)+z(18))**2
         ddz(1)=ddz(1)+(delt/bndprm**2)*(z(1)-z(7)-z(13)+z(18))
         ddz(7)=ddz(7)-(delt/bndprm**2)*(z(1)-z(7)-z(13)+z(18))
         ddz(13)=ddz(13)-(delt/bndprm**2)*(z(1)-z(7)-z(13)+z(18))
         ddz(18)=ddz(18)+(delt/bndprm**2)*(z(1)-z(7)-z(13)+z(18))

         vtors=vtors+(z(1)-z(9)-z(11)+z(16))**2
         ddz(1)=ddz(1)+(delt/bndprm**2)*(z(1)-z(9)-z(11)+z(16))
         ddz(9)=ddz(9)-(delt/bndprm**2)*(z(1)-z(9)-z(11)+z(16))
         ddz(11)=ddz(11)-(delt/bndprm**2)*(z(1)-z(9)-z(11)+z(16))
         ddz(16)=ddz(16)+(delt/bndprm**2)*(z(1)-z(9)-z(11)+z(16))

         vtors=vtors+(z(2)-z(3)-z(8)+z(10))**2
         ddz(2)=ddz(2)+(delt/bndprm**2)*(z(2)-z(3)-z(8)+z(10))
         ddz(3)=ddz(3)-(delt/bndprm**2)*(z(2)-z(3)-z(8)+z(10))
         ddz(8)=ddz(8)-(delt/bndprm**2)*(z(2)-z(3)-z(8)+z(10))
         ddz(10)=ddz(10)+(delt/bndprm**2)*(z(2)-z(3)-z(8)+z(10))

         vtors=vtors+(z(2)-z(11)-z(20)+z(27))**2
         ddz(2)=ddz(2)+(delt/bndprm**2)*(z(2)-z(11)-z(20)+z(27))
         ddz(11)=ddz(11)-(delt/bndprm**2)*(z(2)-z(11)-z(20)+z(27))
         ddz(20)=ddz(20)-(delt/bndprm**2)*(z(2)-z(11)-z(20)+z(27))
         ddz(27)=ddz(27)+(delt/bndprm**2)*(z(2)-z(11)-z(20)+z(27))

         vtors=vtors+(z(2)-z(14)-z(7)+z(21))**2
         ddz(2)=ddz(2)+(delt/bndprm**2)*(z(2)-z(14)-z(7)+z(21))
         ddz(14)=ddz(14)-(delt/bndprm**2)*(z(2)-z(14)-z(7)+z(21))
         ddz(7)=ddz(7)-(delt/bndprm**2)*(z(2)-z(14)-z(7)+z(21))
         ddz(21)=ddz(21)+(delt/bndprm**2)*(z(2)-z(14)-z(7)+z(21))

         vtors=vtors+(z(3)-z(11)-z(22)+z(28))**2
         ddz(3)=ddz(3)+(delt/bndprm**2)*(z(3)-z(11)-z(22)+z(28))
         ddz(11)=ddz(11)-(delt/bndprm**2)*(z(3)-z(11)-z(22)+z(28))
         ddz(22)=ddz(22)-(delt/bndprm**2)*(z(3)-z(11)-z(22)+z(28))
         ddz(28)=ddz(28)+(delt/bndprm**2)*(z(3)-z(11)-z(22)+z(28))

         vtors=vtors+(z(3)-z(16)-z(5)+z(21))**2
         ddz(3)=ddz(3)+(delt/bndprm**2)*(z(3)-z(16)-z(5)+z(21))
         ddz(16)=ddz(16)-(delt/bndprm**2)*(z(3)-z(16)-z(5)+z(21))
         ddz(5)=ddz(5)-(delt/bndprm**2)*(z(3)-z(16)-z(5)+z(21))
         ddz(21)=ddz(21)+(delt/bndprm**2)*(z(3)-z(16)-z(5)+z(21))

         vtors=vtors+(z(3)-z(18)-z(10)+z(24))**2
         ddz(3)=ddz(3)+(delt/bndprm**2)*(z(3)-z(18)-z(10)+z(24))
         ddz(18)=ddz(18)-(delt/bndprm**2)*(z(3)-z(18)-z(10)+z(24))
         ddz(10)=ddz(10)-(delt/bndprm**2)*(z(3)-z(18)-z(10)+z(24))
         ddz(24)=ddz(24)+(delt/bndprm**2)*(z(3)-z(18)-z(10)+z(24))

         vtors=vtors+(z(4)-z(2)-z(9)+z(7))**2
         ddz(4)=ddz(4)+(delt/bndprm**2)*(z(4)-z(2)-z(9)+z(7))
         ddz(2)=ddz(2)-(delt/bndprm**2)*(z(4)-z(2)-z(9)+z(7))
         ddz(9)=ddz(9)-(delt/bndprm**2)*(z(4)-z(2)-z(9)+z(7))
         ddz(7)=ddz(7)+(delt/bndprm**2)*(z(4)-z(2)-z(9)+z(7))

         vtors=vtors+(z(4)-z(3)-z(6)+z(5))**2
         ddz(4)=ddz(4)+(delt/bndprm**2)*(z(4)-z(3)-z(6)+z(5))
         ddz(3)=ddz(3)-(delt/bndprm**2)*(z(4)-z(3)-z(6)+z(5))
         ddz(6)=ddz(6)-(delt/bndprm**2)*(z(4)-z(3)-z(6)+z(5))
         ddz(5)=ddz(5)+(delt/bndprm**2)*(z(4)-z(3)-z(6)+z(5))

         vtors=vtors+(z(4)-z(13)-z(25)+z(31))**2
         ddz(4)=ddz(4)+(delt/bndprm**2)*(z(4)-z(13)-z(25)+z(31))
         ddz(13)=ddz(13)-(delt/bndprm**2)*(z(4)-z(13)-z(25)+z(31))
         ddz(25)=ddz(25)-(delt/bndprm**2)*(z(4)-z(13)-z(25)+z(31))
         ddz(31)=ddz(31)+(delt/bndprm**2)*(z(4)-z(13)-z(25)+z(31))

         vtors=vtors+(z(5)-z(7)-z(34)+z(32))**2
         ddz(5)=ddz(5)+(delt/bndprm**2)*(z(5)-z(7)-z(34)+z(32))
         ddz(7)=ddz(7)-(delt/bndprm**2)*(z(5)-z(7)-z(34)+z(32))
         ddz(34)=ddz(34)-(delt/bndprm**2)*(z(5)-z(7)-z(34)+z(32))
         ddz(32)=ddz(32)+(delt/bndprm**2)*(z(5)-z(7)-z(34)+z(32))

         vtors=vtors+(z(5)-z(27)-z(15)+z(38))**2
         ddz(5)=ddz(5)+(delt/bndprm**2)*(z(5)-z(27)-z(15)+z(38))
         ddz(27)=ddz(27)-(delt/bndprm**2)*(z(5)-z(27)-z(15)+z(38))
         ddz(15)=ddz(15)-(delt/bndprm**2)*(z(5)-z(27)-z(15)+z(38))
         ddz(38)=ddz(38)+(delt/bndprm**2)*(z(5)-z(27)-z(15)+z(38))

         vtors=vtors+(z(6)-z(1)-z(14)+z(11))**2
         ddz(6)=ddz(6)+(delt/bndprm**2)*(z(6)-z(1)-z(14)+z(11))
         ddz(1)=ddz(1)-(delt/bndprm**2)*(z(6)-z(1)-z(14)+z(11))
         ddz(14)=ddz(14)-(delt/bndprm**2)*(z(6)-z(1)-z(14)+z(11))
         ddz(11)=ddz(11)+(delt/bndprm**2)*(z(6)-z(1)-z(14)+z(11))

         vtors=vtors+(z(6)-z(20)-z(35)+z(43))**2
         ddz(6)=ddz(6)+(delt/bndprm**2)*(z(6)-z(20)-z(35)+z(43))
         ddz(20)=ddz(20)-(delt/bndprm**2)*(z(6)-z(20)-z(35)+z(43))
         ddz(35)=ddz(35)-(delt/bndprm**2)*(z(6)-z(20)-z(35)+z(43))
         ddz(43)=ddz(43)+(delt/bndprm**2)*(z(6)-z(20)-z(35)+z(43))

         vtors=vtors+(z(7)-z(22)-z(32)+z(45))**2
         ddz(7)=ddz(7)+(delt/bndprm**2)*(z(7)-z(22)-z(32)+z(45))
         ddz(22)=ddz(22)-(delt/bndprm**2)*(z(7)-z(22)-z(32)+z(45))
         ddz(32)=ddz(32)-(delt/bndprm**2)*(z(7)-z(22)-z(32)+z(45))
         ddz(45)=ddz(45)+(delt/bndprm**2)*(z(7)-z(22)-z(32)+z(45))

         vtors=vtors+(z(7)-z(28)-z(18)+z(39))**2
         ddz(7)=ddz(7)+(delt/bndprm**2)*(z(7)-z(28)-z(18)+z(39))
         ddz(28)=ddz(28)-(delt/bndprm**2)*(z(7)-z(28)-z(18)+z(39))
         ddz(18)=ddz(18)-(delt/bndprm**2)*(z(7)-z(28)-z(18)+z(39))
         ddz(39)=ddz(39)+(delt/bndprm**2)*(z(7)-z(28)-z(18)+z(39))

         vtors=vtors+(z(8)-z(1)-z(17)+z(13))**2
         ddz(8)=ddz(8)+(delt/bndprm**2)*(z(8)-z(1)-z(17)+z(13))
         ddz(1)=ddz(1)-(delt/bndprm**2)*(z(8)-z(1)-z(17)+z(13))
         ddz(17)=ddz(17)-(delt/bndprm**2)*(z(8)-z(1)-z(17)+z(13))
         ddz(13)=ddz(13)+(delt/bndprm**2)*(z(8)-z(1)-z(17)+z(13))

         vtors=vtors+(z(8)-z(6)-z(33)+z(35))**2
         ddz(8)=ddz(8)+(delt/bndprm**2)*(z(8)-z(6)-z(33)+z(35))
         ddz(6)=ddz(6)-(delt/bndprm**2)*(z(8)-z(6)-z(33)+z(35))
         ddz(33)=ddz(33)-(delt/bndprm**2)*(z(8)-z(6)-z(33)+z(35))
         ddz(35)=ddz(35)+(delt/bndprm**2)*(z(8)-z(6)-z(33)+z(35))

         vtors=vtors+(z(9)-z(22)-z(36)+z(46))**2
         ddz(9)=ddz(9)+(delt/bndprm**2)*(z(9)-z(22)-z(36)+z(46))
         ddz(22)=ddz(22)-(delt/bndprm**2)*(z(9)-z(22)-z(36)+z(46))
         ddz(36)=ddz(36)-(delt/bndprm**2)*(z(9)-z(22)-z(36)+z(46))
         ddz(46)=ddz(46)+(delt/bndprm**2)*(z(9)-z(22)-z(36)+z(46))

         vtors=vtors+(z(9)-z(30)-z(16)+z(39))**2
         ddz(9)=ddz(9)+(delt/bndprm**2)*(z(9)-z(30)-z(16)+z(39))
         ddz(30)=ddz(30)-(delt/bndprm**2)*(z(9)-z(30)-z(16)+z(39))
         ddz(16)=ddz(16)-(delt/bndprm**2)*(z(9)-z(30)-z(16)+z(39))
         ddz(39)=ddz(39)+(delt/bndprm**2)*(z(9)-z(30)-z(16)+z(39))

         vtors=vtors+(z(10)-z(1)-z(19)+z(12))**2
         ddz(10)=ddz(10)+(delt/bndprm**2)*(z(10)-z(1)-z(19)+z(12))
         ddz(1)=ddz(1)-(delt/bndprm**2)*(z(10)-z(1)-z(19)+z(12))
         ddz(19)=ddz(19)-(delt/bndprm**2)*(z(10)-z(1)-z(19)+z(12))
         ddz(12)=ddz(12)+(delt/bndprm**2)*(z(10)-z(1)-z(19)+z(12))

         vtors=vtors+(z(10)-z(9)-z(37)+z(36))**2
         ddz(10)=ddz(10)+(delt/bndprm**2)*(z(10)-z(9)-z(37)+z(36))
         ddz(9)=ddz(9)-(delt/bndprm**2)*(z(10)-z(9)-z(37)+z(36))
         ddz(37)=ddz(37)-(delt/bndprm**2)*(z(10)-z(9)-z(37)+z(36))
         ddz(36)=ddz(36)+(delt/bndprm**2)*(z(10)-z(9)-z(37)+z(36))

         vtors=vtors+(z(11)-z(32)-z(27)+z(52))**2
         ddz(11)=ddz(11)+(delt/bndprm**2)*(z(11)-z(32)-z(27)+z(52))
         ddz(32)=ddz(32)-(delt/bndprm**2)*(z(11)-z(32)-z(27)+z(52))
         ddz(27)=ddz(27)-(delt/bndprm**2)*(z(11)-z(32)-z(27)+z(52))
         ddz(52)=ddz(52)+(delt/bndprm**2)*(z(11)-z(32)-z(27)+z(52))

         vtors=vtors+(z(11)-z(34)-z(28)+z(51))**2
         ddz(11)=ddz(11)+(delt/bndprm**2)*(z(11)-z(34)-z(28)+z(51))
         ddz(34)=ddz(34)-(delt/bndprm**2)*(z(11)-z(34)-z(28)+z(51))
         ddz(28)=ddz(28)-(delt/bndprm**2)*(z(11)-z(34)-z(28)+z(51))
         ddz(51)=ddz(51)+(delt/bndprm**2)*(z(11)-z(34)-z(28)+z(51))

         vtors=vtors+(z(12)-z(2)-z(26)+z(20))**2
         ddz(12)=ddz(12)+(delt/bndprm**2)*(z(12)-z(2)-z(26)+z(20))
         ddz(2)=ddz(2)-(delt/bndprm**2)*(z(12)-z(2)-z(26)+z(20))
         ddz(26)=ddz(26)-(delt/bndprm**2)*(z(12)-z(2)-z(26)+z(20))
         ddz(20)=ddz(20)+(delt/bndprm**2)*(z(12)-z(2)-z(26)+z(20))

         vtors=vtors+(z(12)-z(4)-z(29)+z(25))**2
         ddz(12)=ddz(12)+(delt/bndprm**2)*(z(12)-z(4)-z(29)+z(25))
         ddz(4)=ddz(4)-(delt/bndprm**2)*(z(12)-z(4)-z(29)+z(25))
         ddz(29)=ddz(29)-(delt/bndprm**2)*(z(12)-z(4)-z(29)+z(25))
         ddz(25)=ddz(25)+(delt/bndprm**2)*(z(12)-z(4)-z(29)+z(25))

         vtors=vtors+(z(13)-z(3)-z(30)+z(22))**2
         ddz(13)=ddz(13)+(delt/bndprm**2)*(z(13)-z(3)-z(30)+z(22))
         ddz(3)=ddz(3)-(delt/bndprm**2)*(z(13)-z(3)-z(30)+z(22))
         ddz(30)=ddz(30)-(delt/bndprm**2)*(z(13)-z(3)-z(30)+z(22))
         ddz(22)=ddz(22)+(delt/bndprm**2)*(z(13)-z(3)-z(30)+z(22))

         vtors=vtors+(z(13)-z(36)-z(31)+z(58))**2
         ddz(13)=ddz(13)+(delt/bndprm**2)*(z(13)-z(36)-z(31)+z(58))
         ddz(36)=ddz(36)-(delt/bndprm**2)*(z(13)-z(36)-z(31)+z(58))
         ddz(31)=ddz(31)-(delt/bndprm**2)*(z(13)-z(36)-z(31)+z(58))
         ddz(58)=ddz(58)+(delt/bndprm**2)*(z(13)-z(36)-z(31)+z(58))

         vtors=vtors+(z(14)-z(34)-z(47)+z(62))**2
         ddz(14)=ddz(14)+(delt/bndprm**2)*(z(14)-z(34)-z(47)+z(62))
         ddz(34)=ddz(34)-(delt/bndprm**2)*(z(14)-z(34)-z(47)+z(62))
         ddz(47)=ddz(47)-(delt/bndprm**2)*(z(14)-z(34)-z(47)+z(62))
         ddz(62)=ddz(62)+(delt/bndprm**2)*(z(14)-z(34)-z(47)+z(62))

         vtors=vtors+(z(14)-z(42)-z(21)+z(52))**2
         ddz(14)=ddz(14)+(delt/bndprm**2)*(z(14)-z(42)-z(21)+z(52))
         ddz(42)=ddz(42)-(delt/bndprm**2)*(z(14)-z(42)-z(21)+z(52))
         ddz(21)=ddz(21)-(delt/bndprm**2)*(z(14)-z(42)-z(21)+z(52))
         ddz(52)=ddz(52)+(delt/bndprm**2)*(z(14)-z(42)-z(21)+z(52))

         vtors=vtors+(z(15)-z(2)-z(23)+z(8))**2
         ddz(15)=ddz(15)+(delt/bndprm**2)*(z(15)-z(2)-z(23)+z(8))
         ddz(2)=ddz(2)-(delt/bndprm**2)*(z(15)-z(2)-z(23)+z(8))
         ddz(23)=ddz(23)-(delt/bndprm**2)*(z(15)-z(2)-z(23)+z(8))
         ddz(8)=ddz(8)+(delt/bndprm**2)*(z(15)-z(2)-z(23)+z(8))

         vtors=vtors+(z(15)-z(14)-z(48)+z(47))**2
         ddz(15)=ddz(15)+(delt/bndprm**2)*(z(15)-z(14)-z(48)+z(47))
         ddz(14)=ddz(14)-(delt/bndprm**2)*(z(15)-z(14)-z(48)+z(47))
         ddz(48)=ddz(48)-(delt/bndprm**2)*(z(15)-z(14)-z(48)+z(47))
         ddz(47)=ddz(47)+(delt/bndprm**2)*(z(15)-z(14)-z(48)+z(47))

         vtors=vtors+(z(16)-z(32)-z(56)+z(64))**2
         ddz(16)=ddz(16)+(delt/bndprm**2)*(z(16)-z(32)-z(56)+z(64))
         ddz(32)=ddz(32)-(delt/bndprm**2)*(z(16)-z(32)-z(56)+z(64))
         ddz(56)=ddz(56)-(delt/bndprm**2)*(z(16)-z(32)-z(56)+z(64))
         ddz(64)=ddz(64)+(delt/bndprm**2)*(z(16)-z(32)-z(56)+z(64))

         vtors=vtors+(z(16)-z(45)-z(21)+z(51))**2
         ddz(16)=ddz(16)+(delt/bndprm**2)*(z(16)-z(45)-z(21)+z(51))
         ddz(45)=ddz(45)-(delt/bndprm**2)*(z(16)-z(45)-z(21)+z(51))
         ddz(21)=ddz(21)-(delt/bndprm**2)*(z(16)-z(45)-z(21)+z(51))
         ddz(51)=ddz(51)+(delt/bndprm**2)*(z(16)-z(45)-z(21)+z(51))

         vtors=vtors+(z(17)-z(4)-z(24)+z(9))**2
         ddz(17)=ddz(17)+(delt/bndprm**2)*(z(17)-z(4)-z(24)+z(9))
         ddz(4)=ddz(4)-(delt/bndprm**2)*(z(17)-z(4)-z(24)+z(9))
         ddz(24)=ddz(24)-(delt/bndprm**2)*(z(17)-z(4)-z(24)+z(9))
         ddz(9)=ddz(9)+(delt/bndprm**2)*(z(17)-z(4)-z(24)+z(9))

         vtors=vtors+(z(17)-z(37)-z(55)+z(66))**2
         ddz(17)=ddz(17)+(delt/bndprm**2)*(z(17)-z(37)-z(55)+z(66))
         ddz(37)=ddz(37)-(delt/bndprm**2)*(z(17)-z(37)-z(55)+z(66))
         ddz(55)=ddz(55)-(delt/bndprm**2)*(z(17)-z(37)-z(55)+z(66))
         ddz(66)=ddz(66)+(delt/bndprm**2)*(z(17)-z(37)-z(55)+z(66))

         vtors=vtors+(z(18)-z(16)-z(54)+z(56))**2
         ddz(18)=ddz(18)+(delt/bndprm**2)*(z(18)-z(16)-z(54)+z(56))
         ddz(16)=ddz(16)-(delt/bndprm**2)*(z(18)-z(16)-z(54)+z(56))
         ddz(54)=ddz(54)-(delt/bndprm**2)*(z(18)-z(16)-z(54)+z(56))
         ddz(56)=ddz(56)+(delt/bndprm**2)*(z(18)-z(16)-z(54)+z(56))

         vtors=vtors+(z(18)-z(46)-z(24)+z(57))**2
         ddz(18)=ddz(18)+(delt/bndprm**2)*(z(18)-z(46)-z(24)+z(57))
         ddz(46)=ddz(46)-(delt/bndprm**2)*(z(18)-z(46)-z(24)+z(57))
         ddz(24)=ddz(24)-(delt/bndprm**2)*(z(18)-z(46)-z(24)+z(57))
         ddz(57)=ddz(57)+(delt/bndprm**2)*(z(18)-z(46)-z(24)+z(57))

         vtors=vtors+(z(19)-z(4)-z(23)+z(6))**2
         ddz(19)=ddz(19)+(delt/bndprm**2)*(z(19)-z(4)-z(23)+z(6))
         ddz(4)=ddz(4)-(delt/bndprm**2)*(z(19)-z(4)-z(23)+z(6))
         ddz(23)=ddz(23)-(delt/bndprm**2)*(z(19)-z(4)-z(23)+z(6))
         ddz(6)=ddz(6)+(delt/bndprm**2)*(z(19)-z(4)-z(23)+z(6))

         vtors=vtors+(z(19)-z(17)-z(50)+z(55))**2
         ddz(19)=ddz(19)+(delt/bndprm**2)*(z(19)-z(17)-z(50)+z(55))
         ddz(17)=ddz(17)-(delt/bndprm**2)*(z(19)-z(17)-z(50)+z(55))
         ddz(50)=ddz(50)-(delt/bndprm**2)*(z(19)-z(17)-z(50)+z(55))
         ddz(55)=ddz(55)+(delt/bndprm**2)*(z(19)-z(17)-z(50)+z(55))

         vtors=vtors+(z(20)-z(5)-z(42)+z(34))**2
         ddz(20)=ddz(20)+(delt/bndprm**2)*(z(20)-z(5)-z(42)+z(34))
         ddz(5)=ddz(5)-(delt/bndprm**2)*(z(20)-z(5)-z(42)+z(34))
         ddz(42)=ddz(42)-(delt/bndprm**2)*(z(20)-z(5)-z(42)+z(34))
         ddz(34)=ddz(34)+(delt/bndprm**2)*(z(20)-z(5)-z(42)+z(34))

         vtors=vtors+(z(20)-z(47)-z(43)+z(76))**2
         ddz(20)=ddz(20)+(delt/bndprm**2)*(z(20)-z(47)-z(43)+z(76))
         ddz(47)=ddz(47)-(delt/bndprm**2)*(z(20)-z(47)-z(43)+z(76))
         ddz(43)=ddz(43)-(delt/bndprm**2)*(z(20)-z(47)-z(43)+z(76))
         ddz(76)=ddz(76)+(delt/bndprm**2)*(z(20)-z(47)-z(43)+z(76))

         vtors=vtors+(z(21)-z(28)-z(59)+z(70))**2
         ddz(21)=ddz(21)+(delt/bndprm**2)*(z(21)-z(28)-z(59)+z(70))
         ddz(28)=ddz(28)-(delt/bndprm**2)*(z(21)-z(28)-z(59)+z(70))
         ddz(59)=ddz(59)-(delt/bndprm**2)*(z(21)-z(28)-z(59)+z(70))
         ddz(70)=ddz(70)+(delt/bndprm**2)*(z(21)-z(28)-z(59)+z(70))

         vtors=vtors+(z(22)-z(54)-z(45)+z(77))**2
         ddz(22)=ddz(22)+(delt/bndprm**2)*(z(22)-z(54)-z(45)+z(77))
         ddz(54)=ddz(54)-(delt/bndprm**2)*(z(22)-z(54)-z(45)+z(77))
         ddz(45)=ddz(45)-(delt/bndprm**2)*(z(22)-z(54)-z(45)+z(77))
         ddz(77)=ddz(77)+(delt/bndprm**2)*(z(22)-z(54)-z(45)+z(77))

         vtors=vtors+(z(22)-z(56)-z(46)+z(78))**2
         ddz(22)=ddz(22)+(delt/bndprm**2)*(z(22)-z(56)-z(46)+z(78))
         ddz(56)=ddz(56)-(delt/bndprm**2)*(z(22)-z(56)-z(46)+z(78))
         ddz(46)=ddz(46)-(delt/bndprm**2)*(z(22)-z(56)-z(46)+z(78))
         ddz(78)=ddz(78)+(delt/bndprm**2)*(z(22)-z(56)-z(46)+z(78))

         vtors=vtors+(z(23)-z(26)-z(60)+z(69))**2
         ddz(23)=ddz(23)+(delt/bndprm**2)*(z(23)-z(26)-z(60)+z(69))
         ddz(26)=ddz(26)-(delt/bndprm**2)*(z(23)-z(26)-z(60)+z(69))
         ddz(60)=ddz(60)-(delt/bndprm**2)*(z(23)-z(26)-z(60)+z(69))
         ddz(69)=ddz(69)+(delt/bndprm**2)*(z(23)-z(26)-z(60)+z(69))

         vtors=vtors+(z(24)-z(30)-z(61)+z(71))**2
         ddz(24)=ddz(24)+(delt/bndprm**2)*(z(24)-z(30)-z(61)+z(71))
         ddz(30)=ddz(30)-(delt/bndprm**2)*(z(24)-z(30)-z(61)+z(71))
         ddz(61)=ddz(61)-(delt/bndprm**2)*(z(24)-z(30)-z(61)+z(71))
         ddz(71)=ddz(71)+(delt/bndprm**2)*(z(24)-z(30)-z(61)+z(71))

         vtors=vtors+(z(25)-z(8)-z(41)+z(33))**2
         ddz(25)=ddz(25)+(delt/bndprm**2)*(z(25)-z(8)-z(41)+z(33))
         ddz(8)=ddz(8)-(delt/bndprm**2)*(z(25)-z(8)-z(41)+z(33))
         ddz(41)=ddz(41)-(delt/bndprm**2)*(z(25)-z(8)-z(41)+z(33))
         ddz(33)=ddz(33)+(delt/bndprm**2)*(z(25)-z(8)-z(41)+z(33))

         vtors=vtors+(z(25)-z(10)-z(44)+z(37))**2
         ddz(25)=ddz(25)+(delt/bndprm**2)*(z(25)-z(10)-z(44)+z(37))
         ddz(10)=ddz(10)-(delt/bndprm**2)*(z(25)-z(10)-z(44)+z(37))
         ddz(44)=ddz(44)-(delt/bndprm**2)*(z(25)-z(10)-z(44)+z(37))
         ddz(37)=ddz(37)+(delt/bndprm**2)*(z(25)-z(10)-z(44)+z(37))

         vtors=vtors+(z(26)-z(6)-z(38)+z(14))**2
         ddz(26)=ddz(26)+(delt/bndprm**2)*(z(26)-z(6)-z(38)+z(14))
         ddz(6)=ddz(6)-(delt/bndprm**2)*(z(26)-z(6)-z(38)+z(14))
         ddz(38)=ddz(38)-(delt/bndprm**2)*(z(26)-z(6)-z(38)+z(14))
         ddz(14)=ddz(14)+(delt/bndprm**2)*(z(26)-z(6)-z(38)+z(14))

         vtors=vtors+(z(26)-z(48)-z(69)+z(87))**2
         ddz(26)=ddz(26)+(delt/bndprm**2)*(z(26)-z(48)-z(69)+z(87))
         ddz(48)=ddz(48)-(delt/bndprm**2)*(z(26)-z(48)-z(69)+z(87))
         ddz(69)=ddz(69)-(delt/bndprm**2)*(z(26)-z(48)-z(69)+z(87))
         ddz(87)=ddz(87)+(delt/bndprm**2)*(z(26)-z(48)-z(69)+z(87))

         vtors=vtors+(z(27)-z(21)-z(68)+z(59))**2
         ddz(27)=ddz(27)+(delt/bndprm**2)*(z(27)-z(21)-z(68)+z(59))
         ddz(21)=ddz(21)-(delt/bndprm**2)*(z(27)-z(21)-z(68)+z(59))
         ddz(68)=ddz(68)-(delt/bndprm**2)*(z(27)-z(21)-z(68)+z(59))
         ddz(59)=ddz(59)+(delt/bndprm**2)*(z(27)-z(21)-z(68)+z(59))

         vtors=vtors+(z(27)-z(62)-z(38)+z(75))**2
         ddz(27)=ddz(27)+(delt/bndprm**2)*(z(27)-z(62)-z(38)+z(75))
         ddz(62)=ddz(62)-(delt/bndprm**2)*(z(27)-z(62)-z(38)+z(75))
         ddz(38)=ddz(38)-(delt/bndprm**2)*(z(27)-z(62)-z(38)+z(75))
         ddz(75)=ddz(75)+(delt/bndprm**2)*(z(27)-z(62)-z(38)+z(75))

         vtors=vtors+(z(28)-z(56)-z(70)+z(90))**2
         ddz(28)=ddz(28)+(delt/bndprm**2)*(z(28)-z(56)-z(70)+z(90))
         ddz(56)=ddz(56)-(delt/bndprm**2)*(z(28)-z(56)-z(70)+z(90))
         ddz(70)=ddz(70)-(delt/bndprm**2)*(z(28)-z(56)-z(70)+z(90))
         ddz(90)=ddz(90)+(delt/bndprm**2)*(z(28)-z(56)-z(70)+z(90))

         vtors=vtors+(z(28)-z(64)-z(39)+z(77))**2
         ddz(28)=ddz(28)+(delt/bndprm**2)*(z(28)-z(64)-z(39)+z(77))
         ddz(64)=ddz(64)-(delt/bndprm**2)*(z(28)-z(64)-z(39)+z(77))
         ddz(39)=ddz(39)-(delt/bndprm**2)*(z(28)-z(64)-z(39)+z(77))
         ddz(77)=ddz(77)+(delt/bndprm**2)*(z(28)-z(64)-z(39)+z(77))

         vtors=vtors+(z(29)-z(8)-z(40)+z(17))**2
         ddz(29)=ddz(29)+(delt/bndprm**2)*(z(29)-z(8)-z(40)+z(17))
         ddz(8)=ddz(8)-(delt/bndprm**2)*(z(29)-z(8)-z(40)+z(17))
         ddz(40)=ddz(40)-(delt/bndprm**2)*(z(29)-z(8)-z(40)+z(17))
         ddz(17)=ddz(17)+(delt/bndprm**2)*(z(29)-z(8)-z(40)+z(17))

         vtors=vtors+(z(29)-z(23)-z(73)+z(60))**2
         ddz(29)=ddz(29)+(delt/bndprm**2)*(z(29)-z(23)-z(73)+z(60))
         ddz(23)=ddz(23)-(delt/bndprm**2)*(z(29)-z(23)-z(73)+z(60))
         ddz(73)=ddz(73)-(delt/bndprm**2)*(z(29)-z(23)-z(73)+z(60))
         ddz(60)=ddz(60)+(delt/bndprm**2)*(z(29)-z(23)-z(73)+z(60))

         vtors=vtors+(z(30)-z(54)-z(71)+z(88))**2
         ddz(30)=ddz(30)+(delt/bndprm**2)*(z(30)-z(54)-z(71)+z(88))
         ddz(54)=ddz(54)-(delt/bndprm**2)*(z(30)-z(54)-z(71)+z(88))
         ddz(71)=ddz(71)-(delt/bndprm**2)*(z(30)-z(54)-z(71)+z(88))
         ddz(88)=ddz(88)+(delt/bndprm**2)*(z(30)-z(54)-z(71)+z(88))

         vtors=vtors+(z(30)-z(65)-z(39)+z(78))**2
         ddz(30)=ddz(30)+(delt/bndprm**2)*(z(30)-z(65)-z(39)+z(78))
         ddz(65)=ddz(65)-(delt/bndprm**2)*(z(30)-z(65)-z(39)+z(78))
         ddz(39)=ddz(39)-(delt/bndprm**2)*(z(30)-z(65)-z(39)+z(78))
         ddz(78)=ddz(78)+(delt/bndprm**2)*(z(30)-z(65)-z(39)+z(78))

         vtors=vtors+(z(31)-z(10)-z(40)+z(19))**2
         ddz(31)=ddz(31)+(delt/bndprm**2)*(z(31)-z(10)-z(40)+z(19))
         ddz(10)=ddz(10)-(delt/bndprm**2)*(z(31)-z(10)-z(40)+z(19))
         ddz(40)=ddz(40)-(delt/bndprm**2)*(z(31)-z(10)-z(40)+z(19))
         ddz(19)=ddz(19)+(delt/bndprm**2)*(z(31)-z(10)-z(40)+z(19))

         vtors=vtors+(z(31)-z(24)-z(72)+z(61))**2
         ddz(31)=ddz(31)+(delt/bndprm**2)*(z(31)-z(24)-z(72)+z(61))
         ddz(24)=ddz(24)-(delt/bndprm**2)*(z(31)-z(24)-z(72)+z(61))
         ddz(72)=ddz(72)-(delt/bndprm**2)*(z(31)-z(24)-z(72)+z(61))
         ddz(61)=ddz(61)+(delt/bndprm**2)*(z(31)-z(24)-z(72)+z(61))

         vtors=vtors+(z(32)-z(59)-z(64)+z(92))**2
         ddz(32)=ddz(32)+(delt/bndprm**2)*(z(32)-z(59)-z(64)+z(92))
         ddz(59)=ddz(59)-(delt/bndprm**2)*(z(32)-z(59)-z(64)+z(92))
         ddz(64)=ddz(64)-(delt/bndprm**2)*(z(32)-z(59)-z(64)+z(92))
         ddz(92)=ddz(92)+(delt/bndprm**2)*(z(32)-z(59)-z(64)+z(92))

         vtors=vtors+(z(32)-z(70)-z(52)+z(81))**2
         ddz(32)=ddz(32)+(delt/bndprm**2)*(z(32)-z(70)-z(52)+z(81))
         ddz(70)=ddz(70)-(delt/bndprm**2)*(z(32)-z(70)-z(52)+z(81))
         ddz(52)=ddz(52)-(delt/bndprm**2)*(z(32)-z(70)-z(52)+z(81))
         ddz(81)=ddz(81)+(delt/bndprm**2)*(z(32)-z(70)-z(52)+z(81))

         vtors=vtors+(z(33)-z(12)-z(53)+z(26))**2
         ddz(33)=ddz(33)+(delt/bndprm**2)*(z(33)-z(12)-z(53)+z(26))
         ddz(12)=ddz(12)-(delt/bndprm**2)*(z(33)-z(12)-z(53)+z(26))
         ddz(53)=ddz(53)-(delt/bndprm**2)*(z(33)-z(12)-z(53)+z(26))
         ddz(26)=ddz(26)+(delt/bndprm**2)*(z(33)-z(12)-z(53)+z(26))

         vtors=vtors+(z(33)-z(19)-z(67)+z(50))**2
         ddz(33)=ddz(33)+(delt/bndprm**2)*(z(33)-z(19)-z(67)+z(50))
         ddz(19)=ddz(19)-(delt/bndprm**2)*(z(33)-z(19)-z(67)+z(50))
         ddz(67)=ddz(67)-(delt/bndprm**2)*(z(33)-z(19)-z(67)+z(50))
         ddz(50)=ddz(50)+(delt/bndprm**2)*(z(33)-z(19)-z(67)+z(50))

         vtors=vtors+(z(34)-z(59)-z(62)+z(94))**2
         ddz(34)=ddz(34)+(delt/bndprm**2)*(z(34)-z(59)-z(62)+z(94))
         ddz(59)=ddz(59)-(delt/bndprm**2)*(z(34)-z(59)-z(62)+z(94))
         ddz(62)=ddz(62)-(delt/bndprm**2)*(z(34)-z(59)-z(62)+z(94))
         ddz(94)=ddz(94)+(delt/bndprm**2)*(z(34)-z(59)-z(62)+z(94))

         vtors=vtors+(z(34)-z(68)-z(51)+z(81))**2
         ddz(34)=ddz(34)+(delt/bndprm**2)*(z(34)-z(68)-z(51)+z(81))
         ddz(68)=ddz(68)-(delt/bndprm**2)*(z(34)-z(68)-z(51)+z(81))
         ddz(51)=ddz(51)-(delt/bndprm**2)*(z(34)-z(68)-z(51)+z(81))
         ddz(81)=ddz(81)+(delt/bndprm**2)*(z(34)-z(68)-z(51)+z(81))

         vtors=vtors+(z(35)-z(12)-z(49)+z(29))**2
         ddz(35)=ddz(35)+(delt/bndprm**2)*(z(35)-z(12)-z(49)+z(29))
         ddz(12)=ddz(12)-(delt/bndprm**2)*(z(35)-z(12)-z(49)+z(29))
         ddz(49)=ddz(49)-(delt/bndprm**2)*(z(35)-z(12)-z(49)+z(29))
         ddz(29)=ddz(29)+(delt/bndprm**2)*(z(35)-z(12)-z(49)+z(29))

         vtors=vtors+(z(35)-z(15)-z(63)+z(48))**2
         ddz(35)=ddz(35)+(delt/bndprm**2)*(z(35)-z(15)-z(63)+z(48))
         ddz(15)=ddz(15)-(delt/bndprm**2)*(z(35)-z(15)-z(63)+z(48))
         ddz(63)=ddz(63)-(delt/bndprm**2)*(z(35)-z(15)-z(63)+z(48))
         ddz(48)=ddz(48)+(delt/bndprm**2)*(z(35)-z(15)-z(63)+z(48))

         vtors=vtors+(z(36)-z(18)-z(65)+z(54))**2
         ddz(36)=ddz(36)+(delt/bndprm**2)*(z(36)-z(18)-z(65)+z(54))
         ddz(18)=ddz(18)-(delt/bndprm**2)*(z(36)-z(18)-z(65)+z(54))
         ddz(65)=ddz(65)-(delt/bndprm**2)*(z(36)-z(18)-z(65)+z(54))
         ddz(54)=ddz(54)+(delt/bndprm**2)*(z(36)-z(18)-z(65)+z(54))

         vtors=vtors+(z(36)-z(71)-z(58)+z(85))**2
         ddz(36)=ddz(36)+(delt/bndprm**2)*(z(36)-z(71)-z(58)+z(85))
         ddz(71)=ddz(71)-(delt/bndprm**2)*(z(36)-z(71)-z(58)+z(85))
         ddz(58)=ddz(58)-(delt/bndprm**2)*(z(36)-z(71)-z(58)+z(85))
         ddz(85)=ddz(85)+(delt/bndprm**2)*(z(36)-z(71)-z(58)+z(85))

         vtors=vtors+(z(37)-z(13)-z(57)+z(30))**2
         ddz(37)=ddz(37)+(delt/bndprm**2)*(z(37)-z(13)-z(57)+z(30))
         ddz(13)=ddz(13)-(delt/bndprm**2)*(z(37)-z(13)-z(57)+z(30))
         ddz(57)=ddz(57)-(delt/bndprm**2)*(z(37)-z(13)-z(57)+z(30))
         ddz(30)=ddz(30)+(delt/bndprm**2)*(z(37)-z(13)-z(57)+z(30))

         vtors=vtors+(z(37)-z(61)-z(66)+z(103))**2
         ddz(37)=ddz(37)+(delt/bndprm**2)*(z(37)-z(61)-z(66)+z(103))
         ddz(61)=ddz(61)-(delt/bndprm**2)*(z(37)-z(61)-z(66)+z(103))
         ddz(66)=ddz(66)-(delt/bndprm**2)*(z(37)-z(61)-z(66)+z(103))
         ddz(103)=ddz(103)+(delt/bndprm**2)*(z(37)-z(61)-z(66)+z(103))

         vtors=vtors+(z(38)-z(42)-z(82)+z(96))**2
         ddz(38)=ddz(38)+(delt/bndprm**2)*(z(38)-z(42)-z(82)+z(96))
         ddz(42)=ddz(42)-(delt/bndprm**2)*(z(38)-z(42)-z(82)+z(96))
         ddz(82)=ddz(82)-(delt/bndprm**2)*(z(38)-z(42)-z(82)+z(96))
         ddz(96)=ddz(96)+(delt/bndprm**2)*(z(38)-z(42)-z(82)+z(96))

         vtors=vtors+(z(39)-z(45)-z(83)+z(100))**2
         ddz(39)=ddz(39)+(delt/bndprm**2)*(z(39)-z(45)-z(83)+z(100))
         ddz(45)=ddz(45)-(delt/bndprm**2)*(z(39)-z(45)-z(83)+z(100))
         ddz(83)=ddz(83)-(delt/bndprm**2)*(z(39)-z(45)-z(83)+z(100))
         ddz(100)=ddz(100)+(delt/bndprm**2)*(z(39)-z(45)-z(83)+z(100))

         vtors=vtors+(z(40)-z(44)-z(84)+z(99))**2
         ddz(40)=ddz(40)+(delt/bndprm**2)*(z(40)-z(44)-z(84)+z(99))
         ddz(44)=ddz(44)-(delt/bndprm**2)*(z(40)-z(44)-z(84)+z(99))
         ddz(84)=ddz(84)-(delt/bndprm**2)*(z(40)-z(44)-z(84)+z(99))
         ddz(99)=ddz(99)+(delt/bndprm**2)*(z(40)-z(44)-z(84)+z(99))

         vtors=vtors+(z(41)-z(19)-z(49)+z(23))**2
         ddz(41)=ddz(41)+(delt/bndprm**2)*(z(41)-z(19)-z(49)+z(23))
         ddz(19)=ddz(19)-(delt/bndprm**2)*(z(41)-z(19)-z(49)+z(23))
         ddz(49)=ddz(49)-(delt/bndprm**2)*(z(41)-z(19)-z(49)+z(23))
         ddz(23)=ddz(23)+(delt/bndprm**2)*(z(41)-z(19)-z(49)+z(23))

         vtors=vtors+(z(41)-z(40)-z(101)+z(84))**2
         ddz(41)=ddz(41)+(delt/bndprm**2)*(z(41)-z(40)-z(101)+z(84))
         ddz(40)=ddz(40)-(delt/bndprm**2)*(z(41)-z(40)-z(101)+z(84))
         ddz(101)=ddz(101)-(delt/bndprm**2)*(z(41)-z(40)-z(101)+z(84))
         ddz(84)=ddz(84)+(delt/bndprm**2)*(z(41)-z(40)-z(101)+z(84))

         vtors=vtors+(z(42)-z(68)-z(96)+z(112))**2
         ddz(42)=ddz(42)+(delt/bndprm**2)*(z(42)-z(68)-z(96)+z(112))
         ddz(68)=ddz(68)-(delt/bndprm**2)*(z(42)-z(68)-z(96)+z(112))
         ddz(96)=ddz(96)-(delt/bndprm**2)*(z(42)-z(68)-z(96)+z(112))
         ddz(112)=ddz(112)+(delt/bndprm**2)*(z(42)-z(68)-z(96)+z(112))

         vtors=vtors+(z(42)-z(86)-z(52)+z(94))**2
         ddz(42)=ddz(42)+(delt/bndprm**2)*(z(42)-z(86)-z(52)+z(94))
         ddz(86)=ddz(86)-(delt/bndprm**2)*(z(42)-z(86)-z(52)+z(94))
         ddz(52)=ddz(52)-(delt/bndprm**2)*(z(42)-z(86)-z(52)+z(94))
         ddz(94)=ddz(94)+(delt/bndprm**2)*(z(42)-z(86)-z(52)+z(94))

         vtors=vtors+(z(43)-z(15)-z(53)+z(23))**2
         ddz(43)=ddz(43)+(delt/bndprm**2)*(z(43)-z(15)-z(53)+z(23))
         ddz(15)=ddz(15)-(delt/bndprm**2)*(z(43)-z(15)-z(53)+z(23))
         ddz(53)=ddz(53)-(delt/bndprm**2)*(z(43)-z(15)-z(53)+z(23))
         ddz(23)=ddz(23)+(delt/bndprm**2)*(z(43)-z(15)-z(53)+z(23))

         vtors=vtors+(z(43)-z(38)-z(97)+z(82))**2
         ddz(43)=ddz(43)+(delt/bndprm**2)*(z(43)-z(38)-z(97)+z(82))
         ddz(38)=ddz(38)-(delt/bndprm**2)*(z(43)-z(38)-z(97)+z(82))
         ddz(97)=ddz(97)-(delt/bndprm**2)*(z(43)-z(38)-z(97)+z(82))
         ddz(82)=ddz(82)+(delt/bndprm**2)*(z(43)-z(38)-z(97)+z(82))

         vtors=vtors+(z(44)-z(17)-z(58)+z(24))**2
         ddz(44)=ddz(44)+(delt/bndprm**2)*(z(44)-z(17)-z(58)+z(24))
         ddz(17)=ddz(17)-(delt/bndprm**2)*(z(44)-z(17)-z(58)+z(24))
         ddz(58)=ddz(58)-(delt/bndprm**2)*(z(44)-z(17)-z(58)+z(24))
         ddz(24)=ddz(24)+(delt/bndprm**2)*(z(44)-z(17)-z(58)+z(24))

         vtors=vtors+(z(44)-z(72)-z(99)+z(115))**2
         ddz(44)=ddz(44)+(delt/bndprm**2)*(z(44)-z(72)-z(99)+z(115))
         ddz(72)=ddz(72)-(delt/bndprm**2)*(z(44)-z(72)-z(99)+z(115))
         ddz(99)=ddz(99)-(delt/bndprm**2)*(z(44)-z(72)-z(99)+z(115))
         ddz(115)=ddz(115)+(delt/bndprm**2)*(z(44)-z(72)-z(99)+z(115))

         vtors=vtors+(z(45)-z(70)-z(100)+z(110))**2
         ddz(45)=ddz(45)+(delt/bndprm**2)*(z(45)-z(70)-z(100)+z(110))
         ddz(70)=ddz(70)-(delt/bndprm**2)*(z(45)-z(70)-z(100)+z(110))
         ddz(100)=ddz(100)-(delt/bndprm**2)*(z(45)-z(70)-z(100)+z(110))
         ddz(110)=ddz(110)+(delt/bndprm**2)*(z(45)-z(70)-z(100)+z(110))

         vtors=vtors+(z(45)-z(90)-z(51)+z(92))**2
         ddz(45)=ddz(45)+(delt/bndprm**2)*(z(45)-z(90)-z(51)+z(92))
         ddz(90)=ddz(90)-(delt/bndprm**2)*(z(45)-z(90)-z(51)+z(92))
         ddz(51)=ddz(51)-(delt/bndprm**2)*(z(45)-z(90)-z(51)+z(92))
         ddz(92)=ddz(92)+(delt/bndprm**2)*(z(45)-z(90)-z(51)+z(92))

         vtors=vtors+(z(46)-z(39)-z(98)+z(83))**2
         ddz(46)=ddz(46)+(delt/bndprm**2)*(z(46)-z(39)-z(98)+z(83))
         ddz(39)=ddz(39)-(delt/bndprm**2)*(z(46)-z(39)-z(98)+z(83))
         ddz(98)=ddz(98)-(delt/bndprm**2)*(z(46)-z(39)-z(98)+z(83))
         ddz(83)=ddz(83)+(delt/bndprm**2)*(z(46)-z(39)-z(98)+z(83))

         vtors=vtors+(z(46)-z(88)-z(57)+z(102))**2
         ddz(46)=ddz(46)+(delt/bndprm**2)*(z(46)-z(88)-z(57)+z(102))
         ddz(88)=ddz(88)-(delt/bndprm**2)*(z(46)-z(88)-z(57)+z(102))
         ddz(57)=ddz(57)-(delt/bndprm**2)*(z(46)-z(88)-z(57)+z(102))
         ddz(102)=ddz(102)+(delt/bndprm**2)*(z(46)-z(88)-z(57)+z(102))

         vtors=vtors+(z(47)-z(27)-z(86)+z(68))**2
         ddz(47)=ddz(47)+(delt/bndprm**2)*(z(47)-z(27)-z(86)+z(68))
         ddz(27)=ddz(27)-(delt/bndprm**2)*(z(47)-z(27)-z(86)+z(68))
         ddz(86)=ddz(86)-(delt/bndprm**2)*(z(47)-z(27)-z(86)+z(68))
         ddz(68)=ddz(68)+(delt/bndprm**2)*(z(47)-z(27)-z(86)+z(68))

         vtors=vtors+(z(48)-z(20)-z(75)+z(42))**2
         ddz(48)=ddz(48)+(delt/bndprm**2)*(z(48)-z(20)-z(75)+z(42))
         ddz(20)=ddz(20)-(delt/bndprm**2)*(z(48)-z(20)-z(75)+z(42))
         ddz(75)=ddz(75)-(delt/bndprm**2)*(z(48)-z(20)-z(75)+z(42))
         ddz(42)=ddz(42)+(delt/bndprm**2)*(z(48)-z(20)-z(75)+z(42))

         vtors=vtors+(z(49)-z(53)-z(105)+z(107))**2
         ddz(49)=ddz(49)+(delt/bndprm**2)*(z(49)-z(53)-z(105)+z(107))
         ddz(53)=ddz(53)-(delt/bndprm**2)*(z(49)-z(53)-z(105)+z(107))
         ddz(105)=ddz(105)-(delt/bndprm**2)*(z(49)-z(53)-z(105)+z(107))
         ddz(107)=ddz(107)+(delt/bndprm**2)*(z(49)-z(53)-z(105)+z(107))

         vtors=vtors+(z(50)-z(25)-z(79)+z(44))**2
         ddz(50)=ddz(50)+(delt/bndprm**2)*(z(50)-z(25)-z(79)+z(44))
         ddz(25)=ddz(25)-(delt/bndprm**2)*(z(50)-z(25)-z(79)+z(44))
         ddz(79)=ddz(79)-(delt/bndprm**2)*(z(50)-z(25)-z(79)+z(44))
         ddz(44)=ddz(44)+(delt/bndprm**2)*(z(50)-z(25)-z(79)+z(44))

         vtors=vtors+(z(50)-z(29)-z(91)+z(73))**2
         ddz(50)=ddz(50)+(delt/bndprm**2)*(z(50)-z(29)-z(91)+z(73))
         ddz(29)=ddz(29)-(delt/bndprm**2)*(z(50)-z(29)-z(91)+z(73))
         ddz(91)=ddz(91)-(delt/bndprm**2)*(z(50)-z(29)-z(91)+z(73))
         ddz(73)=ddz(73)+(delt/bndprm**2)*(z(50)-z(29)-z(91)+z(73))

         vtors=vtors+(z(51)-z(64)-z(104)+z(120))**2
         ddz(51)=ddz(51)+(delt/bndprm**2)*(z(51)-z(64)-z(104)+z(120))
         ddz(64)=ddz(64)-(delt/bndprm**2)*(z(51)-z(64)-z(104)+z(120))
         ddz(104)=ddz(104)-(delt/bndprm**2)*(z(51)-z(64)-z(104)+z(120))
         ddz(120)=ddz(120)+(delt/bndprm**2)*(z(51)-z(64)-z(104)+z(120))

         vtors=vtors+(z(52)-z(51)-z(106)+z(104))**2
         ddz(52)=ddz(52)+(delt/bndprm**2)*(z(52)-z(51)-z(106)+z(104))
         ddz(51)=ddz(51)-(delt/bndprm**2)*(z(52)-z(51)-z(106)+z(104))
         ddz(106)=ddz(106)-(delt/bndprm**2)*(z(52)-z(51)-z(106)+z(104))
         ddz(104)=ddz(104)+(delt/bndprm**2)*(z(52)-z(51)-z(106)+z(104))

         vtors=vtors+(z(53)-z(63)-z(107)+z(117))**2
         ddz(53)=ddz(53)+(delt/bndprm**2)*(z(53)-z(63)-z(107)+z(117))
         ddz(63)=ddz(63)-(delt/bndprm**2)*(z(53)-z(63)-z(107)+z(117))
         ddz(107)=ddz(107)-(delt/bndprm**2)*(z(53)-z(63)-z(107)+z(117))
         ddz(117)=ddz(117)+(delt/bndprm**2)*(z(53)-z(63)-z(107)+z(117))

         vtors=vtors+(z(55)-z(25)-z(74)+z(41))**2
         ddz(55)=ddz(55)+(delt/bndprm**2)*(z(55)-z(25)-z(74)+z(41))
         ddz(25)=ddz(25)-(delt/bndprm**2)*(z(55)-z(25)-z(74)+z(41))
         ddz(74)=ddz(74)-(delt/bndprm**2)*(z(55)-z(25)-z(74)+z(41))
         ddz(41)=ddz(41)+(delt/bndprm**2)*(z(55)-z(25)-z(74)+z(41))

         vtors=vtors+(z(55)-z(31)-z(89)+z(72))**2
         ddz(55)=ddz(55)+(delt/bndprm**2)*(z(55)-z(31)-z(89)+z(72))
         ddz(31)=ddz(31)-(delt/bndprm**2)*(z(55)-z(31)-z(89)+z(72))
         ddz(89)=ddz(89)-(delt/bndprm**2)*(z(55)-z(31)-z(89)+z(72))
         ddz(72)=ddz(72)+(delt/bndprm**2)*(z(55)-z(31)-z(89)+z(72))

         vtors=vtors+(z(57)-z(65)-z(108)+z(118))**2
         ddz(57)=ddz(57)+(delt/bndprm**2)*(z(57)-z(65)-z(108)+z(118))
         ddz(65)=ddz(65)-(delt/bndprm**2)*(z(57)-z(65)-z(108)+z(118))
         ddz(108)=ddz(108)-(delt/bndprm**2)*(z(57)-z(65)-z(108)+z(118))
         ddz(118)=ddz(118)+(delt/bndprm**2)*(z(57)-z(65)-z(108)+z(118))

         vtors=vtors+(z(58)-z(57)-z(109)+z(108))**2
         ddz(58)=ddz(58)+(delt/bndprm**2)*(z(58)-z(57)-z(109)+z(108))
         ddz(57)=ddz(57)-(delt/bndprm**2)*(z(58)-z(57)-z(109)+z(108))
         ddz(109)=ddz(109)-(delt/bndprm**2)*(z(58)-z(57)-z(109)+z(108))
         ddz(108)=ddz(108)+(delt/bndprm**2)*(z(58)-z(57)-z(109)+z(108))

         vtors=vtors+(z(60)-z(33)-z(93)+z(67))**2
         ddz(60)=ddz(60)+(delt/bndprm**2)*(z(60)-z(33)-z(93)+z(67))
         ddz(33)=ddz(33)-(delt/bndprm**2)*(z(60)-z(33)-z(93)+z(67))
         ddz(93)=ddz(93)-(delt/bndprm**2)*(z(60)-z(33)-z(93)+z(67))
         ddz(67)=ddz(67)+(delt/bndprm**2)*(z(60)-z(33)-z(93)+z(67))

         vtors=vtors+(z(60)-z(35)-z(95)+z(63))**2
         ddz(60)=ddz(60)+(delt/bndprm**2)*(z(60)-z(35)-z(95)+z(63))
         ddz(35)=ddz(35)-(delt/bndprm**2)*(z(60)-z(35)-z(95)+z(63))
         ddz(95)=ddz(95)-(delt/bndprm**2)*(z(60)-z(35)-z(95)+z(63))
         ddz(63)=ddz(63)+(delt/bndprm**2)*(z(60)-z(35)-z(95)+z(63))

         vtors=vtors+(z(61)-z(36)-z(102)+z(65))**2
         ddz(61)=ddz(61)+(delt/bndprm**2)*(z(61)-z(36)-z(102)+z(65))
         ddz(36)=ddz(36)-(delt/bndprm**2)*(z(61)-z(36)-z(102)+z(65))
         ddz(102)=ddz(102)-(delt/bndprm**2)*(z(61)-z(36)-z(102)+z(65))
         ddz(65)=ddz(65)+(delt/bndprm**2)*(z(61)-z(36)-z(102)+z(65))

         vtors=vtors+(z(62)-z(52)-z(116)+z(106))**2
         ddz(62)=ddz(62)+(delt/bndprm**2)*(z(62)-z(52)-z(116)+z(106))
         ddz(52)=ddz(52)-(delt/bndprm**2)*(z(62)-z(52)-z(116)+z(106))
         ddz(116)=ddz(116)-(delt/bndprm**2)*(z(62)-z(52)-z(116)+z(106))
         ddz(106)=ddz(106)+(delt/bndprm**2)*(z(62)-z(52)-z(116)+z(106))

         vtors=vtors+(z(63)-z(26)-z(76)+z(38))**2
         ddz(63)=ddz(63)+(delt/bndprm**2)*(z(63)-z(26)-z(76)+z(38))
         ddz(26)=ddz(26)-(delt/bndprm**2)*(z(63)-z(26)-z(76)+z(38))
         ddz(76)=ddz(76)-(delt/bndprm**2)*(z(63)-z(26)-z(76)+z(38))
         ddz(38)=ddz(38)+(delt/bndprm**2)*(z(63)-z(26)-z(76)+z(38))

         vtors=vtors+(z(66)-z(31)-z(79)+z(40))**2
         ddz(66)=ddz(66)+(delt/bndprm**2)*(z(66)-z(31)-z(79)+z(40))
         ddz(31)=ddz(31)-(delt/bndprm**2)*(z(66)-z(31)-z(79)+z(40))
         ddz(79)=ddz(79)-(delt/bndprm**2)*(z(66)-z(31)-z(79)+z(40))
         ddz(40)=ddz(40)+(delt/bndprm**2)*(z(66)-z(31)-z(79)+z(40))

         vtors=vtors+(z(66)-z(58)-z(119)+z(109))**2
         ddz(66)=ddz(66)+(delt/bndprm**2)*(z(66)-z(58)-z(119)+z(109))
         ddz(58)=ddz(58)-(delt/bndprm**2)*(z(66)-z(58)-z(119)+z(109))
         ddz(119)=ddz(119)-(delt/bndprm**2)*(z(66)-z(58)-z(119)+z(109))
         ddz(109)=ddz(109)+(delt/bndprm**2)*(z(66)-z(58)-z(119)+z(109))

         vtors=vtors+(z(67)-z(29)-z(74)+z(40))**2
         ddz(67)=ddz(67)+(delt/bndprm**2)*(z(67)-z(29)-z(74)+z(40))
         ddz(29)=ddz(29)-(delt/bndprm**2)*(z(67)-z(29)-z(74)+z(40))
         ddz(74)=ddz(74)-(delt/bndprm**2)*(z(67)-z(29)-z(74)+z(40))
         ddz(40)=ddz(40)+(delt/bndprm**2)*(z(67)-z(29)-z(74)+z(40))

         vtors=vtors+(z(67)-z(49)-z(121)+z(105))**2
         ddz(67)=ddz(67)+(delt/bndprm**2)*(z(67)-z(49)-z(121)+z(105))
         ddz(49)=ddz(49)-(delt/bndprm**2)*(z(67)-z(49)-z(121)+z(105))
         ddz(121)=ddz(121)-(delt/bndprm**2)*(z(67)-z(49)-z(121)+z(105))
         ddz(105)=ddz(105)+(delt/bndprm**2)*(z(67)-z(49)-z(121)+z(105))

         vtors=vtors+(z(69)-z(35)-z(80)+z(49))**2
         ddz(69)=ddz(69)+(delt/bndprm**2)*(z(69)-z(35)-z(80)+z(49))
         ddz(35)=ddz(35)-(delt/bndprm**2)*(z(69)-z(35)-z(80)+z(49))
         ddz(80)=ddz(80)-(delt/bndprm**2)*(z(69)-z(35)-z(80)+z(49))
         ddz(49)=ddz(49)+(delt/bndprm**2)*(z(69)-z(35)-z(80)+z(49))

         vtors=vtors+(z(69)-z(43)-z(113)+z(97))**2
         ddz(69)=ddz(69)+(delt/bndprm**2)*(z(69)-z(43)-z(113)+z(97))
         ddz(43)=ddz(43)-(delt/bndprm**2)*(z(69)-z(43)-z(113)+z(97))
         ddz(113)=ddz(113)-(delt/bndprm**2)*(z(69)-z(43)-z(113)+z(97))
         ddz(97)=ddz(97)+(delt/bndprm**2)*(z(69)-z(43)-z(113)+z(97))

         vtors=vtors+(z(71)-z(46)-z(114)+z(98))**2
         ddz(71)=ddz(71)+(delt/bndprm**2)*(z(71)-z(46)-z(114)+z(98))
         ddz(46)=ddz(46)-(delt/bndprm**2)*(z(71)-z(46)-z(114)+z(98))
         ddz(114)=ddz(114)-(delt/bndprm**2)*(z(71)-z(46)-z(114)+z(98))
         ddz(98)=ddz(98)+(delt/bndprm**2)*(z(71)-z(46)-z(114)+z(98))

         vtors=vtors+(z(72)-z(37)-z(85)+z(57))**2
         ddz(72)=ddz(72)+(delt/bndprm**2)*(z(72)-z(37)-z(85)+z(57))
         ddz(37)=ddz(37)-(delt/bndprm**2)*(z(72)-z(37)-z(85)+z(57))
         ddz(85)=ddz(85)-(delt/bndprm**2)*(z(72)-z(37)-z(85)+z(57))
         ddz(57)=ddz(57)+(delt/bndprm**2)*(z(72)-z(37)-z(85)+z(57))

         vtors=vtors+(z(73)-z(33)-z(80)+z(53))**2
         ddz(73)=ddz(73)+(delt/bndprm**2)*(z(73)-z(33)-z(80)+z(53))
         ddz(33)=ddz(33)-(delt/bndprm**2)*(z(73)-z(33)-z(80)+z(53))
         ddz(80)=ddz(80)-(delt/bndprm**2)*(z(73)-z(33)-z(80)+z(53))
         ddz(53)=ddz(53)+(delt/bndprm**2)*(z(73)-z(33)-z(80)+z(53))

         vtors=vtors+(z(73)-z(41)-z(111)+z(101))**2
         ddz(73)=ddz(73)+(delt/bndprm**2)*(z(73)-z(41)-z(111)+z(101))
         ddz(41)=ddz(41)-(delt/bndprm**2)*(z(73)-z(41)-z(111)+z(101))
         ddz(111)=ddz(111)-(delt/bndprm**2)*(z(73)-z(41)-z(111)+z(101))
         ddz(101)=ddz(101)+(delt/bndprm**2)*(z(73)-z(41)-z(111)+z(101))

         vtors=vtors+(z(87)-z(43)-z(95)+z(53))**2
         ddz(87)=ddz(87)+(delt/bndprm**2)*(z(87)-z(43)-z(95)+z(53))
         ddz(43)=ddz(43)-(delt/bndprm**2)*(z(87)-z(43)-z(95)+z(53))
         ddz(95)=ddz(95)-(delt/bndprm**2)*(z(87)-z(43)-z(95)+z(53))
         ddz(53)=ddz(53)+(delt/bndprm**2)*(z(87)-z(43)-z(95)+z(53))

         vtors=vtors+(z(89)-z(44)-z(103)+z(58))**2
         ddz(89)=ddz(89)+(delt/bndprm**2)*(z(89)-z(44)-z(103)+z(58))
         ddz(44)=ddz(44)-(delt/bndprm**2)*(z(89)-z(44)-z(103)+z(58))
         ddz(103)=ddz(103)-(delt/bndprm**2)*(z(89)-z(44)-z(103)+z(58))
         ddz(58)=ddz(58)+(delt/bndprm**2)*(z(89)-z(44)-z(103)+z(58))

         vtors=vtors+(z(91)-z(41)-z(93)+z(49))**2
         ddz(91)=ddz(91)+(delt/bndprm**2)*(z(91)-z(41)-z(93)+z(49))
         ddz(41)=ddz(41)-(delt/bndprm**2)*(z(91)-z(41)-z(93)+z(49))
         ddz(93)=ddz(93)-(delt/bndprm**2)*(z(91)-z(41)-z(93)+z(49))
         ddz(49)=ddz(49)+(delt/bndprm**2)*(z(91)-z(41)-z(93)+z(49))

         vtors=vtors*0.5d0*delt/bndprm**2

         vpuck=vpuck+(z(2)+z(3)+z(4)-3.0d0*z(1))**2
         ddz(1)=ddz(1)-(gam2*9.0d0/bndprm**2)*(z(2)+z(3)+z(4)-3.0d0*z(1))
         ddz(2)=ddz(2)+(gam2*3.0d0/bndprm**2)*(z(2)+z(3)+z(4)-3.0d0*z(1))
         ddz(3)=ddz(3)+(gam2*3.0d0/bndprm**2)*(z(2)+z(3)+z(4)-3.0d0*z(1))
         ddz(4)=ddz(4)+(gam2*3.0d0/bndprm**2)*(z(2)+z(3)+z(4)-3.0d0*z(1))

         vpuck=vpuck+(z(1)+z(5)+z(6)-3.0d0*z(2))**2
         ddz(2)=ddz(2)-(gam2*9.0d0/bndprm**2)*(z(1)+z(5)+z(6)-3.0d0*z(2))
         ddz(1)=ddz(1)+(gam2*3.0d0/bndprm**2)*(z(1)+z(5)+z(6)-3.0d0*z(2))
         ddz(5)=ddz(5)+(gam2*3.0d0/bndprm**2)*(z(1)+z(5)+z(6)-3.0d0*z(2))
         ddz(6)=ddz(6)+(gam2*3.0d0/bndprm**2)*(z(1)+z(5)+z(6)-3.0d0*z(2))

         vpuck=vpuck+(z(1)+z(7)+z(9)-3.0d0*z(3))**2
         ddz(3)=ddz(3)-(gam2*9.0d0/bndprm**2)*(z(1)+z(7)+z(9)-3.0d0*z(3))
         ddz(1)=ddz(1)+(gam2*3.0d0/bndprm**2)*(z(1)+z(7)+z(9)-3.0d0*z(3))
         ddz(7)=ddz(7)+(gam2*3.0d0/bndprm**2)*(z(1)+z(7)+z(9)-3.0d0*z(3))
         ddz(9)=ddz(9)+(gam2*3.0d0/bndprm**2)*(z(1)+z(7)+z(9)-3.0d0*z(3))

         vpuck=vpuck+(z(1)+z(8)+z(10)-3.0d0*z(4))**2
         ddz(4)=ddz(4)-(gam2*9.0d0/bndprm**2)*(z(1)+z(8)+z(10)-3.0d0*z(4))
         ddz(1)=ddz(1)+(gam2*3.0d0/bndprm**2)*(z(1)+z(8)+z(10)-3.0d0*z(4))
         ddz(8)=ddz(8)+(gam2*3.0d0/bndprm**2)*(z(1)+z(8)+z(10)-3.0d0*z(4))
         ddz(10)=ddz(10)+(gam2*3.0d0/bndprm**2)*(z(1)+z(8)+z(10)-3.0d0*z(4))

         vpuck=vpuck+(z(2)+z(11)+z(14)-3.0d0*z(5))**2
         ddz(5)=ddz(5)-(gam2*9.0d0/bndprm**2)*(z(2)+z(11)+z(14)-3.0d0*z(5))
         ddz(2)=ddz(2)+(gam2*3.0d0/bndprm**2)*(z(2)+z(11)+z(14)-3.0d0*z(5))
         ddz(11)=ddz(11)+(gam2*3.0d0/bndprm**2)*(z(2)+z(11)+z(14)-3.0d0*z(5))
         ddz(14)=ddz(14)+(gam2*3.0d0/bndprm**2)*(z(2)+z(11)+z(14)-3.0d0*z(5))

         vpuck=vpuck+(z(2)+z(12)+z(15)-3.0d0*z(6))**2
         ddz(6)=ddz(6)-(gam2*9.0d0/bndprm**2)*(z(2)+z(12)+z(15)-3.0d0*z(6))
         ddz(2)=ddz(2)+(gam2*3.0d0/bndprm**2)*(z(2)+z(12)+z(15)-3.0d0*z(6))
         ddz(12)=ddz(12)+(gam2*3.0d0/bndprm**2)*(z(2)+z(12)+z(15)-3.0d0*z(6))
         ddz(15)=ddz(15)+(gam2*3.0d0/bndprm**2)*(z(2)+z(12)+z(15)-3.0d0*z(6))

         vpuck=vpuck+(z(3)+z(11)+z(16)-3.0d0*z(7))**2
         ddz(7)=ddz(7)-(gam2*9.0d0/bndprm**2)*(z(3)+z(11)+z(16)-3.0d0*z(7))
         ddz(3)=ddz(3)+(gam2*3.0d0/bndprm**2)*(z(3)+z(11)+z(16)-3.0d0*z(7))
         ddz(11)=ddz(11)+(gam2*3.0d0/bndprm**2)*(z(3)+z(11)+z(16)-3.0d0*z(7))
         ddz(16)=ddz(16)+(gam2*3.0d0/bndprm**2)*(z(3)+z(11)+z(16)-3.0d0*z(7))

         vpuck=vpuck+(z(4)+z(12)+z(19)-3.0d0*z(8))**2
         ddz(8)=ddz(8)-(gam2*9.0d0/bndprm**2)*(z(4)+z(12)+z(19)-3.0d0*z(8))
         ddz(4)=ddz(4)+(gam2*3.0d0/bndprm**2)*(z(4)+z(12)+z(19)-3.0d0*z(8))
         ddz(12)=ddz(12)+(gam2*3.0d0/bndprm**2)*(z(4)+z(12)+z(19)-3.0d0*z(8))
         ddz(19)=ddz(19)+(gam2*3.0d0/bndprm**2)*(z(4)+z(12)+z(19)-3.0d0*z(8))

         vpuck=vpuck+(z(3)+z(13)+z(18)-3.0d0*z(9))**2
         ddz(9)=ddz(9)-(gam2*9.0d0/bndprm**2)*(z(3)+z(13)+z(18)-3.0d0*z(9))
         ddz(3)=ddz(3)+(gam2*3.0d0/bndprm**2)*(z(3)+z(13)+z(18)-3.0d0*z(9))
         ddz(13)=ddz(13)+(gam2*3.0d0/bndprm**2)*(z(3)+z(13)+z(18)-3.0d0*z(9))
         ddz(18)=ddz(18)+(gam2*3.0d0/bndprm**2)*(z(3)+z(13)+z(18)-3.0d0*z(9))

         vpuck=vpuck+(z(4)+z(13)+z(17)-3.0d0*z(10))**2
         ddz(10)=ddz(10)-(gam2*9.0d0/bndprm**2)*(z(4)+z(13)+z(17)-3.0d0*z(10))
         ddz(4)=ddz(4)+(gam2*3.0d0/bndprm**2)*(z(4)+z(13)+z(17)-3.0d0*z(10))
         ddz(13)=ddz(13)+(gam2*3.0d0/bndprm**2)*(z(4)+z(13)+z(17)-3.0d0*z(10))
         ddz(17)=ddz(17)+(gam2*3.0d0/bndprm**2)*(z(4)+z(13)+z(17)-3.0d0*z(10))

         vpuck=vpuck+(z(5)+z(7)+z(21)-3.0d0*z(11))**2
         ddz(11)=ddz(11)-(gam2*9.0d0/bndprm**2)*(z(5)+z(7)+z(21)-3.0d0*z(11))
         ddz(5)=ddz(5)+(gam2*3.0d0/bndprm**2)*(z(5)+z(7)+z(21)-3.0d0*z(11))
         ddz(7)=ddz(7)+(gam2*3.0d0/bndprm**2)*(z(5)+z(7)+z(21)-3.0d0*z(11))
         ddz(21)=ddz(21)+(gam2*3.0d0/bndprm**2)*(z(5)+z(7)+z(21)-3.0d0*z(11))

         vpuck=vpuck+(z(6)+z(8)+z(23)-3.0d0*z(12))**2
         ddz(12)=ddz(12)-(gam2*9.0d0/bndprm**2)*(z(6)+z(8)+z(23)-3.0d0*z(12))
         ddz(6)=ddz(6)+(gam2*3.0d0/bndprm**2)*(z(6)+z(8)+z(23)-3.0d0*z(12))
         ddz(8)=ddz(8)+(gam2*3.0d0/bndprm**2)*(z(6)+z(8)+z(23)-3.0d0*z(12))
         ddz(23)=ddz(23)+(gam2*3.0d0/bndprm**2)*(z(6)+z(8)+z(23)-3.0d0*z(12))

         vpuck=vpuck+(z(9)+z(10)+z(24)-3.0d0*z(13))**2
         ddz(13)=ddz(13)-(gam2*9.0d0/bndprm**2)*(z(9)+z(10)+z(24)-3.0d0*z(13))
         ddz(9)=ddz(9)+(gam2*3.0d0/bndprm**2)*(z(9)+z(10)+z(24)-3.0d0*z(13))
         ddz(10)=ddz(10)+(gam2*3.0d0/bndprm**2)*(z(9)+z(10)+z(24)-3.0d0*z(13))
         ddz(24)=ddz(24)+(gam2*3.0d0/bndprm**2)*(z(9)+z(10)+z(24)-3.0d0*z(13))

         vpuck=vpuck+(z(5)+z(20)+z(27)-3.0d0*z(14))**2
         ddz(14)=ddz(14)-(gam2*9.0d0/bndprm**2)*(z(5)+z(20)+z(27)-3.0d0*z(14))
         ddz(5)=ddz(5)+(gam2*3.0d0/bndprm**2)*(z(5)+z(20)+z(27)-3.0d0*z(14))
         ddz(20)=ddz(20)+(gam2*3.0d0/bndprm**2)*(z(5)+z(20)+z(27)-3.0d0*z(14))
         ddz(27)=ddz(27)+(gam2*3.0d0/bndprm**2)*(z(5)+z(20)+z(27)-3.0d0*z(14))

         vpuck=vpuck+(z(6)+z(20)+z(26)-3.0d0*z(15))**2
         ddz(15)=ddz(15)-(gam2*9.0d0/bndprm**2)*(z(6)+z(20)+z(26)-3.0d0*z(15))
         ddz(6)=ddz(6)+(gam2*3.0d0/bndprm**2)*(z(6)+z(20)+z(26)-3.0d0*z(15))
         ddz(20)=ddz(20)+(gam2*3.0d0/bndprm**2)*(z(6)+z(20)+z(26)-3.0d0*z(15))
         ddz(26)=ddz(26)+(gam2*3.0d0/bndprm**2)*(z(6)+z(20)+z(26)-3.0d0*z(15))

         vpuck=vpuck+(z(7)+z(22)+z(28)-3.0d0*z(16))**2
         ddz(16)=ddz(16)-(gam2*9.0d0/bndprm**2)*(z(7)+z(22)+z(28)-3.0d0*z(16))
         ddz(7)=ddz(7)+(gam2*3.0d0/bndprm**2)*(z(7)+z(22)+z(28)-3.0d0*z(16))
         ddz(22)=ddz(22)+(gam2*3.0d0/bndprm**2)*(z(7)+z(22)+z(28)-3.0d0*z(16))
         ddz(28)=ddz(28)+(gam2*3.0d0/bndprm**2)*(z(7)+z(22)+z(28)-3.0d0*z(16))

         vpuck=vpuck+(z(10)+z(25)+z(31)-3.0d0*z(17))**2
         ddz(17)=ddz(17)-(gam2*9.0d0/bndprm**2)*(z(10)+z(25)+z(31)-3.0d0*z(17))
         ddz(10)=ddz(10)+(gam2*3.0d0/bndprm**2)*(z(10)+z(25)+z(31)-3.0d0*z(17))
         ddz(25)=ddz(25)+(gam2*3.0d0/bndprm**2)*(z(10)+z(25)+z(31)-3.0d0*z(17))
         ddz(31)=ddz(31)+(gam2*3.0d0/bndprm**2)*(z(10)+z(25)+z(31)-3.0d0*z(17))

         vpuck=vpuck+(z(9)+z(22)+z(30)-3.0d0*z(18))**2
         ddz(18)=ddz(18)-(gam2*9.0d0/bndprm**2)*(z(9)+z(22)+z(30)-3.0d0*z(18))
         ddz(9)=ddz(9)+(gam2*3.0d0/bndprm**2)*(z(9)+z(22)+z(30)-3.0d0*z(18))
         ddz(22)=ddz(22)+(gam2*3.0d0/bndprm**2)*(z(9)+z(22)+z(30)-3.0d0*z(18))
         ddz(30)=ddz(30)+(gam2*3.0d0/bndprm**2)*(z(9)+z(22)+z(30)-3.0d0*z(18))

         vpuck=vpuck+(z(8)+z(25)+z(29)-3.0d0*z(19))**2
         ddz(19)=ddz(19)-(gam2*9.0d0/bndprm**2)*(z(8)+z(25)+z(29)-3.0d0*z(19))
         ddz(8)=ddz(8)+(gam2*3.0d0/bndprm**2)*(z(8)+z(25)+z(29)-3.0d0*z(19))
         ddz(25)=ddz(25)+(gam2*3.0d0/bndprm**2)*(z(8)+z(25)+z(29)-3.0d0*z(19))
         ddz(29)=ddz(29)+(gam2*3.0d0/bndprm**2)*(z(8)+z(25)+z(29)-3.0d0*z(19))

         vpuck=vpuck+(z(14)+z(15)+z(38)-3.0d0*z(20))**2
         ddz(20)=ddz(20)-(gam2*9.0d0/bndprm**2)*(z(14)+z(15)+z(38)-3.0d0*z(20))
         ddz(14)=ddz(14)+(gam2*3.0d0/bndprm**2)*(z(14)+z(15)+z(38)-3.0d0*z(20))
         ddz(15)=ddz(15)+(gam2*3.0d0/bndprm**2)*(z(14)+z(15)+z(38)-3.0d0*z(20))
         ddz(38)=ddz(38)+(gam2*3.0d0/bndprm**2)*(z(14)+z(15)+z(38)-3.0d0*z(20))

         vpuck=vpuck+(z(11)+z(32)+z(34)-3.0d0*z(21))**2
         ddz(21)=ddz(21)-(gam2*9.0d0/bndprm**2)*(z(11)+z(32)+z(34)-3.0d0*z(21))
         ddz(11)=ddz(11)+(gam2*3.0d0/bndprm**2)*(z(11)+z(32)+z(34)-3.0d0*z(21))
         ddz(32)=ddz(32)+(gam2*3.0d0/bndprm**2)*(z(11)+z(32)+z(34)-3.0d0*z(21))
         ddz(34)=ddz(34)+(gam2*3.0d0/bndprm**2)*(z(11)+z(32)+z(34)-3.0d0*z(21))

         vpuck=vpuck+(z(16)+z(18)+z(39)-3.0d0*z(22))**2
         ddz(22)=ddz(22)-(gam2*9.0d0/bndprm**2)*(z(16)+z(18)+z(39)-3.0d0*z(22))
         ddz(16)=ddz(16)+(gam2*3.0d0/bndprm**2)*(z(16)+z(18)+z(39)-3.0d0*z(22))
         ddz(18)=ddz(18)+(gam2*3.0d0/bndprm**2)*(z(16)+z(18)+z(39)-3.0d0*z(22))
         ddz(39)=ddz(39)+(gam2*3.0d0/bndprm**2)*(z(16)+z(18)+z(39)-3.0d0*z(22))

         vpuck=vpuck+(z(12)+z(33)+z(35)-3.0d0*z(23))**2
         ddz(23)=ddz(23)-(gam2*9.0d0/bndprm**2)*(z(12)+z(33)+z(35)-3.0d0*z(23))
         ddz(12)=ddz(12)+(gam2*3.0d0/bndprm**2)*(z(12)+z(33)+z(35)-3.0d0*z(23))
         ddz(33)=ddz(33)+(gam2*3.0d0/bndprm**2)*(z(12)+z(33)+z(35)-3.0d0*z(23))
         ddz(35)=ddz(35)+(gam2*3.0d0/bndprm**2)*(z(12)+z(33)+z(35)-3.0d0*z(23))

         vpuck=vpuck+(z(13)+z(36)+z(37)-3.0d0*z(24))**2
         ddz(24)=ddz(24)-(gam2*9.0d0/bndprm**2)*(z(13)+z(36)+z(37)-3.0d0*z(24))
         ddz(13)=ddz(13)+(gam2*3.0d0/bndprm**2)*(z(13)+z(36)+z(37)-3.0d0*z(24))
         ddz(36)=ddz(36)+(gam2*3.0d0/bndprm**2)*(z(13)+z(36)+z(37)-3.0d0*z(24))
         ddz(37)=ddz(37)+(gam2*3.0d0/bndprm**2)*(z(13)+z(36)+z(37)-3.0d0*z(24))

         vpuck=vpuck+(z(17)+z(19)+z(40)-3.0d0*z(25))**2
         ddz(25)=ddz(25)-(gam2*9.0d0/bndprm**2)*(z(17)+z(19)+z(40)-3.0d0*z(25))
         ddz(17)=ddz(17)+(gam2*3.0d0/bndprm**2)*(z(17)+z(19)+z(40)-3.0d0*z(25))
         ddz(19)=ddz(19)+(gam2*3.0d0/bndprm**2)*(z(17)+z(19)+z(40)-3.0d0*z(25))
         ddz(40)=ddz(40)+(gam2*3.0d0/bndprm**2)*(z(17)+z(19)+z(40)-3.0d0*z(25))

         vpuck=vpuck+(z(15)+z(35)+z(43)-3.0d0*z(26))**2
         ddz(26)=ddz(26)-(gam2*9.0d0/bndprm**2)*(z(15)+z(35)+z(43)-3.0d0*z(26))
         ddz(15)=ddz(15)+(gam2*3.0d0/bndprm**2)*(z(15)+z(35)+z(43)-3.0d0*z(26))
         ddz(35)=ddz(35)+(gam2*3.0d0/bndprm**2)*(z(15)+z(35)+z(43)-3.0d0*z(26))
         ddz(43)=ddz(43)+(gam2*3.0d0/bndprm**2)*(z(15)+z(35)+z(43)-3.0d0*z(26))

         vpuck=vpuck+(z(14)+z(34)+z(42)-3.0d0*z(27))**2
         ddz(27)=ddz(27)-(gam2*9.0d0/bndprm**2)*(z(14)+z(34)+z(42)-3.0d0*z(27))
         ddz(14)=ddz(14)+(gam2*3.0d0/bndprm**2)*(z(14)+z(34)+z(42)-3.0d0*z(27))
         ddz(34)=ddz(34)+(gam2*3.0d0/bndprm**2)*(z(14)+z(34)+z(42)-3.0d0*z(27))
         ddz(42)=ddz(42)+(gam2*3.0d0/bndprm**2)*(z(14)+z(34)+z(42)-3.0d0*z(27))

         vpuck=vpuck+(z(16)+z(32)+z(45)-3.0d0*z(28))**2
         ddz(28)=ddz(28)-(gam2*9.0d0/bndprm**2)*(z(16)+z(32)+z(45)-3.0d0*z(28))
         ddz(16)=ddz(16)+(gam2*3.0d0/bndprm**2)*(z(16)+z(32)+z(45)-3.0d0*z(28))
         ddz(32)=ddz(32)+(gam2*3.0d0/bndprm**2)*(z(16)+z(32)+z(45)-3.0d0*z(28))
         ddz(45)=ddz(45)+(gam2*3.0d0/bndprm**2)*(z(16)+z(32)+z(45)-3.0d0*z(28))

         vpuck=vpuck+(z(19)+z(33)+z(41)-3.0d0*z(29))**2
         ddz(29)=ddz(29)-(gam2*9.0d0/bndprm**2)*(z(19)+z(33)+z(41)-3.0d0*z(29))
         ddz(19)=ddz(19)+(gam2*3.0d0/bndprm**2)*(z(19)+z(33)+z(41)-3.0d0*z(29))
         ddz(33)=ddz(33)+(gam2*3.0d0/bndprm**2)*(z(19)+z(33)+z(41)-3.0d0*z(29))
         ddz(41)=ddz(41)+(gam2*3.0d0/bndprm**2)*(z(19)+z(33)+z(41)-3.0d0*z(29))

         vpuck=vpuck+(z(18)+z(36)+z(46)-3.0d0*z(30))**2
         ddz(30)=ddz(30)-(gam2*9.0d0/bndprm**2)*(z(18)+z(36)+z(46)-3.0d0*z(30))
         ddz(18)=ddz(18)+(gam2*3.0d0/bndprm**2)*(z(18)+z(36)+z(46)-3.0d0*z(30))
         ddz(36)=ddz(36)+(gam2*3.0d0/bndprm**2)*(z(18)+z(36)+z(46)-3.0d0*z(30))
         ddz(46)=ddz(46)+(gam2*3.0d0/bndprm**2)*(z(18)+z(36)+z(46)-3.0d0*z(30))

         vpuck=vpuck+(z(17)+z(37)+z(44)-3.0d0*z(31))**2
         ddz(31)=ddz(31)-(gam2*9.0d0/bndprm**2)*(z(17)+z(37)+z(44)-3.0d0*z(31))
         ddz(17)=ddz(17)+(gam2*3.0d0/bndprm**2)*(z(17)+z(37)+z(44)-3.0d0*z(31))
         ddz(37)=ddz(37)+(gam2*3.0d0/bndprm**2)*(z(17)+z(37)+z(44)-3.0d0*z(31))
         ddz(44)=ddz(44)+(gam2*3.0d0/bndprm**2)*(z(17)+z(37)+z(44)-3.0d0*z(31))

         vpuck=vpuck+(z(21)+z(28)+z(51)-3.0d0*z(32))**2
         ddz(32)=ddz(32)-(gam2*9.0d0/bndprm**2)*(z(21)+z(28)+z(51)-3.0d0*z(32))
         ddz(21)=ddz(21)+(gam2*3.0d0/bndprm**2)*(z(21)+z(28)+z(51)-3.0d0*z(32))
         ddz(28)=ddz(28)+(gam2*3.0d0/bndprm**2)*(z(21)+z(28)+z(51)-3.0d0*z(32))
         ddz(51)=ddz(51)+(gam2*3.0d0/bndprm**2)*(z(21)+z(28)+z(51)-3.0d0*z(32))

         vpuck=vpuck+(z(23)+z(29)+z(49)-3.0d0*z(33))**2
         ddz(33)=ddz(33)-(gam2*9.0d0/bndprm**2)*(z(23)+z(29)+z(49)-3.0d0*z(33))
         ddz(23)=ddz(23)+(gam2*3.0d0/bndprm**2)*(z(23)+z(29)+z(49)-3.0d0*z(33))
         ddz(29)=ddz(29)+(gam2*3.0d0/bndprm**2)*(z(23)+z(29)+z(49)-3.0d0*z(33))
         ddz(49)=ddz(49)+(gam2*3.0d0/bndprm**2)*(z(23)+z(29)+z(49)-3.0d0*z(33))

         vpuck=vpuck+(z(21)+z(27)+z(52)-3.0d0*z(34))**2
         ddz(34)=ddz(34)-(gam2*9.0d0/bndprm**2)*(z(21)+z(27)+z(52)-3.0d0*z(34))
         ddz(21)=ddz(21)+(gam2*3.0d0/bndprm**2)*(z(21)+z(27)+z(52)-3.0d0*z(34))
         ddz(27)=ddz(27)+(gam2*3.0d0/bndprm**2)*(z(21)+z(27)+z(52)-3.0d0*z(34))
         ddz(52)=ddz(52)+(gam2*3.0d0/bndprm**2)*(z(21)+z(27)+z(52)-3.0d0*z(34))

         vpuck=vpuck+(z(23)+z(26)+z(53)-3.0d0*z(35))**2
         ddz(35)=ddz(35)-(gam2*9.0d0/bndprm**2)*(z(23)+z(26)+z(53)-3.0d0*z(35))
         ddz(23)=ddz(23)+(gam2*3.0d0/bndprm**2)*(z(23)+z(26)+z(53)-3.0d0*z(35))
         ddz(26)=ddz(26)+(gam2*3.0d0/bndprm**2)*(z(23)+z(26)+z(53)-3.0d0*z(35))
         ddz(53)=ddz(53)+(gam2*3.0d0/bndprm**2)*(z(23)+z(26)+z(53)-3.0d0*z(35))

         vpuck=vpuck+(z(24)+z(30)+z(57)-3.0d0*z(36))**2
         ddz(36)=ddz(36)-(gam2*9.0d0/bndprm**2)*(z(24)+z(30)+z(57)-3.0d0*z(36))
         ddz(24)=ddz(24)+(gam2*3.0d0/bndprm**2)*(z(24)+z(30)+z(57)-3.0d0*z(36))
         ddz(30)=ddz(30)+(gam2*3.0d0/bndprm**2)*(z(24)+z(30)+z(57)-3.0d0*z(36))
         ddz(57)=ddz(57)+(gam2*3.0d0/bndprm**2)*(z(24)+z(30)+z(57)-3.0d0*z(36))

         vpuck=vpuck+(z(24)+z(31)+z(58)-3.0d0*z(37))**2
         ddz(37)=ddz(37)-(gam2*9.0d0/bndprm**2)*(z(24)+z(31)+z(58)-3.0d0*z(37))
         ddz(24)=ddz(24)+(gam2*3.0d0/bndprm**2)*(z(24)+z(31)+z(58)-3.0d0*z(37))
         ddz(31)=ddz(31)+(gam2*3.0d0/bndprm**2)*(z(24)+z(31)+z(58)-3.0d0*z(37))
         ddz(58)=ddz(58)+(gam2*3.0d0/bndprm**2)*(z(24)+z(31)+z(58)-3.0d0*z(37))

         vpuck=vpuck+(z(20)+z(47)+z(48)-3.0d0*z(38))**2
         ddz(38)=ddz(38)-(gam2*9.0d0/bndprm**2)*(z(20)+z(47)+z(48)-3.0d0*z(38))
         ddz(20)=ddz(20)+(gam2*3.0d0/bndprm**2)*(z(20)+z(47)+z(48)-3.0d0*z(38))
         ddz(47)=ddz(47)+(gam2*3.0d0/bndprm**2)*(z(20)+z(47)+z(48)-3.0d0*z(38))
         ddz(48)=ddz(48)+(gam2*3.0d0/bndprm**2)*(z(20)+z(47)+z(48)-3.0d0*z(38))

         vpuck=vpuck+(z(22)+z(54)+z(56)-3.0d0*z(39))**2
         ddz(39)=ddz(39)-(gam2*9.0d0/bndprm**2)*(z(22)+z(54)+z(56)-3.0d0*z(39))
         ddz(22)=ddz(22)+(gam2*3.0d0/bndprm**2)*(z(22)+z(54)+z(56)-3.0d0*z(39))
         ddz(54)=ddz(54)+(gam2*3.0d0/bndprm**2)*(z(22)+z(54)+z(56)-3.0d0*z(39))
         ddz(56)=ddz(56)+(gam2*3.0d0/bndprm**2)*(z(22)+z(54)+z(56)-3.0d0*z(39))

         vpuck=vpuck+(z(25)+z(50)+z(55)-3.0d0*z(40))**2
         ddz(40)=ddz(40)-(gam2*9.0d0/bndprm**2)*(z(25)+z(50)+z(55)-3.0d0*z(40))
         ddz(25)=ddz(25)+(gam2*3.0d0/bndprm**2)*(z(25)+z(50)+z(55)-3.0d0*z(40))
         ddz(50)=ddz(50)+(gam2*3.0d0/bndprm**2)*(z(25)+z(50)+z(55)-3.0d0*z(40))
         ddz(55)=ddz(55)+(gam2*3.0d0/bndprm**2)*(z(25)+z(50)+z(55)-3.0d0*z(40))

         vpuck=vpuck+(z(29)+z(50)+z(67)-3.0d0*z(41))**2
         ddz(41)=ddz(41)-(gam2*9.0d0/bndprm**2)*(z(29)+z(50)+z(67)-3.0d0*z(41))
         ddz(29)=ddz(29)+(gam2*3.0d0/bndprm**2)*(z(29)+z(50)+z(67)-3.0d0*z(41))
         ddz(50)=ddz(50)+(gam2*3.0d0/bndprm**2)*(z(29)+z(50)+z(67)-3.0d0*z(41))
         ddz(67)=ddz(67)+(gam2*3.0d0/bndprm**2)*(z(29)+z(50)+z(67)-3.0d0*z(41))

         vpuck=vpuck+(z(27)+z(47)+z(62)-3.0d0*z(42))**2
         ddz(42)=ddz(42)-(gam2*9.0d0/bndprm**2)*(z(27)+z(47)+z(62)-3.0d0*z(42))
         ddz(27)=ddz(27)+(gam2*3.0d0/bndprm**2)*(z(27)+z(47)+z(62)-3.0d0*z(42))
         ddz(47)=ddz(47)+(gam2*3.0d0/bndprm**2)*(z(27)+z(47)+z(62)-3.0d0*z(42))
         ddz(62)=ddz(62)+(gam2*3.0d0/bndprm**2)*(z(27)+z(47)+z(62)-3.0d0*z(42))

         vpuck=vpuck+(z(26)+z(48)+z(63)-3.0d0*z(43))**2
         ddz(43)=ddz(43)-(gam2*9.0d0/bndprm**2)*(z(26)+z(48)+z(63)-3.0d0*z(43))
         ddz(26)=ddz(26)+(gam2*3.0d0/bndprm**2)*(z(26)+z(48)+z(63)-3.0d0*z(43))
         ddz(48)=ddz(48)+(gam2*3.0d0/bndprm**2)*(z(26)+z(48)+z(63)-3.0d0*z(43))
         ddz(63)=ddz(63)+(gam2*3.0d0/bndprm**2)*(z(26)+z(48)+z(63)-3.0d0*z(43))

         vpuck=vpuck+(z(31)+z(55)+z(66)-3.0d0*z(44))**2
         ddz(44)=ddz(44)-(gam2*9.0d0/bndprm**2)*(z(31)+z(55)+z(66)-3.0d0*z(44))
         ddz(31)=ddz(31)+(gam2*3.0d0/bndprm**2)*(z(31)+z(55)+z(66)-3.0d0*z(44))
         ddz(55)=ddz(55)+(gam2*3.0d0/bndprm**2)*(z(31)+z(55)+z(66)-3.0d0*z(44))
         ddz(66)=ddz(66)+(gam2*3.0d0/bndprm**2)*(z(31)+z(55)+z(66)-3.0d0*z(44))

         vpuck=vpuck+(z(28)+z(56)+z(64)-3.0d0*z(45))**2
         ddz(45)=ddz(45)-(gam2*9.0d0/bndprm**2)*(z(28)+z(56)+z(64)-3.0d0*z(45))
         ddz(28)=ddz(28)+(gam2*3.0d0/bndprm**2)*(z(28)+z(56)+z(64)-3.0d0*z(45))
         ddz(56)=ddz(56)+(gam2*3.0d0/bndprm**2)*(z(28)+z(56)+z(64)-3.0d0*z(45))
         ddz(64)=ddz(64)+(gam2*3.0d0/bndprm**2)*(z(28)+z(56)+z(64)-3.0d0*z(45))

         vpuck=vpuck+(z(30)+z(54)+z(65)-3.0d0*z(46))**2
         ddz(46)=ddz(46)-(gam2*9.0d0/bndprm**2)*(z(30)+z(54)+z(65)-3.0d0*z(46))
         ddz(30)=ddz(30)+(gam2*3.0d0/bndprm**2)*(z(30)+z(54)+z(65)-3.0d0*z(46))
         ddz(54)=ddz(54)+(gam2*3.0d0/bndprm**2)*(z(30)+z(54)+z(65)-3.0d0*z(46))
         ddz(65)=ddz(65)+(gam2*3.0d0/bndprm**2)*(z(30)+z(54)+z(65)-3.0d0*z(46))

         vpuck=vpuck+(z(38)+z(42)+z(75)-3.0d0*z(47))**2
         ddz(47)=ddz(47)-(gam2*9.0d0/bndprm**2)*(z(38)+z(42)+z(75)-3.0d0*z(47))
         ddz(38)=ddz(38)+(gam2*3.0d0/bndprm**2)*(z(38)+z(42)+z(75)-3.0d0*z(47))
         ddz(42)=ddz(42)+(gam2*3.0d0/bndprm**2)*(z(38)+z(42)+z(75)-3.0d0*z(47))
         ddz(75)=ddz(75)+(gam2*3.0d0/bndprm**2)*(z(38)+z(42)+z(75)-3.0d0*z(47))

         vpuck=vpuck+(z(38)+z(43)+z(76)-3.0d0*z(48))**2
         ddz(48)=ddz(48)-(gam2*9.0d0/bndprm**2)*(z(38)+z(43)+z(76)-3.0d0*z(48))
         ddz(38)=ddz(38)+(gam2*3.0d0/bndprm**2)*(z(38)+z(43)+z(76)-3.0d0*z(48))
         ddz(43)=ddz(43)+(gam2*3.0d0/bndprm**2)*(z(38)+z(43)+z(76)-3.0d0*z(48))
         ddz(76)=ddz(76)+(gam2*3.0d0/bndprm**2)*(z(38)+z(43)+z(76)-3.0d0*z(48))

         vpuck=vpuck+(z(33)+z(60)+z(73)-3.0d0*z(49))**2
         ddz(49)=ddz(49)-(gam2*9.0d0/bndprm**2)*(z(33)+z(60)+z(73)-3.0d0*z(49))
         ddz(33)=ddz(33)+(gam2*3.0d0/bndprm**2)*(z(33)+z(60)+z(73)-3.0d0*z(49))
         ddz(60)=ddz(60)+(gam2*3.0d0/bndprm**2)*(z(33)+z(60)+z(73)-3.0d0*z(49))
         ddz(73)=ddz(73)+(gam2*3.0d0/bndprm**2)*(z(33)+z(60)+z(73)-3.0d0*z(49))

         vpuck=vpuck+(z(40)+z(41)+z(74)-3.0d0*z(50))**2
         ddz(50)=ddz(50)-(gam2*9.0d0/bndprm**2)*(z(40)+z(41)+z(74)-3.0d0*z(50))
         ddz(40)=ddz(40)+(gam2*3.0d0/bndprm**2)*(z(40)+z(41)+z(74)-3.0d0*z(50))
         ddz(41)=ddz(41)+(gam2*3.0d0/bndprm**2)*(z(40)+z(41)+z(74)-3.0d0*z(50))
         ddz(74)=ddz(74)+(gam2*3.0d0/bndprm**2)*(z(40)+z(41)+z(74)-3.0d0*z(50))

         vpuck=vpuck+(z(32)+z(59)+z(70)-3.0d0*z(51))**2
         ddz(51)=ddz(51)-(gam2*9.0d0/bndprm**2)*(z(32)+z(59)+z(70)-3.0d0*z(51))
         ddz(32)=ddz(32)+(gam2*3.0d0/bndprm**2)*(z(32)+z(59)+z(70)-3.0d0*z(51))
         ddz(59)=ddz(59)+(gam2*3.0d0/bndprm**2)*(z(32)+z(59)+z(70)-3.0d0*z(51))
         ddz(70)=ddz(70)+(gam2*3.0d0/bndprm**2)*(z(32)+z(59)+z(70)-3.0d0*z(51))

         vpuck=vpuck+(z(34)+z(59)+z(68)-3.0d0*z(52))**2
         ddz(52)=ddz(52)-(gam2*9.0d0/bndprm**2)*(z(34)+z(59)+z(68)-3.0d0*z(52))
         ddz(34)=ddz(34)+(gam2*3.0d0/bndprm**2)*(z(34)+z(59)+z(68)-3.0d0*z(52))
         ddz(59)=ddz(59)+(gam2*3.0d0/bndprm**2)*(z(34)+z(59)+z(68)-3.0d0*z(52))
         ddz(68)=ddz(68)+(gam2*3.0d0/bndprm**2)*(z(34)+z(59)+z(68)-3.0d0*z(52))

         vpuck=vpuck+(z(35)+z(60)+z(69)-3.0d0*z(53))**2
         ddz(53)=ddz(53)-(gam2*9.0d0/bndprm**2)*(z(35)+z(60)+z(69)-3.0d0*z(53))
         ddz(35)=ddz(35)+(gam2*3.0d0/bndprm**2)*(z(35)+z(60)+z(69)-3.0d0*z(53))
         ddz(60)=ddz(60)+(gam2*3.0d0/bndprm**2)*(z(35)+z(60)+z(69)-3.0d0*z(53))
         ddz(69)=ddz(69)+(gam2*3.0d0/bndprm**2)*(z(35)+z(60)+z(69)-3.0d0*z(53))

         vpuck=vpuck+(z(39)+z(46)+z(78)-3.0d0*z(54))**2
         ddz(54)=ddz(54)-(gam2*9.0d0/bndprm**2)*(z(39)+z(46)+z(78)-3.0d0*z(54))
         ddz(39)=ddz(39)+(gam2*3.0d0/bndprm**2)*(z(39)+z(46)+z(78)-3.0d0*z(54))
         ddz(46)=ddz(46)+(gam2*3.0d0/bndprm**2)*(z(39)+z(46)+z(78)-3.0d0*z(54))
         ddz(78)=ddz(78)+(gam2*3.0d0/bndprm**2)*(z(39)+z(46)+z(78)-3.0d0*z(54))

         vpuck=vpuck+(z(40)+z(44)+z(79)-3.0d0*z(55))**2
         ddz(55)=ddz(55)-(gam2*9.0d0/bndprm**2)*(z(40)+z(44)+z(79)-3.0d0*z(55))
         ddz(40)=ddz(40)+(gam2*3.0d0/bndprm**2)*(z(40)+z(44)+z(79)-3.0d0*z(55))
         ddz(44)=ddz(44)+(gam2*3.0d0/bndprm**2)*(z(40)+z(44)+z(79)-3.0d0*z(55))
         ddz(79)=ddz(79)+(gam2*3.0d0/bndprm**2)*(z(40)+z(44)+z(79)-3.0d0*z(55))

         vpuck=vpuck+(z(39)+z(45)+z(77)-3.0d0*z(56))**2
         ddz(56)=ddz(56)-(gam2*9.0d0/bndprm**2)*(z(39)+z(45)+z(77)-3.0d0*z(56))
         ddz(39)=ddz(39)+(gam2*3.0d0/bndprm**2)*(z(39)+z(45)+z(77)-3.0d0*z(56))
         ddz(45)=ddz(45)+(gam2*3.0d0/bndprm**2)*(z(39)+z(45)+z(77)-3.0d0*z(56))
         ddz(77)=ddz(77)+(gam2*3.0d0/bndprm**2)*(z(39)+z(45)+z(77)-3.0d0*z(56))

         vpuck=vpuck+(z(36)+z(61)+z(71)-3.0d0*z(57))**2
         ddz(57)=ddz(57)-(gam2*9.0d0/bndprm**2)*(z(36)+z(61)+z(71)-3.0d0*z(57))
         ddz(36)=ddz(36)+(gam2*3.0d0/bndprm**2)*(z(36)+z(61)+z(71)-3.0d0*z(57))
         ddz(61)=ddz(61)+(gam2*3.0d0/bndprm**2)*(z(36)+z(61)+z(71)-3.0d0*z(57))
         ddz(71)=ddz(71)+(gam2*3.0d0/bndprm**2)*(z(36)+z(61)+z(71)-3.0d0*z(57))

         vpuck=vpuck+(z(37)+z(61)+z(72)-3.0d0*z(58))**2
         ddz(58)=ddz(58)-(gam2*9.0d0/bndprm**2)*(z(37)+z(61)+z(72)-3.0d0*z(58))
         ddz(37)=ddz(37)+(gam2*3.0d0/bndprm**2)*(z(37)+z(61)+z(72)-3.0d0*z(58))
         ddz(61)=ddz(61)+(gam2*3.0d0/bndprm**2)*(z(37)+z(61)+z(72)-3.0d0*z(58))
         ddz(72)=ddz(72)+(gam2*3.0d0/bndprm**2)*(z(37)+z(61)+z(72)-3.0d0*z(58))

         vpuck=vpuck+(z(51)+z(52)+z(81)-3.0d0*z(59))**2
         ddz(59)=ddz(59)-(gam2*9.0d0/bndprm**2)*(z(51)+z(52)+z(81)-3.0d0*z(59))
         ddz(51)=ddz(51)+(gam2*3.0d0/bndprm**2)*(z(51)+z(52)+z(81)-3.0d0*z(59))
         ddz(52)=ddz(52)+(gam2*3.0d0/bndprm**2)*(z(51)+z(52)+z(81)-3.0d0*z(59))
         ddz(81)=ddz(81)+(gam2*3.0d0/bndprm**2)*(z(51)+z(52)+z(81)-3.0d0*z(59))

         vpuck=vpuck+(z(49)+z(53)+z(80)-3.0d0*z(60))**2
         ddz(60)=ddz(60)-(gam2*9.0d0/bndprm**2)*(z(49)+z(53)+z(80)-3.0d0*z(60))
         ddz(49)=ddz(49)+(gam2*3.0d0/bndprm**2)*(z(49)+z(53)+z(80)-3.0d0*z(60))
         ddz(53)=ddz(53)+(gam2*3.0d0/bndprm**2)*(z(49)+z(53)+z(80)-3.0d0*z(60))
         ddz(80)=ddz(80)+(gam2*3.0d0/bndprm**2)*(z(49)+z(53)+z(80)-3.0d0*z(60))

         vpuck=vpuck+(z(57)+z(58)+z(85)-3.0d0*z(61))**2
         ddz(61)=ddz(61)-(gam2*9.0d0/bndprm**2)*(z(57)+z(58)+z(85)-3.0d0*z(61))
         ddz(57)=ddz(57)+(gam2*3.0d0/bndprm**2)*(z(57)+z(58)+z(85)-3.0d0*z(61))
         ddz(58)=ddz(58)+(gam2*3.0d0/bndprm**2)*(z(57)+z(58)+z(85)-3.0d0*z(61))
         ddz(85)=ddz(85)+(gam2*3.0d0/bndprm**2)*(z(57)+z(58)+z(85)-3.0d0*z(61))

         vpuck=vpuck+(z(42)+z(68)+z(86)-3.0d0*z(62))**2
         ddz(62)=ddz(62)-(gam2*9.0d0/bndprm**2)*(z(42)+z(68)+z(86)-3.0d0*z(62))
         ddz(42)=ddz(42)+(gam2*3.0d0/bndprm**2)*(z(42)+z(68)+z(86)-3.0d0*z(62))
         ddz(68)=ddz(68)+(gam2*3.0d0/bndprm**2)*(z(42)+z(68)+z(86)-3.0d0*z(62))
         ddz(86)=ddz(86)+(gam2*3.0d0/bndprm**2)*(z(42)+z(68)+z(86)-3.0d0*z(62))

         vpuck=vpuck+(z(43)+z(69)+z(87)-3.0d0*z(63))**2
         ddz(63)=ddz(63)-(gam2*9.0d0/bndprm**2)*(z(43)+z(69)+z(87)-3.0d0*z(63))
         ddz(43)=ddz(43)+(gam2*3.0d0/bndprm**2)*(z(43)+z(69)+z(87)-3.0d0*z(63))
         ddz(69)=ddz(69)+(gam2*3.0d0/bndprm**2)*(z(43)+z(69)+z(87)-3.0d0*z(63))
         ddz(87)=ddz(87)+(gam2*3.0d0/bndprm**2)*(z(43)+z(69)+z(87)-3.0d0*z(63))

         vpuck=vpuck+(z(45)+z(70)+z(90)-3.0d0*z(64))**2
         ddz(64)=ddz(64)-(gam2*9.0d0/bndprm**2)*(z(45)+z(70)+z(90)-3.0d0*z(64))
         ddz(45)=ddz(45)+(gam2*3.0d0/bndprm**2)*(z(45)+z(70)+z(90)-3.0d0*z(64))
         ddz(70)=ddz(70)+(gam2*3.0d0/bndprm**2)*(z(45)+z(70)+z(90)-3.0d0*z(64))
         ddz(90)=ddz(90)+(gam2*3.0d0/bndprm**2)*(z(45)+z(70)+z(90)-3.0d0*z(64))

         vpuck=vpuck+(z(46)+z(71)+z(88)-3.0d0*z(65))**2
         ddz(65)=ddz(65)-(gam2*9.0d0/bndprm**2)*(z(46)+z(71)+z(88)-3.0d0*z(65))
         ddz(46)=ddz(46)+(gam2*3.0d0/bndprm**2)*(z(46)+z(71)+z(88)-3.0d0*z(65))
         ddz(71)=ddz(71)+(gam2*3.0d0/bndprm**2)*(z(46)+z(71)+z(88)-3.0d0*z(65))
         ddz(88)=ddz(88)+(gam2*3.0d0/bndprm**2)*(z(46)+z(71)+z(88)-3.0d0*z(65))

         vpuck=vpuck+(z(44)+z(72)+z(89)-3.0d0*z(66))**2
         ddz(66)=ddz(66)-(gam2*9.0d0/bndprm**2)*(z(44)+z(72)+z(89)-3.0d0*z(66))
         ddz(44)=ddz(44)+(gam2*3.0d0/bndprm**2)*(z(44)+z(72)+z(89)-3.0d0*z(66))
         ddz(72)=ddz(72)+(gam2*3.0d0/bndprm**2)*(z(44)+z(72)+z(89)-3.0d0*z(66))
         ddz(89)=ddz(89)+(gam2*3.0d0/bndprm**2)*(z(44)+z(72)+z(89)-3.0d0*z(66))

         vpuck=vpuck+(z(41)+z(73)+z(91)-3.0d0*z(67))**2
         ddz(67)=ddz(67)-(gam2*9.0d0/bndprm**2)*(z(41)+z(73)+z(91)-3.0d0*z(67))
         ddz(41)=ddz(41)+(gam2*3.0d0/bndprm**2)*(z(41)+z(73)+z(91)-3.0d0*z(67))
         ddz(73)=ddz(73)+(gam2*3.0d0/bndprm**2)*(z(41)+z(73)+z(91)-3.0d0*z(67))
         ddz(91)=ddz(91)+(gam2*3.0d0/bndprm**2)*(z(41)+z(73)+z(91)-3.0d0*z(67))

         vpuck=vpuck+(z(52)+z(62)+z(94)-3.0d0*z(68))**2
         ddz(68)=ddz(68)-(gam2*9.0d0/bndprm**2)*(z(52)+z(62)+z(94)-3.0d0*z(68))
         ddz(52)=ddz(52)+(gam2*3.0d0/bndprm**2)*(z(52)+z(62)+z(94)-3.0d0*z(68))
         ddz(62)=ddz(62)+(gam2*3.0d0/bndprm**2)*(z(52)+z(62)+z(94)-3.0d0*z(68))
         ddz(94)=ddz(94)+(gam2*3.0d0/bndprm**2)*(z(52)+z(62)+z(94)-3.0d0*z(68))

         vpuck=vpuck+(z(53)+z(63)+z(95)-3.0d0*z(69))**2
         ddz(69)=ddz(69)-(gam2*9.0d0/bndprm**2)*(z(53)+z(63)+z(95)-3.0d0*z(69))
         ddz(53)=ddz(53)+(gam2*3.0d0/bndprm**2)*(z(53)+z(63)+z(95)-3.0d0*z(69))
         ddz(63)=ddz(63)+(gam2*3.0d0/bndprm**2)*(z(53)+z(63)+z(95)-3.0d0*z(69))
         ddz(95)=ddz(95)+(gam2*3.0d0/bndprm**2)*(z(53)+z(63)+z(95)-3.0d0*z(69))

         vpuck=vpuck+(z(51)+z(64)+z(92)-3.0d0*z(70))**2
         ddz(70)=ddz(70)-(gam2*9.0d0/bndprm**2)*(z(51)+z(64)+z(92)-3.0d0*z(70))
         ddz(51)=ddz(51)+(gam2*3.0d0/bndprm**2)*(z(51)+z(64)+z(92)-3.0d0*z(70))
         ddz(64)=ddz(64)+(gam2*3.0d0/bndprm**2)*(z(51)+z(64)+z(92)-3.0d0*z(70))
         ddz(92)=ddz(92)+(gam2*3.0d0/bndprm**2)*(z(51)+z(64)+z(92)-3.0d0*z(70))

         vpuck=vpuck+(z(57)+z(65)+z(102)-3.0d0*z(71))**2
         ddz(71)=ddz(71)-(gam2*9.0d0/bndprm**2)*(z(57)+z(65)+z(102)-3.0d0*z(71))
         ddz(57)=ddz(57)+(gam2*3.0d0/bndprm**2)*(z(57)+z(65)+z(102)-3.0d0*z(71))
         ddz(65)=ddz(65)+(gam2*3.0d0/bndprm**2)*(z(57)+z(65)+z(102)-3.0d0*z(71))
         ddz(102)=ddz(102)+(gam2*3.0d0/bndprm**2)*(z(57)+z(65)+z(102)-3.0d0*z(71))

         vpuck=vpuck+(z(58)+z(66)+z(103)-3.0d0*z(72))**2
         ddz(72)=ddz(72)-(gam2*9.0d0/bndprm**2)*(z(58)+z(66)+z(103)-3.0d0*z(72))
         ddz(58)=ddz(58)+(gam2*3.0d0/bndprm**2)*(z(58)+z(66)+z(103)-3.0d0*z(72))
         ddz(66)=ddz(66)+(gam2*3.0d0/bndprm**2)*(z(58)+z(66)+z(103)-3.0d0*z(72))
         ddz(103)=ddz(103)+(gam2*3.0d0/bndprm**2)*(z(58)+z(66)+z(103)-3.0d0*z(72))

         vpuck=vpuck+(z(49)+z(67)+z(93)-3.0d0*z(73))**2
         ddz(73)=ddz(73)-(gam2*9.0d0/bndprm**2)*(z(49)+z(67)+z(93)-3.0d0*z(73))
         ddz(49)=ddz(49)+(gam2*3.0d0/bndprm**2)*(z(49)+z(67)+z(93)-3.0d0*z(73))
         ddz(67)=ddz(67)+(gam2*3.0d0/bndprm**2)*(z(49)+z(67)+z(93)-3.0d0*z(73))
         ddz(93)=ddz(93)+(gam2*3.0d0/bndprm**2)*(z(49)+z(67)+z(93)-3.0d0*z(73))

         vpuck=vpuck+(z(50)+z(84)+z(101)-3.0d0*z(74))**2
         ddz(74)=ddz(74)-(gam2*9.0d0/bndprm**2)*(z(50)+z(84)+z(101)-3.0d0*z(74))
         ddz(50)=ddz(50)+(gam2*3.0d0/bndprm**2)*(z(50)+z(84)+z(101)-3.0d0*z(74))
         ddz(84)=ddz(84)+(gam2*3.0d0/bndprm**2)*(z(50)+z(84)+z(101)-3.0d0*z(74))
         ddz(101)=ddz(101)+(gam2*3.0d0/bndprm**2)*(z(50)+z(84)+z(101)-3.0d0*z(74))

         vpuck=vpuck+(z(47)+z(82)+z(96)-3.0d0*z(75))**2
         ddz(75)=ddz(75)-(gam2*9.0d0/bndprm**2)*(z(47)+z(82)+z(96)-3.0d0*z(75))
         ddz(47)=ddz(47)+(gam2*3.0d0/bndprm**2)*(z(47)+z(82)+z(96)-3.0d0*z(75))
         ddz(82)=ddz(82)+(gam2*3.0d0/bndprm**2)*(z(47)+z(82)+z(96)-3.0d0*z(75))
         ddz(96)=ddz(96)+(gam2*3.0d0/bndprm**2)*(z(47)+z(82)+z(96)-3.0d0*z(75))

         vpuck=vpuck+(z(48)+z(82)+z(97)-3.0d0*z(76))**2
         ddz(76)=ddz(76)-(gam2*9.0d0/bndprm**2)*(z(48)+z(82)+z(97)-3.0d0*z(76))
         ddz(48)=ddz(48)+(gam2*3.0d0/bndprm**2)*(z(48)+z(82)+z(97)-3.0d0*z(76))
         ddz(82)=ddz(82)+(gam2*3.0d0/bndprm**2)*(z(48)+z(82)+z(97)-3.0d0*z(76))
         ddz(97)=ddz(97)+(gam2*3.0d0/bndprm**2)*(z(48)+z(82)+z(97)-3.0d0*z(76))

         vpuck=vpuck+(z(56)+z(83)+z(100)-3.0d0*z(77))**2
         ddz(77)=ddz(77)-(gam2*9.0d0/bndprm**2)*(z(56)+z(83)+z(100)-3.0d0*z(77))
         ddz(56)=ddz(56)+(gam2*3.0d0/bndprm**2)*(z(56)+z(83)+z(100)-3.0d0*z(77))
         ddz(83)=ddz(83)+(gam2*3.0d0/bndprm**2)*(z(56)+z(83)+z(100)-3.0d0*z(77))
         ddz(100)=ddz(100)+(gam2*3.0d0/bndprm**2)*(z(56)+z(83)+z(100)-3.0d0*z(77))

         vpuck=vpuck+(z(54)+z(83)+z(98)-3.0d0*z(78))**2
         ddz(78)=ddz(78)-(gam2*9.0d0/bndprm**2)*(z(54)+z(83)+z(98)-3.0d0*z(78))
         ddz(54)=ddz(54)+(gam2*3.0d0/bndprm**2)*(z(54)+z(83)+z(98)-3.0d0*z(78))
         ddz(83)=ddz(83)+(gam2*3.0d0/bndprm**2)*(z(54)+z(83)+z(98)-3.0d0*z(78))
         ddz(98)=ddz(98)+(gam2*3.0d0/bndprm**2)*(z(54)+z(83)+z(98)-3.0d0*z(78))

         vpuck=vpuck+(z(55)+z(84)+z(99)-3.0d0*z(79))**2
         ddz(79)=ddz(79)-(gam2*9.0d0/bndprm**2)*(z(55)+z(84)+z(99)-3.0d0*z(79))
         ddz(55)=ddz(55)+(gam2*3.0d0/bndprm**2)*(z(55)+z(84)+z(99)-3.0d0*z(79))
         ddz(84)=ddz(84)+(gam2*3.0d0/bndprm**2)*(z(55)+z(84)+z(99)-3.0d0*z(79))
         ddz(99)=ddz(99)+(gam2*3.0d0/bndprm**2)*(z(55)+z(84)+z(99)-3.0d0*z(79))

         vpuck=vpuck+(z(60)+z(105)+z(107)-3.0d0*z(80))**2
         ddz(80)=ddz(80)-(gam2*9.0d0/bndprm**2)*(z(60)+z(105)+z(107)-3.0d0*z(80))
         ddz(60)=ddz(60)+(gam2*3.0d0/bndprm**2)*(z(60)+z(105)+z(107)-3.0d0*z(80))
         ddz(105)=ddz(105)+(gam2*3.0d0/bndprm**2)*(z(60)+z(105)+z(107)-3.0d0*z(80))
         ddz(107)=ddz(107)+(gam2*3.0d0/bndprm**2)*(z(60)+z(105)+z(107)-3.0d0*z(80))

         vpuck=vpuck+(z(59)+z(104)+z(106)-3.0d0*z(81))**2
         ddz(81)=ddz(81)-(gam2*9.0d0/bndprm**2)*(z(59)+z(104)+z(106)-3.0d0*z(81))
         ddz(59)=ddz(59)+(gam2*3.0d0/bndprm**2)*(z(59)+z(104)+z(106)-3.0d0*z(81))
         ddz(104)=ddz(104)+(gam2*3.0d0/bndprm**2)*(z(59)+z(104)+z(106)-3.0d0*z(81))
         ddz(106)=ddz(106)+(gam2*3.0d0/bndprm**2)*(z(59)+z(104)+z(106)-3.0d0*z(81))

         vpuck=vpuck+(z(61)+z(108)+z(109)-3.0d0*z(85))**2
         ddz(85)=ddz(85)-(gam2*9.0d0/bndprm**2)*(z(61)+z(108)+z(109)-3.0d0*z(85))
         ddz(61)=ddz(61)+(gam2*3.0d0/bndprm**2)*(z(61)+z(108)+z(109)-3.0d0*z(85))
         ddz(108)=ddz(108)+(gam2*3.0d0/bndprm**2)*(z(61)+z(108)+z(109)-3.0d0*z(85))
         ddz(109)=ddz(109)+(gam2*3.0d0/bndprm**2)*(z(61)+z(108)+z(109)-3.0d0*z(85))

         vpuck=vpuck+(z(62)+z(96)+z(112)-3.0d0*z(86))**2
         ddz(86)=ddz(86)-(gam2*9.0d0/bndprm**2)*(z(62)+z(96)+z(112)-3.0d0*z(86))
         ddz(62)=ddz(62)+(gam2*3.0d0/bndprm**2)*(z(62)+z(96)+z(112)-3.0d0*z(86))
         ddz(96)=ddz(96)+(gam2*3.0d0/bndprm**2)*(z(62)+z(96)+z(112)-3.0d0*z(86))
         ddz(112)=ddz(112)+(gam2*3.0d0/bndprm**2)*(z(62)+z(96)+z(112)-3.0d0*z(86))

         vpuck=vpuck+(z(63)+z(97)+z(113)-3.0d0*z(87))**2
         ddz(87)=ddz(87)-(gam2*9.0d0/bndprm**2)*(z(63)+z(97)+z(113)-3.0d0*z(87))
         ddz(63)=ddz(63)+(gam2*3.0d0/bndprm**2)*(z(63)+z(97)+z(113)-3.0d0*z(87))
         ddz(97)=ddz(97)+(gam2*3.0d0/bndprm**2)*(z(63)+z(97)+z(113)-3.0d0*z(87))
         ddz(113)=ddz(113)+(gam2*3.0d0/bndprm**2)*(z(63)+z(97)+z(113)-3.0d0*z(87))

         vpuck=vpuck+(z(65)+z(98)+z(114)-3.0d0*z(88))**2
         ddz(88)=ddz(88)-(gam2*9.0d0/bndprm**2)*(z(65)+z(98)+z(114)-3.0d0*z(88))
         ddz(65)=ddz(65)+(gam2*3.0d0/bndprm**2)*(z(65)+z(98)+z(114)-3.0d0*z(88))
         ddz(98)=ddz(98)+(gam2*3.0d0/bndprm**2)*(z(65)+z(98)+z(114)-3.0d0*z(88))
         ddz(114)=ddz(114)+(gam2*3.0d0/bndprm**2)*(z(65)+z(98)+z(114)-3.0d0*z(88))

         vpuck=vpuck+(z(66)+z(99)+z(115)-3.0d0*z(89))**2
         ddz(89)=ddz(89)-(gam2*9.0d0/bndprm**2)*(z(66)+z(99)+z(115)-3.0d0*z(89))
         ddz(66)=ddz(66)+(gam2*3.0d0/bndprm**2)*(z(66)+z(99)+z(115)-3.0d0*z(89))
         ddz(99)=ddz(99)+(gam2*3.0d0/bndprm**2)*(z(66)+z(99)+z(115)-3.0d0*z(89))
         ddz(115)=ddz(115)+(gam2*3.0d0/bndprm**2)*(z(66)+z(99)+z(115)-3.0d0*z(89))

         vpuck=vpuck+(z(64)+z(100)+z(110)-3.0d0*z(90))**2
         ddz(90)=ddz(90)-(gam2*9.0d0/bndprm**2)*(z(64)+z(100)+z(110)-3.0d0*z(90))
         ddz(64)=ddz(64)+(gam2*3.0d0/bndprm**2)*(z(64)+z(100)+z(110)-3.0d0*z(90))
         ddz(100)=ddz(100)+(gam2*3.0d0/bndprm**2)*(z(64)+z(100)+z(110)-3.0d0*z(90))
         ddz(110)=ddz(110)+(gam2*3.0d0/bndprm**2)*(z(64)+z(100)+z(110)-3.0d0*z(90))

         vpuck=vpuck+(z(67)+z(101)+z(111)-3.0d0*z(91))**2
         ddz(91)=ddz(91)-(gam2*9.0d0/bndprm**2)*(z(67)+z(101)+z(111)-3.0d0*z(91))
         ddz(67)=ddz(67)+(gam2*3.0d0/bndprm**2)*(z(67)+z(101)+z(111)-3.0d0*z(91))
         ddz(101)=ddz(101)+(gam2*3.0d0/bndprm**2)*(z(67)+z(101)+z(111)-3.0d0*z(91))
         ddz(111)=ddz(111)+(gam2*3.0d0/bndprm**2)*(z(67)+z(101)+z(111)-3.0d0*z(91))

         vpuck=vpuck+(z(70)+z(104)+z(120)-3.0d0*z(92))**2
         ddz(92)=ddz(92)-(gam2*9.0d0/bndprm**2)*(z(70)+z(104)+z(120)-3.0d0*z(92))
         ddz(70)=ddz(70)+(gam2*3.0d0/bndprm**2)*(z(70)+z(104)+z(120)-3.0d0*z(92))
         ddz(104)=ddz(104)+(gam2*3.0d0/bndprm**2)*(z(70)+z(104)+z(120)-3.0d0*z(92))
         ddz(120)=ddz(120)+(gam2*3.0d0/bndprm**2)*(z(70)+z(104)+z(120)-3.0d0*z(92))

         vpuck=vpuck+(z(73)+z(105)+z(121)-3.0d0*z(93))**2
         ddz(93)=ddz(93)-(gam2*9.0d0/bndprm**2)*(z(73)+z(105)+z(121)-3.0d0*z(93))
         ddz(73)=ddz(73)+(gam2*3.0d0/bndprm**2)*(z(73)+z(105)+z(121)-3.0d0*z(93))
         ddz(105)=ddz(105)+(gam2*3.0d0/bndprm**2)*(z(73)+z(105)+z(121)-3.0d0*z(93))
         ddz(121)=ddz(121)+(gam2*3.0d0/bndprm**2)*(z(73)+z(105)+z(121)-3.0d0*z(93))

         vpuck=vpuck+(z(68)+z(106)+z(116)-3.0d0*z(94))**2
         ddz(94)=ddz(94)-(gam2*9.0d0/bndprm**2)*(z(68)+z(106)+z(116)-3.0d0*z(94))
         ddz(68)=ddz(68)+(gam2*3.0d0/bndprm**2)*(z(68)+z(106)+z(116)-3.0d0*z(94))
         ddz(106)=ddz(106)+(gam2*3.0d0/bndprm**2)*(z(68)+z(106)+z(116)-3.0d0*z(94))
         ddz(116)=ddz(116)+(gam2*3.0d0/bndprm**2)*(z(68)+z(106)+z(116)-3.0d0*z(94))

         vpuck=vpuck+(z(69)+z(107)+z(117)-3.0d0*z(95))**2
         ddz(95)=ddz(95)-(gam2*9.0d0/bndprm**2)*(z(69)+z(107)+z(117)-3.0d0*z(95))
         ddz(69)=ddz(69)+(gam2*3.0d0/bndprm**2)*(z(69)+z(107)+z(117)-3.0d0*z(95))
         ddz(107)=ddz(107)+(gam2*3.0d0/bndprm**2)*(z(69)+z(107)+z(117)-3.0d0*z(95))
         ddz(117)=ddz(117)+(gam2*3.0d0/bndprm**2)*(z(69)+z(107)+z(117)-3.0d0*z(95))

         vpuck=vpuck+(z(71)+z(108)+z(118)-3.0d0*z(102))**2
         ddz(102)=ddz(102)-(gam2*9.0d0/bndprm**2)*(z(71)+z(108)+z(118)-3.0d0*z(102))
         ddz(71)=ddz(71)+(gam2*3.0d0/bndprm**2)*(z(71)+z(108)+z(118)-3.0d0*z(102))
         ddz(108)=ddz(108)+(gam2*3.0d0/bndprm**2)*(z(71)+z(108)+z(118)-3.0d0*z(102))
         ddz(118)=ddz(118)+(gam2*3.0d0/bndprm**2)*(z(71)+z(108)+z(118)-3.0d0*z(102))

         vpuck=vpuck+(z(72)+z(109)+z(119)-3.0d0*z(103))**2
         ddz(103)=ddz(103)-(gam2*9.0d0/bndprm**2)*(z(72)+z(109)+z(119)-3.0d0*z(103))
         ddz(72)=ddz(72)+(gam2*3.0d0/bndprm**2)*(z(72)+z(109)+z(119)-3.0d0*z(103))
         ddz(109)=ddz(109)+(gam2*3.0d0/bndprm**2)*(z(72)+z(109)+z(119)-3.0d0*z(103))
         ddz(119)=ddz(119)+(gam2*3.0d0/bndprm**2)*(z(72)+z(109)+z(119)-3.0d0*z(103))

         vpuck=vpuck*0.5d0*gam2*3.0d0/bndprm**2

         vlatt=vpuck+vtors


   ! **************************************************************************************************
   !            COMBINE LATTICE AND C-H SYSTEM
   ! **************************************************************************************************

         ! Total Potential 
         vv=vt+vlatt-0.5*rkc*(z(1)-qqq)**2

         ! Convert total potential to AU
         vv = vv / MyConsts_Hartree2eV

         ! Forces (negative sign because F = -dV/dx )
         IF ( rho == 0.0 ) THEN
            Forces(1) = 0.0
            Forces(2) = 0.0
         ELSE
            Forces(1) = -dvtdrho*xh/rho
            Forces(2) = -dvtdrho*yh/rho
         ENDIF
         Forces(3)=-dvtds1
         
         DO nn = 4, NrNonFrozen
            IF ( nn == 4 ) THEN
               Forces(4)=-(dvtds2+ddz(1)-rkc*(z(1)-qqq))
!               Forces(4)= 0.0  ! TO FIX EVEN C1 ATOM
            ELSE IF ( nn == 5 ) THEN
               Forces(5)=-(-(dvtds1+dvtds2)/3.0+ddz(2)+rkc*(z(1)-qqq)/3.0)
            ELSE IF ( nn == 6 ) THEN
               Forces(6)=-(-(dvtds1+dvtds2)/3.0+ddz(3)+rkc*(z(1)-qqq)/3.0)
            ELSE IF ( nn == 7 ) THEN
               Forces(7)=-(-(dvtds1+dvtds2)/3.0+ddz(4)+rkc*(z(1)-qqq)/3.0)
            ELSE 
               Forces(nn) = -ddz(nn-3)
            END IF
         END DO

!          ! Accelerations 
!          axh=-dvvdxh/(rmh*ev_tu)
!          ayh=-dvvdyh/(rmh*ev_tu)
!          azh=-dvvdzh/(rmh*ev_tu)
!          DO nn = 1, 121
!             azc(nn)=-dvvdz(nn)/(rmc*ev_tu)
!          END DO 

         ! Transform forces in atomic units (from eV Ang^-1 to Hartree Bohr^-1) 
         Forces(:) = Forces(:) * MyConsts_Bohr2Ang / MyConsts_Hartree2eV

      END FUNCTION VHSticking
      
END MODULE PotentialModule

