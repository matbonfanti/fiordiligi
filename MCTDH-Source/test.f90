PROGRAM TestMCTDHFunction
   IMPLICIT NONE

   REAL*8 :: vC, vH, coup
   REAL*8 :: zh, zc, rho
   INTEGER :: i,j
      
   zc = 0.0d0
   rho = 0.0d0

   DO j = 0, 0
      zc = 0.8 + j*0.1d0

      print*, "# zc = ", zC*0.52917721092d0

      DO i = 0, 100
         zh = 1.0d0 + i*0.1d0

         CALL hstick_carbon( zc, vC ) 
         CALL hstick_hydro( rho, 0.0d0, zh, vH ) 
         CALL hstick_coupling( rho, 0.0d0, zh, zc, coup ) 
         PRINT*, zh, (vC+vH+coup)*27.21138505d0, VHTrapping( 0.0d0, 0.0d0, zh, zc )*27.21138505d0

      END DO

      PRINT*, " "
      PRINT*, " "
      PRINT*, " "

   END DO

   CONTAINS

      REAL*8 FUNCTION VHTrapping( xhbis, yhbis, zhbis, zcbis ) RESULT(vv) 
         IMPLICIT NONE
         REAL*8, INTENT(IN)  :: xhbis, yhbis, zhbis, zcbis 
         REAL*8  :: xh, yh, zh, zc 

      REAL*8, PARAMETER    :: MyConsts_Bohr2Ang       = 0.52917721092d0  !< Conversion factor from Bohr to Ang (from NIST )      
      REAL*8, PARAMETER    :: MyConsts_Hartree2eV     = 27.21138505d0    !< Conversion factor from Hartree to eV (from NIST )

         REAL*8, PARAMETER :: di = 0.00775           !< potential well (eV)
         REAL*8, PARAMETER :: alphi = 0.954          !< curvature (Ang-1)
         REAL*8, PARAMETER :: ai = 4.01              !< equilibrium position (Ang)
         REAL*8, PARAMETER :: alp = 1.4              !< a_f in the text (Ang-1) 
         REAL*8, PARAMETER :: alp2 = 1.65            !< b_f in the text (Ang-1)
         REAL*8, PARAMETER :: ba = 2.9               !< c_f in the text (Ang-1)
         REAL*8, PARAMETER :: rhoa = 1.0             !< rho_f in the text (Ang)
         REAL*8, PARAMETER :: bs = 6.0               !< b_s in the text (Ang-1)
         REAL*8, PARAMETER :: rhos = 0.53            !< rho_s in the text (Ang)
         REAL*8, PARAMETER :: delt = 2.651           !< twisting potential parameter delta (eV)
         REAL*8, PARAMETER :: gam2 = 1.958           !< puckering paramter gamma_2 (eV)
         REAL*8, PARAMETER :: bndprm = 2.4612        !< lattice constant (in Ang)
         REAL*8 :: rkc                               !< force constant for single carbon displacement, setup at runtime


         REAL*8 :: a, b, c1, c2, d0

         REAL*8 :: dbdzc, dc1dzc, dc2dzc, dd0dzc, df2dr, dkrdzc, dkrdzh
         REAL*8 :: dswdr, dvdzc, dvdzh, dvidzc, dvidzh, dvqdr, dvqdzc
         REAL*8 :: dvqdzh, dvtdrho, dvtds1, dvtds2
         REAL*8 :: dzgdzc, dzmdzc

         REAL*8 :: fexp, ff1, ff2, ff3
         REAL*8 :: rho, rkrho, sub1, sub2, sw
         REAL*8 :: v, vi, vq, vt
         REAL*8 :: zg, zm

         INTEGER :: nn

         xh = xhbis * MyConsts_Bohr2Ang
         yh = yhbis * MyConsts_Bohr2Ang
         zh = zhbis * MyConsts_Bohr2Ang
         zc = zcbis * MyConsts_Bohr2Ang

         ! Compute some relevant coordinates for the PES
         rkc = (36.0*gam2+6.0*delt)/(bndprm**2)
         ! positions of H and C1 with respect to zero
         sub1=zh
         sub2=zc
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

         ! Total Potential 
         vv=vt

         ! Upper and lower energy cutoff
         if ( vv < -1.0 ) vv = -1.0
         if ( vv > 20.0 ) vv = 20.0
      vv = vv / MyConsts_Hartree2eV

      END FUNCTION VHTrapping


END PROGRAM