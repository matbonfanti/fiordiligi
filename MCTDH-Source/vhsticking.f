!*******************************************************************************
!> H-Graphene potential by Jackson and coworkers.
!> @ref http://pubs.acs.org/doi/abs/10.1021/jp057136%2B [ J.Phys.Chem.B 2006,110,18811 ]
!> @ref other papers jackson
!*******************************************************************************     
 
      SUBROUTINE vhsticking( rhx, rhy, rhz, rcz ,vv ) 
      IMPLICIT NONE
      REAL*8  rhx, rhy, rhz, rcz, vv

      ! temporary variables to store positions in Ang units
      REAL*8 :: xh, yh, zh
      REAL*8, DIMENSION(121) :: z

      REAL*8, PARAMETER    :: MyConsts_Bohr2Ang       = 0.52917721092               !< Conversion factor from Bohr to Ang (from NIST reference) 
      REAL*8, PARAMETER    :: MyConsts_Hartree2eV     = 27.21138505                 !< Conversion factor from Hartree to eV (from NIST reference)

      REAL*8, PARAMETER :: di = 0.00775           !< potential well (eV)
      REAL*8, PARAMETER :: alphi = 0.954          !< curvature (Ang-1)
      REAL*8, PARAMETER :: ai = 4.01              !< equilibrium position (Ang)
      ! parameters to describe the behaviour along rho at large rho
      REAL*8, PARAMETER :: alp = 1.4              !< a_f in the text (Ang-1) 
      REAL*8, PARAMETER :: alp2 = 1.65            !< b_f in the text (Ang-1)
      REAL*8, PARAMETER :: ba = 2.9               !< c_f in the text (Ang-1)
      REAL*8, PARAMETER :: rhoa = 1.0             !< rho_f in the text (Ang)
      ! switching function from small rho to large rho
      REAL*8, PARAMETER :: bs = 6.0               !< b_s in the text (Ang-1)
      REAL*8, PARAMETER :: rhos = 0.53            !< rho_s in the text (Ang)

      REAL*8, PARAMETER :: delt = 2.651           !< twisting potential parameter delta (eV)
      REAL*8, PARAMETER :: gam2 = 1.958           !< puckering paramter gamma_2 (eV)
      REAL*8, PARAMETER :: bndprm = 2.4612        !< lattice constant (in Ang)

      REAL*8 :: rkc      

      REAL*8 :: a, b, c1, c2, d0

      REAL*8 :: dbdzc, dc1dzc, dc2dzc, dd0dzc, df2dr, dkrdzc, dkrdzh
      REAL*8 :: dswdr, dvdzc, dvdzh, dvidzc, dvidzh, dvqdr, dvqdzc
      REAL*8 :: dvqdzh, dvtdrho, dvtds1, dvtds2
      REAL*8 :: dzgdzc, dzmdzc

      REAL*8 :: fexp, ff1, ff2, ff3
      REAL*8 :: qqq, rho, rkrho, sub1, sub2, sw
      REAL*8 :: v, vi, vlatt, vpuck, vq, vt, vtors
      REAL*8 :: zg, zm

      INTEGER :: nn

      xh = rhx * MyConsts_Bohr2Ang
      yh = rhy * MyConsts_Bohr2Ang
      zh = rhz * MyConsts_Bohr2Ang
      z(1) = rcz * MyConsts_Bohr2Ang
      z(2)= 0.0000000000 * MyConsts_Bohr2Ang 
      z(3)= 0.0000000000 * MyConsts_Bohr2Ang 
      z(4)= 0.0000000000 * MyConsts_Bohr2Ang 
      z(5)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(6)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(7)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(8)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(9)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(10)=-0.1148251303 * MyConsts_Bohr2Ang 
      z(11)=-0.0747453192 * MyConsts_Bohr2Ang 
      z(12)=-0.0747453192 * MyConsts_Bohr2Ang 
      z(13)=-0.0747453192 * MyConsts_Bohr2Ang 
      z(14)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(15)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(16)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(17)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(18)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(19)=-0.2202278756 * MyConsts_Bohr2Ang 
      z(20)=-0.2348303476 * MyConsts_Bohr2Ang 
      z(21)=-0.1549387993 * MyConsts_Bohr2Ang 
      z(22)=-0.2348303476 * MyConsts_Bohr2Ang 
      z(23)=-0.1549387993 * MyConsts_Bohr2Ang 
      z(24)=-0.1549387993 * MyConsts_Bohr2Ang 
      z(25)=-0.2348303476 * MyConsts_Bohr2Ang 
      z(26)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(27)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(28)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(29)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(30)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(31)=-0.2139867766 * MyConsts_Bohr2Ang 
      z(32)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(33)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(34)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(35)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(36)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(37)=-0.1835004879 * MyConsts_Bohr2Ang 
      z(38)=-0.2393535431 * MyConsts_Bohr2Ang 
      z(39)=-0.2393535431 * MyConsts_Bohr2Ang 
      z(40)=-0.2393535431 * MyConsts_Bohr2Ang 
      z(41)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(42)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(43)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(44)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(45)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(46)=-0.2201206755 * MyConsts_Bohr2Ang 
      z(47)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(48)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(49)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(50)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(51)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(52)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(53)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(54)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(55)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(56)=-0.2259041611 * MyConsts_Bohr2Ang 
      z(57)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(58)=-0.1889262549 * MyConsts_Bohr2Ang 
      z(59)=-0.1806026811 * MyConsts_Bohr2Ang 
      z(60)=-0.1806026811 * MyConsts_Bohr2Ang 
      z(61)=-0.1806026811 * MyConsts_Bohr2Ang 
      z(62)=-0.1991261779 * MyConsts_Bohr2Ang 
      z(63)=-0.1991261779 * MyConsts_Bohr2Ang 
      z(64)=-0.1991261779 * MyConsts_Bohr2Ang 
      z(65)=-0.1991261779 * MyConsts_Bohr2Ang 
      z(66)=-0.1991261779 * MyConsts_Bohr2Ang
      z(67)=-0.1991261779 * MyConsts_Bohr2Ang
      z(68)=-0.1912863158 * MyConsts_Bohr2Ang
      z(69)=-0.1912863158 * MyConsts_Bohr2Ang
      z(70)=-0.1912863158 * MyConsts_Bohr2Ang
      z(71)=-0.1912863158 * MyConsts_Bohr2Ang
      z(72)=-0.1912863158 * MyConsts_Bohr2Ang
      z(73)=-0.1912863158 * MyConsts_Bohr2Ang
      z(74)=-0.2042714793 * MyConsts_Bohr2Ang
      z(75)=-0.2042714793 * MyConsts_Bohr2Ang
      z(76)=-0.2042714793 * MyConsts_Bohr2Ang
      z(77)=-0.2042714793 * MyConsts_Bohr2Ang
      z(78)=-0.2042714793 * MyConsts_Bohr2Ang
      z(79)=-0.2042714793 * MyConsts_Bohr2Ang
      z(80)=-0.1755612797 * MyConsts_Bohr2Ang
      z(81)=-0.1755612797 * MyConsts_Bohr2Ang
      z(82)=-0.1973088448 * MyConsts_Bohr2Ang
      z(83)=-0.1973088448 * MyConsts_Bohr2Ang
      z(84)=-0.1973088448 * MyConsts_Bohr2Ang
      z(85)=-0.1755612797 * MyConsts_Bohr2Ang
      z(86)=-0.1849896025 * MyConsts_Bohr2Ang
      z(87)=-0.1849896025 * MyConsts_Bohr2Ang
      z(88)=-0.1849896025 * MyConsts_Bohr2Ang
      z(89)=-0.1849896025 * MyConsts_Bohr2Ang
      z(90)=-0.1849896025 * MyConsts_Bohr2Ang
      z(91)=-0.1849896025 * MyConsts_Bohr2Ang
      z(92)=-0.1801311916 * MyConsts_Bohr2Ang
      z(93)=-0.1801311916 * MyConsts_Bohr2Ang
      z(94)=-0.1801311916 * MyConsts_Bohr2Ang
      z(95)=-0.1801311916 * MyConsts_Bohr2Ang
      z(96)=-0.1849768649 * MyConsts_Bohr2Ang
      z(97)=-0.1849768649 * MyConsts_Bohr2Ang
      z(98)=-0.1849768649 * MyConsts_Bohr2Ang
      z(99)=-0.1849768649 * MyConsts_Bohr2Ang
      z(100)=-0.1849768649 * MyConsts_Bohr2Ang
      z(101)=-0.1849768649 * MyConsts_Bohr2Ang
      z(102)=-0.1801311916 * MyConsts_Bohr2Ang
      z(103)=-0.1801311916 * MyConsts_Bohr2Ang
      z(104)=-0.1719524214 * MyConsts_Bohr2Ang
      z(105)=-0.1719524214 * MyConsts_Bohr2Ang
      z(106)=-0.1719524214 * MyConsts_Bohr2Ang
      z(107)=-0.1719524214 * MyConsts_Bohr2Ang
      z(108)=-0.1719524214 * MyConsts_Bohr2Ang
      z(109)=-0.1719524214 * MyConsts_Bohr2Ang
      z(110)=-0.1659361800 * MyConsts_Bohr2Ang
      z(111)=-0.1659361800 * MyConsts_Bohr2Ang
      z(112)=-0.1659361800 * MyConsts_Bohr2Ang
      z(113)=-0.1659361800 * MyConsts_Bohr2Ang
      z(114)=-0.1659361800 * MyConsts_Bohr2Ang
      z(115)=-0.1659361800 * MyConsts_Bohr2Ang
      z(116)=-0.1774915409 * MyConsts_Bohr2Ang
      z(117)=-0.1774915409 * MyConsts_Bohr2Ang
      z(118)=-0.1774915409 * MyConsts_Bohr2Ang
      z(119)=-0.1774915409 * MyConsts_Bohr2Ang
      z(120)=-0.1774915409 * MyConsts_Bohr2Ang
      z(121)=-0.1774915409 * MyConsts_Bohr2Ang

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

      rkc = (36.0*gam2+6.0*delt)/(bndprm**2)

      ! Compute the parameters for the morse + gaussian 
      ! functional form for the collinear potential

      ! D_0 is the morse depth (eV)
      d0=0.474801+0.9878257*sub2-1.3921499*sub2**2              
     & +0.028278*sub2**3-1.8879928*sub2**4                       
     & +0.11*exp(-8.0*(sub2-0.28)**2)
      dd0dzc=0.9878257-2.7842998*sub2+0.084834*sub2**2-         
     & 7.5519712*sub2**3+(                                      
     & -1.76*(sub2-0.28)*exp(-8.0*(sub2-0.28)**2))

      ! A is the morse curvature (Ang)
      a=2.276211

      ! Z_M is the morse equilibrium distance (Ang)
      zm=0.87447*(sub2)+1.17425
      dzmdzc=0.87447

      ! C_1 is the asympotic harmonic potential for C1 vibrations (eV)
      c1=0.5*rkc*(sub2)**2-0.00326
      dc1dzc=rkc*(sub2)

      ! C_2 is the gaussian height (eV)
      c2= 0.3090344*exp(-2.741813*(sub2-0.2619756))+         
     & 0.03113325*exp(3.1844857*(sub2-0.186741))+            
     & 0.02*exp(-20.0*(sub2-0.1)**2) 
      dc2dzc=-0.8473145*exp(-2.741813*(sub2-0.2619756))+        
     & 0.0991434*exp(3.1844857*(sub2-0.186741))+(                
     & -0.8*(sub2-0.1)*exp(-20.0*(sub2-0.1)**2))

      ! B is the gaussian curvature (Ang-2)
      b=4.00181*exp(1.25965*(sub2-0.58729)) 
      dbdzc=5.0408799*exp(1.25965*(sub2-0.58729))

      ! Z_G is the center of the gaussian (Ang)
      zg=1.99155*(sub2)+1.46095
      dzgdzc=1.99155

      ! Compute the potential and the derivatives of the 
      ! collinear potential V_0

      v=c1+(c1+d0)*(exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+       
     & c2*exp(-b*(sub1-zg)**2)

      dvdzh=-b*c2*(sub1-zg)*2.0*exp(-b*(sub1-zg)**2) +                  
     & (c1+d0)*2.0*(-a*exp(-2.0*a*(sub1-zm))+a*exp(-a*(sub1-zm)))

      dvdzc=dc1dzc+(dc1dzc+dd0dzc)*                             
     & (exp(-2.0*a*(sub1-zm))-2.0*exp(-a*(sub1-zm)))+            
     & (c1+d0)*(a*dzmdzc*2.0*exp(-2.0*a*(sub1-zm))+              
     & (-a*dzmdzc)*2.0*exp(-a*(sub1-zm)))+                       
     & dc2dzc*exp(-b*(sub1-zg)**2)+                              
     & (-c2*dbdzc*(sub1-zg)**2+c2*b*2.0*(sub1-zg)*dzgdzc)*       
     & exp(-b*(sub1-zg)**2)

      ! Compute the force constant (and derivatives) for small rho 
      ! potential rkrho(zh-q,zc-q)

      rkrho=3.866259*exp(-17.038588*(sub2-0.537621)**2+         
     & 0.312355*(sub2-0.537621)*(sub1-2.003753)-                 
     & 4.479864*(sub1-2.003753)**2)+                             
     & 4.317415*exp(-11.931770*(sub2-0.286858)**2+               
     & 18.540974*(sub2-0.286858)*(sub1-1.540947)-                
     & 14.537321*(sub1-1.540947)**2)

      dkrdzc=(-34.077176*(sub2-0.537621)+0.312355*              
     & (sub1-2.003753))                                          
     & *3.866259*exp(-17.038588*(sub2-0.537621)**2+              
     & 0.312355*(sub2-0.537621)*(sub1-2.003753)-                 
     & 4.479864*(sub1-2.003753)**2)+                             
     & (-23.86354*(sub2-0.286858)+18.540974*(sub1-1.540947))     
     & *4.317415*exp(-11.931770*(sub2-0.286858)**2+              
     & 18.540974*(sub2-0.286858)*(sub1-1.540947)-                
     & 14.537321*(sub1-1.540947)**2)

      dkrdzh=(0.312355*(sub2-0.537621)-8.959728*(sub1-2.003753))        
     & *3.866259*exp(-17.038588*(sub2-0.537621)**2+                      
     & 0.312355*(sub2-0.537621)*(sub1-2.003753)-                         
     & 4.479864*(sub1-2.003753)**2)+                                     
     & (18.540974*(sub2-0.286858)-29.074642*(sub1-1.540947))             
     & *4.317415*exp(-11.931770*(sub2-0.286858)**2+                      
     & 18.540974*(sub2-0.286858)*(sub1-1.540947)-                        
     & 14.537321*(sub1-1.540947)**2)

      ! Compute the small rho potential/derivatives

      vq=v+0.5*rkrho*rho**2
      dvqdzh=dvdzh+0.5*dkrdzh*rho**2
      dvqdzc=dvdzc+0.5*dkrdzc*rho**2
      dvqdr=rkrho*rho

      ! Compute the  "infinite" rho potential/derivatives

      vi=0.5*rkc*sub2**2-0.00326+                              
     & di*(exp(-2.0*alphi*(sub1-ai))-2.0*exp(-alphi*(sub1-ai)))
      dvidzh=di*(-2.0*alphi*exp(-2.0*alphi*(sub1-ai))+         
     &       2.0*alphi*exp(-alphi*(sub1-ai)))
      dvidzc=rkc*sub2

      ! Switching function and associated functions

      fexp=exp(-ba*(rho-rhoa))
      ff1=1.0+fexp
      ff2=exp(-2.0*alp*rho)-2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)
      ff3=(vi-v)*ff2+vi
      sw=1.0/(1.0+exp(-bs*(rho-rhos)))

      ! Total H,C1 potential/derivatives

      vt=vq*(1.0-sw)+ff3*sw
   
      df2dr=-2.0*alp*exp(-2.0*alp*rho)-               
     &    2.0*exp(-alp*rho)*exp(-alp2*rho/ff1)*       
     &    (-alp-(alp2/ff1)-alp2*rho*ba*fexp/(ff1**2))
      dswdr=(bs*exp(-bs*(rho-rhos)))/((1.0+exp(-bs*(rho-rhos)))**2)

      dvtds1=dvqdzh*(1.0-sw)+sw*((dvidzh-dvdzh)*ff2+dvidzh)
      dvtds2=dvqdzc*(1.0-sw)+sw*((dvidzc-dvdzc)*ff2+dvidzc)
      dvtdrho=dvqdr*(1.0-sw)+vq*(-dswdr)+sw*(vi-v)*df2dr+ff3*dswdr

! **************************************************************************************************
!                           POTENTIAL FOR THE GRAPHITE LATTICE
! **************************************************************************************************

      vpuck=0.0
      vtors=0.0

      vtors=vtors+(z(1)-z(5)-z(12)+z(15))**2       
      vtors=vtors+(z(1)-z(7)-z(13)+z(18))**2       
      vtors=vtors+(z(1)-z(9)-z(11)+z(16))**2       
      vtors=vtors+(z(2)-z(3)-z(8)+z(10))**2        
      vtors=vtors+(z(2)-z(11)-z(20)+z(27))**2      
      vtors=vtors+(z(2)-z(14)-z(7)+z(21))**2       
      vtors=vtors+(z(3)-z(11)-z(22)+z(28))**2      
      vtors=vtors+(z(3)-z(16)-z(5)+z(21))**2       
      vtors=vtors+(z(3)-z(18)-z(10)+z(24))**2      
      vtors=vtors+(z(4)-z(2)-z(9)+z(7))**2         
      vtors=vtors+(z(4)-z(3)-z(6)+z(5))**2         
      vtors=vtors+(z(4)-z(13)-z(25)+z(31))**2      
      vtors=vtors+(z(5)-z(7)-z(34)+z(32))**2       
      vtors=vtors+(z(5)-z(27)-z(15)+z(38))**2      
      vtors=vtors+(z(6)-z(1)-z(14)+z(11))**2       
      vtors=vtors+(z(6)-z(20)-z(35)+z(43))**2      
      vtors=vtors+(z(7)-z(22)-z(32)+z(45))**2      
      vtors=vtors+(z(7)-z(28)-z(18)+z(39))**2      
      vtors=vtors+(z(8)-z(1)-z(17)+z(13))**2       
      vtors=vtors+(z(8)-z(6)-z(33)+z(35))**2       
      vtors=vtors+(z(9)-z(22)-z(36)+z(46))**2      
      vtors=vtors+(z(9)-z(30)-z(16)+z(39))**2      
      vtors=vtors+(z(10)-z(1)-z(19)+z(12))**2      
      vtors=vtors+(z(10)-z(9)-z(37)+z(36))**2      
      vtors=vtors+(z(11)-z(32)-z(27)+z(52))**2     
      vtors=vtors+(z(11)-z(34)-z(28)+z(51))**2     
      vtors=vtors+(z(12)-z(2)-z(26)+z(20))**2      
      vtors=vtors+(z(12)-z(4)-z(29)+z(25))**2      
      vtors=vtors+(z(13)-z(3)-z(30)+z(22))**2      
      vtors=vtors+(z(13)-z(36)-z(31)+z(58))**2     
      vtors=vtors+(z(14)-z(34)-z(47)+z(62))**2                                                                                                                                   
      vtors=vtors+(z(14)-z(42)-z(21)+z(52))**2                                                                                                                                   
      vtors=vtors+(z(15)-z(2)-z(23)+z(8))**2                                                                                                                                     
      vtors=vtors+(z(15)-z(14)-z(48)+z(47))**2                                                                                                                                   
      vtors=vtors+(z(16)-z(32)-z(56)+z(64))**2                                                                                                                                   
      vtors=vtors+(z(16)-z(45)-z(21)+z(51))**2                                                                                                                                   
      vtors=vtors+(z(17)-z(4)-z(24)+z(9))**2                                                                                                                                     
      vtors=vtors+(z(17)-z(37)-z(55)+z(66))**2                                                                                                                                   
      vtors=vtors+(z(18)-z(16)-z(54)+z(56))**2                                                                                                                                   
      vtors=vtors+(z(18)-z(46)-z(24)+z(57))**2                                                                                                                                   
      vtors=vtors+(z(19)-z(4)-z(23)+z(6))**2                                                                                                                                     
      vtors=vtors+(z(19)-z(17)-z(50)+z(55))**2                                                                                                                                   
      vtors=vtors+(z(20)-z(5)-z(42)+z(34))**2                                                                                                                                    
      vtors=vtors+(z(20)-z(47)-z(43)+z(76))**2                                                                                                                                   
      vtors=vtors+(z(21)-z(28)-z(59)+z(70))**2                                                                                                                                   
      vtors=vtors+(z(22)-z(54)-z(45)+z(77))**2                                                                                                                                   
      vtors=vtors+(z(22)-z(56)-z(46)+z(78))**2                                                                                                                                   
      vtors=vtors+(z(23)-z(26)-z(60)+z(69))**2                                                                                                                                   
      vtors=vtors+(z(24)-z(30)-z(61)+z(71))**2                                                                                                                                   
      vtors=vtors+(z(25)-z(8)-z(41)+z(33))**2                                                                                                                                    
      vtors=vtors+(z(25)-z(10)-z(44)+z(37))**2                                                                                                                                   
      vtors=vtors+(z(26)-z(6)-z(38)+z(14))**2                                                                                                                                    
      vtors=vtors+(z(26)-z(48)-z(69)+z(87))**2                                                                                                                                   
      vtors=vtors+(z(27)-z(21)-z(68)+z(59))**2                                                                                                                                   
      vtors=vtors+(z(27)-z(62)-z(38)+z(75))**2                                                                                                                                   
      vtors=vtors+(z(28)-z(56)-z(70)+z(90))**2                                                                                                                                   
      vtors=vtors+(z(28)-z(64)-z(39)+z(77))**2                                                                                                                                   
      vtors=vtors+(z(29)-z(8)-z(40)+z(17))**2                                                                                                                                    
      vtors=vtors+(z(29)-z(23)-z(73)+z(60))**2                                                                                                                                   
      vtors=vtors+(z(30)-z(54)-z(71)+z(88))**2                                                                                                                                   
      vtors=vtors+(z(30)-z(65)-z(39)+z(78))**2                                                                                                                                   
      vtors=vtors+(z(31)-z(10)-z(40)+z(19))**2                                                                                                                                   
      vtors=vtors+(z(31)-z(24)-z(72)+z(61))**2                                                                                                                                   
      vtors=vtors+(z(32)-z(59)-z(64)+z(92))**2                                                                                                                                   
      vtors=vtors+(z(32)-z(70)-z(52)+z(81))**2                                                                                                                                   
      vtors=vtors+(z(33)-z(12)-z(53)+z(26))**2
      vtors=vtors+(z(33)-z(19)-z(67)+z(50))**2
      vtors=vtors+(z(34)-z(59)-z(62)+z(94))**2
      vtors=vtors+(z(34)-z(68)-z(51)+z(81))**2
      vtors=vtors+(z(35)-z(12)-z(49)+z(29))**2
      vtors=vtors+(z(35)-z(15)-z(63)+z(48))**2
      vtors=vtors+(z(36)-z(18)-z(65)+z(54))**2
      vtors=vtors+(z(36)-z(71)-z(58)+z(85))**2
      vtors=vtors+(z(37)-z(13)-z(57)+z(30))**2
      vtors=vtors+(z(37)-z(61)-z(66)+z(103))**2
      vtors=vtors+(z(38)-z(42)-z(82)+z(96))**2
      vtors=vtors+(z(39)-z(45)-z(83)+z(100))**2
      vtors=vtors+(z(40)-z(44)-z(84)+z(99))**2
      vtors=vtors+(z(41)-z(19)-z(49)+z(23))**2
      vtors=vtors+(z(41)-z(40)-z(101)+z(84))**2
      vtors=vtors+(z(42)-z(68)-z(96)+z(112))**2
      vtors=vtors+(z(42)-z(86)-z(52)+z(94))**2
      vtors=vtors+(z(43)-z(15)-z(53)+z(23))**2
      vtors=vtors+(z(43)-z(38)-z(97)+z(82))**2
      vtors=vtors+(z(44)-z(17)-z(58)+z(24))**2
      vtors=vtors+(z(44)-z(72)-z(99)+z(115))**2
      vtors=vtors+(z(45)-z(70)-z(100)+z(110))**2
      vtors=vtors+(z(45)-z(90)-z(51)+z(92))**2
      vtors=vtors+(z(46)-z(39)-z(98)+z(83))**2
      vtors=vtors+(z(46)-z(88)-z(57)+z(102))**2
      vtors=vtors+(z(47)-z(27)-z(86)+z(68))**2
      vtors=vtors+(z(48)-z(20)-z(75)+z(42))**2
      vtors=vtors+(z(49)-z(53)-z(105)+z(107))**2
      vtors=vtors+(z(50)-z(25)-z(79)+z(44))**2
      vtors=vtors+(z(50)-z(29)-z(91)+z(73))**2
      vtors=vtors+(z(51)-z(64)-z(104)+z(120))**2
      vtors=vtors+(z(52)-z(51)-z(106)+z(104))**2
      vtors=vtors+(z(53)-z(63)-z(107)+z(117))**2
      vtors=vtors+(z(55)-z(25)-z(74)+z(41))**2
      vtors=vtors+(z(55)-z(31)-z(89)+z(72))**2
      vtors=vtors+(z(57)-z(65)-z(108)+z(118))**2
      vtors=vtors+(z(58)-z(57)-z(109)+z(108))**2
      vtors=vtors+(z(60)-z(33)-z(93)+z(67))**2
      vtors=vtors+(z(60)-z(35)-z(95)+z(63))**2
      vtors=vtors+(z(61)-z(36)-z(102)+z(65))**2
      vtors=vtors+(z(62)-z(52)-z(116)+z(106))**2
      vtors=vtors+(z(63)-z(26)-z(76)+z(38))**2
      vtors=vtors+(z(66)-z(31)-z(79)+z(40))**2
      vtors=vtors+(z(66)-z(58)-z(119)+z(109))**2
      vtors=vtors+(z(67)-z(29)-z(74)+z(40))**2
      vtors=vtors+(z(67)-z(49)-z(121)+z(105))**2
      vtors=vtors+(z(69)-z(35)-z(80)+z(49))**2
      vtors=vtors+(z(69)-z(43)-z(113)+z(97))**2
      vtors=vtors+(z(71)-z(46)-z(114)+z(98))**2
      vtors=vtors+(z(72)-z(37)-z(85)+z(57))**2
      vtors=vtors+(z(73)-z(33)-z(80)+z(53))**2
      vtors=vtors+(z(73)-z(41)-z(111)+z(101))**2
      vtors=vtors+(z(87)-z(43)-z(95)+z(53))**2
      vtors=vtors+(z(89)-z(44)-z(103)+z(58))**2
      vtors=vtors+(z(91)-z(41)-z(93)+z(49))**2

      vtors=vtors*0.5d0*delt/bndprm**2

      vpuck=vpuck+(z(2)+z(3)+z(4)-3.0d0*z(1))**2   
      vpuck=vpuck+(z(1)+z(5)+z(6)-3.0d0*z(2))**2   
      vpuck=vpuck+(z(1)+z(7)+z(9)-3.0d0*z(3))**2   
      vpuck=vpuck+(z(1)+z(8)+z(10)-3.0d0*z(4))**2  
      vpuck=vpuck+(z(2)+z(11)+z(14)-3.0d0*z(5))**2 
      vpuck=vpuck+(z(2)+z(12)+z(15)-3.0d0*z(6))**2 
      vpuck=vpuck+(z(3)+z(11)+z(16)-3.0d0*z(7))**2 
      vpuck=vpuck+(z(4)+z(12)+z(19)-3.0d0*z(8))**2 
      vpuck=vpuck+(z(3)+z(13)+z(18)-3.0d0*z(9))**2 
      vpuck=vpuck+(z(4)+z(13)+z(17)-3.0d0*z(10))**2
      vpuck=vpuck+(z(5)+z(7)+z(21)-3.0d0*z(11))**2 
      vpuck=vpuck+(z(6)+z(8)+z(23)-3.0d0*z(12))**2 
      vpuck=vpuck+(z(9)+z(10)+z(24)-3.0d0*z(13))**2
      vpuck=vpuck+(z(5)+z(20)+z(27)-3.0d0*z(14))**2
      vpuck=vpuck+(z(6)+z(20)+z(26)-3.0d0*z(15))**2
      vpuck=vpuck+(z(7)+z(22)+z(28)-3.0d0*z(16))**2
      vpuck=vpuck+(z(10)+z(25)+z(31)-3.0d0*z(17))**2
      vpuck=vpuck+(z(9)+z(22)+z(30)-3.0d0*z(18))**2 
      vpuck=vpuck+(z(8)+z(25)+z(29)-3.0d0*z(19))**2 
      vpuck=vpuck+(z(14)+z(15)+z(38)-3.0d0*z(20))**2
      vpuck=vpuck+(z(11)+z(32)+z(34)-3.0d0*z(21))**2
      vpuck=vpuck+(z(16)+z(18)+z(39)-3.0d0*z(22))**2
      vpuck=vpuck+(z(12)+z(33)+z(35)-3.0d0*z(23))**2
      vpuck=vpuck+(z(13)+z(36)+z(37)-3.0d0*z(24))**2
      vpuck=vpuck+(z(17)+z(19)+z(40)-3.0d0*z(25))**2
      vpuck=vpuck+(z(15)+z(35)+z(43)-3.0d0*z(26))**2
      vpuck=vpuck+(z(14)+z(34)+z(42)-3.0d0*z(27))**2
      vpuck=vpuck+(z(16)+z(32)+z(45)-3.0d0*z(28))**2
      vpuck=vpuck+(z(19)+z(33)+z(41)-3.0d0*z(29))**2
      vpuck=vpuck+(z(18)+z(36)+z(46)-3.0d0*z(30))**2
      vpuck=vpuck+(z(17)+z(37)+z(44)-3.0d0*z(31))**2
      vpuck=vpuck+(z(21)+z(28)+z(51)-3.0d0*z(32))**2
      vpuck=vpuck+(z(23)+z(29)+z(49)-3.0d0*z(33))**2
      vpuck=vpuck+(z(21)+z(27)+z(52)-3.0d0*z(34))**2
      vpuck=vpuck+(z(23)+z(26)+z(53)-3.0d0*z(35))**2
      vpuck=vpuck+(z(24)+z(30)+z(57)-3.0d0*z(36))**2
      vpuck=vpuck+(z(24)+z(31)+z(58)-3.0d0*z(37))**2
      vpuck=vpuck+(z(20)+z(47)+z(48)-3.0d0*z(38))**2
      vpuck=vpuck+(z(22)+z(54)+z(56)-3.0d0*z(39))**2
      vpuck=vpuck+(z(25)+z(50)+z(55)-3.0d0*z(40))**2
      vpuck=vpuck+(z(29)+z(50)+z(67)-3.0d0*z(41))**2
      vpuck=vpuck+(z(27)+z(47)+z(62)-3.0d0*z(42))**2
      vpuck=vpuck+(z(26)+z(48)+z(63)-3.0d0*z(43))**2
      vpuck=vpuck+(z(31)+z(55)+z(66)-3.0d0*z(44))**2
      vpuck=vpuck+(z(28)+z(56)+z(64)-3.0d0*z(45))**2
      vpuck=vpuck+(z(30)+z(54)+z(65)-3.0d0*z(46))**2
      vpuck=vpuck+(z(38)+z(42)+z(75)-3.0d0*z(47))**2
      vpuck=vpuck+(z(38)+z(43)+z(76)-3.0d0*z(48))**2
      vpuck=vpuck+(z(33)+z(60)+z(73)-3.0d0*z(49))**2
      vpuck=vpuck+(z(40)+z(41)+z(74)-3.0d0*z(50))**2
      vpuck=vpuck+(z(32)+z(59)+z(70)-3.0d0*z(51))**2
      vpuck=vpuck+(z(34)+z(59)+z(68)-3.0d0*z(52))**2
      vpuck=vpuck+(z(35)+z(60)+z(69)-3.0d0*z(53))**2
      vpuck=vpuck+(z(39)+z(46)+z(78)-3.0d0*z(54))**2
      vpuck=vpuck+(z(40)+z(44)+z(79)-3.0d0*z(55))**2
      vpuck=vpuck+(z(39)+z(45)+z(77)-3.0d0*z(56))**2
      vpuck=vpuck+(z(36)+z(61)+z(71)-3.0d0*z(57))**2
      vpuck=vpuck+(z(37)+z(61)+z(72)-3.0d0*z(58))**2
      vpuck=vpuck+(z(51)+z(52)+z(81)-3.0d0*z(59))**2
      vpuck=vpuck+(z(49)+z(53)+z(80)-3.0d0*z(60))**2
      vpuck=vpuck+(z(57)+z(58)+z(85)-3.0d0*z(61))**2
      vpuck=vpuck+(z(42)+z(68)+z(86)-3.0d0*z(62))**2
      vpuck=vpuck+(z(43)+z(69)+z(87)-3.0d0*z(63))**2
      vpuck=vpuck+(z(45)+z(70)+z(90)-3.0d0*z(64))**2
      vpuck=vpuck+(z(46)+z(71)+z(88)-3.0d0*z(65))**2
      vpuck=vpuck+(z(44)+z(72)+z(89)-3.0d0*z(66))**2
      vpuck=vpuck+(z(41)+z(73)+z(91)-3.0d0*z(67))**2
      vpuck=vpuck+(z(52)+z(62)+z(94)-3.0d0*z(68))**2
      vpuck=vpuck+(z(53)+z(63)+z(95)-3.0d0*z(69))**2
      vpuck=vpuck+(z(51)+z(64)+z(92)-3.0d0*z(70))**2
      vpuck=vpuck+(z(57)+z(65)+z(102)-3.0d0*z(71))**2
      vpuck=vpuck+(z(58)+z(66)+z(103)-3.0d0*z(72))**2
      vpuck=vpuck+(z(49)+z(67)+z(93)-3.0d0*z(73))**2
      vpuck=vpuck+(z(50)+z(84)+z(101)-3.0d0*z(74))**2
      vpuck=vpuck+(z(47)+z(82)+z(96)-3.0d0*z(75))**2
      vpuck=vpuck+(z(48)+z(82)+z(97)-3.0d0*z(76))**2
      vpuck=vpuck+(z(56)+z(83)+z(100)-3.0d0*z(77))**2
      vpuck=vpuck+(z(54)+z(83)+z(98)-3.0d0*z(78))**2
      vpuck=vpuck+(z(55)+z(84)+z(99)-3.0d0*z(79))**2
      vpuck=vpuck+(z(60)+z(105)+z(107)-3.0d0*z(80))**2
      vpuck=vpuck+(z(59)+z(104)+z(106)-3.0d0*z(81))**2
      vpuck=vpuck+(z(61)+z(108)+z(109)-3.0d0*z(85))**2
      vpuck=vpuck+(z(62)+z(96)+z(112)-3.0d0*z(86))**2
      vpuck=vpuck+(z(63)+z(97)+z(113)-3.0d0*z(87))**2
      vpuck=vpuck+(z(65)+z(98)+z(114)-3.0d0*z(88))**2
      vpuck=vpuck+(z(66)+z(99)+z(115)-3.0d0*z(89))**2
      vpuck=vpuck+(z(64)+z(100)+z(110)-3.0d0*z(90))**2
      vpuck=vpuck+(z(67)+z(101)+z(111)-3.0d0*z(91))**2
      vpuck=vpuck+(z(70)+z(104)+z(120)-3.0d0*z(92))**2
      vpuck=vpuck+(z(73)+z(105)+z(121)-3.0d0*z(93))**2
      vpuck=vpuck+(z(68)+z(106)+z(116)-3.0d0*z(94))**2
      vpuck=vpuck+(z(69)+z(107)+z(117)-3.0d0*z(95))**2
      vpuck=vpuck+(z(71)+z(108)+z(118)-3.0d0*z(102))**2
      vpuck=vpuck+(z(72)+z(109)+z(119)-3.0d0*z(103))**2

      vpuck=vpuck*0.5d0*gam2*3.0d0/bndprm**2

      vlatt=vpuck+vtors

! **************************************************************************************************
!            COMBINE LATTICE AND C-H SYSTEM
! **************************************************************************************************

!     Total Potential 
      vv=vt+vlatt-0.5*rkc*(z(1)-qqq)**2

!     Convert total potential to AU
      vv = vv / MyConsts_Hartree2eV

      END SUBROUTINE vhsticking
