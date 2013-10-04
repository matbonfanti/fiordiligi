!***************************************************************************************
!*                              DifferenceDerivatives
!***************************************************************************************
!*
!*        CONTAINS:          ThreePointDifferenceDerivative( Function, Point, Delta )
!*                           FivePointDifferenceDerivative( Function, Point, Delta )
!*                           FivePointDifferenceSecondOrderDerivative( Function, Point, Delta )
!*                           FivePointDifferenceThirdOrderDerivative( Function, Point, Delta )
!*                           FivePointDifferenceFourthOrderDerivative( Function, Point, Delta )
!*
!*                           DifferenceGradient( Function, DerivativeMethod, Point, Delta, Gradient )
!*
!*                           ThreePointDifferenceHessian( Function, Point, Delta, Hessian )
!*                           FivePointDifferenceHessian( Function, Point, Delta, Hessian )
!*
!***************************************************************************************
!*
!*        AUTHOR(S):          1st: M.F.Somers, dec 2004.
!*
!***************************************************************************************

MODULE DifferenceDerivatives
#include "preprocessoptions.cpp"
IMPLICIT NONE

PRIVATE AbsVector

CONTAINS

!***************************************************************************************
!* Absolute value of the vector ...

REAL FUNCTION AbsVector( Vector )
IMPLICIT NONE
REAL, DIMENSION(:) :: Vector

    AbsVector = SQRT( SUM( Vector(:) * Vector(:) ) )

END FUNCTION AbsVector

!***************************************************************************************
!* Calculate the three point difference derivative...

REAL FUNCTION ThreePointDifferenceDerivative( TheFunction, AtPoint, WithDelta )
IMPLICIT NONE
REAL, DIMENSION(:)         :: AtPoint, WithDelta
REAL                       :: Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(WithDelta),"DIFFERENCEDERIVATIVESMODULE.THREEPOINTDIFFERENCEDERIVATIVE:" // &
                                              " vectors do not match.")

    Factor = 1.0 / ( 2.0 * AbsVector( WithDelta ) )
    ThreePointDifferenceDerivative = Factor * ( TheFunction( AtPoint + WithDelta ) - &
                                                TheFunction( AtPoint - WithDelta ) )

END FUNCTION ThreePointDifferenceDerivative

!***************************************************************************************
!* Calculate the five point difference derivative...

REAL FUNCTION FivePointDifferenceDerivative( TheFunction, AtPoint, WithDelta )
IMPLICIT NONE
REAL, DIMENSION(:)         :: AtPoint, WithDelta
REAL                       :: Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(WithDelta),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCEDERIVATIVE:" // &
                                              " vectors do not match.")

    Factor = 1.0 / ( 12.0 * AbsVector( WithDelta ) )
    FivePointDifferenceDerivative = Factor * ( +1.0 * TheFunction( AtPoint - 2.0 * WithDelta ) +  &
                                               -8.0 * TheFunction( AtPoint - 1.0 * WithDelta ) +  &
                                               +8.0 * TheFunction( AtPoint + 1.0 * WithDelta ) +  &
                                               -1.0 * TheFunction( AtPoint + 2.0 * WithDelta ) )

END FUNCTION FivePointDifferenceDerivative

!***************************************************************************************
!* Calculate the five point difference second order derivative...

REAL FUNCTION FivePointDifferenceSecondOrderDerivative( TheFunction, AtPoint, WithDelta )
IMPLICIT NONE
REAL, DIMENSION(:)         :: AtPoint, WithDelta
REAL                       :: Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(WithDelta),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCESECONDORDERDERIVATIVE:" // &
                                              " vectors do not match.")

    Factor = 1.0 / ( 12.0 * AbsVector( WithDelta ) ** 2.0 )
    FivePointDifferenceSecondOrderDerivative = Factor * ( -1.00 * TheFunction( AtPoint - 2.0 * WithDelta ) +  &
                                                          +16.0 * TheFunction( AtPoint - 1.0 * WithDelta ) +  &
                                                          -30.0 * TheFunction( AtPoint ) +  &
                                                          +16.0 * TheFunction( AtPoint + 1.0 * WithDelta ) +  &
                                                          -1.00 * TheFunction( AtPoint + 2.0 * WithDelta ) )

END FUNCTION FivePointDifferenceSecondOrderDerivative

!***************************************************************************************
!* Calculate the five point difference third order derivative...

REAL FUNCTION FivePointDifferenceThirdOrderDerivative( TheFunction, AtPoint, WithDelta )
IMPLICIT NONE
REAL, DIMENSION(:)         :: AtPoint, WithDelta
REAL                       :: Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(WithDelta),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCETHIRDORDERDERIVATIVE:" // &
                                              " vectors do not match.")

    Factor = 1.0 / ( 2.0 * AbsVector( WithDelta ) ** 3.0 )
    FivePointDifferenceThirdOrderDerivative = Factor * ( -1.0 * TheFunction( AtPoint - 2.0 * WithDelta ) +  &
                                                         +2.0 * TheFunction( AtPoint - 1.0 * WithDelta ) +  &
                                                         -2.0 * TheFunction( AtPoint + 1.0 * WithDelta ) +  &
                                                         +1.0 * TheFunction( AtPoint + 2.0 * WithDelta ) )

END FUNCTION FivePointDifferenceThirdOrderDerivative

!***************************************************************************************
!* Calculate the five point difference fourth order derivative...

REAL FUNCTION FivePointDifferenceFourthOrderDerivative( TheFunction, AtPoint, WithDelta )
IMPLICIT NONE
REAL, DIMENSION(:)         :: AtPoint, WithDelta
REAL                       :: Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(WithDelta),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCEFOURTHORDERDERIVATIVE:" // &
                                              " vectors do not match.")

    Factor = 1.0 / ( AbsVector( WithDelta ) ** 4.0 )
    FivePointDifferenceFourthOrderDerivative = Factor * ( +1.0 * TheFunction( AtPoint - 2.0 * WithDelta ) +  &
                                                          -4.0 * TheFunction( AtPoint - 1.0 * WithDelta ) +  &
                                                          +6.0 * TheFunction( AtPoint ) +  &
                                                          -4.0 * TheFunction( AtPoint + 1.0 * WithDelta ) +  &
                                                          +1.0 * TheFunction( AtPoint + 2.0 * WithDelta ) )

END FUNCTION FivePointDifferenceFourthOrderDerivative

!***************************************************************************************
!* Calculates the gradient at given point...

SUBROUTINE DifferenceGradient( TheFunction, DerivativeMethod, AtPoint, Delta, Gradient )
IMPLICIT NONE
REAL, DIMENSION(:)                  :: Gradient, AtPoint
REAL, DIMENSION(1:SIZE(Gradient))   :: TempVector
INTEGER                             :: N
REAL                                :: Delta

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction

    REAL FUNCTION DerivativeMethod( TheFunction, AtPoint, WithDelta )
         REAL, DIMENSION(:) :: AtPoint, WithDelta
         INTERFACE
               REAL FUNCTION TheFunction( Point )
                   REAL, DIMENSION(:) :: Point
               END FUNCTION TheFunction
         END INTERFACE
    END FUNCTION DerivativeMethod
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(Gradient),"DIFFERENCEDERIVATIVESMODULE.DIFFERENCEGRADIENT:" // &
                                             " vectors do not match.")

    DO N = 1, SIZE( Gradient )
       TempVector = 0.0
       TempVector( N ) = Delta
       Gradient( N ) = DerivativeMethod( TheFunction, AtPoint, TempVector )
    END DO

END SUBROUTINE DifferenceGradient

!***************************************************************************************
!* Calculates the hessian at given point wtih three point diffrence method...

SUBROUTINE ThreePointDifferenceHessian( TheFunction, AtPoint, Delta, Hessian )
IMPLICIT NONE
REAL, DIMENSION(:,:)                :: Hessian
REAL, DIMENSION(:)                  :: AtPoint
REAL, DIMENSION(1:SIZE(Hessian,1))  :: TempVector, Gradient1, Gradient2
INTEGER                             :: N
REAL                                :: Delta, SqrtDelta, Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(Hessian,1),"DIFFERENCEDERIVATIVESMODULE.THREEPOINTDIFFERENCEHESSIAN:" // &
                                              " vectors do not match.")
    CALL ERROR(SIZE(AtPoint)/=SIZE(Hessian,2),"DIFFERENCEDERIVATIVESMODULE.THREEPOINTDIFFERENCEHESSIAN:" // &
                                              " vectors do not match.")

    SqrtDelta = SQRT( Delta )
    Factor = 1.0 / ( 2.0 * Delta )

    DO N = 1, SIZE(Hessian,1)
       TempVector = AtPoint
       TempVector( N ) = TempVector( N ) - Delta
       CALL DifferenceGradient( TheFunction, ThreePointDifferenceDerivative, TempVector, SqrtDelta, Gradient1 )
       TempVector( N ) = TempVector( N ) + Delta + Delta
       CALL DifferenceGradient( TheFunction, ThreePointDifferenceDerivative, TempVector, SqrtDelta, Gradient2 )

       Hessian( N, : ) = Factor * ( Gradient2( : ) - Gradient1( : ) )
    END DO

END SUBROUTINE ThreePointDifferenceHessian

!***************************************************************************************
!* Calculates the hessian at given point with 5 point difference method...

SUBROUTINE FivePointDifferenceHessian( TheFunction, AtPoint, Delta, Hessian )
IMPLICIT NONE
REAL, DIMENSION(:,:)                :: Hessian
REAL, DIMENSION(:)                  :: AtPoint
REAL, DIMENSION(1:SIZE(Hessian,1))  :: TempVector, Gradient1, Gradient2, Gradient3, Gradient4
INTEGER                             :: N
REAL                                :: Delta, SqrtDelta, Factor

INTERFACE
    REAL FUNCTION TheFunction( Point )
         REAL, DIMENSION(:) :: Point
    END FUNCTION TheFunction
END INTERFACE

    CALL ERROR(SIZE(AtPoint)/=SIZE(Hessian,1),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCEHESSIAN:" // &
                                              " vectors do not match.")
    CALL ERROR(SIZE(AtPoint)/=SIZE(Hessian,2),"DIFFERENCEDERIVATIVESMODULE.FIVEPOINTDIFFERENCEHESSIAN:" // &
                                              " vectors do not match.")

    SqrtDelta = SQRT( Delta )
    Factor = 1.0 / ( 12.0 * Delta )

    DO N = 1, SIZE(Hessian,1)
       TempVector = AtPoint
       TempVector( N ) = TempVector( N ) - Delta - Delta
       CALL DifferenceGradient( TheFunction, FivePointDifferenceDerivative, TempVector, SqrtDelta, Gradient1 )
       TempVector( N ) = TempVector( N ) + Delta
       CALL DifferenceGradient( TheFunction, FivePointDifferenceDerivative, TempVector, SqrtDelta, Gradient2 )
       TempVector( N ) = TempVector( N ) + Delta + Delta
       CALL DifferenceGradient( TheFunction, FivePointDifferenceDerivative, TempVector, SqrtDelta, Gradient3 )
       TempVector( N ) = TempVector( N ) + Delta
       CALL DifferenceGradient( TheFunction, FivePointDifferenceDerivative, TempVector, SqrtDelta, Gradient4 )

       Hessian( N, : ) = Factor * ( Gradient1(:) - 8.0 *  Gradient2(:) + 8.0 * Gradient3(:) - Gradient4(:) )
    END DO

END SUBROUTINE FivePointDifferenceHessian

!**************************** END OF FILE **********************************************

END MODULE DifferenceDerivatives
