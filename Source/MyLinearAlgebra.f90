!***************************************************************************************
!*                              MODULE MyLinearAlgebra
!***************************************************************************************
!
!>  \brief            Wrapper for linear algebra operations
!>  \details          This module implements some linear algebra operations 
!>                    by referencing either LAPACK or Numerical Recipes for FORTRAN90.\n
!>                    When the operatios is trivial, it is directly implemented in
!>                    the module.
!
!***************************************************************************************
!
!>  \author           Matteo Bonfanti
!>  \version          1.0
!>  \date             June 2012
!
!***************************************************************************************
!
!>  \par Updates
!>  \arg 9 March 2013: diagonalization implemented with Numerical Recipes
!
!>  \todo          Implement diagonalization with LAPACK
!>  \todo          Implement matrix inversion with Numerical Recipes
!
!***************************************************************************************
!
!>  \remark       The module can be preprocessed with both the Numerical Recipes and 
!>                the LAPACK options, but in this case the LAPACK have priority 
!
!***************************************************************************************
!
!>  \ref         "NUMERICAL RECIPES IN FORTRAN 90:
!>               The Art of PARALLEL Scientific Computing"
!>               Chapter B3. Interpolation and Extrapolation 
!>               ISBN 0-521-57439-0 
!>               Copyright (C) 1986-1996 by Cambridge University Press 
!
!***************************************************************************************
MODULE MyLinearAlgebra
   USE ErrorTrap
   USE MyConsts
#if defined(WITH_NR)
   USE NRUtility
#endif
   
   IMPLICIT NONE

!********************************************************************************************************
   CONTAINS
!********************************************************************************************************



!*******************************************************************************
!          TheOneWithIdentityMatrix
!*******************************************************************************
!> Function giving back the identity matrix of size N in real format
!>
!> @param      N   size of the output matrix 
!> @returns    Identity matrix of size N
!*******************************************************************************
   FUNCTION TheOneWithIdentityMatrix( N ) RESULT( IdentityMatrix )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
      REAL, DIMENSION(N,N) :: IdentityMatrix
      INTEGER :: i

      ! initialize the matrix
      IdentityMatrix = 0.0

      ! set the diagonal element equal to 1
      DO i=1,N
         IdentityMatrix(i,i) = 1.0
      END DO

   END FUNCTION TheOneWithIdentityMatrix

   
!*******************************************************************************
!          TheOneWithInverseMatrix
!*******************************************************************************
!> Function giving back the inverse of a matrix.
!> The inverse is computed by solving a system of linear equations A * X = B
!> with A being the matrix to invert and B being a unit matrix.
!> This code is a wrapper to linear algebra libraries:
!> \arg Lapack routines xGESV for solving linear systems with a general matrix
!> \see http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve.html#ga7815e25af1fb6f8ee8fd5fd0fd1dc245
!>
!> @param      N          Integer size of the input square matrix.
!> @param      Matrix     NxN array with the matrix of which to compute the inverse.
!> @returns    Inverse    NxN array with the inverse of the matrix.
!*******************************************************************************
   FUNCTION TheOneWithInverseMatrix( Matrix, N ) RESULT( Inverse )
      IMPLICIT NONE
      INTEGER                            :: N
      REAL, DIMENSION(N,N), INTENT(IN)   :: Matrix
      REAL, DIMENSION(N,N)               :: Inverse
      
#if defined(WITH_LAPACK)
      INTEGER( SHORT_INTEGER_KIND )                 :: DimShort, Stat
      REAL, DIMENSION( N, N )                       :: Mat
      INTEGER( SHORT_INTEGER_KIND ), DIMENSION( N ) :: Pivot
      CHARACTER(100)                                :: ErrMsg
#endif

      ! Check and define the dimension of the matrices
      CALL ERROR( N /= SIZE(Matrix,2) , " TheOneWithInverseMatrix: input matrix is not square ")

      ! Initialize inverse as identity matrix
      Inverse = TheOneWithIdentityMatrix( N )
      
#if defined(WITH_LAPACK)
      ! Make a copy of the input matrix
      Mat = Matrix
      ! Define the dimension in a lapack compatible integer kind
      DimShort = N

      ! Check kind of real data
      IF ( KIND( Mat(1,1) ) == SINGLE_PRECISION_KIND ) THEN
            ! use lapack routine (single precision, general matrix, linear system solution )
            CALL SGESV( DimShort, DimShort, Mat, DimShort, Pivot, Inverse, DimShort, Stat )
            
      ELSE IF ( KIND( Mat(1,1) ) == DOUBLE_PRECISION_KIND ) THEN
            ! use lapack routine (double precision, general matrix, linear system solution )
            CALL DGESV( DimShort, DimShort, Mat, DimShort, Pivot, Inverse, DimShort, Stat )
      END IF

      ! chech if result is correctly computed
      IF ( Stat < 0 ) THEN
         WRITE(ErrMsg, *) " TheOneWithInverseMatrix: The argument ", -Stat, " had an illegal value."
         CALL AbortWithError( ErrMsg )
      END IF
      IF ( Stat > 0 ) THEN
         WRITE(ErrMsg, *) " TheOneWithInverseMatrix: The factor U is exactly singular, so the solution could not be computed."
         CALL AbortWithError( ErrMsg )
      END IF
#endif
#if !defined(WITH_LAPACK)
      CALL AbortWithError( " TheOneWithInverseMatrix: Matrix inversion implemented only with LAPACK ")
#endif

   END FUNCTION TheOneWithInverseMatrix

!*******************************************************************************
!          TheOneWithMatrixMultiplication
!*******************************************************************************
!> Function giving back the matrix matrix product of two real matrices. \n
!> Apart from a direct implementation of the product, this code is a wrapper
!> to BLAS subroutines xGEMM :
!> \see http://en.wikipedia.org/wiki/General_Matrix_Multiply
!> \see http://www.netlib.org/blas/dgemm.f
!>
!> @param      Matrix1     N1 x M  array with a first matrix.
!> @param      Matrix2     M x N2  array with a second matrix.
!> @returns    ProductM    N1 x N2 array with the product Matrix1 * Matrix2.
!*******************************************************************************
   FUNCTION TheOneWithMatrixMultiplication( Matrix1, Matrix2 ) RESULT( ProductM )
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(IN) :: Matrix1
      REAL, DIMENSION(:,:), INTENT(IN) :: Matrix2
      REAL, DIMENSION(size(Matrix1,1),size(Matrix2,2)) :: ProductM

#if defined(WITH_LAPACK)
      INTEGER( SHORT_INTEGER_KIND )                 :: N1, N2, M
#endif
#if !defined(WITH_LAPACK)
      INTEGER               :: i, j, k
#endif
      ! Check and define the dimension of the matrices
      CALL ERROR( size(Matrix1,2) /= SIZE(Matrix1,1) , &
                       " TheOneWithMatrixMultiplication: mismatch in matrix dimensions ")

#if defined(WITH_LAPACK)
      ! define dimensions
      N1 = size( Matrix1, 1 )
      N2 = size( Matrix2, 2 )
      M  = size( Matrix1, 2 )

      ! Check kind of real data
      IF ( KIND( Matrix1(1,1) ) == SINGLE_PRECISION_KIND ) THEN
            ! use blas routine (single precision, general matrix, matrix matrix product)
            CALL SGEMM ( 'N', 'N', N1, N2, M, 1.0, Matrix1, N1, Matrix2, M, 0.0, ProductM, N1 )
            
      ELSE IF ( KIND( Matrix1(1,1) ) == DOUBLE_PRECISION_KIND ) THEN
            ! use blas routine (double precision, general matrix, matrix matrix product )
            CALL DGEMM ( 'N', 'N', N1, N2, M, 1.0, Matrix1, N1, Matrix2, M, 0.0, ProductM, N1 )
      END IF
#endif
#if !defined(WITH_LAPACK)
      ProductM = 0.0
      DO i = 1, size( Matrix1, 1 )
         DO j  = 1, size( Matrix2, 2 )
            DO k = 1, size( Matrix1, 2 )
               ProductM(i,j) = ProductM(i,j) + Matrix1(i,k)*Matrix2(k,j)
            END DO
         END DO
      END DO
#endif

   END FUNCTION TheOneWithMatrixMultiplication
   
!*******************************************************************************
!          TheOneWithMatrixVectorProduct
!*******************************************************************************
!> Function giving back the real-matrix real-vector product. \n
!> Apart from a direct implementation of the product, this code is a wrapper
!> to BLAS subroutines xGEMM :
!> \see http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
!> \see http://www.netlib.org/blas/dgemv.f
!>
!> @param      Matrix      N x M  real array with a matrix.
!> @param      Vector      M      real array with a vector.
!> @returns    ProductV    N      array with the product Matrix * Vector.
!*******************************************************************************
   FUNCTION TheOneWithMatrixVectorProduct( Matrix, Vector ) RESULT( ProductV )
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(IN) :: Matrix
      REAL, DIMENSION(:), INTENT(IN)   :: Vector
      REAL, DIMENSION(size(Matrix,1)) :: ProductV

#if defined(WITH_LAPACK)
      INTEGER( SHORT_INTEGER_KIND )   :: N, M
#endif
#if !defined(WITH_LAPACK)
      INTEGER                         :: i, j
#endif
      ! Check and define the dimension of the matrices
      CALL ERROR( size(Matrix,2) /= SIZE(Vector) , &
                       " TheOneWithMatrixVectorProduct: mismatch in matrix dimensions ")

#if defined(WITH_LAPACK)
      ! define dimensions
      N = size( Matrix, 1 )
      M = size( Matrix, 2 )

      ! Check kind of real data
      IF ( KIND( Matrix(1,1) ) == SINGLE_PRECISION_KIND ) THEN
            ! use blas routine (single precision, general matrix, matrix vector product)
!             CALL SGEMM ( 'N', 'N', N, 1, M, 1.0, Matrix, N, Vector, M, 0.0, ProductV, N )
             CALL SGEMV ( 'N', N, M, 1.0, Matrix, N, Vector, 1, 0.0, ProductV, 1 )
            
      ELSE IF ( KIND( Matrix(1,1) ) == DOUBLE_PRECISION_KIND ) THEN
            ! use blas routine (double precision, general matrix, matrix vector product )
!             CALL DGEMM ( 'N', 'N', N, 1, M, 1.0, Matrix, N, Vector, M, 0.0, ProductV, N )
             CALL DGEMV ( 'N', N, M, 1.0, Matrix, N, Vector, 1, 0.0, ProductV, 1 )
      END IF
#endif
#if !defined(WITH_LAPACK)
      ProductV = 0.0
      DO j  = 1, size(Matrix,2)
         DO i = 1, size(Matrix,1)
               ProductV(i) = ProductV(i) + Matrix(i,j)*Vector(j)
         END DO
      END DO
#endif

   END FUNCTION TheOneWithMatrixVectorProduct

!*******************************************************************************
!          TheOneWithVectorDotVector
!*******************************************************************************
!> Function giving back the real-matrix real-vector product. \n
!> Apart from a direct implementation of the product, this code is a wrapper
!> to BLAS subroutines xGEMM :
!> \see http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms
!> \see http://www.netlib.org/blas/dgemv.f
!>
!> @param      Vector1      N  real array with a vector.
!> @param      Vector2      N  real array with a vector.
!> @returns    VDotV        the scalar product Vector1 * Vector2.
!*******************************************************************************
   FUNCTION TheOneWithVectorDotVector( Vector1, Vector2 ) RESULT( VDotV )
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN)   :: Vector1, Vector2
      REAL                             :: VDotV
      REAL(kind=SINGLE_PRECISION_KIND) :: SDOT
      REAL(kind=DOUBLE_PRECISION_KIND) :: DDOT

#if defined(WITH_LAPACK)
      INTEGER( SHORT_INTEGER_KIND )   :: N
#endif
#if !defined(WITH_LAPACK)
      INTEGER                         :: i
#endif
      ! Check and define the dimension of the matrices
      CALL ERROR( size(Vector1) /= SIZE(Vector2) , &
                       " TheOneWithVectorDotVector: mismatch in matrix dimensions ")

#if defined(WITH_LAPACK)
      ! define dimensions
      N = size( Vector1 )

      ! Check kind of real data
      IF ( KIND( Vector1(1) ) == SINGLE_PRECISION_KIND ) THEN
            ! use blas routine (single precision, vector scalar product)
             VDotV = SDOT( N, Vector1,1, Vector2,1 ) 
            
      ELSE IF ( KIND( Vector1(1) ) == DOUBLE_PRECISION_KIND ) THEN
            ! use blas routine (double precision, vector scalar product )
             VDotV = DDOT( N, Vector1,1, Vector2,1 ) 
      END IF
#endif
#if !defined(WITH_LAPACK)
      VDotV = 0.0
      DO i = 1, size(Vector1)
            VDotV = VDotV + Vector1(i)*Vector2(i)
      END DO
#endif

   END FUNCTION TheOneWithVectorDotVector


!*******************************************************************************
!          TheOneWithDiagonalization
!*******************************************************************************
!> Function giving the eigenvectors and eigenvalues of a NxN symmetric real 
!> matrix. \n
!> This code is a wrapper to the numerical recipe subroutines tred2 and tqli
!> (householder reduction + QL algorithm )
!> \see http://apps.nrbook.com/fortran/index.html
!>
!> @param    Matrix         N x N  real symmetric matrix.
!> @returns  EigenVectors   N x N  real matrix to store eigenvectors of Matrix.
!> @returns  EigenValues    N      real array to store the eigenvalues of Matrix.
!*******************************************************************************
   SUBROUTINE TheOneWithDiagonalization(Matrix,EigenVectors,EigenValues)
      IMPLICIT NONE
      REAL, DIMENSION(:,:), INTENT(IN)  :: Matrix
      REAL, DIMENSION(:,:), INTENT(OUT) :: EigenVectors
      REAL, DIMENSION(:), INTENT(OUT)   :: EigenValues
      REAL, DIMENSION(:), ALLOCATABLE   :: OffDiagonal
      INTEGER  :: N
      
      ! Check and define the dimension of the matrices
      N = size(Matrix,1)
      CALL ERROR( size(Matrix,2) /= N , " TheOneWithDiagonalization: input matrix is not square ")
      CALL ERROR( size(EigenVectors,1) /= N , " TheOneWithDiagonalization: eigenvector matrix mismatch (1) ")
      CALL ERROR( size(EigenVectors,2) /= N , " TheOneWithDiagonalization: eigenvector matrix mismatch (2) ")
      CALL ERROR( size(EigenValues) /= N , " TheOneWithDiagonalization: eigenvalues vector mismatch ")
      
#if defined(WITH_NR)
      ALLOCATE( OffDiagonal(N) )

      EigenVectors=Matrix
      CALL tred2(EigenVectors,EigenValues,OffDiagonal) ! Householder reduction
      CALL tqli(EigenValues,OffDiagonal,EigenVectors)        ! QL algorithm       

      DEALLOCATE( OffDiagonal )
#endif
#if !defined(WITH_NR)
      CALL AbortWithError( " TheOneWithDiagonalization: Matrix diagonalization implemented only with NR ")
#endif

   END SUBROUTINE TheOneWithDiagonalization
   
#if defined(WITH_NR)
!* * * * * * * * * * * * * * * * * NR SUBROUTINE * * * * * * * * * * * * * * * * * * * *
!* The following subroutines are taken from NR for FORTRAN 90/95

SUBROUTINE tred2(a,d,e,novectors)
IMPLICIT NONE
REAL , DIMENSION(:,:), INTENT(INOUT) :: a
REAL, DIMENSION(:), INTENT(OUT) :: d,e
LOGICAL, OPTIONAL, INTENT(IN) :: novectors
! Householder reduction of a real, symmetric, NxN matrix a. On output, a is replaced by
! an orthogonal matrix Q effecting the trasformation. d returns the diagonal elements of the
! tridiagonal matrix, and e the off-diagonal elements, with e(1)=0. If the optional argument
! novectors is present and .true. only eigenvalues are to be found subsequently, in which
! case a contains no useful information on output
INTEGER(I4B) :: i,j,l,n
REAL :: f,g,h,hh,scale
REAL, DIMENSION(size(a,1)) :: gg
LOGICAL :: yesvec
n = assert_eq(size(a,1), size(a,2), size(d), size(e), 'tred2')
if (present(novectors)) then
   yesvec=.not. novectors
else
   yesvec=.true.
end if
do i=n,2,-1
   l=i-1
   h=0.0
   if (l>1) then
      scale=sum(abs(a(i,1:l)))
      if (scale == 0 ) then      ! skip trasformation
         e(i)=a(i,l)
      else
         a(i,1:l)=a(i,1:l)/scale   ! use scale a's for trasformation
         h=sum(a(i,1:l)**2)   ! form sigma in h
         f=a(i,l)
         g=-sign(sqrt(h),f)
         e(i)=scale*g
         h=h-f*g         ! now h is equation (11.2.4)
         a(i,l)=f-g         ! store u in the ith row of a
         if (yesvec) a(1:l,i)=a(i,1:l)/h   !store u/H in ith column of a
         do j=1,l
            e(j)=(dot_product(a(j,1:j),a(i,1:j)) &   ! unused elements of e
            +dot_product(a(j+1:l,j),a(i,j+1:l)))/h
         end do
         f=dot_product(e(1:l),a(i,1:l))
         hh=f/(h+h)      ! form K eq. (11.2.11)
         e(1:l)=e(1:l)-hh*a(i,1:l)   ! form q and store in e overwriting p
         do j=1,l         ! reduce a, eq (11.2.13)
            a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
         end do
      end if
   else
      e(i)=a(i,l)
   end if
   d(i)=h
end do
if (yesvec) d(1)=0.0
e(1)=0.0
do i=1,n   ! begin accumulation of transformation matrices
   if (yesvec) then
      l=i-1
      if (d(i)/=0.0) then    ! this block skipped when i=1. Use u and u/H stored to form P.Q
         gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
         a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
      end if
      d(i)=a(i,i)
      a(i,i)=1.0   ! Reset row and column of a to identity matrix for next iteration
      a(i,1:l)=0.0
      a(1:l,i)=0.0
   else
      d(i)=a(i,i)
   end if
end do
END SUBROUTINE tred2

SUBROUTINE tqli(d,e,z)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(INOUT) :: d,e
REAL, DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
! QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
! simmetric, tridiagonal matrix, or of a real simmetric matrix previously reduced by tred2
! d is a vector of lenght N. On input, its elements are the diagonal elements of the tridiagonal
! matrix. On output, it returns the eigenvalue. The vector e inputs the subdiagonal elements of the
! tridiagonal matrix, with e(1) arbitrary. On output e is destroyed. When finding only the eigenvalues, the
! optional argument z is omitted. If the eigenvectors of a tridiagonal matrix are desired, the NxN matrix z
! is input as the identity matrix. If the eigenvectors of a matrix that has been reduced by tred2 are
! required, then z is input as the matrix output by tred2. In either case, the kth column of z returns the
! normalized eigenvector corrisponding to d(k).
INTEGER(I4B) :: i,iter,l,m,n,ndum
REAL :: b,c,dd,f,g,p,r,s
REAL, DIMENSION(size(e)) :: ff
n=assert_eq(size(d),size(e),'tqli: n')
if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
e(:)=eoshift(e(:),1)
do l=1,n
    iter=0
   iterate: do
      do m=l,n-1
         dd=abs(d(m))+abs(d(m+1))
         if (abs(e(m))+dd == dd ) exit
      end do
      if (m==l) exit iterate
      if (iter == 30) call nrerror('too many iteration in tqli')
      iter=iter+1
      g=(d(l+1)-d(l))/(2.0*e(l))
      r=pythag(g,1.0)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1.0
      c=1.0
      p=0.0
      do i=m-1,l,-1
         f=s*e(i)
         b=c*e(i)
         r=pythag(f,g)
         e(i+1)=r
         if (r==0.0) then
            d(i+1)=d(i+1)-p
            e(m)=0.0
            cycle iterate
         end if
         s=f/r
         c=g/r
         g=d(i+1)-p
         r=(d(i)-g)*s+2.0*c*b
         p=s*r
         d(i+1)=g+p
         g=c*r-b
         if (present(z)) then
            ff(1:n)=z(1:n,i+1)
            z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
            z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
         end if
      end do
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.0
   end do iterate
end do
END SUBROUTINE tqli

! FUNCTION pythag(a,b)
! IMPLICIT NONE
! REAL, INTENT(IN) :: a,b
! REAL :: pythag
! REAL :: absa,absb
!    absa=abs(a)
!    absb=abs(b)
!    if (absa>absb) then
!       pythag=absa*sqrt(1.0+(absb/absa)**2)
!    else if (absb == 0.0) then
!       pythag=0.0
!    else
!       pythag=absb*sqrt(1.0+(absa/absb)**2)
!    end if
! END FUNCTION pythag

#endif
   
END MODULE MyLinearAlgebra
