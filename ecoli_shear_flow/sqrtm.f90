!===============================================================
!  Compute the principal square root X of an N-by-N real matrix A
!  specialized to the case N=3
!
!  1. N. J. Higham, Computing real square roots of a real  matrix, 
!          Linear Algebra and Appl., 88/89 (1987), pp. 405-430.
!  2. A. Bjorck and S. Hammarling, A Schur method for the square root of a matrix, 
!          Lin Algebra Appl, 52-53:127--140 (1983).
!===============================================================

subroutine sqrtm(AA,X)

  use nrtype
  IMPLICIT NONE

  integer(I4B), parameter :: N=3
  integer(I4B), parameter :: LDA = N
  integer(I4B), parameter :: LDVS = 2 * N
  integer(I4B), parameter :: LWORK = 3 * N

!!! input
  real(dp), dimension(N,N), intent(in) :: AA

!!! ouput
  real(dp), dimension(N,N), intent(out) :: X

  ! local variables

  complex(dpc),dimension(LDA,N) :: A
  complex(dpc),dimension(LDVS,N) :: VS
  complex(dpc),dimension(LWORK) :: WORK
  complex(dpc),dimension(N) :: W

  complex(dpc),dimension(N) :: RWORK
  integer(i4b) :: SDIM,INFO

  logical*8    BWORK(1:N)
  logical      select
  external     select

  integer(i4b) :: I,J,K
  complex(dpc),dimension(N,N) :: Q, QR, QT
  complex(dpc),dimension(N,N) :: CX, R
  complex(dpc) :: summ


  do i=1,N
     do j=1,N
        A(i,j)=dcmplx(AA(i,j))
     end do
  end do

!!!      write(6,*) 'original matrix A'
!!!      PRINT 1010, ((A(I,J), J = 1, N), I = 1, N)
!!!
!!!     Compute the eigenvalues 'W' and Schur vectors of A.
!!!
  call zgees('VECTORS AND EIGENVALUES','Not sorted',SELECT, &
       N, A, LDA, SDIM, W, VS, LDVS, WORK, LWORK,    &
       RWORK, BWORK, INFO)

!!!      write(6,*) 'info', info

!!!  after calling ZGEES, A implies an upper triangular matrix T
!!!  compute square root R of T, a column at a time.

  R=0.0_dp
  summ=0.0_dp

  do j=1,n
     R(j,j) = sqrt(A(j,j))
     do i=j-1,1,-1
        do k = i+1,j-1
           summ = summ + R(i,k)*R(k,j)
        enddo
        R(i,j) = (A(i,j) - summ)/(R(i,i) + R(j,j))
     enddo
  enddo

!!! copy Unitary matrix Q from VS

  do i=1,N
     do j=1,N
        Q(i,j)=VS(i,j)
     end do
  end do

!!! compute square root X of A : X=real(Q*R*Q**H)

  call TransposeMatrix(Q,QT,N)
  call MatrixMul(Q,R,QR,N)
  call MatrixMul(QR,QT,CX,N)

  do j=1,N
     do i=1,N
        X(i,j)=real(CX(i,j))
     end do
  end do

!!!      write(6,*) 'Square root X=real(Q*R*Q**H) of A'
!!!      PRINT 1080, ((X(I,J), J = 1, N), I = 1, N)


!!! 1010 FORMAT (6(3X, F17.14))
!!! 1080 FORMAT (3(3X, F17.14))


  RETURN
END subroutine sqrtm

!==================================================================

LOGICAL FUNCTION SELECT (ARG)

  use nrtype
  IMPLICIT NONE

  complex(dpc) :: ARG

!!!     The value computed is always .TRUE.  It is computed in
!!!     the peculiar way below to avoid compiler messages about
!!!     unused function arguments.

  PRINT 100
  SELECT = (ARG .EQ. ARG)
!!!        SELECT = .TRUE.

100 FORMAT (///1X, '**** ERROR:  ', &
       'SELECT FUNCTION CALLED BUT NOT AVAILABLE. ****')

  STOP
END FUNCTION SELECT

!==================================================================
!  square matrix (3x3) multiplication C=A*B
!==================================================================

Subroutine MatrixMul(A,B,C,N)

  use nrtype
  implicit none

!!!  input
  integer(i4b), intent(in) ::  N
  complex(dpc), dimension(N,N), intent(in) ::  A, B

!!! output
  complex(dpc), dimension(N,N), intent(out) :: C

!!! local variables
  integer(i4b) :: i, j, k

!!! initialize

  C=0.0_dp

!!! compute multiplication

  do j = 1, N
     do i = 1, N
        do k=1,N
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        end do
     end do
  end do

  return
end Subroutine MatrixMul

!==================================================================
!  transpose of a complex square matrix
!==================================================================

subroutine TransposeMatrix(A,AT,N)

  use nrtype
  implicit none

!!! watch out : this is complex conjugate transpose

!!! input
  integer(i4b), intent(in) :: N
  complex(dpc), dimension(N,N), intent(in) :: A
!!! output

  complex(dpc), dimension(N,N), intent(out) :: AT

!!! local variables
  integer(i4b) :: i, j

!!! compute transpose

  do i=1,N
     do j=1,N
        AT(i,j)=dconjg(A(j,i))
     end do
  end do

  return
end subroutine TransposeMatrix

!==================================================================
!==================================================================
