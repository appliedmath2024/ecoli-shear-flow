MODULE interface_rod

  INTERFACE
     SUBROUTINE initial_higdon(L_X,D1,D2,D3,L0,r0,w,kappa,kappa2,tau,tau2)
       USE nrtype; USE commonval
       IMPLICIT NONE

       REAL(DP), DIMENSION(:,:), INTENT(OUT) :: L_X
       REAL(DP), DIMENSION(:,:), INTENT(OUT) :: D1, D2, D3
       REAL(DP), DIMENSION(:), INTENT(OUT) :: L0, kappa, kappa2, tau, tau2
       REAL(DP), INTENT(IN) :: r0, w

!!! local variables
       INTEGER(I4B) :: i, j, iPt
       REAL(DP) :: k, deno, denoD3, s0
       REAL(DP), dimension(nPt) :: x,y,z,dx,dy,dz,r,dr,ddr,dddr, ddx,ddy,ddz,dddx,dddy,dddz,s,TT
       REAL(DP), dimension(3,nPT) :: dTds
       REAL(DP),dimension(3) :: crs
       REAL(DP), DIMENSION(3,nPt) :: D1ds, D2ds, D3ds, DH1, DH2, DH3
       real(dp), dimension(3,3) :: TA

       REAL(DP), DIMENSION(10000*Expmt+1) :: Tspan, Sspan

     END SUBROUTINE initial_higdon
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     subroutine calcForceMoment_mot(iStep,klok,nSkip,coord_pt,D1,D2,D3, &
          moment_pt,force_pt,kappa,tau,kkh,Tan)
       use nrtype; use commonval
       IMPLICIT NONE

       integer(i4b), intent(in) :: iStep, klok, nSkip
       real(dp), dimension(:,:), intent(in) :: coord_pt, D1, D2, D3
       real(dp), dimension(:,:), intent(out) :: moment_pt, force_pt,kkh
!!$real(dp), dimension(:,:), intent(inout) :: allEne ! all energy type
       real(dp), dimension(:), intent(in) :: kappa, tau, Tan

!!$ local
       integer(i4b) :: iPt
       real(dp), dimension(1:3, 1:nPt-1) :: DH1, DH2, DH3

!!! FH and NH should be written in cartesian coordinates
!!!   in subroutine MomentnForce

       real(dp), dimension(1:3, 1:nPt+1) :: FH, NH ! two extra points
       real(dp), dimension(1:3, 1:nPt-1) :: FHT, NHT
       real(dp), dimension(1:3, 1:nPt+1) :: TC ! =dX/ds ! tangent vectors

     END SUBROUTINE calcForceMoment_mot
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
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

     END subroutine sqrtm
  END INTERFACE
!!!!!!!!!!!!!!!
  Interface
     SUBROUTINE initial_body(nref)
       USE nrtype; USE commonval
       IMPLICIT NONE
       INTEGER(I4B), INTENT(INOUT) :: nref
       REAL(DP) :: theta, dheight, t, rr, turn_rate, bdhand
       INTEGER(I4B) :: n, i, j, ii, jj, kk, nvmax, ntmax, nv, nt, nr, nh, nh2
       REAL(DP), DIMENSION(3) :: x1, x2, x3, x4, Rx, Nx, Bx
       INTEGER(I4B), DIMENSION(3) :: next
       INTEGER(I4B), DIMENSION(:,:), Allocatable :: v
       REAL(DP), DIMENSION(:,:), Allocatable :: x

     end SUBROUTINE initial_body
  end interface
!!!!!!!!!!!!!!!
  INTERFACE
     SUBROUTINE frame_CSP(LX,Forc,Torq,YCM,Fram)
       USE nrtype; USE commonval
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:), INTENT(OUT)  :: Fram
       REAL(DP), DIMENSION(:,:), INTENT(IN)  :: LX
       REAL(DP), DIMENSION(:), INTENT(IN)  :: Forc, Torq
       REAL(DP), DIMENSION(:), INTENT(OUT)  :: YCM
       REAL(DP), DIMENSION(3,3)  ::Tran, PSA, IPSA, MatM, SS, Fram0, Fram1, SNST, SOST, Test, NTorq, TranSS
       REAL(DP), DIMENSION(3) :: Diag
       REAL(DP) :: detA, theta
       INTEGER(I4B) :: i, j, LWORK=8, info
       REAL(DP), DIMENSION(8) :: WORK

     END SUBROUTINE frame_CSP
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     subroutine StSolve(data_Pt,moment_Pt,force_Pt,data_vel,LX,LF,LU)
       USE omp_lib
       USE nrtype; USE commonval
       implicit none


       REAL(DP),DIMENSION(:,:,:), INTENT(IN)  :: data_Pt,moment_Pt,force_Pt
       REAL(DP),DIMENSION(:,:),   INTENT(IN)  :: LX, LF
       REAL(DP),DIMENSION(:,:,:), INTENT(OUT) :: data_vel
       REAL(DP),DIMENSION(:,:),   INTENT(OUT) :: LU
       REAL(DP),DIMENSION(nPt) :: r,H1,H2,H3,deno,dx,dy,dz,fdotX,rXdel
       REAL(DP),DIMENSION(nPt) :: uSt,vSt,wSt,uRot,vRot,wRot
       REAL(DP),DIMENSION(nPt) :: xc, yc, zc
       REAL(DP),DIMENSION(Kvmax) :: lr,lH1,lH2,lH3,ldeno,ldx,ldy,ldz,lfdotX,lxc,lyc,lzc,lrXdel
       REAL(DP),DIMENSION(Kvmax) :: LU1, LU2, LU3
       integer(I4B) :: jPt, i, j

     END SUBROUTINE StSolve
  END INTERFACE
  !!!!!!!!!!!!!!!
  INTERFACE
     subroutine StSolve_wall(data_Pt,moment_Pt,force_Pt,data_vel,LX,LF,LU)
       USE omp_lib
       USE nrtype; USE commonval
       implicit none


real(dp),dimension(:,:,:),intent(in) :: data_Pt, moment_Pt,force_Pt
real(dp),dimension(:,:,:),intent(out) :: data_vel
REAL(DP),DIMENSION(:,:),  INTENT(IN)  :: LX, LF
REAL(DP),DIMENSION(:,:),  INTENT(OUT) :: LU
real(dp),dimension(3,nPt,nrod) :: image_pt
real(dp),dimension(3,Kvmax) :: LX_im
real(dp),dimension(nPt) :: r,H1,H2,H3,deno,dx,dy,dz,fdotX,rXdel,rD1,rD2,deno2,gdotX,pdotX,LxX1,LxX2,LxX3,H4
real(dp),dimension(nPt) :: uSt, vSt, wSt, uRot, vRot, wRot
real(dp),dimension(Kvmax) :: lr, lH1, lH2, lH3, ldeno, ldx, ldy, ldz, lfdotX, lrXdel
real(dp),dimension(Kvmax) :: lD1, lD2, ldeno2, lgdotX, lpdotX, lLxX1, lLxX2, lLxX3, lH4
real(dp),dimension(Kvmax) :: LU1, LU2, LU3
integer(I4B) :: jPt, k, i
real(dp) :: hh,lhh

     END SUBROUTINE StSolve_Wall
  END INTERFACE

!!!!!!!!!!!!!!!
  INTERFACE
     subroutine StSolveMk(data_Pt,moment_Pt,force_Pt,marker_pt,data_vel,LX,LF)
       USE omp_lib
       USE nrtype; USE commonval
       implicit none

       real(dp),dimension(:,:,:),intent(in) :: data_Pt,moment_Pt,force_Pt
       real(dp),dimension(:,:),intent(in) :: marker_pt
       REAL(DP),DIMENSION(:,:),   INTENT(IN)  :: LX, LF
       real(dp),dimension(:,:),intent(out) :: data_vel
       real(dp),dimension(3,nPt,nrod) :: image_pt
       real(dp),dimension(3,Kvmax) :: LX_im
       real(dp),dimension(:),allocatable :: r,H1,H2,H3,deno,dx,dy,dz,fdotX,rXdel,rD1,rD2,deno2,gdotX,pdotX,LxX1,LxX2,LxX3,H4
       real(dp),dimension(:),allocatable :: uSt,vSt,wSt,uRot,vRot,wRot
       real(dp),dimension(:),allocatable :: xc,yc,zc
       integer(I4B) :: jPt,N,i,j
       real(dp) :: hh

     END SUBROUTINE StSolveMk
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     subroutine AngSolve(data_Pt,moment_Pt,force_Pt,LX,LF,data_ang)
       USE omp_lib
       USE nrtype; USE commonval
       implicit none

       real(dp),dimension(:,:,:),intent(in) :: data_Pt,moment_Pt,force_Pt
       real(dp),dimension(:,:),intent(in) :: LX,LF
       real(dp),dimension(:,:,:),intent(out) :: data_ang
       real(dp),dimension(nPt) :: r,dx,dy,dz,ndotX,Q1,Q2,Q3, rXdel, Qinv
       real(dp),dimension(nPt) :: uARot,vARot,wARot,uADi,vADi,wADi
       real(dp),dimension(nPt) :: xc, yc, zc
       integer(I4B) :: jPt, i, j

     END SUBROUTINE AngSolve
  END INTERFACE
  !!!!!!!!!!!!!!!
  INTERFACE
     subroutine AngSolve_wall(data_Pt,moment_Pt,force_Pt,LX,LF,data_ang)
       USE omp_lib
       USE nrtype; USE commonval
       implicit none

real(dp),dimension(:,:,:),intent(in) :: data_Pt,moment_Pt,force_Pt
real(dp),dimension(:,:,:),intent(out) :: data_ang
real(dp),dimension(:,:),intent(in) :: LX,LF
real(dp),dimension(3,nPt,nrod) :: image_pt
real(dp),dimension(3,Kvmax) :: LX_im
real(dp),dimension(nPt) :: r, dx, dy, dz, ndotX, Q1, Q2, Q3, rXdel, Qinv
real(dp),dimension(nPt) :: rD1, rD2, deno, deno2, gdotX, H4, H5, qdotX, pdotX
real(dp),dimension(nPt) :: uARot, vARot, wARot, uADi, vADi, wADi
real(dp),dimension(nPt) :: xc, yc, zc
integer(I4B) :: iPt, jPt, k, i
real(dp) :: hh,lhh


     END SUBROUTINE AngSolve_Wall
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     subroutine updateTriad_mot (iStep,ang_pt,Tan,D1, D2, D3)
       USE omp_lib
       use nrtype; use commonval
       implicit none

       integer(i4b), intent(in) :: iStep
       real(dp), dimension(:,:), intent(in) :: ang_pt
       real(dp), dimension(:), intent(in) :: Tan
       real(dp), dimension(:,:), intent(inout) :: D1, D2, D3
       real(dp), dimension(3,nPt) :: DN1, DN2, DN3 ! triads at time step n+1
       real(dp), dimension(3,nPt) :: DT1, DT2, DT3
       real(dp), dimension(3,nPt) :: wxD1, wxD2, wxD3 ! cross product of ang_pt and D1
       integer(i4b) :: I, J, K, iPt, iDim
       real(dp) :: temp1, temp2, temp3
       real(dp), dimension(1:nPt) :: normw, invnormw
       real(dp) :: angw, smooth, cosval, sinval
       real(dp), dimension(3) :: normal

     END SUBROUTINE updateTriad_mot
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     SUBROUTINE rotation(Vec,theta,Tan)
       USE nrtype; USE commonval
       IMPLICIT NONE
       REAL(DP), DIMENSION(:), INTENT(INOUT)  :: Vec
       REAL(DP), INTENT(IN)  :: theta
       REAL(DP), DIMENSION(:), INTENT(In)  :: Tan

       REAL(DP), DIMENSION(3,3)  ::Tran, Iden, PSA
       REAL(DP), DIMENSION(3) :: Omg
       REAL(DP) :: detA
       INTEGER(I4B) :: i, j

     END SUBROUTINE rotation
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     SUBROUTINE selfcontactForce(coord_pt,force_pt)
       USE omp_lib
       USE nrtype; USE commonval
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(in) :: coord_pt
       REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: force_pt

!!! local variables
       INTEGER(I4B) :: i, j, k, jPt
       REAL(DP), DIMENSION(nPt) :: xc, yc, zc, dx, dy, dz, r, ten, fx, fy, fz, dist

     END SUBROUTINE selfcontactForce
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     SUBROUTINE selfcontactForce_Fla_Body(coord_pt,force_pt, LX, LF)
       USE omp_lib
       USE nrtype; USE commonval
       IMPLICIT NONE
       REAL(DP), DIMENSION(:,:,:), INTENT(in) :: coord_pt
       REAL(DP), DIMENSION(:,:), INTENT(in) :: LX
       REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: force_pt
       REAL(DP), DIMENSION(:,:), INTENT(inout) :: LF

!!! local variables
       INTEGER(I4B) :: i, j, k, jPt
       REAL(DP), DIMENSION(nPt) :: xc, yc, zc, dx, dy, dz, r, ten, fx, fy, fz, dist

     END SUBROUTINE selfcontactForce_Fla_Body
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
     SUBROUTINE Force(LX,LF)
       USE omp_lib
       USE nrtype; USE commonval
       IMPLICIT NONE

       REAL(DP), DIMENSION(:,:), INTENT(IN)  :: LX
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: LF

       INTEGER(I4B)  :: i, j, k, ii, jj, kk, i_x, i_y, i_z
       REAL(DP), DIMENSION(3) :: x1, x2, x3, x4, ff
       REAL(DP) :: norm, Ten
       REAL(DP), DIMENSION(3,3,Ktmax) :: PSLF

     END SUBROUTINE Force
  END INTERFACE
!!!!!!!!!!!!!!!
  INTERFACE
      SUBROUTINE selfcontactForce_wall(coord_pt,force_pt, LX, LF)
        USE nrtype; USE commonval
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:,:), INTENT(in) :: coord_pt
        REAL(DP), DIMENSION(:,:), INTENT(in) :: LX
        REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: force_pt
        REAL(DP), DIMENSION(:,:), INTENT(inout) :: LF
      
      !!! local variables
        INTEGER(I4B) :: i, j, k, jPt
        REAL(DP) :: dy, r, ten, fy, dist

     END SUBROUTINE selfcontactForce_wall
  END INTERFACE

!!!!!!!!!!!!!!!
END MODULE interface_rod

