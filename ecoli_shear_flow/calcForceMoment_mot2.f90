!  Calculate forces at fiber points   sook 11/15/08
!========================================================================

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

!!! calculate the tangent vectors


!$OMP PARALLEL DO PRIVATE(iPt)
do iPt=2,nPt
  TC(:,iPt)=(coord_pt(:,iPt)-coord_pt(:,iPt-1))/ds
enddo
!$OMP END PARALLEL DO

TC(:,1)=TC(:,2)         ! X0=2X1-X2
TC(:,nPt+1)=TC(:,nPt)   ! X_{k+1}=2*X_k-X_{k-1}

!!! define the orthonormal triads at half points

  call TriadAtHalfPt(D1,D2,D3,DH1,DH2,DH3)

!!! compute the internal forces at half point

  call InternalForceAtHalfPt (coord_pt,DH1,DH2,DH3,TC,FH,FHT)

!!! compute the internal moments at half point

  call InternalMomentAtHalfPt (DH1,DH2,DH3,D1,D2,D3,NH,NHT,kappa,tau,kkh)

!!!write(*,*) 'test1'
!!!write(*,10) TC
!!!write(*,*) 'test2'
!!!write(*,10) FH
!!!write(*,*) 'test3'
!!!write(*,10) NH
!!!10 format (f20.16,1x,f20.16,1x,f20.16)

!!! compute the applied moment and force to the filament

  call MomentnForce(coord_pt,NH,FH,TC,moment_pt,force_pt)

!!! write the energy
!!$if (mod(iStep,nSkip)==0) then
!!$       call WriteEnergy(iStep, klok, FHT, NHT,allEne)
!!$endif

!!!============================================================================
contains
!!!============================================================================

  subroutine TriadAtHalfPt(D1,D2,D3,DH1,DH2,DH3)
    use nrtype; use commonval
    implicit none

!!! input
!!! D1,D2,D3 : orthonormal triads at boundary points

!!! output
!!! DH1, DH2, DH3 : orthonormal triads at half points

!!! Variables declaration
    real(dp), dimension(:,:), intent(in) :: D1, D2, D3
    real(dp), dimension(:,:), intent(out) :: DH1, DH2, DH3

!!! local variables
    integer(i4b) :: iPt, i, j, iDim
    real(dp), dimension(3,3) :: TA

!!! compute orthonormal triad at half point
!$OMP PARALLEL DO PRIVATE(iPt, i, j, TA)
 do iPt=1,nPt-1
    do j=1,3; do i=1,3
      TA(i,j) = D1(i,iPt+1)*D1(j,iPt)+D2(i,iPt+1)*D2(j,iPt)+D3(i,iPt+1)*D3(j,iPt)
    end do; end do

    call sqrtm(TA,TA)

    DH1(:,iPt) = TA(:,1)*D1(1,iPt)+TA(:,2)*D1(2,iPt)+TA(:,3)*D1(3,iPt)

    DH2(:,iPt) = TA(:,1)*D2(1,iPt)+TA(:,2)*D2(2,iPt)+TA(:,3)*D2(3,iPt)

    DH3(:,iPt) = TA(:,1)*D3(1,iPt)+TA(:,2)*D3(2,iPt)+TA(:,3)*D3(3,iPt)
 end do
!$OMP END PARALLEL DO

  end subroutine TriadAtHalfPt

!!!============================================================================

  subroutine InternalForceAtHalfPt (coord_pt,DH1,DH2,DH3,TC,FH,FHT)
    use nrtype; use commonval
    implicit none

!!! Input arguments

!!! nCon - number of connections in the pc file
!!! nPt - # of boundary points of the filament
!!! DH1, DH2, DH3 - orthonormal triad at half boundary point

!!! Output arguments

!!! FH - internal force vector at each half point
!!! FH will be written out into cartesian coordinates

!!! Variables declaration
    real(dp), dimension(:,:), intent(in) :: coord_pt
    real(dp), dimension(:,:), intent(in) :: DH1, DH2, DH3
    real(dp), dimension(:,:), intent(in) :: TC
    real(dp), dimension(:,:), intent(out) :: FH
    real(dp), dimension(:,:), intent(out) :: FHT !temp for FH in basis(DH1,DH2,DH3)

!!! local variables
    integer(i4b) :: iPt

!!! compute the internal forces at each half point
!$OMP PARALLEL
!$OMP DO PRIVATE(iPt)
  do iPt = 1, nPt-1
     FHT(1,iPt)= shear * sum( DH1(:,iPt)*TC(:,iPt+1) )

     FHT(2,iPt)= shear * sum( DH2(:,iPt)*TC(:,iPt+1) )

     FHT(3,iPt)= stretch * (sum(DH3(:,iPt)*TC(:,iPt+1))-1.0_dp)
  end do
!$OMP END DO

! FHT above are components in the basis (DH1,DH2,DH3)
! they need to be changed into cartesian coordinates to be used
! in calculation of force and torque

!$OMP DO PRIVATE(iPt)
  do iPt = 1,nPt-1
      FH(1,iPt+1) = FHT(1,iPt)*DH1(1,iPt) &
                   +FHT(2,iPt)*DH2(1,iPt) &
                   +FHT(3,iPt)*DH3(1,iPt)

      FH(2,iPt+1) = FHT(1,iPt)*DH1(2,iPt) &
                   +FHT(2,iPt)*DH2(2,iPt) &
                   +FHT(3,iPt)*DH3(2,iPt)

      FH(3,iPt+1) = FHT(1,iPt)*DH1(3,iPt) &
                   +FHT(2,iPt)*DH2(3,iPt) &
                   +FHT(3,iPt)*DH3(3,iPt)
  end do
!$OMP END DO
!$OMP END PARALLEL

!!! by reflection

    FH(1,1:nmot) = 0.0_dp*FH(1,nmot+1) ! tethered motor
    FH(2,1:nmot) = 0.0_dp*FH(2,nmot+1) ! tethered motor
    FH(3,1:nmot) = 0.0_dp*FH(3,nmot+1) ! tethered motor
    FH(1,nPt+1) = 0.0_dp*FH(1,nPt)
    FH(2,nPt+1) = 0.0_dp*FH(2,nPt)
    FH(3,nPt+1) = 0.0_dp*FH(3,nPt)

  end subroutine InternalForceAtHalfPt

!!!============================================================================

  subroutine InternalMomentAtHalfPt (DH1,DH2,DH3,D1,D2,D3,NH,NHT,kappa,tau,kkh)
    use nrtype; use commonval
    implicit none

!!! Input arguments

!!! nPt - # of boundary points of the filament
!!! DH1, DH2, DH3 - orthonormal triad at each half point
!!! D1, D2, D3 - orthonormal triad at each boundary point

!!! Output arguments

!!! NH - internal moment vector at each half point
!!! NH will be written out into cartesian coordinates

!!! Variables declaration
    real(dp), dimension(:,:), intent(in) :: D1, D2, D3
    real(dp), dimension(:,:), intent(in) :: DH1, DH2, DH3
    real(dp), dimension(:,:), intent(out) :: NH
    real(dp), dimension(:,:), intent(out) :: NHT !temp for NH in basis(DH1,DH2,DH3)
    real(dp), dimension(:), intent(in) :: kappa, tau
    real(dp), dimension(:,:), intent(out) :: kkh

!!! local variables
    integer(i4b) :: iPt, iDim
    real(dp), dimension(1:3,1:nPt-1) :: KH !curvature at half point
    real(dp), dimension(1:3,1:nPt-1) :: DD1, DD2, DD3 ! = dD/ds

!!! compute a curvature at half point

!$OMP PARALLEL
!$OMP DO PRIVATE(iPt)
   do iPt = 1, nPt-1
      DD1(:,iPt)=(D1(:,iPt+1)-D1(:,iPt))/ds
      DD2(:,iPt)=(D2(:,iPt+1)-D2(:,iPt))/ds
      DD3(:,iPt)=(D3(:,iPt+1)-D3(:,iPt))/ds

      KH(1,iPt)=   sum(DD2(:,iPt)*DH3(:,iPt))
      KH(2,iPt)=   sum(DD3(:,iPt)*DH1(:,iPt))
      KH(3,iPt)=   sum(DD1(:,iPt)*DH2(:,iPt))

      NHT(1,iPt) = arry_bend1(iPt) * (KH(1,iPt)-curva2)
      NHT(2,iPt) = arry_bend2(iPt) * (KH(2,iPt)-Kappa(iPt))
      NHT(3,iPt) = arry_twist(iPt) * (KH(3,iPt)-tau(iPt)) ! single well energy type
    end do
!$OMP END DO

! NHT above are components in the basis (DH1,DH2,DH3)
! they need to be changed into cartesian coordinates to be used
! in calculation of force and torque

!$OMP DO PRIVATE(iPt)
   do iPt = 1,nPt-1
      NH(1,iPt+1) = NHT(1,iPt)*DH1(1,iPt) &
                   +NHT(2,iPt)*DH2(1,iPt) &
                   +NHT(3,iPt)*DH3(1,iPt)

      NH(2,iPt+1) = NHT(1,iPt)*DH1(2,iPt) &
                   +NHT(2,iPt)*DH2(2,iPt) &
                   +NHT(3,iPt)*DH3(2,iPt)

      NH(3,iPt+1) = NHT(1,iPt)*DH1(3,iPt) &
                   +NHT(2,iPt)*DH2(3,iPt) &
                   +NHT(3,iPt)*DH3(3,iPt)
   end do
!$OMP END DO
!$OMP END Parallel

    kkh=KH

!!! NHT above are components in the basis (DH1,DH2,DH3)
!!! they need to be changed into cartesian coordinates to be used
!!! in calculation of force and torque

!!! by reflection
    NH(:,1)=Tan
    NH(1,nPt+1) = 0.0_dp*NH(1,nPt)
    NH(2,nPt+1) = 0.0_dp*NH(2,nPt)
    NH(3,nPt+1) = 0.0_dp*NH(3,nPt)

  end subroutine InternalMomentAtHalfPt

!!!=======================================================================

  subroutine MomentnForce(coord_pt, NH, FH, TC, moment_pt, force_pt)
    use nrtype; use commonval
    implicit none

!!! Input arguments

!!! coord_pt - vector of boundary points and triads
!!!          - use only first nPt data
!!! FH - force transmitted across a cross section at each half point
!!!    - cartesian coordinates
!!! NH - moment transmitted across a cross section at each half point
!!!    - cartesian coordinates

!!! Output arguments

!!! moment_pt - applied moment
!!! force_pt - applied force

!!! Variables declaration
    real(dp), dimension(:,:), intent(in) :: coord_pt
    real(dp), dimension(:,:), intent(in) :: NH, FH
    real(dp), dimension(:,:), intent(in) :: TC
    real(dp), dimension(:,:), intent(out) :: moment_pt, force_pt

!!! local variables
    integer(i4b) :: iPt, iDim

!!! initialize
    force_pt=0.0_dp
    moment_pt=0.0_dp

!!! compute the appled moment to the filemant


!$OMP Parallel DO PRIVATE(iPt)
    do iPt = 1, nPt
       moment_pt(1,iPt)= (NH(1,iPt+1)-NH(1,iPt))/ds   &
              +( (TC(2,iPt+1)*FH(3,iPt+1)-TC(3,iPt+1)*FH(2,iPt+1)) &
              +(TC(2,iPt)*FH(3,iPt)-TC(3,iPt)*FH(2,iPt)) )/2.0_dp

       moment_pt(2,iPt)= (NH(2,iPt+1)-NH(2,iPt))/ds   &
              +( (TC(3,iPt+1)*FH(1,iPt+1)-TC(1,iPt+1)*FH(3,iPt+1)) &
              +(TC(3,iPt)*FH(1,iPt)-TC(1,iPt)*FH(3,iPt)) )/2.0_dp

       moment_pt(3,iPt)= (NH(3,iPt+1)-NH(3,iPt))/ds   &
              +( (TC(1,iPt+1)*FH(2,iPt+1)-TC(2,iPt+1)*FH(1,iPt+1)) &
              +(TC(1,iPt)*FH(2,iPt)-TC(2,iPt)*FH(1,iPt)) )/2.0_dp

      force_pt(:,iPt)=(FH(:,iPt+1)-FH(:,iPt))/ds
     end do
!$OMP END Parallel DO

1010 FORMAT (3(3X, F20.16))

  end subroutine MomentnForce

!!!=======================================================================
!!!=======================================================================

!!$subroutine WriteEnergy(iStep, klok, FHT, NHT,allEne)
!!$use nrtype; use commonval
!!$implicit none
!!$integer(i4b), intent(in) :: iStep, klok
!!$real(dp), dimension(:,:), intent(in) :: FHT, NHT
!!$real(dp), dimension(:,:), intent(inout) :: allEne
!!$real(dp) :: bendE, twistE, shearE, stretchE
!!$integer(i4b) :: i
!!$
!!$bendE=ds*(1.0_dp/2.0_dp)*(sum(NHT(1,:)**2/arry_bend1)+sum(NHT(2,:)**2/arry_bend2))
!!$twistE=ds*(1.0_dp/2.0_dp)*(sum(NHT(3,:)**2/arry_twist))
!!$shearE=ds*(1.0_dp/2.0_dp)*(sum(FHT(1,:)**2/shear)+sum(FHT(2,:)**2/shear))
!!$stretchE=ds*(1.0_dp/2.0_dp)*(sum(FHT(3,:)**2/stretch))
!!$
!!$!--------------------------
!!$! write the energy
!!$!-------------------------
!!$
!!$  i=klok+1
!!$  allEne(i,1)=iStep*dt
!!$  allEne(i,2)= bendE
!!$  allEne(i,3)= twistE
!!$  allEne(i,4)= shearE
!!$  allEne(i,5)= stretchE
!!$
!!$end subroutine WriteEnergy

!!!============================================================
end subroutine calcForceMoment_mot
