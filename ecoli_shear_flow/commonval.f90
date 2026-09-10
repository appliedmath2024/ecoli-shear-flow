MODULE commonval
!!! Global  parameters and variables
!!! Units: microns, gram, second
  USE nrtype
  IMPLICIT NONE
  REAL(DP), PARAMETER :: rho=1.0_dp*10.0_dp**(-12), mu=0.01_dp*10.0_dp**(-4)
  REAL(DP), PARAMETER :: Domain_x=10.0_dp, Length=8.5_dp
  INTEGER(I4B), PARAMETER :: nPt=374

  INTEGER(I4B), PARAMETER :: Expmt=2, nrod=4 !! symmetric=4+1
  INTEGER(I4B), PARAMETER :: N1=64*Expmt, N2=64*Expmt, N3=64*Expmt
  REAL(DP), PARAMETER :: H=Domain_x/N3               !mesh width

  REAL(DP) :: ds, del, R_cell, H_cell, A_cell ! Warning: should be nPt-1

  REAL(DP) ::  rot, rot2, hand, dt, freq, freq_body, bend1, bend2, twist, shear, stretch, sigma, &
               stiff_contact, Dmin, inv_Dmin, stiff_contact_wall, Dmin_wall, inv_Dmin_wall
  REAL(DP) :: curva1, curva2, curva3, eta, alpha, beta     ! intrinsic properties
  REAL(DP), DIMENSION(nPt):: d2tau ! intrinsic properties
  INTEGER(I4B), DIMENSION(nPt):: switch

  REAL(DP), DIMENSION(nPt-1):: arry_bend1, arry_bend2, arry_twist
  REAL(DP), PARAMETER :: epsil = 0.001_dp
  INTEGER(I4B) :: nmot, nhook  ! CCW(+1) or CW(-1)

  INTEGER(I4B), DIMENSION(nrod) :: Index_rod
  REAL(DP), DIMENSION(:,:), Allocatable  :: LX0, Coeff, Length0
  REAL(DP), DIMENSION(:), Allocatable  :: Radius0, Area0
  INTEGER(I4B), DIMENSION(:,:), Allocatable  :: Index0, IndexM
  INTEGER(I4B) :: Kvmax, Ktmax, KvmaxM, KtmaxM
  REAL(DP) :: TMss, c_mss, c_vel, axis_x, axis_y, axis_z
  REAL(DP), DIMENSION(:), Allocatable :: Mss
  REAL(DP), DIMENSION(3,3) :: MI, IMI

END MODULE commonval
