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

  angw = 2.0_dp*pi*freq!!*(1.0_dp-exp(-dt*iStep/0.0001_dp))
!!!smooth=0.1_dp

!!!----------------------------------------------------------------
!!! update triads of motors
!!!----------------------------------------------------------------

!!!----------------------------------------------------------------
!!! update triads of non-motors
!!!----------------------------------------------------------------
!!! calculate tangent vectors
  !$OMP PARALLEL
  !$OMP DO PRIVATE(iPt,temp1,temp2,temp3)
  do iPt = 1, nPt

     normw(iPt) = sqrt( ang_pt(1,iPt)**2+ang_pt(2,iPt)**2 &
          +ang_pt(3,iPt)**2 )
     invnormw(iPt)=1.0_dp/normw(iPt)

     temp1 = D1(1,iPt)*ang_pt(1,iPt)+D1(2,iPt)*ang_pt(2,iPt) &
          +D1(3,iPt)*ang_pt(3,iPt)

     DT1(:,iPt)=temp1*ang_pt(:,iPt)*(invnormw(iPt)**2)	

     temp2 = D2(1,iPt)*ang_pt(1,iPt)+D2(2,iPt)*ang_pt(2,iPt) &
          +D2(3,iPt)*ang_pt(3,iPt)

     DT2(:,iPt)=temp2*ang_pt(:,iPt)*(invnormw(iPt)**2)

     temp3 = D3(1,iPt)*ang_pt(1,iPt)+D3(2,iPt)*ang_pt(2,iPt) &
          +D3(3,iPt)*ang_pt(3,iPt)

     DT3(:,iPt)=temp3*ang_pt(:,iPt)*(invnormw(iPt)**2)

  end do
  !$OMP END DO

  !$OMP DO PRIVATE(iPt)
  do  iPt = 1, nPt
     wxD1(1,iPt)=ang_pt(2,iPt)*D1(3,iPt)-ang_pt(3,iPt)*D1(2,iPt)
     wxD1(2,iPt)=ang_pt(3,iPt)*D1(1,iPt)-ang_pt(1,iPt)*D1(3,iPt)
     wxD1(3,iPt)=ang_pt(1,iPt)*D1(2,iPt)-ang_pt(2,iPt)*D1(1,iPt)

     wxD2(1,iPt)=ang_pt(2,iPt)*D2(3,iPt)-ang_pt(3,iPt)*D2(2,iPt)
     wxD2(2,iPt)=ang_pt(3,iPt)*D2(1,iPt)-ang_pt(1,iPt)*D2(3,iPt)
     wxD2(3,iPt)=ang_pt(1,iPt)*D2(2,iPt)-ang_pt(2,iPt)*D2(1,iPt)

     wxD3(1,iPt)=ang_pt(2,iPt)*D3(3,iPt)-ang_pt(3,iPt)*D3(2,iPt)
     wxD3(2,iPt)=ang_pt(3,iPt)*D3(1,iPt)-ang_pt(1,iPt)*D3(3,iPt)
     wxD3(3,iPt)=ang_pt(1,iPt)*D3(2,iPt)-ang_pt(2,iPt)*D3(1,iPt)
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(iPt,cosval,sinval)
  do  iPt = 1, nPt; 	
     cosval=cos(normw(iPt)*dt); sinval=sin(normw(iPt)*dt)
     DN1(:,iPt) = cosval*D1(:,iPt)+(1.0_dp-cosval)*DT1(:,iPt)  &
          + sinval*wxD1(:,iPt)*invnormw(iPt)

     DN2(:,iPt) = cosval*D2(:,iPt)+(1.0_dp-cosval)*DT2(:,iPt) &
          +sinval*wxD2(:,iPt)*invnormw(iPt)

     DN3(:,iPt) = cosval*D3(:,iPt)+(1.0_dp-cosval)*DT3(:,iPt) &
          +sinval*wxD3(:,iPt)*invnormw(iPt)
  end do;
  !$OMP END DO

  !$OMP END PARALLEL

  D1=DN1; D2=DN2; D3=DN3

end subroutine updateTriad_mot
