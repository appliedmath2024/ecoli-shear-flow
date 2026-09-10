!   nohup ./Stokes >& m.out &

PROGRAM rod_body
  USE omp_lib
  USE nrtype; USE interface_rod; USE commonval
  IMPLICIT NONE

  REAL(DP) :: Simulation_Time = 3.0_dp
  INTEGER(I4B) :: Iter_S, Iter_N, nSkip, frame, nSkip2, nSkip3
  INTEGER(I4B) :: i, j, k, iStep, kink_speed, oneP_num
  REAL(DP) :: tSwitch1, tSwitch2, tSwitch3, tSwitch4, tSwitch5, tSwitch6, temp, hook_rate
  INTEGER(I4B) :: switch1, switch2, switch3, switch4, switch5, switch6
  REAL(DP), DIMENSION(nrod) :: rod_height, rod_theta, body_height

  REAL(DP), DIMENSION(3,nPt) :: coord_pt0, D10, D20, D30
  REAL(DP), DIMENSION(3*nrod) :: temp3
  REAL(DP), DIMENSION(2*nrod) :: temp2
  REAL(DP), DIMENSION(nrod) :: temp1
  REAL(DP), DIMENSION(3,nPt,nrod) :: coord_pt, D1, D2, D3
  REAL(DP), DIMENSION(3,nPt,nrod) :: moment_pt, force_pt, vel_pt, ang_pt
  REAL(DP), DIMENSION(nPt):: KAPPA, Tau, KAPPA2, Tau2
  REAL(DP), DIMENSION(3,nPt-1,nrod):: kkh


  REAL(DP), DIMENSION(nPt-1) :: L0, ang
!!$  REAL(DP), DIMENSION(5) :: allEne

!!!! time bend twist shear stretch !!!!
  INTEGER(I4B) :: klok
  real(dp) :: start_time, end_time, z0, stiff, radius, wave, pitch, stiff_tan, angw, dist1
  real(DP) :: elapsed, avg_per_iter, remaining
  character(len=100) :: data_dir
!!!!!!!!!!!!!!!!!!!for body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(I4B) :: nref, kk, wall_switch
  REAL(DP), DIMENSION(:,:), Allocatable  :: LX, LF, LY, LU
  REAL(DP), DIMENSION(3) :: Tan, Tan2, Kforc, PSY, Torq, TTorq, Forc, FForc, &
       YCM, VCM, YCMold, Omega, normal, YCMn, YCMn2
  REAL(DP), DIMENSION(22+3*nrod) :: AD
  REAL(DP), DIMENSION(13+nrod) :: AD2
  REAL(DP), DIMENSION(3,nrod) :: axis_pt, axis_ct, axis_ct_move, Torq_ct
  REAL(DP), DIMENSION(3,3)  :: Fram, Framn, Framn2, Rtt, RotMat, PsMat
  REAL(DP), DIMENSION(3,3,nrod)  :: Trid, Tridn, Tridn2
  REAL(DP) :: detA, theta, walldist

  ! background flow
  REAL(DP) :: tube_R, tube_L, sconst, PoisU, Uflow

  !!!!!!!!!! for velofield !!!!!!!!!!
  INTEGER(I4B), PARAMETER :: nmk1=41, nmk2=41;
  INTEGER(I4B), PARAMETER :: nMk = nmk1*nmk2
  REAL(DP), DIMENSION(3) :: vecmk1, vecmk2, oript
  REAL(DP) :: a1, a2, a3, b1, b2, b3, width, interval
  REAL(DP), DIMENSION(3,nMk,3) :: coord_mk, vel_mk
  
!!!====================================================
!!!====================================================
  call OMP_set_NUM_THREADS(16);
!!!====================================================

  wall_switch=1; !!0:off, 1:on
  
  klok = 0
!ds=Length/(nPt-1.0_dp)  ! Warning: should be nPt-1
  hand= -1.0_dp        ! left-handed (-1), right-handed (+1)
  radius=0.35_dp   ! helical radius
  pitch=2.5_dp
  wave=hand*2.0_dp*pi/pitch    ! wave number: 2*pi/|w| - wave length
  rot = 1.0_dp   ! direction of rotation: CCW (+1), CW (-1)
  rot2 = rot        ! direction of rotation for the last flagellum

!!!====================================================
!!! parameter values
!!!====================================================
  dt=0.00000004_dp

  nmot=1   ! number of fixed points/motor points
  nhook=3
  freq=100.0_dp/1.0_dp
  freq_body=0.0_dp
  bend1=1.0_dp*10.0_dp**(-3)*3.5_dp
  bend2=bend1
  twist=bend1*1.0_dp
  shear=1.0_dp
  stretch=shear*1.0_dp

  curva1=0.0_dp      !radius / (radius**2+(pitch/2.0_dp/pi)**2)
  curva2=0.0_dp      ! intrinsic curvature of D^2
  curva3=hand*pitch/2.0_dp/pi / (radius**2 + (pitch/2.0_dp/pi)**2)

  stiff=20.0_dp ; stiff_tan=10.0_dp;  c_mss=2.0_dp;
  stiff_contact=100.0_dp; Dmin=0.09_dp; inv_Dmin=1.0_dp/Dmin
  stiff_contact_wall=100.0_dp; Dmin_wall=0.09_dp; inv_Dmin_wall=1.0_dp/Dmin_wall
  alpha=0.0_dp*10.0_dp**1; beta=alpha; ! alpha=1/(alpha*ds)
  hook_rate=20.0_dp

  uflow=0.0_dp
  walldist=0.8_dp; sconst=10.0_dp;
  PoisU=sconst*50.0_dp/4.0_dp; !! radius=50/2
!!! ==============================================================
!!! Initial Configuration and the rest length
!!!===============================================================

!!!====================================================
!!! cell Body
!!!====================================================
  R_cell=0.5_dp; H_cell=1.4_dp; A_cell=0.0_dp !! cap length, the other length (pitch), helical radius
  axis_x=R_cell;  axis_z=axis_x; axis_y=axis_x
!!$  body_height=R_cell + 0.5_dp*H_cell
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! symmetry setting test
!!$  rod_height(1)=0.5_dp*H_cell+R_cell
!!$  rod_height(2:nrod)=0.5_dp*H_cell-0.05_dp !+R_cell*cos(89.0_dp*pi/180_dp)
!!$
!!$  rod_theta(1)=0.0_dp;
!!$  do i=2,nrod;
!!$     rod_theta(i)=2.0_dp*(i-2.0_dp)*pi/(1.0_dp*(nrod-1));
!!$  end do

  rod_height(1)=0.5_dp*H_cell-0.05_dp !+R_cell*cos(89.0_dp*pi/180_dp)
  rod_height(2)=0.5_dp*H_cell+R_cell*cos(15.0_dp*pi/180_dp);
  rod_height(3)=0.5_dp*H_cell-0.05_dp !+R_cell*cos(89.0_dp*pi/180_dp);
  rod_height(4)=0.5_dp*H_cell+R_cell*cos(15.0_dp*pi/180_dp);

  do i=1,nrod;
     rod_theta(i)=2.0_dp*(i-1.0_dp)*pi/(1.0_dp*nrod);
  end do

!!$  theta=30
!!$  rod_height=0.5_dp*H_cell+cos(theta*pi/180_dp)*R_cell
!!$  
!!$  do i=1,nrod;
!!$     rod_theta(i)=2.0_dp*(i-1.0_dp)*pi/(1.0_dp*nrod);
!!$  end do
  
  oneP_num = nint(0.01_dp/dt)
  nSkip = nint(oneP_num/4.0_dp)
  nSkip2 = nint(oneP_num/20.0_dp)
  nSkip3 = nint(oneP_num/20.0_dp)
  
  kink_speed = nint(10.0_dp*oneP_num/npt)
  do i=1,nPt
     switch(i)=(i-1)*kink_speed
  end do

  data_dir = '../Data/restart_data/'
  frame=80  !!40=0.1s
!!$  Iter_S = frame*nSkip
  Iter_S = 0
  Iter_N = nint(Simulation_Time/dt)

  call initial_higdon(coord_pt0,D10,D20,D30,L0,radius,wave,kappa,kappa2,tau,tau2); del=3.0_dp*ds
!!Bend array
  arry_bend1=bend1; arry_bend2=bend2; arry_twist=twist;
  do i=1,nmot+nhook-1
        arry_bend1(i)=bend1/hook_rate
        arry_bend2(i)=bend2/hook_rate
        arry_twist(i)=twist*1.0_dp
  end do
  arry_bend1(nmot+nhook-1)=(arry_bend1(nmot+nhook-1)+arry_bend1(nmot+nhook))/2.0_dp
  arry_bend2(nmot+nhook-1)=(arry_bend2(nmot+nhook-1)+arry_bend2(nmot+nhook))/2.0_dp
  arry_twist(nmot+nhook-1)=(arry_twist(nmot+nhook-1)+arry_twist(nmot+nhook))/2.0_dp

  print*, 'nPt, ds, H, del', nPt, ds, H, del
  print*, 'bend, twist, shear, stiff', bend1, twist, shear, stiff
  print*, 'freq, Lenght, mu, rho', freq, Length, mu, rho
  print*, 'nhook, rot', nhook, rot


  nref=3;
!!$  if (H_cell.eq.0.0_dp) then
!!$     Call initial_body_sphere(nref);
!!$  else
     Call initial_body(nref);
!!$  end if
  KtmaxM=Ktmax; KvmaxM=Kvmax ! for movie
  Allocate(IndexM(KtmaxM,3)); ! for movie
  IndexM=Index0; Deallocate(Area0);
  Deallocate(Index0); Deallocate(LX0); Deallocate(Length0); Deallocate(Radius0) ! for movie
  nref=3;
!!$  if (H_cell.eq.0.0_dp) then
!!$     Call initial_body_sphere(nref);
!!$  else
     Call initial_body(nref);
!!$  end if
  Allocate(LX(3,Kvmax)); Allocate(LY(3,Kvmax));  Allocate(LU(3,Kvmax))
  Allocate(Coeff(3,Kvmax)); Allocate(LF(3,Kvmax)); Allocate(Mss(Kvmax))
  do i=1,3; LX(i,1:Kvmax)=LX0(:,i);
  end do

  YCM=0.0_dp
  Fram=0.0_dp; Fram(1,1)=1.0_dp; Fram(2,2)=1.0_dp; Fram(3,3)=1.0_dp; !!initial frame
  Coeff=LX

  do i=1,Kvmax; LY(:,i)=YCM+matmul(Fram,Coeff(:,i));
  end do
     
!!==========================================================================

!!! attach the filaments to the body;
     !  do j=1,3; coord_pt(j,:,1)=coord_pt0(j,:)-coord_pt0(j,1)+LX(j,1)+ds/2.0_dp*Fram(j,3); enddo
     !  D1(:,:,1)=D10; D2(:,:,1)=D20; D3(:,:,1)=D30
     !  axis_ct=0.0_dp ! (0,0,0) is center

  do i=1,nrod;
     if (rod_height(i).ge.0.5_dp*H_cell) then
        axis_pt(1,i)=sqrt(axis_z**2-(rod_height(i)-0.5_dp*H_cell)**2)*cos(rod_theta(i))
        axis_pt(2,i)=sqrt(axis_z**2-(rod_height(i)-0.5_dp*H_cell)**2)*sin(rod_theta(i))
        axis_pt(3,i)=rod_height(i)-0.5_dp*H_cell
        Tan=(/axis_pt(1,i)/axis_x**2,axis_pt(2,i)/axis_x**2,axis_pt(3,i)/axis_z**2/)
     elseif (rod_height(i).le.-0.5_dp*H_cell) then
        axis_pt(1,i)=sqrt(axis_z**2-(rod_height(i)+0.5_dp*H_cell)**2)*cos(rod_theta(i))
        axis_pt(2,i)=sqrt(axis_z**2-(rod_height(i)+0.5_dp*H_cell)**2)*sin(rod_theta(i))
        axis_pt(3,i)=rod_height(i)+0.5_dp*H_cell
        Tan=(/axis_pt(1,i)/axis_x**2,axis_pt(2,i)/axis_x**2,axis_pt(3,i)/axis_z**2/)
     else
        axis_pt(1,i)=axis_x*cos(rod_theta(i))
        axis_pt(2,i)=axis_x*sin(rod_theta(i))
        axis_pt(3,i)=0.0_dp
        Tan=(/axis_pt(1,i)/axis_x**2,axis_pt(2,i)/axis_x**2,axis_pt(3,i)/axis_z**2/)
     end if
     theta=sqrt(Tan(1)**2+Tan(2)**2+Tan(3)**2);
     axis_ct(:,i)=Tan/theta; print*, i, axis_pt(:,i), axis_ct(:,i)
  end do

  do i=1,nrod;
     !     Tan2=0.0_dp; Tan2(3)=1.0_dp; theta=0.0_dp
     D1(:,:,i)=D10; !do j=1,nPt; call rotation(D1(:,j,i),theta,Tan2); end do
     D2(:,:,i)=D20; !do j=1,nPt; call rotation(D2(:,j,i),theta,Tan2); end do
     D3(:,:,i)=D30; !do j=1,nPt; call rotation(D3(:,j,i),theta,Tan2); end do

     if (rod_height(i).le.-0.5_dp*H_cell) then
        coord_pt(1,:,i)=-coord_pt0(1,:)+coord_pt0(1,1);
        coord_pt(2,:,i)= coord_pt0(2,:)-coord_pt0(2,1);
        coord_pt(3,:,i)=-coord_pt0(3,:)+coord_pt0(3,1);
        do j=1,nPt; coord_pt(:,j,i)= coord_pt(:,j,i)-coord_pt(:,1,i)
        end do
        D1(1,:,i)=-D1(1,:,i); D1(3,:,i)=-D1(3,:,i)
        D2(1,:,i)=-D2(1,:,i); D2(3,:,i)=-D2(3,:,i)
        D3(1,:,i)=-D3(1,:,i); D3(3,:,i)=-D3(3,:,i)
        Tan=-axis_ct(:,i)
     else
        do j=1,nPt; coord_pt(:,j,i)= coord_pt0(:,j)-coord_pt0(:,1)
        end do
        Tan=axis_ct(:,i)
     end if

     if (abs(rod_height(i)).ne.0.5_dp*H_cell+R_cell) then
        Tan2(1)=Tan(3)*Fram(2,3)-Tan(2)*Fram(3,3)
        Tan2(2)=Tan(1)*Fram(3,3)-Tan(3)*Fram(1,3)
        Tan2(3)=Tan(2)*Fram(1,3)-Tan(1)*Fram(2,3)
        
        theta=sqrt(Tan2(1)**2+Tan2(2)**2+Tan2(3)**2)
        Tan2=Tan2/theta; theta=asin(theta);
        
        do j=1,nPt; call rotation(D1(:,j,i),theta,Tan2);
        end do
        do j=1,nPt; call rotation(D2(:,j,i),theta,Tan2);
        end do
        do j=1,nPt; call rotation(D3(:,j,i),theta,Tan2);
        end do
        do j=1,nPt; call rotation(coord_pt(:,j,i),theta,Tan2);
        end do;
     end if

     do j=1,3; coord_pt(j,:,i)=coord_pt(j,:,i)+axis_pt(j,i);
     end do
        
     if (rod_height(i).ge.0.5_dp*H_cell) then
        coord_pt(3,:,i)=coord_pt(3,:,i) + 0.5_dp*H_cell
        axis_pt(3,i)=axis_pt(3,i) + 0.5_dp*H_cell
     elseif (rod_height(i).le.-0.5_dp*H_cell) then
        coord_pt(3,:,i)=coord_pt(3,:,i) - 0.5_dp*H_cell
        axis_pt(3,i)=axis_pt(3,i) - 0.5_dp*H_cell
     else
        coord_pt(3,:,i)=coord_pt(3,:,i) + rod_height(i)
        axis_pt(3,i)=axis_pt(3,i) + rod_height(i)
     end if
     
     print*,'Dists=0 all;',axis_pt(:,i)-coord_pt(:,1,i)
  end do

  print*,'Dists=0 all;',maxval(abs(LX-LY)),LX(:,1)-coord_pt(:,1,1)
  print*,'L and YCM=', maxval(length0),minval(length0),sum(Coeff(1,:)),sum(Coeff(2,:)),sum(Coeff(3,:))

  ! rotation with tilted body with theta
!!$  theta=-10.0_dp*pi/180.0_dp; Tan=Fram(:,1)
!!$  do i=2,3; call rotation(Fram(:,i),theta,Tan);
!!$  end do
!!$  do k=1,nrod
!!$     do i=1,nPt; call rotation(D1(:,i,k),theta,Tan);
!!$     end do
!!$     do i=1,nPt; call rotation(D2(:,i,k),theta,Tan);
!!$     end do
!!$     do i=1,nPt; call rotation(D3(:,i,k),theta,Tan);
!!$     end do
!!$     do i=1,nPt; coord_pt(:,i,k)= coord_pt(:,i,k)-YCM
!!$        call rotation(coord_pt(:,i,k),theta,Tan);
!!$        coord_pt(:,i,k)= coord_pt(:,i,k)+YCM
!!$     end do
!!$  end do
!!$  do i=1,Kvmax; LX(:,i)=LX(:,i)-YCM
!!$     call rotation(LX(:,i),theta,Tan);
!!$     LX(:,i)=LX(:,i)+YCM
!!$  end do

  !! transliation due to wall
  if (wall_switch.eq.1) then
     LX(2,:)=LX(2,:)+walldist
     coord_pt(2,:,:)=coord_pt(2,:,:)+walldist
     YCM(2)=YCM(2)+walldist
  end if

  print*,'Dists=0 all;',maxval(abs(LX-LY)),LX(:,1)-coord_pt(:,1,1)
  print*,'L and YCM=', maxval(length0),minval(length0),sum(Coeff(1,:)),sum(Coeff(2,:)),sum(Coeff(3,:))

!!!====================================================
!!! open files
!!!====================================================
  if (frame .eq. 0) then
     open(20, file = 'coord_pt.m',status='replace')
     open(21, file = 'triad1.m',status='replace')
     open(22, file = 'triad2.m',status='replace')
     open(23, file = 'triad3.m',status='replace')
     open(30, file = 'body_pt.m',status='replace')
     open(31, file = 'Fram.m',status='replace')
     open(32, file = 'YCM.m',status='replace')
     open(33, file = 'VCM.m',status='replace')
     open(40, file = 'curv_tors.m',status='replace')
     open(35, file = 'AD.m',status='replace')
     open(36, file = 'AD2.m',status='replace')
  else
     open(20, file = trim(data_dir)//'coord_pt.m', status='old')
     open(21, file = trim(data_dir)//'triad1.m', status='old')
     open(22, file = trim(data_dir)//'triad2.m', status='old')
     open(23, file = trim(data_dir)//'triad3.m', status='old')
     open(30, file = trim(data_dir)//'body_pt.m', status='old')
     open(31, file = trim(data_dir)//'Fram.m', status='old')
     open(32, file = trim(data_dir)//'YCM.m', status='old')
     open(33, file = trim(data_dir)//'VCM.m', status='old')
     open(40, file = trim(data_dir)//'curv_tors.m', status='old')
     open(35, file = trim(data_dir)//'AD.m', status='old')
     open(36, file = trim(data_dir)//'AD2.m', status='old')
  end if

  open(34, file = 'body_face.m',status='replace')

!!$  open(50, file = 'coord_mk.m',status='replace')
!!$  open(51, file = 'vel_mk.m',status='replace')
!!$  open(52, file = 'vel_Fram.m',status='replace')

  open(60, file = 'force_pt.m',status='replace')
  
  open(97, file = 'CFL_condition.m',status='replace')
  open(99, file = 'status.m',status='replace')

  moment_pt=0.0_dp; force_pt=0.0_dp;
  vel_pt=0.0_dp; ang_pt=0.0_dp;
  VCM=0.0_dp; LU=0.0_dp;

  vel_pt(3,:,:) = Uflow
  vel_pt(3,:,:) = vel_pt(3,:,:) + sconst*coord_pt(2,:,:);
  LU(3,:) = LU(3,:) + sconst*LX(2,:);
  ang_pt(1,:,:) = ang_pt(1,:,:) + sconst/2.0_dp
  
  AD=0.0_dp;
  AD(5:7) = YCM;
  AD2=0.0_dp;

  if (frame .ne. 0) then
     do k=1,frame+1
        do i=1,nPt
           read(20,150) temp3
           do j=1,nrod; coord_pt(:,i,j)=temp3((j-1)*3+1:j*3)
           end do
           read(21,150) temp3
           do j=1,nrod; D1(:,i,j)=temp3((j-1)*3+1:j*3)
           end do
           read(22,150) temp3
           do j=1,nrod; D2(:,i,j)=temp3((j-1)*3+1:j*3)
           end do
           read(23,150) temp3
           do j=1,nrod; D3(:,i,j)=temp3((j-1)*3+1:j*3)
           end do
        end do

        do i=1,kvmaxM
           read(30,33) LX(:,i)
        end do

        read(31,99) Fram(:,1),Fram(:,2),Fram(:,3)
        read(32,33) YCM
        read(33,33) VCM
     end do
!!$  end if
!!$
!!$  if (frame.eq.80) then
     if (wall_switch.eq.1) then
  !! transliation due to wall
        LX(2,:)=LX(2,:)+walldist
        coord_pt(2,:,:)=coord_pt(2,:,:)+walldist
        YCM(2)=YCM(2)+walldist
     end if

     LX(1,:)=LX(1,:)-YCM(1)
     coord_pt(1,:,:)=coord_pt(1,:,:)-YCM(1)
     YCM(1)=YCM(1)-YCM(1)

     LX(3,:)=LX(3,:)-YCM(3)
     coord_pt(3,:,:)=coord_pt(3,:,:)-YCM(3)
     YCM(3)=YCM(3)-YCM(3)
  end if

  close(20); close(21); close(22); close(23);
  close(30); close(31); close(32); close(33);
  close(35); close(36); close(40);
  open(20, file = 'coord_pt.m', status='replace')
  open(21, file = 'triad1.m', status='replace')
  open(22, file = 'triad2.m', status='replace')
  open(23, file = 'triad3.m', status='replace')
  open(30, file = 'body_pt.m', status='replace')
  open(31, file = 'Fram.m', status='replace')
  open(32, file = 'YCM.m', status='replace')
  open(33, file = 'VCM.m', status='replace')
  open(40, file = 'curv_tors.m', status='replace')
  open(35, file = 'AD.m', status='replace')
  open(36, file = 'AD2.m', status='replace')

  Omega=0.0_dp;
  YCMn=YCM; Framn=Fram;
  YCMn2=YCM; Framn2=Fram;
  Trid=0.0_dp;
  do i=1,nrod; Trid(:,1,i)=D1(:,1,i); Trid(:,2,i)=D2(:,1,i); Trid(:,3,i)=D3(:,1,i);
  end do
  Tridn=Trid
  Tridn2=Trid

 
!!!====================================================
!!!  write initial data & setting
!!!====================================================
  write(99,*) 'bend1   = ',bend1,';'
  write(99,*) 'bend2   = ',bend2,';'
  write(99,*) 'twist   = ',twist,';'
  write(99,*) 'shear   = ',shear,';'
  write(99,*) 'stretch = ',stretch,';'
  write(99,*) 'freq    = ', freq,';'
  write(99,*) 'nPt   = ', nPt,';'
  write(99,*) 'nBd   = ', kvmaxM,';'
  write(99,*) 'nhook = ', nhook,';'
  write(99,*) 'nmot = ' , nmot,';'
  write(99,*) 'ds    = ', ds,';'
  write(99,*) 'h     = ', h,';'
  write(99,*) 'del   = ', del,';'
  write(99,*) 'tFinal   = ', Simulation_Time,';'
  write(99,*) 'nSkip    = ', nSkip,';'
  write(99,*) 'tnSkip    = ', nSkip*dt,';'
  write(99,*) 'nSkip2    = ', nSkip2,';'
  write(99,*) 'tnSkip2    = ', nSkip2*dt,';'
  write(99,*) 'tSwitch1 = ', tSwitch1,';'
  write(99,*) 'tSwitch2 = ', tSwitch2,';'
  write(99,*) 'tSwitch3 = ', tSwitch3,';'
  write(99,*) 'dt       = ', dt,';'
  write(99,*) 'radius       = ', radius,';'
  write(99,*) 'pitch       = ', pitch,';'
  write(99,*) 'curva1       = ', curva1,';'
  write(99,*) 'curva3       = ', curva3,';'
  write(99,*) 'Height = ' , Length,';'
  write(99,*) 'Length = ' , sum(L0(:)),';'
  write(99,*) 'fla_s = ' , coord_pt(3,1,nrod),';'
  write(99,*) 'fla_e = ' , coord_pt(3,nPt,nrod),';'
  write(99,*) 'nflag = ' , nrod,';'
  write(99,*) 'Domain = ' , Domain_x,';'
  write(99,*) 'Iter_N=', Iter_N,';'
  write(99,*) 'oneP_num=', oneP_num,';'
  write(99,*) 'kink_speed=', kink_speed,';'
  write(99,*) 'wall_switch=', wall_switch,';'
  if (wall_switch.eq.1) write(99,*) 'walldist=', walldist,';'
  
  do i=1,nrod
     write(99,*) 'rod_height=', rod_height(i),';'
     write(99,*) 'rod_theta=', rod_theta(i),';'
  end do

  do i=1,nrod
     write(99,*) 'block_angle=', acos(axis_ct(3,i))*180/pi,';'
  end do

!!$  if (Iter_S .eq. 0) then
     do i=1,nPt
        do j=1,nrod; temp3((j-1)*3+1:j*3)=coord_pt(:,i,j)
        end do
        write(20,150) temp3
        do j=1,nrod; temp3((j-1)*3+1:j*3)=D1(:,i,j)
        end do
        write(21,150) temp3
        do j=1,nrod; temp3((j-1)*3+1:j*3)=D2(:,i,j)
        end do
        write(22,150) temp3
        do j=1,nrod; temp3((j-1)*3+1:j*3)=D3(:,i,j)
        end do
        write(23,150) temp3
     end do

     do i=1,kvmaxM
        write(30,33) LX(:,i)
     end do

     write(31,99) Fram(:,1),Fram(:,2),Fram(:,3)
     write(32,33) YCM
     write(33,33) VCM

     do i=1,nPt-1
        do j=1,nrod
           temp2((j-1)*2+1)=kappa(i)
           temp2(j*2)=tau(i)
        end do
        write(40,150) temp2
     end do
     write(35,21) AD
     write(36,21) AD2
!!$  else
!!$     close(20); close(21); close(22); close(23);
!!$     close(30); close(31); close(32); close(33);
!!$     close(35); close(36); close(40);
!!$     open(20, file = 'coord_pt.m', status='old',position='append')
!!$     open(21, file = 'triad1.m', status='old',position='append')
!!$     open(22, file = 'triad2.m', status='old',position='append')
!!$     open(23, file = 'triad3.m', status='old',position='append')
!!$     open(30, file = 'body_pt.m', status='old',position='append')
!!$     open(31, file = 'Fram.m', status='old',position='append')
!!$     open(32, file = 'YCM.m', status='old',position='append')
!!$     open(33, file = 'VCM.m', status='old',position='append')
!!$     open(40, file = 'curv_tors.m', status='old',position='append')
!!$     open(35, file = 'AD.m', status='old',position='append')
!!$     open(36, file = 'AD2.m', status='old',position='append')
!!$  end if
  write(34,60) IndexM

  start_time = OMP_GET_WTIME()

!!!!======================================================
!!!!Iteration begins
!!!!======================================================

  do iStep = Iter_S+1, Iter_N
!!$  do iStep = 1,1
!!$     print *, iStep

     
     do i=1,nrod
        do kk=1,3; axis_ct_move(kk,i)=sum(Fram(kk,:)*axis_ct(:,i));
        end do
        Torq_ct(:,i)=-0.002_dp*axis_ct_move(:,i)
     end do



     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,1),D1(:,:,1),D2(:,:,1),D3(:,:,1), &
          moment_pt(:,:,1),force_pt(:,:,1),kappa,tau,kkh(:,:,1),Torq_ct(:,1))
     
     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,2),D1(:,:,2),D2(:,:,2),D3(:,:,2), &
          moment_pt(:,:,2),force_pt(:,:,2),kappa2,tau2,kkh(:,:,2),Torq_ct(:,2))
     
     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,3),D1(:,:,3),D2(:,:,3),D3(:,:,3), &
          moment_pt(:,:,3),force_pt(:,:,3),kappa2,tau2,kkh(:,:,3),Torq_ct(:,3))
 
     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,4),D1(:,:,4),D2(:,:,4),D3(:,:,4), &
          moment_pt(:,:,4),force_pt(:,:,4),kappa2,tau2,kkh(:,:,4),Torq_ct(:,4))
!!$
!!$     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,5),D1(:,:,5),D2(:,:,5),D3(:,:,5), &
!!$          moment_pt(:,:,5),force_pt(:,:,5),kappa2,tau2,kkh(:,:,5),Torq_ct(:,5))
!!$
!!$     call calcForceMoment_mot(iStep,klok,nSkip,coord_pt(:,:,6),D1(:,:,6),D2(:,:,6),D3(:,:,6), &
!!$          moment_pt(:,:,6),force_pt(:,:,6),kappa2,tau2,kkh(:,:,6),Torq_ct(:,6))

!!! cell body force
     Forc=0.0_dp; Torq=0.0_dp; LF=0.0_dp
     do i=1,Kvmax
        do kk=1,3; PSY(kk)=sum(Fram(kk,:)*Coeff(:,i));
        end do
        Kforc=c_mss*(LX(:,i)-PSY-YCM)
        LF(:,i)=-Kforc

        if (wall_switch.eq.1) then
           if (LX(2,i)<Dmin_wall) then;
              LF(2,i)=LF(2,i)+stiff_contact_wall*(1.0_dp-LX(2,i)*inv_Dmin_wall)
           end if
        end if
        
     end do

!!XXXXXXXXXXXXXX ADD these 10 lines with stiff=100.0 (may change) and use Forc and Torq
!!$     Forc=0.0_dp; Torq=0.0_dp;
!!$     do i=1,nrod
!!$        do kk=1,3; PSY(kk)=sum(Fram(kk,:)*axis_pt(:,i)); end do
!!$        Kforc=c_mss*(coord_pt(:,1,i)-(YCM+PSY))
!!$        Forc=Forc+Kforc
!!$        force_pt(:,1,i)=force_pt(:,1,i)-Kforc/ds
!!$        Torq(1)=Torq(1)+PSY(2)*kforc(3)-PSY(3)*kforc(2)
!!$        Torq(2)=Torq(2)+PSY(3)*kforc(1)-PSY(1)*kforc(3)
!!$        Torq(3)=Torq(3)+PSY(1)*kforc(2)-PSY(2)*kforc(1)
!!$     end do
     
     Call selfcontactForce(coord_pt, force_pt)
     Call selfcontactForce_Fla_Body(coord_pt, force_pt, LX, LF)
!!$      Call selfcontactForce_wall(coord_pt, force_pt, LX, LF)

!!! Need to have proper scale for moments and forces to send to Velocity Solves
     moment_Pt=moment_Pt*ds ; force_Pt=force_Pt*ds;
     do i=1,nrod
        Tan=axis_ct_move(:,i)
        Tan2=stiff_tan*(sum(D2(:,1,i)*Tan)*D1(:,1,i)-sum(D1(:,1,i)*Tan)*D2(:,1,i))*ds
        Torq=Torq+Tan2+Torq_ct(:,i)
        moment_Pt(:,1,i)=-Tan2+moment_Pt(:,1,i)
     end do

     if (wall_switch.eq.1) then
        do i=1,nrod; do j=1,nPt; !! repulsive force
           dist1=coord_pt(2,j,i);
           if (dist1<Dmin_wall) then;
              force_Pt(2,j,i)=force_Pt(2,j,i)+stiff_contact_wall*(1.0_dp-dist1*inv_Dmin_wall)
           end if
        end do; end do
     end if


     if (wall_switch.eq.0) then
        call StSolve(Coord_Pt,moment_Pt,force_Pt,vel_pt,LX,LF,LU)
        call AngSolve(Coord_Pt,moment_Pt,force_Pt,LX,LF,ang_Pt)
     else
        call StSolve_wall(Coord_Pt,moment_Pt,force_Pt,vel_pt,LX,LF,LU)
        call AngSolve_wall(Coord_Pt,moment_Pt,force_Pt,LX,LF,ang_Pt)
     end if

!!! background flow
     vel_pt(3,:,:) = vel_pt(3,:,:) + Uflow

!!! shear flow
     vel_pt(3,:,:) = vel_pt(3,:,:) + sconst*coord_pt(2,:,:);
     LU(3,:) = LU(3,:) + sconst*LX(2,:);
     ang_pt(1,:,:) = ang_pt(1,:,:) + sconst/2.0_dp
     
!!! check CFL condition
     if ( dt*maxval( abs(vel_pt(1,:,:))+abs(vel_pt(2,:,:))+abs(vel_pt(3,:,:)) ) > H ) then
        write(97,*) 'CFL violation'
        print *,  'CFL violation'
        write(97,*) 'dt*maxval(abs(U)+abs(V)+abs(W))',dt*maxval( abs(vel_pt(1,:,:))+abs(vel_pt(2,:,:))+abs(vel_pt(3,:,:)) )
        write(97,*) 'iStep', iStep
        write(97,*) '***************************************'
        stop
     end if

     coord_pt=coord_pt+dt*(vel_pt + force_pt*alpha)  !! slip
     LX = LX + dt*(LU+LF*alpha)

     Call frame_CSP(LX,Forc,Torq,YCM,Fram)
!!$   Fram=Framn;
!!$   YCM=YCMn;

!!XXXXXXXXXXXXXX remove these 5 lines
     do i=1,nrod;
        do kk=1,3; PSY(kk)=sum(Fram(kk,:)*axis_pt(:,i));
        end do
        coord_pt(:,1,i)=YCM+PSY;
     end do

     do i=1,nrod
!!!-----------------------------------
!!! add torque slip
!!!-----------------------------------
        do j = 1, nPt
           ang_pt(:,j,i)= ang_pt(:,j,i)+ moment_pt(3,j,i)*D3(:,j,i)*beta
        enddo
        
        call updateTriad_mot (iStep,ang_pt(:,:,i),Tan,D1(:,:,i),D2(:,:,i),D3(:,:,i))
     end do

!!! write files
     if (mod(iStep,nSkip)==0) then
        klok = klok + 1

        AD(1)=dt*iStep
        do i=1,Kvmax
           do kk=1,3
              PSY(kk)=sum(Fram(kk,:)*Coeff(:,i))
           end do
           LY(:,i)=YCM+PSY
        end do
        AD(2) =  maxval(abs(LX-LY))
        AD(3) = 0.0_dp; AD(4) =  0.0_dp
        do i=1,nrod
           do kk=1,3; PSY(kk)=sum(Fram(kk,:)*axis_pt(:,i)); end do
           dist1=maxval(abs(coord_pt(:,1,i)-(YCM+PSY)))
           AD(3) = max(dist1, AD(3))
           Tan=axis_ct_move(:,i)
           dist1=maxval(abs(Tan-D3(:,1,i)))
           AD(4) = max(dist1, AD(4))
        end do
        AD(5:7) = YCM;
        VCM=(YCM-YCMn)/(dt*nSkip);
        AD(8:10) = VCM; YCMn=YCM

        do i=1,3;
           do j=1,3; Rtt(i,j)=sum(Fram(i,:)*Framn(j,:)); enddo;
        enddo
        theta=acos((Rtt(1,1)+Rtt(2,2)+Rtt(3,3)-1.0_dp)/2.0_dp)
        detA=theta/(dt*nSkip);
        Omega(1)=(Rtt(3,2)-Rtt(2,3))*detA/(2.0_dp*sin(theta))
        Omega(2)=(Rtt(1,3)-Rtt(3,1))*detA/(2.0_dp*sin(theta))
        Omega(3)=(Rtt(2,1)-Rtt(1,2))*detA/(2.0_dp*sin(theta))
        AD(11:13) = Omega;
        AD(14)=sum(Fram(:,3)*Torq);
        AD(15)=sum(Fram(:,3)*Forc);
        AD(16)=sum(Fram(:,3)*Omega)/(2.0_dp*pi)
        AD(17)=sum(Fram(:,3)*VCM)

        AD(18:20)=Fram(:,3)

        do i=1,3; Omega(i)=sum(LF(i,:)); end do
        AD(21)=sum(Fram(:,3)*Omega);
        AD(22)=90.0_dp-acos(-Fram(2,3))*180.0_dp/pi

        do i=1,nrod; Trid(:,1,i)=D1(:,1,i); Trid(:,2,i)=D2(:,1,i); Trid(:,3,i)=D3(:,1,i); end do
!!        do i=1,3; PsMat(:,i)=Framn(i,:); end do
!!        RotMat=matmul(Fram, PsMat)
        do k=1,nrod
!!           PsMat=matmul(RotMat,Tridn(:,:,k)); Tridn(:,:,k)=PsMat
           do i=1,3;
              do j=1,3; Rtt(i,j)=sum(Trid(i,:,k)*Tridn(j,:,k));
              enddo;
           enddo
           theta=acos((Rtt(1,1)+Rtt(2,2)+Rtt(3,3)-1.0_dp)/2.0_dp)
           detA=theta/(dt*nSkip);
           Omega(1)=(Rtt(3,2)-Rtt(2,3))*detA/(2.0_dp*sin(theta))
           Omega(2)=(Rtt(1,3)-Rtt(3,1))*detA/(2.0_dp*sin(theta))
           Omega(3)=(Rtt(2,1)-Rtt(1,2))*detA/(2.0_dp*sin(theta))
           AD(22+k)=sum(Omega*Tridn(:,3,k))/(2.0_dp*pi)
        end do
        Framn=Fram; Tridn=Trid
        do k=1,nrod
          do i=1,3; Omega(i)=sum(force_pt(i,:,k)); end do
          AD(22+nrod+k)=sum(Fram(:,3)*Omega)
          Omega=-force_Pt(:,1,k)
          AD(22+2*nrod+k)=sum(Fram(:,3)*Omega)
        end do

        write(35,21) AD
!!$        print*, 'AD', iStep, AD,rot
!!!!!!!!!!!!!!!!!!!!!!

        do i=1,nPt
           do j=1,nrod; temp3((j-1)*3+1:j*3)=coord_pt(:,i,j)
           end do
           write(20,150) temp3
           do j=1,nrod; temp3((j-1)*3+1:j*3)=D1(:,i,j)
           end do
           write(21,150) temp3
           do j=1,nrod; temp3((j-1)*3+1:j*3)=D2(:,i,j)
           end do
           write(22,150) temp3
           do j=1,nrod; temp3((j-1)*3+1:j*3)=D3(:,i,j)
           end do
           write(23,150) temp3
        end do

        do i=1,kvmaxM
           write(30,33) LX(:,i)
        end do

        write(31,99) Fram(:,1),Fram(:,2),Fram(:,3)
        write(32,33) YCM
        write(33,33) VCM

        do i=1,nPt-1
           do j=1,nrod
              temp2((j-1)*2+1)=sqrt(kkh(1,i,j)**2+kkh(2,i,j)**2)
              temp2(j*2)=kkh(3,i,j)
           end do
           write(40,150) temp2
        end do
       
     end if

     if (mod(iStep,1*nSkip)==0) then
        print*, 'AD', iStep, AD;

        end_time = OMP_GET_WTIME()
        elapsed = end_time-start_time
        avg_per_iter = elapsed / iStep
        remaining = avg_per_iter * (Iter_N - iStep)

        print *, "remaining time: ", remaining/60., '(m)'
        print *, "remaining time: ", remaining/3600., '(h)'
     endif
  

!!!======================================================
  end do ! iStep
!!!======================================================
21 format(50(e23.16,1x))
33 format(3(e23.16,1x))
55 format(8(e23.16,1x))
99 format(9(e23.16,1x))
100 format(10(e23.16,1x))
150 format(50(e23.16,1x))
60 format(i7)

  !call cpu_time(end)
  end_time = OMP_GET_WTIME()
  print *, "Elapsed time = ", end_time-start_time, '(s)'
  print *, "Elapsed time = ", (end_time-start_time)/60., '(m)'
  print *, "Elapsed time = ", (end_time-start_time)/3600., '(h)'

  !!======================================================
!!!!Iteration ends
  !!======================================================

END PROGRAM
