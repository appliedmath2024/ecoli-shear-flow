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

  k=2.0_dp

  Tspan=(/ (Length*(i-1.0_dp)/(10000.0_dp*Expmt), i=1,10000*Expmt+1 )/)
  Sspan(1)=0.0_dp;
  do i=1,10000*Expmt
     Sspan(i+1)=Sspan(i)+(0.5_dp)*(Length/(10000.0_dp*Expmt))&
          *(sqrt( (2.0_dp*r0*k*Tspan(i)*exp(-k*Tspan(i)*Tspan(i)))**2     &
          +(r0*w*(1.0_dp-exp(-k*Tspan(i)*Tspan(i))))**2+1.0_dp)    &
          + sqrt( (2.0_dp*r0*k*Tspan(i+1)*exp(-k*Tspan(i+1)*Tspan(i+1)))**2     &
          +(r0*w*(1.0_dp-exp(-k*Tspan(i+1)*Tspan(i+1))))**2+1.0_dp) )
  end do
  print*,Sspan(10000*Expmt+1)

  ds=Sspan(10000*Expmt+1)/(1.0_dp*(nPt-nmot-nhook))
  do i=1, nPt; TT(i)=ds*(i-1) ;
  end do

  s(nmot+nhook)=0.0_dp; s(nPt)=Length
  do i=2,nPt-nmot-nhook;
     do j=1,10000*Expmt
        if ((Sspan(j)<=TT(i)) .AND. (TT(i)<Sspan(j+1))) then
           s(nmot+nhook+i-1)=(Tspan(j+1)-Tspan(j))*(TT(i)-Sspan(j))/&
                (Sspan(j+1)-Sspan(j))+Tspan(j)
        end if;
     end do;
  end do
  
!!! initialize
  dx=0.0_dp; dy=0.0_dp; dz=0.0_dp; r=0.0_dp; dr=0.0_dp; ddr=0.0_dp; dddr=0.0_dp
  ddx=0.0_dp; ddy=0.0_dp; ddz=0.0_dp; dddx=0.0_dp; dddy=0.0_dp; dddz=0.0_dp
  L_X=0.0_dp
  
!!! Lagrangian parameter
  do i=1, nmot+nhook-1
     L_X(3,i)=ds*(i-(nmot+nhook))
  end do
  
  do i = nmot+nhook, nPt
     L_X(3,i)=s(i)
     r(i)=r0*(1.0_dp-exp(-k*s(i)*s(i)))
     L_X(1,i)=r(i)*cos(w*s(i-nmot+1))
     L_X(2,i)=r(i)*sin(w*s(i-nmot+1))
  end do
  
  do i=1,nPt-1
     L0(i)=sqrt( (L_X(1,i+1)-L_X(1,i))**2 &
          +(L_X(2,i+1)-L_X(2,i))**2 &
          +(L_X(3,i+1)-L_X(3,i))**2 )
  end do
!!!print*,L0, 10000

  do i=nmot+nhook,nPt
     dr(i)=2.0_dp*r0*k*s(i)*exp(-k*s(i)*s(i))
     
     dx(i)= dr(i)*cos(w*s(i))-r(i)*w*sin(w*s(i))
     dy(i)= dr(i)*sin(w*s(i))+r(i)*w*cos(w*s(i))
     dz(i)= 1.0_dp;
     
     deno=sqrt(dx(i)**2+dy(i)**2+dz(i)**2)
     
     D3(1,i)=dx(i)/deno
     D3(2,i)=dy(i)/deno
     D3(3,i)=dz(i)/deno
     
     ddr(i)=2.0_dp*r0*k*exp(-k*s(i)*s(i))*(1.0_dp-2.0_dp*k*s(i)*s(i))
     
     ddx(i)=ddr(i)*cos(w*s(i))-2.0_dp*w*dr(i)*sin(w*s(i))-r(i)*w**2*cos(w*s(i))
     ddy(i)=ddr(i)*sin(w*s(i))+2.0_dp*w*dr(i)*cos(w*s(i))-r(i)*w**2*sin(w*s(i))
     ddz(i)=0.0_dp
     
     TT(i)=dx(i)*ddx(i)+dy(i)*ddy(i)+dz(i)*ddz(i)
     
     dTds(1,i)=(ddx(i)*deno**2-dx(i)*TT(i))/deno**3
     dTds(2,i)=(ddy(i)*deno**2-dy(i)*TT(i))/deno**3
     dTds(3,i)=(ddz(i)*deno**2-dz(i)*TT(i))/deno**3

     denoD3=sqrt(dTds(1,i)**2+dTds(2,i)**2+dTds(3,i)**2)

     D1(1,i)=dTds(1,i)/denoD3
     D1(2,i)=dTds(2,i)/denoD3
     D1(3,i)=dTds(3,i)/denoD3

     D2(1,i)=D3(2,i)*D1(3,i)-D3(3,i)*D1(2,i)
     D2(2,i)=D3(3,i)*D1(1,i)-D3(1,i)*D1(3,i)
     D2(3,i)=D3(1,i)*D1(2,i)-D3(2,i)*D1(1,i)

     dddr(i)=2.0_dp*r0*k*(4.0_dp*k**2*s(i)**3-6.0_dp*k*s(i))*exp(-k*s(i)*s(i))

     dddx(i)=dddr(i)*cos(w*s(i))-3.0_dp*w*ddr(i)*sin(w*s(i))-3.0_dp*dr(i)*w**2*cos(w*s(i))+w**3*r(i)*sin(w*s(i))
     dddy(i)=dddr(i)*sin(w*s(i))+3.0_dp*w*ddr(i)*cos(w*s(i))-3.0_dp*dr(i)*w**2*sin(w*s(i))-w**3*r(i)*cos(w*s(i))
     dddz(i)=0.0_dp

     crs(1)=dy(i)*ddz(i)-dz(i)*ddy(i)
     crs(2)=dz(i)*ddx(i)-dx(i)*ddz(i);
     crs(3)=dx(i)*ddy(i)-dy(i)*ddx(i);
     KAPPA(i)=sqrt(crs(1)**2+crs(2)**2+crs(3)**2)/deno**3;
     TAU(i)=(crs(1)*dddx(i)+crs(2)*dddy(i)+crs(3)*dddz(i))/(crs(1)**2+crs(2)**2+crs(3)**2);

  end do
!!!print*,D1(:,1),D1(:,2),D2(:,1),D2(:,2),D3(:,1),D3(:,2), Kappa(1), Kappa(2),Tau(1),Tau(2)
!!!dr=KAPPA; ddr=TAU; do i=nmot,nPt-1; KAPPA(i)=(dr(i-nmot+1)+dr(i-nmot+2))/2.0_dp;
!!!                                    TAU(i)=(ddr(i-nmot+1)+ddr(i-nmot+2))/2.0_dp; enddo
  KAPPA(1:nmot+nhook-1)=0.0_dp; TAU(1:nmot+nhook-1)=0.0_dp
  dr=KAPPA; ddr=TAU;
  do i=1,nPt-1;
     KAPPA(i)=(dr(i)+dr(i+1))/2.0_dp;
     TAU(i)=(ddr(i)+ddr(i+1))/2.0_dp;
  end do
  D1(1,1:nmot+nhook-1)=1.0_dp
  D1(2,1:nmot+nhook-1)=0.0_dp
  D1(3,1:nmot+nhook-1)=0.0_dp
  D2(1,1:nmot+nhook-1)=0.0_dp
  D2(2,1:nmot+nhook-1)=1.0_dp
  D2(3,1:nmot+nhook-1)=0.0_dp
  D3(1,1:nmot+nhook-1)=0.0_dp
  D3(2,1:nmot+nhook-1)=0.0_dp
  D3(3,1:nmot+nhook-1)=1.0_dp

  kappa2=kappa
  tau2=tau

!!!----------------------------------
end subroutine initial_higdon
