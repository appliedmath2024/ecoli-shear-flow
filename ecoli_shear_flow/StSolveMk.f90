!**************************************************************
!!StSolve - Solving Stoke's Eqn with Stokeslets and Rotlets     
!**************************************************************

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

  !-----------------------------------------

  N = size(marker_pt,2)

  ALLOCATE(r(N),H1(N),H2(N),H3(N),deno(N),dx(N),dy(N),dz(N),fdotX(N),rXdel(N))
  ALLOCATE(rD1(N),rD2(N),deno2(N),gdotX(N),pdotX(N),LxX1(N),LxX2(N),LxX3(N),H4(N))
  ALLOCATE(uSt(N),vSt(N),wSt(N),uRot(N),vRot(N),wRot(N))
  ALLOCATE(xc(N),yc(N),zc(N))

  !!! Make Image Points!!!!!!!!!!!!!!!!!
  image_pt(1,:,:)=data_pt(1,:,:); image_pt(2,:,:)=-data_pt(2,:,:); image_pt(3,:,:)=data_pt(3,:,:);
  LX_im(1,:)=LX(1,:); LX_im(2,:)=-LX(2,:); LX_im(3,:)=LX(3,:);

  !Zeroing out u, v, w components of velocity
  ! from forces due to Stokeslets
  uSt=0.0_dp; vSt=0.0_dp; wSt=0.0_dp

  !Zeroing out u, v, w components of velocity
  ! from torques due to Rotlets
  uRot=0.0_dp; vRot=0.0_dp; wRot=0.0_dp

  xc=marker_pt(1,:)
  yc=marker_pt(2,:)
  zc=marker_pt(3,:)

  !$OMP PARALLEL PRIVATE(dx,dy,dz,r,rxdel,deno,fdotX,H1,H2,H3)
  !$OMP DO REDUCTION (+:uSt,vSt,wSt,uRot,vRot,wRot)

  do jPt=1,nPt
     do j=1,nrod
        !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
        dx=(xc-data_Pt(1,jPt,j))
        dy=(yc-data_Pt(2,jPt,j))	
        dz=(zc-data_Pt(3,jPt,j)) 				       
        r=dsqrt(dx**2+dy**2+dz**2); rXdel=r*r+del*del
        deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
        H2=1.0_dp/deno
        H1=(rXdel+del*del)*H2
        fdotX=(force_Pt(1,jPt,j)*dx+force_Pt(2,jPt,j)*dy+force_Pt(3,jPt,j)*dz)
        uSt=uSt+force_Pt(1,jPt,j)*H1+fdotX*dx*H2
        vSt=vSt+force_Pt(2,jPt,j)*H1+fdotX*dy*H2
        wSt=wSt+force_Pt(3,jPt,j)*H1+fdotX*dz*H2
!!!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
!!! H3=(1/2)*G'/r
        H3=(3.0_dp*del*del+2.0_dp*rXdel)/(2.0_dp*rXdel*deno)
        uRot=uRot+H3*(moment_Pt(2,jPt,j)*dz-moment_Pt(3,jPt,j)*dy)
        vRot=vRot+H3*(moment_Pt(3,jPt,j)*dx-moment_Pt(1,jPt,j)*dz)
        wRot=wRot+H3*(moment_Pt(1,jPt,j)*dy-moment_Pt(2,jPt,j)*dx)
     end do
  end do
  !$OMP END DO

!!!!! image system for rod <-- rod
  !$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,deno,H1,H2,H3,H4,deno2,rD1,rD2,fdotX,gdotX,rXdel,hh,pdotX,LxX1,LxX2,LxX3) REDUCTION(+:uSt,vSt,wSt,uRot,vRot,wRot)
  do jPt=1,nPt; do i=1,nrod
     !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
     dx=(xc-image_Pt(1,jPt,i))
     dy=(yc-image_Pt(2,jPt,i))
     dz=(zc-image_Pt(3,jPt,i))
     hh=abs(image_Pt(2,jPt,i))
     r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
     deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
     H2=1.0_dp/deno
     H1=(rXdel+del*del)*H2
     fdotX=(force_Pt(1,jPt,i)*dx+force_Pt(2,jPt,i)*dy+force_Pt(3,jPt,i)*dz)
     deno2=1.0_dp/(0.5_dp*deno*(rXdel))
     rD2=-3.0_dp*deno2
     rD1=(rXdel-(3.0_dp*del*del))*deno2
     gdotX=(-force_pt(1,jPt,i)*dx+force_pt(2,jPt,i)*dy-force_pt(3,jPt,i)*dz)

     uSt=uSt-force_Pt(1,jPt,i)*H1-fdotX*dx*H2-(hh*hh*(-force_pt(1,jPt,i)*rD1+gdotX*dx*rD2)) &
          +(2.0_dp*hh*((dx*force_pt(2,jPt,i)-dy*force_pt(1,jPt,i))*H2+dy*dx*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(-dy*force_pt(1,jPt,i))))
     vSt=vSt-force_Pt(2,jPt,i)*H1-fdotX*dy*H2-(hh*hh*(force_pt(2,jPt,i)*rD1+gdotX*dy*rD2)) &
          +(2.0_dp*hh*((2.0_dp*dy*force_pt(2,jPt,i))*H2+gdotX*(-0.5_dp*rD1+(del*del*rD2))+dy*dy*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(dx*force_pt(1,jPt,i)+dz*force_pt(3,jPt,i))))
     wSt=wSt-force_Pt(3,jPt,i)*H1-fdotX*dz*H2-(hh*hh*(-force_pt(3,jPt,i)*rD1+gdotX*dz*rD2)) &
          +(2.0_dp*hh*((dz*force_pt(2,jPt,i)-dy*force_pt(3,jPt,i))*H2+dy*dz*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(-dy*force_pt(3,jPt,i))))
     ! !!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
     ! !! H3=(1/2)*G'/r
     H3=(3.0_dp*del*del+2.0_dp*rXdel)/(2.0_dp*rXdel*deno)
     LxX1 =dz*moment_Pt(2,jPt,i)-dy*moment_Pt(3,jPt,i)
     LxX2 =dx*moment_Pt(3,jPt,i)-dz*moment_Pt(1,jPt,i)
     LxX3 =dy*moment_Pt(1,jPt,i)-dx*moment_Pt(2,jPt,i)
     pdotX=-LxX2
     H4=15.0_dp*del*del/(deno*rXdel*rXdel)

     uRot=uRot-H3*LxX1 &
          + dy*moment_Pt(3,jPt,i)*(del*del*rD2)    - dx*dy*pdotX*rD2  &
          - hh*moment_Pt(3,jPt,i)*(0.5_dp*rD1+H2)  + dx*hh*pdotX*rD2 &
          + (- hh*dy + hh*hh)*LxX1*H4
     vRot=vRot-H3*LxX2 &
          - dy*dy*pdotX*rD2 &
          + hh*dy*pdotX*rD2 &
          + (- hh*dy + hh*hh)*LxX2*H4
     wRot=wRot-H3*LxX3 &
          - dy*moment_Pt(1,jPt,i)*(del*del*rD2)    - dz*dy*pdotX*rD2  &
          + hh*moment_Pt(1,jPt,i)*(0.5_dp*rD1+H2)  + dz*hh*pdotX*rD2 &
          + (- hh*dy + hh*hh)*LxX3*H4

  end do; end do
  !$OMP END DO
 
  !$OMP DO PRIVATE(jPt,dx,dy,dz,r,rXdel,deno,H2,H1,fdotX) REDUCTION(+:uSt,vSt,wSt)
  do jPt=1,Kvmax
!!!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
     dx=(xc-LX(1,jPt))
     dy=(yc-LX(2,jPt))
     dz=(zc-LX(3,jPt))
     r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
     deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
     H2=1.0_dp/deno
     H1=(rXdel+del*del)*H2
     fdotX=(LF(1,jPt)*dx+LF(2,jPt)*dy+LF(3,jPt)*dz)
     uSt=uSt+LF(1,jPt)*H1+fdotX*dx*H2
     vSt=vSt+LF(2,jPt)*H1+fdotX*dy*H2
     wSt=wSt+LF(3,jPt)*H1+fdotX*dz*H2
  end do
  !$OMP END DO

!!!! image system for rod <-- body
  !$OMP DO PRIVATE(jPt,dx,dy,dz,r,deno,H1,H2,deno2,rD1,rD2,fdotX,gdotX,rXdel,hh) REDUCTION(+:uSt,vSt,wSt)
  do jPt=1,Kvmax;
     !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
     dx=(xc-LX_im(1,jPt))
     dy=(yc-LX_im(2,jPt))
     dz=(zc-LX_im(3,jPt))
     hh=abs(LX_im(2,jPt))
     r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
     deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
     H2=1.0_dp/deno
     H1=(rXdel+del*del)*H2
     fdotX=(LF(1,jPt)*dx+LF(2,jPt)*dy+LF(3,jPt)*dz)
     deno2=1.0_dp/(0.5_dp*deno*(rXdel))
     rD2=-3.0_dp*deno2
     rD1=(rXdel-(3.0_dp*del*del))*deno2
     gdotX=(-LF(1,jPt)*dx+LF(2,jPt)*dy-LF(3,jPt)*dz)

     uSt=uSt-LF(1,jPt)*H1-fdotX*dx*H2-(hh*hh*(-LF(1,jPt)*rD1+gdotX*dx*rD2)) &
          +(2.0_dp*hh*((dx*LF(2,jPt)-dy*LF(1,jPt))*H2+dy*dx*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(-dy*LF(1,jPt))))
     vSt=vSt-LF(2,jPt)*H1-fdotX*dy*H2-(hh*hh*(LF(2,jPt)*rD1+gdotX*dy*rD2)) &
          +(2.0_dp*hh*((2.0_dp*dy*LF(2,jPt))*H2+gdotX*(-0.5_dp*rD1+(del*del*rD2))+dy*dy*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(dx*LF(1,jPt)+dz*LF(3,jPt))))
     wSt=wSt-LF(3,jPt)*H1-fdotX*dz*H2-(hh*hh*(-LF(3,jPt)*rD1+gdotX*dz*rD2)) &
          +(2.0_dp*hh*((dz*LF(2,jPt)-dy*LF(3,jPt))*H2+dy*dz*gdotX*0.5_dp*rD2)) &
          +(2.0_dp*hh*(0.5_dp*del*del*rD2*(-dy*LF(3,jPt))))

  end do
  !$OMP END DO
  !$OMP END PARALLEL

  data_vel(1,:)=(uSt+uRot)/mu
  data_vel(2,:)=(vSt+vRot)/mu
  data_vel(3,:)=(wSt+wRot)/mu

  !open(14,file="uvwSt.dat")
  !do iPt=1,nPt
  !write(14,30) uSt(ipt), vSt(ipt),wSt(ipt)
  !end do
  !close(14)

  !30 format(3(E23.16),1x)

  DEALLOCATE(r,H1,H2,H3,deno,dx,dy,dz,fdotX,rXdel)
  DEALLOCATE(rD1,rD2,deno2,gdotX,pdotX,LxX1,LxX2,LxX3,H4)
  DEALLOCATE(uSt,vSt,wSt,uRot,vRot,wRot)
  DEALLOCATE(xc,yc,zc)

END SUBROUTINe StSolveMk
