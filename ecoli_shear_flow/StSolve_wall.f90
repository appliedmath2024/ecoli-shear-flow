!**************************************************************
!!StSolve - Solving Stoke's Eqn with Stokeslets and Rotlets
!**************************************************************

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

!-----------------------------------------
!!! Make Image Points!!!!!!!!!!!!!!!!!
image_pt(1,:,:)=data_pt(1,:,:); image_pt(2,:,:)=-data_pt(2,:,:); image_pt(3,:,:)=data_pt(3,:,:);
LX_im(1,:)=LX(1,:); LX_im(2,:)=-LX(2,:); LX_im(3,:)=LX(3,:);

do k=1,nrod

! from forces due to Stokeslets
uSt=0.0_dp; vSt=0.0_dp; wSt=0.0_dp
uRot=0.0_dp; vRot=0.0_dp; wRot=0.0_dp
LU1=0.0_dp; LU2=0.0_dp; LU3=0.0_dp

!$OMP PARALLEL

!!!!! rod <-- rod
!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,deno,H1,H2,H3,fdotX,rXdel) REDUCTION(+:uSt,vSt,wSt,uRot,vRot,wRot)
do jPt=1,nPt; do i=1,nrod
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    dx=(data_Pt(1,:,k)-data_Pt(1,jPt,i))
    dy=(data_Pt(2,:,k)-data_Pt(2,jPt,i))
    dz=(data_Pt(3,:,k)-data_Pt(3,jPt,i))
    r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
    deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
    H2=1.0_dp/deno
    H1=(rXdel+del*del)*H2
    fdotX=(force_Pt(1,jPt,i)*dx+force_Pt(2,jPt,i)*dy+force_Pt(3,jPt,i)*dz)
    uSt=uSt+force_Pt(1,jPt,i)*H1+fdotX*dx*H2
    vSt=vSt+force_Pt(2,jPt,i)*H1+fdotX*dy*H2
    wSt=wSt+force_Pt(3,jPt,i)*H1+fdotX*dz*H2
 !!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
 !! H3=(1/2)*G'/r
    H3=(3.0_dp*del*del+2.0_dp*rXdel)/(2.0_dp*rXdel*deno)
    uRot=uRot+H3*(moment_Pt(2,jPt,i)*dz-moment_Pt(3,jPt,i)*dy)
    vRot=vRot+H3*(moment_Pt(3,jPt,i)*dx-moment_Pt(1,jPt,i)*dz)
    wRot=wRot+H3*(moment_Pt(1,jPt,i)*dy-moment_Pt(2,jPt,i)*dx)
end do; end do
!$OMP END DO

!!!!! image system for rod <-- rod
!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,deno,H1,H2,H3,H4,deno2,rD1,rD2,fdotX,gdotX,rXdel,hh,pdotX,LxX1,LxX2,LxX3) REDUCTION(+:uSt,vSt,wSt,uRot,vRot,wRot)
do jPt=1,nPt; do i=1,nrod
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    dx=(data_Pt(1,:,k)-image_Pt(1,jPt,i))
    dy=(data_Pt(2,:,k)-image_Pt(2,jPt,i))
    dz=(data_Pt(3,:,k)-image_Pt(3,jPt,i))
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

!!!! rod <-- body
!$OMP DO PRIVATE(jPt,dx,dy,dz,r,deno,H1,H2,fdotX,rXdel) REDUCTION(+:uSt,vSt,wSt)
do jPt=1,Kvmax;
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    dx=(data_Pt(1,:,k)-LX(1,jPt))
    dy=(data_Pt(2,:,k)-LX(2,jPt))
    dz=(data_Pt(3,:,k)-LX(3,jPt))
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
    dx=(data_Pt(1,:,k)-LX_im(1,jPt))
    dy=(data_Pt(2,:,k)-LX_im(2,jPt))
    dz=(data_Pt(3,:,k)-LX_im(3,jPt))
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

data_vel(1,:,k)=(uSt+uRot)/mu; data_vel(2,:,k)=(vSt+vRot)/mu; data_vel(3,:,k)=(wSt+wRot)/mu
END DO

!$OMP PARALLEL
!!!! CELL BODY <-- rod
!$OMP DO PRIVATE(jPt,i,ldx,ldy,ldz,lr,ldeno,lH1,lH2,lH3,lfdotX,lrXdel) REDUCTION(+:LU1,LU2,LU3)
do jPt=1,nPt; do i=1,nrod
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    ldx=(LX(1,:)-data_Pt(1,jPt,i))
    ldy=(LX(2,:)-data_Pt(2,jPt,i))
    ldz=(LX(3,:)-data_Pt(3,jPt,i))
    lr=dsqrt(ldx*ldx+ldy*ldy+ldz*ldz);  lrXdel=lr*lr+del*del
    ldeno=8.0_dp*pi*(lrXdel)*sqrt(lrXdel)
    lH2=1.0_dp/ldeno
    lH1=(lrXdel+del*del)*lH2
    lfdotX=(force_Pt(1,jPt,i)*ldx+force_Pt(2,jPt,i)*ldy+force_Pt(3,jPt,i)*ldz)
    LU1=LU1+force_Pt(1,jPt,i)*lH1+lfdotX*ldx*lH2
    LU2=LU2+force_Pt(2,jPt,i)*lH1+lfdotX*ldy*lH2
    LU3=LU3+force_Pt(3,jPt,i)*lH1+lfdotX*ldz*lH2
 !!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
 !! H3=(1/2)*G'/r
    lH3=(3.0_dp*del*del+2.0_dp*lrXdel)/(2.0_dp*ldeno*lrXdel)
    LU1=LU1+lH3*(moment_Pt(2,jPt,i)*ldz-moment_Pt(3,jPt,i)*ldy)
    LU2=LU2+lH3*(moment_Pt(3,jPt,i)*ldx-moment_Pt(1,jPt,i)*ldz)
    LU3=LU3+lH3*(moment_Pt(1,jPt,i)*ldy-moment_Pt(2,jPt,i)*ldx)
end do; end do
!$OMP END DO

!!!! image system for body <-- rod
!$OMP DO PRIVATE(jPt,i,ldx,ldy,ldz,lr,ldeno,lH1,lH2,lH3,lH4,ldeno2,lD1,lD2,lfdotX,lgdotX,lrXdel,lhh,lpdotX,lLxX1,lLxX2,lLxX3) REDUCTION(+:LU1,LU2,LU3)
do jPt=1,nPt; do i=1,nrod
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    ldx=(LX(1,:)-image_Pt(1,jPt,i))
    ldy=(LX(2,:)-image_Pt(2,jPt,i))
    ldz=(LX(3,:)-image_Pt(3,jPt,i))
    lhh=abs(image_Pt(2,jPt,i))
    lr=dsqrt(ldx*ldx+ldy*ldy+ldz*ldz); lrXdel=lr*lr+del*del
    ldeno=8.0_dp*pi*(lrXdel)*sqrt(lrXdel)
    lH2=1.0_dp/ldeno
    lH1=(lrXdel+del*del)*lH2
    lfdotX=(force_Pt(1,jPt,i)*ldx+force_Pt(2,jPt,i)*ldy+force_Pt(3,jPt,i)*ldz)
    ldeno2=1.0_dp/(0.5_dp*ldeno*(lrXdel))
    lD2=-3.0_dp*ldeno2
    lD1=(lrXdel-(3.0_dp*del*del))*ldeno2
    lgdotX=(-force_pt(1,jPt,i)*ldx+force_pt(2,jPt,i)*ldy-force_pt(3,jPt,i)*ldz)

    LU1=LU1-force_Pt(1,jPt,i)*lH1-lfdotX*ldx*lH2-(lhh*lhh*(-force_pt(1,jPt,i)*lD1+lgdotX*ldx*lD2)) &
        +(2.0_dp*lhh*((ldx*force_pt(2,jPt,i)-ldy*force_pt(1,jPt,i))*lH2+ldy*ldx*lgdotX*0.5_dp*lD2)) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(-ldy*force_pt(1,jPt,i))))
    LU2=LU2-force_Pt(2,jPt,i)*lH1-lfdotX*ldy*lH2-(lhh*lhh*(force_pt(2,jPt,i)*lD1+lgdotX*ldy*lD2)) &
        +(2.0_dp*lhh*((2.0_dp*ldy*force_pt(2,jPt,i))*lH2+lgdotX*(-0.5_dp*lD1+(del*del*lD2))+ldy*ldy*lgdotX*0.5_dp*lD2)) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(ldx*force_pt(1,jPt,i)+ldz*force_pt(3,jPt,i))))
    LU3=LU3-force_Pt(3,jPt,i)*lH1-lfdotX*ldz*lH2-(lhh*lhh*(-force_pt(3,jPt,i)*lD1+lgdotX*ldz*lD2)) &
        +(2.0_dp*lhh*((ldz*force_pt(2,jPt,i)-ldy*force_pt(3,jPt,i))*lH2+ldy*ldz*lgdotX*0.5_dp*lD2)) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(-ldy*force_pt(3,jPt,i))))
!!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
 !! H3=(1/2)*G'/r
    lH3=(3.0_dp*del*del+2.0_dp*lrXdel)/(2.0_dp*lrXdel*ldeno)
    lLxX1 =ldz*moment_Pt(2,jPt,i)-ldy*moment_Pt(3,jPt,i)
    lLxX2 =ldx*moment_Pt(3,jPt,i)-ldz*moment_Pt(1,jPt,i)
    lLxX3 =ldy*moment_Pt(1,jPt,i)-ldx*moment_Pt(2,jPt,i)
    lpdotX=-lLxX2
    lH4=15.0_dp*del*del/(ldeno*lrXdel*lrXdel)
    LU1=LU1-lH3*(moment_Pt(2,jPt,i)*ldz-moment_Pt(3,jPt,i)*ldy) &
        + ldy*moment_Pt(3,jPt,i)*(del*del*lD2)    - ldx*ldy*lpdotX*lD2  &
        - lhh*moment_Pt(3,jPt,i)*(0.5_dp*lD1+lH2) + ldx*lhh*lpdotX*lD2 &
        - lhh*ldy*lLxX1*lH4 + lhh*lhh*lLxX1*lH4
    LU2=LU2-lH3*(moment_Pt(3,jPt,i)*ldx-moment_Pt(1,jPt,i)*ldz) &
        - lpdotX*(0.5_dp*del*del*lD2)  - ldy*ldy*lpdotX*lD2 &
        + lpdotX*(0.5_dp*del*del*lD2)  + lhh*ldy*lpdotX*lD2 &
        - lhh*ldy*lLxX2*lH4 + lhh*lhh*lLxX2*lH4
    LU3=LU3-lH3*(moment_Pt(1,jPt,i)*ldy-moment_Pt(2,jPt,i)*ldx) &
        - ldy*moment_Pt(1,jPt,i)*(del*del*lD2)    - ldz*ldy*lpdotX*lD2  &
        + lhh*moment_Pt(1,jPt,i)*(0.5_dp*lD1+lH2)  + ldz*lhh*lpdotX*lD2 &
        - lhh*ldy*lLxX3*lH4 + lhh*lhh*lLxX3*lH4
end do; end do
!$OMP END DO

!!!!! body <-- body
!$OMP DO PRIVATE(jPt,i,ldx,ldy,ldz,lr,ldeno,lH1,lH2,lfdotX,lrXdel) REDUCTION(+:LU1,LU2,LU3)
do jPt=1,Kvmax
 !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    ldx=(LX(1,:)-LX(1,jPt))
    ldy=(LX(2,:)-LX(2,jPt))
    ldz=(LX(3,:)-LX(3,jPt))
    lr=dsqrt(ldx*ldx+ldy*ldy+ldz*ldz); lrXdel=lr*lr+del*del
    ldeno=8.0_dp*pi*(lrXdel)*sqrt(lrXdel)
    lH2=1.0_dp/ldeno
    lH1=(lrXdel+del*del)*lH2
    lfdotX=(LF(1,jPt)*ldx+LF(2,jPt)*ldy+LF(3,jPt)*ldz)
    LU1=LU1+LF(1,jPt)*lH1+lfdotX*ldx*lH2
    LU2=LU2+LF(2,jPt)*lH1+lfdotX*ldy*lH2
    LU3=LU3+LF(3,jPt)*lH1+lfdotX*ldz*lH2
end do
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!! image system for cell body <-- body
!$OMP DO PRIVATE(jPt,i,ldx,ldy,ldz,lr,ldeno,lH1,lH2,ldeno2,lD1,lD2,lfdotX,lgdotX,lrXdel,lhh) REDUCTION(+:LU1,LU2,LU3)
do jPt=1,Kvmax
! !!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
    ldx=(LX(1,:)-LX_im(1,jPt))
    ldy=(LX(2,:)-LX_im(2,jPt))
    ldz=(LX(3,:)-LX_im(3,jPt))
    lhh=abs(LX_im(2,jPt))
    lr=dsqrt(ldx*ldx+ldy*ldy+ldz*ldz); lrXdel=lr*lr+del*del
    ldeno=8.0_dp*pi*(lrXdel)*sqrt(lrXdel)
    lH2=1.0_dp/ldeno
    lH1=(lrXdel+del*del)*lH2
    lfdotX=(LF(1,jPt)*ldx+LF(2,jPt)*ldy+LF(3,jPt)*ldz)
    ldeno2=1.0_dp/(0.5_dp*ldeno*(lrXdel))
    lD2=-3.0_dp*ldeno2
    lD1=(lrXdel-(3.0_dp*del*del))*ldeno2
    lgdotX=(-LF(1,jPt)*ldx+LF(2,jPt)*ldy-LF(3,jPt)*ldz)

    LU1=LU1-LF(1,jPt)*lH1-lfdotX*ldx*lH2-(lhh*lhh*(-LF(1,jPt)*lD1+lgdotX*ldx*lD2)) &
        +(2.0_dp*lhh*((ldx*LF(2,jPt)-ldy*LF(1,jPt))*lH2+(ldy*ldx*lgdotX*0.5_dp*lD2))) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(-ldy*LF(1,jPt))))
    LU2=LU2-LF(2,jPt)*lH1-lfdotX*ldy*lH2-(lhh*lhh*(LF(2,jPt)*lD1+lgdotX*ldy*lD2)) &
        +(2.0_dp*lhh*((2.0_dp*ldy*LF(2,jPt))*lH2+lgdotX*(-0.5_dp*lD1+(del*del*lD2))+ldy*ldy*lgdotX*0.5_dp*lD2)) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(ldx*LF(1,jPt)+ldz*LF(3,jPt))))
    LU3=LU3-LF(3,jPt)*lH1-lfdotX*ldz*lH2-(lhh*lhh*(-LF(3,jPt)*lD1+lgdotX*ldz*lD2)) &
        +(2.0_dp*lhh*((ldz*LF(2,jPt)-ldy*LF(3,jPt))*lH2+(ldy*ldz*lgdotX*0.5_dp*lD2))) &
        +(2.0_dp*lhh*(0.5_dp*del*del*lD2*(-ldy*LF(3,jPt))))

end do
!$OMP END DO

!$OMP END PARALLEL

LU(1,:)=LU1/mu; LU(2,:)=LU2/mu; LU(3,:)=LU3/mu

END SUBROUTINE
