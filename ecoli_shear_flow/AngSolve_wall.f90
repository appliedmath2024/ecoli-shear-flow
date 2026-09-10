!!**************************************************************
!!AngSolve - Solving for Angular Velocity with
!!           Rotlet like term from the force (force_pt) and
!!           Dipole like term from the moment/torque (moment_pt)
!**************************************************************

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
!-----------------------------------------
!!! Make Image Points!!!!!!!!!!!!!!!!!
image_pt(1,:,:)=data_pt(1,:,:); image_pt(2,:,:)=-data_pt(2,:,:); image_pt(3,:,:)=data_pt(3,:,:);
LX_im(1,:)=LX(1,:); LX_im(2,:)=-LX(2,:); LX_im(3,:)=LX(3,:);

do k=1,nrod

uARot=0.0_dp; vARot=0.0_dp; wARot=0.0_dp; uADi=0.0_dp; vADi=0.0_dp; wADi=0.0_dp
xc=data_Pt(1,:,k); yc=data_Pt(2,:,k); zc=data_Pt(3,:,k)

!$OMP PARALLEL

!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,Q1,rXdel) REDUCTION(+:uARot,vARot,wARot)
do jPt=1,nPt; do i=1,nrod
  !!Rotlet due to the force
  !!blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
  !! Q1=(1/2)*G'/r
   dx=(xc-data_Pt(1,jPt,i))
   dy=(yc-data_Pt(2,jPt,i))	
   dz=(zc-data_Pt(3,jPt,i)) 	
   r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
   Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))			
   uARot=uARot+Q1*(force_Pt(2,jPt,i)*dz-force_Pt(3,jPt,i)*dy)
   vARot=vARot+Q1*(force_Pt(3,jPt,i)*dx-force_Pt(1,jPt,i)*dz)
   wARot=wARot+Q1*(force_Pt(1,jPt,i)*dy-force_Pt(2,jPt,i)*dx)
end do; end do
!$OMP END DO

!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,Q1,rXdel,hh,deno,deno2,H4,rD1,rD2,gdotX) REDUCTION(+:uARot,vARot,wARot)
do jPt=1,nPt; do i=1,nrod
  !!Rotlet due to the force
  !!blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
  !! Q1=(1/2)*G'/r
    dx=(xc-image_Pt(1,jPt,i))
    dy=(yc-image_Pt(2,jPt,i))
    dz=(zc-image_Pt(3,jPt,i))
    hh=abs(image_Pt(2,jPt,i))
    r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
    deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
    Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))
    H4=15.0_dp*del*del/(deno*rXdel*rXdel)
    deno2=1.0_dp/(0.5_dp*deno*(rXdel))
    rD2=-3.0_dp*deno2
    rD1=(rXdel-(3.0_dp*del*del))*deno2
    gdotX=(-force_pt(1,jPt,i)*dx+force_pt(2,jPt,i)*dy-force_pt(3,jPt,i)*dz)
			
    uARot=uARot-Q1*(force_Pt(2,jPt,i)*dz-force_Pt(3,jPt,i)*dy) &
	  + hh*hh*H4*(force_Pt(2,jPt,i)*dz+force_Pt(3,jPt,i)*dy) &
	  + hh*( -force_Pt(3,jPt,i)*(r*r*H4+del*del*rD2) - H4*(-dx*dx*force_Pt(3,jPt,i)+dx*dz*force_Pt(1,jPt,i)) ) &
	  - 2.0_dp*hh*Q1*force_Pt(3,jPt,i) + hh*gdotX*dz*(rD2-H4)

    vARot=vARot-Q1*(force_Pt(3,jPt,i)*dx-force_Pt(1,jPt,i)*dz) &
	  + hh*hh*H4*(-force_Pt(3,jPt,i)*dx+force_Pt(1,jPt,i)*dz) &
	  - hh*H4*(-dx*dy*force_Pt(3,jPt,i)+dy*dz*force_Pt(1,jPt,i))
	
    wARot=wARot-Q1*(force_Pt(1,jPt,i)*dy-force_Pt(2,jPt,i)*dx) &
	  + hh*hh*H4*(-force_Pt(1,jPt,i)*dy-force_Pt(2,jPt,i)*dx) &
	  + hh*( force_Pt(1,jPt,i)*(r*r*H4+del*del*rD2) - H4*(-dx*dz*force_Pt(3,jPt,i)+dz*dz*force_Pt(1,jPt,i)) ) &
	  + 2.0_dp*hh*Q1*force_Pt(1,jPt,i) - hh*gdotX*dx*(rD2-H4)
end do; end do
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$OMP DO PRIVATE(jPt,dx,dy,dz,r,Q1,rXdel) REDUCTION(+:uARot,vARot,wARot)
do jPt=1,Kvmax;
    dx=(xc-LX(1,jPt))
    dy=(yc-LX(2,jPt))
    dz=(zc-LX(3,jPt))
    r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
    Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))
    uARot=uARot+Q1*(LF(2,jPt)*dz-LF(3,jPt)*dy)
    vARot=vARot+Q1*(LF(3,jPt)*dx-LF(1,jPt)*dz)
    wARot=wARot+Q1*(LF(1,jPt)*dy-LF(2,jPt)*dx)
end do
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!! image system for body point
!$OMP DO PRIVATE(jPt,dx,dy,dz,r,Q1,rXdel,hh,deno,deno2,H4,rD1,rD2,gdotX) REDUCTION(+:uARot,vARot,wARot)
do jPt=1,Kvmax;
  !!Rotlet due to the force
  !!blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
  !! Q1=(1/2)*G'/r
    dx=(xc-LX_im(1,jPt))
    dy=(yc-LX_im(2,jPt))
    dz=(zc-LX_im(3,jPt))
    hh=abs(LX_im(2,jPt))
    r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
    deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
    Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))
    H4=15.0_dp*del*del/(deno*rXdel*rXdel)
    deno2=1.0_dp/(0.5_dp*deno*(rXdel))
    rD2=-3.0_dp*deno2
    rD1=(rXdel-(3.0_dp*del*del))*deno2
    gdotX=(-LF(1,jPt)*dx+LF(2,jPt)*dy-LF(3,jPt)*dz)

    uARot= uARot-Q1*(LF(2,jPt)*dz-LF(3,jPt)*dy) &
	     + hh*hh*H4*(LF(2,jPt)*dz+LF(3,jPt)*dy) &
	     + hh*dx*H4*(LF(3,jPt)*dx-LF(1,jPt)*dz)  - hh*LF(3,jPt)*(r*r*H4 + del*del*rD2)   &
	     - 2.0_dp*hh*Q1*LF(3,jPt) + hh*gdotX*dz*(rD2-H4)

    vARot= vARot-Q1*(LF(3,jPt)*dx-LF(1,jPt)*dz)  &
         + hh*hh*H4*(-LF(3,jPt)*dx+LF(1,jPt)*dz) &
	     - hh*dy*H4*(-LF(3,jPt)*dx+LF(1,jPt)*dz)
	
    wARot= wARot-Q1*( LF(1,jPt)*dy-LF(2,jPt)*dx) &
         + hh*hh*H4*(-LF(1,jPt)*dy-LF(2,jPt)*dx) &
	     + hh*dz*H4*( LF(3,jPt)*dx-LF(1,jPt)*dz) + hh*LF(1,jPt)*(r*r*H4 + del*del*rD2)   &
    	 + 2.0_dp*hh*Q1*LF(1,jPt) - hh*gdotX*dx*(rD2-H4)
end do
!$OMP END DO Nowait

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,ndotX,Q2,Q3,rXdel,Qinv) REDUCTION(+:uADi,vADi,wADi)
do jPt=1,nPt; do i=1,nrod
 !!Dipole due to the moment
 !!Blob is: (15*del^4)/(8*pi*(del^2+r^2)^(7/2))
 !!Q2=(-1/4)*(G''/r^2-G'/r^3)
 !!Q3=(-1/4)*(G'/r-blob)
   dx=(xc-data_Pt(1,jPt,i))
   dy=(yc-data_Pt(2,jPt,i))
   dz=(zc-data_Pt(3,jPt,i))
   r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
   Qinv=1.0_dp/(32.0_dp*pi*rXdel*rXdel*rXdel*sqrt(rXdel))
   ndotX=dx*moment_Pt(1,jPt,i)+dy*moment_Pt(2,jPt,i)+dz*moment_Pt(3,jPt,i)
   Q2=(15.0_dp*del*del+6.0_dp*rXdel)*Qinv
   Q3=(15.0_dp*del**(4.0_dp)-3.0_dp*del*del*rXdel-2.0_dp*rXdel*rXdel)*Qinv
   uADi=uADi+Q2*dx*ndotX+Q3*moment_Pt(1,jPt,i)
   vADi=vADi+Q2*dy*ndotX+Q3*moment_Pt(2,jPt,i)
   wADi=wADi+Q2*dz*ndotX+Q3*moment_Pt(3,jPt,i)
end do; end do
!$OMP END DO

!!!!!!!!!!!!!!!!!!!!!! image system of Dipole due to the moment
!$OMP DO PRIVATE(jPt,i,dx,dy,dz,r,ndotX,Q1,Q2,Q3,rXdel,Qinv,deno,deno2,H4,rD2,hh,H5,qdotX,pdotX) REDUCTION(+:uADi,vADi,wADi)
do jPt=1,nPt; do i=1,nrod
 !!Dipole due to the moment
 !!Blob is: (15*del^4)/(8*pi*(del^2+r^2)^(7/2))
 !!Q2=(-1/4)*(G''/r^2-G'/r^3)
 !!Q3=(-1/4)*(G'/r-blob)
   dx=(xc-image_Pt(1,jPt,i))
   dy=(yc-image_Pt(2,jPt,i))
   dz=(zc-image_Pt(3,jPt,i))
   hh=abs(image_Pt(2,jPt,i))

   r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
   Qinv=1.0_dp/(32.0_dp*pi*rXdel*rXdel*rXdel*sqrt(rXdel))
   ndotX=dx*moment_Pt(1,jPt,i)+dy*moment_Pt(2,jPt,i)+dz*moment_Pt(3,jPt,i)
   Q2=(15.0_dp*del*del+6.0_dp*rXdel)*Qinv
   Q3=(15.0_dp*del**(4.0_dp)-3.0_dp*del*del*rXdel-2.0_dp*rXdel*rXdel)*Qinv
   deno=8.0_dp*pi*(rXdel)*sqrt(rXdel)
   Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))
   H4=15.0_dp*del*del/(deno*rXdel*rXdel)
   deno2=1.0_dp/(0.5_dp*deno*(rXdel))
   rD2=-3.0_dp*deno2
   H5=105.0_dp*del*del/(deno*rXdel*rXdel*rXdel)
   qdotX=dx*moment_pt(1,jPt,i) + dz*moment_pt(3,jPt,i)
   pdotX=-dx*moment_Pt(3,jPt,i)+dz*moment_Pt(1,jPt,i)

   uADi=uADi - Q2*dx*ndotX - Q3*moment_Pt(1,jPt,i) &
		 + 0.5_dp*(H4-rD2)*( dz*pdotX-dy*dy*moment_pt(1,jPt,i)) &
		 + hh*H4*dy*moment_pt(1,jPt,i) - 0.5_dp*(H4*r*r+rD2*del*del)*moment_pt(1,jPt,i) &
		 + 0.5_dp*H4*dx*qdotX  &
		 - 0.5_dp*hh*( H4*(4.0_dp*dy*moment_pt(1,jPt,i) - dx*moment_pt(2,jPt,i))  &
                     - H5*dy*(r*r*moment_pt(1,jPt,i) - ndotX*dx) ) &
		 + 0.5_dp*hh*hh*( (2.0_dp*H4 - H5*r*r)*moment_pt(1,jPt,i) + H5*ndotX*dx )


   vADi=vADi - Q2*dy*ndotX  - Q3*moment_Pt(2,jPt,i) &
    	 + 0.5_dp*(H4-rD2)*dy*qdotX &
	     + H4*qdotX *(-hh + 0.5_dp*dy)  &
		 - 0.5_dp*hh*( H4*(2.0_dp*dy*moment_pt(2,jPt,i) - qdotX) &
                     - H5*dy*(r*r*moment_pt(2,jPt,i) - ndotX*dy) ) &
        + 0.5_dp*hh*hh*( (2.0_dp*H4 - H5*r*r)*moment_pt(2,jPt,i) + H5*ndotX*dy )


   wADi=wADi - Q2*dz*ndotX - Q3*moment_Pt(3,jPt,i) &
		 + 0.5_dp*(H4-rD2)*(-dx*pdotX-dy*dy*moment_pt(3,jPt,i)) &
		 + hh*H4*dy*moment_pt(3,jPt,i) - 0.5_dp*(H4*r*r+rD2*del*del)*moment_pt(3,jPt,i) &
		 + 0.5_dp*H4*dz*qdotx  &
		 - 0.5_dp*hh*( H4*(4.0_dp*dy*moment_pt(3,jPt,i) - dz*moment_pt(2,jPt,i))  &
                     - H5*dy*(r*r*moment_pt(3,jPt,i) - ndotX*dz) ) &
		 + 0.5_dp*hh*hh*( (2.0_dp*H4 - H5*r*r)*moment_pt(3,jPt,i) + H5*ndotX*dz )

end do; end do
!$OMP END DO

!$OMP END PARALLEL

data_ang(1,:,k)=(uADi+uARot)/mu; data_ang(2,:,k)=(vADi+vARot)/mu
data_ang(3,:,k)=(wADi+wARot)/mu

END DO

END SUBROUTINE
