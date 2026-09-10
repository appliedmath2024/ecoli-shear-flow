!**************************************************************
!!StSolve - Solving Stoke's Eqn with Stokeslets and Rotlets
!**************************************************************

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

  do i=1,nrod;

!!!Zeroing out u, v, w components of velocity from forces due to Stokeslets
     uSt=0.0_dp; vSt=0.0_dp; wSt=0.0_dp
!!!Zeroing out u, v, w components of velocity from torques due to Rotlets
     uRot=0.0_dp; vRot=0.0_dp; wRot=0.0_dp
     xc=data_Pt(1,:,i); yc=data_Pt(2,:,i); zc=data_Pt(3,:,i)

     !$OMP PARALLEL
     !$OMP DO PRIVATE(j,jPt,dx,dy,dz,r,rXdel,deno,H2,H1,fdotX,H3) REDUCTION(+:uSt,vSt,wSt,uRot,vRot,wRot)
     do jPt=1,nPt;
        do j=1,nrod
!!!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
           dx=(xc-data_Pt(1,jPt,j))
           dy=(yc-data_Pt(2,jPt,j))
           dz=(zc-data_Pt(3,jPt,j))
           r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=r*r+del*del
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
        end do;
     end do
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
     !$OMP END PARALLEL
     data_vel(1,:,i)=(uSt+uRot)/mu
     data_vel(2,:,i)=(vSt+vRot)/mu
     data_vel(3,:,i)=(wSt+wRot)/mu
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LU=0.0_dp; LU1=0.0_dp; LU2=0.0_dp; LU3=0.0_dp
  !$OMP PARALLEL
  !$OMP DO PRIVATE(jPt,ldx,ldy,ldz,lr,lrXdel,ldeno,lH2,lH1,lfdotX) REDUCTION(+:LU1,LU2,LU3)
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

  !$OMP DO PRIVATE(j,jPt,ldx,ldy,ldz,lr,lrXdel,ldeno,lH2,lH1,lfdotX,lH3) REDUCTION(+:LU1,LU2,LU3)
  do jPt=1,nPt;
     do j=1,nrod
!!!Stokeslet/force blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
        ldx=(LX(1,:)-data_Pt(1,jPt,j))
        ldy=(LX(2,:)-data_Pt(2,jPt,j))
        ldz=(LX(3,:)-data_Pt(3,jPt,j))
        lr=dsqrt(ldx*ldx+ldy*ldy+ldz*ldz);  lrXdel=lr*lr+del*del
        ldeno=8.0_dp*pi*(lrXdel)*sqrt(lrXdel)
        lH2=1.0_dp/ldeno
        lH1=(lrXdel+del*del)*lH2
        lfdotX=(force_Pt(1,jPt,j)*ldx+force_Pt(2,jPt,j)*ldy+force_Pt(3,jPt,j)*ldz)
        LU1=LU1+force_Pt(1,jPt,j)*lH1+lfdotX*ldx*lH2
        LU2=LU2+force_Pt(2,jPt,j)*lH1+lfdotX*ldy*lH2
        LU3=LU3+force_Pt(3,jPt,j)*lH1+lfdotX*ldz*lH2
!!!Rotlet/moment blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
!!! H3=(1/2)*G'/r
        lH3=(3.0_dp*del*del+2.0_dp*lrXdel)/(2.0_dp*ldeno*lrXdel)
        LU1=LU1+lH3*(moment_Pt(2,jPt,j)*ldz-moment_Pt(3,jPt,j)*ldy)
        LU2=LU2+lH3*(moment_Pt(3,jPt,j)*ldx-moment_Pt(1,jPt,j)*ldz)
        LU3=LU3+lH3*(moment_Pt(1,jPt,j)*ldy-moment_Pt(2,jPt,j)*ldx)
     end do;
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  LU(1,:)=LU1; LU(2,:)=LU2; LU(3,:)=LU3; LU=LU/mu

END SUBROUTINE StSolve
