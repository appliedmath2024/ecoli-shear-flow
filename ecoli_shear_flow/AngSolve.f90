!!**************************************************************
!!AngSolve - Solving for Angular Velocity with
!!           Rotlet like term from the force (force_pt) and
!!           Dipole like term from the moment/torque (moment_pt)
!**************************************************************

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

  do i=1,nrod
     uARot=0.0_dp; vARot=0.0_dp; wARot=0.0_dp
     uADi=0.0_dp; vADi=0.0_dp; wADi=0.0_dp
     xc=data_Pt(1,:,i); yc=data_Pt(2,:,i); zc=data_Pt(3,:,i)

     !$OMP PARALLEL
     !$OMP DO PRIVATE(j,jPt,dx,dy,dz,r,rXdel,Q1,Qinv,ndotX,Q2,Q3) REDUCTION(+:uARot,vARot,wARot,uADi,vADi,wADi)
     do jPt=1,nPt;
        do j=1,nrod;
!!!Rotlet due to the force
!!!blob is: (15*del^4)/(8*pi*(r^2+del^2)^7/2)
!!! Q1=(1/2)*G'/r
           dx=(xc-data_Pt(1,jPt,j))
           dy=(yc-data_Pt(2,jPt,j))	
           dz=(zc-data_Pt(3,jPt,j)) 				
           r=dsqrt(dx*dx+dy*dy+dz*dz); rXdel=(r*r+del*del)
           Q1=(3.0_dp*del*del+2.0_dp*rXdel)/(16.0_dp*pi*rXdel*rXdel*sqrt(rXdel))
           uARot=uARot+Q1*(force_Pt(2,jPt,j)*dz-force_Pt(3,jPt,j)*dy)
           vARot=vARot+Q1*(force_Pt(3,jPt,j)*dx-force_Pt(1,jPt,j)*dz)
           wARot=wARot+Q1*(force_Pt(1,jPt,j)*dy-force_Pt(2,jPt,j)*dx)

           Qinv=1.0_dp/(32.0_dp*pi*rXdel*rXdel*rXdel*sqrt(rXdel))
           ndotX=dx*moment_Pt(1,jPt,j)+dy*moment_Pt(2,jPt,j)+dz*moment_Pt(3,jPt,j)
           Q2=(15.0_dp*del*del+6.0_dp*rXdel)*Qinv
           Q3=(15.0_dp*del**(4.0_dp)-3.0_dp*del*del*rXdel-2.0_dp*rXdel*rXdel)*Qinv
           uADi=uADi+Q2*dx*ndotX+Q3*moment_Pt(1,jPt,j)
           vADi=vADi+Q2*dy*ndotX+Q3*moment_Pt(2,jPt,j)
           wADi=wADi+Q2*dz*ndotX+Q3*moment_Pt(3,jPt,j)
        end do;
     end do
     !$OMP END DO

     !$OMP DO PRIVATE(jPt,dx,dy,dz,r,rXdel,Q1) REDUCTION(+:uARot,vARot,wARot)
     do jPt=1,size(LX,2)
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
     !$OMP END PARALLEL

     data_ang(1,:,i)=(uADi+uARot)/mu
     data_ang(2,:,i)=(vADi+vARot)/mu
     data_ang(3,:,i)=(wADi+wARot)/mu
  end do

END SUBROUTINE AngSolve
