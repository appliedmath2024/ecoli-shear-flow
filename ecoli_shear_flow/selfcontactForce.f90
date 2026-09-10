SUBROUTINE selfcontactForce(coord_pt,force_pt)
  USE omp_lib
  USE nrtype; USE commonval
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:,:), INTENT(in) :: coord_pt
  REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: force_pt

!!! local variables
  INTEGER(I4B) :: i, j, k, jPt
  REAL(DP), DIMENSION(nPt) :: xc, yc, zc, dx, dy, dz, r, ten, fx, fy, fz, dist
  REAL(DP) :: max_dist

  do i = 1, nrod-1

     xc=coord_pt(1,:,i); yc=coord_pt(2,:,i); zc=coord_pt(3,:,i)
     fx=0.0_dp; fy=0.0_dp; fz=0.0_dp

     !$OMP PARALLEL
     !$OMP DO PRIVATE(j,jPt,k,dx,dy,dz,r,dist,ten,max_dist) REDUCTION(+:fx,fy,fz)
     do jPt=1,nPt;
        do j=i+1,nrod;
           dx=(xc-coord_pt(1,jPt,j))
           dy=(yc-coord_pt(2,jPt,j))
           dz=(zc-coord_pt(3,jPt,j))

           dist=dsqrt(dx*dx+dy*dy+dz*dz);
           r=1.0_dp/dist-inv_Dmin;
           max_dist=maxval(r)
           if (max_dist > 0.0_dp) then
             ten=0.0_dp;
             do k=1,nPt;
                if (r(k)>0.0_dp) then;
                   ten(k)=stiff_contact * r(k);
!!!               ten(k)=stiff_contact * (tan(r(k))-0.0_dp);
                endif;
             end do
             dx = ten*dx; dy = ten*dy; dz = ten*dz
             fx = fx + dx
             fy = fy + dy
             fz = fz + dz
             force_pt(1,jPt,j)=force_pt(1,jPt,j) -  sum( dx )
             force_pt(2,jPt,j)=force_pt(2,jPt,j) -  sum( dy )
             force_pt(3,jPt,j)=force_pt(3,jPt,j) -  sum( dz )
           end if
        end do;
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     force_pt(1,:,i)=force_pt(1,:,i)+fx;
     force_pt(2,:,i)=force_pt(2,:,i)+fy;
     force_pt(3,:,i)=force_pt(3,:,i)+fz;
  end do

end subroutine selfcontactForce
