SUBROUTINE selfcontactForce_Fla_Body(coord_pt,force_pt, LX, LF)
  USE omp_lib
  USE nrtype; USE commonval
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:,:), INTENT(in) :: coord_pt
  REAL(DP), DIMENSION(:,:), INTENT(in) :: LX
  REAL(DP), DIMENSION(:,:,:), INTENT(inout) :: force_pt
  REAL(DP), DIMENSION(:,:), INTENT(inout) :: LF

!!! local variables
  INTEGER(I4B) :: i, j, k, jPt
  REAL(DP), DIMENSION(nPt) :: xc, yc, zc, dx, dy, dz, r, ten, fx, fy, fz, dist
  REAL(DP) :: max_dist

  do i = 1, nrod

     xc=coord_pt(1,:,i); yc=coord_pt(2,:,i); zc=coord_pt(3,:,i)
     fx=0.0_dp; fy=0.0_dp; fz=0.0_dp

     !$OMP PARALLEL
     !$OMP DO PRIVATE(jPt,k,dx,dy,dz,r,dist,ten,max_dist) REDUCTION(+:fx,fy,fz)
     do jPt=1,Kvmax
        dx = (xc-LX(1,jPt))
        dy = (yc-LX(2,jPt))
        dz = (zc-LX(3,jPt))

        dist = dsqrt(dx*dx+dy*dy+dz*dz);
        r(nhook+2:nPt) = 1.0_dp/dist(nhook+2:nPt) - inv_Dmin;
        max_dist = maxval(r(nhook+2:nPt))

        if (max_dist > 0.0_dp) then
           ten=0.0_dp;
           do k=nhook+2,nPt
              if (r(k)>0.0_dp) then
                 ten(k)=stiff_contact * r(k)
              endif
           enddo
           dx = ten*dx; dy = ten*dy; dz = ten*dz
           fx = fx + dx; fy = fy + dy; fz = fz + dz
           LF(1,jPt)=LF(1,jPt) -  sum( dx )*ds
           LF(2,jPt)=LF(2,jPt) -  sum( dy )*ds
           LF(3,jPt)=LF(3,jPt) -  sum( dz )*ds
        end if
     end do
     !$OMP END DO
     !$OMP END PARALLEL

     force_pt(1,:,i)=force_pt(1,:,i)+fx
     force_pt(2,:,i)=force_pt(2,:,i)+fy
     force_pt(3,:,i)=force_pt(3,:,i)+fz
  end do

end subroutine selfcontactForce_Fla_Body
