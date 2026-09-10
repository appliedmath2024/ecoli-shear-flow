SUBROUTINE initial_body(nref)
USE nrtype; USE commonval
IMPLICIT NONE
INTEGER(I4B), INTENT(INOUT) :: nref
REAL(DP) :: theta, psi, t, rr, ratio
INTEGER(I4B) :: n, i, j, ii, jj, kk, nvmax, ntmax, nv, nt, nr
REAL(DP), DIMENSION(3) :: x1, x2, x3, x4, Rx
INTEGER(I4B), DIMENSION(3) :: next
INTEGER(I4B), DIMENSION(40000,3) :: v
REAL(DP), DIMENSION(40000,3)  :: x

x=0.0_dp; v=0
ratio=R_cell/H_cell; theta=pi/2.0_dp; psi=1.0_dp/sqrt(2.0_dp); print*, ratio

if (ratio < 0.375_dp) then
   x(1,:)=(/0.0_dp,0.0_dp, R_cell+H_cell/2.0_dp/)
   x(2,:)=(/0.0_dp,0.0_dp,-R_cell-H_cell/2.0_dp/)
   do j=1,4; ii=2+j; jj=6+j; kk=10+j
      x(ii,:)=(/R_cell*cos((j-1.0_dp)*theta),R_cell*sin((j-1.0_dp)*theta),  H_cell/2.0_dp/)
      x(jj,:)=(/R_cell*cos((j-0.5_dp)*theta),R_cell*sin((j-0.5_dp)*theta),  0.0_dp/)
      x(kk,:)=(/R_cell*cos((j-1.0_dp)*theta),R_cell*sin((j-1.0_dp)*theta), -H_cell/2.0_dp/)
   end do

   v(1 ,:)=(/1 ,3 , 4/); v(2 ,:)=(/1, 4, 5/); v(3 ,:)=(/1, 5, 6/);  v(4 ,:)=(/1, 6 , 3/);
   v(5 ,:)=(/3 ,7 , 4/); v(6 ,:)=(/4, 8, 5/); v(7 ,:)=(/5, 9, 6/);  v(8 ,:)=(/6, 10 ,3/);
   v(9 ,:)=(/3 ,10, 7/); v(10,:)=(/4, 7, 8/); v(11,:)=(/5, 8, 9/);  v(12,:)=(/6, 9 ,10/);
   v(13,:)=(/11, 7, 10/); v(14,:)=(/12,8, 7/); v(15,:)=(/13,9, 8/);  v(16,:)=(/14,10 ,9/);
   v(17,:)=(/11,12, 7/); v(18,:)=(/12,13,8/); v(19,:)=(/13,14,9/);  v(20,:)=(/14,11,10/);
   v(21,:)=(/11,2, 12/); v(22,:)=(/12,2,13/); v(23,:)=(/13,2,14/);  v(24,:)=(/14,2,11/);
   
   nv=14; nt=24;  Ktmax=nt*4**nref; do nr=1,nref; call refine(x,v,nv,nt); enddo
   
elseif ( (ratio >=0.375_dp) .AND. (ratio < 0.75_dp) ) then
   x(1,:)=(/0.0_dp,0.0_dp, R_cell+H_cell/2.0_dp/)
   x(2,:)=(/0.0_dp,0.0_dp,-R_cell-H_cell/2.0_dp/)
   do j=1,4; ii=2+j; jj=6+j
      x(ii,:)=(/R_cell*cos((j-1.0_dp)*theta),R_cell*sin((j-1.0_dp)*theta),  H_cell/2.0_dp/)
      x(jj,:)=(/R_cell*cos((j-0.5_dp)*theta),R_cell*sin((j-0.5_dp)*theta), -H_cell/2.0_dp/)
      
   end do
  
   v(1 ,:)=(/1 ,3 , 4/); v(2 ,:)=(/1, 4, 5/); v(3 ,:)=(/1, 5, 6/);  v(4 ,:)=(/1, 6 , 3/);
   v(5 ,:)=(/3 ,7 , 4/); v(6 ,:)=(/4, 8, 5/); v(7 ,:)=(/5, 9, 6/);  v(8 ,:)=(/6, 10 ,3/);
   v(9 ,:)=(/3 ,10, 7/); v(10,:)=(/4, 7, 8/); v(11,:)=(/5, 8, 9/);  v(12,:)=(/6, 9 ,10/);
   v(13,:)=(/7 ,2 , 8/); v(14,:)=(/8, 2, 9/); v(15,:)=(/9,2, 10/);  v(16,:)=(/10,2 , 7/);
   
   nv=10; nt=16;  Ktmax=nt*4**nref; do nr=1,nref; call refine(x,v,nv,nt); enddo
   
elseif ( (ratio >=0.75_dp) .AND. (ratio < 1.5_dp) ) then
   x(1,:)=(/0.0_dp,0.0_dp, R_cell+H_cell/2.0_dp/)
   x(2,:)=(/0.0_dp,0.0_dp,-R_cell-H_cell/2.0_dp/)
   do j=1,4; ii=2+j; t=j-0.5_dp
      x(ii,:)=(/psi*R_cell*cos(t*theta),psi*R_cell*sin(t*theta),psi*R_cell+H_cell/2.0_dp/)
   end do
   do j=1,8; ii=6+j; t=0.5_dp*(j-1.0_dp)
      x(ii,:)=(/R_cell*cos(t*theta),R_cell*sin(t*theta), H_cell/2.0_dp/)
   end do
   do j=1,4; ii=14+j; t=j-0.5_dp+0.25_dp
      x(ii,:)=(/psi*R_cell*cos(t*theta),psi*R_cell*sin(t*theta),-psi*R_cell-H_cell/2.0_dp/)
   end do
   do j=1,8; ii=18+j; t=0.5_dp*(j-1.0_dp)+0.25_dp
      x(ii,:)=(/R_cell*cos(t*theta),R_cell*sin(t*theta), -H_cell/2.0_dp/)
   end do
   v(1 ,:)=(/1 ,3 , 4/); v(2 ,:)=(/1, 4, 5/); v(3 ,:)=(/1, 5, 6/);  v(4 ,:)=(/1, 6, 3/);
   
   v(5 ,:)=(/3 ,9 , 4/); v(6 ,:)=(/4, 11,5/); v(7 ,:)=(/5, 13,6/);  v(8 ,:)=(/6, 7, 3/);
   
   v(9 ,:)=(/3 ,8 , 9/); v(10,:)=(/4, 9,10/); v(11,:)=(/4,10,11/);  v(12,:)=(/5,11,12/);
   v(13,:)=(/5 ,12,13/); v(14,:)=(/6,13,14/); v(15,:)=(/6,14, 7/);  v(16,:)=(/3, 7, 8/);
   
   v(17,:)=(/15,21,20/); v(18,:)=(/16,22,21/); v(19,:)=(/16,23,22/);  v(20,:)=(/17,24,23/);
   v(21,:)=(/17,25,24/); v(22,:)=(/18,26,25/); v(23,:)=(/18,19,26/);  v(24,:)=(/15,20,19/);

   v(25,:)=(/15,16,21/); v(26,:)=(/16,17,23/); v(27,:)=(/17,18,25/);  v(28,:)=(/18,15,19/);

      
   v(29,:)=(/2 ,16,15/); v(30,:)=(/2 ,17,16/); v(31,:)=(/2 ,18,17/);  v(32,:)=(/2 ,15,18/);

   
   do j=1,7; kk=32+j; ii=6+j; jj=18+j; v(kk,:)=(/ii,jj,ii+1/);
   end do; v(40,:)=(/14,26,7/);
   do j=1,7; kk=40+j; ii=18+j; jj=7+j; v(kk,:)=(/ii,ii+1,jj/);
   end do; v(48,:)=(/26,19,7/);
         
   nv=26; nt=48; Ktmax=nt*4**nref; do nr=1,nref; call refine(x,v,nv,nt); enddo
   
end if

do i=1,40000; if (sum(abs(x(i,:)))==0) exit ; enddo
   Kvmax=i-1;
   print*,'Kvtmax=',Kvmax,Ktmax, v(Ktmax,:), v(Ktmax+1,:)
   
   Allocate(Index0(Ktmax,3)); Allocate(LX0(Kvmax,3));
   Allocate(Length0(Ktmax,3)); Allocate(Radius0(Kvmax)); Allocate(Area0(Kvmax))
   Index0=v(1:Ktmax,:); do i=1,Kvmax; LX0(i,:)=x(i,:); end do
   
   if (A_cell > 0.000001_dp) then
      do i=1,Kvmax
         Rx=0.0_dp; x1=x(i,:); x1(3)=0.0_dp
         x2(1)=A_cell*sin(2.0_dp*pi*x(i,3)/H_cell); x2(2)=0.0_dp; x2(3)=x(i,3)
         if (x(i,3) > H_cell/2.0_dp) then
            Rx(1)=-2.0_dp*pi*A_cell/H_cell; Rx(3)=1.0_dp
         else if (x(i,3) < -H_cell/2.0_dp) then
            Rx(1)=-2.0_dp*pi*A_cell/H_cell; Rx(3)=1.0_dp
         else
            Rx(1)=2.0_dp*pi*A_cell/H_cell*cos(2.0_dp*pi*x(i,3)/H_cell); Rx(3)=1.0_dp
         end if
         t=sqrt(Rx(1)**2+Rx(3)**2); Rx=Rx/t
         LX0(i,1)= Rx(3)*x(i,1)+x2(1)
         LX0(i,2)=x(i,2)
         LX0(i,3)=-Rx(1)*x(i,1)+x2(3)
      end do
      !! rotation to make the z-diretional ~~~
      Rx(1)=-2.0_dp*pi*A_cell/H_cell; Rx(2)=0.0_dp; Rx(3)=1.0_dp
      t=sqrt(Rx(1)**2+Rx(3)**2); Rx=Rx/t; t=acos(Rx(3)); print*,'theta', t*180/pi
      do i=1,Kvmax
         x1(1)= cos(t)*LX0(i,1)+sin(t)*LX0(i,3)
         x1(3)=-sin(t)*LX0(i,1)+cos(t)*LX0(i,3)
         LX0(i,1)=x1(1); LX0(i,3)=x1(3)
      end do
   end if
   
   do i=1,Ktmax
      ii=Index0(i,1); jj=Index0(i,2); kk=Index0(i,3)
      x1=LX0(ii,:); x2=LX0(jj,:); x3=LX0(kk,:)
      x4=x2-x1; Length0(i,1)=sqrt(x4(1)**2+x4(2)**2+x4(3)**2);
      x4=x3-x2; Length0(i,2)=sqrt(x4(1)**2+x4(2)**2+x4(3)**2);
      x4=x1-x3; Length0(i,3)=sqrt(x4(1)**2+x4(2)**2+x4(3)**2);
   end do
   
   do i=1,Kvmax
      x4=LX0(i,:); Radius0(i)=sqrt(x4(1)**2+x4(2)**2+x4(3)**2);
   end do
   
   Area0=0.0_dp
   do i=1,Ktmax
      ii=Index0(i,1); jj=Index0(i,2); kk=Index0(i,3)
      x1=LX0(ii,:); x2=LX0(jj,:); x3=LX0(kk,:)
      x4(1)=(x2(2)-x1(2))*(x3(3)-x1(3))-(x2(3)-x1(3))*(x3(2)-x1(2))
      x4(2)=(x2(3)-x1(3))*(x3(1)-x1(1))-(x2(1)-x1(1))*(x3(3)-x1(3))
      x4(3)=(x2(1)-x1(1))*(x3(2)-x1(2))-(x2(2)-x1(2))*(x3(1)-x1(1))
      rr=sqrt(x4(1)**2+x4(2)**2+x4(3)**2)/2.0_dp
      Area0(ii)=Area0(ii)+rr/3.0_dp;
      Area0(jj)=Area0(jj)+rr/3.0_dp
      Area0(kk)=Area0(kk)+rr/3.0_dp
   enddo


CONTAINS
  SUBROUTINE refine(x,v,nv,nt)
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: x
  INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT) :: v
  INTEGER(I4B), INTENT(INOUT) :: nv, nt
  INTEGER(I4B), DIMENSION(3)  :: next3, vnew
  REAL(DP), DIMENSION(3)  :: xout
  INTEGER(I4B), DIMENSION(:,:), Allocatable  :: a
  INTEGER(I4B) :: k, j, v1, v2
  vnew=0; next3=(/2,3,1/); Allocate(a(nv,nv)); a=0;
  do k=1,nt;  do j=1,3
      v1=v(k,next3(j)); v2=v(k,next3(next3(j)));
      if(a(v1,v2)==0) then
        nv=nv+1; vnew(j)=nv
        call locate(v1,v2,x,xout); x(nv,:)=xout;
        a(v1,v2)=nv; a(v2,v1)=nv
      else
        vnew(j)=a(v1,v2)
      end if
    end do
     v(1*nt+k,:)=(/v(k,1),vnew(3),vnew(2)/)
     v(2*nt+k,:)=(/v(k,2),vnew(1),vnew(3)/)
     v(3*nt+k,:)=(/v(k,3),vnew(2),vnew(1)/)
     v(k,:)=vnew;
  end do
  nt=nt*4;deallocate(a)
  END SUBROUTINE refine

  SUBROUTINE locate(v1,v2,x,xout)
  INTEGER(I4B), INTENT(IN) :: v1, v2
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: x
  REAL(DP), DIMENSION(:), INTENT(OUT) :: xout
  REAL(DP) :: rr,  r1, r2, rs
  xout=(x(v1,:)+x(v2,:))/2.0_dp
  r1=sqrt( x(v1,1)**2 + x(v1,2)**2 )
  r2=sqrt( x(v2,1)**2 + x(v2,2)**2 )
  if ( (r1 < (R_cell - 10.0_dp**(-10))) .OR. (r2 < (R_cell - 10.0_dp**(-10))) ) then
     rs=sign(1.0_dp,xout(3))
     xout(3)=xout(3)-rs*H_cell/2.0_dp
     rr=sqrt(xout(1)**2+xout(2)**2+xout(3)**2)
     xout=R_cell*xout/rr; xout(3)=xout(3)+rs*H_cell/2.0_dp
  else
     rr=sqrt(xout(1)**2+xout(2)**2)
     xout(1:2)=R_cell*xout(1:2)/rr
  end if
  END SUBROUTINE locate

END SUBROUTINE
