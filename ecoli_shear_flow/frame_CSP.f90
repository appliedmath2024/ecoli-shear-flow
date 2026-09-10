SUBROUTINE frame_CSP(LX,Forc,Torq,YCM,Fram)
  USE nrtype; USE commonval
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(OUT)  :: Fram
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: LX
  REAL(DP), DIMENSION(:), INTENT(IN)  :: Forc, Torq
  REAL(DP), DIMENSION(:), INTENT(OUT)  :: YCM
  REAL(DP), DIMENSION(3,3)  ::Tran, PSA, IPSA, MatM, SS, Fram0, Fram1, SNST, SOST, Test, NTorq, TranSS
  REAL(DP), DIMENSION(3) :: Diag
  REAL(DP) :: detA, theta
  INTEGER(I4B) :: i, j, LWORK=8, info
  REAL(DP), DIMENSION(8) :: WORK

  do i=1,3; YCM(i)=Forc(i)/(c_mss*Kvmax)+sum(LX(i,:))/Kvmax;
  end do

  do i=1,3;
     do j=1,3
        PSA(i,j)=sum(Coeff(i,:)*LX(j,:)) !! PSA=A
     end do;
  end do

  do i=1,3;
     do j=1,3;
        Tran(i,j)=sum(PSA(:,i)*PSA(:,j)) !! Tran=A^{t}A
     end do;
  end do

!!! sqrt of A =========================================================
  call sqrtm(Tran,MatM)  !! M=MatM

!!!Test=matmul(MatM,matM); print*,1,maxval(abs(Test-Tran))

!!! inverse of A ===================================================================
  detA=PSA(1,1)*PSA(2,2)*PSA(3,3)+PSA(1,2)*PSA(2,3)*PSA(3,1)+PSA(1,3)*PSA(2,1)*PSA(3,2) &
       -PSA(1,3)*PSA(2,2)*PSA(3,1)-PSA(1,2)*PSA(2,1)*PSA(3,3)-PSA(1,1)*PSA(2,3)*PSA(3,2)
  IPSA(1,1)=PSA(2,2)*PSA(3,3)-PSA(2,3)*PSA(3,2); IPSA(1,2)=-(PSA(1,2)*PSA(3,3)-PSA(1,3)*PSA(3,2));
  IPSA(1,3)=PSA(1,2)*PSA(2,3)-PSA(1,3)*PSA(2,2);
  IPSA(2,1)=-(PSA(2,1)*PSA(3,3)-PSA(2,3)*PSA(3,1)); IPSA(2,2)=(PSA(1,1)*PSA(3,3)-PSA(1,3)*PSA(3,1));
  IPSA(2,3)=-(PSA(1,1)*PSA(2,3)-PSA(1,3)*PSA(2,1))
  IPSA(3,1)=(PSA(2,1)*PSA(3,2)-PSA(2,2)*PSA(3,1)); IPSA(3,2)=-(PSA(1,1)*PSA(3,2)-PSA(1,2)*PSA(3,1));
  IPSA(3,3)=(PSA(1,1)*PSA(2,2)-PSA(1,2)*PSA(2,1))
  IPSA=IPSA/detA; !! IPSA= inv(A)

!!!Tran=matmul(PSA,IPSA); print*,2, Tran

!!! similar decompose of MatM S^{T}DiagS=MatM ==============================================
  TranSS=MatM
  Call DSYEV( 'V', 'U', 3, TranSS, 3, Diag, WORK, LWORK, INFO )
  do i=1,3;
     do j=1,3;
        SS(i,j)=TranSS(j,i);
     enddo;
  enddo !! SS=S, TranSS=tran(S)

!!!print*,3, TranSS
!!!print*,4,SS
!!!print*,5,Diag
!!!Tran=0.0_dp; Tran(1,1)=Diag(1); Tran(2,2)=Diag(2); Tran(3,3)=Diag(3)
!!!Test=matmul(Tran,SS); Tran=matmul(tranSS,Test);
!!!print*,6,maxval(abs(MatM-Tran))

  NTorq=0.0_dp; NTorq(1,2)=-Torq(3); NTorq(1,3)=Torq(2)
  NTorq(2,1)=Torq(3); NTorq(2,3)=-Torq(1); NTorq(3,1)=-Torq(2); NTorq(3,2)=Torq(1)
  Tran=matmul(NTorq,TranSS); SNST=matmul(SS,Tran)
  do i=1,3;
     do j=1,3;
        SOST(i,j)=SNST(i,j)/(Diag(i)+Diag(j));
     end do;
  end do; !!! SOST=Omg2

!!!
  Tran=matmul(SOST,SS);
  Fram1=matmul(TranSS,Tran); !! Fram1=Omg1
  Fram1=Fram1/c_mss; do i=1,3; Fram1(i,i)=Fram1(i,i)+1.0_dp; enddo; !! Fram1=R1
  Fram0=matmul(MatM,IPSA)
  Fram=matmul(Fram1,Fram0);
  
!!!do i=1,3; do j=1,3; Tran(i,j)=Fram(j,i); enddo; enddo; TranSS=matmul(Fram,Tran); print*,TranSS
  
END SUBROUTINE frame_CSP
