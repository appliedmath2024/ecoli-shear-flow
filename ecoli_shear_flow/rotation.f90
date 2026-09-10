SUBROUTINE rotation(Vec,theta,Tan)
  USE nrtype; USE commonval
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(INOUT)  :: Vec
  REAL(DP), INTENT(IN)  :: theta
  REAL(DP), DIMENSION(:), INTENT(In)  :: Tan

  REAL(DP), DIMENSION(3,3)  ::Tran, Iden, PSA
  REAL(DP), DIMENSION(3) :: Omg
  REAL(DP) :: detA
  INTEGER(I4B) :: i, j

  Iden=0.0_dp; Iden(1,1)=1.0_dp; Iden(2,2)=1.0_dp; Iden(3,3)=1.0_dp

  PSA(1,1)=Tan(1)*Tan(1); PSA(1,2)=Tan(1)*Tan(2); PSA(1,3)=Tan(1)*Tan(3) 
  PSA(2,1)=Tan(2)*Tan(1); PSA(2,2)=Tan(2)*Tan(2); PSA(2,3)=Tan(2)*Tan(3)
  PSA(3,1)=Tan(3)*Tan(1); PSA(3,2)=Tan(3)*Tan(2); PSA(3,3)=Tan(3)*Tan(3)
  Tran=0.0_dp; Tran(1,2)=-Tan(3); Tran(1,3)= Tan(2); Tran(2,1)=Tan(3);
  Tran(2,3)=-Tan(1); Tran(3,1)=-Tan(2); Tran(3,2)=Tan(1) 
  Tran=cos(theta)*Iden+(1.0_dp-cos(theta))*PSA+sin(theta)*Tran
  Omg=Matmul(Tran,Vec); Vec=Omg

END SUBROUTINE rotation
