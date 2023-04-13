SUBROUTINE known_term
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k,r
!Definition of b - known terms array
DO k=1,Nz
 DO j=1,Ny
  DO i=1,Nx
   
   r=i+Nx*(j-1)+Nx*Ny*(k-1)
   b_Phi(r)= -div_V(i,j,k)
   
  END DO
 END DO
END DO
RETURN 
END SUBROUTINE known_term
