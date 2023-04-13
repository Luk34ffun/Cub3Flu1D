SUBROUTINE associate_var
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k,r

DO k=1,Nz
 DO j=1,Ny
  DO i=1,Nx
  r=i+Nx*(j-1)+Nx*Ny*(k-1)
  Phi(i,j,k)=x_Phi(r,0)
  X_vx(r)= ux(I,j,k)
  X_vy(r)= uy(i,J,k)
  X_vz(r)= uz(i,j,K)
  END DO  
 END DO
END DO


RETURN
END SUBROUTINE associate_var
