SUBROUTINE grad_Phi()
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k

DO k=1,Nz
  DO j=1,Ny
    DO i=0,Nx
    
     grad_Phi_x(I,j,k) = (Phi(i+1,j,k)-Phi(i,j,k))/dx 
    
    END DO
  END DO
END DO

DO k=1,Nz
  DO j=0,Ny
    DO i=1,Nx

     grad_Phi_y(i,J,k) = (Phi(i,j+1,k)-Phi(i,j,k))/dy
    
    END DO
  END DO
END DO

DO k=0,Nz!k=0,Nz
  DO j=1,Ny
    DO i=1,Nx
    
     grad_Phi_z(i,j,K) = (Phi(i,j,k+1)-Phi(i,j,k))/dz 
    
    END DO
  END DO
END DO

RETURN
END SUBROUTINE grad_Phi
