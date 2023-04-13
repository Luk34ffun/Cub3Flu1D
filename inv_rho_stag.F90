SUBROUTINE inv_rho_stag
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k

DO k=1,Nz
  DO j=1,Ny
    DO i=0,Nx
     
     inv_rho_x(I,j,k) = 2/(rho(i,j,k)+rho(i+1,j,k))
    
    END DO
  END DO
END DO

DO k=1,Nz
  DO j=0,Ny
    DO i=1,Nx
     
     inv_rho_y(i,J,k) = 2/(rho(i,j,k)+rho(i,j+1,k))
    
    END DO
  END DO
END DO

DO k=0,Nz
  DO j=1,Ny
    DO i=1,Nx
     
     inv_rho_z(i,j,K) = 2/(rho(i,j,k)+rho(i,j,k+1))
    
    END DO
  END DO
END DO
