SUBROUTINE U_field()
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k

DO k=1,Nz
  DO j=1,Ny
    DO i=0,Nx

    ux(I,j,k) = ux(I,j,k)+grad_Phi_x(I,j,k) 
    
    END DO
  END DO
END DO


DO k=1,Nz
  DO j=0,Ny
    DO i=1,Nx
  
    uy(i,J,k) = uy(i,J,k)+grad_Phi_y(i,J,k)
 
    END DO
  END DO
END DO


DO k=0,Nz
  DO j=1,Ny
    DO i=1,Nx

    uz(i,j,K) = uz(i,j,K)+grad_Phi_z(i,j,K)

    END DO
  END DO
END DO

RETURN
END SUBROUTINE U_field
