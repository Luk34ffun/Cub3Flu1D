SUBROUTINE vol_flowrate
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k
REAL(myk), DIMENSION (Nx*Ny) :: tmp
k=1
DO j=1,Ny
  DO i=1,Nx
    tmp(k) = uz(i,j,0) *dx*dy
    k = k+1
  END DO
END DO
q_dot(1,1) = sum(tmp)

k=1
DO j=1,Ny
  DO i=1,Nx
    tmp(k) = uz(i,j,Nz) *dx*dy
    k = k+1
  END DO
END DO
q_dot(1,2) = sum(tmp)


RETURN
END SUBROUTINE vol_flowrate
