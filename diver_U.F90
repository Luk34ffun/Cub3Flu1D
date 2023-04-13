SUBROUTINE diver_U
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
     div_U(i,j,k) = (ux(I,j,k)-ux(I-1,j,k))/dx+(uy(i,J,k)-uy(i,J-1,k))/dy+(uz(i,j,K)-uz(i,j,K-1))/dz
    END DO
  END DO
END DO
maxdivU=MAXVAL(ABS(div_U))
!posmax=MAXLOC(div_U,MASK=ABS(div_U))
RETURN
END SUBROUTINE diver_U
