SUBROUTINE vzmaxbound
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j

DO j=1,Ny
  DO i=1,Nx
  vz0(i,j,1)=uz(i,j,1)
  END DO
END DO
vz0max=maxval(abs(vz0))
DO j=1,Ny
  DO i=1,Nx
   IF (abs(vz0(i,j,1)) == vz0max) THEN
   i_max=i
   j_max=j
   GOTO 100
   END IF
  END DO
END DO
100 CONTINUE
IF (MOD(it,itout)==0) THEN
CALL stampa(3)
END IF
RETURN
END SUBROUTINE vzmaxbound
