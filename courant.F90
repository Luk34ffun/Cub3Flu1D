SUBROUTINE courant
USE mod_fluid
IMPLICIT NONE
ni0=1.5E-5 ![m^2/s] viscosità cinematica aria @ 20°C
L0=0.2E-2  ![m]     diametro opening
U0 = Re*ni0/L0
ux_max = maxval(abs(ux))
uy_max = maxval(abs(uy))
uz_max = maxval(abs(uz))

Co = dt*(ux_max/dx+uy_max/dy+uz_max/dz)*U0
RETURN
END SUBROUTINE courant
