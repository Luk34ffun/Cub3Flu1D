SUBROUTINE clean
USE mod_fluid
IMPLICIT NONE
!grid variables
deallocate(X,Y,Z,mu,ux,uy,uz,H_ux,H_uy,H_uz,P_ux,P_uy,P_uz)
DEALLOCATE(Txx,Tyy,Tzz,Txy,Txz,Tyz,diff_x,diff_y,diff_z,hx,hy,hz)
DEALLOCATE(A_Phi,b_Phi,X0,X_vx,X_vy,X_vz,div_V)
DEALLOCATE(grad_Phi_x,grad_Phi_y,grad_Phi_z)
DEALLOCATE(th,H_th,P_th)
RETURN
END SUBROUTINE clean
