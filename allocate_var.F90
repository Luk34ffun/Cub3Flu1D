subroutine allocate_var
use mod_fluid
!grid variables
allocate(X(0:Nx+1))
allocate(Y(0:Ny+1))
allocate(Z(0:Nz+1))
!Physical variables
allocate(mu(0:Nx+1,0:Ny+1,0:Nz+1))
!Velocity field var
allocate(ux(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(uy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(uz(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(H_ux(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(H_uy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(H_uz(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(P_ux(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(P_uy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(P_uz(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(diff_x(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(diff_y(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(diff_z(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(hx(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(hy(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(hz(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(div_U(Nx,Ny,Nz))
ALLOCATE(vz(1:Nx,1:Ny,1))
!Viscous tensor var
allocate(Txx(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(Tyy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(Tzz(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(Txy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(Txz(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(Tyz(0:Nx+1,0:Ny+1,0:Nz+1))
!Force field var
ALLOCATE(fz(0:Nx+1,0:Ny+1,0:Nz+1))
!Temperature field var
ALLOCATE(th(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(H_th(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(P_th(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(l(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(diff_th(0:Nx+1,0:Ny+1,0:Nz+1))
!Iterative method variables
ALLOCATE(A_Phi(Nt,Nt))
ALLOCATE(b_Phi(Nt))
ALLOCATE(X0(Nt))
ALLOCATE(X_Phi(Nt,0:1))
ALLOCATE(div_V(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(grad_Phi_x(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(grad_Phi_y(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(grad_Phi_z(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(Phi(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(X_vx(Nt))
ALLOCATE(X_vy(Nt))
ALLOCATE(X_vz(Nt))
!
ALLOCATE(vz0(Nt,Nt,1))

return
endsubroutine allocate_var
