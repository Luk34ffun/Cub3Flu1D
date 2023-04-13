subroutine rhs
use mod_fluid

implicit none
integer::i,j,k

!-----------Components of the viscous stress tensor T-----------------------------------------------------------------
!! Here I suppose that the extradiagonal components are placed at the edges of the cell
!! while the diagonal ones are placed at the center of the cell
do k=1,Nz
  do j=1,Ny
    do i=1,Nx

Txy(I,J,k)=0.25_myk*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j+1,k)+mu(i+1,j+1,k))*((ux(I,j+1,k)-ux(I,j,k))/dy+(uy(i+1,J,k)-&
uy(i,J,k))/dx)

Txz(I,j,K)=0.25_myk*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k+1)+mu(i+1,j,k+1))*((ux(I,j,k+1)-ux(I,j,k))/dz+(uz(i+1,j,K)-&
uz(i,j,K))/dx)

Tyz(i,J,K)=0.25_myk*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*((uy(i,J,k+1)-uy(i,J,k))/dz+(uz(i,j+1,K)-&
uz(i,j,K))/dy)


    Txx(i,j,k)=2.0_myk*mu(i,j,k)*((ux(I,j,k)-ux(I-1,j,k))/dx)

    Tyy(i,j,k)=2.0_myk*mu(i,j,k)*((uy(i,J,k)-uy(i,J-1,k))/dy)

    Tzz(i,j,k)=2.0_myk*mu(i,j,k)*((uz(i,j,K)-uz(i,j,K-1))/dz)

    enddo
  enddo
enddo
WRITE(*,'(6x,A)') '- stress tensor components written'
!----------- Divergence of T ---------------
!!<div(T),i> = d/dx(Txx)+d/dy(Txy)+d/dz(Txz)
!!<div(T),j> = d/dx(Tyx)+d/dy(Tyy)+d/dz(Tyz)
!!<div(T),k> = d/dx(Tzx)+d/dy(Tzy)+d/dz(Tzz)
do k=1,Nz
  do j=1,Ny
    do i=1,Nx
!Diffusive term (dimensionless) : diff = 1/Re * Laplac(U)
diff_x(I,j,k) = 0.5_myk*inv_re*((Txx(i+1,j,k)-Txx(i,j,k))/dx+(Txy(I,J,k)-Txy(I,J-1,k))/dy+(Txz(I,j,K)-Txz(I,j,K-1))/dz)!*2._myk !il prof moltiplica per 0.5
diff_y(i,J,k) = 0.5_myk*inv_re*((Txy(I,J,k)-Txy(I-1,J,k))/dx+(Tyy(i,j+1,k)-Tyy(i,j,k))/dy+(Tyz(i,J,K)-Tyz(i,J,K-1))/dz)!*2._myk
diff_z(i,j,K) = 0.5_myk*inv_re*((Txz(I,j,K)-Txz(I-1,j,K))/dx+(Tyz(i,J,K)-Tyz(i,J-1,K))/dy+(Tzz(i,j,k+1)-Tzz(i,j,k))/dz)!*2._myk
    END DO
  END DO
END DO
WRITE(*,'(6x,A)') '- diffusive part written'
!______________________________________________________________________
!Convective term (dimensionless) : h = u_j * d(u_i)/dx_j= (hx,hy,hz)^t

DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx


hx(I,j,k) = (ux(I,j,k)*(ux(I+1,j,k)-ux(I-1,j,k))*0.5_myk/dx&
+0.25_myk*(uy(i,J,k)+uy(i+1,J,k)+uy(i+1,J-1,k)+uy(i,J-1,k))*(ux(I,j+1,k)-ux(I,j-1,k))*0.5_myk/dy&
+0.25_myk*(uz(i,j,K)+uz(i+1,j,K)+uz(i+1,j,K-1)+uz(i,j,K-1))*(ux(I,j,k+1)-ux(I,j,k-1))*0.5_myk/dz)




hy(i,J,k) = (0.25_myk*(ux(I,j,k)+ux(I-1,j,k)+ux(I-1,j+1,k)+ux(I,j+1,k))*(uy(i+1,J,k)-uy(i-1,J,k))*0.5_myk/dx&
+uy(i,J,k)*(uy(i,J+1,k)-uy(i,J-1,k))*0.5_myk/dy&
+0.25_myk*(uz(i,j,K)+uz(i,j+1,K)+uz(i,j+1,K-1)+uz(i,j,K-1))*(uy(i,J,k+1)-uy(i,J,k-1))*0.5_myk/dz)




hz(i,j,K) = (0.25_myk*(ux(I,j,k)+ux(I,j,k+1)+ux(I-1,j,k+1)+ux(I-1,j,k))*(uz(i+1,j,K)-uz(i-1,j,K))*0.5_myk/dx&
+0.25_myk*(uy(i,J,k)+uy(i,J,k+1)+uy(i,J-1,k+1)+uy(i,J-1,k))*(uz(i,j+1,K)-uz(i,j-1,K))*0.5_myk/dy&
+uz(i,j,K)*(uz(i,j,K+1)-uz(i,j,K-1))*0.5_myk/dz)


    END DO
  END DO
END DO

IF (mod(it,itout) == 0) THEN
CALL stampa(7)
END IF
WRITE(*,'(6x,A)') '- convective part written'

!Forcing term
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx

fz(i,j,k) = Ri*(th(i,j,k)*Th_0-Th_s)/dTh0
    END DO
  END DO
END DO
WRITE(*,'(6x,A)') '- forcing term written'
!Resulting right hand side term
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx

    H_ux(I,j,k)=diff_x(I,j,k)-hx(I,j,k)
    H_uy(i,J,k)=diff_y(i,J,k)-hy(i,J,k)
    H_uz(i,j,K)=diff_y(i,j,K)-hz(i,j,K) + fz(i,j,k)

    END DO
  END DO
END DO
WRITE(*,'(6x,A)') '- RHS of the v-field eq. written'
!---- Temperature field (theta)----------------------------------------------------
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx

!---- convective term l=u_i*d(th)/dx_i
    l(i,j,k) = 0.5_myk*((ux(I,j,k)*(th(i+1,j,k)-th(i,j,k))/dx+ux(I-1,j,k)*(th(i,j,k)-th(i-1,j,k))/dx)+&
                      (uy(i,J,k)*(th(i,j+1,k)-th(i,j,k))/dy+uy(i,J-1,k)*(th(i,j,k)-th(i,j-1,k))/dy)+&
                      (uz(i,j,K)*(th(i,j,k+1)-th(i,j,k))/dz+uz(i,j,K-1)*(th(i,j,k)-th(i,j,k-1))/dz))

!---- diffusive term 1/(Re*Pr) * lapl(th) 
     diff_th(i,j,k) = 1.0_myk/(Re*Pr)*((th(i+1,j,k)-2*th(i,j,k)+th(i-1,j,k))/dx**2+(th(i,j+1,k)-&
                      2*th(i,j,k)+th(i,j-1,k))/dy**2+(th(i,j,k+1)-2*th(i,j,k)+th(i,j,k-1))/dz**2)   
     
     !--- RHS temperature term
     H_th(i,j,k) = -l(i,j,k)+diff_th(i,j,k)
    
    END DO
  END DO
END DO
WRITE(*,'(6x,A)') '- RHS of the temp eq. written'
return
endsubroutine rhs
