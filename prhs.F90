      subroutine prhs()
      use mod_fluid
      implicit none
!____._________________________________________________________________
      integer:: i,j,k
      real(myk)::  b_n
!____._________________________________________________________________

      b_n=dt*Brk(i_step,i_kutta)
!____._________________________________________________________________
      do k=1,Nz
      do j=1,Ny-1
      do i=1,Nx
         !auxiliary velocity field (V) equation
         P_ux(I,j,k)=ux(I,j,k)+b_n*H_ux(I,j,k)
         P_uy(i,J,k)=uy(i,J,k)+b_n*H_uy(i,J,k)
         
      enddo
      enddo
!____.j=Nr
      j=Ny
      do i=1,Nx
         P_ux(I,j,k)=ux(I,j,k)+b_n*H_ux(I,j,k)
      enddo
      enddo
!____.
      do k=1,Nz
      do j=1,Ny
      do i=1,Nx

         P_uz(i,j,K)=uz(i,j,K)+b_n*H_uz(i,j,K)

         !temperature equation
        P_th(i,j,k)=th(i,j,k)+b_n*H_th(i,j,k)
      enddo
      enddo
      enddo
!____._________________________________________________________________
      return
      end


