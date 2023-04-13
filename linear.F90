      subroutine linear()
      use mod_fluid
      implicit none
      !In questa subroutine u = V (campo ausiliario)
!____._________________________________________________________________
      integer:: i,j,k
      real(myk):: a_n
!____._________________________________________________________________
      a_n=dt*Ark(i_step,i_kutta)
!____._________________________________________________________________
      do k=1,Nz
      do j=1,Ny-1
      do i=1,Nx


         ux(I,j,k)=P_ux(I,j,k)+a_n*H_ux(I,j,k)         
         uy(i,J,k)=P_uy(i,J,k)+a_n*H_uy(i,J,k)
      enddo
      enddo
!____._______________________________________________________________
      j=Ny
      do i=1,Nx
         ux(I,j,k)=P_ux(I,j,k)+a_n*H_ux(I,j,k)
      enddo
      enddo
!____._______________________________________________________________
      do k=1,Nz
      do j=1,Ny
      do i=1,Nx
         uz(i,j,K)=P_uz(i,j,K)+a_n*H_uz(i,j,K)


         !---- temperature field -----------------------------------
         th(i,j,k) = P_th(i,j,k)+a_n*H_th(i,j,k)
      enddo
      enddo
      enddo

      !Vz
      DO j=1,Ny
        DO i=1,Nx
        
         vz(i,j,1) = uz(i,j,Nz)!store the value of auxiliary field 
                                !that will be used to correct the B.C.
                                !on top surface (opening) ,referred to grad(Phi)
        END DO
      END DO
!____._________________________________________________________________
      return
      end


