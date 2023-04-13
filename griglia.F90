subroutine griglia()
use mod_fluid
implicit none
integer:: i,j,k
!____.__________________________________________________________________
dx = Xmax/float(Nx)
X(0) = -dx/2.0_myk
do i=1,Nx+1
   X(i) = X(i-1) + dx
enddo
!____.__________________________________________________________________
dy = Ymax/float(Ny)
Y(0) = -dy/2.0_myk
do j=1,Ny+1
   Y(j) = Y(j-1) + dy
enddo
!____.__________________________________________________________________
dz = Zmax/float(Nz)
Z(0) = -dz/2.0_myk
do k=1,Nz+1
   Z(k) = Z(k-1) + dz
enddo
WRITE(*,'(A,/)') ' X-Y-Z divisions / spatial discretization steps'
WRITE(*,'(A,I3,A,E13.8)') ' [Nx] = ',Nx,' -----> dx = ',dx
WRITE(*,'(A,I3,A,E13.8)') ' [Ny] = ',Ny,' -----> dy = ',dy
WRITE(*,'(A,I3,A,E13.8,//)') ' [Nz] = ',Nz,' -----> dz = ',dz
!____.__________________________________________________________________
return
end
