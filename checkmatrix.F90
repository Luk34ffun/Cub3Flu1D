SUBROUTINE checkmatrix   !check the nonzero terms of the sparse matrix A_Phi
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k,r
INTEGER :: ii,jj,kk,s
REAL(myk) :: nzero

REAL(myk) :: a0,a_1,a_Nx,a_NxNy !coefficients
REAL(myk) :: d_1  !inv_den

d_1 = 1._myk/(dx*dy*dz)**2 
!define the coefficients of the non-zero terms
a0     = d_1*-2._myk*((dy*dz)**2+(dx*dz)**2+(dx*dy)**2)
a_1    = d_1*(dy*dz)**2
a_Nx   = d_1*(dx*dz)**2
a_NxNy = d_1*(dx*dy)**2

!Writing the non-zero terms of the sparse Matrix A_Phi & the relative nodes which they are referred to
OPEN(UNIT=1,FILE='checkmat.dat',STATUS='unknown',ACTION='WRITE')
WRITE(1,'(//,2x,A,//)') '---------------------------- NON-ZERO TERMS OF THE SPARSE MATRIX A_Phi ----------------------------'
WRITE(1,'(2x,A,F15.8,5x,A)') 'a0     = ', a0,' (if the i,j,k node is not a boundary node) '
WRITE(1,'(2x,A,F15.8)') 'a_1    = ', a_1
WRITE(1,'(2x,A,F15.8)') 'a_Nx   = ', a_Nx
WRITE(1,'(2x,A,F15.8,//)') 'a_NxNy = ', a_NxNy
DO k=1,Nz
 DO j=1,Ny
  DO i=1,Nx
   r=i+(j-1)*Nx+(k-1)*Nx*Ny
   DO kk=1,Nz
    DO jj=1,Ny
     DO ii=1,Nx
      s=ii+(jj-1)*Nx+(kk-1)*Nx*Ny
      IF (abs(A_Phi(r,s))>1E-10) THEN
       nzero = A_Phi(r,s)
       IF (r.EQ.s) THEN
        WRITE(1,'(1x,A)') ' --------------------------------------------------------------------'
        WRITE(1,'(1x,A,I4,A,I4,A,F18.9,10x,A,I3,A,I3,A,I3,A)')'| A_Phi (',r,',',s,') =',nzero,' node (',ii,',',jj,',',kk,') |' 
        WRITE(1,'(1x,A)') ' --------------------------------------------------------------------'
       ELSE
        WRITE(1,'(1x,A,I4,A,I4,A,F18.9,10x,A,I3,A,I3,A,I3,A)')'  A_Phi (',r,',',s,') =',nzero,' node (',ii,',',jj,',',kk,')' 
       END IF
     END IF
     END DO
    END DO
   END DO
   WRITE(1,'(///)')
  END DO
 END DO
END DO
CLOSE(1)
RETURN
END SUBROUTINE checkmatrix
