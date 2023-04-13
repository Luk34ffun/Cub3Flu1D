SUBROUTINE Phi_iter(metho)
USE mod_fluid
IMPLICIT NONE
INTEGER :: i
INTEGER,INTENT(in) :: metho
REAL(myk),DIMENSION(Nt) :: tmp_Phi
REAL(myk) :: err_iter
!CALL buildmatrix

last_i=0
err_iter=1e-1 !guess on the starting error : initializing variable

SELECT CASE (metho)

!JACOBI ITERATION (method==1)
CASE(1)

DO i=0,max_iter-1
!While loop
 IF(err_iter<max_err) EXIT
 DO m=1,Nt

  DO n=1,Nt
  tmp_Phi(n)=A_Phi(m,n)*x_Phi(n,0)
  END DO

 x_Phi(m,1)=(1._myk/A_Phi(m,m))*(b_Phi(m)-(sum(tmp_Phi)-A_Phi(m,m)*x_Phi(m,0)))
 END DO
 err_iter= MAXVAL(ABS(x_Phi(:,1)-x_Phi(:,0)))
 x_Phi(:,0)=x_Phi(:,1)
last_i=i+1
END DO

!GAUSS - SEIDEL ITERATION (method==2)
CASE(2)

DO i=0,max_iter-1
!While loop
 IF(err_iter<max_err) EXIT
 DO m=1,Nt

  DO n=1,Nt
  IF (m>n) THEN
  tmp_Phi(n)=A_Phi(m,n)*x_Phi(n,1)
  ELSE
  tmp_Phi(n)=A_Phi(m,n)*x_Phi(n,0)
  END IF
  END DO

 x_Phi(m,1)=(1._myk/A_Phi(m,m))*(b_Phi(m)-(sum(tmp_Phi)-A_Phi(m,m)*x_Phi(m,0)))
 END DO
 err_iter= MAXVAL(ABS(x_Phi(:,1)-x_Phi(:,0)))
 x_Phi(:,0)=x_Phi(:,1)
 last_i=i+1
END DO

!SUCCESSIVE OVER RELAXATION ITERATION (method==3)
CASE(3)
!WRITE(*,*) 'inizio sor'
DO i=0,max_iter-1
!While loop
 IF(err_iter<max_err) EXIT
 
 DO m=1,Nt
!WRITE(*,*) 'incr_righe'
  DO n=1,Nt
!WRITE(*,*) 'incr_col'
  IF (m>n) THEN
 !WRITE(*,*) 'se m>n allora'  
  tmp_Phi(n)=A_Phi(m,n)*x_Phi(n,1)
 !WRITE(*,*) 'scrivo vett temp'  
  ELSE
  tmp_Phi(n)=A_Phi(m,n)*x_Phi(n,0)
  END IF
  END DO

 x_Phi(m,1)=(1._myk/A_Phi(m,m))*(b_Phi(m)-(sum(tmp_Phi)-A_Phi(m,m)*x_Phi(m,0)))
 x_Phi(m,1)=omg*x_Phi(m,1)+(1._myk-omg)*x_Phi(m,0)
! WRITE(*,*) ' scrivo x(',m,',1)'
 END DO
 err_iter= MAXVAL(ABS(x_Phi(:,1)-x_Phi(:,0)))
 !WRITE(*,*) 'aggiorno l errore',m,',',n
 x_Phi(:,0)=x_Phi(:,1)
 !WRITE(*,*) ' soluzione :',m
 last_i=i+1
END DO
IF (it==itmax) THEN
OPEN(UNIT=1,FILE=path_fold//'iter_sol.bin',STATUS='unknown',FORM='unformatted',ACTION='WRITE')
WRITE(1) x_Phi(:,0)
WRITE(*,*) it
CLOSE(1)
END IF
END SELECT

RETURN
END SUBROUTINE Phi_iter
