SUBROUTINE stampa(cs)
USE mod_fluid
IMPLICIT NONE
INTEGER, INTENT(IN) :: cs
INTEGER :: i,j,k,r,num
REAL(myk) :: time 
time=it*dt
!-------------------------------------------------------------------------
SELECT CASE (cs)

CASE(1)
!------------------- Printing Binary output (vel/temp field) -----------
write(fileres(1:1),'(A)') ver        !a
write(fileres(4:9),'(I6.6)') it      ! ve000it.bin
opfile = path_fold//'bin/'//trim(fileres)
write(*,*) 'Writing the file ',opfile
open(unit=1,file=opfile,form='unformatted')
write(1) ux,uy,uz
CLOSE(1)

!---------------------------- Printing Temperature Field ---------------------------------
!temp bin
WRITE(fileth(1:2),'(A)') 'th'
WRITE(fileth(3:8),'(I6.6)') it
opfileth=path_fold//'bin/'//fileth//'.bin'
opfileth=TRIM(opfileth)
WRITE(*,*)'Writing the file ',opfileth
OPEN(UNIT=9,FILE=opfileth,FORM='unformatted')
WRITE(9) th
CLOSE(9)



resth=path_fold//'temperature/'//fileth//extns
resth=TRIM(resth)
WRITE(*,*) 'Writing the file ' ,resth

OPEN(UNIT=9,FILE=resth,STATUS='unknown',ACTION='WRITE')
WRITE(9,'(/,5x,A,//)')'======= TEMPERATURE FIELD (dimensionless) ======='

DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
    WRITE(9,220) i,j,k, th(i,j,k)
220 FORMAT (1x,'Theta*(',I2,',',I2,',',I2,') = ',F18.9)

    END DO
  END DO
END DO
CLOSE(9)

CASE(2)

!----------------------- Printing V ------------------------------

!WRITE(filename(1:2),'(A)') 'vf'
!WRITE(filename(3:8),'(I6.6)') it
!resfile=path_fold//'V/'//TRIM(filename)//extns
!WRITE(*,*) 'Writing the file ' ,resfile
!OPEN(UNIT=9,FILE=resfile,STATUS='unknown',ACTION='WRITE')
!WRITE(9,'(/,15x,A,//)')'======= V-FIELD (dimensionless auxiliary velocity field) ======='
!WRITE(9,'(1x,A,I3,A,I3,A,6x,A,F6.3)') ' approx position of bottom entry center = (',int(x_1/dx),',',int(y_1/dy),', 0 )' &
!&,'r_bottom = ',r_1
!WRITE(9,'(1x,A,I3,A,I3,A,I3,A,6x,A,F6.3,//)') ' approx position of top entry center    = (',int(x_2/dx),',',int(y_2/dy),',',Nz,')' &
!&,'r_top    = ',r_2
!DO k=1,Nz
!  DO j=1,Ny
!    DO i=1,Nx
!     r=i+Nx*(j-1)+Nx*Ny*(k-1)
!     WRITE(9,140) i,j,k, X_vx(r),i,j,k,X_vy(r),i,j,k,X_vz(r)
!     140 FORMAT (1x,'Vx(',I2,',',I2,',',I2,') = ',F18.9,1x,15x, &
!     'Vy(',I2,',',I2,',',I2,') = ',F18.9,15x, &
!     'Vz(',I2,',',I2,',',I2,') = ',F18.9)
!     END DO
!   END DO
!END DO
!CLOSE(9)

!---------------------Printing Phi-----------------------------------

SELECT CASE (method)
!JACOBI (method==1)
CASE(1)
!Phi
WRITE(filePhi(1:3),'(A)') jac
WRITE(filePhi(4:9),'(I6.6)') it
resPhi=path_fold//'Phi/jacobi/'//filePhi//extns
WRITE(*,*) 'Writing the file ' ,resPhi
OPEN(UNIT=9,FILE=resPhi,STATUS='unknown',ACTION='WRITE')
WRITE(9,*) '======= JACOBI METHOD ======='
WRITE(9,210) last_i
210 FORMAT (1x,//,'------ iteration n°(',I5,') of solution search :  ------ ',//)
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
    r=i+Nx*(j-1)+Nx*Ny*(k-1)

WRITE(9,150) i,j,k, x_Phi(r,0)
150 FORMAT (1x,'Phi(',I2,',',I2,',',I2,') = ',F18.9)
    END DO
  END DO
END DO
CLOSE(9)

!GAUSS SEIDEL (method==2)
CASE(2)
!Phi
WRITE(filePhi(1:3),'(A)') gs
WRITE(filePhi(4:9),'(I6.6)') it
resPhi=path_fold//'Phi/gs/'//filePhi//extns
WRITE(*,*) 'Writing the file ' ,resPhi
OPEN(UNIT=9,FILE=resPhi,STATUS='unknown',ACTION='WRITE')
WRITE(9,*) '======= GAUSS-SEIDEL METHOD ======='
WRITE(9,210) last_i
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
    r=i+Nx*(j-1)+Nx*Ny*(k-1)

WRITE(9,150) i,j,k, x_Phi(r,0)
    END DO
  END DO
END DO
CLOSE(9)

!S.O.R. (method==3)
CASE(3)
!Phi
WRITE(filePhi(1:3),'(A)') sor
WRITE(filePhi(4:9),'(I6.6)') it
resPhi=path_fold//'Phi/sor/'//filePhi//extns
WRITE(*,*) 'Writing the file ' ,resPhi
OPEN(UNIT=9,FILE=resPhi,STATUS='unknown',ACTION='WRITE')
WRITE(9,*) '======= S.O.R. METHOD ======='
WRITE(9,210) last_i
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
    r=i+Nx*(j-1)+Nx*Ny*(k-1)

WRITE(9,150) i,j,k, x_Phi(r,0)

    END DO
  END DO
END DO
CLOSE(9)

END SELECT !end select case (method)


CASE(3)
!Print max val of vz @ entry

resfile=path_fold//'V/vzmax.txt'
WRITE(*,*) 'Writing the file ' ,resfile
OPEN(UNIT=9,FILE=resfile,POSITION="append",STATUS='unknown',ACTION='WRITE')
WRITE(9,'(A,F13.8,A,/)') '====== Max Vz @ k=1 (t=',time,') ======'

WRITE(9,135)  i_max, j_max,vz0max
     135 FORMAT (1x,'Vz_max(',I2,',',I2,', 1) = ',F18.9)
WRITE(9,'(A,//)') '-------------------------------------'     
CLOSE(9)

CASE(4)
!Print div(U) -> it must be of O(10^-7)
WRITE(divu(1:4),'(A)') 'divU'
WRITE(divU(5:10),'(I6.6)') it
resfile=path_fold//'U/'//divu//'.txt'
WRITE(*,*) 'Writing the file ',resfile
OPEN(UNIT=9,FILE=resfile,STATUS='unknown',ACTION='WRITE')
WRITE(9,'(/,5x,A,//)')'======= DIVERGENCE of U* (dimensionless)  ======='
WRITE(9,'(1x,A,I3,A,I3,A,6x,A,F9.5)') ' position of bottom entry = (',int(x_1/dx),',',int(y_1/dy),', 0 )',&
&'r_bottom = ',r_1
WRITE(9,'(1x,A,I3,A,I3,A,I3,A,6x,A,F9.5,//)') ' position of top entry    = (',int(x_2/dx),',',int(y_2/dy),',',Nz,')',&
&'r_top    = ',r_2

DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
     IF (abs(div_U(i,j,k)).GE.1e-3) THEN
      WRITE(9,'(1x,A,I2,A,I2,A,I2,A,F18.9,5x,A)') 'div_U(',i,',',j,',',k,') = ',div_U(i,j,k),' <=====(too high)==== '
     ELSE
      WRITE(9,222) i,j,k, div_U(i,j,k)
  222 FORMAT (1x,'div_U(',I2,',',I2,',',I2,') = ',F18.9)
     END IF
    END DO
  END DO
END DO

WRITE(9,'(//,A,F18.9)') 'Max(div_U) = ', maxdivU
CLOSE(9)
!write the binary for visualization pipeline
resfile=path_fold//'bin/'//divu//'.bin'   !OCCHIO
WRITE(*,*) 'Writing the file ',resfile
OPEN(UNIT=1,FILE=resfile,FORM='unformatted')
WRITE(1) div_U
CLOSE(1)



CASE(5)
WRITE(filename(1:2),'(A)') 'uf'
WRITE(filename(3:8),'(I6.6)') it
resfile=path_fold//'U/'//TRIM(filename)//'.txt'
WRITE(*,*) 'Writing the file ' ,resfile
OPEN(UNIT=9,FILE=resfile,STATUS='unknown',ACTION='WRITE')
WRITE(9,'(/,15x,A,//)')'======= U-FIELD  ======='
WRITE(9,'(1x,A,I3,A,I3,A,6x,A,F9.5)') ' approx position of bottom entry center = (',int(x_1/dx),',',int(y_1/dy),', 0 )',&
'r_bottom = ',r_1
WRITE(9,'(1x,A,I3,A,I3,A,I3,A,6x,A,F9.5,//)') ' approx position of top entry center    = (',int(x_2/dx),',',int(y_2/dy),',',Nz,')',&
'r_top    = ',r_2

DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
     
     
     WRITE(9,145) i,j,k, ux(I,j,k),i,j,k,uy(i,J,k),i,j,k,uz(i,j,K)
     145 FORMAT (1x,'Ux(',I2,',',I2,',',I2,') = ',F18.9,15x, &
     'Uy(',I2,',',I2,',',I2,') = ',F18.9,15x, &
     'Uz(',I2,',',I2,',',I2,') = ',F18.9)
     END DO
   END DO
END DO
CLOSE(9)




CASE(6)

!Print the Volumetric Flow Rate q_dot(1,1) = inflow vol rate 
                             !  q_dot(1,2) = outflow vol rate

resfile=path_fold//'U/qdot.dat'
WRITE(*,*) 'Writing the file ' ,resfile
OPEN(UNIT=9,FILE=resfile,STATUS='unknown',POSITION='append',ACTION='write')
num= INT(it/itout)
IF (num .eq. 1) THEN
WRITE(9,'(/,3x,A,///)') '========= Net Volume flow rate (-Q_inlet + Q_outlet) (dimensionless =========='
WRITE(9,'(A,F13.8,A,I6,A,F13.8,A,//)') 'Q_dot_net_wrt_q1 (t = ',time,'--> n_iter = ',it,' ) =',&
& abs(q_dot(1,2)-q_dot(1,1))/q_dot(1,1)*100._myk,' %'
WRITE(9,'(A,F13.8,A,F13.8,//,A)')'Q_dot(1) = ', q_dot(1,1), '  \\\\\\\  Q_dot(2) = ', q_dot(1,2),&
'--------------------------------------------------------------'

ELSE
WRITE(9,'(A,F13.8,A,I6,A,F13.8,A,//)')  'Q_dot_net_wrt_q1 (t = ',time,'--> n_iter = ',it,' ) =',&
& abs(q_dot(1,2)-q_dot(1,1))/q_dot(1,1)*100._myk,' %'
WRITE(9,'(A,F13.8,A,F13.8,//,A)') ' Q_dot(1) = ', q_dot(1,1),'  \\\\\\\  Q_dot(2) = ', q_dot(1,2),&
'--------------------------------------------------------------'
END IF
CLOSE(9)


CASE(7)
!Print the convective term
WRITE(filere(1:4),'(A)')'conv'
WRITE(filere(5:10),'(I6.6)') it
opfile1=path_fold//'bin/'//TRIM(filere)//'.bin'
opfile1=TRIM(opfile1)
!filename1 =path_fold//'U/conv/'//TRIM(filere)//'.txt'
!filename1 = TRIM(filename1)
write(*,*) 'Writing the file ',opfile1
!WRITE(*,*) 'Writing the file ',filename1

open(unit=1,file=opfile1,form='unformatted')
write(1) hx,hy,hz
CLOSE(1)
!OPEN(1,FILE=filename1,STATUS='unknown',ACTION='write')
!WRITE(1,'(/,15x,A,//)')'======= Convective part H (dimensionless) ======='
!WRITE(1,'(1x,A,I3,A,I3,A)') ' position of bottom entry = (',int(x_1/dx),',',int(y_1/dy),', 0 )'
!WRITE(1,'(1x,A,I3,A,I3,A,I3,A//)') ' position of top entry    = (',int(x_2/dx),',',int(y_2/dy),',',Nz,')'

!DO k=1,Nz
!  DO j=1,Ny
!    DO i=1,Nx
          
     !WRITE(1,141) i,j,k, hx(i,j,k),i,j,k,hy(i,j,k),i,j,k,hz(i,j,k)
     !141 FORMAT (1x,'Hx(',I2,',',I2,',',I2,') = ',F18.9,1x,15x, &
     !'Hy(',I2,',',I2,',',I2,') = ',F18.9,15x, &
     !'Hz(',I2,',',I2,',',I2,') = ',F18.9)
     !END DO
   !END DO
!END DO
!CLOSE(1)
END SELECT
RETURN
END SUBROUTINE stampa
