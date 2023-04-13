SUBROUTINE sim_spec
!This subroutine create the directory and sons of the actual simulation and writes a file
!sim_spec.txt that shows the simulation specs
USE mod_fluid
IMPLICIT NONE
INTEGER :: valori(1:8)
INTEGER :: num_res


OPEN(1,FILE='restart.dat',STATUS='OLD',ACTION='READ')
READ(1,*) num_res
CLOSE(1)


OPEN(UNIT=1,FILE="sim.dat",STATUS='OLD',ACTION='READ')
READ(1,*) num_sim
CLOSE(1)

CALL DATE_AND_TIME (VALUES=valori)
IF (num_res.EQ.0) THEN
 !predispongo folders per nuova simulazione
 num_sim = num_sim + 1

 OPEN(UNIT=1,FILE='sim.dat',STATUS='REPLACE',ACTION="WRITE")
 WRITE(1,'(I4)') num_sim
 CLOSE(1)

 fold_num(1:3) ='sim'
 WRITE(fold_num(4:7),'(I4.4)') num_sim
 path_fold = './DATA/'//fold_num//'/'
 path_fold = TRIM(path_fold) !'./DATA/sim000#/


 CALL SYSTEM ('mkdir -p '//path_fold//'temperature/')
 CALL SYSTEM ('mkdir -p '//path_fold//'bin/; mkdir -p '//path_fold//'U/conv/')
 CALL SYSTEM ('mkdir -p '//path_fold//'V/')
 CALL SYSTEM ('mkdir -p '//path_fold//'Phi/gs; mkdir '//path_fold//'Phi/jacobi;mkdir '&
 &//path_fold//'Phi/sor')

 OPEN(UNIT=1,FILE=path_fold//'sim_spec.txt',STATUS='new',ACTION="WRITE")
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(/,A,I4,A,I2.2,A,I2.2,/,A,13x,A,I2.2,A,I2.2,A,I2.2,/)') ' DATE (yyyy/mm/dd) = ', valori(1),'/', valori(2),'/', valori(3),&
 & ' TIME', ' = ', valori(5),':',valori(6),':',valori(7)
 WRITE(1,'(A,A,/)') ' Comment : ',comment
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(A,/)') ' X-Y-Z divisions / spatial discretization steps'
 WRITE(1,'(A,I3,A,E13.8)') ' [Nx] = ',Nx,' -----> dx = ',dx
 WRITE(1,'(A,I3,A,E13.8)') ' [Ny] = ',Ny,' -----> dy = ',dy
 WRITE(1,'(A,I3,A,E13.8,//)') ' [Nz] = ',Nz,' -----> dz = ',dz
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(A)') ' Time Specs '
 WRITE(1,'(/,4x,A,F13.6,/,A,I6,/,A,I6,/)') ' dt = ', dt, ' it_max = ',itmax,' it_out = ',itout
 CLOSE(1)

ELSE
 !se non ho resettato (num_res.NE.0), aggiorno il file sim.dat con la continuazione della simulazione
 fold_num(1:3) ='sim'
 WRITE(fold_num(4:7),'(I4.4)') num_sim
 path_fold = './DATA/'//fold_num//'/'
 path_fold = TRIM(path_fold) !'./DATA/sim000#/

 OPEN(UNIT=1,FILE=path_fold//'sim_spec.txt',STATUS='old',POSITION='APPEND',ACTION="WRITE")
 WRITE(1,'(A,/)') ' CONTINUATION OF THE PREVIOUS SIMULATION '
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(/,A,I4,A,I2.2,A,I2.2,/,A,13x,A,I2.2,A,I2.2,A,I2.2,/)') ' DATE (yyyy/mm/dd) = ', valori(1),'/', valori(2),'/', valori(3),&
 & ' TIME',' = ', valori(5),':',valori(6),':',valori(7)
 WRITE(1,'(A,A,/)') ' Comment : ',comment
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(A,/)') ' X-Y-Z divisions / spatial discretization steps'
 WRITE(1,'(A,I3,A,E13.8)') ' [Nx] = ',Nx,' -----> dx = ',dx
 WRITE(1,'(A,I3,A,E13.8)') ' [Ny] = ',Ny,' -----> dy = ',dy
 WRITE(1,'(A,I3,A,E13.8,//)') ' [Nz] = ',Nz,' -----> dz = ',dz 
 WRITE(1,'(/,A,//)') ' ======================================================== '
 WRITE(1,'(A)') ' Time Specs '
 WRITE(1,'(/,4x,A,F13.6,/,A,I6,/,A,I6,/)') ' dt = ', dt, ' it_max = ',itmax,' it_out = ',itout
 CLOSE(1)
END IF
RETURN
END SUBROUTINE sim_spec
