PROGRAM reset
!This program allows to clean up (eventually) the folders where the previous simulation data are stored and set to 0 the value in restart.dat
!to perform a new simulation
IMPLICIT NONE
CHARACTER(1) :: choice
CHARACTER(8) :: fold_num ='sim4567/'
CHARACTER(30) :: path_fold
INTEGER :: num_sim


OPEN(UNIT=2,FILE='restart.dat',STATUS='replace', ACTION='WRITE')
WRITE(2,*) '0'
WRITE(2,*) '0'
CLOSE(2)

OPEN(UNIT=1,FILE="sim.dat",STATUS='OLD',ACTION='READ')
READ(1,*) num_sim
CLOSE(1)

WRITE(*,'(A,/,20x,A)') ' Would you also to delete the previous simulation data ? ',' ### Type y/n ### '
READ(*,*) choice
IF (choice == 'y') THEN 
 IF (num_sim==0) THEN
  WRITE(fold_num(4:7),'(I4.4)') num_sim+1
  WRITE(*,'(A,/,A,A,A)') ' ==== A new Simulation can now be performed ==== ',' ==== Target Folder : ',fold_num,' (empty) ===='
 ELSE !(num_sim>0)
  WRITE(fold_num(4:7),'(I4.4)') num_sim
  path_fold = './DATA/'//fold_num
  path_fold = TRIM(path_fold) !'./DATA/sim000#/
  CALL SYSTEM ('rm -r '//path_fold)
  WRITE(*,'(A,A,A)') ' The folder ',fold_num,' content has been deleted' 
  WRITE(*,'(A,/,A,A,A)') ' ==== A new Simulation can now be performed ==== ',' ==== Target Folder : ',fold_num,' (empty) ===='
  OPEN(UNIT=1,FILE='sim.dat',STATUS='REPLACE',ACTION="WRITE")
  WRITE(1,*) num_sim-1
  CLOSE(1)
 END IF
ELSE !(if choice =/ y)
 IF (num_sim==0) THEN
  WRITE(fold_num(4:7),'(I4.4)') num_sim+1
  WRITE(*,'(A,/,A,A,A)') ' ==== A new Simulation can now be performed ==== ',' ==== Target Folder : ',fold_num,' (empty) ===='
 ELSE
  WRITE(fold_num(4:7),'(I4.4)') num_sim
  path_fold = './DATA/'//fold_num
  path_fold = TRIM(path_fold) !'./DATA/sim000#/
  WRITE(*,'(A)') ' Previous simulation folder '//fold_num//' has not been deleted'
  num_sim = num_sim+1
  WRITE(fold_num(4:7),'(I4.4)') num_sim
  WRITE(*,'(A,A)') ' The current Target folder is :',fold_num
 END IF
END IF

END PROGRAM reset
