subroutine init_field
use mod_fluid
implicit none
logical:: ext
real(myk)::r2,sigma
integer::i,j,k

mu = 1._myk

IF (itmin.NE.0) THEN  !se (itmin.NE.0) voglio proseguire la prec simul

 OPEN(UNIT=1,FILE="sim.dat",STATUS='OLD',ACTION='READ')
 READ(1,*) num_sim
 CLOSE(1)

 fold_num(1:3) ='sim'
 WRITE(fold_num(4:7),'(I4.4)') num_sim
 path_fold = './DATA/'//fold_num//'/'
 path_fold = TRIM(path_fold) !'./DATA/sim000#/


write(fileres(1:1),'(A)') ver
write(fileres(4:9),'(I6.6)') itmin
opfile = path_fold//'bin/'//trim(fileres)

WRITE(fileth(1:2),'(A)') 'th'
WRITE(fileth(3:8),'(I6.6)') itmin
opfileth=path_fold//'bin/'//fileth//'.bin'
opfileth=TRIM(opfile)



!potrei toglierlo perchè superfluo
inquire(FILE=opfile,exist=ext) 
!hp1: ext=true (esiste il file opfile)=> CALL leggi()
!hp2: ext=false => definisco I.C. da comandi (1)
  IF (.NOT.ext) THEN 
   99 CONTINUE
   WRITE(*,'(1x,A)') ' - The file related to the I.C. was not found '
   WRITE(*,'(1x,A)') ' - I.C. will be initialized from command '
   !(1) Initial field definition

   !I.C. velocity field
   DO k=1,Nz
    DO j=1,Ny
     DO i=1,Nx
      r2 = (0.5_myk*(x(i)+x(i+1))-0.5_myk)**2+(y(j)-0.5_myk)**2+(z(k)-0.5)**2
      sigma = 0.2_myk
     !ux(i,j,k)=(1.0_myk*exp(-r2/(2*(sigma**2))))
    
     !I.C. temperature field
     !th(i,j,k) = 1.0_myk      !(273.15_myk + 20.0_myk)/Th_0
 
     END DO
    END DO
   END DO
    ux=0.0_myk
    uy=0.0_myk
    uz=0.0_myk
    th = 1.0_myk
   
  ELSE
   call leggi()
  END IF
  GOTO 100
 ELSE !(itmin.EQ.0)
 GOTO 99
 END IF
100 CONTINUE 
CALL boundary_cond(1)
CALL boundary_cond(2)
WRITE(*,'(1x,A)') ' - From [init_field] CALL [boundary_cond(1-2)]-> related to U-field + theta-field '
RETURN
END SUBROUTINE init_field
