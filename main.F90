program main
use mod_fluid
!______________________________________________________________________
!     Dichiarazione variabili locali diretta e globali tramite
!     include/common
!_____
!      real t1,t2,dtime,dummy,global_timef,t3,t4
!      real int_Tau_rz,pi,dt_rk,U_bulk
!      real Wi,r_sp_0,r_sp_max
!      integer part_flag
!____.__________________________________________________________________
!____.
implicit none

!local variables character
CHARACTER(2),PARAMETER :: iterazione='it', delta_t='dt'
CHARACTER(6),PARAMETER :: iteraz_step='i_step'
CHARACTER(1),PARAMETER :: time='t'
LOGICAL :: ext

WRITE(*,'(/,2x,A,/)') ' (Optional) Insert a comment and press ENTER (max 100 char) : '
READ(*,'(A)') comment

WRITE(*,'(////,1x,A,//)')                       '=========================================='
WRITE(*,'(4x,A,11x,A,8x,A)')                     '___',        '|'     , '___'
WRITE(*,'(3x,A,6x,A,3x,A,3x,A,4x,A)')           '|'   ,'|','|','|___'  ,'|'
WRITE(*,'(3x,A,6x,A,3x,A,3x,A,3x,A)')           '|'   ,'|','|','|   |' ,'|--' 
WRITE(*,'(3x,A,3x,A,3x,A,3x,A,3x,A,//)')        '|___','|___|','|___|' ,'|___', 'FLU1D'
WRITE(*,'(1x,A,//)')                            '=========================================='
!Subroutines Calling

call input()
WRITE(*,'(A,//)') '(01) check [input]'
!Unknown terms array dimension
Nt=Nx*Ny*Nz

!!Geometrical definition of the  openings

!Bottom opening Specs
 x_1=0.25_myk*Xmax
 y_1=Ymax*0.5_myk
 r_1= 0.1_myk
 sigma_u = 0.03_myk
!Top opening Specs
 y_2=y_1
 x_2=0.6_myk*Xmax
 r_2=r_1!*0.5_myk

WRITE(*,'(A,/)') '===================================================='
WRITE(*,'(6x,A,//)')   ' Writing the openings specs : '
WRITE(*,'(2x,A,6x,A,/)') ' (1) Bottom Opening',' (2) Top   Opening'
WRITE(*,'(6x,A,F6.4,5x,A,5x,A,F6.4)')   ' x_1 = ',x_1,'||',' x_2 = ', x_2
WRITE(*,'(6x,A,F6.4,5x,A,5x,A,F6.4)')   ' y_1 = ',y_1,'||',' y_2 = ', y_2
WRITE(*,'(6x,A,F6.4,5x,A,5x,A,F6.4,/)') ' r_1 = ',r_1,'||',' r_2 = ', r_2
WRITE(*,'(A,/)') '===================================================='

call allocate_var()
WRITE(*,'(1x,A,/)') '(02) check [allocate_var]'
WRITE(*,'(A,//)') ' ==============================================='
call griglia()
WRITE(*,'(A,//)') ' ==============================================='
WRITE(*,'(1x,A,/)') '(03) check [griglia]'

CALL sim_spec 

call init_rkcoef()
WRITE(*,'(/,1x,A)') '(04) check [init_rkcoef]'
call init_field()
WRITE(*,*) '(05) check [init_field]'
CALL buildmatrix()
WRITE(*,*) '(06) check [buildmatrix]'
call stampa(1)
WRITE(*,'(1x,A,//)') '(07) check [stampa(1)]'

!10 WRITE(*,'(A,/,A,/,A,/,A)')' insert the method that should be used :',' [1] for Jacobi',' [2] for Gauss-Seidel',' [3] for S.O.R.'
!READ(*,*)  method
!290 FORMAT (' the method was chosen is ',I1,1x, '- ',A,' method')

!IF (method==1) THEN
!WRITE(*,290) method, jac
!ELSE IF (method==2) THEN
!WRITE(*,290) method, gs

!ELSE IF (method==3) THEN
!WRITE(*,290) method, sor

!ELSE 
!WRITE(*,*) 'ERROR ~ The value entered must be chosen between the values 1(Jacobi) - 2(G-S) - 3(S.O.R)'
!GO TO 10
!END IF

!Selector : S.O.R. method (forced)
method=3
!Define first guess to solution (initializing iterative method array)
INQUIRE(FILE=path_fold//'iter_sol.bin',EXIST=ext) !If the sim is restarting, it allows to use the last solution found,
IF (ext) THEN
 OPEN(UNIT=1,FILE=path_fold//'iter_sol.bin',FORM='unformatted',STATUS='old',ACTION='read')
 READ(1) X_Phi(:,0)
 CLOSE(1)
 WRITE(*,'(1x,A,/)') ' - initializing the first guess solution using the last found in the previous simulation'
ELSE                                              !otherwise initializing the first guess as zero array.
 x_Phi(:,0) = 0.0_myk
 WRITE(*,'(1x,A,/)') ' - initializing the first guess solution using the 0-array'
END IF
!Simulation CPU time(start)
call cpu_time(start)
!____.__________________________________________________________________
!     Inizio ciclo temporale
!____
WRITE(*,'(A,//,3x,A,//,A,//)') '================================',' Time Cycle Starting ','================================'

write(*,*)'------------------------------------------------------------'
!write(*,*)' it    i_step              t                  dt'
write(*,23)iterazione, iteraz_step, time, delta_t
23 FORMAT (4x,A,7x,A,12x,A,20x,A)
write(*,*)'------------------------------------------------------------'


time_iter: do it=itmin+1,itmax
 IF(it>itmin+1) THEN
  IF(Co<1) THEN
   WRITE(*,*)'-------------------------------------------------------------------------'
   WRITE(*,'(1x,A,I5,A,F13.7,A)')'|~ Courant N° of the previous time step (',it-1,') = ',Co,' (<-OK) ~|'
   WRITE(*,*)'-------------------------------------------------------------------------'

  ELSE
   WRITE(*,*)'----------------------------------------------------------------------------------'
   WRITE(*,'(1x,A,I5,A,F13.7,A)')'|~ Courant N° of the previous time step (',it-1,') = ',Co,' (<-TOO HIGH!!!) ~|'
   WRITE(*,*)'----------------------------------------------------------------------------------'
  END IF
 END IF
!Runge-Kutta low storage 3/4 steps
    rkutta_iter: do i_step=1,n_step
        t=t+dt*(ark(i_step,i_kutta)+brk(i_step,i_kutta))
        write(*,*)'------------------------------------------------------------'
        write(*,'(I7,7x,I1,8x,f13.7,8x,f13.7)')it,i_step,t,dt
        write(*,*)'------------------------------------------------------------'

         call prhs()
         WRITE(*,*)'(08) check [prhs]'
         call rhs()
         WRITE(*,*) '(09) check [rhs]'
         call linear()
         WRITE(*,*) '(10) check [linear]'
    END DO rkutta_iter

    CALL boundary_cond(1) !B.C. on the auxiliary vel. field (V)
    WRITE(*,*) '(11) check [boundary_cond(1)] -> releted to V-field'
    CALL vzmaxbound
    WRITE(*,*) '(12) check [vzmaxbound]'
    CALL diver_V
    WRITE(*,*) '(13) check [diver_V]'
    CALL known_term
    WRITE(*,*) '(14) check [known_term]'
    CALL Phi_iter(method)
    WRITE(*,*) '(15) check [Phi_iter]'
    CALL associate_var
    WRITE(*,*) '(16) check [associate_var]'
    CALL boundary_cond(3) !B.C. on the irrotational vel. field (grad(Phi))
    WRITE(*,*) '(17) check [boundary_cond(3)] -> related to Phi-field'
    CALL grad_Phi()
    WRITE(*,*) '(18) check [grad_Phi]'
    CALL U_field()
    WRITE(*,*) '(19) check [U_field]'
    CALL vol_flowrate
    WRITE(*,*) '(20) check [vol_flowrate]'
    CALL diver_U
    WRITE(*,*) '(21) check [diver_U]'
    if (mod(it,itout)==0) then
     CALL vol_flowrate
     WRITE(*,*) '(22) check [vol_flowrate]'
     CALL stampa(1)  
     CALL stampa(2)
     CALL stampa(4)
     !CALL stampa(5)
     CALL stampa(6)
     open(unit=1,file='restart.dat',status="unknown")
     write(1,*) it
     WRITE(1,*) t
     close(1)
    endif
    CALL courant
   enddo time_iter
   CALL courant
   IF(Co<1) THEN
   WRITE(*,*)'-------------------------------------------------------------------------'
   WRITE(*,'(1x,A,I5,A,F13.7,A)')'|~ Courant N° of the previous time step (',it-1,') = ',Co,' (<-OK) ~|'
   WRITE(*,*)'-------------------------------------------------------------------------'

  ELSE
   WRITE(*,*)'----------------------------------------------------------------------------------'
   WRITE(*,'(1x,A,I5,A,F13.7,A)')'|~ Courant N° of the previous time step (',it-1,') = ',Co,' (<-TOO HIGH!!!) ~|'
   WRITE(*,*)'----------------------------------------------------------------------------------'
  END IF

 CALL clean
 WRITE(*,'(1x,A,//)') '(23) check [clean]'
WRITE(*,*)           ' |---------------------------------------------------------| ' 
WRITE(*,*)           ' |----------> SIMULATION SUCCESSFULLY COMPLETED <----------| '
WRITE(*,'(1x,A,//)') ' |---------------------------------------------------------| ' 
!Simulation CPU time(finish)
call cpu_time(finish)
OPEN(UNIT=1,FILE=path_fold//'sim_spec.txt',STATUS='old',POSITION="APPEND",ACTION="WRITE")
WRITE(1,'(/,A,F13.6,A)') ' CPU elapsed time = ', finish-start,' [sec]'
WRITE(1,'(/,A,//)') ' ======================================================== '
CLOSE(1)
!!!!DEBUG!!! 
!!!!DEBUG!!!        call prhs_pipe(i_kutta,n)
!!!!DEBUG!!!        call rhs_pipe()
!!!!DEBUG!!!        call linear_pipe(i_kutta,n)
!!!!DEBUG!!! 
!!!!DEBUG!!! !____.__________________________________________________________________
!!!!DEBUG!!! !     Serve a tagliare le componenti di velocita' ad alta frequenza vicino
!!!!DEBUG!!! !     l''asse per problemi di stabilita'
!!!!DEBUG!!! !_____
!!!!DEBUG!!! 
!!!!DEBUG!!!       if(Nt.gt.1)  call taglia_alto_pipe()
!!!!DEBUG!!! 
!!!!DEBUG!!!       call cornici_pipe()
!!!!DEBUG!!! 
!!!!DEBUG!!! !_____
!!!!DEBUG!!!       call Poisson_pipe()
!!!!DEBUG!!! !____.__________________________________________________________________
!!!!DEBUG!!!       call lapl_u()
!!!!DEBUG!!! !____.__________________________________________________________________
!!!!DEBUG!!!       print*, 'controllo ut',rho_ut(nt/2,Nr+1,nz/2),rho_ut(nt/2,Nr,nz/2)
!!!!DEBUG!!!      & ,0.5*(rho_ut(nt/2,Nr+1,nz/2)+rho_ut(nt/2,Nr,nz/2))
!!!!DEBUG!!! !     Aggiornamento dt e t
!!!!DEBUG!!!           t=t+dt*(Ark(n,i_kutta)+Brk(n,i_kutta))
!!!!DEBUG!!!           write(*,*)'+----------------------------------------+'
!!!!DEBUG!!!           write(*,*)'|',it/n_step,t,dt
!!!!DEBUG!!!           write(*,*)'+----------------------------------------+'
!!!!DEBUG!!! !____.__________________________________________________________________
!!!!DEBUG!!! !     Stampa i campi se la soluzione �fisica e l''iterata �quella di
!!!!DEBUG!!! !     stampa.
!!!!DEBUG!!! !_____
!!!!DEBUG!!!        if(mod(it,itout).eq.0) then
!!!!DEBUG!!!          call stampa_pipe()
!!!!DEBUG!!!          write(*,*)'stampa'
!!!!DEBUG!!! !____.
!!!!DEBUG!!!          open(unit=2,file='DAT/restart.dat',status="unknown")
!!!!DEBUG!!!          write(2,*)it/n_step
!!!!DEBUG!!!          close(2)
!!!!DEBUG!!!       endif
!!!!DEBUG!!! !____.
!!!!DEBUG!!!       !t2=global_timef()
!!!!DEBUG!!!       if((n.eq.n_step).and.((t2-t1).gt.28200000)) then
!!!!DEBUG!!!          call stampa_pipe()
!!!!DEBUG!!!          write(*,*)'stampa di restart: ',t2-t1
!!!!DEBUG!!!          open(unit=2,file='DAT/restart.dat',status="unknown")
!!!!DEBUG!!!          write(2,*)it/n_step
!!!!DEBUG!!!          close(2)
!!!!DEBUG!!! 	 stop
!!!!DEBUG!!!       endif
!!!!DEBUG!!! 
!!!!DEBUG!!! 
!!!!DEBUG!!!       t2=global_timef()
!!!!DEBUG!!!       write(*,*)'tempo tot',t2-t1
!!!!DEBUG!!! !____.__________________________________________________________________
!!!!DEBUG!!! !     Fine iterata
!!!!DEBUG!!! !_____
!____.__________________________________________________________________
end

!!!!____.__________________________________________________________________
!!!      function global_timef()
!!!      real global_timef
!!!      real t1
!!!      integer tiarray(8)
!!!      call date_and_time(VALUES=tiarray)
!!!      !write(*,*)'int',tiarray(6),tiarray(7),tiarray(8)
!!!      t1=float(tiarray(8))*0.001+float(tiarray(7))+float(tiarray(6))*60.
!!!      global_timef=t1
!!!      return
!!!      end
!!!      function omp_get_thread_num()
!!!      integer t1,omp_get_thread_num
!!!      t1=0
!!!      omp_get_thread_num=t1
!!!      return
!!!      end
!!!
