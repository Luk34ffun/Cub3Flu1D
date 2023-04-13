      subroutine input()
      use mod_fluid
!______________________________________________________________________
      implicit none
!______________________________________________________________________
      open(unit=1,file='input.dat',status="unknown")
          read(1,*) ver
          read(1,*) fileres
          read(1,*) !salto path1
          read(1,*) !salto foldername
          read(1,*) n_step
          READ(1,*) !salto itmin_postproc
          read(1,*) itmin
          read(1,*) itmax
          read(1,*) itout
          read(1,*) dt
          read(1,*) Xmax
          read(1,*) Ymax
          read(1,*) Zmax
          read(1,*) Nx
          read(1,*) Ny
          read(1,*) Nz
          !input adimensional numbers
          read(1,*) Re
          READ(1,*) Pr
          READ(1,*) Ri
          !input variables for iteration method
          READ(1,*) max_iter
          READ(1,*) max_err
          !for S.O.R. method
          READ(1,*) omg
!--------------------------------------------------!
          !Buoyancy source temperature (dimensional)
           READ(1,*) Th_s
      close(1)
!______________________________________________________________________
      if(itmin.eq.999999) then
        open(unit=1,file='restart.dat',status="unknown")
        read(1,*) itmin !itmin passa da 999999 al valore letto dentro restart.dat
        READ(1,*) t
        close(1)
      endif
!______________________________________________________________________
      itmin=itmin
      itmax=itmin+itmax
      itout=itout

      if (n_step.eq.4) then
          i_kutta=2
      else
      n_step=3
      i_kutta=1
      endif
!______________________________________________________________________
      Inv_Re=1./Re
      Inv_RePr = Inv_Re/Pr
      end

