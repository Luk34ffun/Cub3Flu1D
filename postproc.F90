module mod_fluid

implicit none
!space
integer,parameter:: myk=kind(0.0D0)
integer:: Nx,Ny,Nz
real(myk),dimension(:),allocatable:: X,Y,Z
real(myk):: Xmax,Ymax,Zmax
real(myk):: dx,dy,dz

!time
integer:: itmin,itmax,itout,it,foldnu
real(myk):: dt,t

!file names
character(1):: ver
character(128)::fileres,fileres1,fileres2,fileres3
character(128)::path1
CHARACTER(8) :: foldername
character(256)::opfile,opfileh,opfileth,opfilediv

!variables
real(myk),dimension(:,:,:)  ,allocatable:: rho_ux,rho_uy,rho_uz,th,divU
real(myk),dimension(:,:,:)  ,allocatable:: hx,hy,hz
real(myk),dimension(:,:,:,:),allocatable:: grid

!physical parameters
real(myk):: Re,Inv_Re

contains
!______________________________________________________________________
!______________________________________________________________________
!______________________________________________________________________

!______________________________________________________________________
!______________________________________________________________________
subroutine input()
implicit none
integer:: n_step
open(unit=1,file='../../test/input.dat',status="old",ACTION='read')
    read(1,*) ver
    read(1,*) fileres
    read(1,*) path1
    read(1,*) foldername
    read(1,*) n_step
    read(1,*) itmin
    read(1,*) !salto itmin_prog_main =999999
    read(1,*) !itmax
    read(1,*) itout
    read(1,*) dt
    read(1,*) Xmax
    read(1,*) Ymax
    read(1,*) Zmax
    read(1,*) Nx
    read(1,*) Ny
    read(1,*) Nz
    read(1,*) Re
close(1)
Inv_Re=1./Re
OPEN(UNIT=1,FILE='../../test/sim.dat',STATUS='old',ACTION='read')
READ(1,'(I4.4)') foldnu
CLOSE(1)
OPEN(UNIT=1,FILE='../../test/restart.dat',STATUS='old',ACTION='read')
READ(1,*) itmax
WRITE(*,*) itmax
CLOSE(1)
endsubroutine input
!______________________________________________________________________
!______________________________________________________________________
subroutine allocate_var
allocate(X(0:Nx+1))
allocate(Y(0:Ny+1))
allocate(Z(0:Nz+1))
allocate(grid(3,1:Nx,1:Ny,1:Nz))
allocate(rho_ux(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(rho_uy(0:Nx+1,0:Ny+1,0:Nz+1))
allocate(rho_uz(0:Nx+1,0:Ny+1,0:Nz+1))

ALLOCATE(hx(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(hy(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(hz(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(th(0:Nx+1,0:Ny+1,0:Nz+1))
ALLOCATE(divU(0:Nx,0:Ny,0:Nz))
return
endsubroutine allocate_var
!______________________________________________________________________
!______________________________________________________________________
subroutine griglia()
implicit none
integer:: i,j,k
dx = Xmax/float(Nx)
X(0) = -dx/2.0_myk
do i=1,Nx+1
   X(i) = X(i-1) + dx
enddo
dy = Ymax/float(Ny)
Y(0) = -dy/2.0_myk
do i=1,Ny+1
   Y(i) = Y(i-1) + dy
enddo
dz = Zmax/float(Nz)
Z(0) = -dz/2.0_myk
do i=1,Nz+1
   Z(i) = Z(i-1) + dz
enddo
write(*,*) 'dx =',dx,' Nx =',Nx
write(*,*) 'dy =',dy,' Ny =',Ny
write(*,*) 'dz =',dz,' Nz =',Nz
return
end
!______________________________________________________________________
!______________________________________________________________________
SUBROUTINE leggi()
IMPLICIT NONE
INTEGER :: ios

WRITE(foldername(4:7),'(I4.4)') foldnu
!leggo campo di velocità U
write(fileres(1:1),'(A)') ver
write(fileres(4:9),'(I6.6)') it
write(fileres(11:13),'(A)') 'bin'
opfile = TRIM(path1)//foldername//'bin/'//trim(fileres)
opfile = TRIM(opfile)
write(*,'(1x,A,/)') '~ reading the .bin file '//opfile(1:45)
open(unit=1,file=opfile,form='unformatted')
read(1) rho_ux,rho_uy,rho_uz
close(1)
!leggo parte convettiva
WRITE(fileres1(1:4),'(A)')'conv'
WRITE(fileres1(5:10),'(I6.6)') it
WRITE(fileres1(11:14),'(A)')'.bin'
opfileh=TRIM(path1)//foldername//'bin/'//TRIM(fileres1)
opfileh=TRIM(opfileh) !-------------------------------------(*)
WRITE(*,*) '~ reading the .bin file ',opfileh
OPEN(UNIT=1,FILE=opfileh,FORM='unformatted',ACTION='read')
READ(1,IOSTAT=ios) hx,hy,hz
CLOSE(1)
!leggo temperatura
WRITE(fileres2(1:2),'(A)') 'th'
WRITE(fileres2(3:8),'(I6.6)') it
WRITE(fileres2(9:12),'(A)') '.bin'
opfileth=TRIM(path1)//foldername//'bin/'//fileres2
opfileth=TRIM(opfileth)
WRITE(*,*)'~ reading the .bin file ',opfileth
OPEN(UNIT=1,FILE=opfileth,FORM='unformatted',ACTION='READ')
READ(1,IOSTAT=ios) th
CLOSE(1)
!leggo divergenza di U
WRITE(fileres3(1:4),'(A)') 'divU'
WRITE(fileres3(5:10),'(I6.6)') it
WRITE(fileres3(11:14),'(A)') '.bin'
opfilediv=TRIM(path1)//foldername//'bin/'//fileres3
opfilediv=TRIM(opfilediv)
WRITE(*,*)'~ reading the .bin file ',opfilediv
OPEN(UNIT=1,FILE=opfilediv,FORM='unformatted',ACTION='READ')
READ(1,IOSTAT=ios) divU
CLOSE(1)
RETURN
ENDSUBROUTINE leggi
!______________________________________________________________________
!______________________________________________________________________
SUBROUTINE copy
IMPLICIT NONE
CALL SYSTEM ('mkdir '//TRIM(path1)//foldername//'vtk/ ; cp -r ./vtk/* '//TRIM(path1)//foldername//'vtk/ ; rm ./vtk/*')
WRITE(*,*) ' === vtk files have been copied to the folder : ',TRIM(path1)//foldername//'vtk ==='
RETURN
END SUBROUTINE copy
!______________________________________________________________________
!_____________________________________________________________________
SUBROUTINE erase
IMPLICIT NONE
CHARACTER(1) :: chr
WRITE(*,'(A,/,5x,A)') ' would you to delete the bin DATA ? ',' === type y/n ==='
READ(*,'(A)') chr
IF (chr =='y') THEN
 CALL SYSTEM ('rm '//TRIM(path1)//foldername//'bin/*.bin')
 WRITE(*,*) ' === bin files have been deleted from the folder : ',TRIM(path1)//foldername//'bin ==='
ELSE
 WRITE(*,*) ' === bin files have been left in the folder : ',TRIM(path1)//foldername//'bin ==='
END IF
RETURN
END SUBROUTINE erase
!______________________________________________________________________
!_______________________________________________________________________
end module mod_fluid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module paraview_module
use mod_fluid
implicit none
character(1)::c
contains
!______________________________________________________________________
!______________________________________________________________________
!______________________________________________________________________
!______________________________________________________________________
SUBROUTINE stampa_vtk()
IMPLICIT NONE
REAL(myk):: time
WRITE(fileres(1:1),'(A)') ver
WRITE(fileres(4:9),'(I6.6)') it
WRITE(fileres(11:13),'(A)') 'vtk'
opfile  = './vtk/'//trim(fileres)

IF (it.NE.itmax) THEN
 WRITE(*,'(A,A,//)') ' ->  the .vtk file has been printed : ',opfile(1:20)
ELSE
 WRITE(*,'(A,A,///,7x,A,//)') ' ->  the .vtk file has been printed : ',opfile(1:20),&
 &' ==== .vtk files are ready to be opened with ParaView ==== '
END IF
!open(unit=1,file=opfile,form='unformatted',access='stream',convert='big_endian')
OPEN(UNIT=1,FILE=opfile,FORM='FORMATTED',ACTION="WRITE")
time = it*dt
call write_head(1,time)
call write_grid(1)
call write_velocity(1)
!CALL write_convective(1)
!CALL write_temperature(1)
!CALL write_divU(1)
close(1)

return
ENDSUBROUTINE stampa_vtk
!______________________________________________________________________
!______________________________________________________________________
SUBROUTINE write_head(nfile,time)
IMPLICIT NONE
INTEGER nfile
REAL(myk)time
CHARACTER(10):: time_strg


WRITE(nfile,'(A)') "# vtk DataFile Version 3.0"
WRITE(nfile,'(A)') "Fluid t= 0.0"
WRITE(nfile,'(A)') "ASCII"

RETURN
ENDSUBROUTINE write_head
!______________________________________________________________________
!______________________________________________________________________
SUBROUTINE write_grid(nfile)
IMPLICIT NONE
INTEGER i,j,k,nfile,ii
REAL xx,yy,zz,sx,cx
INTEGER:: appo_nx,appo_ny,appo_nz,appo_tot
CHARACTER(128):: dimens_string, size_string, tot_dimens_string, point_data_string


appo_nx = Nx; appo_ny = Ny; appo_Nz=Nz
appo_tot = appo_nx*appo_ny*appo_nz

dimens_string = 'DIMENSIONS ' 
WRITE(size_string,'(I5)') appo_nx
dimens_string = TRIM(dimens_string)//TRIM(size_string)
WRITE(size_string,'(I5)') appo_ny
dimens_string = TRIM(dimens_string)//TRIM(size_string)
WRITE(size_string,'(I5)') appo_nz
dimens_string = TRIM(dimens_string)//TRIM(size_string)

tot_dimens_string = 'POINTS ' 
WRITE(size_string,'(I13)') appo_tot
tot_dimens_string = TRIM(tot_dimens_string)//TRIM(size_string)//' double'

point_data_string = 'POINT_DATA '//TRIM(size_string)


!Type=STRUCTURED_POINTS,STRUCTURED_GRID,UNSTRUCTURED_GRID,POLYDATA,RECTILINEAR_GRID,FIELD

WRITE(nfile,'(A)') "DATASET STRUCTURED_GRID"
WRITE(nfile,'(A)') TRIM(dimens_string)
WRITE(nfile,'(A)') TRIM(tot_dimens_string)


DO k=1,Nz   
 DO j=1,Ny
  DO i=1,Nx
    WRITE(nfile,15) X(i),Y(j),Z(k)
  END DO
 END DO
END DO

write(nfile,'(A)') trim(point_data_string)

15    format(12(E13.6,1x))
return
endsubroutine write_grid
!______________________________________________________________________
!______________________________________________________________________
subroutine write_velocity(nfile)
implicit none
integer i,j,k,nfile
real rr,r_u,r_v,r_w,ra,time
real r_r,r_t,sx,cx,xx,yy,zz


write(nfile,'(A)') "VECTORS U double"

!!! parte centrale
do k=1,Nz
 do j=1,Ny
  do i=1,Nx
   r_t=0.5*(rho_ux(I,j,k)+rho_ux(I-1,j,k))
   r_r=0.5*(rho_uy(i,J,k)+rho_uy(i,J-1,k))
   r_w=0.5*(rho_uz(i,j,K)+rho_uz(i,j,K-1))
   write(nfile,15)r_t,r_r,r_w
  enddo
 enddo
enddo

15    format(12(E13.6,1x))
return
endsubroutine write_velocity
!______________________________________________________________________
!______________________________________________________________________

subroutine write_convective(nfile)
implicit none
integer i,j,k,nfile
real rr,r_u,r_v,r_w,ra,time
real r_x,r_y,r_z,sx,cx,xx,yy,zz

!open(nfile,file=opfile,form='formatted',position='append')

write(nfile,'(A)') "VECTORS H double"

!!! parte centrale
do k=1,Nz
 do j=1,Ny
  do i=1,Nx
   r_x=0.5*(hx(I,j,k)+hx(I-1,j,k))
   r_y=0.5*(hy(i,J,k)+hy(i,J-1,k))
   r_z=0.5*(hz(i,j,K)+hz(i,j,K-1))
   write(nfile,15)r_x,r_y,r_z
  enddo
 enddo
enddo

15    format(12(E13.6,1x))  !perchè mette 12 e non 3?
return
endsubroutine write_convective
!______________________________________________________________________
!______________________________________________________________________

subroutine write_temperature(nfile)
implicit none
integer i,j,k,nfile
REAL :: r_theta

!open(nfile,file=opfile,form='formatted',position='append')

write(nfile,'(A)') "SCALARS Theta double"
WRITE(nfile,'(A,1x,A)') "LOOKUP_TABLE","default"
!!! parte centrale
do k=1,Nz
 do j=1,Ny
  do i=1,Nx
   r_theta = th(i,j,k)
   write(nfile,'(E13.6)') r_theta
  enddo
 enddo
enddo


return
endsubroutine write_temperature

!______________________________________________________________________
!______________________________________________________________________
subroutine write_divU(nfile)
implicit none
integer i,j,k,nfile
REAL :: r_divU

!open(nfile,file=opfile,form='formatted',position='append')

write(nfile,'(A)') "SCALARS div_U double"
WRITE(nfile,'(A,1x,A)') "LOOKUP_TABLE","default"
!!! parte centrale
do k=1,Nz
 do j=1,Ny
  do i=1,Nx
   r_divU = divU(i,j,k)
   write(nfile,'(E13.6)') r_divU
  enddo
 enddo
enddo


return
endsubroutine write_divU

!!!subroutine write_scalar(T,var_name,nfile)
!!!implicit none
!!!integer i,j,k,nfile
!!!character(*):: var_name
!!!real:: sca
!!!!real,dimension(1:Nr+1,Nz):: T
!!!real(myk),dimension(:,:,:):: T
!!!character(128) :: var_name_string
!!!
!!!var_name_string = 'SCALARS '//var_name//' double 1'
!!!
!!!write(nfile) trim(var_name_string),new_line(c)
!!!write(nfile) "LOOKUP_TABLE default",new_line(c)
!!!
!!!!do k=1,Nz
!!!! do j=1,Nr+1
!!!!  write(nfile) T(j,k)
!!!! enddo
!!!!enddo
!!!  write(nfile) T
!!!
!!!15    format(12(E21.14,1x))
!!!return
!!!endsubroutine write_scalar
!___________________________________________________________________
!___________________________________________________________________
endmodule paraview_module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

program postproc
use mod_fluid
use paraview_module
implicit none

call input()
call allocate_var()
call griglia()

DO it = itmin,itmax,itout

   call leggi()
   call stampa_vtk()

END DO
CALL copy
CALL erase
END PROGRAM postproc
