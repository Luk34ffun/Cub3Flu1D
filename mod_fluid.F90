module mod_fluid
implicit none
save
!functional variables (for storage)
INTEGER :: num_sim
CHARACTER (7) :: fold_num
CHARACTER (15) :: path_fold
CHARACTER(100) :: comment
REAL :: start, finish


!space variables
integer,parameter:: myk=kind(0.0D0)
integer:: Nx,Ny,Nz
INTEGER :: Nt
INTEGER :: m,n
real(myk),dimension(:),allocatable:: X,Y,Z
real(myk) :: Xmax,Ymax,Zmax
real(myk) :: dx,dy,dz
!REAL(myk) :: delta

!time variables
real(myk):: ark(4,2),brk(4,2)
integer:: n_step,i_kutta,i_step
integer:: ik(0:3)
integer:: itmin,itmax,itout,it
real(myk):: dt,t

!Stability variables
REAL(myk) :: Co
REAL(myk) :: ux_max,uy_max,uz_max
REAL(myk) :: U0,L0,Ni0

!file names
character(1):: ver
character(128)::fileres,fileres1,fileres2,fileres3,fileb,fileresb
character(128)::path1,path2,path3,path4,path5,path6,path7
character(256)::opfile,opfile1,opfileth
CHARACTER(len=46) :: matA='A - matrix of the coefficients of the unknowns'
CHARACTER(len=18) :: vecx='x - solution array'
CHARACTER(len=20) :: vecb='b - known term array'
CHARACTER(1) :: b_vec='b'
CHARACTER(len=3) :: jac='jac'
CHARACTER(len=3) :: gs='g_s'
CHARACTER(len=3) :: sor='sor'
CHARACTER(len=4) :: extns='.dat'
CHARACTER(9) :: filename
CHARACTER(128) :: filename1
CHARACTER(50) :: resfile,resfile1
CHARACTER(9) :: filePhi
CHARACTER(8) :: fileth
CHARACTER(50) :: resPhi
CHARACTER(39) :: resth
CHARACTER(12) :: path8='temperature/'
CHARACTER(10) :: divU,filere
!variables
real(myk),dimension(:,:,:),allocatable:: mu

!Alia
REAL(myk),PARAMETER :: Pi = 3.1415927_myk

!Physical variables
real(myk),dimension (:,:,:), allocatable:: ux,uy,uz
real(myk),dimension (:,:,:), allocatable:: H_ux,H_uy,H_uz,H_th
real(myk),dimension (:,:,:), allocatable:: P_ux,P_uy,P_uz,P_th
real(myk),dimension (:,:,:), allocatable:: Txx,Tyy,Tzz,Txy,Txz,Tyz
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: diff_x, diff_y, diff_z, diff_th
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: hx, hy, hz, fz, l
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: th !dimensionless temperature
REAL(myk), DIMENSION (1,2) :: q_dot !Volumetric flow rate
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: div_V, Phi, div_U
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: grad_Phi_x,grad_Phi_y,grad_Phi_z
REAL(myk), DIMENSION (:,:,:), ALLOCATABLE :: vz
REAL(myk) :: maxdivU

!!Opening geometrical variables
REAL(myk) :: x_1,y_1,r_1,  x_2,y_2,r_2
! handle variable to tuning the amplitude of bivariate distrib of vel in bottom entry
REAL(myk) :: sigma_u 
!!Iteration method variables
REAL(myk), DIMENSION (:,:), ALLOCATABLE :: A_Phi 
REAL(myk), DIMENSION (:), ALLOCATABLE :: X0,b_Phi,X_vx,X_vy,X_vz
REAL(myk), DIMENSION (:,:), ALLOCATABLE :: X_Phi
INTEGER :: method,max_iter,last_i
REAL(myk) :: max_err,omg
!physical parameters
real(myk):: Re,Inv_Re,Pr,Ri,Inv_RePr

!temp variables
REAL(myk) :: Th_s                                ![k] max temp @ the buoyancy source
REAL(myk),PARAMETER :: Th_0 = 273.15_myk+20.0_myk![k] temp of the box 
                                                 !box considered as an adiabatic container
                                                 !and as an ideal temperature source

REAL(myk) :: dTh0

!Variable for vzmaxbound subroutine
REAL(myk),DIMENSION(:,:,:), ALLOCATABLE ::vz0
INTEGER :: i_max,j_max
REAL(myk) :: vz0max

contains
end
