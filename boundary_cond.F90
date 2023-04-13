SUBROUTINE boundary_cond(f)
USE mod_fluid
IMPLICIT NONE
INTEGER :: i,j,k
INTEGER,INTENT(IN) :: f
REAL(myk) :: x_th,y_th,sigma_th,r_sq,r
REAL(myk) :: tau
tau=itmax*dt
!!Geometrical definition of the  openings
!Bottom opening Specs
 !x_1=0.25_myk
 !y_1=Ymax*0.5_myk
 !r_1= 0.05_myk
 !sigma_u = 0.03_myk
!Top opening Specs
 !y_2=y_1
 !x_2=0.6_myk
 !r_2=r_1

!Characteristic Delta_temperature
 dTh0 = Th_0 - Th_s

!buoyancy source specs
x_th=Xmax*0.5_myk
y_th=Ymax*0.5_myk
sigma_th=0.03

!Reference value for vel & temp
!Th0 = 273.15D0 + 40 ![K] max temp @ the buoyancy source (@ center of it)
!th_b = 273.15_myk + 20.0_myk  ![k] temp of the box surface

SELECT CASE (f)

!CASE(f==1) -> B.C. on the initial velocity (u-field) & on the ausiliary velocity (v-field)
CASE(1)
!z=0 (bottom surface with opening)
DO j=1,Ny
  DO i=1,Nx
   r = sqrt((X(i) -x_1)**2+(Y(j)-y_1)**2)
    IF (r<r_1) THEN

     !inflow conditions (ux=0 & uy=0 & uz=uz(r) -> inflow has only z component of vel.)
     ux(I,j,0) = -ux(I,j,1)
     uy(i,J,0) = -uy(i,J,1)
     uz(i,j,0) = (1._myk-(r/r_1)**2)*exp(-t/tau)
     !uz(i,j,0) = exp(-(r**2)/(1.4*sigma_u**2))!*exp(-t/tau)
    
    ELSE
     
     !no slip
     ux(I,j,0) = -ux(I,j,1)
     uy(i,J,0) = -uy(i,J,1)
     !impermeability
     uz(i,j,0) = 0.0_myk
     
     END IF
  END DO
END DO

!z=1 (top surface with opening)
DO j=1,Ny
  DO i=1,Nx
   r = sqrt((X(i)-x_2)**2+(Y(j)-y_2)**2)
    IF (r<r_2) THEN


     !Memo: la condizione sul flusso in uscita è imposta assumendo noto il valore di Phi
     !      ciò implica che il valore sulla velocità non può essere assegnato e si impone
     !      alternativamente una condizione sul valore della sua derivata.

     !outflow condition (dux/dz=0 & duy/dz=0 & duz/dz=0)
     ux(I,j,Nz+1) = ux(I,j,Nz  )
     uy(i,J,Nz+1) = uy(i,J,Nz  )                                      !  *********( 1 )*********
     uz(i,j,Nz+1) = uz(i,j,Nz-1) 
     uz(i,j,Nz  ) = uz(i,j,Nz-1)     



    ELSE
     
     !no slip
     ux(I,j,Nz+1) = -ux(I,j,Nz)
     uy(i,J,Nz+1) = -uy(i,J,Nz)
     !impermeability
     uz(i,j,Nz) = 0.0_myk
     
    END IF
  END DO
END DO


!y=0 & y=1 (solid surfaces)
DO k=1,Nz !k=1,Nz cambio estremi
  DO i=1,Nx
    
    !impermeability
     !y=0
     uy(i,0,k) = 0.0_myk
     !y=1
     uy(i,Ny,k) = 0.0_myk
 
    !no slip
     !y=0
     ux(I,0,k) = -ux(I,1,k)
     uz(i,0,K) = -uz(i,1,K) 
    
     !y=1
     ux(I,Ny+1,k) = -ux(I,Ny,k)
     uz(i,Ny+1,K) = -uz(i,Ny,K)
  
  ENDDO 
ENDDO


!x=0 & x=1 (solid surfaces)
DO k=1,Nz
  DO j=1,Ny
     
     !impermeability
      !x=0
      ux(0,j,k) = 0.0_myk
      !x=1
      ux(Nx,j,k) = 0.0_myk

     !no slip 
      !x=0
      uy(0,J,k) = -uy(1,J,k)
      uz(0,j,K) = -uz(1,j,K) 
      !x=1
      uy(Nx+1,J,k) = -uy(Nx,J,k)
      uz(Nx+1,j,K) = -uz(Nx,j,K)
     
     
  END DO
END DO

CASE(2)
!------ Temperature field --------------------------------------

!z=0 (low surface of the box - where the bouyancy source is placed)
  DO j=1,Ny
    DO i=1,Nx
     r_sq=(X(i)-x_th)**2+(Y(j)-y_th)**2
     IF (sqrt(r_sq)>sigma_th) THEN

     !imposing the temperature of the box
      th(i,j,0) = 2._myk - th(i,j,1)
     
     ELSE
     
     !buoyancy source
      th(i,j,0) = 2._myk*((Th_s-Th_0)/Th_0*exp(-0.5_myk*r_sq/(sigma_th**2))+1) - th(i,j,1)
     END IF

!z=1 (top surface of the box)
     th(i,j,Nz+1) = 2._myk - th(i,j,Nz)

  END DO
END DO

!x=0 & x=1
DO k=1,Nz
  DO j=1,Ny
    
    th(   0,j,k) = 2._myk - th( 1,j,k)
    th(Nx+1,j,k) = 2._myk - th(Nx,j,k)
  
  END DO
END DO

!y=0 & y=1
DO k=1,Nz
  DO i=1,Nx

   th(i,   0,k) = 2._myk - th(i, 1,k)
   th(i,Ny+1,k) = 2._myk - th(i,Ny,k)
  
  END DO
END DO
!---------------------------------------------------------------

!CASE(field==3) ->B.C. on the irrotational field (grad(Phi))

!Anche se le equazioni relative alle b.c. sono già state usate  per risolvere il sistema Ax=b per Phi,
!formalmente sono solamente stati risolti i nodi interni e non i ghost nodes !
!queste b.c. mi permettono di assegnare i valori ai ghost nodes con i valori del campo
!Phi all'interno del dominio fluido già noti!!! 

CASE(3)
!z=0 (bottom surface with openings) 
DO j=1,Ny
  DO i=1,Nx
    r = sqrt((X(i)-x_1)**2+(Y(j)-y_1)**2)
    IF (r<r_1) THEN
    
    !inflow conditions ( (d(Phi)/dz = 0 @ z=0) (Neumann b.c.)

    Phi(i,j,0) =  Phi(i,j,1)

    ELSE!solid part of surface

    !impermeability (d(Phi)/dz =0)
    Phi(i,j,0) = Phi(i,j,1)

    END IF

  END DO
END DO

!z=1 (top surface with openings)
DO j=1,Ny
  DO i=1,Nx
    r = sqrt((X(i)-x_2)**2+(Y(j)-y_2)**2)
    IF (r<r_2) THEN
    
    !outflow conditions (Phi = 0 @ z=1)  (Dirichlet b.c.)
    Phi(i,j,Nz+1) = -Phi(i,j,Nz)!+(Phi(i,j,1)+Phi(i,j,0))/2

    ELSE!solid part of surface

    !impermeability (d(Phi)/dz=0)
    Phi(i,j,Nz+1) = Phi(i,j,Nz)

    END IF

  END DO
END DO

!y=0 & y=1 (solid surfaces)
DO k=1,Nz
  DO i=1,Nx
    
    !impermeability(d(Phi)/dy=0)
     !y=0
     Phi(i,0,k) = Phi(i,1,k)
     !y=1
     Phi(i,Ny+1,k) = Phi(i,Ny,k)
   
   ENDDO
ENDDO


!x=0 & x=1 (solid surface)
DO k=1,Nz
  DO j=1,Ny
       
       !impermeability(d(Phi)/dx=0)
        !x=0
        Phi(0,j,k) = Phi(1,j,k)
        !x=1
        Phi(Nx+1,j,k) = Phi(Nx,j,k)
  END DO
END DO


END SELECT
RETURN
END SUBROUTINE boundary_cond
