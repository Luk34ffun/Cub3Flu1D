SUBROUTINE buildmatrix
!!This Subroutine build the matrices A of the unknown term coefficients of the linear system Ax=b
!! the matrix is defined starting from its principal diagonal terms and wrt these.
USE mod_fluid 
IMPLICIT NONE
INTEGER   :: i,j,k,r !rows 
INTEGER   :: ii,jj,kk,c !columns
REAL(myk) :: w_1, w_Nx, w_NxNy !weights
REAL(myk) :: a0,a_1,a_Nx,a_NxNy !coefficients
REAL(myk) :: rb,rt
REAL(myk) :: d_1  !inv_den

d_1 = 1._myk/(dx*dy*dz)**2 
!define the coefficients of the non-zero terms
a0     = d_1*-2._myk*((dy*dz)**2+(dx*dz)**2+(dx*dy)**2)
a_1    = d_1*(dy*dz)**2
a_Nx   = d_1*(dx*dz)**2
a_NxNy = d_1*(dx*dy)**2
!Inizializzo la matrice con tutti zeri; seguitamente partendo dai termini sulla diagonale principale,
!spostandomi sulle colonne, assegno i valori ai termini delle altre 6 diagonali minori.
A_phi=0._myk
!------------Definition unknown term coefficients Matrix A  ----------
DO k=1,Nz
  DO j=1,Ny
    DO i=1,Nx
     r=i+(j-1)*Nx+(k-1)*Nx*Ny   

      DO kk=1,Nz
        DO jj=1,Ny
          DO ii=1,Nx
           c=ii+(jj-1)*Nx+(kk-1)*Nx*Ny
    ext1:    IF (c.EQ.r) THEN !mi posiziono sulla diag principale
            
            !conditions block
    cblock:  IF (kk.EQ.1) THEN                                                                                    ![IF]
                                                                          
              IF ( (jj.EQ.1).or.(jj.EQ.Ny) ) THEN                                            ![IF]
                 
               IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN                           ![IF] 
                 w_NxNy = 1._myk                   
                 IF (Ny.EQ.1) THEN !caso 2D             ![IF]
                  w_Nx   = 2._myk
                 ELSE              !caso 3D             ![ELSE]
                  w_Nx   = 1._myk
                 END IF                                 ![ENDIF]
                 w_1    = 1._myk                   
                ELSE !(ii.ne.1 & ii.ne.Nx)                                   ![ELSE]                          
                 w_NxNy = 1._myk
                 IF (Ny.EQ.1) THEN !caso 2D     ![IF]
                  w_Nx  = 2._myk
                 ELSE              !caso 3D     ![ELSE]
                  w_Nx = 1._myk
                 END IF                         ![ENDIF]
                 w_1    = 0._myk                   
                END IF ![ic1]                                                ![ENDIF]

               ELSE !(jj.ne.1 & jj.ne.Ny)                                                   ![ELSE]
                IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN            ![IF]
                 w_NxNy = 1._myk
                 w_Nx   = 0._myk
                 w_1    = 1._myk
                ELSE !(ii.ne.1 & ii.ne.Nx)                     ![ELSE]
                 rb=sqrt((X(ii)-x_1)**2+(Y(jj)-y_1)**2)
                  IF (rb>r_1) THEN    ![IF]
                   w_NxNy =  1._myk
                   w_Nx   =  0._myk
                   w_1    =  0._myk
                  ELSE !(rb<r_1)      ![ELSE]
                   w_NxNy =  1._myk   !sempre Neumann b.c. ( d(Phi)/dz = 0 -> Phi(i,j,0) = Phi(i,j,1) )
                   w_Nx   =  0._myk
                   w_1    =  0._myk
                  END IF              ![ENDIF]
                END IF                                         ![ENDIF]
               END IF ![jc1]                                                                 ![ENDIF]

              !proseguo col cblock
              ELSE IF (kk.EQ.Nz) THEN !-------------------------------------                                      ![ELSEIF]
               
               ! [jc2:]  
               IF ( (jj.EQ.1).or.(jj.EQ.Ny) ) THEN                                ![IF]
                
                ! [ic2:]  
                IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN                    ![IF]
                 w_NxNy = 1._myk
                 IF (Ny.EQ.1) THEN  !caso 2D           ![IF]
                  w_Nx  = 2._myk
                 ELSE               !caso 3D           ![ELSE]
                  w_Nx   = 1._myk
                 END IF                                ![ENDIF]
                 w_1    = 1._myk
                ELSE !(ii.ne.1 & ii.ne.Nx)                             ![ELSE]
                 w_NxNy = 1._myk
                 IF (Ny.EQ.1) THEN  !caso 2D        ![IF]
                  w_Nx  = 2._myk
                 ELSE                               ![ELSE]
                  w_Nx   = 1._myk   !caso 3D
                 END IF                             ![ENDIF]
                 w_1    = 0._myk
                END IF! [ic2]                                          ![ENDIF]           

               ELSE !(jj.ne.1 & jj.ne.Ny)                                         ![ELSE] 
                IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN        ![IF]
                 w_NxNy = 1._myk
                 w_Nx   = 0._myk
                 w_1    = 1._myk
                ELSE !(ii.ne.1 & ii.ne.Nx)                 ![ELSE]
                 rt=sqrt((X(ii)-x_2)**2+(Y(jj)-y_2)**2)
                  IF (rt>r_2) THEN       ![IF]
                   w_NxNy =  1._myk !(impermeability->Neumann b.c.)
                   w_Nx   =  0._myk
                   w_1    =  0._myk
                  ELSE !(rt<r_2)         ![ELSE]
                   w_NxNy = -1._myk !(opening->Dirichlet b.c.)
                   w_Nx   =  0._myk
                   w_1    =  0._myk
                  END IF                 ![ENDIF]
                END IF                                      ![ENDIF]
               END IF ![jc2]                                                      ![ENDIF]
             
             ELSE !(kk.ne.1 & kk.ne.Nz)                                                                           ![ELSE]
              IF ( (jj.EQ.1).or.(jj.EQ.Ny) ) THEN                          ![IF]

               IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN              ![IF]
                w_NxNy =  0._myk
                IF (Ny.EQ.1) THEN !caso 2D         ![IF]
                 w_Nx   =  2._myk
                ELSE              !caso 3D         ![ELSE]
                 w_Nx   =  1._myk
                END IF                             ![ENDIF]
                w_1    =  1._myk
               ELSE !(ii.ne.1 & ii.ne.Nx)                       ![ELSE]
                w_NxNy =  0._myk
                IF (Ny.EQ.1) THEN  !caso 2D    ![IF]
                 w_Nx   =  2._myk
                ELSE               !caso 3D    ![ELSE]
                 w_Nx   =  1._myk
                END IF                         ![ENDIF]
                w_1    =  0._myk                                
               END IF                                           ![ENDIF]

              ELSE !(jj.ne.1 & jj.ne.Ny)                                    ![ELSE]
               IF ( (ii.EQ.1).or.(ii.EQ.Nx) ) THEN  ![IF]
                w_NxNy =  0._myk
                w_Nx   =  0._myk
                w_1    =  1._myk
               ELSE !(ii.ne.1 & ii.ne.Nx)           ![ELSE]
                w_NxNy =  0._myk
                w_Nx   =  0._myk
                w_1    =  0._myk
               END IF                               ![ENDIF]
              
              END IF                                                        ![ENDIF]
             
             END IF cblock                                                                                        ![ENDIF]

           A_Phi(r,c) = a0+w_1*a_1+w_Nx*a_Nx+w_NxNy*a_NxNy !cond ext1
          
          !Characterization of the terms belonging to the diagonals wrt the principal diagonal
          IF (c-1>0) THEN
           IF (ii.NE.1) THEN ! (if ii == 1 -> nodo boundary)
            A_Phi(r,c-1) = a_1
           END IF
          END IF

          IF (c+1 .LE. Nt) THEN
           IF (ii.NE.Nx) THEN ! (if ii == Nx -> nodo boundary)
            A_Phi(r,c+1) = a_1
           END IF
          END IF

          IF ((c-Nx >0).and.(jj.NE.1)) THEN  !(if jj == 1 -> nodo boundary)
           A_Phi(r,c-Nx      ) = a_Nx
          END IF

          IF ((c+Nx .LE. Nt).AND.(jj.NE.Ny)) THEN !(if jj == Ny -> nodo boundary)
           A_Phi(r,c+Nx      ) = a_Nx 
          END IF

          IF ((c-Nx*Ny>0).AND.(kk.NE.1)) THEN !(if kk == 1 -> nodo boundary)
           A_Phi(c,c-Nx*Ny) = a_NxNy  
          END IF

          IF ((c+Nx*Ny .LE. Nt).AND.(kk.NE.Nz)) THEN !(if kk == 1 -> nodo boundary) 
           A_Phi(c,c+Nx*Ny)  = a_NxNy
          END IF
           
        END IF ext1 
        END DO
     END DO
    END DO
    
  END DO
 END DO
END DO
CALL checkmatrix
RETURN
END SUBROUTINE buildmatrix

