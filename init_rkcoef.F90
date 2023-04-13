      subroutine init_rkcoef()
      use mod_fluid
!______________________________________________________________________
!
!     assegna i coefficienti di Runge Kutta
!     e calcolo il puntatore ciclico ik che mi permette
!     di fare R-K su un indice it non ciclico
!
!______________________________________________________________________
      implicit none
!e______________________________________________________________________
!     Assegno i coefficienti di R-K
      Ark(1,1)=8./15.
      Ark(2,1)=5./12.
      Ark(3,1)=3./4.
      Brk(1,1)=0.
      Brk(2,1)=-17./60.
      Brk(3,1)=-5./12.
      Ark(1,2)=8./17.
      Ark(2,2)=17./60.
      Ark(3,2)=5./12.
      Ark(4,2)=3./4.
      Brk(1,2)=0.
      Brk(2,2)=-15./68.
      Brk(3,2)=-17./60.
      Brk(4,2)=-5./12.
!______________________________________________________________________
!     Calcolo il contatore ciclico
      ik(0)=n_step
      ik(1)=1
      ik(2)=2
      ik(3)=3
!______________________________________________________________________
      end
                                           
