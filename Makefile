FC=gfortran
#FFLAGS= -O3 -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow

FFLAGS= -O3 -Wall -Wextra -fcheck=all -ffpe-trap=zero,overflow

SRC=mod_fluid.F90 input.F90 allocate_var.F90 griglia.F90 init_rkcoef.F90  leggi.F90 init_field.F90 boundary_cond.F90 prhs.F90 rhs.F90 linear.F90 stampa.F90 main.F90 buildmatrix.F90 diver_V.F90 inv_rho_stag.F90 known_term.F90 Phi_iter.F90 associate_var.F90 grad_Phi.F90 U_field.F90 clean.F90 vzmaxbound.F90 diver_U.F90 vol_flowrate.F90 checkmatrix.F90 sim_spec.F90 courant.F90

OBJ=${SRC:.F90=.o}

%.o: %.F90

			$(FC) $(FFLAGS) -o $@ -c $<
fluid.exe: $(OBJ)
			$(FC) $(FFLAGS) -o $@  $(OBJ) 

clean:
			@rm -f *.o *.mod fluid

