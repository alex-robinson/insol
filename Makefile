.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make insol      : compiles the main program test_insol.x"
	@echo " make clean      : cleans object files"
	@echo ""

objdir = .obj

ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I/home/robinson/apps/netcdf/netcdf/include
	LFLAGS		 = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w 
	else
	    DFLAGS   = -vec-report0 -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I/opt/local/include
	LFLAGS		 = -L/opt/local/lib -lnetcdff -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
	else
	    DFLAGS   = -O3
	endif
endif

## Cluster ##
# FC			 = ifort
# FLAGS          = -module $(objdir) -L$(objdir) -I/home/robinson/apps/netcdf/netcdf/include
# DFLAGS         = -w -C -traceback -ftrapuv -fpe0 -check all -vec-report0
# RELEASEFLAGS   = -vec-report0 -O3
# LFLAGS		 = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

## Individual libraries or modules ##
$(objdir)/ncio3.o: ../ncio/ncio3.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: ../coord/coordinates.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: ../coord/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/insolation.o: insolation.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/sinsol_orbit.o: sinsol_orbit.f
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

insol: $(objdir)/ncio3.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

insol0: $(objdir)/ncio3.o $(objdir)/sinsol_orbit.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

clean:
	rm -f test_insol.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
