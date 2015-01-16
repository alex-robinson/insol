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

# PATH options
objdir = .obj

# Command-line options at make call
ifort ?= 0
debug ?= 0 

## GFORTRAN OPTIONS (default) ##
FC = gfortran
#LIB = /usr/lib
#INC = /usr/include
LIB = /opt/local/lib
INC = /opt/local/include

FLAGS  = -I$(objdir) -J$(objdir) -I$(INC)
LFLAGS = -L$(LIB) -lnetcdff -lnetcdf

DFLAGS = -O3
ifeq ($(debug), 1)
    DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
endif

ifeq ($(ifort),1) 
	## IFORT OPTIONS ##
    FC = ifort 
    LIB = /home/robinson/apps/netcdf/netcdf/lib
    INC = /home/robinson/apps/netcdf/netcdf/include

	FLAGS        = -module $(objdir) -L$(objdir) -I$(INC)
	LFLAGS		 = -L$(LIB) -lnetcdf

	DFLAGS   = -vec-report0 -O3
	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w 
	endif
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/insolation.o: insolation.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/sinsol_orbit.o: sinsol_orbit.f
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

test: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

test_insol0: $(objdir)/ncio.o $(objdir)/sinsol_orbit.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol0.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol0.x is ready."
	@echo " "

insol65N: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol_65N.x $^ test_insol_65N.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol_65N.x is ready."
	@echo " "

clean:
	rm -f test_insol.x $(objdir)/*.o $(objdir)/*.mod

