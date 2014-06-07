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
libdir = ..

# netcdf_inc = /usr/include
# netcdf_lib = /usr/lib
netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib
netcdf_inc_ifort = /home/robinson/apps/netcdf/netcdf/include
netcdf_lib_ifort = /home/robinson/apps/netcdf/netcdf/lib

# Command-line options at make call
ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I$(netcdf_inc_ifort)
	LFLAGS		 = -L$(netcdf_lib_ifort) -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w 
	else
	    DFLAGS   = -vec-report0 -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I$(netcdf_inc)
	LFLAGS		 = -L$(netcdf_lib) -lnetcdff -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow \
	               -fbacktrace -fcheck=all -fbackslash
	else
	    DFLAGS   = -O3 -fbackslash
	endif
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: $(libdir)/ncio/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(libdir)/coord/coordinates.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/coord/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/insolation.o: insolation.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/sinsol_orbit.o: sinsol_orbit.f
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

insol65N: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol_65N.x $^ test_insol_65N.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol_65N.x is ready."
	@echo " "

insol: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

insol0: $(objdir)/ncio.o $(objdir)/sinsol_orbit.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_insol.x $^ test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

clean:
	rm -f test_insol.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
