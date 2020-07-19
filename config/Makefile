.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = ./
objdir = .obj
bindir = ./
testdir = ./

# Command-line options at make call
debug    ?= 0
parallel ?= 0 

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS = $(DFLAGS_PROFILE)
endif

###############################################
##							
## List of rules and source files
##
###############################################

## Individual libraries or modules ##
$(objdir)/ncio.o: ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/interp1D.o: interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: insolation.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/sinsol_orbit.o: sinsol_orbit.f
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

test: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_insol.x $^ $(srcdir)/test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol.x is ready."
	@echo " "

test_insol0: $(objdir)/ncio.o $(objdir)/sinsol_orbit.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_insol0.x $^ $(srcdir)/test_insol.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol0.x is ready."
	@echo " "

insol65N: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/insolation.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -o $(bindir)/test_insol_65N.x $^ $(srcdir)/test_insol_65N.f90 $(LFLAGS)
	@echo " "
	@echo "    test_insol_65N.x is ready."
	@echo " "

clean:
	rm -f test_insol.x $(objdir)/*.o $(objdir)/*.mod

