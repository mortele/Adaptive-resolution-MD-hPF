GFORTRAN 			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -Wextra -Wconversion -std=f2008 -pedantic-errors
OPTIMIZATION_FLAGS  = #-O2 

# TODO: Check if -frepack-arrays actually is faster.
AGRESSIVE_OPTIMIZE  = #-faggressive-function-elimination -frepack-arrays  

DEBUG_FLAGS 		= -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all
FORTRAN_COMPIILER 	= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(AGRESSIVE_OPTIMIZE) $(DEBUG_FLAGS)
PROGRAM 			= AdapResoMD-hPF

SRC = 	app/main.f90 			\
		src/system.f90 			\
		src/particles.f90		\
		src/potential.f90		\
		src/parameters.f90		\
		src/field.f90			\
		src/randomGenerator.f90	\
		src/initialStates.f90 	\
		src/fileWriter.f90		\
		src/integrator.f90		\
		src/sampler.f90			\

OBJ = 	main.o 					\
		system.o 				\
		particles.o 			\
		potential.o 			\
		parameters.o 			\
		field.o 				\
		randomGenerator.o 		\
		initialStates.o 		\
		fileWriter.o 			\
		integrator.o 			\
		sampler.o 				\

MOD = 	system.mod 				\
		particles.mod 			\
		potential.mod 			\
		parameters.mod 			\
		field.mod				\
		randomGenerator.mod 	\
		initialStates.mod 		\
		fileWriter.mod 			\
		integrator.mod 			\
		sampler.mod 			\

all: $(OBJ)
	$(FORTRAN_COMPIILER) $(OBJ) -o app/$(PROGRAM).app

main.o: app/main.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

# Make sure parameters.f90 is kompiled first, so parameters.mod is available to
# the other modules at compile time.
parameters.o parameters.mod: src/parameters.f90
	$(FORTRAN_COMPIILER) -c $<

# Secondly, we compile system.f90 so that other modules which depend on it can 
# use system.mod.
system.o system.mod: src/system.f90 parameters.mod 
	$(FORTRAN_COMPIILER) -c $<

# Thirdly, make sure particles.f90 is compiled so other modules can use it.
particles.o particles.mod: src/particles.f90 
	$(FORTRAN_COMPIILER) -c $<

# The rest of the modules are compiled using the % wild card. The %.o %.mod: 
# possibly only works with the gfortran.
%.o %.mod: src/%.f90 parameters.mod system.mod 
	$(FORTRAN_COMPIILER) -c $<

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o


# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ xxx
#		$< code.1
#		$^ code.1 code.2