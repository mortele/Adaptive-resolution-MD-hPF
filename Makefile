GFORTRAN 			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -std=f2008ts
FORTRAN_COMPIILER 	= $(GFORTRAN) $(COMPILER_FLAGS)
PROGRAM 			= AdapResoMD-hPF

SRC = 	app/main.f90 			\
		src/system.f90 			\
		src/particles.f90		\
		src/potential.f90		\
		src/parameters.f90		\
		src/field.f90			\
		src/randomGenerator.f90	\
		src/initialStates.f90 	\

OBJ = 	main.o 					\
		system.o 				\
		particles.o 			\
		potential.o 			\
		parameters.o 			\
		field.o 				\
		randomGenerator.o 		\
		initialStates.o 		\

MOD = 	system.mod 				\
		particles.mod 			\
		potential.mod 			\
		parameters.mod 			\
		field.mod				\
		randomGenerator.mod 	\
		initialStates.mod 		\

all: $(OBJ)
	$(FORTRAN_COMPIILER) $(OBJ) -o app/$(PROGRAM).app

main.o: app/main.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

# Make sure parameters.f90 is kompiled first, so parameters.mod is available to
# the other modules at compile time.
parameters.o parameters.mod: src/parameters.f90
	$(FORTRAN_COMPIILER) -c $^

# Secondly, we compile system.f90 so that other modules which depend on it can 
# use system.mod.
system.o system.mod: src/system.f90 parameters.mod 
	$(FORTRAN_COMPIILER) -c $^

# The rest of the modules are compiled using the % wild card. The %.o %.mod: 
# possibly only works with the gfortran.
%.o %.mod: src/%.f90 parameters.mod system.mod 
	$(FORTRAN_COMPIILER) -c $^

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o


# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ xxx
#		$< code.1
#		$^ code.1 code.2