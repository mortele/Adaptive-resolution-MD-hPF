.DEFAULT_GOAL := all

GFORTRAN 			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -Wextra -Wconversion -std=f2008 -pedantic-errors -ffree-line-length-200
OPTIMIZATION_FLAGS  = -Og 	# All default optimizations which dont interfer with 
							# debugging flags and functionality.

# TODO: Check if -frepack-arrays actually is faster.
# TODO: Check that -Ofast doesnt break something by enabling unsafe math flags.
AGRESSIVE_OPTIMIZE  = -Ofast -faggressive-function-elimination -frepack-arrays -march=native -flto

# TODO: Profile guided optimization; compile with -fprofile-generate and run the 
# code, then compile with -fprofile-use to let gfortran learn from the run time
# execution of the code to optimize it further.

DEBUG_FLAGS 				= -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all
FORTRAN_COMPIILER 			= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS)
FORTRAN_COMPIILER_OPTIMIZE 	= $(GFORTRAN) $(COMPILER_FLAGS) $(AGRESSIVE_OPTIMIZE) 
PROGRAM 					= AdapResoMD-hPF
TEST_PROGRAM 				= UnitTests
FRUIT_SRC 					= ext/FRUIT/src/fruit.f90
FRUIT_MOD 					= fruit.mod fruit_util.mod
FRUIT_OBJ					= fruit.o

SRC = 	app/main.f90 			\
		test/unitTests.f90		\
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

OBJ = 	system.o 				\
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

TEST_SRC = test/system_test.f90

TEST_OBJ = system_test.o

TEST_MOD = system_test.mod

all: $(PROGRAM)

$(PROGRAM): main.o $(OBJ)
	$(FORTRAN_COMPIILER) $(OBJ) main.o -o app/$(PROGRAM).app

main.o: app/main.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

# Make sure parameters.f90 is compiled first, so parameters.mod is available to
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

# Fourth, initialStates.f90 depends on randomGenerator.f90 so we compile that.
randomGenerator.o randomGenerator.mod: src/randomGenerator.f90
	$(FORTRAN_COMPIILER) -c $<

# integrator.f90 depends on potential.f90, so lets do that one next.
potential.o potential.mod: src/potential.f90 parameters.mod
	$(FORTRAN_COMPIILER) -c $<

# The rest of the modules are compiled using the % wild card. The %.o %.mod: 
# possibly only works with the gfortran.
%.o %.mod: src/%.f90 parameters.mod system.mod particles.mod potential.mod particles.mod
	$(FORTRAN_COMPIILER) -c $<

%.o %.mod: test/%.f90 $(MOD) 
	$(FORTRAN_COMPIILER) -c $<

unitTests.o: test/unitTests.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

fruit.o fruit.mod fruit_util.mod: $(FRUIT_SRC)
	$(FORTRAN_COMPIILER) -c $<

test: $(OBJ) $(FRUIT_OBJ) $(TEST_OBJ) unitTests.o
	$(FORTRAN_COMPIILER) $(OBJ) $(FRUIT_OBJ) $(TEST_OBJ) unitTests.o -o test/$(TEST_PROGRAM).app

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o


# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ = xxx
#		$< = code.1
#		$^ = code.1 code.2