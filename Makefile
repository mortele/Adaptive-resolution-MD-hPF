.DEFAULT_GOAL := all

GFORTRAN  			= gfortran
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
TEST_FLAGS 					= -fprofile-arcs #-ftest-coverage
FORTRAN_COMPILER 			= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS)
FORTRAN_COMPILER_OPTIMIZE 	= $(GFORTRAN) $(COMPILER_FLAGS) $(AGRESSIVE_OPTIMIZE) 
PROGRAM 					= AdapResoMD-hPF
TEST_PROGRAM 				= UnitTests
FRUIT_SRC 					= ext/FRUIT/src/fruit.f90
FRUIT_MOD 					= fruit.mod fruit_util.mod
FRUIT_OBJ					= fruit.o

VPATH = src test build bin

SRC  = $(wildcard src/*.f90)
OBJ  = $(patsubst src/%.f90, %.o,   $(SRC))
MODD = $(patsubst src/%.f90, %.mod, $(SRC))

# Main does not produce a module, just the .o file, so we remove it from the 
# list of module files.
MOD  = $(filter-out main.mod, $(MODD)) 

TEST_SRC = $(wildcard test/*.f90)
TEST_OBJ = $(patsubst test/%.f90, %.o,   $(TEST_SRC))
TEST_MOD = $(patsubst test/%.f90, %.mod, $(TEST_SRC))

all: $(OBJ)
	$(FORTRAN_COMPILER) $(OBJ) -o $(PROGRAM) 

test: $(filter-out main.o, $(OBJ)) $(TEST_OBJ) fruit.o
	$(FORTRAN_COMPILER) $(filter-out main.o, $(OBJ)) fruit.o $(TEST_OBJ) -o $(TEST_PROGRAM) 

field.o field.mod: src/field.f90 system.o parameters.o 
	$(FORTRAN_COMPILER) -c $<
file_writer.o file_writer.mod: src/file_writer.f90 parameters.o system.o 
	$(FORTRAN_COMPILER) -c $<
initial_states.o initial_states.mod: src/initial_states.f90 parameters.o system.o random_generator.o 
	$(FORTRAN_COMPILER) -c $<
random_generator.o random_generator.mod: src/random_generator.f90 
	$(FORTRAN_COMPILER) -c $<
integrator.o integrator.mod: src/integrator.f90 parameters.o system.o particles.o potential.o 
	$(FORTRAN_COMPILER) -c $<
main.o: src/main.f90 parameters.o sampler.o potential.o integrator.o file_writer.o initial_states.o particles.o
	$(FORTRAN_COMPILER) -c $<
parameters.o parameters.mod: src/parameters.f90 
	$(FORTRAN_COMPILER) -c $<
particles.o particles.mod: src/particles.f90 
	$(FORTRAN_COMPILER) -c $<
potential.o potential.mod: src/potential.f90 parameters.o particles.o system.o 
	$(FORTRAN_COMPILER) -c $<
sampler.o sampler.mod: src/sampler.f90 potential.o particles.o parameters.o 
	$(FORTRAN_COMPILER) -c $<
system.o system.mod: src/system.f90 parameters.o 
	$(FORTRAN_COMPILER) -c $<
system_test.o: test/system_test.f90 fruit.o system.o particles.o parameters.o 
	$(FORTRAN_COMPILER) -c $<
unit_tests.o: unit_tests.f90 system_test.o 
	$(FORTRAN_COMPILER) -c $<
fruit.o fruit.mod: ext/FRUIT/src/fruit.f90
	$(FORTRAN_COMPILER) -c $<

#%.o: src/%.f90
#	$(FORTRAN_COMPILER) -c $<

depend depend.mk:
	./ext/makedepf90/makedepf90 src/*.f90 test/*.f90 -m '%m.mod' -b '' > depend.mk

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o *.gcda *.gcno	build/*.mod build/*.o build/*.gcda build/*.gcno




# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ = xxx
#		$< = code.1
#		$^ = code.1 code.2