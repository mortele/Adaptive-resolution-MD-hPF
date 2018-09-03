.DEFAULT_GOAL := all

GFORTRAN  			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -Wextra -Wconversion -std=f2008 -pedantic-errors -ffree-line-length-0 -Jbuild -ftest-coverage -fprofile-arcs
OPTIMIZATION_FLAGS  = -Og 	# All default optimizations which dont interfer with 
							# debugging flags and functionality.

# TODO: Check if -frepack-arrays actually is faster.
# TODO: Check that -Ofast doesnt break something by enabling unsafe math flags.
AGRESSIVE_OPTIMIZE  = -Ofast -faggressive-function-elimination -frepack-arrays -march=native -flto

# TODO: Profile guided optimization; compile with -fprofile-generate and run the 
# code, then compile with -fprofile-use to let gfortran learn from the run time
# execution of the code to optimize it further.

DEBUG_FLAGS 				= -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all -g
TEST_FLAGS 					= -fprofile-arcs #-ftest-coverage
FORTRAN_COMPILER 			= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS)
FORTRAN_COMPILER_OPTIMIZE 	= $(GFORTRAN) $(COMPILER_FLAGS) $(AGRESSIVE_OPTIMIZE) 
PROGRAM 					= AdapResoMD-hPF.app
TEST_PROGRAM 				= UnitTests.app
FRUIT_SRC 					= ext/FRUIT/src/fruit.f90
FRUIT_MOD 					= fruit.mod fruit_util.mod
FRUIT_OBJ					= fruit.o

VPATH = src:test:build:bin

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
	$(FORTRAN_COMPILER) $(patsubst %.o, build/%.o, $(OBJ)) -o bin/$(PROGRAM) 

test: $(filter-out main.o, $(OBJ)) $(TEST_OBJ) fruit.o
	$(FORTRAN_COMPILER) $(patsubst %.o, build/%.o, $(filter-out main.o, $(OBJ))) $(patsubst %.o, build/%.o, $(TEST_OBJ) fruit.o) -o bin/$(TEST_PROGRAM) 

field.o field.mod: src/field.f90 system.o parameters.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
file_writer.o file_writer.mod: src/file_writer.f90 parameters.o system.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
initial_states.o initial_states.mod: src/initial_states.f90 parameters.o system.o random_generator.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
random_generator.o random_generator.mod: src/random_generator.f90 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
integrator.o integrator.mod: src/integrator.f90 parameters.o system.o particles.o potential.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
main.o: src/main.f90 parameters.o sampler.o potential.o integrator.o file_writer.o initial_states.o particles.o
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
parameters.o parameters.mod: src/parameters.f90 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
particles.o particles.mod: src/particles.f90 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
potential.o potential.mod: src/potential.f90 parameters.o particles.o system.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
sampler.o sampler.mod: src/sampler.f90 potential.o particles.o parameters.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
system.o system.mod: src/system.f90 parameters.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
system_test.o: test/system_test.f90 fruit.o system.o particles.o parameters.o 
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
random_generator_test.o: test/random_generator_test.f90 random_generator.o fruit.o
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@) 
potential_test.o: test/potential_test.f90 potential.o system.o parameters.o fruit.o
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@) 
unit_tests.o: unit_tests.f90 $(patsubst build/%.o, %.o, $(wildcard build/*_test.o))
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)
fruit.o fruit.mod: ext/FRUIT/src/fruit.f90
	$(FORTRAN_COMPILER) -c $< -o build/$(notdir $@)


#%.o: src/%.f90
#	$(FORTRAN_COMPILER) -c $<

depend depend.mk:
	./ext/makedepf90/makedepf90 src/*.f90 test/*.f90 -m '%m.mod' -b '' > depend.mk

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o test/*.mod test/*.o4 *.gcda *.gcno	build/*.mod build/*.o build/*.gcda build/*.gcno




# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ = xxx
#		$< = code.1
#		$^ = code.1 code.2