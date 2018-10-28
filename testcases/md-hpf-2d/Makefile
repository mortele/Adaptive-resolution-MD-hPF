.DEFAULT_GOAL := all

GFORTRAN  			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -Wextra -Wconversion -std=f2008 -pedantic-errors -ffree-line-length-0 -J../../build -ftest-coverage -fprofile-arcs
OPTIMIZATION_FLAGS  = #-Og 	# All default optimizations which dont interfer with 
							# debugging flags and functionality.

# TODO: Check if -frepack-arrays actually is faster.
# TODO: Check that -Ofast doesnt break something by enabling unsafe math flags.
AGRESSIVE_OPTIMIZE  = -Ofast -faggressive-function-elimination -frepack-arrays -march=native -flto

# TODO: Profile guided optimization; compile with -fprofile-generate and run the 
# code, then compile with -fprofile-use to let gfortran learn from the run time
# execution of the code to optimize it further.

DEBUG_FLAGS 				= -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all -g
TEST_FLAGS 					= -fprofile-arcs #-ftest-coverage
FORTRAN_COMPILER_DEBUG 		= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS)
#FORTRAN_COMPILER 			= $(GFORTRAN) -J../../build -ffree-line-length-0
FORTRAN_COMPILER        	= $(GFORTRAN) $(COMPILER_FLAGS) $(AGRESSIVE_OPTIMIZE) 
PROGRAM 					= md-hpf-2d.app

VPATH = ../../src:../../test:../../build:../../bin

SRC  = $(wildcard ../../src/*.f90)
OBJ  = $(patsubst ../../src/%.f90, %.o,   $(SRC))
MODD = $(patsubst ../../src/%.f90, %.mod, $(SRC))

# Main does not produce a module, just the .o file, so we remove it from the 
# list of module files.
MOD  = $(filter-out main.mod, $(MODD))
OBJJ = $(filter-out main.o,   $(OBJ))

all: md-hpf-2d.o $(OBJJ)
	$(FORTRAN_COMPILER) ../../build/md-hpf-2d.o $(patsubst %.o, ../../build/%.o, $(OBJJ)) -o $(PROGRAM) 

field.o field.mod: ../../src/field.f90 system.o parameters.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
file_writer.o file_writer.mod: ../../src/file_writer.f90 parameters.o system.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
initial_states.o initial_states.mod: ../../src/initial_states.f90 parameters.o system.o random_generator.o particles.o
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
random_generator.o random_generator.mod: ../../src/random_generator.f90 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
integrator.o integrator.mod: ../../src/integrator.f90 parameters.o system.o particles.o potential.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
md-hpf-2d.o: md-hpf-2d.f90 parameters.o sampler.o potential.o integrator.o file_writer.o initial_states.o particles.o field.o
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
parameters.o parameters.mod: ../../src/parameters.f90 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
particles.o particles.mod: ../../src/particles.f90 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
potential.o potential.mod: ../../src/potential.f90 parameters.o particles.o system.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
sampler.o sampler.mod: ../../src/sampler.f90 potential.o particles.o parameters.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)
system.o system.mod: ../../src/system.f90 parameters.o 
	$(FORTRAN_COMPILER) -c $< -o ../../build/$(notdir $@)

clean: 
	@/bin/rm -f *.mod *.o ../../src/*.mod ../../src/*.o ../../app/*.mod ../../app/*.o ../../test/*.mod ../../test/*.o4 *.gcda *.gcno	../../build/*.mod ../../build/*.o ../../build/*.gcda ../../build/*.gcno


# Make macro cheat sheet
# 
# xxx: code.1 code.2
#		$@ = xxx
#		$< = code.1
#		$^ = code.1 code.2