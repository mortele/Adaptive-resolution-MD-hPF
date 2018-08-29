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

include depend.mk

SRC = $(wildcard src/*.f90)
OBJ = $(patsubst src/%.f90, %.o,   $(SRC))
MOD = $(patsubst src/%.f90, %.mod, $(SRC))

TEST_SRC = $(wildcard test/*.f90)
TEST_OBJ = $(patsubst test/%.f90, %.o,   $(SRC))
TEST_MOD = $(patsubst test/%.f90, %.mod, $(SRC))

all: $(OBJ)
	$(FORTRAN_COMPILER) $(OBJ) $(MOD) -o $(PROGRAM) 

%.o %.mod: src/%.f90
	$@
	#$(FORTRAN_COMPILER) -J. -c $< -o $@

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