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
TEST_FLAGS 					= -fprofile-arcs -ftest-coverage
FORTRAN_COMPIILER 			= $(GFORTRAN) $(COMPILER_FLAGS) $(OPTIMIZATION_FLAGS) $(DEBUG_FLAGS) $(TEST_FLAGS)
FORTRAN_COMPIILER_OPTIMIZE 	= $(GFORTRAN) $(COMPILER_FLAGS) $(AGRESSIVE_OPTIMIZE) 
PROGRAM 					= AdapResoMD-hPF
TEST_PROGRAM 				= UnitTests
FRUIT_SRC 					= ext/FRUIT/src/fruit.f90
FRUIT_MOD 					= fruit.mod fruit_util.mod


SRC = $(wildcard src/*.f90)
OBJ = $(patsubst src/%.f90, build/%.o,   $(SRC))
MOD = $(patsubst src/%.f90, build/%.mod, $(SRC))

# the names of every source file, separated by space;
# you can use a wildcard, or specify the explicitly
F90SRC=$(SRC)

# directory in which the object files are located (optional);
# it must end with a slash
OBJPREFIX=build/

# note: it's probably a bad idea to put this snippet
#       before the first target rule in the makefile

# extract the dependencies of the project
# (the '@' before 'for' silences the output of this command)
#
# note: this is kind of a hack, as all it does is to search the
#       sources for a line that looks like "use <module>"
dependencies.mk: $(F90SRC)
	@for f in $^; do \
	    printf "%s:" "$(OBJPREFIX)$${f%.f90}.o"; \
	    awk -v p="$(OBJPREFIX)" \
	        '$$1 == "use" && NF == 2 { printf " %s%s.o",p,$$2 }' "$$f"; \
	    echo; \
	done >$@.tmp; \
	mv $@.tmp $@

include dependencies.mk

all: $(OBJ)
  $(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $(OBJ) -o prog 

obj/%.o: src/%.cpp
  $(CC) $(CFLAGS) -c $< -o $@
