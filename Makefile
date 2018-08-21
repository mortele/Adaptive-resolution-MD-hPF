GFORTRAN 			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -std=f2008ts
FORTRAN_COMPIILER 	= $(GFORTRAN) $(COMPILER_FLAGS)

SRC = 	main.f90 	\
		system.f90 	\

OBJ = 	main.o 		\
		system.o 	\

MOD = 	system.mod 	\


all: $(OBJ)
	$(FORTRAN_COMPIILER) $(OBJ) -o main.app

main.o: main.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

%.o %.mod: %.f90
	$(FORTRAN_COMPIILER) -c $^

clean: 
	@/bin/rm -f *.mod *.o


