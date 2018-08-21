GFORTRAN 			= gfortran
COMPILER_FLAGS 		= -fimplicit-none -fmodule-private -Wall -std=f2008ts
FORTRAN_COMPIILER 	= $(GFORTRAN) $(COMPILER_FLAGS)
PROGRAM 			= AdapResoMD-hPF

SRC = 	app/main.f90 	\
		src/system.f90 	\

OBJ = 	main.o 		\
		system.o 	\

MOD = 	system.mod 	\


all: $(OBJ)
	$(FORTRAN_COMPIILER) $(OBJ) -o app/$(PROGRAM).app

main.o: app/main.f90 $(MOD)
	$(FORTRAN_COMPIILER) -c $<

%.o %.mod: src/%.f90
	$(FORTRAN_COMPIILER) -c $^

clean: 
	@/bin/rm -f *.mod *.o src/*.mod src/*.o app/*.mod app/*.o


