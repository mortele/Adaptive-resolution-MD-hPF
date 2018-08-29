field.o: src/field.f90 system.o parameters.o 
fileWriter.o: src/fileWriter.f90 parameters.o system.o 
initialStates.o: src/initialStates.f90 parameters.o system.o random_generator.o 
integrator.o: src/integrator.f90 parameters.o system.o particles.o potential.o 
main.o: src/main.f90 parameters.o sampler.o potential.o integrator.o file_writer.o initial_states.o particles.o
parameters.o: src/parameters.f90 
particles.o: src/particles.f90 
potential.o: src/potential.f90 parameters.o particles.o system.o 
randomGenerator.o: src/randomGenerator.f90 
sampler.o: src/sampler.f90 potential.o particles.o parameters.o 
system.o: src/system.f90 parameters.o 
system_test.o: test/system_test.f90 system.o particles.o parameters.o 
unitTests.o: unitTests.f90 system_test.o 
