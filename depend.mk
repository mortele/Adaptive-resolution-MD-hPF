field.o field.mod : src/field.f90 system.mod parameters.mod 
fileWriter.o file_writer.mod : src/fileWriter.f90 parameters.mod system.mod 
initialStates.o initial_states.mod : src/initialStates.f90 parameters.mod system.mod random_generator.mod 
integrator.o integrator.mod : src/integrator.f90 parameters.mod system.mod particles.mod potential.mod 
main.o : src/main.f90 parameters.mod sampler.mod potential.mod integrator.mod file_writer.mod initial_states.mod particles.mod 
parameters.o parameters.mod : src/parameters.f90 
particles.o particles.mod : src/particles.f90 
potential.o potential.mod : src/potential.f90 parameters.mod particles.mod system.mod 
randomGenerator.o random_generator.mod : src/randomGenerator.f90 
sampler.o sampler.mod : src/sampler.f90 potential.mod particles.mod parameters.mod 
system.o system.mod : src/system.f90 parameters.mod 
system_test.o system_test.mod : test/system_test.f90 system.mod particles.mod parameters.mod 
unitTests.o : test/unitTests.f90 system_test.mod 
