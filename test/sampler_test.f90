module sampler_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,          only:   assert_equals,                      &
                                assert_true
    use parameters,     only:   number_of_time_steps,               &
                                number_of_particles,                &
                                number_of_dimensions
    use particles,      only:   positions,                          &
                                velocities,                         &
                                forces,                             &
                                masses,                             &
                                types
    use sampler,        only:   kinetic_energy,                     &
                                potential_energy,                   &
                                total_energy,                       &
                                store_energy                          ! Subroutine
    use potential,      only:   Ek, V, E
    implicit none
    
    private 

    public ::   setup,                      &
                teardown,                   &
                test_store_energy
contains

    subroutine setup()
    end subroutine setup 

    subroutine teardown()
    end subroutine teardown

    subroutine test_store_energy()
        number_of_particles  = 5
        number_of_dimensions = 3
        number_of_time_steps = 10
        allocate(positions (number_of_dimensions, number_of_particles))
        allocate(velocities(number_of_dimensions, number_of_particles))
        allocate(forces    (number_of_dimensions, number_of_particles))
        allocate(masses    (number_of_particles))
        allocate(types     (number_of_particles))

        call random_number(positions)
        call random_number(velocities)
        forces = 0.0
        masses = 1.0
        types  = 1

        allocate(potential_energy(number_of_time_steps + 1))
        allocate(kinetic_energy  (number_of_time_steps + 1))
        allocate(total_energy    (number_of_time_steps + 1))
        call store_energy(kinetic_energy, potential_energy, total_energy)

        call assert_equals(Ek, kinetic_energy(1),   "1  test_store_energy : Stored kinetic energy does not equal computed kinetic energy in the potential module")
        call assert_equals(V,  potential_energy(1), "2  test_store_energy : Stored potential energy does not equal computed potential energy in the potential module")
        call assert_equals(E,  total_energy(1),     "3  test_store_energy : Stored total energy does not equal computed total energy in the potential module")

        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(types)
        deallocate(potential_energy)
        deallocate(kinetic_energy)
        deallocate(total_energy)
    end subroutine test_store_energy
end module sampler_test