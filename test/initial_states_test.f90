module initial_states_test
    use, intrinsic :: iso_fortran_env, only: real64, int32, output_unit
    use fruit,          only:   assert_equals,                  &
                                assert_true
    use system,         only:   system_size
    use particles,      only:   positions,                      &
                                velocities,                     &
                                forces,                         &
                                masses,                         &
                                types
    use initial_states, only:   random_initial_state,           & ! Subroutine
                                fcc_initial_state,              & ! Subroutine
                                setup_initial_state,            & ! Subroutine
                                allocate_arrays                   ! Subroutine
    use parameters,     only:   number_of_dimensions,           &
                                number_of_particles,            &
                                fcc_lattice_constant,           &
                                fcc_number_of_unit_cells,       &
                                temperature,                    &
                                system_size_x,                  &
                                system_size_y,                  &
                                system_size_z,                  &
                                initial_configuration,          &
                                lennard_jones_epsilon,          &
                                lennard_jones_sigma
    use potential,      only:   compute_forces,                 & ! Subroutine
                                V


    implicit none
    private 

    public ::   setup,                      &
                teardown,                   &
                test_allocate_arrays,       &
                test_random_initial_state,  &
                test_fcc_initial_state

contains

    subroutine setup()
    end subroutine setup 

    subroutine teardown()
    end subroutine teardown

    subroutine test_allocate_arrays()
        number_of_dimensions = 2
        number_of_particles  = 8
        call allocate_arrays(positions, velocities, forces, masses, types)
        call assert_true(allocated(positions),    "1  test_allocate_arrays : Positions was not allocated by allocate_arrays()")
        call assert_true(allocated(velocities),   "2  test_allocate_arrays : Velocities was not allocated by allocate_arrays()")
        call assert_true(allocated(forces),       "3  test_allocate_arrays : Forces was not allocated by allocate_arrays()")
        call assert_true(allocated(masses),       "4  test_allocate_arrays : Masses was not allocated by allocate_arrays()")
        call assert_true(allocated(types),        "5  test_allocate_arrays : Types was not allocated by allocate_arrays()")
        call assert_equals(2, size(positions, 1), "6  test_allocate_arrays : positions array allocated to wrong size")
        call assert_equals(8, size(positions, 2), "7  test_allocate_arrays : positions array allocated to wrong size")
        call assert_equals(2, size(velocities,1), "8  test_allocate_arrays : velocities array allocated to wrong size")
        call assert_equals(8, size(velocities,2), "9 test_allocate_arrays : velocities array allocated to wrong size")
        call assert_equals(2, size(forces,    1), "10 test_allocate_arrays : forces array allocated to wrong size")
        call assert_equals(8, size(forces,    2), "11 test_allocate_arrays : forces array allocated to wrong size")
        call assert_equals(8, size(masses),       "12 test_allocate_arrays : masses array allocated to wrong size")
        call assert_equals(8, size(types),        "13 test_allocate_arrays : types array allocated to wrong size")
        
        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(types)
    end subroutine test_allocate_arrays


    subroutine test_random_initial_state()
        logical :: silent = .true.

        number_of_dimensions  = 3
        number_of_particles   = 10
        initial_configuration = "random"
        call setup_initial_state(positions, velocities, forces, masses, types, silent)
        call assert_true(allocated(positions),                       "1  test_random_initial_state : Arrays not allocated by random_initial_state()")
        call assert_equals(number_of_dimensions, size(positions, 1), "2  test_random_initial_state : Arrays not allocated to correct size random_initial_state()")
        call assert_equals(number_of_particles,  size(positions, 2), "3  test_random_initial_state : Arrays not allocated to correct size random_initial_state()")

        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(types)
    end subroutine test_random_initial_state


    subroutine test_fcc_initial_state()
        integer (int32) :: i, j
        real (real64)   :: FCC_potential_energy, potential_energy
        real (real64)   :: dx = 0.02
        logical :: silent = .true.

        fcc_lattice_constant     = 5.26
        lennard_jones_epsilon    = 1.0
        lennard_jones_sigma      = 3.405
        number_of_dimensions     = 3
        fcc_number_of_unit_cells = 2
        initial_configuration    = "fcc"
        call setup_initial_state(positions, velocities, forces, masses, types, silent)

        call assert_equals(4*2**3, number_of_particles, "1  test_fcc_initial_state : Number of particles not correctly assigned in test_fcc_initial_state")

        ! Test that the FCC lattice is the minium of the potential energy. We 
        ! move all particles in every direction and check if the potential 
        ! energy ever decreases in any other configuration.
        call compute_forces(positions, forces)
        FCC_potential_energy = V

        do i = 1, number_of_particles
            do j = 1, number_of_dimensions
                
                positions(j,i) = positions(j,i) + dx
                call compute_forces(positions, forces)
                potential_energy = V
                positions(j,i) = positions(j,i) - dx
                ! write (output_unit, fmt="( f15.10, f15.10 )") potential_energy, FCC_potential_energy
                call assert_true(FCC_potential_energy < potential_energy, char(i*number_of_dimensions+i)//" test_fcc_initial_state : The FCC state is not the minimum of the potential energy")
            end do
        end do
    end subroutine test_fcc_initial_state

end module initial_states_test