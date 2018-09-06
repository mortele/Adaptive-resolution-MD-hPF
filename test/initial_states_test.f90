module initial_states_test
    use, intrinsic :: iso_fortran_env, only: real64, int32, output_unit
    use fruit,          only:   assert_equals,                      &
                                assert_true,                        &
                                assert_false
    use system,         only:   system_size,                        &
                                apply_periodic_boundary_conditions    ! Subroutine
    use particles,      only:   positions,                          &
                                velocities,                         &
                                forces,                             &
                                masses,                             &
                                types
    use initial_states, only:   random_initial_state,               & ! Subroutine
                                fcc_initial_state,                  & ! Subroutine
                                setup_initial_state,                & ! Subroutine
                                allocate_arrays                       ! Subroutine
    use parameters,     only:   number_of_dimensions,               &
                                number_of_particles,                &
                                fcc_lattice_constant,               &
                                fcc_number_of_unit_cells,           &
                                temperature,                        &
                                system_size_x,                      &
                                system_size_y,                      &
                                system_size_z,                      &
                                initial_configuration,              &
                                lennard_jones_epsilon,              &
                                lennard_jones_sigma
    use potential,      only:   compute_forces,                     & ! Subroutine
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
        
        if (allocated(positions)) then
            deallocate(positions)
        end if
        if (allocated(velocities)) then
            deallocate(velocities)
        end if
        if (allocated(forces)) then
            deallocate(forces)
        end if
        if (allocated(masses)) then
            deallocate(masses)
        end if
        if (allocated(types)) then
            deallocate(types)
        end if
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
        call assert_false(all(abs(positions - 0.0_real64) < 1e-10_real64), "4  test_random_initial_state : The random initial state did not position particles randomly (all particles left in [0,0,0]")
        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(types)

        ! Test the *default to random* if the initial_configuration string isnt
        ! recognized.
        number_of_dimensions  = 3
        number_of_particles   = 10
        initial_configuration = "whatev"
        call setup_initial_state(positions, velocities, forces, masses, types, silent)
        call assert_true(allocated(positions),                       "5  test_random_initial_state : Arrays not allocated by random_initial_state()")
        call assert_equals(number_of_dimensions, size(positions, 1), "6  test_random_initial_state : Arrays not allocated to correct size random_initial_state()")
        call assert_equals(number_of_particles,  size(positions, 2), "7  test_random_initial_state : Arrays not allocated to correct size random_initial_state()")
        call assert_false(all(abs(positions - 0.0_real64) < 1e-10_real64), "8  test_random_initial_state : The random initial state did not position particles randomly (all particles left in [0,0,0]")

        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(types)
    end subroutine test_random_initial_state


    subroutine test_fcc_initial_state()
        integer         :: i, j
        real (real64)   :: FCC_potential_energy, potential_energy
        real (real64)   :: dx = 1e-5_real64
        logical :: silent = .true.

        lennard_jones_epsilon    = 1.0
        lennard_jones_sigma      = 3.405
        number_of_dimensions     = 3
        fcc_lattice_constant     = dsqrt(2.0_real64) * (2.0_real64)**(1.0_real64 / 6.0_real64) * lennard_jones_sigma
        fcc_number_of_unit_cells = 1
        initial_configuration    = "fcc"
        call setup_initial_state(positions, velocities, forces, masses, types, silent)

        call assert_equals(4*fcc_number_of_unit_cells**3, number_of_particles, "1  test_fcc_initial_state : Number of particles not correctly assigned in test_fcc_initial_state")

        ! Test that the FCC lattice is the minium of the potential energy. We 
        ! move all particles in every direction and check if the potential 
        ! energy ever decreases in any other configuration.
        call compute_forces(positions, forces)
        FCC_potential_energy = V
        do i = 1, number_of_particles
            do j = 1, number_of_dimensions
                
                positions(j,i) = positions(j,i) + dx
                call apply_periodic_boundary_conditions(positions)
                call compute_forces(positions, forces)
                potential_energy = V
                positions(j,i) = positions(j,i) - dx
                call apply_periodic_boundary_conditions(positions)
                call assert_true(FCC_potential_energy < potential_energy, char(i)//" test_fcc_initial_state : The FCC state is not the minimum of the potential energy")
            end do
        end do

        ! If the FCC lattice is the minimum energy configuration, then all 
        ! forces on every atom should vanish in this configuration.
        call compute_forces(positions, forces)
        do i = 1, number_of_particles
            call assert_equals([0.0_real64, 0.0_real64, 0.0_real64], forces(:,i), 3, 1e-14_real64, char(i)//" test_fcc_initial_state : Particle net force should be zero in the FCC lattice")
        end do

        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(types)
        deallocate(masses)
    end subroutine test_fcc_initial_state

end module initial_states_test