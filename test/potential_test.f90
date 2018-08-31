module potential_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use potential,  only:   lennard_jones_force,        &
                            lennard_jones_potential,    &
                            compute_forces
    use fruit,      only:   assert_equals,              & ! Subroutine
                            assert_false                  ! Subroutine
    use parameters, only:   lennard_jones_sigma,        &
                            lennard_jones_epsilon,      &
                            lennard_jones_cutoff,       &
                            number_of_particles,        &
                            number_of_dimensions
    use system,     only:   system_size
    implicit none
    private

    public ::   setup,                          &
                teardown,                       &
                test_lennard_jones_force,       &
                test_lennard_jones_potential


contains
    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_lennard_jones_force
    end subroutine test_lennard_jones_force

    subroutine test_lennard_jones_potential
    end subroutine test_lennard_jones_potential

end module potential_test