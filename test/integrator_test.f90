module integrator_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,              &   ! Subroutine
                                    assert_true                     ! Subroutine
    use integrator,         only:   integrate_one_step              ! Subroutine
    use potential,          only:   compute_forces                  ! Subroutine
    use particles,          only:   masses,                     &
                                    positions,                  &
                                    velocities,                 &
                                    forces
    use parameters,         only:   number_of_particles,        &
                                    number_of_dimensions,       &
                                    number_of_time_steps,       &
                                    lennard_jones_sigma,        &
                                    lennard_jones_epsilon,      &
                                    time_step
    use system,             only:   apply_periodic_boundary_conditions
    implicit none
    
    private

    public  ::  setup,                      &
                teardown,                   &
                test_integrate_one_step,    &
                test_half_move,             &
                test_move
contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_half_move()
    end subroutine test_half_move

    subroutine test_move()
    end subroutine test_move

    subroutine test_integrate_one_step()
    end subroutine test_integrate_one_step

end module integrator_test