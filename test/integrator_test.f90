module integrator_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,              &   ! Subroutine
                                    assert_true                     ! Subroutine
    use integrator,         only:   integrate_one_step,         &   ! Subroutine
                                    half_move,                  &   ! Subroutine
                                    move                            ! Subroutine
    use potential,          only:   compute_forces                  ! Subroutine
    use particles,          only:   masses,                     &
                                    positions,                  &
                                    velocities,                 &
                                    forces,                     &
                                    types
    use parameters,         only:   number_of_particles,        &
                                    number_of_dimensions,       &
                                    number_of_time_steps,       &
                                    lennard_jones_sigma,        &
                                    lennard_jones_epsilon,      &
                                    time_step
    use system,             only:   apply_periodic_boundary_conditions ! Subroutine
    implicit none
    
    private ::  zero_force_calculator

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
        real (real64), allocatable, dimension(:,:) :: old_velocities
        integer (int32) :: i

        ! 1D tests
        number_of_dimensions  = 1
        number_of_particles   = 1
        time_step             = 0.01
        allocate(velocities (number_of_dimensions, number_of_particles))
        allocate(forces     (number_of_dimensions, number_of_particles))
        allocate(masses     (number_of_particles))
        velocities  = 0.0
        forces      = 0.0
        masses      = 97.2982149

        call half_move(velocities, forces)
        call assert_equals(0.0_real64, velocities(1,1), "1  test_half_move : A half move with no forces shouldnt change the velocity")
        
        forces = 214.135985
        call half_move(velocities, forces)

        call assert_equals(0.5_real64 * forces(1,1) * time_step / masses(1), velocities(1,1), "2  test_half_move : The half move didnt change the particles velocity by the correct amount")

        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)

        ! 3D tests
        number_of_dimensions  = 3
        number_of_particles   = 10
        time_step             = 0.01
        allocate(old_velocities (number_of_dimensions, number_of_particles))
        allocate(velocities     (number_of_dimensions, number_of_particles))
        allocate(forces         (number_of_dimensions, number_of_particles))
        allocate(masses         (number_of_particles))
        velocities  = 0.0
        forces      = 0.0
        
        call random_number(velocities)
        call random_number(forces)
        call random_number(masses)
        velocities  = velocities * 20.0 - 10.0 ! Uniform random in [-10, 10)
        forces      = forces     * 20.0 - 10.0 ! Uniform random in [-10, 10)
        masses      = masses     * 50.0 + 1.0  ! Uniform random in [1, 51)
        old_velocities = velocities
        
        call half_move(velocities, forces)

        do i = 1, number_of_particles
            call assert_equals(old_velocities(:,i) + 0.5_real64 * forces(:,i) * time_step / masses(i), velocities(:,i), 3, "10  test_half_move : The half move didnt change the particles velocity by the correct amount")
        end do

        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)
        deallocate(old_velocities) 
    end subroutine test_half_move

    subroutine test_move()
        real (real64), allocatable, dimension(:,:) :: old_positions
        integer (int32) :: i

        ! 1D tests
        number_of_dimensions  = 1
        number_of_particles   = 1
        time_step             = 0.01
        allocate(positions  (number_of_dimensions, number_of_particles))
        allocate(velocities (number_of_dimensions, number_of_particles))
        positions   = 0.0
        velocities  = 0.0

        call move(positions, velocities)
        call assert_equals(0.0_real64, positions(1,1), "1  test_move : A move with no velocity shouldnt change the position")
        
        velocities = 314.2849502
        call move(positions, velocities)

        call assert_equals(velocities(1,1) * time_step, positions(1,1), "2  test_move : The move didnt change the particles position by the correct amount")

        deallocate(positions)
        deallocate(velocities)

                ! 3D tests
        number_of_dimensions  = 3
        number_of_particles   = 10
        time_step             = 0.01
        allocate(old_positions  (number_of_dimensions, number_of_particles))
        allocate(positions      (number_of_dimensions, number_of_particles))
        allocate(velocities     (number_of_dimensions, number_of_particles))
        positions   = 0.0
        velocities  = 0.0
        
        call random_number(positions)
        call random_number(velocities)
        positions   = positions  * 10.0        ! Uniform random in [0, 10)
        velocities  = velocities * 20.0 - 10.0 ! Uniform random in [-10, 10)
        old_positions = positions
        
        call move(positions, velocities)

        do i = 1, number_of_particles
            call assert_equals(old_positions(:,i) + velocities(:,i) * time_step, positions(:,i), 3, "10  test_move : The move didnt change the particles position by the correct amount")
        end do

        deallocate(positions)
        deallocate(velocities)
        deallocate(old_positions) 
    end subroutine test_move

    subroutine test_integrate_one_step()
        ! First we test with no forces, one particle, in 1D.
        number_of_dimensions  = 1
        number_of_particles   = 1
        time_step             = 0.01
        allocate(positions  (number_of_dimensions, number_of_particles))
        allocate(velocities (number_of_dimensions, number_of_particles))
        allocate(forces     (number_of_dimensions, number_of_particles))
        allocate(masses     (number_of_particles))

        positions   = 1.0
        velocities  = 0.0
        forces      = 0.0
        masses      = 1.0
        call integrate_one_step(positions, velocities, forces, zero_force_calculator)
        
        call assert_equals(1.0_real64, positions(1,1),  "1  test_integrate_one_step : The one step integration shouldnt change the position of a particle when the velocities and forces are both zero")
        call assert_equals(0.0_real64, velocities(1,1), "2  test_integrate_one_step : The one step integration shouldnt change the velocity of a particle when the forces are set to zero")


        lennard_jones_sigma     = 3.405
        lennard_jones_epsilon   = 1.0
    end subroutine test_integrate_one_step

    subroutine zero_force_calculator(positions, forces)
        implicit none
        real (real64), dimension(:,:), intent(in)     :: positions
        real (real64), dimension(:,:), intent(in out) :: forces

        ! forces = 0.0

        ! To make the compiler with -Wall -Wextra -pedantic-errors not complain
        ! about unused parameter positions.
        forces = 0.0 * positions 
    end subroutine
end module integrator_test