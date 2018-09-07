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
    use system,             only:   system_size,                &
                                    apply_periodic_boundary_conditions ! Subroutine
    implicit none
    
    real (real64), dimension(:), allocatable :: constant_force

    private ::  zero_force_calculator,      &
                constant_force_calculator

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
        real (real64), allocatable, dimension(:) :: initial_velocity
        real (real64), allocatable, dimension(:) :: initial_position
        real (real64), allocatable, dimension(:) :: expected_velocity
        real (real64), allocatable, dimension(:) :: expected_positions
        integer (int32) :: i

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

        ! Lets make sure that also the default force calculation produce the 
        ! same result if we arrange for compute_forces to give no forces.
        lennard_jones_sigma     = 3.405
        lennard_jones_epsilon   = 1.0
        call integrate_one_step(positions, velocities, forces) ! Default call

        call assert_equals(1.0_real64, positions(1,1),  "1  test_integrate_one_step : The one step integration shouldnt change the position of a particle when the velocities and forces are both zero")
        call assert_equals(0.0_real64, velocities(1,1), "2  test_integrate_one_step : The one step integration shouldnt change the velocity of a particle when the forces are set to zero")

        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)

        ! We next test a known case which we can compute by hand:
        !
        ! If we set the force constant, equal to F0 = (a0,a0,a0), then the 
        ! *velocity* obeys the recurrence relation
        !
        !    v       =   v      +  Δt  a  
        !      i+1         i             0
        ! 
        ! with solution
        ! 
        !    v       =   v      +  i Δt a
        !      i           0              0
        !
        ! The *position* obeys the recurrence relation 
        !                             ╭     ╮       1    2
        !    p       =   p      +  Δt │	v   │   +  ─── Δt  a
        !      i+1         i          ╰   i ╯       2        0
        !                             ╭                 ╮       1    2  
        !            =   p      +  Δt │	v   +  i Δt a   │   +  ─── Δt  a    
        !                  i          ╰   0           0 ╯       2        0 
        !
        ! with solution
        !               i Δt ╭                      ╮
        !    p       = ───── │ a   Δt i   +   2 v   │  + p
        !      i         2   ╰   0                0 ╯      0
        !
        !
        ! We test now if the velocity verlet algorithm satisfies this. Note 
        ! that this should be satisfied *to machine precision*, i.e. ~1e-15.
        number_of_particles  = 1
        number_of_dimensions = 3
        time_step            = 0.01
        system_size          = [10.0, 10.0, 10.0]
        allocate(constant_force   (number_of_dimensions))
        allocate(expected_velocity(number_of_dimensions))
        allocate(initial_velocity (number_of_dimensions))
        allocate(initial_position (number_of_dimensions))
        allocate(positions        (number_of_dimensions, number_of_particles))
        allocate(velocities       (number_of_dimensions, number_of_particles))
        allocate(forces           (number_of_dimensions, number_of_particles))
        constant_force  = [ 0.85_real64,  2.0_real64,  -1.98_real64]
        velocities(:,1) = [-3.3_real64,   0.0_real64,   1.5_real64]
        positions (:,1) = [ 5.0_real64,   5.0_real64,   5.0_real64]
        masses          = 1.0

        expected_velocity  = velocities(:,1)
        expected_positions = positions (:,1) 
        initial_velocity   = velocities(:,1)
        initial_position   = positions (:,1)

        call constant_force_calculator(positions, forces)
        do i = 1, 10
            call integrate_one_step(positions, velocities, forces, constant_force_calculator)
            expected_positions = 0.5_real64 * i * time_step * (constant_force * time_step * i + 2.0_real64 * initial_velocity) + initial_position
            expected_velocity  = initial_velocity + time_step * i * constant_force
            call assert_equals(expected_positions, positions(:,1),  3, 1e-14_real64, "20 test_integrate_one_step : Integrating one step with constant force does not reproduce the known closed form difference equation solution (position)")
            call assert_equals(expected_velocity,  velocities(:,1), 3, 1e-14_real64, "10 test_integrate_one_step : Integrating one step with constant force does not reproduce the known closed form difference equation solution (velocity)")
        end do


        

    end subroutine test_integrate_one_step

    subroutine zero_force_calculator(positions, forces)
        implicit none
        real (real64), dimension(:,:), intent(in)     :: positions
        real (real64), dimension(:,:), intent(in out) :: forces

        ! To make the compiler with -Wall -Wextra -pedantic-errors not complain
        ! about unused parameter positions.
        forces = 0.0 * positions 
    end subroutine

    subroutine constant_force_calculator(positions, forces)
        implicit none
        real (real64), dimension(:,:), intent(in)     :: positions
        real (real64), dimension(:,:), intent(in out) :: forces
        integer (int32) :: i

        do i = 1, number_of_particles
            forces(:,i) = constant_force
        end do

        ! To make the compiler with -Wall -Wextra -pedantic-errors not complain
        ! about unused parameter positions.
        forces = forces + 0.0 * positions 
    end subroutine constant_force_calculator
end module integrator_test