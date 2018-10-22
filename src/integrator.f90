module integrator
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use potential,      only: compute_forces                      ! Subroutine
    use particles,      only: masses
    use system,         only: apply_periodic_boundary_conditions  ! Subroutine
    use parameters,     only: time_step,                          &
                              number_of_particles,                &
                              number_of_dimensions
    implicit none
    private

    public ::   integrate_one_step,         &
                half_move,                  &
                move

    ! To overload the integrate_one_step function in order to facilitate 
    ! testing with different force calculations, we write an explicit interface
    ! using the almost identical subroutines. 
    ! 
    ! TODO: Figure out if this is horrible, awful, godforbidden bad practice to
    !       basically just copy paste 95% of the code in these subroutines. Is 
    !       there a cleaner way to do this?
    interface integrate_one_step
        ! Takes arguments (position, velocities, forces) 
        module procedure integrate_one_step_default_force
        ! Takes arguments (position, velocities, forces, force_calculator) 
        module procedure integrate_one_step_argument_force 
    end interface integrate_one_step
contains

    ! Integrate a single step using the velocity Verlet algorithm. Uses the 
    ! default force calculator (potential module :: compute_forces).
    subroutine integrate_one_step_default_force(positions, velocities, forces)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: positions
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:,:), intent(in out) :: forces

        logical, save :: first_step = .true.

        if (first_step) then
            forces = 0.0_real64
            call compute_forces(positions, forces)
            first_step = .false.
        end if
        
        call half_move     (velocities, forces)
        call move          (positions,  velocities)

        call apply_periodic_boundary_conditions(positions)
        
        call compute_forces(positions,  forces)
        call half_move     (velocities, forces)
    end subroutine integrate_one_step_default_force


    ! Integrate a single step using the velocity Verlet algorithm. Uses the 
    ! specified force calculator given as argument.
    subroutine integrate_one_step_argument_force(positions, velocities, forces, force_calculator)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: positions
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:,:), intent(in out) :: forces
        interface 
            subroutine force_calculator(pos, for)
                use, intrinsic :: iso_fortran_env, only: real64
                real (real64), dimension(:,:), intent(in)     :: pos
                real (real64), dimension(:,:), intent(in out) :: for
            end subroutine
        end interface
        
        logical, save :: first_step = .true.

        if (first_step) then
            forces = 0.0_real64
            call force_calculator(positions, forces)
            first_step = .false.
        end if
        
        call half_move     (velocities, forces)
        call move          (positions,  velocities)

        call apply_periodic_boundary_conditions(positions)
        
        call force_calculator(positions,  forces)
        call half_move       (velocities, forces)
    end subroutine integrate_one_step_argument_force

    subroutine half_move(velocities, forces)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:,:), intent(in)     :: forces
        integer (int32) :: i
        
        do i = 1, number_of_particles
            velocities(:, i) = velocities(:, i) + 0.5_real64 * time_step / masses(i) * forces(:, i)
        end do
    end subroutine half_move

    subroutine move(positions, velocities)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: positions
        real (real64), dimension(:,:), intent(in)     :: velocities

        integer (int32) :: i
        do i = 1, number_of_particles
            positions(:, i) = positions(:, i) + velocities(:, i) * time_step
        end do
    end subroutine move
end module integrator
