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

    public :: integrate_one_step

contains

    ! Integrate a single step using the velocity Verlet algorithm.
    subroutine integrate_one_step(positions, velocities, forces)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: positions
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:,:), intent(in out) :: forces
        logical, save :: first_step = .true.

        if (first_step) then
            forces = 0
            call compute_forces(positions, forces)
            first_step = .false.
        end if
        
        call half_move     (velocities, forces)
        call move          (positions,  velocities)
        call compute_forces(positions,  forces)
        call half_move     (velocities, forces)
    end subroutine integrate_one_step

    subroutine half_move(velocities, forces)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:,:), intent(in)     :: forces
        integer (int32) :: i
        
        do i = 1, number_of_particles
            velocities(:, i) = velocities(:, i) + 0.5 * time_step / masses(i) * forces(:, i)
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
