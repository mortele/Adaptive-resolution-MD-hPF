module initial_states
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,           only: system_size
    use particles,        only: masses
    use random_generator, only: random_normal
    use parameters,       only: number_of_particles,  &
                                number_of_dimensions, &
                                temperature
    implicit none
    private 
    public :: random_initial_state

contains

    subroutine random_initial_state(positions, velocities)
        real (real64), dimension(:,:), intent(in out) :: positions
        real (real64), dimension(:,:), intent(in out) :: velocities
        integer (int32)  :: i, j
        real    (real64) :: T = temperature

        do i = 1, number_of_particles
            ! Assign a uniform random position to each particle, then scale each
            ! dimension with the system size to get a uniform distribution in 
            ! the complete simulation box.
            call random_number(positions(:, i))
            
            do j = 1, number_of_dimensions
                positions(j, i) = positions(j, i) * system_size(j)
            end do
        end do

        do i = 1, number_of_particles
            ! Assign the velocity of each particle to a normal distribution, 
            ! with zero mean and unit variance.
            call random_normal(velocities(:, i))
            
            ! Scale the standard deviation to be (kB T / m) according to the 
            ! Boltzmann distribution at temperature T.
            velocities(:, i) = velocities(:, i) * T / masses(i)
        end do
    end subroutine random_initial_state
end module initial_states