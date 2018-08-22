module initial_states
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,          only: system_size
    use parameters,      only: number_of_particles
    implicit none
    private 


contains

    subroutine random_initial_state(positions, velocities)
        real (real64), dimension(:), intent(in)  :: positions
        real (real64), dimension(:), intent(out) :: velocities
        
        
    end subroutine random_initial_state
end module initial_states