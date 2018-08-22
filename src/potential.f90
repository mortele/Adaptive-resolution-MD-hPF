module potential
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use particles,  only: masses
    implicit none
    private 

    public :: compute_forces

contains

    subroutine compute_forces(positions, forces)
        real (real64), dimension(:), intent(in)     :: positions
        real (real64), dimension(:), intent(in out) :: forces
        forces = 0
    end subroutine compute_forces
end module potential