module particles
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private 

    integer (int32), public, dimension(:), allocatable :: types
    real   (real64), public, dimension(:), allocatable :: positions,  &
                                                          velocities, &
                                                          forces,     &
                                                          masses



contains
end module particles