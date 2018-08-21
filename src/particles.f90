module particle
    use iso_fortran_env, only: real64, int32
    implicit none
    private 

    real (real64), public, dimension(:), allocatable :: positions,  &
                                                        velocities, &
                                                        forces
    real (real64), public, dimension(:), allocatable :: masses


contains
end module particle