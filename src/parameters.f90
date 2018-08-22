module parameters
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private 

    real (real64), public, parameter :: system_size_x = 10.0
    real (real64), public, parameter :: system_size_y = 10.0
    real (real64), public, parameter :: system_size_z = 10.0

    integer (int32), public, parameter :: number_of_particles = int(1e3)
    
    ! Number of field vertices per dimension, total number of vertices is 
    ! field_nodes^3.
    integer (int32), public, parameter :: field_nodes         = 25



contains

end module parameters