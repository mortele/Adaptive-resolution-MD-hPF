module parameters
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private 

    real (real64), public, parameter :: system_size_x = 10.0
    real (real64), public, parameter :: system_size_y = 10.0
    real (real64), public, parameter :: system_size_z = 10.0

    integer (int32), public, parameter :: number_of_particles  = int(1e3)
    integer (int32), public, parameter :: number_of_dimensions = 3
    
    ! Number of field vertices per dimension, total number of vertices is 
    ! field_nodes^3.
    integer (int32), public, parameter :: field_nodes = 25
    real (real64), public, parameter :: temperature = 1.0
    character (*), public, parameter :: out_file_name = "positions.xyz"

    ! Time step used in the velocity Verlet integration.
    real (real64), public, parameter :: time_step = 0.01

    ! Step length used in the numerical derivatives of the density field.
    real (real64), public, parameter :: step_length = 0.01

    integer (int32), public, parameter :: number_of_time_steps = int(1e3)




contains

end module parameters