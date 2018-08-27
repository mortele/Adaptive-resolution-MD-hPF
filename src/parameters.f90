module parameters
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    private 

    ! Pre-set fabricated initial states. 
    ! Current options: 
    !    - "random"         ! Random uniform positions, velocities according to 
    !                       ! Maxwell-Boltzmann distribution at given 
    !                       ! temperature. Default case used.
    !    - "fcc"            ! Solid Ar configuration, face centered cubic 
    !                       ! lattice of given lattice constant, b. Velocities
    !                       ! According to Maxwell-Boltzmann distribution at 
    !                       ! given temperature.
    character (*), public, parameter :: initial_configuration = "fcc"

    ! FCC parameters. 
    real (real64), public, parameter :: fcc_lattice_constant = 5.26
    integer (int32), public, parameter :: fcc_number_of_unit_cells = 5

    real (real64), public :: system_size_x = 10.0
    real (real64), public :: system_size_y = 10.0
    real (real64), public :: system_size_z = 10.0

    integer (int32), public            :: number_of_particles  = int(1e2)
    integer (int32), public, parameter :: number_of_dimensions = 3
    
    ! Number of field vertices per dimension, total number of vertices is 
    ! field_nodes^3.
    integer (int32), public, parameter :: field_nodes = 25
    real (real64),   public, parameter :: temperature = 1.0

    character (*), public, parameter :: out_file_name  = "positions.xyz"
    character (*), public, parameter :: info_file_name = "systeminfo.out"

    ! Time step used in the velocity Verlet integration.
    real (real64), public, parameter :: time_step = 0.001

    ! Step length used in the numerical derivatives of the density field.
    real (real64), public, parameter :: step_length = 0.01

    integer (int32), public, parameter :: number_of_time_steps = 3

    ! Potential parameters.
    real (real64), public, parameter :: lennard_jones_epsilon = 1.0
    real (real64), public, parameter :: lennard_jones_sigma   = 3.405
    real (real64), public, parameter :: lennard_jones_cutoff  = 100.0


contains

end module parameters