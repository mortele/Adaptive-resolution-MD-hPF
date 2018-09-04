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
    character (6), public :: initial_configuration = "random"
    !character (3), public :: initial_configuration = "fcc"

    ! FCC parameters. 
    real (real64), public :: fcc_lattice_constant = 5.26
    integer (int32), public :: fcc_number_of_unit_cells = 5

    real (real64), public :: system_size_x = 10.0
    real (real64), public :: system_size_y = 10.0
    real (real64), public :: system_size_z = 10.0

    integer (int32), public :: number_of_particles  = int(1e1)
    integer (int32), public :: number_of_dimensions = 3
    
    ! Number of field vertices per dimension, total number of vertices is 
    ! field_nodes^3.
    integer (int32), public :: number_of_field_nodes = 25
    real (real64),   public :: temperature = 1.0

    character (13),  public :: out_file_name      = "positions.xyz"
    character (14),  public :: info_file_name     = "systeminfo.out"
    character (16),  public :: silent_output_file = "silentoutput.out"
    integer (int32), public :: silent_output_ID

    ! Time step used in the velocity Verlet integration.
    real (real64), public :: time_step = 0.001

    ! Step length used in the numerical derivatives of the density field.
    real (real64), public :: step_length = 0.01

    integer (int32), public :: number_of_time_steps = 3

    ! Potential parameters.
    real (real64), public :: lennard_jones_epsilon = 1.0
    real (real64), public :: lennard_jones_sigma   = 3.405
    real (real64), public :: lennard_jones_cutoff  = 100.0


contains

end module parameters