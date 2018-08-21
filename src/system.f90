module system
    use iso_fortran_env, only: real64, int32
    implicit none
    private 

    real (real64), public, parameter :: system_size_x = 10.0
    real (real64), public, parameter :: system_size_y = 10.0
    real (real64), public, parameter :: system_size_z = 10.0

    real (real64), public, parameter, dimension(3) :: system_size =            &
                                                            [system_size_x,    &
                                                             system_size_y,    &
                                                             system_size_z]

    integer (int32), public, parameter :: number_of_particles = int(1e3)

contains
end module system
