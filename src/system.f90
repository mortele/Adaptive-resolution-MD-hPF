module system
    use, intrinsic :: iso_fortran_env, only:  real64, int32
    use parameters,      only:  system_size_x, &
                                system_size_y, &
                                system_size_z
    implicit none
    private 

    real (real64), public, parameter, dimension(3) :: system_size =            &
                                                            [system_size_x,    &
                                                             system_size_y,    &
                                                             system_size_z]




contains
end module system
