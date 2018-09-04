module field_test
    use, intrinsic :: iso_fortran_env,  only: int32, real64
    use fruit,      only:   assert_equals,                       &
                            assert_true
    use system,     only:   system_size
    use parameters, only:   number_of_field_nodes,              &
                            number_of_dimensions,               &
                            number_of_particles,                &
                            system_size_x,                      &
                            system_size_y,                      &
                            system_size_z
    use field,      only:   allocate_field_arrays,              & ! Submodule
                            compute_density_field,              & ! Submodule
                            density_field,                      &
                            position_of_density_nodes,          &
                            density_gradient
    
    use particles,  only:   positions
    implicit none
    private 

    public  ::  setup,                          &
                teardown,                       &
                test_compute_density_field,     &
                test_allocate_field_arrays
contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_compute_density_field()
        integer (int32) :: i, j, k

        ! We start off with a 1D test of the density field computation.
        number_of_field_nodes = 10
        number_of_particles   = 10
        number_of_dimensions  = 3
        system_size_x         = 10.0_real64
        system_size_y         = 10.0_real64
        system_size_z         = 10.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        allocate(positions(number_of_dimensions, number_of_particles))

        do i = 1, number_of_particles
            positions(1,i) = i-1    ! A line of equidistant particles, separated 
                                    ! by a distance of 1.0.
            positions(2:,i) = [0.0_real64, 0.0_real64]
        end do


        call compute_density_field(positions)

        
    end subroutine test_compute_density_field

    subroutine test_allocate_field_arrays()

    end subroutine test_allocate_field_arrays

end module field_test