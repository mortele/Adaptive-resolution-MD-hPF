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
        !integer (int32) :: i, j, k

        ! We start off with a 1D test of the density field computation, with a 
        ! single particle.
        number_of_field_nodes = 3
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 3.0_real64
        system_size_y         = 3.0_real64
        system_size_z         = 3.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        allocate(positions(number_of_dimensions, number_of_particles))
        positions(:,1) = [0.5,0.0,0.0]
        call compute_density_field(positions)
        call assert_equals(0.5_real64, density_field(1,1,1), "1  test_compute_density_field : A particle placed in the middle between field nodes 1 and 2 should distribute its mass evenly across the nodes")
        call assert_equals(0.5_real64, density_field(2,1,1), "2  test_compute_density_field : A particle placed in the middle between field nodes 1 and 2 should distribute its mass evenly across the nodes")
        
        
    end subroutine test_compute_density_field

    subroutine test_allocate_field_arrays()

    end subroutine test_allocate_field_arrays

end module field_test