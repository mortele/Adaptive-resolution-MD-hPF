module system_test
    use, intrinsic :: iso_fortran_env, only: real64
    use fruit
    use parameters, only:   system_size_x,                      &
                            system_size_y,                      &
                            system_size_z,                      &
                            number_of_particles,                &
                            number_of_dimensions               
    use particles,  only:   positions,                          &
                            velocities,                         &
                            masses
    use system,     only:   remove_linear_momentum,             & ! Subroutine
                            apply_periodic_boundary_conditions, & ! Subroutine
                            compute_distance_minimum_image,     & ! Subroutine
                            system_size
    implicit none
    private

    public ::   test_periodic_boundary_conditions,      &
                test_distance_minimum_image
contains

    subroutine test_periodic_boundary_conditions()
        real (real64), parameter :: test_tollerance = 1e-10
        system_size_x = 10.0
        system_size_y = 10.0
        system_size_z = 10.0
        system_size = [system_size_x, system_size_y, system_size_z]
        number_of_dimensions = 3
        number_of_particles  = 11
        
        allocate(positions(number_of_dimensions, number_of_particles))
        positions(:,1)  = [1.0,   2.0,   3.0]
        positions(:,2)  = [11.0,  4.0,   5.0]
        positions(:,3)  = [12.0,  13.0,  6.0]
        positions(:,4)  = [7.0,   14.0, 15.0]
        positions(:,5)  = [16.0,  17.0, 18.0]
        positions(:,6)  = [-1.0,  8.0,  19.0]
        positions(:,7)  = [-2.0,  11.5,  9.0]
        positions(:,8)  = [12.5,  13.5, -3.0]
        positions(:,9)  = [-4.0,  -5.0, -6.0]
        positions(:,10) = [-7.0,  14.5, -8.0]
        positions(:,11) = [-9.0,  1.5,   2.5]
        call apply_periodic_boundary_conditions(positions)

        call assert_equals(positions(:,1),  [1.0_real64,   2.0_real64,   3.0_real64], 3, "No changes expected")
        call assert_equals(positions(:,2),  [1.0_real64,   4.0_real64,   5.0_real64], 3, "Change expected in first component")
        call assert_equals(positions(:,3),  [2.0_real64,   3.0_real64,   6.0_real64], 3, "Change expected in first and second components")
        call assert_equals(positions(:,4),  [7.0_real64,   4.0_real64,   5.0_real64], 3, "Change expected in second and third components")
        call assert_equals(positions(:,4),  [7.0_real64,   4.0_real64,   5.0_real64], 3, "Change expected in second and third components")
        call assert_equals(positions(:,5),  [6.0_real64,   7.0_real64,   8.0_real64], 3, "Change expected in all components")
        call assert_equals(positions(:,6),  [9.0_real64,   8.0_real64,   9.0_real64], 3, "Change expected in first and third components")
        call assert_equals(positions(:,7),  [8.0_real64,   1.5_real64,   9.0_real64], 3, "Change expected in first and second components")        
        call assert_equals(positions(:,8),  [2.5_real64,   3.5_real64,   7.0_real64], 3, "Change expected in all components")                
        call assert_equals(positions(:,9),  [6.0_real64,   5.0_real64,   4.0_real64], 3, "Change expected in all components")        
        call assert_equals(positions(:,10), [3.0_real64,   4.5_real64,   2.0_real64], 3, "Change expected in all components")        
        call assert_equals(positions(:,11), [1.0_real64,   1.5_real64,   2.5_real64], 3, "Change expected in first component")        
        
        deallocate(positions)
    end subroutine test_periodic_boundary_conditions

    subroutine test_distance_minimum_image()
        implicit none
        
    end subroutine test_distance_minimum_image
end module system_test