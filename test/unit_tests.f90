program unit_tests
    use fruit
    use system_test
    implicit none

    call init_fruit
    
    ! Battery of tests
    call test_periodic_boundary_conditions
    call test_distance_minimum_image
    
    call fruit_summary

end program unit_tests