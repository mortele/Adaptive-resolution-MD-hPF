program unit_tests
    use fruit
    implicit none

    call init_fruit
    
    ! Battery of tests
    call system_all_tests()
    call random_generator_all_tests()
    call potential_all_tests()
    ! ================
    
    call fruit_summary
    
contains
    
    subroutine system_all_tests()
        use system_test
        
        call setup
        print *, " "
        print *, " ..running test: test_periodic_boundary_conditions"
        call set_unit_name ('test_periodic_boundary_conditions')    
        call run_test_case(test_periodic_boundary_conditions, "test_periodic_boundary_conditions")
        call teardown

        call setup
        print *, " "
        print *, " ..running test: test_distance_minimum_image"
        call set_unit_name ('test_distance_minimum_image')    
        call run_test_case(test_distance_minimum_image, "test_distance_minimum_image") 
        call teardown

        call setup
        print *, " "
        print *, " ..running test: test_remove_linear_momentum"
        call set_unit_name ('test_remove_linear_momentum')    
        call run_test_case(test_remove_linear_momentum, "test_remove_linear_momentum")
        call teardown
    end subroutine system_all_tests

    subroutine random_generator_all_tests()
        use random_generator_test

        call setup
        print *, " "
        print *, " ..running test: test_random_normal"
        call set_unit_name ('test_random_normal')    
        call run_test_case(test_random_normal, "test_random_normal")
        call teardown        
    end subroutine random_generator_all_tests

    subroutine potential_all_tests()
        use potential_test

        call setup
        print *, " "
        print *, " ..running test: test_lennard_jones_force"
        call set_unit_name ('test_lennard_jones_force')    
        call run_test_case(test_lennard_jones_force, "test_lennard_jones_force")
        call teardown      
        
        call setup
        print *, " "
        print *, " ..running test: test_lennard_jones_potential"
        call set_unit_name ('test_lennard_jones_potential')    
        call run_test_case(test_lennard_jones_potential, "test_lennard_jones_potential")
        call teardown     
    end subroutine potential_all_tests
end program unit_tests