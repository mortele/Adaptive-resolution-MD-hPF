program unit_tests
    use fruit
    use parameters,     only:   silent_output_file,     &
                                silent_output_ID
    implicit none

    call init_fruit()
    call open_alternative_output_file()
    
    ! Battery of tests
    call system_all_tests()
    call random_generator_all_tests()
    call potential_all_tests()
    call initial_states_all_tests()
    ! ================
    
    call fruit_summary()
    call close_alternative_output_file()
contains

    subroutine open_alternative_output_file()
        logical :: file_exists

        ! Check if the file exists.
        inquire(file = silent_output_file, exist = file_exists)

        if (file_exists) then
            ! Replace the existing file.
            open(   newunit = silent_output_ID,     &
                    file    = silent_output_file,   &
                    status  = "replace",            &
                    action  = "write")
        else 

            ! Create a new file.
            open(   newunit = silent_output_ID,     &
                    file    = silent_output_file,   &
                    status  = "new",                &
                    action  = "write")
        end if
    end subroutine open_alternative_output_file

    subroutine close_alternative_output_file()
        logical :: is_open

        inquire(file = silent_output_file, opened = is_open)

        if (is_open) then
            close(silent_output_ID, status = 'delete')
        end if
    end subroutine close_alternative_output_file
    

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

        call setup
        print *, " "
        print *, " ..running test: test_compute_forces"
        call set_unit_name ('test_compute_forces')    
        call run_test_case(test_compute_forces, "test_compute_forces")
        call teardown     

    end subroutine potential_all_tests

    subroutine initial_states_all_tests()
        use initial_states_test

        call setup
        print *, " "
        print *, " ..running test: test_allocate_arrays"
        call set_unit_name ('test_allocate_arrays')    
        call run_test_case(test_allocate_arrays, "test_allocate_arrays")
        call teardown      
        
        call setup
        print *, " "
        print *, " ..running test: test_random_initial_state"
        call set_unit_name ('test_random_initial_state')    
        call run_test_case(test_random_initial_state, "test_random_initial_state")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_fcc_initial_state"
        call set_unit_name ('test_fcc_initial_state')    
        call run_test_case(test_fcc_initial_state, "test_fcc_initial_state")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_setup_initial_state"
        call set_unit_name ('test_setup_initial_state')    
        call run_test_case(test_setup_initial_state, "test_setup_initial_state")
        call teardown

    end subroutine initial_states_all_tests

end program unit_tests