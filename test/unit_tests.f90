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
    call field_all_tests()
    call sampler_all_tests()
    call integrator_all_tests()
    call file_writer_all_tests()
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
        print *, " ..running test: test_compute_forces_md"
        call set_unit_name ('test_compute_forces_md')    
        call run_test_case(test_compute_forces_md, "test_compute_forces_md")
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
    end subroutine initial_states_all_tests

    subroutine field_all_tests()
        use field_test

        call setup
        print *, " "
        print *, " ..running test: test_allocate_field_arrays"
        call set_unit_name ('test_allocate_field_arrays')    
        call run_test_case(test_allocate_field_arrays, "test_allocate_field_arrays")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_compute_density_field"
        call set_unit_name ('test_compute_density_field')    
        call run_test_case(test_compute_density_field, "test_compute_density_field")
        call teardown        

        call setup
        print *, " "
        print *, " ..running test: test_compute_density_field_periodic_boundaries"
        call set_unit_name ('test_compute_density_field_periodic_boundaries')    
        call run_test_case(test_compute_density_field_periodic_boundaries, "test_compute_density_field_periodic_boundaries")
        call teardown        
        
        call setup
        print *, " "
        print *, " ..running test: test_compute_density_gradient"
        call set_unit_name ('test_compute_density_gradient')    
        call run_test_case(test_compute_density_gradient, "test_compute_density_gradient")
        call teardown  
        
        call setup
        print *, " "
        print *, " ..running test: test_interpolate_density_field"
        call set_unit_name ('test_interpolate_density_field')    
        call run_test_case(test_interpolate_density_field, "test_interpolate_density_field")
        call teardown  

        call setup
        print *, " "
        print *, " ..running test: test_interpolate_density_gradient"
        call set_unit_name ('test_interpolate_density_gradient')    
        call run_test_case(test_interpolate_density_gradient, "test_interptest_interpolate_density_gradientolate_density_field")
        call teardown  

        
    end subroutine field_all_tests

    subroutine sampler_all_tests()
        use sampler_test

        call setup
        print *, " "
        print *, " ..running test: test_store_energy"
        call set_unit_name ('test_store_energy')    
        call run_test_case(test_store_energy, "test_store_energy")
        call teardown     

    end subroutine sampler_all_tests

    subroutine integrator_all_tests()
        use integrator_test

        call setup
        print *, " "
        print *, " ..running test: test_half_move"
        call set_unit_name ('test_half_move')    
        call run_test_case(test_half_move, "test_half_move")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_move"
        call set_unit_name ('test_move')    
        call run_test_case(test_move, "test_move")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_integrate_one_step"
        call set_unit_name ('test_integrate_one_step')    
        call run_test_case(test_integrate_one_step, "test_integrate_one_step")
        call teardown     
        
    end subroutine integrator_all_tests

    subroutine file_writer_all_tests()
        use file_writer_test

        call setup
        print *, " "
        print *, " ..running test: test_write_state"
        call set_unit_name ('test_write_state')    
        call run_test_case(test_write_state, "test_write_state")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_read_state"
        call set_unit_name ('test_read_state')    
        call run_test_case(test_read_state, "test_read_state")
        call teardown   

        call setup
        print *, " "
        print *, " ..running test: test_write_info"
        call set_unit_name ('test_write_info')    
        call run_test_case(test_write_info, "test_write_info")
        call teardown     

        call setup
        print *, " "
        print *, " ..running test: test_read_state_lammps"
        call set_unit_name ('test_read_state_lammps')    
        call run_test_case(test_read_state_lammps, "test_read_state_lammps")
        call teardown     
        
    end subroutine file_writer_all_tests

end program unit_tests