module file_writer_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,              & ! Subroutine
                                    assert_true                   ! Subroutine
    use file_writer,        only:   write_state,                & ! Subroutine
                                    write_info                    ! Subroutine
    use system,             only:   system_size
    use particles,          only:   positions,                  &
                                    types                               
    use parameters,         only:   number_of_particles,        &
                                    number_of_dimensions,       &
                                    out_file_name,              &
                                    number_of_particles,        &
                                    number_of_dimensions,       &
                                    temperature,                &
                                    time_step,                  &
                                    step_length,                &
                                    number_of_time_steps,       &
                                    info_file_name,             &
                                    system_size_x,              &
                                    system_size_y,              &
                                    system_size_z
    implicit none
    private 

    public  ::  setup,                      &
                teardown,                   &
                test_write_state,           &
                test_write_info

contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_write_state()
        logical :: file_exists
        integer :: out_file_ID, info_file_ID


        number_of_dimensions = 3
        number_of_particles  = 10
        system_size_x        = 10.0
        system_size_y        = 10.0
        system_size_z        = 10.0
        system_size          = [system_size_x, system_size_y, system_size_z]
        out_file_name        = "testingfw.xyz"
        info_file_name       = "testinfofw.out"
        allocate(positions  (number_of_dimensions, number_of_particles))
        allocate(types      (number_of_particles))
        types = 1
        call random_number(positions)
        positions = positions * system_size(1)  ! Uniform random in range [0, 10)

        ! If the file exists, we delete it.
        inquire(file = out_file_name, exist = file_exists)
        if (file_exists) then
            open(   newunit = out_file_ID,      &
                    file    = out_file_name,    &
                    status  = "replace",        &
                    action  = "write")
            close(out_file_ID, status = 'delete')
        end if

        call write_state(positions, 1, types)
        inquire(file = out_file_name, exist = file_exists)
        call assert_true(file_exists, "1  test_write_state : write_state did not create a file with the correct name for the output")

        deallocate(positions)
        deallocate(types)
    end subroutine

    subroutine test_write_info()
    end subroutine

end module file_writer_test