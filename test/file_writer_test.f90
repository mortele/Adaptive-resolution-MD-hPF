module file_writer_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,              & ! Subroutine
                                    assert_true                   ! Subroutine
    use file_writer,        only:   write_state,                & ! Subroutine
                                    write_info                    ! Subroutine
    use system,             only:   system_size                                 
    use parameters,         only:   number_of_particles,        &
                                    number_of_dimensions,       &
                                    out_file_name,              &
                                    number_of_particles,        &
                                    number_of_dimensions,       &
                                    temperature,                &
                                    time_step,                  &
                                    step_length,                &
                                    number_of_time_steps,       &
                                    info_file_name
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
    end subroutine

    subroutine test_write_info()
    end subroutine

end module file_writer_test