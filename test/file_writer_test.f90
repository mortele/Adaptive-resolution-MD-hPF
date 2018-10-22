module file_writer_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,              & ! Subroutine
                                    assert_true                   ! Subroutine
    use file_writer,        only:   write_state,                & ! Subroutine
                                    read_state,                 & ! Subroutine
                                    close_outfile,              & ! Subroutine
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
                test_read_state,            &
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

    subroutine test_read_state()
        implicit none
        real (real64),   dimension(:,:), allocatable :: positions_read
        real (real64),   dimension(:),   allocatable :: types_real
        integer (int32), dimension(:),   allocatable :: types_read
        integer (int32) :: i
        
        number_of_dimensions = 1
        number_of_particles  = 3
        out_file_name        = "testingfw.xyz"


        allocate(positions (number_of_dimensions, number_of_particles))
        allocate(types     (number_of_particles))

        positions(1,1) = 25985.25_real64
        positions(1,2) = 95224.10_real64
        positions(1,3) = 0.259827_real64
        types(1) = 1
        types(2) = 2
        types(3) = 30

        call close_outfile()
        call write_state(positions, 1, types)
        call close_outfile()
        call read_state(out_file_name, positions, types)
        call close_outfile()

        call assert_equals([25985.25_real64, 95224.10_real64, 0.259827_real64], positions(1,:), 3, 1e-15_real64, "1  test_read_state : read positions are not equal to the positions just written to file")
        call assert_equals([1,2,30], types, 3, "2  test_read_state : read types are not equal to the types just written to file")
        deallocate(positions)
        deallocate(types)

        number_of_particles  = 10
        number_of_dimensions = 3
        allocate(positions     (number_of_dimensions, number_of_particles))
        allocate(positions_read(number_of_dimensions, number_of_particles))
        allocate(types         (number_of_particles))
        allocate(types_real    (number_of_particles))
        allocate(types_read    (number_of_particles))
        call random_number(positions)
        call random_number(types_real)
        types       = floor(10 * types_real)
        positions   = 10_real64 * positions
        
        call close_outfile()
        call write_state(positions, 1, types)
        call close_outfile()
        call read_state(out_file_name, positions_read, types_read)
        call close_outfile()

        do i = 1, number_of_particles
            call assert_equals(positions(:,i), positions_read(:,i), 3, "3  test_read_state : read positions are not equal to the positions just written to file")
            call assert_equals(types(i),       types_read(i),          "4  test_read_state : read types are not equal to the types just written to file")
        end do


        deallocate(positions)
        deallocate(positions_read)
        deallocate(types)
        deallocate(types_real)
        deallocate(types_read)


    end subroutine

    subroutine test_write_info()
    end subroutine


end module file_writer_test
