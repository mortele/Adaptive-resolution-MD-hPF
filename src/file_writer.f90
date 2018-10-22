module file_writer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,     only:   system_size 
    use parameters, only:   number_of_particles,    &
                            number_of_dimensions,   &
                            out_file_name,          &
                            number_of_particles,    &
                            number_of_dimensions,   &
                            temperature,            &
                            time_step,              &
                            step_length,            &
                            number_of_time_steps,   &
                            info_file_name

    
    implicit none
    private

    public ::   write_state,            &
                read_state,             &
                close_outfile,          &
                write_info             

contains

    subroutine write_state(positions, current_time_step, types)
        implicit none
        real    (real64), dimension(:,:), intent(in)           :: positions
        integer (int32),  dimension(:),   optional, intent(in) :: types
        integer (int32),                  optional, intent(in) :: current_time_step
        integer (int32)       :: i, j
        integer (int32), save :: file_ID 
        logical :: file_exists
        logical :: file_is_open
        
        ! If the file is not open, we have to open it before we continue.
        inquire(file = out_file_name, exist = file_exists)
        if (file_exists) then
            inquire(file = out_file_name, opened = file_is_open)
        else 
            file_is_open = .false.
        end if

        if (.not. file_is_open)  then

            ! Check if the file exists.
            inquire(file = out_file_name, exist = file_exists)

            if (file_exists) then

                ! Replace the existing file.
                open(   newunit = file_ID,          &
                        file    = out_file_name,    &
                        status  = "replace",        &
                        action  = "write")
                file_is_open = .true.
            else 

                ! Create a new file.
                open(   newunit = file_ID,          &
                        file    = out_file_name,    &
                        status  = "new",            &
                        action  = "write")
                file_is_open = .true.
            end if
        end if
    
        ! Write positions of all particles to out_file_name in the .xyz 
        ! file format. The .xyz format starts with a line containing the 
        ! total number of particles, followed by an optional line. We use 
        ! this second line to write the time step.
        write(file_ID, *) number_of_particles

        if (present(current_time_step)) then
            write(file_ID, *) "#current time step ", current_time_step
        else 
            write(file_ID, *) "#"
        end if
        
        do i = 1, number_of_particles  

            ! The start of each row is the particle identifier, contained in
            ! the type(:) array.
            if (present(types)) then    
                write(file_ID, fmt = "( i3 )", advance = "no") types(i)
            else
                write(file_ID, fmt = "( i3 )", advance = "no") 1
            end if


            ! Write the x, y, and z positions, supressing the newline with
            ! advance = "no".
            do j = 1, number_of_dimensions
                write(  file_ID,                 &
                        fmt     = "( f30.18 )",  &
                        advance = "no")             positions(j, i)
            end do

            ! Write a newline after the position of each particle.
            write(file_ID, *) " " 
        end do
    end subroutine write_state

    subroutine read_state(file_name, positions, types)
        implicit none
        character  (len=*),                  intent(in)     :: file_name
        real       (real64), dimension(:,:), intent(in out) :: positions
        integer    (int32),  dimension(:),   intent(in out) :: types
        
        integer (int32) ::  number_of_particles_in_file,        &
                            time_step_in_file,                  &
                            file_ID,                            &
                            i
        real (real64), dimension(:), allocatable :: type_position_line
        character (len=200) :: time_step_line
        logical :: file_exists

        inquire(file = file_name, exist = file_exists)
        if (file_exists) then
            open(   newunit = file_ID,      &
                    file    = file_name,    &
                    status  = "old",        &
                    action  = "read")
        else 
            print *, "Attempted to open file ", file_name
            error stop "Error: The file does not exist."
        end if 

        read(file_ID, *) number_of_particles_in_file
        read(file_ID, *) time_step_line

        allocate(type_position_line (number_of_dimensions+1))

        do i = 1, number_of_particles_in_file
            read(file_ID, *) type_position_line
            types(i) = int(type_position_line(1))
            positions(:,i) = type_position_line(2:number_of_dimensions+1)
        end do

        number_of_particles = number_of_particles_in_file
    end subroutine read_state

    subroutine write_info
        implicit none
        integer (int32), save :: file_ID 
        logical               :: file_exists

        ! Check if the file exists.
        inquire(file = info_file_name, exist = file_exists)

        if (file_exists) then

            ! Replace the existing file.
            open(   newunit = file_ID,          &
                    file    = info_file_name,   &
                    status  = "replace",        &
                    action  = "write")
        else 

            ! Create a new file.
            open(   newunit = file_ID,          &
                    file    = info_file_name,   &
                    status  = "new",            &
                    action  = "write")
        end if

        write(file_ID, *) "number_of_time_steps ", number_of_time_steps
        write(file_ID, *) "number_of_particles ", number_of_particles
        write(file_ID, *) "number_of_dimensions ", number_of_dimensions
        write(file_ID, *) "temperature ", temperature
        write(file_ID, *) "time_step ", time_step
        write(file_ID, *) "step_length ", step_length
        close(file_ID)
    end subroutine write_info


    subroutine close_outfile() 
        logical :: is_open
        integer :: file_ID

        inquire(file = out_file_name, opened = is_open, number = file_ID)

        if (is_open) then
            close(file_ID, status = "keep")
        end if 
    end subroutine
end module file_writer