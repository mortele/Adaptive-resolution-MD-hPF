module file_writer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use parameters, only:   number_of_particles,  &
                            number_of_dimensions, &
                            out_file_name
    
    implicit none
    private

    public :: write_state

contains

    subroutine write_state(positions, time_step, types)
        implicit none
        real    (real64), dimension(:,:), intent(in)           :: positions
        integer (int32),  dimension(:),   optional, intent(in) :: types
        integer (int32),                  optional, intent(in) :: time_step
        integer (int32)       :: i, j
        integer (int32), save :: file_ID 
        logical,         save :: file_is_open = .false.
        logical               :: file_exists
        

        ! If the file is not open, we have to open it before we continue.
        if (.not. file_is_open)  then

            ! Check if the file exists.
            inquire(file = out_file_name, exist = file_exists)

            if (file_exists) then

                ! Replace the existing file.
                open(   newunit = file_ID,          &
                        file    = out_file_name,    &
                        status  = "replace",        &
                        action  = "write")
            else 

                ! Create a new file.
                open(   newunit = file_ID,          &
                         file   = out_file_name,    &
                        status  = "new",            &
                        action  = "write")
            end if

            ! Write positions of all particles to out_file_name in the .xyz 
            ! file format. The .xyz format starts with a line containing the 
            ! total number of particles, followed by an optional line. We use 
            ! this second line to write the time step.
            write(file_ID, *) number_of_particles

            if (present(time_step)) then
                write(file_ID, *) "#time step ", time_step
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
                    write(  file_ID,                &
                            fmt     = "( f15.5 )",  &
                            advance = "no")             positions(j, i)
                end do

                ! Write a newline after the position of each particle.
                write(file_ID, *) " " 
            end do
        end if
        
    end subroutine write_state

end module file_writer