module file_writer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use parameters, only:   number_of_particles,  &
                            number_of_dimensions, &
                            out_file_name
    implicit none
    private

    public :: write_state

contains

    subroutine write_state(positions, time_step)
        implicit none
        real    (real64), dimension(:,:), intent(in) :: positions
        integer (int32),  optional,       intent(in) :: time_step
        integer (int32)       :: i, j
        integer (int32), save :: file_ID 
        logical,         save :: file_is_open = .false.
        

        ! If the file is not open, we have to open it before we continue.
        if (.not. file_is_open)  then

            

        end if
        
    end subroutine write_state

end module file_writer