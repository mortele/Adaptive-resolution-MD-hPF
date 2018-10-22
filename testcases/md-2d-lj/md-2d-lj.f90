program md2dlj
    use, intrinsic :: iso_fortran_env,    only: real64, int32
    use particles,          only: positions,                &
                                  velocities,               &
                                  forces,                   &
                                  masses,                   &
                                  types
    use file_writer,        only: write_state,              & ! subroutine
                                  write_info,               & ! subroutine
                                  read_state_lammps           ! subroutine
    use integrator,         only: integrate_one_step          ! subroutine
    use potential,          only: Ek, V, E,                 &
                                  compute_forces              ! subroutine           
    use system,             only: system_size
    use sampler,            only: kinetic_energy,           &
                                  potential_energy,         &
                                  total_energy,             &
                                  store_energy                ! subroutine
    use parameters,         only: number_of_particles,      &
                                  number_of_time_steps,     &
                                  system_size_x,            &
                                  system_size_y,            &
                                  system_size_z,            &
                                  time_step,                &
                                  lennard_jones_sigma,      &
                                  lennard_jones_epsilon,    &
                                  lennard_jones_cutoff,     &
                                  temperature,              &
                                  number_of_dimensions,     &
                                  number_of_particles
    implicit none
    character (len=*), parameter :: lammps_file = "lammps.dump"
    integer (int32) :: i
    
    time_step               = 0.001_real64
    lennard_jones_sigma     = 1.0_real64
    lennard_jones_epsilon   = 1.0_real64
    lennard_jones_cutoff    = 4.0_real64

    call read_state_lammps(lammps_file, positions, velocities, forces, types, masses)
    
    call compute_forces(positions, forces)
    call write_array_to_file("forces.dump",     forces)
    call write_array_to_file("positions.dump",  positions)

    contains 
        subroutine open_file(file_name, file_ID)
            implicit none
            character (len=*), intent(in)     :: file_name
            integer (int32),   intent(in out) :: file_ID
            logical :: file_exists

            inquire(file = file_name, exist = file_exists)   
            if (file_exists) then
                ! Replace the existing file.
                open(   newunit = file_ID,          &
                        file    = file_name,        &
                        status  = "replace",        &
                        action  = "write")
            else 
    
                ! Create a new file.
                open(   newunit = file_ID,          &
                        file    = file_name,        &
                        status  = "new",            &
                        action  = "write")
            end if
        end subroutine


        subroutine write_array_to_file(file_name, array)
            implicit none
            character (len=*), intent(in) :: file_name
            real (real64), dimension(:,:), allocatable, intent(in) :: array
            integer (int32) :: i, file_ID
            
            call open_file(file_name, file_ID)

            do i = 1, number_of_particles
                write(file_ID, '(F20.15,F20.15,F20.15)') array(:,i)
            end do
        end subroutine

end program