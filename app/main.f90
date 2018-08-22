program main
    use, intrinsic :: iso_fortran_env,    only: real64, int32
    use particles,          only: positions, velocities, forces, masses, types
    use system,             only: system_size
    use initial_states,     only: random_initial_state
    use file_writer,        only: write_state
    use integrator,         only: integrate_one_step
    use parameters,         only: number_of_particles,  &
                                  number_of_dimensions, &
                                  number_of_time_steps
    implicit none
    
    integer (int32) :: time_step = 0

    ! Allocate particle arrays and generate the initial state of the system.
    allocate(positions (number_of_dimensions, number_of_particles))
    allocate(velocities(number_of_dimensions, number_of_particles))
    allocate(forces    (number_of_dimensions, number_of_particles))
    allocate(masses    (number_of_particles))
    allocate(types     (number_of_particles))
    call random_initial_state(positions, velocities, types)

    ! Write the initial state to a file for analysis.
    call write_state(positions, time_step, types)

    ! Integrate the equations of motion using the velocity Verlet algorithm.
    do time_step = 0, number_of_time_steps
        call integrate_one_step(positions, velocities, forces)
    end do

    ! Deallocate all arrays. This is probably not neccessary, but ...
    deallocate(positions)
    deallocate(velocities)
    deallocate(forces)
    deallocate(masses)
    deallocate(types)

    ! Since the automatic Visual Studio Code terminal doesnt show exit codes, we
    ! add a simple output showing the program exited with code 0.
    print *, "AdapResoMD-hPF exiting with code 0."
end program main