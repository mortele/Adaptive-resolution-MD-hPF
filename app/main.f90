program main
    use, intrinsic :: iso_fortran_env,    only: real64, int32
    use particles,          only: positions, velocities, forces, masses, types
    use system,             only: system_size
    use initial_states,     only: random_initial_state      ! subroutine
    use file_writer,        only: write_state               ! subroutine
    use integrator,         only: integrate_one_step        ! subroutine
    use sampler,            only: kinetic_energy,       &
                                  potential_energy,     &
                                  total_energy,         &
                                  store_energy              ! subroutine
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

    ! Allocate statistics arrays and compute the initial energy.
    allocate(potential_energy(number_of_time_steps))
    allocate(kinetic_energy  (number_of_time_steps))
    allocate(total_energy    (number_of_time_steps))
    call store_energy(kinetic_energy, potential_energy, total_energy)

    ! Integrate the equations of motion using the velocity Verlet algorithm.
    do time_step = 0, number_of_time_steps
        call integrate_one_step(positions, velocities, forces)

        ! Dump the energies computed in the integration step (in the potential
        ! module) into the appropriate arrays.
        call store_energy(kinetic_energy, potential_energy, total_energy)
    end do

    ! Since the automatic Visual Studio Code terminal doesnt show exit codes, we
    ! add a simple output showing the program exited with code 0.
    print *, "AdapResoMD-hPF exiting with code 0."
end program main