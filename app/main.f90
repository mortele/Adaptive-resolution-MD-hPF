program main
    use, intrinsic :: iso_fortran_env,    only: real64, int32
    use particles,          only: positions, velocities, forces, masses, types
    use initial_states,     only: setup_initial_state       ! subroutine
    use file_writer,        only: write_state,          &   ! subroutine
                                  write_info                ! subroutine
    use integrator,         only: integrate_one_step        ! subroutine
    use potential,          only: Ek, V, E
    use sampler,            only: kinetic_energy,       &
                                  potential_energy,     &
                                  total_energy,         &
                                  store_energy              ! subroutine
    use parameters,         only: number_of_particles,  &
                                  number_of_time_steps
    implicit none
    
    integer (int32) :: time_step = 0

    ! Setup the initial state of the particles. The initial_states module also
    ! handles allocating of the arrays to the appropriate sizes.
    call setup_initial_state(positions, velocities, forces, masses, types)

    ! Write the initial state to a file for analysis.
    call write_state(positions, time_step, types)

    ! Allocate statistics arrays and compute the initial energy.
    allocate(potential_energy(number_of_time_steps + 1))
    allocate(kinetic_energy  (number_of_time_steps + 1))
    allocate(total_energy    (number_of_time_steps + 1))
    call store_energy(kinetic_energy, potential_energy, total_energy)

    ! Integrate the equations of motion using the velocity Verlet algorithm.
    do time_step = 1, number_of_time_steps - 1
        call integrate_one_step(positions, velocities, forces)

        ! Dump the energies computed in the integration step (in the potential
        ! module) into the appropriate arrays.
        call store_energy(kinetic_energy, potential_energy, total_energy)

        ! Write the state to file.
        call write_state(positions, time_step, types)

        print *, E  / number_of_particles,  &
                 Ek / number_of_particles,  &
                 V  / number_of_particles
    end do

    call write_info

    ! Since the automatic Visual Studio Code terminal doesnt show exit codes, we
    ! add a simple output showing the program exited with code 0.
    print *, "AdapResoMD-hPF exiting with code 0."
end program main