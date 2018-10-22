program md2dlj
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

    print *, "hello"
    allocate(positions (3, 2))
    call random_number(positions)
    print *, positions


end program