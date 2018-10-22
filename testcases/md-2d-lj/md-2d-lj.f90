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


    call read_state_lammps(lammps_file, positions, velocities, forces, types, masses)
    
    call compute_forces(positions, forces)

    do i = 1, 10
        print '(F20.15,F20.15,F20.15)', positions(:,i)
    end do


end program