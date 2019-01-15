program hpftoy
    use, intrinsic :: iso_fortran_env,    only: real64, int32, output_unit
    use particles,          only: positions,                &
                                  velocities,               &
                                  forces,                   &
                                  masses,                   &
                                  types
    use file_writer,        only: write_state,              & ! subroutine
                                  write_info,               & ! subroutine
                                  read_state_lammps           ! subroutine
    use integrator,         only: integrate_one_step          ! subroutine
    use potential,          only: Ek, V_md, E,              &
                                  compute_forces_md           ! subroutine           
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
                                  number_of_particles,      &
                                  number_of_field_nodes_x,  &
                                  number_of_field_nodes_y,  &
                                  number_of_field_nodes_z
    use field,              only: allocate_field_arrays,    & ! Subroutine
                                  compute_density_field,    & ! Subroutine
                                  compute_density_gradient, & ! Subroutine
                                  position_of_density_nodes,&
                                  density_field,            &
                                  density_gradient,         &
                                  interpolate_density_gradient! Subroutine
    use initial_states,     only: allocate_arrays             ! Subroutine

    implicit none
    number_of_dimensions = 3
    number_of_particles  = 10
    
    
    
end program