program hpftoy
    use, intrinsic :: iso_fortran_env,    only: real64, int32, output_unit
    use particles,          only: positions,                &
                                  velocities,               &
                                  forces,                   &
                                  forces_hpf,               &
                                  masses,                   &
                                  types
    use file_writer,        only: write_state,              & ! subroutine
                                  write_info,               & ! subroutine
                                  read_state_lammps           ! subroutine
    use integrator,         only: integrate_one_step          ! subroutine
    use potential,          only: compute_forces_hpf!,       & ! subroutine 
    !                              Ek, V_md, E              
    use system,             only: system_size
    !use sampler,            only: kinetic_energy,           &
    !                              potential_energy,         &
    !                              total_energy,             &
    !                              store_energy                ! subroutine
    use parameters,         only: number_of_particles,      &
                                  number_of_time_steps,     &
                                  system_size_x,            &
                                  system_size_y,            &
                                  system_size_z,            &
                                  time_step,                &
                                  number_of_dimensions,     &
                                  number_of_particles,      &
                                  number_of_field_nodes_x,  &
                                  number_of_field_nodes_y,  &
                                  number_of_field_nodes_z,  &
                                  chi,                      &
                                  kappa
    use field,              only: allocate_field_arrays,    & ! Subroutine
                                  compute_density_field,    & ! Subroutine
                                  compute_density_gradient, & ! Subroutine
                                  position_of_density_nodes,&
                                  density_field,            &
                                  density_gradient,         &
                                  interpolate_density_gradient! Subroutine
    use initial_states,     only: allocate_arrays             ! Subroutine

    implicit none
    integer (int32) :: i, j, k


    ! Setup parameters to mimic the OCCAM run
    number_of_dimensions    = 3
    number_of_particles     = 10
    system_size             = [10.0, 10.0, 10.0]
    system_size_x           = 10.0
    system_size_y           = 10.0
    system_size_z           = 10.0
    number_of_field_nodes_x = 2
    number_of_field_nodes_y = 2
    number_of_field_nodes_z = 2
    time_step               = 0.03_real64
    number_of_time_steps    = 100
    kappa                   = 0.05
    chi                     = 0.0

    call allocate_arrays(positions,                         &
                         velocities,                        &
                         forces,                            &
                         masses,                            &
                         types)
    call allocate_field_arrays(density_field,               &
                               density_gradient,            &
                               position_of_density_nodes)
    velocities = 0.0_real64
    forces     = 0.0_real64
    forces_hpf = 0.0_real64
    masses     = 1.0_real64
    
    positions(:, 1 ) = [0.05815856656600038_real64,0.025342305369190687_real64,0.19873935016981903_real64 ]
    positions(:, 2 ) = [0.07215982897620887_real64,0.2504571863305104_real64,0.38362945814487526_real64 ]
    positions(:, 3 ) = [0.49720339870241104_real64,0.8872456338742205_real64,0.9632329182983882_real64 ]
    positions(:, 4 ) = [0.5905505257877675_real64,0.33492055184731184_real64,0.8214094476851357_real64 ]
    positions(:, 5 ) = [0.5156196908805613_real64,0.7799713140086558_real64,0.007959114457485872_real64 ]
    positions(:, 6 ) = [0.8737540746260515_real64,0.12371839911797189_real64,0.6094189843776412_real64 ]
    positions(:, 7 ) = [0.7794755801413096_real64,0.44018012706879395_real64,0.10548340557947777_real64 ]
    positions(:, 8 ) = [0.9931810993676221_real64,0.26069746565869656_real64,0.7264137166754773_real64 ]
    positions(:, 9 ) = [0.2898648348703551_real64,0.5064354936235451_real64,0.8455429407939724_real64 ]
    positions(:, 10 ) = [0.38696815531713213_real64,0.898330488912219_real64,0.3928309115958717_real64 ]
    call compute_density_field(positions, masses)

    do i = 1, number_of_field_nodes_x
        do j = 1, number_of_field_nodes_y
            do k = 1, number_of_field_nodes_z
                print *, i,j,k, density_field(i,j,k)
            end do
        end do
    end do

    
end program