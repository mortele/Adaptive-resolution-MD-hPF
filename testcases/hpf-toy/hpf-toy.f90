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
    number_of_field_nodes_x = 5
    number_of_field_nodes_y = 5
    number_of_field_nodes_z = 5
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
    
    positions(:, 1 )  = [ 9.861513516680217_real64,  5.854609220329756_real64,  4.090152684461406_real64   ]
    positions(:, 2 )  = [ 2.9071826059146924_real64, 8.743543693640504_real64,  9.004829883778097_real64   ]
    positions(:, 3 )  = [ 2.461251549282577_real64,  6.473427427545419_real64,  0.40419124737701484_real64 ]
    positions(:, 4 )  = [ 2.873699847732466_real64,  4.562054499100865_real64,  0.8750826555311375_real64  ]
    positions(:, 5 )  = [ 8.61517700020657_real64,   7.5947341584529235_real64, 7.817588697757908_real64   ]
    positions(:, 6 )  = [ 8.329323455469527_real64,  6.48491076670158_real64,   1.0595951799017522_real64  ]
    positions(:, 7 )  = [ 7.539162628881729_real64,  6.8390933240322624_real64, 7.783573997831711_real64   ]
    positions(:, 8 )  = [ 5.502628575319229_real64,  5.086691204907025_real64,  9.508740311058814_real64   ]
    positions(:, 9 )  = [ 4.68875984367639_real64,   0.7111235962240725_real64, 7.499310540550107_real64   ]
    positions(:, 10 ) = [ 5.628789714994156_real64,  0.5249494827666557_real64, 6.096779978340795_real64   ]

    call compute_density_field(positions, masses)

    do i = 1, number_of_field_nodes_x
        do j = 1, number_of_field_nodes_y
            do k = 1, number_of_field_nodes_z
                print *, i,j,k, density_field(i,j,k)
            end do
        end do
    end do

    
end program