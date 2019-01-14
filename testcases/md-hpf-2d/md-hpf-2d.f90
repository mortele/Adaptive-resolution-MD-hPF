program mdhpf2d
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
                                  number_of_particles,      &
                                  number_of_field_nodes_x,  &
                                  number_of_field_nodes_y,  &
                                  number_of_field_nodes_z
    use field,              only: allocate_field_arrays,    & ! Subroutine
                                  compute_density_field,    & ! Subroutine
                                  compute_density_gradient, & ! Subroutine
                                  position_of_density_nodes,&
                                  density_field,            &
                                  density_gradient

    implicit none
    real (real64), allocatable, dimension(:,:) :: field_visualization_positions
    character (len=*), parameter :: lammps_file = "lammps.dump"
    integer (int32) :: i, number_of_field_nodes_total
    
    time_step               = 0.001_real64
    lennard_jones_sigma     = 1.0_real64
    lennard_jones_epsilon   = 1.0_real64
    lennard_jones_cutoff    = 4.0_real64

    number_of_field_nodes_x = 10
    number_of_field_nodes_y = 5
    number_of_field_nodes_z = 3
    ! Start on the second LAMMPS step and compare positions, velocities, and
    ! energies.
    !call read_state_lammps(lammps_file, 3, positions, velocities, forces, types, masses)
    !forces = 0.0_real64
    !call compute_forces_md(positions, forces)
    !call write_array_to_file("positions.dump",   positions)
    !call write_array_to_file("velocities.dump",  velocities)
    !call write_array_to_file("forces.dump",      forces)
    
    ! Integrate a single step, and compare again.
    !call integrate_one_step(positions, velocities, forces)
    !call write_array_to_file("positions2.dump",  positions)
    !call write_array_to_file("velocities2.dump", velocities)
    !call write_array_to_file("forces2.dump",     forces)
    
    ! Reload on initial LAMMPS step and integrate.
    call read_state_lammps(lammps_file, 1, positions, velocities, forces, types, masses)
    call compute_forces_md(positions, forces)
    call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
    call compute_density_field(positions, masses)

    number_of_field_nodes_total =   number_of_field_nodes_x *       &
                                    number_of_field_nodes_y *       &
                                    number_of_field_nodes_z
    allocate(field_visualization_positions(number_of_dimensions, number_of_particles+number_of_field_nodes_total))

    do i = 1, 1000
        if (modulo(i,10) == 0) then
            write(output_unit, fmt="(f10.1,a1)", advance="no") (i / 1000.0)*100, "%"
            write(output_unit, fmt="(a2)", advance="no") char(13) ! carriage return, \r
        end if
        !call write_state(positions, i, types)
        !call write_state_density_visualization(positions, field_visualization_positions, i, types)
        !call write_energy_to_file(Ek, V)
        call integrate_one_step(positions, velocities, forces)
        !call compute_density_field(positions, masses)
    end do
    print *, "thermalization done"

    do i = 1, 10000
        if (modulo(i,10) == 0) then
            write(output_unit, fmt="(f10.1,a1)", advance="no") (i / 1000.0)*100, "%"
            write(output_unit, fmt="(a2)", advance="no") char(13) ! carriage return, \r
        end if

        !call write_state(positions, i, types)
        call write_state_density_visualization(positions, field_visualization_positions, i, types)
        !call write_energy_to_file(Ek, V)
        call integrate_one_step(positions, velocities, forces)
        call compute_density_field(positions, masses)
    end do

    print *, "md-2d-lj.f90 exited with exit code 0"

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

        subroutine write_energy_to_file(kinetic_energy, potential_energy)
            implicit none
            real (real64), intent(in) :: kinetic_energy, potential_energy
            character (len=*), parameter :: file_name = "energy.dump"
            integer (int32), save :: file_ID
            logical, save :: file_open = .false.

            if (.not. file_open) then
                call open_file(file_name, file_ID)
                file_open = .true.
            end if

            write(file_ID, '(F20.15,F20.15)')                                &
                            kinetic_energy/number_of_particles,                     &
                            potential_energy/number_of_particles
        end subroutine

        subroutine write_array_to_file(file_name, array)
            implicit none
            character (len=*), intent(in) :: file_name
            real (real64), dimension(:,:), allocatable, intent(in) :: array
            integer (int32) :: i, file_ID
            
            call open_file(file_name, file_ID)

            do i = 1, number_of_particles
                write(file_ID, '(F25.15,F25.15,F25.15)') array(:,i)
            end do
        end subroutine

        subroutine compute_density_field_visualization( density_field,                  &
                                                        position_of_density_nodes,      &
                                                        field_visualization_positions)
            implicit none
            real (real64), dimension(:,:,:),   intent(in)     :: density_field
            real (real64), dimension(:,:,:,:), intent(in)     :: position_of_density_nodes
            real (real64), dimension(:,:), allocatable, intent(in out) :: field_visualization_positions

            integer (int32) :: i, j, k, n_tot, index 
            real (real64)   :: sum_z

            if (.not. allocated(field_visualization_positions)) then
                n_tot = number_of_field_nodes_x * number_of_field_nodes_y
                allocate(field_visualization_positions(number_of_dimensions, n_tot))
            end if

            field_visualization_positions = 0.0_real64
            index = 1
            do i = 1, number_of_field_nodes_x
                do j = 1, number_of_field_nodes_y
                    field_visualization_positions(:,index) = position_of_density_nodes(:,i,j,1)
                    index = index + 1
                end do
            end do

            index = 1
            do i = 1, number_of_field_nodes_x
                do j = 1, number_of_field_nodes_y
                    sum_z = 0.0_real64
                    do k = 1, number_of_field_nodes_z
                        sum_z = sum_z + density_field(i,j,k)
                    end do
                    field_visualization_positions(3, index) = ((field_visualization_positions(3, index) + sum_z)**(1./3.))**3 - 15
                    index = index + 1
                    if (field_visualization_positions(1, index) > system_size_x / 2.0) then
                        field_visualization_positions(3, index) = field_visualization_positions(3, index) + 50
                    end if
                end do
            end do
        end subroutine

        subroutine write_state_density_visualization(   positions,                      &
                                                        field_visualization_positions,  &
                                                        step,                              &
                                                        types)
            implicit none
            real (real64),   dimension(:,:), intent(in) :: positions
            real (real64),   dimension(:,:), allocatable, intent(in out) :: field_visualization_positions
            integer (int32), dimension(:),   intent(in) :: types
            integer (int32),               intent(in) :: step

            real (real64), dimension(:,:), allocatable, save :: new_positions
            integer (int32), dimension(:),   allocatable, save :: new_types
            integer (int32) :: i, j, k, n_part, n_tot, index

            call compute_density_field_visualization(density_field, position_of_density_nodes, field_visualization_positions)

            n_tot = number_of_field_nodes_x * number_of_field_nodes_y
            if (.not. allocated(new_positions)) then
                allocate(new_positions(number_of_dimensions, number_of_particles + n_tot))
                allocate(new_types    (number_of_particles + n_tot))
            end if

            do i = 1, number_of_particles
                new_positions(:,i) = positions(:,i)
                new_types(i)       = types(i)
            end do

            index = 1
            do j = number_of_particles+1, number_of_particles+n_tot
                new_positions(:,j) = field_visualization_positions(:,index)
                new_types(j) = 99
                index = index + 1
            end do

            n_part = number_of_particles
            number_of_particles = n_part + n_tot
            call write_state(new_positions, step, new_types)
            number_of_particles = n_part

        end subroutine

end program