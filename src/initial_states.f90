module initial_states
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use random_generator, only: random_normal                   ! Function
    use system,           only: system_size,                &
                                remove_linear_momentum          ! Subroutine
    use parameters,       only: number_of_particles,        &
                                number_of_dimensions,       &
                                temperature,                &
                                fcc_lattice_constant,       &
                                fcc_number_of_unit_cells,   &
                                initial_configuration,      &
                                system_size_x,              &
                                system_size_y,              &
                                system_size_z              
    implicit none
    private 

    private ::  random_initial_state,    &
                fcc_initial_state,       &
                allocate_arrays,         &
                print_fcc_warning
    public  :: setup_initial_state

contains

    subroutine random_initial_state(positions, velocities, masses, types)
        implicit none
        real    (real64), dimension(:,:), intent(in out) :: positions
        real    (real64), dimension(:,:), intent(in out) :: velocities
        real    (real64), dimension(:),   intent(in out) :: masses
        integer (int32),  dimension(:),   intent(in out) :: types

        integer (int32)  :: i, j
        real    (real64) :: T
        real    (real64), dimension(:) :: S(number_of_dimensions)
        T = temperature

        ! Assign each particle a mass of 1.0 and type 1.
        masses = 1.0
        types  = 1

        do i = 1, number_of_particles
            ! Assign a uniform random position to each particle, then scale each
            ! dimension with the system size to get a uniform distribution in 
            ! the complete simulation box.
            call random_number(positions(:, i))
            
            do j = 1, number_of_dimensions
                positions(j, i) = positions(j, i) * system_size(j)
            end do
        end do

        do i = 1, number_of_particles
            ! Assign the velocity of each particle to a normal distribution, 
            ! with zero mean and unit variance.
            call random_normal(velocities(:, i))
            
            ! Scale the standard deviation to be (kB T / m) according to the 
            ! Boltzmann distribution at temperature T.
            velocities(:, i) = velocities(:, i) * T / masses(i)
        end do

        do i = 1, number_of_particles
            ! Assign every particle the type 1.
            types(i) = 1
        end do

        ! Aliasing system_size for brevity while printing to terminal.
        S = system_size

        print *, "╔════════════════════════════════════════════════════╗"
        print *, "║ Initializing a random state.                       ║"
        print *, "╚════════════════════════════════════════════════════╝"
        print *, "   Positions uniformely distributed inside the simulation box"
        print *, "   of size:"
        print *, "   [", S(1), S(2), S(3), "]"
        call remove_linear_momentum(velocities, masses)

    end subroutine random_initial_state

    subroutine fcc_initial_state(positions, velocities, masses, types)
        implicit none
        real    (real64), dimension(:,:), intent(in out) :: positions
        real    (real64), dimension(:,:), intent(in out) :: velocities
        real    (real64), dimension(:),   intent(in out) :: masses
        integer (int32),  dimension(:),   intent(in out) :: types        
        
        integer (int32) :: atoms_per_unit_cell = 4
        real (real64), dimension(:) :: lattice_vector(number_of_dimensions),   &
                                       unit_vector_x (number_of_dimensions),   &
                                       unit_vector_y (number_of_dimensions),   &
                                       unit_vector_z (number_of_dimensions)
        real (real64), dimension(:,:), allocatable :: atom_relative_position
        real (real64)   :: argon_mass = 39.948
        real (real64)   :: b
        integer (int32) :: i, j, k, l, atom_counter
    
        allocate(atom_relative_position(number_of_dimensions, atoms_per_unit_cell))

        ! Argon atomic mass.
        masses = argon_mass
        
        ! All particles are argon atoms.
        types = 1

        ! Set up the Cartesian unit vectors to help with creating the lattice.
        unit_vector_x = [1.0,   0.0,   0.0]
        unit_vector_y = [0.0,   1.0,   0.0]
        unit_vector_z = [0.0,   0.0,   1.0]

        ! Positions of the four atoms w.r.t. the unit cell origin.
        atom_relative_position(:, 1) = [0.0,  0.0,  0.0]
        atom_relative_position(:, 2) = [0.5,  0.5,  0.0]
        atom_relative_position(:, 3) = [0.5,  0.0,  0.5]
        atom_relative_position(:, 4) = [0.0,  0.5,  0.5]
        
        atom_relative_position = atom_relative_position * fcc_lattice_constant

        ! Aliasing the lattice constant for brevity in the expression inside the 
        ! loop.
        b = fcc_lattice_constant

        atom_counter = 1
        do i = 1, fcc_number_of_unit_cells
            do j = 1, fcc_number_of_unit_cells
                do k = 1, fcc_number_of_unit_cells
                    lattice_vector = i * unit_vector_x * b  +       &
                                     j * unit_vector_y * b  +       &
                                     k * unit_vector_z * b

                    do l = 1, atoms_per_unit_cell
                        positions(:, atom_counter) = lattice_vector + atom_relative_position(:, l)
                    end do
                end do
            end do
        end do

        call random_normal(velocities)
        velocities = velocities * temperature / argon_mass
    end subroutine fcc_initial_state

    subroutine setup_initial_state(positions, velocities, forces, masses, types)
        implicit none
        real    (real64), dimension(:,:), allocatable, intent(in out) :: positions
        real    (real64), dimension(:,:), allocatable, intent(in out) :: velocities
        real    (real64), dimension(:,:), allocatable, intent(in out) :: forces
        real    (real64), dimension(:),   allocatable, intent(in out) :: masses
        integer (int32),  dimension(:),   allocatable, intent(in out) :: types

        select case (initial_configuration)
        case ("fcc")             
            call print_fcc_warning()
            number_of_particles = 4 * fcc_number_of_unit_cells**3
            system_size_x = fcc_number_of_unit_cells * fcc_lattice_constant
            system_size_y = fcc_number_of_unit_cells * fcc_lattice_constant
            system_size_z = fcc_number_of_unit_cells * fcc_lattice_constant
            system_size   = [system_size_x, system_size_y, system_size_z]
            call allocate_arrays(positions, velocities, forces, masses, types)
            call fcc_initial_state(positions, velocities, masses, types)

        case ("random") 
            call allocate_arrays(positions, velocities, forces, masses, types)
            call random_initial_state(positions, velocities, masses, types)

        case default 
            print *, "╔════════════════════════════════════════════════════╗"
            print *, "║ Unknown initial configuration.                     ║"
            print *, "╚════════════════════════════════════════════════════╝"
            print *, "   The initial configuration"
            print *, "   <", initial_configuration, ">"
            print *, "   was not recognized. Defaulting to a random system."
            call random_initial_state(positions, velocities, masses, types)

        end select 
    end subroutine setup_initial_state

    subroutine allocate_arrays(positions, velocities, forces, masses, types)
        implicit none
        real    (real64), dimension(:,:), allocatable, intent(in out) :: positions
        real    (real64), dimension(:,:), allocatable, intent(in out) :: velocities
        real    (real64), dimension(:,:), allocatable, intent(in out) :: forces
        real    (real64), dimension(:),   allocatable, intent(in out) :: masses
        integer (int32),  dimension(:),   allocatable, intent(in out) :: types

        allocate(positions (number_of_dimensions, number_of_particles))
        allocate(velocities(number_of_dimensions, number_of_particles))
        allocate(forces    (number_of_dimensions, number_of_particles))
        allocate(masses    (number_of_particles))
        allocate(types     (number_of_particles))
    end subroutine allocate_arrays

    subroutine print_fcc_warning()
        implicit none 
        real (real64), dimension(:) :: S        (number_of_dimensions)

        ! If the FCC initial condition is chosen, we need to override the 
        ! number_of_particles in the input parameters (parameters module) 
        ! and use 4 (the number of atoms per FCC unit cell) times the number  
        ! of unit cells (cells in each dimension, cubed) instead.
        ! 
        ! The system size also needs to be set to the lattice constant times
        ! the number of unit cells in each dimension. Aliasing the new 
        ! system size to S for brevity when printing.
        S = [1.0, 1.0, 1.0] * fcc_number_of_unit_cells * fcc_lattice_constant 
        print *, "╔════════════════════════════════════════════════════╗"
        print *, "║                                _                   ║"
        print *, "║                               (_)                  ║"
        print *, "║      __      ____ _ _ __ _ __  _ _ __   __ _       ║"
        print *, "║      \ \ /\ / / _` | '__| '_ \| | '_ \ / _` |      ║"
        print *, "║       \ V  V / (_| | |  | | | | | | | | (_| |      ║"
        print *, "║        \_/\_/ \__,_|_|  |_| |_|_|_| |_|\__, |      ║"
        print *, "║                                         __/ |      ║"
        print *, "║                                        |___/       ║"
        print *, "╚════════════════════════════════════════════════════╝"
        print *, "   The specified number of particles in the input is"
        print *, "   ignored when using the FCC lattice initial state."
        print *, "   "
        print *, "   Specified number of atoms (disregarded):", number_of_particles
        print *, "   Number of atoms used:                   ", 4*fcc_number_of_unit_cells**3
        print *, "   "
        print *, "   The specified system size is also ignored."
        print *, "   "
        print *, "   Specified system size (disregarded):"
        print *, "   [", system_size_x, system_size_y, system_size_z, "]"
        print *, "   System size in use:"
        print *, "   [", S(1), S(2), S(3), "]"
    end subroutine
end module initial_states