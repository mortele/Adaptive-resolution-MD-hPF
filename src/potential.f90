module potential
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,     only: compute_distance_minimum_image,   & ! Function 
                          system_size
    use particles,  only: masses, velocities
    use parameters, only: number_of_particles,              &
                          lennard_jones_epsilon,            &
                          lennard_jones_sigma,              &
                          lennard_jones_cutoff,             &
                          number_of_dimensions,             &
                          chi,                              &
                          kappa
    use field,      only: density_field,                    &
                          density_gradient,                 &
                          position_of_density_nodes,        &
                          interpolate_density_field,        & ! Function
                          interpolate_density_gradient,     & ! Subroutine
                          compute_density_field,            & ! Subroutine
                          compute_density_gradient            ! Subroutine
    implicit none
    private 

    real (real64), public :: sigma6, sigma12, cutoff_squared, potential_at_cutoff
    real (real64), public  :: Ek, V_md, V_hpf, E  ! Kinetic, potential, and total energies.

    public  :: lennard_jones_force,         &
               lennard_jones_potential,     &
               hpf_force,                   &
               hpf_potential,               &
               compute_forces_md,           &
               compute_forces_hpf

contains

    subroutine compute_forces_md(positions, forces)
        implicit none
        real (real64), dimension(:,:), intent(in)      :: positions
        real (real64), dimension(:,:), intent(in out)  :: forces
        real (real64), dimension(number_of_dimensions) :: distance_vector,  &
                                                          force
        real (real64)   :: dr_squared, r2, r6, r12, pi, volume, Vtail
        integer (int32) :: i, j
        
        sigma6  = lennard_jones_sigma**6
        sigma12 = lennard_jones_sigma**12
        
        ! Calculate the cutoff potential energy.
        cutoff_squared = lennard_jones_cutoff**2
        r2 = 1.0_real64 / cutoff_squared
        r6 = r2*r2*r2
        r12 = r6*r6
        potential_at_cutoff = 4.0_real64 * lennard_jones_epsilon * sigma6 * r6 * (sigma6 * r6 - 1.0_real64)

        forces  = 0.0_real64 ! All elements of dimension(3,:) array is set to zero. 
        force   = 0.0_real64
        V_md    = 0.0_real64
        Ek      = 0.0_real64

        do i = 1, number_of_particles
            do j = i + 1, number_of_particles
                distance_vector = positions(:,j) - positions(:,i)

                ! Apply the minimum image convention to compute the periodic 
                ! boundary condition distance between two particles.
                distance_vector = compute_distance_minimum_image(distance_vector)
                
                !dr_squared      = dot_product(distance_vector, distance_vector)
                dr_squared = distance_vector(1)**2 + distance_vector(2)**2 + distance_vector(3)**2
                if (dr_squared < cutoff_squared) then
                    force           = lennard_jones_force(dr_squared)
                    forces(:,i)     = forces(:,i) + force * distance_vector
                    forces(:,j)     = forces(:,j) - force * distance_vector
                    
                    V_md = V_md + lennard_jones_potential(dr_squared)
                end if
            end do
        end do

        do i = 1, number_of_particles
            Ek = Ek + 0.5 * masses(i) * dot_product(velocities(:,i), velocities(:, i))
        end do

        ! Tail corrections for the Lennard-Jones cutoff
        pi = 3.141592653589793238462643383279502884197169399375105820974_real64
        volume = system_size(1) * system_size(2) * system_size(3)
        Vtail = (8.0_real64 / 3.0_real64) * pi * lennard_jones_epsilon * lennard_jones_sigma**3 * number_of_particles / volume * ((1.0_real64 / 3.0_real64) * (lennard_jones_sigma/lennard_jones_cutoff)**9 - (lennard_jones_sigma/lennard_jones_cutoff)**3)
        !print *, Vtail
        !V = V + Vtail
        E = V_md + Ek
    end subroutine compute_forces_md

    subroutine compute_forces_hpf(positions, forces_hpf)
        implicit none
        real (real64), dimension(:,:), intent(in)      :: positions
        real (real64), dimension(:,:), intent(in out)  :: forces_hpf

        integer (int32) :: i

        do i = 1, number_of_particles
            forces_hpf(:,i) = hpf_force(positions(:,i))
        end do

        do i = 1, number_of_particles
            Ek = Ek + 0.5 * masses(i) * dot_product(velocities(:,i), velocities(:, i))
        end do

    end subroutine compute_forces_hpf

    pure function lennard_jones_force(dr_squared) result(F_divided_by_r)
        implicit none
        real (real64), intent(in) :: dr_squared
        real (real64) :: F_divided_by_r
        real (real64) :: r2, r6, eps
        eps = lennard_jones_epsilon  ! Simply aliasing the variable for brevity.
        r2  = 1.0_real64 / dr_squared
        r6  = r2 * r2 * r2
        F_divided_by_r = 24.0_real64 * sigma6 * eps * r6  * (1.0_real64 - 2.0_real64 * r6 * sigma6) * r2
    end function lennard_jones_force

    pure function lennard_jones_potential(dr_squared) result(potential)
        implicit none
        real (real64), intent(in) :: dr_squared
        real (real64) :: potential
        real (real64) :: r2, r6, eps
        eps = lennard_jones_epsilon  ! Simply aliasing the variable for brevity.
        r2  = 1.0_real64 / dr_squared
        r6  = r2 * r2 * r2
        potential = 4.0_real64 * eps * sigma6 * r6 * (sigma6 * r6 - 1.0_real64) - potential_at_cutoff
    end function lennard_jones_potential

    function hpf_force(position) result(force)
        implicit none
        real (real64), intent(in), dimension(:) :: position
        real (real64), dimension(number_of_dimensions) :: force
        
        call interpolate_density_gradient(density_gradient,          &
                                          position_of_density_nodes, &
                                          position,                  &
                                          force)
        force = -(1.0_real64 / kappa) * force
    end function hpf_force

    function hpf_potential(position) result(potential)
        implicit none
        real (real64), intent(in), dimension(:) :: position
        real (real64) :: potential

        potential = interpolate_density_field(density_field,                &
                                              position_of_density_nodes,    &
                                              position)
    end function hpf_potential

end module potential