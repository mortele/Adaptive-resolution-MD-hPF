module potential
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use particles,  only: masses, velocities
    use parameters, only: number_of_particles,      &
                          lennard_jones_epsilon,    &
                          lennard_jones_sigma,      &
                          lennard_jones_cutoff,     &
                          number_of_dimensions
    implicit none
    private 

    real (real64), private, parameter :: sigma6  = lennard_jones_sigma**6
    real (real64), private, parameter :: sigma12 = lennard_jones_sigma**12
    real (real64), public :: Ek, V, E  ! Kinetic, potential, and total energies.

    private :: lennard_jones_force, lennard_jones_potential
    public  :: compute_forces

contains

    subroutine compute_forces(positions, forces)
        implicit none
        real (real64), dimension(:,:), intent(in)      :: positions
        real (real64), dimension(:,:), intent(in out)  :: forces
        real (real64), dimension(number_of_dimensions) :: distance_vector,  &
                                                          force
        real (real64)   :: dr_squared
        integer (int32) :: i, j

        forces  = 0  ! All elements of dimension(3,:) array is set to zero. 
        V       = 0
        Ek      = 0

        do i = 1, number_of_particles
            do j = i + 1, number_of_particles
                distance_vector = positions(:,j) - positions(:,i)
                dr_squared      = dot_product(distance_vector, distance_vector)
                force           = lennard_jones_force(dr_squared)
                forces(:,i)     = forces(:,i) + force * distance_vector
                forces(:,j)     = forces(:,j) - force * distance_vector
                
                V = V + lennard_jones_potential(dr_squared)
            end do
        end do

        do i = 1, number_of_particles
            Ek = Ek + 0.5 * masses(i) * dot_product(velocities(:,i), velocities(:, i))
        end do

        E = V + Ek
    end subroutine compute_forces

    pure function lennard_jones_force(dr_squared) result(F_divided_by_r)
        implicit none
        real (real64), intent(in) :: dr_squared
        real (real64) :: F_divided_by_r
        real (real64) :: r2, r6, eps
        eps = lennard_jones_epsilon  ! Simply aliasing the variable for brevity.
        r2  = 1.0 / dr_squared
        r6  = r2 * r2 * r2
        F_divided_by_r = 24.0 * sigma6 * eps * r6  * (1.0 - 2.0 * r6 * sigma6) * r2
    end function lennard_jones_force

    pure function lennard_jones_potential(dr_squared) result(potential)
        implicit none
        real (real64), intent(in) :: dr_squared
        real (real64) :: potential
        real (real64) :: r2, r6, eps
        eps = lennard_jones_epsilon  ! Simply aliasing the variable for brevity.
        r2  = 1.0 / dr_squared
        r6  = r2 * r2 * r2
        potential = 4 * eps * sigma6 * r6 * (sigma6 * r6 - 1.0)
    end function lennard_jones_potential

end module potential