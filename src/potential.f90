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

    real (real64), public :: Ek, V, E  ! Kinetic, potential, and total energies.

    private :: lennard_jones_force, lennard_jones_potential
    public  :: compute_forces

contains

    subroutine compute_forces(positions, forces)
        real (real64), dimension(:,:), intent(in)      :: positions
        real (real64), dimension(:,:), intent(in out)  :: forces
        real (real64), dimension(number_of_dimensions) :: distance_vector,  &
                                                          force
        real (real64)   :: dr
        integer (int32) :: i, j

        forces  = 0  ! All elements of dimension(3,:) array is set to zero. 
        V       = 0
        Ek      = 0

        do i = 1, number_of_particles
            do j = i + 1, number_of_particles
                distance_vector = positions(:,j) - positions(:,i)
                dr              = dot_product(distance_vector, distance_vector)
                force           = lennard_jones_force(dr)
                forces(:,i)     = forces(:,i) + force
                forces(:,j)     = forces(:,j) - force
                
                V  = V + lennard_jones_potential(dr)
            end do
        end do

        do i = 1, number_of_particles
            Ek = Ek + 0.5 * masses(i) * dot_product(velocities(:,i), velocities(:, i))
        end do

        E = V + Ek
    end subroutine compute_forces

    pure function lennard_jones_force(r) result(F_divided_by_r)
        implicit none
        real (real64), intent(in)  :: r
        real (real64) :: F_divided_by_r
        real (real64) :: r2, r6

        F_divided_by_r = 0
    end function lennard_jones_force

    pure function lennard_jones_potential(r) result(potential)
        implicit none
        real (real64), intent(in)  :: r
        real (real64) :: potential
        real (real64) :: r2, r6
        
        potential = 0
    end function lennard_jones_potential

end module potential