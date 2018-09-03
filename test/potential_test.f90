module potential_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use potential,  only:   lennard_jones_force,        &
                            lennard_jones_potential,    &
                            compute_forces,             &
                            sigma6,                     &   
                            sigma12,                    &
                            Ek,                         &
                            V,                          &
                            E
    use fruit,      only:   assert_equals,              & ! Subroutine
                            assert_false                  ! Subroutine
    use parameters, only:   lennard_jones_sigma,        &
                            lennard_jones_epsilon,      &
                            lennard_jones_cutoff,       &
                            number_of_particles,        &
                            number_of_dimensions
    use particles,  only:   positions,                  &
                            velocities,                 &
                            forces,                     &
                            masses
    use system,     only:   system_size
    implicit none
    private

    public ::   setup,                          &
                teardown,                       &
                test_lennard_jones_force,       &
                test_lennard_jones_potential,   &
                test_compute_forces


contains
    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_lennard_jones_force()
        implicit none
        real (real64)   :: r, expected_force, computed_force, tollerance
        integer (int32) :: i
        character (10)  :: r_string

        tollerance            = 1e-10
        lennard_jones_epsilon = 1.25875
        lennard_jones_sigma   = 3.29852
        sigma6  = lennard_jones_sigma**6
        sigma12 = lennard_jones_sigma**12


        do i = 1, 20
            call random_number(r)
            r = r * 2.0 * lennard_jones_sigma + 0.75 * lennard_jones_sigma*2.0**(1.0/6.0)   ! Uniform random in range [0.75σ, 2σ)
            expected_force = 24.0 * lennard_jones_epsilon * lennard_jones_sigma**6 / r**7 - 48.0 * lennard_jones_epsilon * lennard_jones_sigma**12 / r**13
            computed_force = lennard_jones_force(r**2) * r
            write(r_string,fmt="(f10.3)") r
            call assert_equals(expected_force, computed_force, tollerance, char(i)//" test_lennard_jones_potential : Computed LJ force is not correct for distance r = "//r_string)
        end do

    end subroutine test_lennard_jones_force

    subroutine test_lennard_jones_potential()
        implicit none
        real (real64)   :: r, expected_potential, computed_potential, tollerance
        integer (int32) :: i
        character (10)  :: r_string

        tollerance            = 1e-10
        lennard_jones_epsilon = 1.25875
        lennard_jones_sigma   = 3.29852
        sigma6  = lennard_jones_sigma**6
        sigma12 = lennard_jones_sigma**12


        do i = 1, 20
            call random_number(r)
            r = r * 2.0 * lennard_jones_sigma + 0.75 * lennard_jones_sigma*2.0**(1.0/6.0)   ! Uniform random in range [0.75σ, 2σ)
            expected_potential = 4.0 * lennard_jones_epsilon * ((lennard_jones_sigma/r)**12 - (lennard_jones_sigma/r)**6)
            computed_potential = lennard_jones_potential(r**2)
            write(r_string,fmt="(f10.3)") r
            call assert_equals(expected_potential, computed_potential, tollerance, char(i)//" test_lennard_jones_potential : Computed LJ potential is not correct for distance r = "//r_string)
        end do

    end subroutine test_lennard_jones_potential

    subroutine test_compute_forces()
        implicit none
        real (real64), dimension(number_of_dimensions) :: unit_vector
        real (real64)   :: tollerance

        tollerance            = 1e-10
        lennard_jones_epsilon = 1.25875
        lennard_jones_sigma   = 3.29852

        ! Tests involving only a single particle.
        number_of_particles = 1
        allocate(positions (number_of_dimensions,number_of_particles))
        allocate(velocities(number_of_dimensions,number_of_particles))
        allocate(forces    (number_of_dimensions,number_of_particles))
        allocate(masses    (number_of_particles))
        masses = 1.0

        positions (:,1) = [3.0, 3.0, 3.0]
        velocities(:,1) = [1.0, 2.0, 3.0]
        call compute_forces(positions, forces)
        call assert_equals([0.0_real64, 0.0_real64, 0.0_real64], forces(:,1), 3, tollerance, "1  test_compute_forces : Computed net force on a system of a single particle is not zero")
        call assert_equals(0.0_real64, V, tollerance, "2  test_compute_forces : Compute net potential on a system of a single particle is not zero")
        call assert_equals(0.5_real64*((1.0)**2 + (2.0)**2 + (3.0)**2), Ek, tollerance, "3  test_compute_forces : Computed kinetic energy of a single particle is not correct")
        call assert_equals(E, Ek, tollerance, "4  test_compute_forces : Computed total energy is not equal to the kinetic energy for a system of a single particle.")
        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)

        ! Three particles placed on a line, the forces from the outer two should
        ! exactly cancel and make the net force on the middle one zero.
        number_of_particles = 3
        allocate(positions (number_of_dimensions,number_of_particles))
        allocate(velocities(number_of_dimensions,number_of_particles))
        allocate(forces    (number_of_dimensions,number_of_particles))
        allocate(masses    (number_of_particles))
        masses = 1.0

        positions(:,1) = [2.0, 0.0, 0.0]
        positions(:,2) = [4.0, 0.0, 0.0]
        positions(:,3) = [6.0, 0.0, 0.0]
        call random_number(velocities)
        call compute_forces(positions, forces)
        call assert_equals(0.0_real64, forces(1,2), tollerance, "1  test_compute_forces : Computed net force should be zero when affected by only two particles on opposite sides at identical distance")


        ! Three particles placed in a triangle, the middle particle should feel
        ! a net force in exactly two directions.
        positions(:,1) = [2.0, 2.0, 2.0]
        positions(:,2) = [4.0, 2.0, 2.0]
        positions(:,3) = [2.0, 4.0, 2.0]
        call compute_forces(positions, forces)
        call assert_equals(0.0_real64, forces(3,1), tollerance, "2  test_compute_forces : Computed net force in the z direction should be zero when all particles have the same z coordinate")
        call assert_equals(lennard_jones_force((2.0_real64)**2)*2.0_real64, forces(1,1), tollerance, "3  test_compute_forces : Computed net force in the x direction should only be the force from the second particle, since the third has the same x coordinate as the first")
        call assert_equals(lennard_jones_force((2.0_real64)**2)*2.0_real64, forces(2,1), tollerance, "4  test_compute_forces : Computed net force in the y direction should only be the force from the third particle, since the second has the same y coordinate as the first")
        unit_vector = [-1.0, -1.0, 0.0]
        unit_vector = unit_vector / norm2(unit_vector)
        call assert_equals(unit_vector, forces(:,1)/norm2(forces(:,1)), 3, tollerance, "5  test_compute_forces : The computed net force does not have the correct overall direction")
        deallocate(positions)
        deallocate(velocities)
        deallocate(forces)
        deallocate(masses)

    end subroutine test_compute_forces
end module potential_test