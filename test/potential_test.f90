module potential_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use potential,  only:   lennard_jones_force,        &
                            lennard_jones_potential,    &
                            compute_forces,             &
                            sigma6,                     &   
                            sigma12
    use fruit,      only:   assert_equals,              & ! Subroutine
                            assert_false                  ! Subroutine
    use parameters, only:   lennard_jones_sigma,        &
                            lennard_jones_epsilon,      &
                            lennard_jones_cutoff,       &
                            number_of_particles,        &
                            number_of_dimensions
    use system,     only:   system_size
    implicit none
    private

    public ::   setup,                          &
                teardown,                       &
                test_lennard_jones_force,       &
                test_lennard_jones_potential


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

end module potential_test