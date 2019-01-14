module sampler
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use parameters, only:   number_of_time_steps
    use particles,  only:   positions,             &
                            forces
    use potential,  only:   compute_forces_md,     &
                            Ek,                    & ! Current total kinetic energy.
                            V_md,                  & ! Current total potential energy. 
                            E                        ! Current total energy.
    implicit none
    private

    real (real64), public, dimension(:), allocatable ::  kinetic_energy,     &
                                                         potential_energy,   &
                                                         total_energy

    public :: store_energy

contains

    subroutine store_energy(kinetic_energy, potential_energy, total_energy)
        implicit none
        real (real64), dimension(:), intent(in out) :: kinetic_energy
        real (real64), dimension(:), intent(in out) :: potential_energy
        real (real64), dimension(:), intent(in out) :: total_energy
        integer (int32), save :: array_index = 1
        logical,         save :: first_call  = .true.

        ! If this is the first time the subroutine is called, we have to 
        ! compute the energy. In subsequent calls, this is already done in the 
        ! potential module, so we just fetch it from there and avoid re-
        ! computing it.
        if (first_call) then
            call compute_forces_md(positions, forces)
        end if

        kinetic_energy  (array_index) = Ek
        potential_energy(array_index) = V_md
        total_energy    (array_index) = E

        array_index = array_index + 1
    end subroutine store_energy

end module sampler