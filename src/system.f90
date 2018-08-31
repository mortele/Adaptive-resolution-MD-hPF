module system
    use, intrinsic :: iso_fortran_env, only:  real64, int32
    use parameters,      only:  system_size_x,          &
                                system_size_y,          &
                                system_size_z,          &
                                number_of_particles,    &
                                number_of_dimensions
    implicit none
    private 

    real (real64), public, dimension(3) :: system_size = [-1.0, -1.0, -1.0]
    public ::   remove_linear_momentum,                 &
                apply_periodic_boundary_conditions,     &
                compute_distance_minimum_image

contains

    subroutine remove_linear_momentum(velocities, masses, silent)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: velocities
        real (real64), dimension(:),   intent(in)     :: masses
        logical, optional,             intent(in)     :: silent
        real (real64), dimension(:) :: total_momentum(number_of_dimensions), &
                                       P             (number_of_dimensions)
        integer (int32) :: i

        total_momentum = 0
        do i = 1, number_of_particles
            total_momentum = total_momentum + velocities(:, i) * masses(i)
        end do
        
        total_momentum = total_momentum / number_of_particles
        do i = 1, number_of_particles
            velocities(:, i) = velocities(:, i) - total_momentum / masses(i)
        end do
        
        ! Aliasing the total momentum of the center of mass for brevity when 
        ! priting to terminal.
        P = total_momentum * number_of_particles
        if (.not. silent) then
            print *, "╔════════════════════════════════════════════════════╗"
            print *, "║ Removing linear momentum.                          ║"
            print *, "╚════════════════════════════════════════════════════╝"
            print *, "   Removed total (center of mass momentum):"
            print *, "   [", P(1), P(2), P(3), "]"
        end if
    end subroutine remove_linear_momentum

    subroutine apply_periodic_boundary_conditions(positions)
        implicit none
        real (real64), dimension(:,:), intent(in out) :: positions
        integer (int32) :: i, j
        
        do i = 1, number_of_particles
            do j = 1, number_of_dimensions
                if (positions(j, i) >= system_size(j)) then
                    positions(j, i) = positions(j, i) - system_size(j)

                else if (positions(j, i) < 0) then
                    positions(j, i) = positions(j, i) + system_size(j)
                end if
            end do
        end do
    end subroutine apply_periodic_boundary_conditions

    pure function compute_distance_minimum_image(vector) result(minimum_image_vector)
        implicit none
        real (real64), dimension(:), intent(in) :: vector(number_of_dimensions)
        real (real64), dimension(:)             :: minimum_image_vector(number_of_dimensions)
        
        minimum_image_vector = vector - system_size * nint(vector / system_size)
    end function compute_distance_minimum_image

end module system
