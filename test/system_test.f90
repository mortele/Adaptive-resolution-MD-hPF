module system_test
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use fruit,              only:   assert_equals,                      & ! Subroutine
                                    assert_false                          ! Subroutine
    use parameters,         only:   system_size_x,                      &
                                    system_size_y,                      &
                                    system_size_z,                      &
                                    number_of_particles,                &
                                    number_of_dimensions               
    use particles,          only:   positions,                          &
                                    velocities,                         &
                                    masses
    use system,             only:   remove_linear_momentum,             & ! Subroutine
                                    apply_periodic_boundary_conditions, & ! Subroutine
                                    compute_distance_minimum_image,     & ! Subroutine
                                    system_size
    use random_generator,   only:   random_normal
    implicit none
    private

    private ::  setup_test_system
    public  ::  setup,                                  &
                teardown,                               &
                test_periodic_boundary_conditions,      &
                test_distance_minimum_image,            &
                test_remove_linear_momentum
contains

subroutine setup
end subroutine setup

subroutine teardown
end subroutine teardown


subroutine test_distance_minimum_image()
    logical, allocatable, dimension(:,:,:) :: distances_too_long    
    real (real64), allocatable, dimension(:,:,:) :: distances
    logical :: any_distance_too_long
    integer (int32) :: i, j, test_number

    call setup_test_system()
    allocate(distances (number_of_dimensions, number_of_particles,  number_of_particles))
    allocate(positions (number_of_dimensions, number_of_particles))
    positions(:,1) = [1.0,  5.0,  6.0]
    positions(:,2) = [9.0,  5.0,  6.0]
    positions(:,3) = [9.0,  5.0,  0.5]
    positions(:,4) = [9.0,  9.0,  9.0]
    positions(:,5) = [1.0,  1.0,  1.0]
    positions(:,6) = [6.0,  2.0,  0.5]
    call random_number(positions(:,7:))
    positions(:,7:) = positions(:,7:) * system_size_x


    do i = 1, number_of_particles
        do j = 1, number_of_particles
            if (i /= j) then
                distances(:,i,j) = compute_distance_minimum_image(positions(:,i) - positions(:,j))
            else
                distances(:,i,j) = 0.0
            end if
        end do
    end do

    ! positions(:,1) - positions(:,2) = [1,5,6] - [9,5,6] = [-8,0,0] -> [2,0,0]
    call assert_equals([ 2.0_real64,  0.0_real64,  0.0_real64], distances(:,1,2), 3, "1  test_distance_minimum_image : Computed distance across boundary of the simulation box incorrect")
    ! positions(:,2) - positions(:,3) = [9,5,6] - [9,5,0.5] = [0,0,5.5] -> [0,0,-4.5]
    call assert_equals([ 0.0_real64,  0.0_real64, -4.5_real64], distances(:,2,3), 3, "2  test_distance_minimum_image : Computed distance across boundary of the simulation box incorrect")
    ! positions(:,4) - positions(:,5) = [9,9,9] - [1,1,1] = [8,8,8] -> [-2,-2,-2]
    call assert_equals([-2.0_real64, -2.0_real64, -2.0_real64], distances(:,4,5), 3, "3  test_distance_minimum_image : Computed distance across boundary of the simulation box incorrect")
    ! positions(:,4) - positions(:,6) = [9,9,9] - [6,2,0.5] = [3,7,8.5] -> [3,-2,-3.5]
    call assert_equals([ 3.0_real64, -3.0_real64, -1.5_real64], distances(:,4,6), 3, "4  test_distance_minimum_image : Computed distance across boundary of the simulation box incorrect")

    test_number = 5
    do i = 1, number_of_particles
        do j = 1, number_of_particles
            if (i /= j) then
                call assert_equals(distances(:,i,j), -distances(:,j,i), 3, char(test_number)//"  test_distance_minimum_image : Distance between i and j is not equal to negative distance between j and i")
                test_number = test_number + 1
            end if
        end do
    end do

    allocate(distances_too_long(number_of_dimensions, number_of_particles, number_of_particles))

    distances_too_long = distances > 5.0
    any_distance_too_long = any(distances_too_long)
    call assert_false(any_distance_too_long, char(test_number)//" test_distance_minimum_image : Some distances were computed to be larger than half the system size")
    test_number = test_number + 1

    distances_too_long = distances < -5.0
    any_distance_too_long = any(distances_too_long)
    call assert_false(any_distance_too_long, char(test_number)//" test_distance_minimum_image : Some distances were computed to be smaller than minus half the system size")
    
    deallocate(positions)
    deallocate(distances_too_long)
end subroutine test_distance_minimum_image

subroutine test_periodic_boundary_conditions()
    logical, allocatable, dimension(:,:) :: outside_box
    logical :: any_outside
    
    call setup_test_system()
    
    allocate(positions (number_of_dimensions, number_of_particles))
    allocate(outside_box(number_of_dimensions, number_of_particles))
    positions(:,1)  = [1.0,   2.0,   3.0]
    positions(:,2)  = [11.0,  4.0,   5.0]
    positions(:,3)  = [12.0,  13.0,  6.0]
    positions(:,4)  = [7.0,   14.0, 15.0]
    positions(:,5)  = [16.0,  17.0, 18.0]
    positions(:,6)  = [-1.0,  8.0,  19.0]
    positions(:,7)  = [-2.0,  11.5,  9.0]
    positions(:,8)  = [12.5,  13.5, -3.0]
    positions(:,9)  = [-4.0,  -5.0, -6.0]
    positions(:,10) = [-7.0,  14.5, -8.0]
    positions(:,11) = [-9.0,  1.5,   2.5]
    call apply_periodic_boundary_conditions(positions)
    
    call assert_equals([1.0_real64,   2.0_real64,   3.0_real64], positions(:,1),  3, "1  test_periodic_boundary_conditions : No changes expected")
    call assert_equals([1.0_real64,   4.0_real64,   5.0_real64], positions(:,2),  3, "2  test_periodic_boundary_conditions : Change expected in first component")
    call assert_equals([2.0_real64,   3.0_real64,   6.0_real64], positions(:,3),  3, "3  test_periodic_boundary_conditions : Change expected in first and second components")
    call assert_equals([7.0_real64,   4.0_real64,   5.0_real64], positions(:,4),  3, "4  test_periodic_boundary_conditions : Change expected in second and third components")
    call assert_equals([7.0_real64,   4.0_real64,   5.0_real64], positions(:,4),  3, "5  test_periodic_boundary_conditions : Change expected in second and third components")
    call assert_equals([6.0_real64,   7.0_real64,   8.0_real64], positions(:,5),  3, "6  test_periodic_boundary_conditions : Change expected in all components")
    call assert_equals([9.0_real64,   8.0_real64,   9.0_real64], positions(:,6),  3, "6  test_periodic_boundary_conditions : Change expected in first and third components")
    call assert_equals([8.0_real64,   1.5_real64,   9.0_real64], positions(:,7),  3, "7  test_periodic_boundary_conditions : Change expected in first and second components")        
    call assert_equals([2.5_real64,   3.5_real64,   7.0_real64], positions(:,8),  3, "8  test_periodic_boundary_conditions : Change expected in all components")                
    call assert_equals([6.0_real64,   5.0_real64,   4.0_real64], positions(:,9),  3, "9  test_periodic_boundary_conditions : Change expected in all components")        
    call assert_equals([3.0_real64,   4.5_real64,   2.0_real64], positions(:,10), 3, "10 test_periodic_boundary_conditions : Change expected in all components")        
    call assert_equals([1.0_real64,   1.5_real64,   2.5_real64], positions(:,11), 3, "11 test_periodic_boundary_conditions : Change expected in first component")        
    
    call random_number(positions)
    positions  = positions * 30.0 - 10.0 ! Uniform distribution in [-10, 20)
    call apply_periodic_boundary_conditions(positions)
    outside_box = positions > 10.0
    any_outside = any(outside_box)
    call assert_false(any_outside, "12 test_periodic_boundary_conditions : One or more of the random positions ended up outside the simulation box")
    
    outside_box = positions < 0.0
    any_outside = any(outside_box)
    call assert_false(any_outside, "13 test_periodic_boundary_conditions : One or more of the random positions ended up outside the simulation box")
    
    deallocate(positions)
    deallocate(outside_box)
end subroutine test_periodic_boundary_conditions

subroutine test_remove_linear_momentum()
    real (real64), dimension(:), allocatable :: total_momentum
    real (real64)   :: test_tollerance = 1e-10
    integer (int32) :: i, j, k, test_number

    call setup_test_system()
    allocate(velocities     (number_of_dimensions, number_of_particles))
    allocate(masses         (number_of_particles))
    allocate(total_momentum (number_of_dimensions))

    test_number = 1
    do k = 1, 2
        do i = 1, 10
            call random_number(masses)
            masses = masses * 9.9 + 0.1  ! Uniform distribution in [0.1, 10.0)
            
            if (k == 1) then
                call random_normal(velocities)
            else 
                call random_number(velocities)
                velocities = velocities * 20 - 10
            end if
            
            call remove_linear_momentum(velocities, masses, .true.)
            total_momentum = 0
            do j = 1, number_of_particles
                total_momentum = total_momentum + masses(j) * velocities(:,j)
            end do

            call assert_equals([0.0_real64, 0.0_real64, 0.0_real64], total_momentum, 3, test_tollerance, char(test_number)//"  test_remove_linear_momentum : The total computed linear momentum is not zero")
            test_number = test_number + 1
        end do
    end do
end subroutine test_remove_linear_momentum

subroutine setup_test_system
    system_size_x = 10.0
    system_size_y = 10.0
    system_size_z = 10.0
    system_size = [system_size_x, system_size_y, system_size_z]
    number_of_dimensions = 3
    number_of_particles  = 11
end subroutine setup_test_system

end module system_test