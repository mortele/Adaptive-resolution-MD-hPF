module field_test
    use, intrinsic :: iso_fortran_env,  only: int32, real64
    use fruit,      only:   assert_equals,                      &
                            assert_true,                        &
                            assert_false
    use system,     only:   system_size
    use parameters, only:   number_of_field_nodes_x,            &
                            number_of_field_nodes_y,            &
                            number_of_field_nodes_z,            &
                            number_of_dimensions,               &
                            number_of_particles,                &
                            system_size_x,                      &
                            system_size_y,                      &
                            system_size_z
    use field,      only:   allocate_field_arrays,              & ! Submodule
                            compute_density_field,              & ! Submodule
                            density_field,                      &
                            position_of_density_nodes,          &
                            density_gradient,                   &
                            number_of_field_nodes,              &
                            compute_density_gradient              ! Subroutine
    
    use particles,  only:   positions,                          &
                            masses
    implicit none
    private 

    public  ::  setup,                          &
                teardown,                       &
                test_compute_density_field,     &
                test_allocate_field_arrays,     &
                test_compute_density_gradient
contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_compute_density_field()
        integer (int32) :: i, j, k

        ! We start off with a 1D test of the density field computation, with a 
        ! single particle.
        number_of_field_nodes_x = 3
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 3.0_real64
        system_size_y         = 3.0_real64
        system_size_z         = 3.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        ! Note: The density nodes we calculate the density at are now positioned
        ! at:
        !
        ! 1   [0, 0, 0]
        ! 2   [1, 0, 0]
        ! 3   [2, 0, 0]

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses   (number_of_particles))
        masses         = 1.0
        positions(:,1) = [0.5,0.0,0.0]
        call compute_density_field(positions, masses)
        
        ! Placing the particle in the middle of two field nodes must allocate 
        ! half of the mass of the particle to each node.
        call assert_equals(0.5_real64, density_field(1,1,1), "1  test_compute_density_field : A particle placed in the middle between field nodes 1 and 2 should distribute its mass evenly across the nodes")
        call assert_equals(0.5_real64, density_field(2,1,1), "2  test_compute_density_field : A particle placed in the middle between field nodes 1 and 2 should distribute its mass evenly across the nodes")
        
        ! Checking the same thing, but at the right hand side edge of the 
        ! system box, making sure the periodic boundary conditions are working
        ! correctly.
        positions(:,1) = [2.5,0.0,0.0]
        call compute_density_field(positions, masses)
        call assert_equals(0.5_real64, density_field(1,1,1), "3  test_compute_density_field : A particle placed in the middle between field node 3 and the right hand side edge of the system box should distribute its mass evenly between node 3 and node 1")
        call assert_equals(0.5_real64, density_field(3,1,1), "4  test_compute_density_field : A particle placed in the middle between field node 3 and the right hand side edge of the system box should distribute its mass evenly between node 3 and node 1")
        
        ! Next, we make sure that a single particle of mass 1.0 gives a grand 
        ! total contribution to the overall density field of exactly 1.0.
        do i = 1, 20
            call random_number(positions)
            positions(2:,1) = [0.0, 0.0]
            positions = positions * system_size_x
            call compute_density_field(positions, masses)
            call assert_equals(1.0_real64, sum(density_field), char(4+i)//" test_compute_density_field : The contribution of a single particle of mass 1.0 should be exactly 1.0 to the overall density field")
        end do

        ! Check next that a particle placed exactly on top of a density field 
        ! vertex should contribute all its mass to *only* this vertex, and give
        ! no contribution to any other vertices.
        do i = 1, 3
            positions(:,1) = [real(i-1), 0.0, 0.0]
            call compute_density_field(positions, masses)
            call assert_equals(1.0_real64, density_field(i,1,1), char(i)//"  test_compute_density_field : A particle placed on top of a density field vertex should contribute all its mass to that vertex, and nothing to any other vertex in the field")
        end do

        ! We now consider a 2D case with a single particle.
        number_of_field_nodes_x = 3
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 3.0_real64
        system_size_y         = 3.0_real64
        system_size_z         = 3.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        ! When placed right in between the four nodes at 
        !
        ! 1   [0, 0, 0]
        ! 2   [1, 0, 0]
        ! 3   [0, 1, 0]
        ! 4   [1, 1, 0]
        !
        ! the mass of our single particle should be distributed evenly across 
        ! all of them.
        positions(:,1) = [0.5, 0.5, 0.0]
        call compute_density_field(positions, masses)
        do i = 1, 2
            do j = 1,2
                call assert_equals(0.25_real64, density_field(i,j,1), "30 test_compute_density_field : A particle in 2D placed in the middle of four density vertices should evenly distribute 25% of its mass onto each of these vertices, and contribute nothing to any others")
            end do
        end do

        ! Next, check that the same holds at the very right hand side of the 
        ! simulation box, i.e. that the periodic boundary conditions also work
        ! correctly in 2D.
        positions(:,1) = [2.5, 2.5, 0.0]
        call compute_density_field(positions, masses)
        call assert_equals(0.25_real64, density_field(1,1,1), "40 test_compute_density_field : A particle in 2D placed in the middle of four density vertices on the right hand side of the simulation box should evenly distribute 25% of its mass onto each of these vertices, two of which should be the vertices with index 1 because of periodic boundary conditions")
        call assert_equals(0.25_real64, density_field(3,1,1), "41 test_compute_density_field : A particle in 2D placed in the middle of four density vertices on the right hand side of the simulation box should evenly distribute 25% of its mass onto each of these vertices, two of which should be the vertices with index 1 because of periodic boundary conditions")
        call assert_equals(0.25_real64, density_field(1,3,1), "42 test_compute_density_field : A particle in 2D placed in the middle of four density vertices on the right hand side of the simulation box should evenly distribute 25% of its mass onto each of these vertices, two of which should be the vertices with index 1 because of periodic boundary conditions")
        call assert_equals(0.25_real64, density_field(3,3,1), "43 test_compute_density_field : A particle in 2D placed in the middle of four density vertices on the right hand side of the simulation box should evenly distribute 25% of its mass onto each of these vertices, two of which should be the vertices with index 1 because of periodic boundary conditions")

        ! Now we move on to full 3D tests (still with a single particle). 
        number_of_field_nodes_x = 3
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 3.0_real64
        system_size_y         = 3.0_real64
        system_size_z         = 3.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        ! Lets check that a particle placed in the center of the cube consisting
        ! of the 8 closest density field vertices distribute its mass evenly
        ! onto all 8.
        positions(:,1) = [0.5, 0.5, 0.5]
        call compute_density_field(positions, masses)
        do i = 1, 2
            do j = 1,2
                do k = 1,2
                    call assert_equals(1.0_real64 / 8.0_real64, density_field(i,j,k), "50 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices")
                end do
            end do
        end do

        ! Redo the same test, but now on the right hand side edge of the 
        ! simulation box. In this way we check that the periodic boundary 
        ! conditions on the density field works in the full 3D configuration.
        positions(:,1) = [2.5, 2.5, 2.5]
        call compute_density_field(positions, masses)
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,1), "60 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,1), "61 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,1), "62 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,3), "63 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,1), "64 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,3), "65 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,3), "66 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,3), "67 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box")

        ! We check that a particle placed on top of a density vertex contributes
        ! its mass *only* to that single vertex in 3D.
        do i = 1, number_of_field_nodes_x
            do j = 1, number_of_field_nodes_y
                do k = 1, number_of_field_nodes_z
                    positions(:,1) = [real(i-1), real(j-1), real(k-1)]
                    call compute_density_field(positions, masses)
                    call assert_equals(1.0_real64, density_field(i,j,k), "70 test_compute_density_field : A particle placed on top of a density field vertex should contribute its mass to that single vertex and not to any other vertices")
                    call assert_equals(0.0_real64, sum(density_field(i+1:,j+1:,k+1:))+sum(density_field(:i-1,:j-1,:k-1)), "70 test_compute_density_field : A particle placed on top of a density field vertex should contribute its mass to that single vertex and not to any other vertices")
                end do
            end do
        end do
        
        deallocate(positions)
        deallocate(masses)
        deallocate(density_field)
        deallocate(density_gradient)
        deallocate(position_of_density_nodes)
    end subroutine test_compute_density_field

    subroutine test_compute_density_gradient()
        integer (int32) :: i, j, k
        real (real64)   :: pi, x, real_i, real_number_of_field_nodes
        real (real64)   :: tollerance, h
        logical         :: any_nonzero

        number_of_dimensions  = 3
        number_of_particles   = 1
        number_of_field_nodes_x = 100
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        system_size_x         = 1.0_real64
        system_size_y         = 1.0_real64
        system_size_z         = 1.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]
        allocate(density_field      (                     number_of_field_nodes_x,number_of_field_nodes_y,number_of_field_nodes_z))
        allocate(density_gradient   (number_of_dimensions,number_of_field_nodes_x,number_of_field_nodes_y,number_of_field_nodes_z))

        ! We first test on a constant density field configuration. The gradient
        ! should vanish everywhere.
        density_field = 3.5
        call compute_density_gradient(density_field)
        any_nonzero = any(density_gradient /= 0.0)
        call assert_false(any_nonzero, "1  test_compute_density_gradient : At least one element of density_gradient was computed to be non-zero for a constant density field configuration")

        ! Secondly, we see if a 1D linear density is differentiated correctly.
        real_number_of_field_nodes  = real(number_of_field_nodes_x)
        do i = 1, number_of_field_nodes_x
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            density_field(i,:,:) = 1.25_real64 * x
        end do
        call compute_density_gradient(density_field)

        ! Excluding the boundary points, since the density field has a 
        ! discontinuity at the x=0 and x=1 edges.
        do i = 2, number_of_field_nodes_x-1, 10
            call assert_equals(1.25_real64, density_gradient(1,i,1,1), 1e-13_real64, "2  test_compute_density_gradient : The gradient of a linear density field configuration was not calculated correctly.")
        end do

        ! We now try a second order polynomial.
        system_size_x         = 3.0_real64
        system_size_y         = 3.0_real64
        system_size_z         = 3.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]
        do i = 1, number_of_field_nodes_x
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            density_field(i,:,:) = (x - 1.5_real64)**2 + 0.75_real64 * x
        end do
        call compute_density_gradient(density_field)

        ! Excluding the boundary points, since the density field has a 
        ! discontinuity at the x=0 and x=1 edges.
        do i = 2, number_of_field_nodes_x-1, 10
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            call assert_equals(2.0_real64*(x-1.5_real64) + 0.75_real64, density_gradient(1,i,1,1), 1e-13_real64, "3  test_compute_density_gradient : The gradient of a linear density field configuration was not calculated correctly.")
        end do


        ! Next we test if compute_density_gradient() can differentiate a 
        ! sinusoidal function in 1D.
        pi = 3.141592653589793238462643383279502884197169399375105820974
        system_size_x         = 2.0_real64*pi
        system_size_y         = 2.0_real64*pi
        system_size_z         = 2.0_real64*pi
        system_size           = [system_size_x, system_size_y, system_size_z]
        do i = 1, number_of_field_nodes_x
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            density_field(i,:,:) = sin(x)
        end do
        call compute_density_gradient(density_field)

        ! The error in the central difference approximation used is bounded by
        !     1                  1           2
        !    ─── R (x ± h)   =  ─── f'''(ξ) h,      x-h < ξ < x+h .
        !    2 h  2              6
        !
        ! Since the derivatives of sines and cosines are themselves sines and 
        ! consines, the f'''(ξ) term is bounded in absolute value by 1. This 
        ! means that the error in the approximation should be less than 
        !                                           
        !                  f(x+h) - f(x-h)           
        !    f'(x)   =    ─────────────────  +  ε,   
        !                        2 h                
        !                2           
        !               h     1  ╭ 2 π ╮2
        !          ε ∝ ─── = ─── │ ─── │  ≈ 0.00065797 ≈ 6.6e-4
        !               6     6  ╰ 100 ╯ 
        !
        tollerance = 6.6e-4_real64
        do i = 1, number_of_field_nodes_x,10
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            call assert_equals(cos(x), density_gradient(1,i,1,1), tollerance, "4  test_compute_density_gradient : Gradient of the 1D sinusoidal density field configuration not calculated correctly")
        end do

        ! Next, lets test an exponential.
        system_size_x         = 1.0_real64
        system_size_y         = 1.0_real64
        system_size_z         = 1.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]
        do i = 1, number_of_field_nodes_x
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            density_field(i,:,:) = exp(x)
        end do
        call compute_density_gradient(density_field)

        ! In this case the f'''(x) term is just an exponential, so the error 
        ! term takes the form 
        !            2
        !           h   ξ       1  ╭  1  ╮2  ξ                    ξ
        !    ε  ∝  ─── e    =  ─── │ ─── │  e  ≈ 1.6666667e-05 * e,
        !           6           6  ╰ 100 ╯ 
        !
        ! with x-h < ξ < x+h. Since exp(x) is monotonically increasing, the 
        ! maximum value of exp(ξ) is exp(x+h). We verify now that the error is 
        ! always smaller than 
        !                           x+h
        !    ε  =  1.6666667e-05 * e
        !
        tollerance  = 1.666666667e-5
        h           = real(system_size_x / number_of_field_nodes_x) ! Step size
        do i = 2, number_of_field_nodes_x-1,10
            real_i = real(i)
            x = system_size_x * (real_i / real_number_of_field_nodes)
            call assert_equals(exp(x), density_gradient(1,i,1,1), tollerance*exp(x+h), "5  test_compute_density_gradient : Gradient of the 1D exponential density field configuration not calculated correctly")
        end do


        deallocate(density_field)
        deallocate(density_gradient)

    end subroutine test_compute_density_gradient

    subroutine test_allocate_field_arrays()
        number_of_field_nodes_x = 7
        number_of_field_nodes_y = 11
        number_of_field_nodes_z = 16
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_dimensions  = 3
        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        call assert_true(allocated(density_field), "1  test_allocate_field_arrays : density_field array not allocated")
        call assert_equals(number_of_field_nodes_x, size(density_field,1), "2  test_allocate_field_arrays : density_field x array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_y, size(density_field,2), "3  test_allocate_field_arrays : density_field y array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_z, size(density_field,3), "4  test_allocate_field_arrays : density_field z array not allocated to the correct size")

        call assert_true(allocated(density_gradient), "5  test_allocate_field_arrays : density_field array not allocated")
        call assert_equals(number_of_dimensions,  size(density_gradient,1), "6  test_allocate_field_arrays : density_gradient array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_x, size(density_gradient,2), "7  test_allocate_field_arrays : density_gradient x array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_y, size(density_gradient,3), "8  test_allocate_field_arrays : density_gradient y array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_z, size(density_gradient,4), "9  test_allocate_field_arrays : density_gradient z array not allocated to the correct size")
        
        call assert_true(allocated(position_of_density_nodes), "10 test_allocate_field_arrays : density_field array not allocated")
        call assert_equals(number_of_dimensions,  size(position_of_density_nodes,1), "11 test_allocate_field_arrays : position_of_density_nodes array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_x, size(position_of_density_nodes,2), "12 test_allocate_field_arrays : position_of_density_nodes x array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_y, size(position_of_density_nodes,3), "13 test_allocate_field_arrays : position_of_density_nodes y array not allocated to the correct size")
        call assert_equals(number_of_field_nodes_z, size(position_of_density_nodes,4), "14 test_allocate_field_arrays : position_of_density_nodes z array not allocated to the correct size")
    
        if (allocated(density_field)) then
            deallocate(density_field)
        end if
        if (allocated(density_gradient)) then
            deallocate(density_gradient)
        end if
        if (allocated(position_of_density_nodes)) then
            deallocate(position_of_density_nodes)
        end if
    end subroutine test_allocate_field_arrays

end module field_test