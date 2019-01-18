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
                            compute_density_gradient,           & ! Subroutine
                            interpolate_density_field,          & ! Subroutine
                            interpolate_density_gradient          ! Subroutine
    
    use particles,  only:   positions,                          &
                            masses
    implicit none
    private 

    public  ::  setup,                                                  &
                teardown,                                               &
                test_compute_density_field,                             &
                test_compute_density_field_periodic_boundaries,         &
                test_allocate_field_arrays,                             &
                test_compute_density_gradient,                          &
                test_interpolate_density_field,                         &
                test_interpolate_density_gradient   
contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_compute_density_field()
        integer (int32) :: i, j, k
        real (real64), allocatable, dimension(:,:,:) :: occam_density

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

        ! Check that the density field values correspond exactly to those found
        ! by OCCAM in a simple system of 10 particles and 2 vertices in each
        ! dimension. 
        number_of_dimensions = 3
        number_of_particles  = 10
        number_of_field_nodes_x = 2
        number_of_field_nodes_y = 2
        number_of_field_nodes_z = 2
        system_size_x         = 10.0_real64
        system_size_y         = 10.0_real64
        system_size_z         = 10.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses   (number_of_particles))
        allocate(occam_density(number_of_field_nodes_x,         &
                               number_of_field_nodes_y,         &
                               number_of_field_nodes_z))
        masses     = 1.0_real64
        
        positions(:,  1 ) = [0.05815856656600038_real64, 0.025342305369190687_real64, 0.19873935016981903_real64   ]
        positions(:,  2 ) = [0.07215982897620887_real64, 0.2504571863305104_real64,   0.38362945814487526_real64   ]
        positions(:,  3 ) = [0.49720339870241104_real64, 0.8872456338742205_real64,   0.9632329182983882_real64    ]
        positions(:,  4 ) = [0.5905505257877675_real64,  0.33492055184731184_real64,  0.8214094476851357_real64    ]
        positions(:,  5 ) = [0.5156196908805613_real64,  0.7799713140086558_real64,   0.007959114457485872_real64  ]
        positions(:,  6 ) = [0.8737540746260515_real64,  0.12371839911797189_real64,  0.6094189843776412_real64    ]
        positions(:,  7 ) = [0.7794755801413096_real64,  0.44018012706879395_real64,  0.10548340557947777_real64   ]
        positions(:,  8 ) = [0.9931810993676221_real64,  0.26069746565869656_real64,  0.7264137166754773_real64    ]
        positions(:,  9 ) = [0.2898648348703551_real64,  0.5064354936235451_real64,   0.8455429407939724_real64    ]
        positions(:, 10 ) = [0.38696815531713213_real64, 0.898330488912219_real64,    0.3928309115958717_real64    ]
        
        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        call compute_density_field(positions, masses)

        occam_density(1,1,1) = 7.360439175029531_real64
        occam_density(2,1,1) = 0.820343600335064_real64
        occam_density(1,2,1) = 0.726866505217675_real64
        occam_density(2,2,1) = 0.081418669862102_real64
        occam_density(1,1,2) = 0.817321740122457_real64
        occam_density(2,1,2) = 0.100435691350726_real64
        occam_density(1,2,2) = 0.083985428583254_real64
        occam_density(2,2,2) = 0.009189189499192_real64

        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    call assert_equals(occam_density(i,j,k), density_field(i,j,k), 1.0e-15_real64, "80 test_compute_density_field : Calculated density field values are not equal to the corresponding values calculated in OCCAM")
                end do
            end do
        end do

        deallocate(positions)
        deallocate(masses)
        deallocate(occam_density)
        deallocate(density_field)
        deallocate(density_gradient)
        deallocate(position_of_density_nodes)


    end subroutine test_compute_density_field

    subroutine test_compute_density_gradient()
        integer (int32) :: i, j, N
        real (real64)   :: pi, x, real_i, real_number_of_field_nodes
        real (real64)   :: tollerance, h
        real (real64), allocatable, dimension(:) :: convergence_rate,       &
                                                    step_lengths,           &
                                                    errors
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
        any_nonzero = any(abs(density_gradient) > 1.0e-14)
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
        pi = 3.141592653589793238462643383279502884197169399375105820974_real64
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

        ! Let us redo the exponential test, while doubling the number of nodes
        ! and try to figure out the error scaling w.r.t. the step length, h.
        ! 
        ! See: 
        !
        ! http://hplgit.github.io/INF5620/doc/notes/lecture_decay.html#___sec17
        !
        system_size_x         = 1.0_real64
        system_size_y         = 1.0_real64
        system_size_z         = 1.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]
        number_of_field_nodes_x = 5
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        
        N = 12
        allocate(convergence_rate(N))
        allocate(step_lengths(N))
        allocate(errors(N))
        convergence_rate = 0
        step_lengths     = 0
        errors           = 0
        
        do j = 1, N
            ! Double the number of nodes each iteration, halving the step 
            ! length.
            if (j /= 1) then
                number_of_field_nodes_x     = number_of_field_nodes_x * 2
                number_of_field_nodes       = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
                real_number_of_field_nodes  = real(number_of_field_nodes_x)
            end if
            step_lengths(j) = system_size_x / real_number_of_field_nodes
                        
            deallocate(density_field)
            deallocate(density_gradient)
            allocate(density_field                         (number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))
            allocate(density_gradient(number_of_dimensions, number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))            

            do i = 1, number_of_field_nodes_x
                real_i = real(i)
                x = system_size_x * (real_i / real_number_of_field_nodes)
                density_field(i,:,:) = exp(x)
            end do
            
            call compute_density_gradient(density_field)

            ! Since density_field(:,1,1) holds the values of exp(x), we just use
            ! it as the known exact value that density_gradient *should* have.
            errors(j) = maxval(abs(density_gradient(1,2:number_of_field_nodes_x-1,1,1) - density_field(2:number_of_field_nodes_x-1,1,1)))
        end do

        do j = 2, N
            convergence_rate(j) = log(errors(j-1) / errors(j)) / log(step_lengths(j-1) / step_lengths(j))
        end do
        call assert_equals(2.0_real64, convergence_rate(N), 0.001_real64, "6  test_compute_density_gradient : The computed convergence rate of the central finite difference scheme used in compute_density_gradient is not sufficiently close to 2.0")
    
        
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

    subroutine test_interpolate_density_field()
        implicit none
        integer (int32) :: i, j, k
        real (real64), allocatable, dimension(:,:) :: points
        real (real64), allocatable, dimension(:)   :: point, values
        real (real64) :: interpolated_density

        ! We start with a single cube (8 density lattice points) with the 
        ! density at (0,0,0) being 1 and vanishing at all other lattice points.
        ! The interpolated density at points on the face of the cube *should* be
        ! 1-d, with d being the distance between the point and the (0,0,0) 
        ! lattice point.
        number_of_field_nodes_x = 2
        number_of_field_nodes_y = 2
        number_of_field_nodes_z = 2
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 2.0_real64
        system_size_y         = 2.0_real64
        system_size_z         = 2.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        if (allocated(positions)) then
            deallocate(positions)
        end if
        if (allocated(masses)) then
            deallocate(masses) 
        end if

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        density_field = 0.0_real64
        allocate(point(number_of_dimensions))
        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses   (number_of_particles))
        masses         = 1.0
        positions(:,1) = [0.0, 0.0, 0.0]
        
        ! This call sets up position_of_density_nodes, so it is neccessary even
        ! though we throw away the computed density values.
        call compute_density_field(positions, masses)
        density_field = 0.0_real64
        density_field(1,1,1) = 1.0_real64
        point = [0.5, 0.5, 0.5]
        interpolated_density = interpolate_density_field(density_field, position_of_density_nodes, point)
        
        ! A single field vertex with value 1, and all others vanishing, should
        ! result in an interpolated value of 1 / 8 evaluated in the middle of 
        ! the cube.
        call assert_equals(1.0_real64 / 8.0_real64, interpolated_density, "11 test_interpolate_density_field : interpolated density in the middle of the cube was not calculated correctly")
        
        ! Repeat with all the other single vertices.
        do i = 1, 2
            do j = 1, 2
                do k = 1, 2
                    density_field = 0.0_real64
                    density_field(i,j,k) = 1.0_real64
                    interpolated_density = interpolate_density_field(density_field, position_of_density_nodes, point)
                    call assert_equals(1.0_real64 / 8.0_real64, interpolated_density, "12 test_interpolate_density_field : interpolated density in the middle of the cube was not calculated correctly")
                end do
            end do
        end do

        ! Next we test on some random data generated and interpolated using 
        ! the scipy.interpolate package in python3 (script in /tools folder).
        density_field( 1 , 1 , 1 ) = 0.47919988065709096_real64
        density_field( 1 , 1 , 2 ) = 0.4432262866925848_real64
        density_field( 1 , 2 , 1 ) = 0.14775650101900195_real64
        density_field( 1 , 2 , 2 ) = 0.3311914156552239_real64
        density_field( 2 , 1 , 1 ) = 0.544945955559647_real64
        density_field( 2 , 1 , 2 ) = 0.33486182332430936_real64
        density_field( 2 , 2 , 1 ) = 0.28336184861496694_real64
        density_field( 2 , 2 , 2 ) = 0.7670556317983035_real64

        allocate(points(number_of_dimensions, 10))
        points(:, 1 )  = [ 0.3226243761072245_real64,  0.16000928172147622_real64, 0.07474857607741736_real64  ]
        points(:, 2 )  = [ 0.5092014137865475_real64,  0.1184896077950004_real64,  0.25241336618190535_real64  ]
        points(:, 3 )  = [ 0.35834120119087876_real64, 0.2809403053449454_real64,  0.0003234148574184914_real64]
        points(:, 4 )  = [ 0.23480781215213598_real64, 0.7843375272202351_real64,  0.9271719737772645_real64   ]
        points(:, 5 )  = [ 0.5604832314057514_real64,  0.5011060615131893_real64,  0.36617926974252146_real64  ]
        points(:, 6 )  = [ 0.9363242421592629_real64,  0.7424379445058036_real64,  0.9412666162250206_real64   ]
        points(:, 7 )  = [ 0.8282735698757658_real64,  0.9659206984116183_real64,  0.2605764978122821_real64   ]
        points(:, 8 )  = [ 0.9899606876712709_real64,  0.7170472284509276_real64,  0.22977023945233688_real64  ]
        points(:, 9 )  = [ 0.8239128542674138_real64,  0.7495958518805627_real64,  0.18555519519490316_real64  ]
        points(:, 10 ) = [ 0.43629733663334413_real64, 0.85447599895436_real64,    0.4369060699031969_real64   ]

        allocate(values(10))
        values( 1 ) = 0.448550408370945_real64
        values( 2 ) = 0.4599482323769342_real64
        values( 3 ) = 0.41668007845152427_real64
        values( 4 ) = 0.41683986940539497_real64
        values( 5 ) = 0.40972231149742805_real64
        values( 6 ) = 0.6196626661445901_real64
        values( 7 ) = 0.37656402561364555_real64
        values( 8 ) = 0.421864931428537_real64
        values( 9 ) = 0.3796523928562688_real64
        values( 10 ) = 0.3609886588271966_real64

        do i = 1, 10
            interpolated_density = interpolate_density_field(density_field, position_of_density_nodes, points(:,i))
            call assert_equals(values(i), interpolated_density, 1e-14_real64, "21 test_interpolate_density_field : interpolated density does not equal the scipy.interpolate values for the tested reference points / field values")
        end do

        ! Check that the interpolated density field values correspond exactly to
        ! those found by OCCAM in a simple system of 10 particles and 2 
        ! vertices in each dimension. 
        number_of_field_nodes_x = 2
        number_of_field_nodes_y = 2
        number_of_field_nodes_z = 2
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 10
        number_of_dimensions  = 3
        system_size_x         = 10.0_real64
        system_size_y         = 10.0_real64
        system_size_z         = 10.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        ! Reallocate positions and masses to enable calling 
        ! compute_density_field. This is only done as a lazy way to compute 
        ! position_of_density_nodes--the density values are disregarded.
        deallocate(positions)
        deallocate(masses)
        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses(number_of_particles))
        positions = 0.0_real64
        masses = 0.0_real64
        call compute_density_field(positions, masses)

        density_field(1,1,1) = 7.360439175029531_real64
        density_field(2,1,1) = 0.820343600335064_real64
        density_field(1,2,1) = 0.726866505217675_real64
        density_field(2,2,1) = 0.081418669862102_real64
        density_field(1,1,2) = 0.817321740122457_real64
        density_field(2,1,2) = 0.100435691350726_real64
        density_field(1,2,2) = 0.083985428583254_real64
        density_field(2,2,2) = 0.009189189499192_real64

        points(:,  1 ) = [0.05815856656600038_real64, 0.025342305369190687_real64, 0.19873935016981903_real64   ]
        points(:,  2 ) = [0.07215982897620887_real64, 0.2504571863305104_real64,   0.38362945814487526_real64   ]
        points(:,  3 ) = [0.49720339870241104_real64, 0.8872456338742205_real64,   0.9632329182983882_real64    ]
        points(:,  4 ) = [0.5905505257877675_real64,  0.33492055184731184_real64,  0.8214094476851357_real64    ]
        points(:,  5 ) = [0.5156196908805613_real64,  0.7799713140086558_real64,   0.007959114457485872_real64  ]
        points(:,  6 ) = [0.8737540746260515_real64,  0.12371839911797189_real64,  0.6094189843776412_real64    ]
        points(:,  7 ) = [0.7794755801413096_real64,  0.44018012706879395_real64,  0.10548340557947777_real64   ]
        points(:,  8 ) = [0.9931810993676221_real64,  0.26069746565869656_real64,  0.7264137166754773_real64    ]
        points(:,  9 ) = [0.2898648348703551_real64,  0.5064354936235451_real64,   0.8455429407939724_real64    ]
        points(:, 10 ) = [0.38696815531713213_real64, 0.898330488912219_real64,    0.3928309115958717_real64    ]
    
        values(1 ) = 6.9948858249501731_real64
        values(2 ) = 6.4648355463455784_real64
        values(3 ) = 4.6718703421534755_real64
        values(4 ) = 5.2864544226407126_real64
        values(5 ) = 5.7378962066850665_real64
        values(6 ) = 5.4204577534819203_real64
        values(7 ) = 5.7283284816729338_real64
        values(8 ) = 5.0307551561218400_real64
        values(9 ) = 5.3904463286467665_real64
        values(10) = 5.3433024738229600_real64

        do i = 1, 10
            interpolated_density = interpolate_density_field(density_field, position_of_density_nodes, points(:,i))
            call assert_equals(values(i), interpolated_density, 1e-15_real64, "30 test_interpolate_density_field : Interpolated density values are not equal to the corresponding interpolated density values calculated in OCCAM")
        end do

        deallocate(points)
        deallocate(values)
        deallocate(positions)
        deallocate(masses)
        deallocate(density_field)
        deallocate(density_gradient)
        deallocate(position_of_density_nodes)
    end subroutine test_interpolate_density_field

    subroutine test_interpolate_density_gradient()
        implicit none
        integer (int32) :: i
        real (real64), allocatable, dimension(:,:) :: points, values
        real (real64), allocatable, dimension(:)   :: interpolated_gradient

        number_of_field_nodes_x = 2
        number_of_field_nodes_y = 2
        number_of_field_nodes_z = 2
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 2.0_real64
        system_size_y         = 2.0_real64
        system_size_z         = 2.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        if (allocated(positions)) then
            deallocate(positions)
        end if
        if (allocated(masses)) then
            deallocate(masses) 
        end if

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        density_field = 0.0_real64
        allocate(interpolated_gradient(number_of_dimensions))
        allocate(points   (number_of_dimensions,10))
        allocate(values   (number_of_dimensions, 10))
        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses   (number_of_particles))
        masses         = 1.0
        positions(:,1) = [0.0, 0.0, 0.0]
        
        ! This call sets up position_of_density_nodes, so it is neccessary even
        ! though we throw away the computed density values.
        call compute_density_field   (positions, masses)
        call compute_density_gradient(density_field)
        density_field    = 0.0_real64
        density_gradient = 0.0_real64
        
        density_gradient(:, 1 , 1 , 1 ) = [ 0.5540271911850455_real64,  0.3961947008639597_real64, 0.3709350754350955_real64  ]
        density_gradient(:, 1 , 1 , 2 ) = [ 0.29989066985724067_real64, 0.5256264968923157_real64, 0.9835075358585761_real64  ]
        density_gradient(:, 1 , 2 , 1 ) = [ 0.5977836288708848_real64,  0.5695607469442536_real64, 0.5911496386331424_real64  ]
        density_gradient(:, 1 , 2 , 2 ) = [ 0.08186290845694477_real64, 0.6297744649628887_real64, 0.2847708818027924_real64  ]
        density_gradient(:, 2 , 1 , 1 ) = [ 0.4147543677628468_real64,  0.7465105563078386_real64, 0.3990646219350965_real64  ]
        density_gradient(:, 2 , 1 , 2 ) = [ 0.6606266575728618_real64,  0.53231087874991_real64,   0.08359344200872987_real64 ]
        density_gradient(:, 2 , 2 , 1 ) = [ 0.6443874169895922_real64,  0.2935243867752405_real64, 0.2603839953354746_real64  ]
        density_gradient(:, 2 , 2 , 2 ) = [ 0.42547562683826445_real64, 0.5521127426044589_real64, 0.8916053404339211_real64  ]

        points(:, 1 )  = [ 0.5775948748056732_real64,   0.18410288727491742_real64, 0.011502436296286556_real64  ]
        points(:, 2 )  = [ 0.16892202728673078_real64,  0.29878402069780086_real64, 0.1594028147281893_real64    ]
        points(:, 3 )  = [ 0.6987832468901581_real64,   0.7371310391437007_real64,  0.0015824144694810416_real64 ]
        points(:, 4 )  = [ 0.9530700783586349_real64,   0.808497523595109_real64,   0.03912110333739105_real64   ]
        points(:, 5 )  = [ 0.924003466783018_real64,    0.7564474589497775_real64,  0.3056221180804195_real64    ]
        points(:, 6 )  = [ 0.3445125008638579_real64,   0.10130670707437384_real64, 0.363201051304211_real64     ]
        points(:, 7 )  = [ 0.5007985545190358_real64,   0.9880660518539662_real64,  0.06302118081541042_real64   ]
        points(:, 8 )  = [ 0.5123913896979235_real64,   0.09895123684176588_real64, 0.847060755818077_real64     ]
        points(:, 9 )  = [ 0.016671988827459905_real64, 0.567731274356605_real64,   0.068584064535122_real64     ]
        points(:, 10 ) = [ 0.17486109033920194_real64,  0.9450808716033482_real64,  0.7137522199688752_real64    ]

        values(:, 1 )  = [ 0.5010012334846653_real64,  0.5635703374727197_real64,  0.390777216015084_real64   ]
        values(:, 2 )  = [ 0.5118083828970728_real64,  0.48799989474168715_real64, 0.46726883870339936_real64 ]
        values(:, 3 )  = [ 0.5843839053101522_real64,  0.4463381487722325_real64,  0.36844537952444056_real64 ]
        values(:, 4 )  = [ 0.5941968437894098_real64,  0.3939981659320265_real64,  0.3157752312771956_real64  ]
        values(:, 5 )  = [ 0.5480036081000176_real64,  0.4555317324527523_real64,  0.42387899302020093_real64 ]
        values(:, 6 )  = [ 0.47502269821515997_real64, 0.520921272713153_real64,   0.4866094799359585_real64  ]
        values(:, 7 )  = [ 0.5966229450301203_real64,  0.44290013053647614_real64, 0.4352912067237929_real64  ]
        values(:, 8 )  = [ 0.46750660652443143_real64, 0.539058626463431_real64,   0.5081330256110408_real64  ]
        values(:, 9 )  = [ 0.5511235516756882_real64,  0.5006726571356255_real64,  0.499410890950402_real64   ]
        values(:, 10 ) = [ 0.2822271319148202_real64,  0.5845211015987357_real64,  0.44623908382011673_real64 ]

        do i = 1, 10
            call interpolate_density_gradient(density_gradient,                &
                                              position_of_density_nodes,       &
                                              points(:,i),                     &
                                              interpolated_gradient)
            call assert_equals(values(:,i), interpolated_gradient, 3, 1e-14_real64, "11 test_interpolate_density_gradient : interpolated density gradient does not equal the scipy.interpolate values for the tested reference points / field values")
        end do

        deallocate(positions)
        deallocate(masses)
        deallocate(density_field)
        deallocate(density_gradient)
        deallocate(position_of_density_nodes)
    end subroutine test_interpolate_density_gradient
    
    subroutine test_compute_density_field_periodic_boundaries() 
        implicit none
        integer (int32) :: i, j, k

        number_of_field_nodes_x = 3
        number_of_field_nodes_y = 3
        number_of_field_nodes_z = 3
        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]
        number_of_particles   = 1
        number_of_dimensions  = 3
        system_size_x         = 1.0_real64
        system_size_y         = 1.0_real64
        system_size_z         = 1.0_real64
        system_size           = [system_size_x, system_size_y, system_size_z]

        allocate(positions(number_of_dimensions, number_of_particles))
        allocate(masses(number_of_particles))

        call allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        positions = 0.0_real64
        masses = 1.0_real64
        
        ! Particle in the center of the (nx,ny,nz) cell should contribute its 
        ! number density equally between the first and the last density 
        ! vertices, i.e. density_field(i,j,k) for i,j,k=1 and 
        ! i,j,k=number_of_field_nodes. Nx, ny, nz is short for 
        ! number_of_field_nodes_x, etc. 
        positions(:,1) = [0.83333333333333333333_real64,    &
                          0.83333333333333333333_real64,    &
                          0.83333333333333333333_real64]
        call compute_density_field(positions, masses)
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,1), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,1), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,1), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,3), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,1), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,3), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,3), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,3), 1.0e-15_real64, "11 test_compute_density_field_periodic_boundaries : computed number density is not distributed correctly w.r.t. the periodic boundaries")
        deallocate(density_field)
        deallocate(density_gradient)
        deallocate(position_of_density_nodes)
        deallocate(positions)
        deallocate(masses)

    end subroutine test_compute_density_field_periodic_boundaries

end module field_test