module field_test
    use, intrinsic :: iso_fortran_env,  only: int32, real64
    use fruit,      only:   assert_equals,                       &
                            assert_true
    use system,     only:   system_size
    use parameters, only:   number_of_field_nodes,              &
                            number_of_dimensions,               &
                            number_of_particles,                &
                            system_size_x,                      &
                            system_size_y,                      &
                            system_size_z
    use field,      only:   allocate_field_arrays,              & ! Submodule
                            compute_density_field,              & ! Submodule
                            density_field,                      &
                            position_of_density_nodes,          &
                            density_gradient
    
    use particles,  only:   positions,                          &
                            masses
    implicit none
    private 

    public  ::  setup,                          &
                teardown,                       &
                test_compute_density_field,     &
                test_allocate_field_arrays
contains

    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_compute_density_field()
        integer (int32) :: i, j, k

        ! We start off with a 1D test of the density field computation, with a 
        ! single particle.
        number_of_field_nodes = 3
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
        number_of_field_nodes = 3
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
        number_of_field_nodes = 3
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
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,1), "60 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,1), "61 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,1), "62 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,1,3), "63 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,1), "64 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,1,3), "65 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(1,3,3), "66 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")
        call assert_equals(1.0_real64 / 8.0_real64, density_field(3,3,3), "67 test_compute_density_field : A particle in 3D placed in the center of a cube of 8 density vertices should evenly distribute 12.5% of its mass onto each of these vertices, periodic boundary version (particle on the right hand side edge of the simulation box)")


    end subroutine test_compute_density_field

    subroutine test_allocate_field_arrays()

    end subroutine test_allocate_field_arrays

end module field_test