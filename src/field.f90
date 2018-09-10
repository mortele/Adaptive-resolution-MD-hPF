module field
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,          only:  system_size
    use parameters,      only:  number_of_field_nodes,              &
                                number_of_dimensions,               &
                                number_of_particles
    implicit none
    private 

    ! The distance between density vertices.
    real (real64), private :: l

    real (real64), public, dimension(:,:,:),   allocatable :: density_field
    real (real64), public, dimension(:,:,:,:), allocatable :: position_of_density_nodes
    real (real64), public, dimension(:,:,:,:), allocatable :: density_gradient


    public  ::  compute_density_field,                      &
                allocate_field_arrays,                      &
                compute_density_gradient

    private ::  compute_density_gradient_2,                 &
                compute_density_gradient_higher_order

    interface compute_density_gradient
        module procedure compute_density_gradient_2
        module procedure compute_density_gradient_higher_order
    end interface

contains

    subroutine compute_density_field(positions, masses)
        implicit none
        real (real64), dimension(:,:), intent(in) :: positions
        real (real64), dimension(:),   intent(in) :: masses
        
        ! Tensor for holding the contribution to each of the density nodes 
        ! closest to the particle in question. There are 8 such points, 
        ! enclosing each particle in a cube of side lengths l.
        real (real64),   dimension(0:1,0:1,0:1) :: vertex_contributions

        ! Indices in the density field tensor corresponding to the vertices 
        ! closest to the particle with positions x_vertex < x, y_vertex < y, and 
        ! z_vertex < z.
        integer (int32), dimension(:), allocatable :: node_vector

        ! The cooresponding indices for the closest vertices with x_vertex > x,
        ! y_vertex > y, and z_vertex > z. Accounts for the periodic boundary 
        ! conditions, so if node_vector(i) equals number_of_field_nodes, then
        ! node_vector_next(i) will be 1 and not number_of_field_nodes+1.
        integer (int32), dimension(:), allocatable :: node_vector_next

        ! Positions of the density field vertices on which the field is 
        ! computed.
        real (real64),   dimension(:,:,:,:), allocatable :: node_positions

        ! Internal position inside the cube of 8 nearest density vertex points.
        real (real64),   dimension(:), allocatable :: p 

        integer (int32) :: i, j, k

        density_field = 0

        allocate(node_vector     (number_of_dimensions))
        allocate(node_vector_next(number_of_dimensions))
        allocate(node_positions  (number_of_dimensions, number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
        allocate(p               (number_of_dimensions))
        
        ! Assuming equal number of field vertices in each spatial direction for 
        ! now, so l = system_size(i) / number_of_field_nodes, i = 1,2,3 (x,y,z).
        l = system_size(1) / number_of_field_nodes

        ! Compute the positions of the vertices.
        do i = 1, number_of_field_nodes 
            do j = 1, number_of_field_nodes 
                do k = 1, number_of_field_nodes 
                    node_positions(:,i,j,k) = [i-1,j-1,k-1] * l
                end do
            end do
        end do

        do i = 1, number_of_particles
            do j = 1, number_of_dimensions 
                node_vector(j) = floor(positions(j,i) / l)
            end do
            do j = 1, number_of_dimensions
                p(j) = positions(j,i) - node_vector(j) * system_size(j) / number_of_field_nodes
            end do

            ! Computing the contribution to each nearest neighbor density vertex
            ! from this particle.
            vertex_contributions(0,0,0) = (l - p(1)) * (l - p(2)) * (l - p(3)) / l**3
            vertex_contributions(1,0,0) =    p(1)    * (l - p(2)) * (l - p(3)) / l**3
            vertex_contributions(0,1,0) = (l - p(1)) *    p(2)    * (l - p(3)) / l**3
            vertex_contributions(0,0,1) = (l - p(1)) * (l - p(2)) *    p(3)    / l**3            
            vertex_contributions(1,1,0) =    p(1)    *    p(2)    * (l - p(3)) / l**3
            vertex_contributions(1,0,1) =    p(1)    * (l - p(2)) *    p(3)    / l**3
            vertex_contributions(0,1,1) = (l - p(1)) *    p(2)    *    p(3)    / l**3
            vertex_contributions(1,1,1) =    p(1)    *    p(2)    *    p(3)    / l**3

            ! Since node_vector indices start at 1, we need to add 1 here
            ! because of the floor function used. It is convenient to use the 
            ! 0:number_of_field_nodes-1 indexing until this point.
            do j = 1, number_of_dimensions
                node_vector(j) = node_vector(j) + 1

                ! Handle periodic boundary conditions for the density field.
                if (node_vector(j) == number_of_field_nodes) then
                    node_vector_next(j) = 1
                else 
                    node_vector_next(j) = node_vector(j) + 1
                end if
            end do
            
            ! Scale the total contribution to the overall density field by the 
            ! mass of the particle in question.
            vertex_contributions = vertex_contributions * masses(i)
            
            ! The tensor vertex_contributions holds the contribution to the 
            ! density field at the 8 nodes closest to the particle in question.
            ! We now use the node_vector to map this onto the global indices of 
            ! the density_field array and add the relevant contributions.
            density_field(node_vector(1),      node_vector(2),      node_vector(3))      = density_field(node_vector(1),      node_vector(2),      node_vector(3))      + vertex_contributions(0,0,0)
            density_field(node_vector_next(1), node_vector(2),      node_vector(3))      = density_field(node_vector_next(1), node_vector(2),      node_vector(3))      + vertex_contributions(1,0,0)
            density_field(node_vector(1),      node_vector_next(2), node_vector(3))      = density_field(node_vector(1),      node_vector_next(2), node_vector(3))      + vertex_contributions(0,1,0)
            density_field(node_vector(1),      node_vector(2),      node_vector_next(3)) = density_field(node_vector(1),      node_vector(2),      node_vector_next(3)) + vertex_contributions(0,0,1)
            density_field(node_vector_next(1), node_vector_next(2), node_vector(3))      = density_field(node_vector_next(1), node_vector_next(2), node_vector(3))      + vertex_contributions(1,1,0)
            density_field(node_vector_next(1), node_vector(2),      node_vector_next(3)) = density_field(node_vector_next(1), node_vector(2),      node_vector_next(3)) + vertex_contributions(1,0,1)
            density_field(node_vector(1),      node_vector_next(2), node_vector_next(3)) = density_field(node_vector(1),      node_vector_next(2), node_vector_next(3)) + vertex_contributions(0,1,1)
            density_field(node_vector_next(1), node_vector_next(2), node_vector_next(3)) = density_field(node_vector_next(1), node_vector_next(2), node_vector_next(3)) + vertex_contributions(1,1,1)
        end do
    end subroutine compute_density_field

    subroutine compute_density_gradient_2(density_field)
        real (real64), dimension(:,:,:), intent(in) :: density_field
        integer (int32) :: i, j, k
        integer (int32) :: i_next,     j_next,     k_next
        integer (int32) :: i_previous, j_previous, k_previous

        ! Assuming equal number of field vertices in each spatial direction for 
        ! now, so l = system_size(i) / number_of_field_nodes, i = 1,2,3 (x,y,z).
        l = system_size(1) / number_of_field_nodes

        density_gradient = 0

        do i = 1, number_of_field_nodes
            do j = 1, number_of_field_nodes
                do k = 1, number_of_field_nodes
                    i_next      = i + 1
                    i_previous  = i - 1

                    j_next      = j + 1
                    j_previous  = j - 1

                    k_next      = k + 1
                    k_previous  = k - 1

                    ! Periodic boundary conditions.
                    if (i == 1) then
                        i_previous = number_of_field_nodes
                    end if
                    if (i == number_of_field_nodes) then
                        i_next = 1
                    end if
                    if (j == 1) then
                        j_previous = number_of_field_nodes
                    end if
                    if (j == number_of_field_nodes) then
                        j_next = 1
                    end if
                    if (k == 1) then
                        k_previous = number_of_field_nodes
                    end if
                    if (k == number_of_field_nodes) then
                        k_next = 1
                    end if

                    density_gradient(1,i,j,k) = (density_field(i_next, j,      k     ) - density_field(i_previous, j,          k         )) / (2.0_real64 * l)
                    density_gradient(2,i,j,k) = (density_field(i,      j_next, k     ) - density_field(i,          j_previous, k         )) / (2.0_real64 * l)
                    density_gradient(3,i,j,k) = (density_field(i,      j,      k_next) - density_field(i,          j,          k_previous)) / (2.0_real64 * l)
                end do
            end do
        end do
    end subroutine compute_density_gradient_2

    subroutine compute_density_gradient_higher_order(density_field, order)
        real (real64), dimension(:,:,:), intent(in) :: density_field
        integer, intent(in) :: order

        !integer (int32) :: i, j, k
        !integer (int32) :: i_next,     j_next,     k_next
        !integer (int32) :: i_previous, j_previous, k_previous

        ! Assuming equal number of field vertices in each spatial direction for 
        ! now, so l = system_size(i) / number_of_field_nodes, i = 1,2,3 (x,y,z).
        l = system_size(1) / number_of_field_nodes

        density_gradient = 0

    end subroutine compute_density_gradient_higher_order

    subroutine allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        implicit none
        real (real64), allocatable, dimension(:,:,:),   intent(in out) :: density_field
        real (real64), allocatable, dimension(:,:,:,:), intent(in out) :: position_of_density_nodes
        real (real64), allocatable, dimension(:,:,:,:), intent(in out) :: density_gradient

        allocate(density_field            (number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
        allocate(position_of_density_nodes(number_of_dimensions,  number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
        allocate(density_gradient         (number_of_dimensions,  number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
    end subroutine allocate_field_arrays

end module field