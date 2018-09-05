module field
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,          only:  system_size
    use parameters,      only:  number_of_field_nodes,              &
                                number_of_dimensions,               &
                                number_of_particles
    implicit none
    private 

    real (real64), public, dimension(:,:,:),   allocatable :: density_field
    real (real64), public, dimension(:,:,:),   allocatable :: position_of_density_nodes
    real (real64), public, dimension(:,:,:,:), allocatable :: density_gradient

    public ::   compute_density_field,              &
                allocate_field_arrays

contains

    subroutine compute_density_field(positions)
        implicit none
        real (real64), dimension(:,:), intent(in) :: positions
        
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

        ! Internal position in the cube of nearest density vertex points.
        real (real64),   dimension(:), allocatable :: p 

        real (real64)   :: l
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

        print *,"=====----=====----====="
        print *,"nodes_pos:  ", node_positions(1,:,1,1)
        print *,"=====----=====----====="

        do i = 1, number_of_particles
            do j = 1, number_of_dimensions 
                node_vector(j) = floor(positions(j,i) / l)
            end do
            do j = 1, number_of_dimensions
                p(j) = positions(j,i) - node_vector(j) * system_size(j) / number_of_field_nodes
            end do
            print *, " ====== "
            print *, "pos/l:       ", positions(1,i) / l
            print *, "l:           ", l
            print *, "pos:         ", positions(:,i)
            print *, "nodes:       ", (node_vector) * l
            print *, "node_vector: ", (node_vector)
            print *, "p:           ", p
            print *, " ====== "
            ! Sum of all the distances |p - vertex(i,j,k)| for all the 8 nearest
            ! vertices.
            ! sum_of_vertex_distances = sqrt(5.0 * ( p(1)**2 + p(2)**2 + p(3)**2 ) + 3 * l**2 - 2 * l * (p(1) + p(2) + p(3)))
            
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
            
            ! print *, vertex_contributions(0,0,0) + vertex_contributions(1,0,0) + vertex_contributions(0,1,0) + vertex_contributions(0,0,1) + vertex_contributions(1,1,0) + vertex_contributions(1,0,1) + vertex_contributions(0,1,1) + vertex_contributions(1,1,1)
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

    subroutine allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        implicit none
        real (real64), allocatable, dimension(:,:,:),   intent(in out) :: density_field
        real (real64), allocatable, dimension(:,:,:),   intent(in out) :: position_of_density_nodes
        real (real64), allocatable, dimension(:,:,:,:), intent(in out) :: density_gradient

        allocate(density_field            (number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
        allocate(position_of_density_nodes(number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
        allocate(density_gradient         (number_of_dimensions,  number_of_field_nodes, number_of_field_nodes, number_of_field_nodes))
    end subroutine allocate_field_arrays

end module field