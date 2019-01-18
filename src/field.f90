module field
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use system,          only:  system_size
    use parameters,      only:  number_of_field_nodes_x,            &
                                number_of_field_nodes_y,            &
                                number_of_field_nodes_z,            &
                                number_of_dimensions,               &
                                number_of_particles
    implicit none
    private 

    ! The distance between density vertices.
    real (real64), private :: lx, ly, lz
    real (real64), private, dimension(:), allocatable :: l

    real (real64),   public, dimension(:,:,:),   allocatable :: density_field
    real (real64),   public, dimension(:,:,:,:), allocatable :: position_of_density_nodes
    real (real64),   public, dimension(:,:,:,:), allocatable :: density_gradient
    integer (int32), public, dimension(:),       allocatable :: number_of_field_nodes

    public  ::  compute_density_field,                      &
                allocate_field_arrays,                      &
                compute_density_gradient,                   &
                interpolate_density_field,                  &
                interpolate_density_gradient

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
        allocate(node_positions  (number_of_dimensions, number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))
        allocate(p               (number_of_dimensions))
        
        ! Distance between each vertex in x, y, and z directions.
        lx = system_size(1) / number_of_field_nodes_x
        ly = system_size(2) / number_of_field_nodes_y
        lz = system_size(3) / number_of_field_nodes_z

        if (.not. allocated(l)) then
            allocate(l(number_of_dimensions))
        end if
        l  = [lx, ly, lz]

        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]

        ! Compute the positions of the vertices.
        do i = 1, number_of_field_nodes_x
            do j = 1, number_of_field_nodes_y 
                do k = 1, number_of_field_nodes_z 
                    node_positions(:,i,j,k) = [i-1,j-1,k-1] * l
                    position_of_density_nodes(:,i,j,k) = node_positions(:,i,j,k)
                end do
            end do
        end do

        do i = 1, number_of_particles
            do j = 1, number_of_dimensions 
                node_vector(j) = floor(positions(j,i) / l(j))
            end do
            do j = 1, number_of_dimensions
                p(j) = positions(j,i) - node_vector(j) * system_size(j) / number_of_field_nodes(j)
            end do

            ! Computing the contribution to each nearest neighbor density vertex
            ! from this particle.
            vertex_contributions(0,0,0) = (lx - p(1)) * (ly - p(2)) * (lz - p(3)) / (lx * ly * lz)
            vertex_contributions(1,0,0) =    p(1)     * (ly - p(2)) * (lz - p(3)) / (lx * ly * lz)
            vertex_contributions(0,1,0) = (lx - p(1)) *    p(2)     * (lz - p(3)) / (lx * ly * lz)
            vertex_contributions(0,0,1) = (lx - p(1)) * (ly - p(2)) *    p(3)     / (lx * ly * lz)
            vertex_contributions(1,1,0) =    p(1)     *    p(2)     * (lz - p(3)) / (lx * ly * lz)
            vertex_contributions(1,0,1) =    p(1)     * (ly - p(2)) *    p(3)     / (lx * ly * lz)
            vertex_contributions(0,1,1) = (lx - p(1)) *    p(2)     *    p(3)     / (lx * ly * lz)
            vertex_contributions(1,1,1) =    p(1)     *    p(2)     *    p(3)     / (lx * ly * lz)

            ! Since node_vector indices start at 1, we need to add 1 here
            ! because of the floor function used. It is convenient to use the 
            ! 0:number_of_field_nodes-1 indexing until this point.
            do j = 1, number_of_dimensions
                node_vector(j) = node_vector(j) + 1

                ! Handle periodic boundary conditions for the density field.
                if (node_vector(j) == number_of_field_nodes(j)) then
                    node_vector_next(j) = 1
                else 
                    node_vector_next(j) = node_vector(j) + 1
                end if
            end do
            
            ! Scale the total contribution to the overall (number density) 
            ! field by the contributions from this particle. Note: Not mass 
            ! density.
            vertex_contributions = vertex_contributions ! * masses(i)
            
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

        ! Distance between each density vertex.
        lx = system_size(1) / number_of_field_nodes_x
        ly = system_size(2) / number_of_field_nodes_y
        lz = system_size(3) / number_of_field_nodes_z

        if (.not. allocated(l)) then
            allocate(l(number_of_dimensions))
        end if
        l  = [lx, ly, lz]

        density_gradient = 0

        do i = 1, number_of_field_nodes_x
            do j = 1, number_of_field_nodes_y
                do k = 1, number_of_field_nodes_z
                    i_next      = i + 1
                    i_previous  = i - 1

                    j_next      = j + 1
                    j_previous  = j - 1

                    k_next      = k + 1
                    k_previous  = k - 1

                    ! Periodic boundary conditions.
                    if (i == 1) then
                        i_previous = number_of_field_nodes_x
                    end if
                    if (i == number_of_field_nodes_x) then
                        i_next = 1
                    end if
                    if (j == 1) then
                        j_previous = number_of_field_nodes_y
                    end if
                    if (j == number_of_field_nodes_y) then
                        j_next = 1
                    end if
                    if (k == 1) then
                        k_previous = number_of_field_nodes_z
                    end if
                    if (k == number_of_field_nodes_z ) then
                        k_next = 1
                    end if

                    density_gradient(1,i,j,k) = (density_field(i_next, j,      k     ) - density_field(i_previous, j,          k         )) / (2.0_real64 * lx)
                    density_gradient(2,i,j,k) = (density_field(i,      j_next, k     ) - density_field(i,          j_previous, k         )) / (2.0_real64 * ly)
                    density_gradient(3,i,j,k) = (density_field(i,      j,      k_next) - density_field(i,          j,          k_previous)) / (2.0_real64 * lz)
                end do
            end do
        end do
    end subroutine compute_density_gradient_2

    subroutine compute_density_gradient_higher_order(density_field, order)
        real (real64), dimension(:,:,:), intent(in) :: density_field
        integer, intent(in) :: order
        real (real64) :: dummy

        !integer (int32) :: i, j, k
        !integer (int32) :: i_next,     j_next,     k_next
        !integer (int32) :: i_previous, j_previous, k_previous

        ! Assuming equal number of field vertices in each spatial direction for 
        ! now, so l = system_size(i) / number_of_field_nodes, i = 1,2,3 (x,y,z).
        lx = system_size(1) / number_of_field_nodes_x
        ly = system_size(2) / number_of_field_nodes_y
        lz = system_size(3) / number_of_field_nodes_z

        if (.not. allocated(l)) then
            allocate(l(number_of_dimensions))
        end if
        l  = [lx, ly, lz]

        number_of_field_nodes = [number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z]

        density_gradient = 0

        ! Supress compiler warnings about unused arguments: order, density_field
        dummy = order * density_field(1,1,1)

    end subroutine compute_density_gradient_higher_order

    real (real64) function interpolate_density_field(density_field,             &
                                                     position_of_density_nodes, &
                                                     point)                     &
                result(interpolated_density)
        implicit none
        real (real64), intent(in), dimension(:,:,:)   :: density_field
        real (real64), intent(in), dimension(:,:,:,:) :: position_of_density_nodes
        real (real64), intent(in), dimension(:)       :: point

        integer :: j, i, x_index, y_index, z_index
        integer, allocatable, dimension(:,:) :: offset
        integer, allocatable, dimension(:)   :: node_vector, opposite
        real (real64), allocatable, dimension(:)   :: partial_volume, cube
        real (real64) :: lattice_edge, point_edge, denominator, c, x0, x1

        ! There are 2ⁿ nearest neighbor lattice points in a system of n 
        ! dimensions, i.e. a hyper-cube in n dimensions has 2ⁿ corners.
        allocate(partial_volume(2**number_of_dimensions))
        allocate(node_vector(number_of_dimensions))
        allocate(cube(number_of_dimensions))
        allocate(offset(2**number_of_dimensions, number_of_dimensions))
        allocate(opposite(2**number_of_dimensions))

                              ! Opposite corners
        offset(1,:) = [0,0,0] ! 8 
        offset(2,:) = [1,0,0] ! 7
        offset(3,:) = [0,1,0] ! 6
        offset(4,:) = [0,0,1] ! 5 
        offset(5,:) = [1,1,0] ! 4
        offset(6,:) = [1,0,1] ! 3
        offset(7,:) = [0,1,1] ! 2
        offset(8,:) = [1,1,1] ! 1
        opposite = [8,7,6,5,4,3,2,1]

        ! Compute which field densities are closest, i.e. which (hyper-)cube we 
        ! are in.
        do j = 1, number_of_dimensions 
            node_vector(j) = floor(point(j) / l(j)) + 1
        end do
        
        ! Compute the partial volumes enclosed by the (hyper-)cubes with corners 
        ! in the interpolation point, and the field density vertices.
        do j = 1, 2**number_of_dimensions
            
            ! Compute the cube sides of the partial volume corresponding to 
            ! lattice vertex j.
            do i = 1, number_of_dimensions
                x_index = node_vector(1) + offset(j, 1)
                y_index = node_vector(2) + offset(j, 2)
                z_index = node_vector(3) + offset(j, 3)
                if (      (node_vector(1) == number_of_field_nodes_x)     &
                    .and. (offset(j,1)    == 1)) then
                    x_index = 1
                end if
                if (      (node_vector(2) == number_of_field_nodes_y)     &
                    .and. (offset(j,2)    == 1)) then
                    y_index = 1
                end if
                if (      (node_vector(3) == number_of_field_nodes_z)     &
                    .and. (offset(j,3)    == 1)) then
                    z_index = 1
                end if

                lattice_edge = position_of_density_nodes(i,               &
                                                         x_index,         & 
                                                         y_index,         &
                                                         z_index)
                point_edge   = point(i)
                cube(i) = abs(lattice_edge - point_edge)
            end do
            
            partial_volume(j) = cube(1)
            do i = 2, number_of_dimensions
                partial_volume(j) = partial_volume(j) * cube(i)
            end do
        end do
        
        interpolated_density = 0.0_real64
        do j = 1, 2**number_of_dimensions
            x_index = node_vector(1) + offset(j, 1)
            y_index = node_vector(2) + offset(j, 2)
            z_index = node_vector(3) + offset(j, 3)
            if (      (node_vector(1) == number_of_field_nodes_x)     &
                .and. (offset(j,1)    == 1)) then
                x_index = 1
            end if
            if (      (node_vector(2) == number_of_field_nodes_y)     &
                .and. (offset(j,2)    == 1)) then
                y_index = 1
            end if
            if (      (node_vector(3) == number_of_field_nodes_z)     &
                .and. (offset(j,3)    == 1)) then
                z_index = 1
            end if
            c = partial_volume(opposite(j)) * density_field(x_index,  & 
                                                            y_index,  &
                                                            z_index)
            interpolated_density = interpolated_density + c
        end do 
        denominator = 1.0_real64
        do i = 1, number_of_dimensions
            ! The offset(8,:) subarray contains [1,1,1].
            x_index = node_vector(1) + offset(8, 1)
            y_index = node_vector(2) + offset(8, 2)
            z_index = node_vector(3) + offset(8, 3)
            !print *, "node_vector, number_of_field_nodes_x : ", node_vector(3), number_of_field_nodes_x
            if (      (node_vector(1) == number_of_field_nodes_x)     &
                .and. (offset(8,1)    == 1)) then
                x_index = 1
            end if
            if (      (node_vector(2) == number_of_field_nodes_y)     &
                .and. (offset(8,2)    == 1)) then
                y_index = 1
            end if
            if (      (node_vector(3) == number_of_field_nodes_z)     &
                .and. (offset(8,3)    == 1)) then
                z_index = 1
            end if
            x0 = position_of_density_nodes(i,                  &
                                           node_vector(1),     &
                                           node_vector(2),     &
                                           node_vector(3))
            x1 = position_of_density_nodes(i,                  &
                                           x_index,            &
                                           y_index,            &
                                           z_index)
            denominator = denominator * (x1 - x0)
        end do
        interpolated_density = interpolated_density / denominator   
    end function interpolate_density_field

    subroutine interpolate_density_gradient(density_gradient,          &
                                            position_of_density_nodes, &
                                            point,                     &
                                            interpolated_gradient)
        implicit none
        real (real64), dimension(:,:,:,:), intent(in)     :: density_gradient
        real (real64), dimension(:,:,:,:), intent(in)     :: position_of_density_nodes
        real (real64), dimension(:),       intent(in)     :: point
        real (real64), dimension(:),       intent(in out) :: interpolated_gradient

        integer :: j, i, x_index, y_index, z_index
        integer, allocatable, dimension(:,:) :: offset
        integer, allocatable, dimension(:)   :: node_vector, opposite
        real (real64), allocatable, dimension(:)   :: partial_volume, cube
        real (real64) :: lattice_edge, point_edge, denominator, c, x0, x1

        ! There are 2ⁿ nearest neighbor lattice points in a system of n 
        ! dimensions, i.e. a hyper-cube in n dimensions has 2ⁿ corners.
        allocate(partial_volume(2**number_of_dimensions))
        allocate(node_vector(number_of_dimensions))
        allocate(cube(number_of_dimensions))
        allocate(offset(2**number_of_dimensions, number_of_dimensions))
        allocate(opposite(2**number_of_dimensions))

                              ! Opposite corners
        offset(1,:) = [0,0,0] ! 8 
        offset(2,:) = [1,0,0] ! 7
        offset(3,:) = [0,1,0] ! 6
        offset(4,:) = [0,0,1] ! 5 
        offset(5,:) = [1,1,0] ! 4
        offset(6,:) = [1,0,1] ! 3
        offset(7,:) = [0,1,1] ! 2
        offset(8,:) = [1,1,1] ! 1
        opposite = [8,7,6,5,4,3,2,1]

        ! Compute which field densities are closest, i.e. which (hyper-)cube we 
        ! are in.
        do j = 1, number_of_dimensions 
            node_vector(j) = int(floor(point(j) / l(j))) + 1
            print *, "l(j) : ", l(j) 
            print *, "point(j)  /  l(j) ", point(j)/l(j)
            print *, "node_vector: ", node_vector
        end do
        
        ! Compute the partial volumes enclosed by the (hyper-)cubes with corners 
        ! in the interpolation point, and the field density vertices.
        do j = 1, 2**number_of_dimensions
            
            ! Compute the cube sides of the partial volume corresponding to 
            ! lattice vertex j.
            do i = 1, number_of_dimensions
                x_index = node_vector(1) + offset(j, 1)
                y_index = node_vector(2) + offset(j, 2)
                z_index = node_vector(3) + offset(j, 3)
                if (      (node_vector(1) == number_of_field_nodes_x)     &
                    .and. (offset(j,1)    == 1)) then
                    x_index = 1
                end if
                if (      (node_vector(2) == number_of_field_nodes_y)     &
                    .and. (offset(j,2)    == 1)) then
                    y_index = 1
                end if
                if (      (node_vector(3) == number_of_field_nodes_z)     &
                    .and. (offset(j,3)    == 1)) then
                    z_index = 1
                end if
                lattice_edge = position_of_density_nodes(i,                             &
                                                         x_index, & 
                                                         y_index, &
                                                         z_index)
                point_edge   = point(i)
                cube(i) = abs(lattice_edge - point_edge)
            end do
            
            partial_volume(j) = cube(1)
            do i = 2, number_of_dimensions
                partial_volume(j) = partial_volume(j) * cube(i)
            end do
        end do
        
        interpolated_gradient = 0.0_real64
        do j = 1, 2**number_of_dimensions
            do i = 1, number_of_dimensions
                x_index = node_vector(1) + offset(j, 1)
                y_index = node_vector(2) + offset(j, 2)
                z_index = node_vector(3) + offset(j, 3)
                if (      (node_vector(1) == number_of_field_nodes_x)     &
                    .and. (offset(j,1)    == 1)) then
                    x_index = 1
                end if
                if (      (node_vector(2) == number_of_field_nodes_y)     &
                    .and. (offset(j,2)    == 1)) then
                    y_index = 1
                end if
                if (      (node_vector(3) == number_of_field_nodes_z)     &
                    .and. (offset(j,3)    == 1)) then
                    z_index = 1
                end if

                c = partial_volume(opposite(j)) * density_gradient(i,       &
                                                                   x_index, & 
                                                                   y_index, &
                                                                   z_index)
                interpolated_gradient(i) = interpolated_gradient(i) + c
            end do
        end do 
        denominator = 1.0_real64
        do i = 1, number_of_dimensions
            ! The offset(8,:) subarray contains [1,1,1].
            x_index = node_vector(1) + offset(8, 1)
            y_index = node_vector(2) + offset(8, 2)
            z_index = node_vector(3) + offset(8, 3)
            !print *, "node_vector, number_of_field_nodes_x : ", node_vector(3), number_of_field_nodes_x
            if (      (node_vector(1) == number_of_field_nodes_x)     &
                .and. (offset(8,1)    == 1)) then
                x_index = 1
            end if
            if (      (node_vector(2) == number_of_field_nodes_y)     &
                .and. (offset(8,2)    == 1)) then
                y_index = 1
            end if
            if (      (node_vector(3) == number_of_field_nodes_z)     &
                .and. (offset(8,3)    == 1)) then
                z_index = 1
            end if

            x0 = position_of_density_nodes(i,                  &
                                           node_vector(1),     &
                                           node_vector(2),     &
                                           node_vector(3))
            x1 = position_of_density_nodes(i,                  &
                                           x_index,            &
                                           y_index,            &
                                           z_index)
            denominator = denominator * (x1 - x0)
        end do
        interpolated_gradient = interpolated_gradient / denominator   

    end subroutine interpolate_density_gradient 

    subroutine allocate_field_arrays(density_field, density_gradient, position_of_density_nodes)
        implicit none
        real (real64), allocatable, dimension(:,:,:),   intent(in out) :: density_field
        real (real64), allocatable, dimension(:,:,:,:), intent(in out) :: position_of_density_nodes
        real (real64), allocatable, dimension(:,:,:,:), intent(in out) :: density_gradient

        allocate(density_field                                  (number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))
        allocate(position_of_density_nodes(number_of_dimensions, number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))
        allocate(density_gradient         (number_of_dimensions, number_of_field_nodes_x, number_of_field_nodes_y, number_of_field_nodes_z))
    end subroutine allocate_field_arrays

end module field