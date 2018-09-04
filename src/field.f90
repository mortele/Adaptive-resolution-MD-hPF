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
        
        integer (int32), dimension(:), allocatable :: node_vector

        ! Internal position in the cube of nearest density vertex points.
        real (real64),   dimension(:), allocatable :: p 

        real (real64)   :: sum_of_vertex_distances, l
        integer (int32) :: i, j

        allocate(node_vector(number_of_dimensions))
        allocate(p          (number_of_dimensions))
        l = system_size(1) / number_of_field_nodes

        do i = 1, number_of_particles
            do j = 1, number_of_dimensions 
                node_vector(j) = floor(positions(j,i) / number_of_field_nodes * system_size(j))
            end do
            do j = 1, number_of_dimensions
                p(j) = positions(j,i) - node_vector(j) * system_size(j) / number_of_field_nodes
            end do
            
            ! Sum of all the distances |p - vertex(i,j,k)| for all the 8 nearest
            ! vertices.
            sum_of_vertex_distances = 5 * ( p(1)**2 + p(2)**2 + p(3)**2 ) + 3 * l**2 - 2 * l * (p(1) + p(2) + p(3))

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