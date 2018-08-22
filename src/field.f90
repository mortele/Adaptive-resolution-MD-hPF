module field
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use parameters,      only: field_nodes
    use system,          only: system_size
    implicit none
    private 

    real (real64), public, dimension(:,:,:),   allocatable :: density_field
    real (real64), public, dimension(:,:,:,:), allocatable :: density_gradient




contains
end module field