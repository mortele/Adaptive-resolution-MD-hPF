module random_generator
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    private :: gasdev_s, normalScalar, normalArray, normalMatrix
    public  :: random_normal

    ! Interface for the gasdev_s(harvest) subroutine, which enables calling it 
    ! with scalar, or array dimension(:) or dimension(:,:) arguments.
    interface random_normal
        module procedure normalScalar, normalArray, normalMatrix
    end interface

contains

    ! gasdev_s(harvest) function from numerical recipes in Fortran, 
    ! http://www.elch.chem.msu.ru/tch/group/FortranBooks/NumericalRecipesinF90.pdf
    ! 
    ! Uses the intrinsic fortran random_number(harvest) subroutine to generate
    ! two uniform random numbers, which are transformed into two Gaussian 
    ! distributed random numbers with zero mean and unit variance.
    subroutine gasdev_s(harvest)
        implicit none

        real (real64), intent(out)  :: harvest
        real (real64)               :: rsq, v1, v2
        real (real64), save         :: g
        logical, save               :: gausStored = .false.

        if (gausStored) then
            harvest = g
            gausStored = .false.
        else
            do
                call random_number(v1)
                call random_number(v2)
                v1  = 2_real64 * v1 - 1_real64
                v2  = 2_real64 * v2 - 1_real64
                rsq = v1**2 + v2**2
                if (rsq > 0.0 .and. rsq < 1.0) exit
            end do
            rsq     = sqrt(-2.0_real64 * log(rsq) / rsq)
            harvest = v1 * rsq
            g       = v2 * rsq
            gausStored = .true.
        end if
    end subroutine

    subroutine normalScalar(scalar)
        implicit none
        real (real64), intent(in out) :: scalar
        call gasdev_s(scalar)
    end subroutine

    subroutine normalArray(array)
        implicit none
        real    (real64), dimension(:), intent(in out) :: array
        integer (int32) :: i
        do i = 1, size(array,1)
            call normalScalar(array(i))
        end do
    end subroutine

    subroutine normalMatrix(matrix)
        implicit none
        real    (real64), dimension(:,:), intent(in out) :: matrix
        integer (int32) :: i
        do i = 1, size(matrix,2)
            call normalArray(matrix(:,i))
        end do
    end subroutine

  end module random_generator
