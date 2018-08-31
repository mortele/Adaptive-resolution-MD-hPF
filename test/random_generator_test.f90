module random_generator_test
    use, intrinsic :: iso_fortran_env,  only: real64, int32
    use fruit,                          only: assert_equals    ! Subroutine
    use random_generator,               only: random_normal
    implicit none

    private

    private ::  compute_mean,               &
                compute_variance,           &
                compute_skewness,           &
                compute_kurtosis
    public  ::  setup,                      &
                teardown,                   &
                test_random_normal
contains
    subroutine setup
    end subroutine setup

    subroutine teardown
    end subroutine teardown

    subroutine test_random_normal()
        real (real64), dimension(:), allocatable :: values
        real (real64)   :: mean, variance, skewness, kurtosis
        real (real64)   :: mean_tollerance     = 0.01
        real (real64)   :: variance_tollerance = 0.1
        real (real64)   :: skewness_tollerance = 1.0
        real (real64)   :: kurtosis_tollerance = 1.0
        integer (int32) :: i
        
        allocate(values(int(1e5)))

        do i = 1, 20
            call random_normal(values)
            mean        = compute_mean    (values)
            variance    = compute_variance(values, mean)
            skewness    = compute_skewness(values, mean, variance)
            kurtosis    = compute_kurtosis(values, mean, variance)

            call assert_equals(0.0_real64, mean,     mean_tollerance,     char(i)//"  test_random_normal : Calculated mean is not sufficiently close to zero")
            call assert_equals(1.0_real64, variance, variance_tollerance, char(i)//"  test_random_normal : Calculated variance is not sufficiently close to zero")
            call assert_equals(0.0_real64, skewness, skewness_tollerance, char(i)//"  test_random_normal : Calculated skewness is not sufficiently close to zero")
            call assert_equals(0.0_real64, kurtosis, kurtosis_tollerance, char(i)//"  test_random_normal : Calculated kurtosis is not sufficiently close to zero")
        end do
    end subroutine test_random_normal

    pure function compute_mean(values) result(mean)
        real (real64), dimension(:), intent(in) :: values
        real (real64)   :: mean
        integer (int32) :: i
        
        mean = 0
        do i = 1, size(values)
            mean = mean + values(i)
        end do
        mean = mean / size(values)
    end function compute_mean

    pure function compute_variance(values, mean) result(variance)
        real (real64), dimension(:), intent(in) :: values
        real (real64),               intent(in) :: mean
        real (real64)   :: variance
        integer (int32) :: i
    
        variance = 0
        do i = 1, size(values)
            variance = variance + (values(i) - mean)**2
        end do

        ! If the mean is an _estimate_, we divide by N-1. If the mean was 
        ! known a priori, the denominator would be N.
        variance = variance / (size(values)) 
    end function compute_variance

    pure function compute_skewness(values, mean, variance) result(skewness)
        real (real64), dimension(:), intent(in) :: values
        real (real64),               intent(in) :: mean, variance
        real (real64)   :: skewness, standard_deviation
        integer (int32) :: i

        standard_deviation  = sqrt(variance)
        skewness            = 0
        do i = 1, size(values)
            skewness = skewness + ((values(i) - mean) / standard_deviation)**3
        end do
        skewness = skewness / size(values)
    end function compute_skewness

    function compute_kurtosis(values, mean, variance) result(kurtosis)
        real (real64), dimension(:), intent(in) :: values
        real (real64),               intent(in) :: variance, mean
        real (real64)   :: kurtosis, standard_deviation
        integer (int32) :: i
    
        standard_deviation  = sqrt(variance)
        kurtosis            = 0
        do i = 1, size(values)
            kurtosis = kurtosis + ((values(i) - mean) / standard_deviation)**4
        end do
        kurtosis = kurtosis / size(values) - 3.0
    end function compute_kurtosis
! https://www.researchgate.net/post/What_is_the_acceptable_range_of_skewness_and_kurtosis_for_normal_distribution_of_data
end module random_generator_test