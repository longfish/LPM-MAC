INCLUDE 'mkl_pardiso.f90'

Program test
    implicit none
    real totals
    integer I
    !$omp parallel do private (I) &
    !$omp& reduction(+:totals)
        do I=1,100
            totals = totals + I
        enddo
    !$omp end parallel do
    print *, 'The calculated sum total is', totals

END PROGRAM test