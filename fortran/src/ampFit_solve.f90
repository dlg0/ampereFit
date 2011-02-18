module ampFit_solve

use la_precision, only: WP=>DP
use f95_lapack, only: la_gelss, la_lange
use constants
use basisFns, only: nBFns

real(kind=DBL), allocatable :: fit_bTheta_GEI(:), fit_bPhi_GEI(:)

contains

subroutine ampFit_solve_svd ( la_A, dataIn ) 

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    real(kind=WP), intent(inout) :: la_A(:,:)

    real(kind=WP), allocatable :: la_B(:), la_s(:)
    real(kind=WP) :: r1
    integer :: la_rank, la_info, nObs, M, N
    real(kind=DBL), allocatable :: fitted_B(:), tmpA(:,:), tmpB(:)
    
    write(*,*) 'Solving Ax=B system for x ...'

    ! M = number of data points
    ! N = number of basis functions

    nObs = size ( dataIn )

    N = size ( la_A, 2 ) 
    M = size ( la_A, 1 )
    
    write(*,*) '    N: ', N
    write(*,*) '    M: ', M

    checkShape: &
    if (nObs*2 /= M ) then 

        write(*,*) '*** ERROR: Shape mismatch in ampFit_solve_svd.'
        write(*,*) '         : ', shape(la_A), shape(dataIn)
        write(*,*) '         : ', 'M: ', M, ' N: ', N

        stop

    endif checkShape

    allocate ( la_B(M), la_s(N) )

    la_B(1:nObs) = dataIn%bTheta_GEI
    la_B(nObs+1:nObs*2) = dataIn%bPhi_GEI

    write(*,*) '    Shape A: ', shape(la_A)
    write(*,*) '    Shape B: ', shape(la_B)

    allocate ( tmpA(M,N), tmpB(M) )
    tmpA = la_A
    tmpB = la_B

    if(size(la_B)/=max(size(la_A,1),size(la_A,2))) then
        write(*,*) '*** ERROR: size B not right'
        stop 
    endif

    if(size(la_s)/=min(size(la_A,1),size(la_A,2))) then
        write(*,*) '*** ERROR: size S not right'
        stop
    endif

    write(*,*) '    Calling la_gelss ...'
    call la_gelss ( la_A, la_B, la_rank, la_s, info = la_info ) 

    write(*,*) '    Rank: ', la_rank
    write(*,*) '    Info: ', la_info

    allocate ( fitted_B(M) )

    fitted_B = matMul ( tmpA, la_B(1:N) )

    write(*,*) fitted_B(1:20)
    write(*,*) tmpB(1:20)
    
    allocate ( fit_bTheta_GEI(nObs), fit_bPhi_GEI(nObs) )

    fit_bTheta_GEI  = fitted_B(1:nObs)
    fit_bPhi_GEI    = fitted_B(nObs+1:nObs*2)

    write(*,*) 'DONE' 

end subroutine ampFit_solve_svd

end module ampFit_solve
