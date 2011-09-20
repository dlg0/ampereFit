module ampFit_solve

use la_precision, only: WP=>DP
use f95_lapack, only: la_gelss, la_lange
use constants
use basisFns, only: nBFns

!real(kind=DBL), allocatable :: fit_bTheta_GEI(:), fit_bPhi_GEI(:)

contains

subroutine ampFit_solve_svd ( basis, dataIn, coeffs ) 

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampBasis), intent(in) :: basis(:,:)
    real(DBL), allocatable, intent(inout) ::  coeffs(:)

    real(kind=WP), allocatable :: la_A(:,:)
    real(kind=WP), allocatable :: la_B(:), la_s(:)
    real(kind=WP) :: r1
    integer :: la_rank, la_info, nObs, M, N
    real(kind=DBL), allocatable :: fitted_B(:), tmpA(:,:), tmpB(:)

    nObs = size ( dataIn )

    write(*,*) '    Size of basis array: ', nBFns*2d0*nObs*2d0*16/(1024d0**2), 'MB'

    allocate ( la_A(nObs*2,nBFns*2) )
    write(*,*) __FILE__, __LINE__
    write(*,*) "nObs: ", nObs, " nBFns: ", nBFns


    la_A(1:nObs,1:nBFns) = basis%bTh
    la_A(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    la_A(1:nObs,nBFns+1:2*nBFns) = -basis%bPh
    la_A(nObs+1:2*nObs,1:nBFns) = basis%bPh

    write(*,*) 'DONE'

    write(*,*) 'Solving Ax=B system for x ...'

    ! M = number of data points
    ! N = number of basis functions

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

    la_B(1:nObs) = dataIn%bT
    la_B(nObs+1:nObs*2) = dataIn%bP

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

    allocate(coeffs(N))
    coeffs  = la_B(1:N)

    !allocate ( fitted_B(M) )

    !fitted_B = matMul ( tmpA, ampFit_solve_svd ) !la_B(1:N) )

    !allocate ( fit_bTheta_GEI(nObs), fit_bPhi_GEI(nObs) )

    !fit_bTheta_GEI  = fitted_B(1:nObs)
    !fit_bPhi_GEI    = fitted_B(nObs+1:nObs*2)

    write(*,*) 'DONE' 

end subroutine ampFit_solve_svd


subroutine ampFit_sumBasis ( basis, dataIn, coeffs )

    type(ampBasis), intent(in) :: basis(:,:)
    type(ampData), intent(inout) :: dataIn(:)
    real(DBL), intent(in) :: coeffs(:)

    real(DBL), allocatable :: tmpA(:,:)
    integer :: nObs, M, N
    real(kind=DBL), allocatable :: fitted_B(:)

    nObs = size ( basis, 1 )

    if(nObs /= size(dataIn)) then

        write(*,*) __FILE__, __LINE__
        stop

    endif

    if(size(coeffs) /= size(basis,2)*2) then

        write(*,*) __FILE__, __LINE__
        write(*,*) "ERROR:"
        write(*,*) "    size(coeffs): ", size(coeffs)
        write(*,*) "    size(basis,2): ", size(basis,2)
        stop

    endif

    allocate ( tmpA(nObs*2,nBFns*2) )

    tmpA(1:nObs,1:nBFns) = basis%bTh
    tmpA(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    tmpA(1:nObs,nBFns+1:2*nBFns) = -basis%bPh
    tmpA(nObs+1:2*nObs,1:nBFns) = basis%bPh

    N = size ( tmpA, 2 ) 
    M = size ( tmpA, 1 )
    
    allocate(fitted_B(M))
    fitted_B = matMul ( tmpA, coeffs )

    dataIn%bT  = fitted_B(1:nObs)
    dataIn%bP  = fitted_B(nObs+1:nObs*2)

end subroutine ampFit_sumBasis

end module ampFit_solve
