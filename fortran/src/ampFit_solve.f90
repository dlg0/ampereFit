module ampFit_solve
!
! subroutine ampFit_solve_svd ( basis, dataIn, coeffs ) 
! subroutine ampFit_sumBasis ( basis, dataIn, coeffs )
!
! D.L. Green and C.L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Sep 2011
!
! ver : 201202
!
use la_precision, only: WP=>DP
use f95_lapack, only: la_gelss, la_lange, la_gesdd
use constants
use basisFns, only: nBFns, nArr
use ampFit_rotate

contains

subroutine ampFit_solve_svd ( basis, dataIn, coeffs ) 

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampBasis), intent(in) :: basis(:,:)
    real(DBL), allocatable, intent(inout) ::  coeffs(:)

    real(kind=WP), allocatable :: la_A(:,:), U(:,:), VT(:,:)
    real(kind=WP), allocatable :: V(:,:), alpha(:,:)
    real(kind=WP), allocatable :: la_B(:), la_s(:), W(:), beta(:)
    real(kind=WP) :: r1, ss, TOL
    integer :: la_rank, la_info, nObs, i, j, jj
    real(kind=DBL), allocatable :: fitted_B(:), tmpA(:,:), tmp(:)

    nObs = size ( dataIn )
    TOL = 0.1              ! SVD rejection level on singular values

    write(*,*) '    Size of basis array: ', nBFns*2d0*nObs*2d0*16/(1024d0**2), 'MB'

    allocate ( U(2*nBFns, 2*nBFns))
    allocate (VT(2*nBFns, 2*nBFns))
    allocate ( W(2*nBFns))
    allocate (la_s(2*nBFns))

    allocate ( la_A(nObs*2, nBFns*2))
!    write(*,*) __FILE__, __LINE__
    write(*,*) "nObs: ", nObs, " nBFns: ", nBFns

! la_A holds the basis function values at all input data locations
! la_A should be 2 lots of [nObs by nBFns]
    la_A(1:nObs,1:nBFns) = basis%bTh
    la_A(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    la_A(1:nObs,nBFns+1:2*nBFns) = basis%bPh
    la_A(nObs+1:2*nObs,1:nBFns) = -basis%bPh

    write(*,*) 'Solving Ax=B system ...'

    ! M = number of data points
    ! N = number of basis functions

    allocate (alpha(2*nBFns, 2*nBFns))
    alpha = matmul ( transpose(la_A),la_A )
!    N = min(M,N)     ! alpha makes things square
!    M = N

!    N = size ( la_A, 2 ) 
!    M = size ( la_A, 1 )
    
!    write(*,*) '    N: ', N
!    write(*,*) '    M: ', M

    call la_gesdd(alpha, la_s, u, vt, w, 'N', la_info)

    if (la_info .ne. 0) then
      err_stat = 10
      err_msg = 'Error in SVD solver'
      write(*,*) err_msg
      print*,'la_info=',la_info
      return
    endif

! singular values, W are in la_s
! print*,'la_info=',la_info

    deallocate (alpha)

! Weed out the w_i
    do i=1,2*nBFns
      if (la_s(i) .le. TOL) la_s(i) = 0.0
    enddo

! Load in the dB data - Apply quality flag to data here
    allocate ( la_B(2*nObs) )
    do i=1,nObs
      la_B(i) = dataIn(i)%bT/dataIn(i)%qual
      la_B(nObs+i) = dataIn(i)%bP/dataIn(i)%qual

!      write(10,*) dataIn(i)%T*180.0/pi, dataIn(i)%P*180.0/pi, &
!                  dataIn(i)%bT, dataIn(i)%bP, &
!                  la_B(i), la_B(nObs+i)           ! for compare with IDL

    enddo

    allocate (beta(2*nBFns))
    beta = matmul(transpose(la_A),la_B)

    deallocate (la_A)
    allocate (V(2*nBFns, 2*nBFns))
    V = transpose(VT)

    checkShape: &
    if (nObs*2 /= size(la_A, 1) ) then 
      err_stat = 11
      err_msg = 'ERROR: Array shape mismatch in ampFit_solve_svd'
      write(*,*) err_msg
      write(*,*) ' : ', shape(la_A), shape(dataIn)
      return 
    endif checkShape

! load in the dB data

!    allocate ( tmpA(M,N), tmpB(M) )
!    tmpA = la_A
!    tmpB = la_B
!
!    if(size(la_B)/=max(size(la_A,1),size(la_A,2))) then
!        write(*,*) '*** ERROR: size B not right'
!        stop 
!    endif
!
!    if(size(la_s)/=min(size(la_A,1),size(la_A,2))) then
!        write(*,*) '*** ERROR: size S not right'
!        stop
!    endif
!
!    write(*,*) '    Calling la_gelss ...'
!    call la_gelss ( la_A, la_B, la_rank, la_s, info = la_info ) 

!    write(*,*) '    Rank: ', la_rank
!    write(*,*) '    Info: ', la_info

    allocate (coeffs(2*nBFns))
    allocate (tmp(2*nObs))
!    coeffs  = la_B(1:N)
! Back-substitution (assumes square matrix)
    do j=1,2*nBFns
      ss = 0.0
      if (la_s(j) .ne. 0.0) then
        do i=1,2*nBFns
          ss = ss+u(i,j)*beta(i)
        enddo
        ss = ss/la_s(j)
      endif
      tmp(j) = ss
    enddo
    do j=1,2*nBFns
      ss = 0.0
      do jj=1,2*nBFns
        ss = ss+v(j,jj)*tmp(jj)
      enddo
      coeffs(j) = ss
    enddo

!print*,'nBF,N=',nBFns,N
!    do i=1,N
!      write(3,*) coeffs(i)
!    enddo
!print*,'Stopeed in ampFit_solve_svd'
! stop

!    write(*,*) 'DONE' 

end subroutine ampFit_solve_svd
!
! ------------------------------------------------------------------------------------
!
subroutine ampFit_sumBasis ( basis, dataIn, coeffs )

    type(ampBasis), intent(in) :: basis(:,:)
    type(ampData), intent(inout) :: dataIn(:)
    real(DBL), intent(in) :: coeffs(:)

    real(DBL), allocatable :: tmpA(:,:)
    integer :: nObs, M, N, j
    real(kind=DBL), allocatable :: fitted_B(:), jPar(:)

    nObs = size ( basis, 1 )

    if(nObs /= size(dataIn)) then
      err_stat = 12
      err_msg = 'Incorrect Array size in sumBasis'
      write(*,*) err_msg
      print*,'nObs=',nObs
      return
    endif

    if(size(coeffs) /= size(basis,2)*2) then
      err_stat = 13
      err_msg = 'Incorrect Array size in sumBasis'
      write(*,*) err_msg
      print*,'Coeffs array size=',size(coeffs)
      return
    endif

    allocate ( tmpA(nObs*2,nBFns*2) )

    tmpA(1:nObs,1:nBFns) = basis%bTh
    tmpA(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    tmpA(1:nObs,nBFns+1:2*nBFns) = basis%bPh
    tmpA(nObs+1:2*nObs,1:nBFns) = -basis%bPh

    N = size ( tmpA, 2 ) 
    M = size ( tmpA, 1 )
    
    allocate(fitted_B(M))
    fitted_B = matMul ( tmpA, coeffs )

    allocate(jPar(size(dataIn)))

    jPar = matMul ( basis%Y, coeffs(nBFns+1:)*(-nArr*(nArr+1.0)) ) / (u0_*(rE+rSat)) * 1e-9 * 1e6! uAm^-2
    !jPar = matMul ( basis%Y, coeffs(1:nBFns)*(-nArr*(nArr+1.0)) ) / (u0_*(rE+rSat)) * 1e-9 * 1e6! uAm^-2

    dataIn%jPar = jPar
    dataIn%bT  = fitted_B(1:nObs)
    dataIn%bP  = fitted_B(nObs+1:nObs*2)

! Update the cartesian components from these new spherical ones
    call rtp_to_xyz_vec ( dataIn )

  end subroutine ampFit_sumBasis

end module ampFit_solve
!
! =====================================================================================
