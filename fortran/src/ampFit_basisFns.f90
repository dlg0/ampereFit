module basisFns
!
!  function numberBFns
!  subroutine create_bFns_at_data ( dataIn, basis )
!
! D.L. Green and C.L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Ver : 201202
! 

use constants
use ampFit_legendre
!use fgsl
!use la_precision, only: WP=>DP

implicit none

real(kind=DBL), allocatable :: nArr(:)
integer, allocatable :: mArr(:), kArr(:)
integer :: nBFns

contains

  integer function numberBFns (MaxK, MaxM)
! assumes full sphere
! Modified from IDL code by CLW

    implicit none

    integer, intent(in) :: MaxK, MaxM
    integer :: i, k, n_left

    numberBFns  = 0

    do i = 0, maxM
      numberBFns = numberBFns + (2*i + 1)
    enddo
    n_left = maxK - maxM
    numberBFns = numberBFns + n_left*(2*maxM+1) - 1

    return

  end function numberBFns
!
! ---------------------------------------------------------------------------------------
!

  subroutine create_bFns_at_data ( dataIn, basis, maxK, maxM )

    use ampFit_data

    implicit none
  
    type(ampData), intent(in) :: dataIn(:)
    type(ampBasis), allocatable, intent (inout) :: basis(:,:)
    integer, intent(in) :: maxK, maxM
 
    integer :: k, l, m, n, s
    integer :: i, j, cnt, nObs
    real(kind=DBL) :: x, lon, coLat, r, dYdr_rDep, Plm, dPlm, rDep

    nObs    = size ( dataIn )
    nBFns = numberBFns (maxK, maxM)

    allocate ( basis(nObs,nBFns) )
    if (.not. allocated ( mArr )) then 
        allocate ( mArr(nBFns), nArr(nBfns), kArr(nBfns) )
    endif

    r = 1.0

    do i=1,nObs       ! for every data point

      coLat = dataIn(i)%T
      lon = dataIn(i)%P

      cnt = 1         ! basis ftn counter
      dYdR_rDep = 1d0

      do m = 0,maxM

        do j = m, maxK
          if (j .gt. 0) then
            call dLegendre(j,abs(m),coLat,Plm,dPlm)
            basis(i,cnt)%PLM = Plm 
            basis(i,cnt)%dPLM = dPlm 

            mArr(cnt) = m
            nArr(cnt) = j 
            kArr(cnt) = j 

            rDep = r**nArr(cnt)*(1.0 + nArr(cnt)/(nArr(cnt) + 1.0 ) &
                    *(1.0 / r )**(2.0 * nArr(cnt) + 1.0 ) )

            basis(i,cnt)%Y = rDep * Plm * cos ( m * lon )
            basis(i,cnt)%br = dYdr_rDep * PLM * cos (m * lon )
            basis(i,cnt)%bTh = 1d0 / r * rDep * dPlm * cos ( m * lon )
            basis(i,cnt)%bPh = -m * rDep / ( r * sin ( coLat ) ) &
                * Plm * sin ( m * lon )

            mArr(cnt) = m
            nArr(cnt) = j 
            kArr(cnt) = j 

            cnt = cnt + 1

            if ( m .gt. 0) then   ! exclude m=0
              basis(i,cnt)%PLM = Plm 
              basis(i,cnt)%dPLM = dPlm 
              basis(i,cnt)%Y = rDep * Plm * sin ( m * lon )
              basis(i,cnt)%br = dYdr_rDep * Plm * sin ( m * lon )
              basis(i,cnt)%bTh = 1d0 / r * rDep * dPlm * sin ( m * lon )
              basis(i,cnt)%bPh = m * rDep / ( r * sin ( coLat ) ) &
                * Plm * cos ( m * lon )
              mArr(cnt) = -m
              nArr(cnt) =  j 
              kArr(cnt) =  j 

              cnt = cnt + 1
            endif

          endif  ! do not count the L=0,M=0 term

        enddo    ! L(j) loop for spherical harmonics
      enddo      ! M loop

    enddo        ! input data loop

  end subroutine create_bFns_at_data

end module basisFns
!
! ====================================================================
