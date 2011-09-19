module basisFns

use constants
use fgsl
use ampFit_namelist
use la_precision, only: WP=>DP

implicit none

real(kind=DBL), allocatable :: &
!    data_brBFnArr(:,:), &
!    data_bThBFnArr(:,:), &
!    data_bPhBFnArr(:,:), &
!    data_YBFnArr(:,:), &
    nArr(:)

type(ampBasis), allocatable :: dataBFn(:,:), gridBFn(:,:)

integer, allocatable :: mArr(:), kArr(:)
real(kind=WP), allocatable :: la_data_basisArr(:,:)
integer :: nBFns

contains

integer function numberBFns ()

    use ampFit_nameList
    implicit none
    integer :: m, k

    write(*,*) 'Calculating number of basis funcitons ...'
    write(*,*) '    basisType: ', basisType

    numberBFns  = 0

    do m = -maxM, maxM
        do k = 0, maxK


            select case ( basisType )
            case ('FULL')

                if ( k >= abs(m) ) numberBFns = numberBFns + 1

            case ('CAP0')

                if ( mod ( k - abs ( m ), 2 ) == 0  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1

            case ('CAP1')

                if ( mod ( k - abs ( m ), 2 ) == 1  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1

            case ('BOTH')

                if ( mod ( k - abs ( m ), 2 ) == 0  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1
                if ( mod ( k - abs ( m ), 2 ) == 1  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1
 
            end select

        enddo
    enddo

    write(*,*) '    nBfns: ', numberBFns
    write(*,*) 'DONE'

    return

end function numberBFns


subroutine create_bFns_at_data ( dataIn, basis )

    use fgsl
    use ampFit_nameList
    use ampFit_data

    implicit none
  
    type(ampData), intent(in) :: dataIn(:)
    type(ampBasis), allocatable, intent (inout) :: basis(:,:)
 
    integer :: k, l, m, n
    integer :: i, cnt, span, nObs
    real(kind=DBL) :: x, lon, coLat, r, rDep, dYdr_rDep

    !real(kind=DBL), allocatable :: &
    !    dataIn_PLMBFnArr(:,:), &
    !    dataIn_dPLMBFnArr(:,:)
 
    ! FGSL 

    integer(fgsl_int) :: fgsl_stat

    !write(*,*) 'Creating basis functions at dataIn locations ...'

    nObs    = size ( dataIn )

    !allocate ( &
    !    dataIn_PLMBFnArr(nObs,nBFns), &
    !    dataIn_dPLMBfnArr(nObs,nBFns), &
    !    dataIn_YBfnArr(nObs,nBFns), &
    !    dataIn_brBfnArr(nObs,nBFns), &
    !    dataIn_bThBfnArr(nObs,nBFns), &
    !    dataIn_bPhBfnArr(nObs,nBFns) )

    allocate ( basis(nObs,nBFns) )
    if (.not. allocated ( mArr )) then 
        allocate ( mArr(nBFns), nArr(nBfns), kArr(nBfns) )
    endif

    r = rE

    do i=1,nObs

        cnt = 1
        do m=-maxM,maxM

            ! l, m, x
            span = fgsl_sf_legendre_array_size ( maxK, abs(m) ) 

            coLat = dataIn(i)%T
            lon = dataIn(i)%P
            x = cos(coLat)

            ! The GSL routines assume you want the l=0,m=0 term, so
            ! to get this to work I am running with this term included
            ! in the fit. It could be removed, but right now it should 
            ! be fine as its coefficient will be small.

            fgsl_stat = fgsl_sf_legendre_sphplm_deriv_array ( &
                    maxK, abs(m), x, &
                    basis(i,cnt:cnt+span-1)%PLM, & 
                    basis(i,cnt:cnt+span-1)%dPLM )

            rDep = 1d0
            dYdR_rDep = 1d0

            basis(i,cnt:cnt+span-1)%Y = rDep * basis(i,cnt:cnt+span-1)%PLM * cos ( abs(m) * lon )
            basis(i,cnt:cnt+span-1)%br = dYdr_rDep * basis(i,cnt:cnt+span-1)%PLM * cos (abs(m) * lon )
            basis(i,cnt:cnt+span-1)%bTh = 1d0 / r * rDep * basis(i,cnt:cnt+span-1)%dPLM * cos ( abs(m) * lon )
            basis(i,cnt:cnt+span-1)%bPh = -m * rDep / ( r * sin ( coLat ) ) &
                * basis(i,cnt:cnt+span-1)%PLM * sin ( m * lon )

            mArr(cnt:cnt+span-1) = m
            nArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 
            kArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 

            cnt = cnt + span

        enddo

    enddo

    !deallocate ( dataIn_PLMBFnArr, dataIn_dPLMBFnArr )

end subroutine create_bFns_at_data

end module basisFns
