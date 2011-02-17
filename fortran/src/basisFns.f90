module basisFns

use constants
use fgsl
use ampFit_namelist
use la_precision, only: WP=>DP

implicit none

real(kind=DBL), allocatable :: &
    data_brBFnArr(:,:), &
    data_bThBFnArr(:,:), &
    data_bPhBFnArr(:,:), &
    data_YBFnArr(:,:), &
    nArr(:)
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


subroutine create_bFns_at_data ( dataIn )

    use fgsl
    use ampFit_nameList
    use ampFit_data, only: nObs=>nSubSet

    implicit none
  
    type(ampData), intent(in) :: dataIn(:)
 
    integer :: k, l, m, n
    integer :: i, cnt, span
    real(kind=DBL) :: x, lon, coLat, r, rDep, dYdr_rDep

    real(kind=DBL), allocatable :: &
        data_PLMBFnArr(:,:), &
        data_dPLMBFnArr(:,:)
 
    ! FGSL 

    integer(fgsl_int) :: fgsl_stat

    write(*,*) 'Creating basis functions at data locations ...'

    allocate ( &
        data_PLMBFnArr(nObs,nBFns), &
        data_dPLMBfnArr(nObs,nBFns), &
        data_YBfnArr(nObs,nBFns), &
        data_brBfnArr(nObs,nBFns), &
        data_bThBfnArr(nObs,nBFns), &
        data_bPhBfnArr(nObs,nBFns) )

    allocate ( mArr(nBFns), nArr(nBfns), kArr(nBfns) )

    r = rE

    do i=1,nObs

        cnt = 1
        do m=-maxM,maxM

            ! l, m, x
            span = fgsl_sf_legendre_array_size ( maxK, abs(m) ) 

            coLat = dataIn(i)%GEI_coLat_rad
            lon = dataIn(i)%GEI_lon_rad
            x = cos(coLat)

            !write(*,*) maxK, m, x, cnt, span, i 

            ! The GSL routines assume you want the l=0,m=0 term, so
            ! to get this to work I am running with this term included
            ! in the fit. It could be removed, but right now it should 
            ! be fine as its coefficient will be small.

            fgsl_stat = fgsl_sf_legendre_sphplm_deriv_array ( &
                    maxK, abs(m), x, &
                    data_PLMBFnArr(i,cnt:cnt+span-1), & 
                    data_dPLMBFnArr(i,cnt:cnt+span-1) )

            rDep = 1d0
            dYdR_rDep = 1d0

            data_YBFnArr(i,cnt:cnt+span-1) = rDep * data_PLMBFnArr(i,cnt:cnt+span-1) * cos ( abs(m) * lon )
            data_brBFnArr(i,cnt:cnt+span-1) = dYdr_rDep * data_PLMBFnArr(i,cnt:cnt+span-1) * cos (abs(m) * lon )
            data_bThBFnArr(i,cnt:cnt+span-1) = 1d0 / r * rDep * data_dPLMBFnArr(i,cnt:cnt+span-1) * cos ( abs(m) * lon )
            data_bPhBFnArr(i,cnt:cnt+span-1) = -m * rDep / ( r * sin ( coLat ) ) &
                * data_PLMBFnArr(i,cnt:cnt+span-1) * sin ( m * lon )

            mArr(cnt:cnt+span-1) = m
            nArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 
            kArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 

            cnt = cnt + span

        enddo

    enddo

    deallocate ( data_PLMBFnArr, data_dPLMBFnArr )

    allocate ( la_data_basisArr(nObs*2,nBFns*2) )

    la_data_basisArr(1:nObs,1:nBFns) = data_bThBFnArr
    la_data_basisArr(nObs+1:2*nObs,nBFns+1:2*nBFns) = data_bThBFnArr
    la_data_basisArr(1:nObs,nBFns+1:2*nBFns) = -data_bPhBFnArr
    la_data_basisArr(nObs+1:2*nObs,1:nBFns) = data_bPhBFnArr

    !write(*,*) '    mArr: ', mArr
    !write(*,*) '    kArr: ', kArr

    write(*,*) '    Size of basis array: ', nBFns*2d0*nObs*2d0*16/(1024d0**2), 'MB'

    write(*,*) 'DONE'

end subroutine create_bFns_at_data

end module basisFns
