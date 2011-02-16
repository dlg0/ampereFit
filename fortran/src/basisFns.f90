module basisFns

use constants
use fgsl
use ampFit_namelist
use ampFit_data

implicit none

real(kind=DBL), allocatable :: &
    data_brBFnArr(:,:), &
    data_bThBFnArr(:,:), &
    data_bPhBFnArr(:,:), &
    data_YBFnArr(:,:), &
    nArr(:)
integer, allocatable :: mArr(:), kArr(:)

real(kind=DBL), allocatable :: basis(:,:)


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
        data_PLMBFnArr(nBFns,nObs), &
        data_dPLMBfnArr(nBFns,nObs), &
        data_YBfnArr(nBFns,nObs), &
        data_brBfnArr(nBFns,nObs), &
        data_bThBfnArr(nBFns,nObs), &
        data_bPhBfnArr(nBFns,nObs) )

    allocate ( mArr(nBFns), nArr(nBfns), kArr(nBfns) )

    r = rE

    do i=1,nSubSet

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
                    data_PLMBFnArr(cnt:cnt+span-1,i), & 
                    data_dPLMBFnArr(cnt:cnt+span-1,i) )

            rDep = 1d0
            dYdR_rDep = 1d0

            data_YBFnArr(cnt:cnt+span-1,i) = rDep * data_PLMBFnArr(cnt:cnt+span-1,i) * cos ( abs(m) * lon )
            data_brBFnArr(cnt:cnt+span-1,i) = dYdr_rDep * data_PLMBFnArr(cnt:cnt+span-1,i) * cos (abs(m) * lon )
            data_bThBFnArr(cnt:cnt+span-1,i) = 1d0 / r * rDep * data_dPLMBFnArr(cnt:cnt+span-1,i) * cos ( abs(m) * lon )
            data_bPhBFnArr(cnt:cnt+span-1,i) = -m * rDep / ( r * sin ( coLat ) ) &
                * data_PLMBFnArr(cnt:cnt+span-1,i) * sin ( m * lon )
            mArr(cnt:cnt+span-1) = m
            nArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 
            kArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 

            cnt = cnt + span

        enddo

    enddo

    deallocate ( data_PLMBFnArr, data_dPLMBFnArr )

    allocate ( basis(nBFns*2,nObs*2) )

    basis(1:nBFns,1:nObs) = data_bThBFnArr
    basis(nBFns+1:2*nBFns,nObs+1:2*nObs) = data_bThBFnArr
    basis(nBFns+1:2*nBFns,1:nObs) = -data_bPhBFnArr
    basis(1:nBFns,nObs+1:2*nObs) = data_bPhBFnArr

    !write(*,*) '    mArr: ', mArr
    !write(*,*) '    kArr: ', kArr

    write(*,*) '    Size of basis array: ', nBFns*2d0*nObs*2d0*16/(1024d0**2), 'MB'

    write(*,*) 'DONE'

end subroutine create_bFns_at_data

end module basisFns
