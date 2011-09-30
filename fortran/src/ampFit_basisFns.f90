module basisFns

use constants
use fgsl
use ampFit_namelist
use la_precision, only: WP=>DP

implicit none

real(kind=DBL), allocatable :: nArr(:)
integer, allocatable :: mArr(:), kArr(:)
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
 
    integer :: k, l, m, n, s
    integer :: i, cnt, span, nObs
    real(kind=DBL) :: x, lon, coLat, r, dYdr_rDep, factorial1, factorial2
    real(DBL), allocatable :: rDep(:), norm(:)

    ! FGSL 

    integer(fgsl_int) :: fgsl_stat

    !write(*,*) 'Creating basis functions at dataIn locations ...'

    nObs    = size ( dataIn )

    allocate ( basis(nObs,nBFns) )
    if (.not. allocated ( mArr )) then 
        allocate ( mArr(nBFns), nArr(nBfns), kArr(nBfns) )
    endif

    r = 1.0

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

            !! Normalised
            !fgsl_stat = fgsl_sf_legendre_sphplm_deriv_array ( &
            !        maxK, abs(m), x, &
            !        basis(i,cnt:cnt+span-1)%PLM, & 
            !        basis(i,cnt:cnt+span-1)%dPLM )

            ! Not normalised
            fgsl_stat = fgsl_sf_legendre_plm_deriv_array ( &
                    maxK, abs(m), x, &
                    basis(i,cnt:cnt+span-1)%PLM, & 
                    basis(i,cnt:cnt+span-1)%dPLM )

            mArr(cnt:cnt+span-1) = m
            nArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 
            kArr(cnt:cnt+span-1) = (/ (i,i=abs(m),maxK) /) 

            !rDep = 1d0
            allocate(rDep(span),norm(span))
            dYdR_rDep = 1d0

            do s=1,span
                rDep(s) = 1.0!r**nArr(cnt+s-1) * ( 1.0 + nArr(cnt+s-1) / ( nArr(cnt+s-1) + 1.0 ) &
                    !* ( 1.0 / r )**( 2.0 * nArr(cnt+s-1) + 1.0 ) )
                !l = nArr(cnt+s-1)
                !factorial1 = fgsl_sf_fact ( l-m )
                !factorial2 = fgsl_sf_fact ( l+m )
                !norm(s) = sqrt ( (2*l+1)/(4*pi) ) * sqrt ( factorial1/factorial2 )
                !!write(*,*) l, m, factorial1, factorial2, norm(s), basis(i,cnt+s-1)%Y
            enddo

            if(m>=0) then
                basis(i,cnt:cnt+span-1)%Y = rDep * basis(i,cnt:cnt+span-1)%PLM * cos ( abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%br = dYdr_rDep * basis(i,cnt:cnt+span-1)%PLM * cos (abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%bTh = 1d0 / r * rDep * basis(i,cnt:cnt+span-1)%dPLM * cos ( abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%bPh = -abs(m) * rDep / ( r * sin ( coLat ) ) &
                    * basis(i,cnt:cnt+span-1)%PLM * sin ( abs(m) * lon )! / norm
            else
                basis(i,cnt:cnt+span-1)%Y = rDep * basis(i,cnt:cnt+span-1)%PLM * sin ( abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%br = dYdr_rDep * basis(i,cnt:cnt+span-1)%PLM * sin (abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%bTh = 1d0 / r * rDep * basis(i,cnt:cnt+span-1)%dPLM * sin ( abs(m) * lon )! / norm
                basis(i,cnt:cnt+span-1)%bPh = -abs(m) * rDep / ( r * sin ( coLat ) ) &
                    * basis(i,cnt:cnt+span-1)%PLM * cos ( abs(m) * lon )! / norm
            endif


            deallocate(rDep,norm)

            cnt = cnt + span

        enddo

    enddo

    !deallocate ( dataIn_PLMBFnArr, dataIn_dPLMBFnArr )

end subroutine create_bFns_at_data

end module basisFns
