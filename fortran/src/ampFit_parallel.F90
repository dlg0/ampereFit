module parallel

#ifdef _parallel_

use constants, only: DBL
use basisFns, only: nBFns

implicit none

! MPI / ScaLAPACK variables

integer :: nRow, nCol, nObs
integer :: iContext,npRow,npCol,myRow,myCol, &
    rowBlockSize,colBlockSize,nRowLocal,nColLocal, &
    rowStartProc,colStartProc
integer :: iAm, nProcs
integer :: descA( 9 ), descB( 9 )


real(kind=DBL), allocatable :: la_A(:,:), la_B(:)

contains

subroutine initParallelEnv ()

    implicit none

    integer :: info, lld
    integer, external :: numRoc

    write(*,*) 'Setting up MPI ...' 

    nRow = nObs*2
    nCol = nBFns*2

    npRow = 2
    npCol = 2

    rowBlockSize = 64 ! blocking factors
    colBlockSize = 64

    rowStartProc = 0
    colStartProc = 0

    info = 0

    !call sl_init(iContext,npRow,npCol)
    !call blacs_gridInfo(iContext,npRow,npCol,myRow,myCol)

    call blacs_pInfo ( iAm, nProcs )

    !write(*,*) 'proc id: ', iAm, nProcs

    if ( nProcs <= 1 ) then 

        write(*,*) 'blacs_pInfo backup needed'

        if ( iAm == 0 ) nProcs = npRow * npCol
        call blacs_setup ( iAm, nProcs )

    endif

    if ( nProcs /= npRow*npCol ) then 

        write(*,*)
        write(*,*)
        write(*,*) 'CONFIG ERROR:'
        write(*,*) '-------------'
        write(*,*) '    nProcs /= npRow * npCol'
        write(*,*) '    Please correct and re-run'
        write(*,*) '    Have a nice day :)'
        write(*,*)
        write(*,*)
 
        return  
                
    endif

    call blacs_get ( -1, 0, iContext )
    call blacs_gridInit ( iContext, 'Row-major', npRow, npCol )
    call blacs_gridInfo ( iContext, npRow, npCol, myRow, myCol )

    nRowLocal = numRoc(nRow,rowBlockSize,myRow,rowStartProc,npRow)
    nColLocal = numRoc(nCol,colBlockSize,myCol,colStartProc,npCol)

    lld = max ( 1, nRowLocal)

    write(*,*) nRow, nCol, rowBlockSize, colBlocksize, rowStartProc, colStartProc, iContext, lld, info
    call descInit( descA, nRow, nCol, rowBlockSize, colBlockSize, &
        rowStartProc, colStartProc, iContext, lld, info )

    write(*,*) nRow, nCol, rowBlockSize, colBlocksize, rowStartProc, colStartProc, iContext, lld, info
    call descInit( descB, nRow, 1   , rowBlockSize, colBlockSize, &
        rowStartProc, colStartProc, iContext, lld, info )

    allocate ( la_A(nRowLocal,nColLocal), la_B(nRowLocal) )


end subroutine initParallelEnv

#endif

end module parallel
