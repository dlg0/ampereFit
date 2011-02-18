program ampereFit

    use constants
    use ampFit_nameList
    use dlg
    use spherHarmFns
    use ampFit_data, &
    only: dataOriginal, nObs=>nSubSet, nVec, &
        ampFit_read_data, ampFit_fill_structures
    use basisFns
    use ampFit_solve
    use write_output

    implicit none

    call init_nameList ()

    call ampFit_read_data ()

    call ampFit_fill_structures ()

    nBFns = numberBFns ()

    call create_bFns_at_data ( dataOriginal )

    ! *** check for track labelling mistake code ***

    ! *** create shifted data structure *** 

    call ampFit_solve_svd ( la_data_basisArr, dataOriginal )

    call write_FAC ()

!!   Setup the basis functions at the observation locations
!!   ------------------------------------------------------
!
!    write(*,*) 'Generating basis functions ..'
!    nBFns   = numberBFns ()
!
!    write(*,*) 'nObs: ', nObs, 'nBFns: ', nBFns
!
!    allocate ( brBFnArr(nBFns,nObs), &
!                bThBFnArr(nBFns,nObs), &
!                bPhBFnArr(nBFns,nObs), &
!                YBFnArr(nBFns,nObs), &
!                alpha(nBFns,nBFns), &
!                beta_(nBFns,nBFns), &
!                coeffs(nBFns,1) )
!
!    call setupSHFns ( nObs, nBFns, rI,  &
!        coLat, mlt * 15.0, brBFnArr, bThBFnArr, bPhBFnArr, YBFnArr, &
!            minCoLat, maxCoLat ) 
!    write(*,*) 'DONE'
!    read(*,*)


     
end program ampereFit
