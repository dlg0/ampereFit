program ampereFit

    use constants
    use ampFit_nameList
    use dlg
    !use spherHarmFns
    use ampFit_data
    use basisFns
    use ampFit_solve
    use write_output
    use ampFit_shift
    use ampFit_rotate

    implicit none

    call init_nameList ()

    call ampFit_read_data ()

    call ampFit_fill_structures ()

    nBFns = numberBFns ()

    call XYZ_to_SPH ( dataOriginal )

    !! *** check for track labelling mistake code ***

    call calculate_intersection_point ( dataOriginal )

    call calculate_shift_rotation_matrix ()

    call create_dataShifted ( dataOriginal, dataShifted )

    call XYZ_to_SPH ( dataShifted )

    call create_dataHalfSphere ( dataShifted )

    call create_bFns_at_data ( dataHalfSphere )

    call ampFit_solve_svd ( la_data_basisArr, dataHalfSphere )

    call write_data ( dataOriginal, fileName = 'ampData_original.nc' )
    call write_data ( dataHalfSphere, fileName = 'ampData_shifted.nc' )

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
