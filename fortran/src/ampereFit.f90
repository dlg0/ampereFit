program ampereFit

    use constants
    use ampFit_nameList
    use dlg
    !use spherHarmFns
    use ampFit_data
    use ampFit_sort
    use basisFns
    use ampFit_solve
    use write_output
    use ampFit_shift
    use ampFit_rotate
    use ISO_C_BINDING
    use f_aacgm

    implicit none

        real(DBL) :: iLat, iLon, hgt
        real(DBL), target :: oLat, oLon, r
        integer :: yr, flg, s, yrsec
        real(DBL) :: mlon, mlt

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

! Try to call the aacgm C library

        iLat = 85.0
        iLon = 45.0
        hgt = 150.0

        flg = 0
        yr = 1990

        s = f_AACGMInit(yr)
        s = f_AACGMConvert(iLat, iLon, hgt, C_LOC(oLat), C_LOC(oLon), C_LOC(r), flg)

        yrsec = 3*24*3600
        mlon = 45.0

        mlt = f_MLTConvertYrsec(yr, yrsec, mlon)

        if(s/=0) then
                write(*,*) 'Error in f_AACGMConvert'
        else
                write(*,*) 'Input: ', iLat, iLon, hgt
                write(*,*) 'Output: ', oLat, oLon
                write(*,*) 'mlon, mlt: ', mlon, mlt
        endif

end program ampereFit
