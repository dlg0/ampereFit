program test_aacgm

        use ISO_C_BINDING
        use f_aacgm
        use constants

        implicit none

        real(DBL) :: iLat, iLon, hgt
        real(DBL), target :: oLat, oLon, r
        integer :: yr, flg, s, yrsec
        real(DBL) :: mlon, mlt

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

end program test_aacgm
