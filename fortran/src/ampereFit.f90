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

    integer :: trk_order(6)
! structure to hold modified (ghost added) data
    type(ampData), allocatable :: &
        dataHalfSphere_ghosted(:), &
        dataHalfSphere_ghosted_fit(:)

    real(kind=DBL) :: clat_lim

        ! AACGM Grid vars
        integer :: nLatGrid, nLonGrid
        real(DBL), allocatable :: aacgmColatGrid_deg(:,:), &
                aacgmLonGrid_deg(:,:), aacgmHgtGrid_km(:,:)
        real(DBL), allocatable, target :: &
                geogCoLat_deg(:,:), &
                geogLat_deg(:,:), &
                geogLon_deg(:,:), &
                geogHgt_km(:,:) 
        type(ampData), allocatable :: dataGrid(:), dataGridShifted(:)

        real :: gridCoLat_deg, gridHgt_km, latStep, lonStep
        integer :: i, j, idx

        real(DBL) :: iLat, iLon, hgt
        real(DBL), target :: oLat, oLon, r
        integer :: yr, flg, s, yrsec
        real(DBL) :: mlon, mlt

        real(DBL), allocatable :: coeffs(:)

    call init_nameList ()

    call ampFit_read_data ()

    call ampFit_fill_structures ()

    nBFns = numberBFns ()

    call XYZ_to_SPH ( dataOriginal )

    call calculate_intersection_point ( dataOriginal )

    call calculate_shift_rotation_matrix ()

    call create_dataShifted ( dataOriginal, dataShifted )

    call XYZ_to_SPH ( dataShifted )

    call sort_struc ( dataShifted, trk_order )

    call tag_lon_strays ( dataShifted, south )

    call create_dataHalfSphere ( dataShifted )

    clat_lim = 25.0
    call ampFit_ghost ( dataHalfSphere, dataHalfSphere_ghosted, clat_lim, trk_order )

    call create_bFns_at_data ( dataHalfSphere_ghosted, dataBFn )

    allocate(coeffs(nBFns*2))
    coeffs = ampFit_solve_svd ( dataBFn, dataHalfSphere_ghosted )

    allocate(dataHalfSphere_ghosted_fit(size(dataHalfSphere_ghosted)))
    dataHalfSphere_ghosted_fit = dataHalfSphere_ghosted 
    call ampFit_sumBasis ( dataBFn, dataHalfSphere_ghosted, coeffs )

    call write_data ( dataOriginal, fileName = 'output/ampData_original.nc' )
    call write_data ( dataHalfSphere, fileName = 'output/ampData_shifted.nc' )
    call write_data ( dataHalfSphere_ghosted, fileName = 'output/ampData_ghosted.nc' )
    call write_data ( dataHalfSphere_ghosted_fit, fileName = 'output/ampData_ghosted_fit.nc' )

! Try to call the aacgm C library

        write(*,*) 'Creating AACGM Grid ...'

        nLatGrid = 10
        nLonGrid = 6

        allocate ( aacgmCoLatGrid_deg(nLatGrid,nLonGrid), &
                aacgmLonGrid_deg(nLatGrid,nLonGrid), &
                aacgmHgtGrid_km(nLatGrid,nLonGrid), &
                geogCoLat_deg(nLatGrid,nLonGrid), &
                geogLat_deg(nLatGrid,nLonGrid), &
                geogLon_deg(nLatGrid,nLonGrid), &
                geogHgt_km(nLatGrid,nLonGrid) )

        allocate ( dataGrid(nLatGrid*nLonGrid) )

        write(*,*) 'nGrid: ', nLatGrid*nLonGrid

        gridColat_deg = 40.0
        latStep = gridCoLat_deg / (nLatGrid+1)
        lonStep = 360.0 / (nLonGrid)
        gridHgt_km = 780.0

        flg = 0 ! 0 -> to aacgm, 1 -> geographic
        yr = 1990

        s = f_AACGMInit(yr)

        do i=1,nLatGrid
                do j=1,nLonGrid
                        aacgmCoLatGrid_deg(i,j) = i * latStep
                        aacgmLonGrid_deg(i,j) = j * lonStep
                        aacgmHgtGrid_km(i,j) = gridHgt_km

                        s = f_AACGMConvert(90.0-aacgmCoLatGrid_deg(i,j), &
                                aacgmLonGrid_deg(i,j), aacgmHgtGrid_km(i,j),&
                                C_LOC(geogLat_deg(i,j)), &
                                C_LOC(geogLon_deg(i,j)), &
                                C_LOC(geogHgt_km(i,j)), flg)

                        if(s/=0) then
                                write(*,*) 'Error in f_AACGMConvert'
                                write(*,*) 'Input: ', 90.0-aacgmCoLatGrid_deg(i,j), &
                                        aacgmLonGrid_deg(i,j), aacgmHgtGrid_km(i,j), flg
                        endif

                        idx = (i-1)*nLonGrid+j

                        dataGrid(idx)%T = (90.0-geogLat_deg(i,j))* degToRad
                        dataGrid(idx)%P = (geogLon_deg(i,j)) * degToRad

                        dataGrid(idx)%R = aacgmHgtGrid_km(i,j)

                        dataGrid(idx)%bX = 0.0
                        dataGrid(idx)%bY = 0.0
                        dataGrid(idx)%bZ = 0.0

                        !write(*,*) dataGrid(idx)%GEI_coLat_rad, &
                        !        dataGrid(idx)%GEI_lon_rad, &
                        !        dataGrid(idx)%GEI_coLat_deg, &
                        !        dataGrid(idx)%GEI_lon_deg, &
                        !        dataGrid(idx)%GEI_R_km


                enddo
        enddo

        call XYZ_from_SPH ( dataGrid )
        allocate ( dataGridShifted (size(dataGrid)) )
        call create_dataShifted ( dataGrid, dataGridShifted )
        call create_bFns_at_data ( dataGridShifted, gridBFn )

        iLat = 85.0
        iLon = 45.0
        hgt = 150.0

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
