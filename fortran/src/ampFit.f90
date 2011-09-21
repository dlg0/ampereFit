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
    use ampFit_geivec_aacgmvec

    implicit none

    integer :: trk_order(6)

    type(ampData), allocatable :: &
        dataHalfSphere(:), &
        dataHalfSphere_ghosted(:), &
        dataHalfSphere_ghosted_fit(:), &
        dataHalfSphere_unshifted_fit(:)

    type(ampBasis), allocatable :: &
        dataBFn(:,:), gridBFn(:,:)

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
    type(ampData), allocatable :: &
            gridAACGM(:), &
            gridGEI(:), &
            gridShiftedGEI(:)

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

    call ShiftData ( dataOriginal, dataShifted )

    call sort_struc ( dataShifted, trk_order )

    call tag_lon_strays ( dataShifted, south )

    call create_dataHalfSphere ( dataShifted, dataHalfSphere )

    clat_lim = 45.0
    call ampFit_ghost ( dataHalfSphere, dataHalfSphere_ghosted, clat_lim, trk_order )

    call create_bFns_at_data ( dataHalfSphere_ghosted, dataBFn )

    call ampFit_solve_svd ( dataBFn, dataHalfSphere_ghosted, coeffs )

    allocate(dataHalfSphere_ghosted_fit(size(dataHalfSphere_ghosted)))
    dataHalfSphere_ghosted_fit = dataHalfSphere_ghosted 
    call ampFit_sumBasis ( dataBFn, dataHalfSphere_ghosted_fit, coeffs )

    allocate(dataHalfSphere_unshifted_fit(size(dataHalfSphere_ghosted_fit)))
    dataHalfSphere_unshifted_fit = dataHalfSphere_ghosted_fit
    call UnShiftData ( dataHalfSphere_ghosted_fit, dataHalfSphere_unshifted_fit )

    call write_data ( dataOriginal, fileName = 'output/ampData_original.nc' )
    call write_data ( dataHalfSphere, fileName = 'output/ampData_shifted.nc' )
    call write_data ( dataHalfSphere_ghosted, fileName = 'output/ampData_ghosted.nc' )
    call write_data ( dataHalfSphere_ghosted_fit, fileName = 'output/ampData_ghosted_fit.nc' )
    call write_data ( dataHalfSphere_unshifted_fit, fileName = 'output/ampData_unshifted_fit.nc' )

    ! Try to call the aacgm C library

    write(*,*) 'Creating AACGM Grid ...'

    nLatGrid = 70
    nLonGrid = 12

    allocate ( aacgmCoLatGrid_deg(nLatGrid,nLonGrid), &
            aacgmLonGrid_deg(nLatGrid,nLonGrid), &
            aacgmHgtGrid_km(nLatGrid,nLonGrid), &
            geogCoLat_deg(nLatGrid,nLonGrid), &
            geogLat_deg(nLatGrid,nLonGrid), &
            geogLon_deg(nLatGrid,nLonGrid), &
            geogHgt_km(nLatGrid,nLonGrid) )

    allocate ( gridGEI(nLatGrid*nLonGrid), gridAACGM(nLatGrid*nLonGrid) )

    write(*,*) 'nGrid: ', nLatGrid*nLonGrid

    gridColat_deg = 40.0
    latStep = gridCoLat_deg / (nLatGrid+1)
    lonStep = 360.0 / (nLonGrid)
    gridHgt_km = (rE + rSat) * 1e-3

    flg = 1 ! 0 -> to aacgm, 1 -> geographic
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
                            stop
                    endif

                    idx = (i-1)*nLonGrid+j

                    gridAACGM(idx)%bR = 0.0
                    gridAACGM(idx)%bT = 0.0
                    gridAACGM(idx)%bP = 0.0
 
                    gridAACGM(idx)%T = aacgmCoLatGrid_deg(i,j)* degToRad
                    gridAACGM(idx)%P = aacgmLonGrid_deg(i,j) * degToRad
                    gridAACGM(idx)%R = aacgmHgtGrid_km(i,j)

                    gridGEI(idx)%bR = 0.0
                    gridGEI(idx)%bT = 0.0
                    gridGEI(idx)%bP = 0.0

                    gridGEI(idx)%T = (90.0-geogLat_deg(i,j))* degToRad
                    gridGEI(idx)%P = (geogLon_deg(i,j)) * degToRad
                    gridGEI(idx)%R = aacgmHgtGrid_km(i,j)

            enddo
    enddo

    call SPH_to_XYZ ( gridAACGM )
    call SPH_to_XYZ ( gridGEI )

    call write_data ( gridAACGM, fileName = 'output/ampData_gridAACGM.nc' )
    call write_data ( gridGEI, fileName = 'output/ampData_gridGEI.nc' )

    allocate ( gridShiftedGEI (size(gridGEI)) )
    gridShiftedGEI = gridGEI
    call ShiftData ( gridGEI, gridShiftedGEI )
    call create_bFns_at_data ( gridShiftedGEI, gridBFn )
    call ampFit_sumBasis ( gridBFn, gridShiftedGEI, coeffs )
    call SPH_to_XYZ ( gridShiftedGEI )
    call UnShiftData ( gridShiftedGEI, gridGEI )

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


    call write_data ( gridGEI, fileName = 'output/ampData_gridUnShiftedGEI.nc' )
    call write_data ( gridShiftedGEI, fileName = 'output/ampData_gridShiftedGEI.nc' )

end program ampereFit
