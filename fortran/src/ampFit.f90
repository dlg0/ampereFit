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
    integer :: nLatGrid, nLonGrid, jGEOPACK
    real(DBL), allocatable :: aacgmColatGrid_deg(:,:), &
            aacgmLonGrid_deg(:,:), aacgmHgtGrid_km(:,:)
    real(DBL), allocatable, target :: &
            geogCoLat_deg(:,:), &
            geogLat_deg(:,:), &
            geogLon_deg(:,:), &
            geogHgt_km(:,:) 
    type(ampData), allocatable :: &
            gridAACGM(:), &
            gridGEO(:), &
            gridGEI(:), &
            gridShiftedGEI(:)

    real :: gridCoLat_deg, gridHgt_km, latStep, lonStep
    integer :: i, j, idx

    real(DBL) :: iLat, iLon, hgt
    real(DBL), target :: oLat, oLon, r
    integer :: flg, s, yearSecond
    integer :: year,dayOfYear,month,day,hour,minute
    real(DBL) :: mlon, mlt, second
    real(DBL) vGSEX, vGSEY, vGSEZ
    real(DBL) :: bT_AACGM_T, bP_AACGM_T, bT_AACGM_P, bP_AACGM_P, &
                thisT_AACGM_deg, thisP_AACGM_deg

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

    call write_data ( dataOriginal, fileName = 'output/ampData_dataGEIRaw.nc' )
    call write_data ( dataHalfSphere, fileName = 'output/ampData_dataGEIShifted.nc' )
    call write_data ( dataHalfSphere_ghosted, fileName = 'output/ampData_dataGEIShiftedGhosted.nc' )
    call write_data ( dataHalfSphere_ghosted_fit, fileName = 'output/ampData_dataFitGEIShiftedGhosted.nc' )
    call write_data ( dataHalfSphere_unshifted_fit, fileName = 'output/ampData_dataFitGEIUnShiftedGhosted.nc' )

    ! Try to call the aacgm C library

    write(*,*) 'Creating AACGM Grid ...'

    nLatGrid = 50
    nLonGrid = 24

    allocate ( aacgmCoLatGrid_deg(nLatGrid,nLonGrid), &
            aacgmLonGrid_deg(nLatGrid,nLonGrid), &
            aacgmHgtGrid_km(nLatGrid,nLonGrid), &
            geogCoLat_deg(nLatGrid,nLonGrid), &
            geogLat_deg(nLatGrid,nLonGrid), &
            geogLon_deg(nLatGrid,nLonGrid), &
            geogHgt_km(nLatGrid,nLonGrid) )

    allocate ( gridGEO(nLatGrid*nLonGrid), gridGEI(nLatGrid*nLonGrid), gridAACGM(nLatGrid*nLonGrid) )

    write(*,*) 'nGrid: ', nLatGrid*nLonGrid

    gridColat_deg = 40.0
    latStep = gridCoLat_deg / (nLatGrid+1)
    lonStep = 360.0 / (nLonGrid)
    gridHgt_km = rSat * 1e-3


    write(*,*) 'Calculating and initialising the time/date variables ...'

    year    = 2010
    month   = 08
    day     = 24
    hour    = (sHr+eHr)/2.0
    minute  = 0
    second  = 0.0

    ! Default solar wind vector (see GEOPACK documentation)
    vGSEX   =  -400.0
    vGSEY   =  0.0
    vGSEZ   =  0.0 

    yearSecond = f_TimeYMDHMSToYrsec ( year,month,day,hour,minute,second )
    dayOfYear = yearSecond / (60*60*24)

    write(*,*) '    year: ', year
    write(*,*) '    month: ', month
    write(*,*) '    day: ', day
    write(*,*) '    hour: ', hour
    write(*,*) '    minute: ', minute
    write(*,*) '    second: ', second
    write(*,*) '    yearSecond: ', yearSecond
    write(*,*) '    dayOfYear: ', dayOfYear

    call RECALC_08 (year,dayOfYear,hour,minute,second,vGSEX,vGSEY,vGSEZ)

    s = f_AACGMInit(year)

    write(*,*) 'DONE'


    flg = 1 ! 0 -> to aacgm, 1 -> geographic
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
                    gridAACGM(idx)%R = aacgmHgtGrid_km(i,j) + rE*1e-3

                    gridGEO(idx)%bR = 0.0
                    gridGEO(idx)%bT = 0.0
                    gridGEO(idx)%bP = 0.0

                    gridGEO(idx)%T = (90.0-geogLat_deg(i,j))* degToRad
                    gridGEO(idx)%P = (geogLon_deg(i,j)) * degToRad
                    gridGEO(idx)%R = geogHgt_km(i,j) + rE*1e-3


            enddo
    enddo

    call SPH_to_XYZ ( gridGEO )
    call write_data ( gridGEO, fileName = 'output/ampData_gridGEO.nc' )

    jGEOPACK = -1
    do i=1,nLatGrid
            do j=1,nLonGrid
     
                ! jGEOPACK<1 GEO->GEI, jGEOAPCK>1 GEI->GEO

                idx = (i-1)*nLonGrid+j

                call GEIGEO_08 ( gridGEI(idx)%X,gridGEI(idx)%Y,gridGEI(idx)%Z,&
                                    gridGEO(idx)%X,gridGEO(idx)%Y,gridGEO(idx)%Z,jGEOPACK )

                gridGEI(idx)%bX = 0.0
                gridGEI(idx)%bY = 0.0
                gridGEI(idx)%bZ = 0.0

            enddo
    enddo

    call XYZ_to_SPH ( gridGEI )
 
    call SPH_to_XYZ ( gridAACGM )
    call SPH_to_XYZ ( gridGEI )

    call write_data ( gridAACGM, fileName = 'output/ampData_gridAACGM.nc' )
    call write_data ( gridGEI, fileName = 'output/ampData_gridGEI.nc' )

    allocate ( gridShiftedGEI (size(gridGEI)) )
    gridShiftedGEI = gridGEI
    call ShiftData ( gridGEI, gridShiftedGEI )
    
    call write_data ( gridShiftedGEI, fileName = 'output/ampData_gridShiftedGEI.nc' )

    call create_bFns_at_data ( gridShiftedGEI, gridBFn )
    call ampFit_sumBasis ( gridBFn, gridShiftedGEI, coeffs )
    call SPH_to_XYZ ( gridShiftedGEI )
    call UnShiftData ( gridShiftedGEI, gridGEI )

    ! Rotate vectors from GEI to AACGM
    do i=1,nLatGrid
            do j=1,nLonGrid

                idx = (i-1)*nLonGrid+j

                call geisph_to_aacgmvec (gridGEI(idx)%R, gridGEI(idx)%T, gridGEI(idx)%P, &
                        gridGEI(idx)%bT, gridGEI(idx)%bP, &
                        bT_AACGM_T, bP_AACGM_T, bT_AACGM_P, bP_AACGM_P, &
                        thisT_AACGM_deg, thisP_AACGM_deg )

                write(*,*) 'AACGM coord comparision: ', i, j, idx
                write(*,*) gridAACGM(idx)%T*radToDeg, gridAACGM(idx)%P*radToDeg
                write(*,*) 90-thisT_AACGM_deg, thisP_AACGM_deg
                write(*,*)

                gridAACGM(idx)%bT = bT_AACGM_T
                gridAACGM(idx)%bP = bP_AACGM_P

           enddo
    enddo


    iLat = 85.0
    iLon = 45.0
    hgt = 150.0

    !yrsec = 3*24*3600
    mlon = 45.0

    mlt = f_MLTConvertYrsec(year, yearSecond, mlon)

    if(s/=0) then
            write(*,*) 'Error in f_AACGMConvert'
    else
            write(*,*) 'Input: ', iLat, iLon, hgt
            write(*,*) 'Output: ', oLat, oLon
            write(*,*) 'mlon, mlt: ', mlon, mlt
    endif


    call write_data ( gridGEI, fileName = 'output/ampData_gridFitUnShiftedGEI.nc' )
    call write_data ( gridShiftedGEI, fileName = 'output/ampData_gridFitShiftedGEI.nc' )
    call write_data ( gridAACGM, fileName = 'output/ampData_gridFitAACGM.nc' )

end program ampereFit
