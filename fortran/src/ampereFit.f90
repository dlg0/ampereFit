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

    integer :: south
    integer :: trk_order(6)
! structure to hold input data + ghosts + extra ghosts
    type(ampData), allocatable :: all_data(:)
    real(kind=DBL) :: clat_lim

    integer :: i, j     ! for testing - clw

    call init_nameList ()

    call ampFit_read_data ()

    call ampFit_fill_structures ()

    nBFns = numberBFns ()

    call XYZ_to_SPH ( dataOriginal )
! output check with IDL - CLW
!    open(unit=1,file='test_sph.dat',form='formatted')
!    do i=1,size(dataOriginal)
!      write(1,10),dataOriginal(i)%GEI_R_km, &
!              dataOriginal(i)%GEI_coLat_deg, &
!              dataOriginal(i)%GEI_lon_deg, &
!              dataOriginal(i)%br_GEI, &
!              dataOriginal(i)%bTheta_GEI, &
!              dataOriginal(i)%bPhi_GEI
!    enddo
!10  format(6f12.3)
!    close(unit=1)

    !! *** check for track labelling mistake code ***

    call calculate_intersection_point ( dataOriginal )

    call calculate_shift_rotation_matrix ()

    call create_dataShifted ( dataOriginal, dataShifted )

    call XYZ_to_SPH ( dataShifted )

! added CLW 23 aug
    call sort_struc ( dataShifted, trk_order )
    south = 0
    call tag_lon_strays ( dataShifted, south )
!    south = 1
!    call tag_lon_strays ( dataShifted, south )
! - - - - - - - - - - - - - - - - - - - - - - - - 
    
    call create_dataHalfSphere ( dataShifted )

! add CLW aug 2011
    clat_lim = 25.0
    call ampFit_ghost(dataHalfSphere, all_data, clat_lim, trk_order)

!    call create_bFns_at_data ( dataHalfSphere )
    call create_bFns_at_data ( all_data )

!    call ampFit_solve_svd ( la_data_basisArr, dataHalfSphere )
    call ampFit_solve_svd ( la_data_basisArr, all_data )

    call write_data ( dataOriginal, fileName = 'ampData_original.nc' )
!    call write_data ( dataHalfSphere, fileName = 'ampData_shifted.nc' )
    call write_data ( all_data, fileName = 'ampData_shifted.nc' )

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
