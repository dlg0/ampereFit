program ampereFit

    use constants
    use ampFit_nameList
    use netcdf
    use dlg
    use spherHarmFns
    use ampFit_data
    use basisFns
!    use lin_sol_svd_int
!    use imsl_libraries

    implicit none

    real(kind=DBL), allocatable :: brBFnArr(:,:), bThBFnArr(:,:), &
        bPhBFnArr(:,:), YBFnArr(:,:)

    real(kind=DBL), parameter :: zero = 1.0d0
    !type(d_options) :: iOpt_(1) = d_options(0,zero)
    real(kind=DBL), allocatable :: S_(:)
    real(kind=DBL), parameter :: small = 1.0d-34
    real(kind=DBL), allocatable :: alpha(:,:), beta_(:,:), &
        coeffs(:,:)
    character(len=100) :: outFileName

    call init_nameList ()

    call ampFit_read_data ()

    call ampFit_fill_structures ()

    nBFns = numberBFns ()

    call create_bFns_at_data ( dataOriginal )
    

    ! *** check for track labelling mistake code ***

    ! *** create shifted data structure *** 


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


    write(*,*) '************************'
    write(*,*) 'ERROR: Need to replace the IMSL svd with LAPACK'
!!   Do the fit
!!   ----------
!
!    write(*,*) 'Doing the fit ...'
!    alpha   = matMul ( bPhBFnArr, transpose ( bPhBFnArr ) )
!    beta_   = transpose ( matMul ( reshape ( bPh, (/1,nObs/) ), transpose ( bPhBFnArr ) ) )
!
!    iOpt_(1)    = d_options ( d_lin_sol_svd_set_small, small )
!
!    call lin_sol_svd ( alpha, beta_, coeffs, rank = nSingular, S = S_, &
!        iOpt = iOpt_)
!    write(*,*) 'DONE'

!   Write the ouput file
!   --------------------

    write(*,*) 'Writing sh_debug.nc ...'
    outFileName = 'jParIrid.nc'
    
    call dlg_check ( nf90_create ( outFileName, nf90_clobber, nc_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nVec', nVec, nVec_id ) )
   
    call dlg_check ( nf90_def_var ( nc_id, 'time', NF90_DOUBLE, (/ nObs_id /), time_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'b_eci', NF90_DOUBLE, (/ nObs_id, nVec_id /), b_eci_id ) )
   
    call dlg_check ( nf90_enddef ( nc_id ) )
    
    call dlg_check ( nf90_put_var ( nc_id, time_id, time ) )
    call dlg_check ( nf90_put_var ( nc_id, b_eci_id, b_eci ) )

    call dlg_check ( nf90_close ( nc_id ) )
    write(*,*) 'DONE'
      
end program ampereFit
