program ampereFit

    use constants
    use read_nameList
    use netcdf
    use dlg
    use spherHarmFns
    use lin_sol_svd_int
    use imsl_libraries
    implicit none

    integer :: nc_id, nObs_id, coLat_id, mlt_id, &
        bTh_id, bPh_id
    integer :: nObs, nBFns
    real(kind=DBL), allocatable :: coLat(:), mlt(:), bTh(:), bPh(:)
    real(kind=DBL), allocatable :: brBFnArr(:,:), bThBFnArr(:,:), &
        bPhBFnArr(:,:), YBFnArr(:,:)

    integer :: nSingular, dim_ids(1)
    real(kind=DBL), parameter :: zero = 1.0d0
    type(d_options) :: iOpt_(1) = d_options(0,zero)
    real(kind=DBL), allocatable :: S_(:)
    real(kind=DBL), parameter :: small = 1.0d-34
    real(kind=DBL), allocatable :: alpha(:,:), beta_(:,:), &
        coeffs(:,:)
    character(len=100) :: outFileName

!   Read the namelist variables
!   ---------------------------

    call init_nameList ()

!   Read in the delta b data from netCDF file
!   -----------------------------------------

    write(*,*) 'Reading ', trim ( deltab_fileName )

    call dlg_check ( nf90_open ( path = deltab_fileName, &
        mode = nf90_nowrite, ncid = nc_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'coLat', coLat_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'mlt', mlt_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'bTh', bTh_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'bPh', bPh_id ) )

    call dlg_check ( nf90_inquire_variable ( nc_id, coLat_id, &
        dimIds = dim_ids ) )
    call dlg_check ( nf90_inquire_dimension ( nc_id, dim_ids(1), &
        len = nObs ) )

    allocate ( coLat(nObs), mlt(nObs), bTh(nObs), bPh(nObs) )

    call dlg_check ( nf90_get_var ( nc_id, coLat_id, coLat ) )
    call dlg_check ( nf90_get_var ( nc_id, mlt_id, mlt ) )
    call dlg_check ( nf90_get_var ( nc_id, bTh_id, bTh ) )
    call dlg_check ( nf90_get_var ( nc_id, bPh_id, bPh ) )

    call dlg_check ( nf90_close ( nc_id ) )
    write (*,*) 'DONE'

!   Setup the basis functions at the observation locations
!   ------------------------------------------------------

    write(*,*) 'Generating basis functions ...'
    nBFns   = numberBFns ()

    write(*,*) 'nObs: ', nObs, 'nBFns: ', nBFns

    allocate ( brBFnArr(nBFns,nObs), &
                bThBFnArr(nBFns,nObs), &
                bPhBFnArr(nBFns,nObs), &
                YBFnArr(nBFns,nObs), &
                alpha(nBFns,nBFns), &
                beta_(nBFns,nBFns), &
                coeffs(nBFns,1) )

    call setupSHFns ( nObs, nBFns, rI,  &
        coLat, mlt * 15.0, brBFnArr, bThBFnArr, bPhBFnArr, YBFnArr, &
            minCoLat, maxCoLat ) 
    write(*,*) 'DONE'
    read(*,*)
!   Do the fit
!   ----------

    write(*,*) 'Doing the fit ...'
    alpha   = matMul ( bPhBFnArr, transpose ( bPhBFnArr ) )
    beta_   = transpose ( matMul ( reshape ( bPh, (/1,nObs/) ), transpose ( bPhBFnArr ) ) )

    iOpt_(1)    = d_options ( d_lin_sol_svd_set_small, small )

    call lin_sol_svd ( alpha, beta_, coeffs, rank = nSingular, S = S_, &
        iOpt = iOpt_)
    write(*,*) 'DONE'

!   Write the ouput file
!   --------------------

    write(*,*) 'Writing sh_debug.nc ...'
    outFileName = 'jParIrid.nc'
    
    call dlg_check ( nf90_create ( outFileName, nf90_clobber, nc_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )
    
    call dlg_check ( nf90_def_var ( nc_id, 'coLat', NF90_DOUBLE, &
        (/ nObs_id /), coLat_id ) )
    
    call dlg_check ( nf90_enddef ( nc_id ) )
    
    call dlg_check ( nf90_put_var ( nc_id, coLat_id, coLat ) )
    call dlg_check ( nf90_close ( nc_id ) )
    write(*,*) 'DONE'
      
end program ampereFit
