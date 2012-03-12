! Output netCDF routines for AMPERE data
!
! D.L. Green and C.L. Waters
! Centre for Space Physics
! University of NEwcastle
! Australia
!
! Ver : 201202
!
module ampFit_write_output

  use netcdf
  use constants
  use ampFit_data
 
  contains

  subroutine write_data ( dataIn, fileName )

    implicit none

    type(ampData) :: dataIn(:)

    character(len=*), intent(in), optional :: fileName
    character(len=100)  :: useFileName

    integer :: utc_id, bR_id, bT_id, bP_id, &
        nObs_id, nc_id, ipln_id, &
        T_id, P_id, R_id, &
        bX_id, bY_id, bZ_id, &
        X_id, Y_id, Z_id, sv_id, jPar_id

    integer :: nObs 

    nObs = size ( dataIn )

    if( present(fileName) ) then
        useFileName = fileName
    else
        useFileName = 'ampData.nc'
    endif
    
    write(*,*) 'Writing ', trim(useFileName), ' ...'

    call check_cdf ( nf90_create ( useFileName, nf90_clobber, nc_id ) )
    call check_cdf ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )
   
    call check_cdf ( nf90_def_var ( nc_id, 'utc', NF90_DOUBLE, (/ nObs_id /), utc_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'sv_number', NF90_DOUBLE, (/ nObs_id /), sv_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'iPln', NF90_DOUBLE, (/ nObs_id /), ipln_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'bR', NF90_DOUBLE, (/ nObs_id /), bR_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'bT', NF90_DOUBLE, (/ nObs_id /), bT_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'bP', NF90_DOUBLE, (/ nObs_id /), bP_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'bX', NF90_DOUBLE, (/ nObs_id /), bx_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'bY', NF90_DOUBLE, (/ nObs_id /), by_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'bZ', NF90_DOUBLE, (/ nObs_id /), bz_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'T', NF90_DOUBLE, (/ nObs_id /), T_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'P', NF90_DOUBLE, (/ nObs_id /), P_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'R', NF90_DOUBLE, (/ nObs_id /), R_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'X', NF90_DOUBLE, (/ nObs_id /), X_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'Y', NF90_DOUBLE, (/ nObs_id /), Y_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'Z', NF90_DOUBLE, (/ nObs_id /), Z_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'jPar', NF90_DOUBLE, (/ nObs_id /), jPar_id ) )
  
    call check_cdf ( nf90_enddef ( nc_id ) )
    
    call check_cdf ( nf90_put_var ( nc_id, utc_id, dataIn%utc ) )
    call check_cdf ( nf90_put_var ( nc_id, sv_id, dataIn%isat ) )
    call check_cdf ( nf90_put_var ( nc_id, ipln_id, dataIn%ipln ) )

    call check_cdf ( nf90_put_var ( nc_id, bR_id, dataIn%bR) )
    call check_cdf ( nf90_put_var ( nc_id, bT_id, dataIn%bT) )
    call check_cdf ( nf90_put_var ( nc_id, bP_id, dataIn%bP) )

    call check_cdf ( nf90_put_var ( nc_id, bx_id, dataIn%bX ) )
    call check_cdf ( nf90_put_var ( nc_id, by_id, dataIn%bY ) )
    call check_cdf ( nf90_put_var ( nc_id, bz_id, dataIn%bZ ) )

    call check_cdf ( nf90_put_var ( nc_id, T_id, dataIn%T ) )
    call check_cdf ( nf90_put_var ( nc_id, P_id, dataIn%P ) )
    call check_cdf ( nf90_put_var ( nc_id, R_id, dataIn%R ) )

    call check_cdf ( nf90_put_var ( nc_id, X_id, dataIn%x ) )
    call check_cdf ( nf90_put_var ( nc_id, Y_id, dataIn%y ) )
    call check_cdf ( nf90_put_var ( nc_id, Z_id, dataIn%z ) )

    call check_cdf ( nf90_put_var ( nc_id, jPar_id, dataIn%jPar ) )

    call check_cdf ( nf90_close ( nc_id ) )

  end subroutine write_data
!
! ------------------------------------------------------------------------------
!
  subroutine write_griddata ( fname, &
                            yy, mnth, dd, hh, mn, ss, avg_int, &
                            Kmax, Mmax, resol_deg, nLat, nLon, &
                            coLat_rad, lon_rad, &
                            dBthet_th, dbPhi_th,&
                            dBthet_ph, dBPhi_ph, jPar )

    implicit none

    character(len=*), intent(in) :: fname

    integer, intent(in) :: yy, mnth, dd, hh, mn, ss, avg_int, &
                           KMax, Mmax, nLat, nLon

    real(kind=DBL), intent(in) :: coLat_rad(:), lon_rad(:), &
                                  dBthet_th(:), dBPhi_th(:), &
                                  dBthet_ph(:), dBPhi_ph(:), jPar(:)

    real(kind=DBL), intent(in) :: resol_deg

    real(kind=DBL), allocatable :: clat_deg(:), lon_hr(:)
 
    integer :: yy_id, mnth_id, dd_id, hh_id, mn_id, ss_id, nc_id
    integer :: avg_int_id, kmax_id, mmax_id, resol_id, nLat_id, nLon_id
    integer :: nObs_id, clat_deg_id, lon_hr_id
    integer :: dbth_th_id, dbph_th_id, dbth_ph_id, dbph_ph_id, j_id

    integer :: nObs 

    nObs = size ( coLat_rad )

! Create file
    call check_cdf ( nf90_create ( fname, nf90_clobber, nc_id ) )

! Define mode
    call check_cdf ( nf90_def_var ( nc_id, 'year', NF90_INT, yy_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'month', NF90_INT, mnth_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'day', NF90_INT, dd_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'hour', NF90_INT, hh_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'minut', NF90_INT, mn_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'sec', NF90_INT, ss_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'avgint', NF90_INT, avg_int_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'kmax', NF90_INT, kmax_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'mmax', NF90_INT, mmax_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'res_deg', NF90_DOUBLE, resol_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'nLatGrid', NF90_INT, nLat_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'nLonGrid', NF90_INT, nLon_id ) )
   
    call check_cdf ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'cLat_deg', NF90_DOUBLE, (/ nObs_id /), clat_deg_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'lon_hr', NF90_DOUBLE, (/ nObs_id /), lon_hr_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'dbth_th', NF90_DOUBLE, (/ nObs_id /), dbth_th_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbph_th', NF90_DOUBLE, (/ nObs_id /), dbph_th_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbth_ph', NF90_DOUBLE, (/ nObs_id /), dbth_ph_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbph_ph', NF90_DOUBLE, (/ nObs_id /), dbph_ph_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'jPar', NF90_DOUBLE, (/ nObs_id /), j_id ) )
 
    call check_cdf ( nf90_enddef ( nc_id ) )

! Data mode
    call check_cdf ( nf90_put_var ( nc_id, yy_id, yy) )
    call check_cdf ( nf90_put_var ( nc_id, mnth_id, mnth) )
    call check_cdf ( nf90_put_var ( nc_id, dd_id, dd) )
    call check_cdf ( nf90_put_var ( nc_id, hh_id, hh) )
    call check_cdf ( nf90_put_var ( nc_id, mn_id, mn) )
    call check_cdf ( nf90_put_var ( nc_id, ss_id, ss) )
    call check_cdf ( nf90_put_var ( nc_id, avg_int_id, avg_int) )
    call check_cdf ( nf90_put_var ( nc_id, kmax_id, Kmax) )
    call check_cdf ( nf90_put_var ( nc_id, mmax_id, Mmax) )
    call check_cdf ( nf90_put_var ( nc_id, resol_id, resol_deg) )
    call check_cdf ( nf90_put_var ( nc_id, nLat_id, nLat) )
    call check_cdf ( nf90_put_var ( nc_id, nLon_id, nLon) )

    allocate ( clat_deg(nObs), lon_hr(nObs) )
! convert coLat rad to degrees
    cLat_deg = coLat_rad*180.0/pi
! convert phi (rad) to decimal hours)
    lon_hr = (lon_rad*180.0/pi)/15.0

    call check_cdf ( nf90_put_var ( nc_id, clat_deg_id, cLat_deg) )
    call check_cdf ( nf90_put_var ( nc_id, lon_hr_id, lon_hr) )
    call check_cdf ( nf90_put_var ( nc_id, dbth_th_id, dBthet_th) )
    call check_cdf ( nf90_put_var ( nc_id, dbph_th_id, dBPhi_th) )
    call check_cdf ( nf90_put_var ( nc_id, dbth_ph_id, dBthet_ph) )
    call check_cdf ( nf90_put_var ( nc_id, dbph_ph_id, dBPhi_ph) )

    call check_cdf ( nf90_put_var ( nc_id, j_id, jPar ) )

    call check_cdf ( nf90_close ( nc_id ) )

  end subroutine write_griddata
!
! --------------------------------------------------------------------------
!
  subroutine write_irddata ( fname, &
                            yy, mnth, dd, hh, mn, ss, avg_int, &
                            coLat_rad, lon_rad, &
                            dBthet_th, dbPhi_th,&
                            dBthet_ph, dBPhi_ph, &
                            ipln_a, isat_a, &
                            qual_a, splice_a )

    implicit none

    character(len=*), intent(in) :: fname

    integer, intent(in) :: yy, mnth, dd, hh, mn, ss, avg_int

    real(kind=DBL), intent(in) :: coLat_rad(:), lon_rad(:), &
                                  dBthet_th(:), dBPhi_th(:), &
                                  dBthet_ph(:), dBPhi_ph(:)

    integer, intent(in) :: ipln_a(:), isat_a(:), splice_a(:)
    real, intent(in) :: qual_a(:)

    real(kind=DBL), allocatable :: clat_deg(:), lon_hr(:)
 
    integer :: yy_id, mnth_id, dd_id, hh_id, mn_id, ss_id, nc_id
    integer :: avg_int_id, ipln_id, isat_id, qual_id, splice_id
    integer :: nObs_id, clat_deg_id, lon_hr_id
    integer :: dbth_th_id, dbph_th_id, dbth_ph_id, dbph_ph_id

    integer :: nObs 

    nObs = size ( coLat_rad )

! Create file
    call check_cdf ( nf90_create ( fname, nf90_clobber, nc_id ) )

! Define mode
    call check_cdf ( nf90_def_var ( nc_id, 'year', NF90_INT, yy_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'month', NF90_INT, mnth_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'day', NF90_INT, dd_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'hour', NF90_INT, hh_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'minut', NF90_INT, mn_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'sec', NF90_INT, ss_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'avgint', NF90_INT, avg_int_id ) )
   
    call check_cdf ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'cLat_deg', NF90_DOUBLE, (/ nObs_id /), clat_deg_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'lon_hr', NF90_DOUBLE, (/ nObs_id /), lon_hr_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'dbth_th', NF90_DOUBLE, (/ nObs_id /), dbth_th_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbph_th', NF90_DOUBLE, (/ nObs_id /), dbph_th_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbth_ph', NF90_DOUBLE, (/ nObs_id /), dbth_ph_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'dbph_ph', NF90_DOUBLE, (/ nObs_id /), dbph_ph_id ) )

    call check_cdf ( nf90_def_var ( nc_id, 'SatPln', NF90_INT, (/ nObs_id /), ipln_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'iSat', NF90_INT, (/ nObs_id /), isat_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'qual', NF90_DOUBLE, (/ nObs_id /), qual_id ) )
    call check_cdf ( nf90_def_var ( nc_id, 'splice', NF90_INT, (/ nObs_id /), splice_id ) )
 
    call check_cdf ( nf90_enddef ( nc_id ) )

! Data mode
    call check_cdf ( nf90_put_var ( nc_id, yy_id, yy) )
    call check_cdf ( nf90_put_var ( nc_id, mnth_id, mnth) )
    call check_cdf ( nf90_put_var ( nc_id, dd_id, dd) )
    call check_cdf ( nf90_put_var ( nc_id, hh_id, hh) )
    call check_cdf ( nf90_put_var ( nc_id, mn_id, mn) )
    call check_cdf ( nf90_put_var ( nc_id, ss_id, ss) )
    call check_cdf ( nf90_put_var ( nc_id, avg_int_id, avg_int) )

    allocate ( clat_deg(nObs), lon_hr(nObs) )
! convert coLat rad to degrees
    cLat_deg = coLat_rad*180.0/pi
! convert phi (rad) to decimal hours)
    lon_hr = (lon_rad*180.0/pi)/15.0

    call check_cdf ( nf90_put_var ( nc_id, clat_deg_id, cLat_deg) )
    call check_cdf ( nf90_put_var ( nc_id, lon_hr_id, lon_hr) )
    call check_cdf ( nf90_put_var ( nc_id, dbth_th_id, dBthet_th) )
    call check_cdf ( nf90_put_var ( nc_id, dbph_th_id, dBPhi_th) )
    call check_cdf ( nf90_put_var ( nc_id, dbth_ph_id, dBthet_ph) )
    call check_cdf ( nf90_put_var ( nc_id, dbph_ph_id, dBPhi_ph) )

    call check_cdf ( nf90_put_var ( nc_id, ipln_id, ipln_a ) )
    call check_cdf ( nf90_put_var ( nc_id, isat_id, isat_a ) )
    call check_cdf ( nf90_put_var ( nc_id, qual_id, qual_a ) )
    call check_cdf ( nf90_put_var ( nc_id, splice_id, splice_a ) )

    call check_cdf ( nf90_close ( nc_id ) )

  end subroutine write_irddata

end module ampFit_write_output
!
! ===================================================================
