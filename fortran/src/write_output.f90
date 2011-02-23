module write_output

use netcdf
use ampFit_solve
use dlg
use constants
 
contains

subroutine write_data ( dataIn, fileName )

    implicit none

    type(ampData) :: dataIn(:)

    character(len=*), intent(in), optional :: fileName
    character(len=100)  :: useFileName

    integer :: utc_id, bR_GEI_id, bT_GEI_id, bP_GEI_id, &
        nObs_id, nc_id, &
        fit_bT_GEI_id, fit_bP_GEI_id, &
        data_GEI_coLat_rad_id, data_GEI_lon_rad_id, data_GEI_R_km_id, &
        bx_GEI_id, by_GEI_id, bz_GEI_id, &
        data_GEI_px_id, data_GEI_py_id, data_GEI_pz_id

    integer :: nObs 

    nObs = size ( dataIn )

    if( present(fileName) ) then
        useFileName = fileName
    else
        useFileName = 'ampData.nc'
    endif
    
    write(*,*) 'Writing ', trim(useFileName), ' ...'

    call dlg_check ( nf90_create ( useFileName, nf90_clobber, nc_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )
   
    call dlg_check ( nf90_def_var ( nc_id, 'utc', NF90_DOUBLE, (/ nObs_id /), utc_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_bR_GEI', NF90_DOUBLE, (/ nObs_id /), bR_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_bT_GEI', NF90_DOUBLE, (/ nObs_id /), bT_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_bP_GEI', NF90_DOUBLE, (/ nObs_id /), bP_GEI_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_bx_GEI', NF90_DOUBLE, (/ nObs_id /), bx_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_by_GEI', NF90_DOUBLE, (/ nObs_id /), by_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_bz_GEI', NF90_DOUBLE, (/ nObs_id /), bz_GEI_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'fit_bT_GEI', NF90_DOUBLE, (/ nObs_id /), fit_bT_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'fit_bP_GEI', NF90_DOUBLE, (/ nObs_id /), fit_bP_GEI_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_coLat_rad', &
            NF90_DOUBLE, (/ nObs_id /), data_GEI_coLat_rad_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_lon_rad', &
            NF90_DOUBLE, (/ nObs_id /), data_GEI_lon_rad_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_R_km', &
            NF90_DOUBLE, (/ nObs_id /), data_GEI_R_km_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_px', NF90_DOUBLE, (/ nObs_id /), data_GEI_px_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_py', NF90_DOUBLE, (/ nObs_id /), data_GEI_py_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_pz', NF90_DOUBLE, (/ nObs_id /), data_GEI_pz_id ) )
  
    call dlg_check ( nf90_enddef ( nc_id ) )
    
    call dlg_check ( nf90_put_var ( nc_id, utc_id, dataIn%utc ) )

    call dlg_check ( nf90_put_var ( nc_id, bR_GEI_id, dataIn%bR_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, bT_GEI_id, dataIn%bTheta_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, bP_GEI_id, dataIn%bPhi_GEI ) )

    call dlg_check ( nf90_put_var ( nc_id, bx_GEI_id, dataIn%dbx ) )
    call dlg_check ( nf90_put_var ( nc_id, by_GEI_id, dataIn%dby ) )
    call dlg_check ( nf90_put_var ( nc_id, bz_GEI_id, dataIn%dbz ) )

    call dlg_check ( nf90_put_var ( nc_id, fit_bT_GEI_id, fit_bTheta_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, fit_bP_GEI_id, fit_bPhi_GEI ) )

    call dlg_check ( nf90_put_var ( nc_id, data_GEI_coLat_rad_id, dataIn%GEI_coLat_rad ) )
    call dlg_check ( nf90_put_var ( nc_id, data_GEI_lon_rad_id, dataIn%GEI_lon_rad ) )
    call dlg_check ( nf90_put_var ( nc_id, data_GEI_R_km_id, dataIn%GEI_R_km ) )

    call dlg_check ( nf90_put_var ( nc_id, data_GEI_px_id, dataIn%px ) )
    call dlg_check ( nf90_put_var ( nc_id, data_GEI_py_id, dataIn%py ) )
    call dlg_check ( nf90_put_var ( nc_id, data_GEI_pz_id, dataIn%pz ) )

    call dlg_check ( nf90_close ( nc_id ) )

    write(*,*) 'DONE'

end subroutine write_data

end module write_output
