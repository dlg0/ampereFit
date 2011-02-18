module write_output

use netcdf
use ampFit_data, &
only: dataOriginal, nObs=>nSubSet, nVec
use ampFit_solve
use dlg
 
contains

subroutine write_FAC ()

    implicit none

    character(len=100) :: outFileName

    integer :: utc_id, bR_GEI_id, bT_GEI_id, bP_GEI_id, &
        nObs_id, nVec_id, nc_id, &
        fit_bT_GEI_id, fit_bP_GEI_id, &
        data_GEI_coLat_rad_id, data_GEI_lon_rad_id

    outFileName = 'jParIrid.nc'
    
    write(*,*) 'Writing ', trim(outFileName), ' ...'

    call dlg_check ( nf90_create ( outFileName, nf90_clobber, nc_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nObs', nObs, nObs_id ) )
    call dlg_check ( nf90_def_dim ( nc_id, 'nVec', nVec, nVec_id ) )
   
    call dlg_check ( nf90_def_var ( nc_id, 'utc', NF90_DOUBLE, (/ nObs_id /), utc_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_bR_GEI', NF90_DOUBLE, (/ nObs_id /), bR_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_bT_GEI', NF90_DOUBLE, (/ nObs_id /), bT_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_bP_GEI', NF90_DOUBLE, (/ nObs_id /), bP_GEI_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'fit_bT_GEI', NF90_DOUBLE, (/ nObs_id /), fit_bT_GEI_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'fit_bP_GEI', NF90_DOUBLE, (/ nObs_id /), fit_bP_GEI_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_coLat_rad', &
            NF90_DOUBLE, (/ nObs_id /), data_GEI_coLat_rad_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'data_GEI_lon_rad', &
            NF90_DOUBLE, (/ nObs_id /), data_GEI_lon_rad_id ) )
   
    call dlg_check ( nf90_enddef ( nc_id ) )
    
    call dlg_check ( nf90_put_var ( nc_id, utc_id, dataOriginal%utc ) )

    call dlg_check ( nf90_put_var ( nc_id, bR_GEI_id, dataOriginal%bR_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, bT_GEI_id, dataOriginal%bTheta_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, bP_GEI_id, dataOriginal%bPhi_GEI ) )

    call dlg_check ( nf90_put_var ( nc_id, fit_bT_GEI_id, fit_bTheta_GEI ) )
    call dlg_check ( nf90_put_var ( nc_id, fit_bP_GEI_id, fit_bPhi_GEI ) )

    call dlg_check ( nf90_put_var ( nc_id, data_GEI_coLat_rad_id, dataOriginal%GEI_coLat_rad ) )
    call dlg_check ( nf90_put_var ( nc_id, data_GEI_lon_rad_id, dataOriginal%GEI_lon_rad ) )
 
    call dlg_check ( nf90_close ( nc_id ) )

    write(*,*) 'DONE'

end subroutine write_fac

end module write_output
