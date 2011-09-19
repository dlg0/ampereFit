module write_output

use netcdf
!use ampFit_solve
use dlg
use constants
 
contains

subroutine write_data ( dataIn, fileName )

    implicit none

    type(ampData) :: dataIn(:)

    character(len=*), intent(in), optional :: fileName
    character(len=100)  :: useFileName

    integer :: utc_id, bR_id, bT_id, bP_id, &
        nObs_id, nc_id, &
        T_id, P_id, R_id, &
        bX_id, bY_id, bZ_id, &
        X_id, Y_id, Z_id

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

    call dlg_check ( nf90_def_var ( nc_id, 'bR', NF90_DOUBLE, (/ nObs_id /), bR_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'bT', NF90_DOUBLE, (/ nObs_id /), bT_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'bP', NF90_DOUBLE, (/ nObs_id /), bP_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'bX', NF90_DOUBLE, (/ nObs_id /), bx_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'bY', NF90_DOUBLE, (/ nObs_id /), by_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'bZ', NF90_DOUBLE, (/ nObs_id /), bz_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'T', NF90_DOUBLE, (/ nObs_id /), T_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'P', NF90_DOUBLE, (/ nObs_id /), P_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'R', NF90_DOUBLE, (/ nObs_id /), R_id ) )

    call dlg_check ( nf90_def_var ( nc_id, 'X', NF90_DOUBLE, (/ nObs_id /), X_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'Y', NF90_DOUBLE, (/ nObs_id /), Y_id ) )
    call dlg_check ( nf90_def_var ( nc_id, 'Z', NF90_DOUBLE, (/ nObs_id /), Z_id ) )
  
    call dlg_check ( nf90_enddef ( nc_id ) )
    
    call dlg_check ( nf90_put_var ( nc_id, utc_id, dataIn%utc ) )

    call dlg_check ( nf90_put_var ( nc_id, bR_id, dataIn%bR) )
    call dlg_check ( nf90_put_var ( nc_id, bT_id, dataIn%bT) )
    call dlg_check ( nf90_put_var ( nc_id, bP_id, dataIn%bP) )

    call dlg_check ( nf90_put_var ( nc_id, bx_id, dataIn%bX ) )
    call dlg_check ( nf90_put_var ( nc_id, by_id, dataIn%bY ) )
    call dlg_check ( nf90_put_var ( nc_id, bz_id, dataIn%bZ ) )

    call dlg_check ( nf90_put_var ( nc_id, T_id, dataIn%T ) )
    call dlg_check ( nf90_put_var ( nc_id, P_id, dataIn%P ) )
    call dlg_check ( nf90_put_var ( nc_id, R_id, dataIn%R ) )

    call dlg_check ( nf90_put_var ( nc_id, X_id, dataIn%x ) )
    call dlg_check ( nf90_put_var ( nc_id, Y_id, dataIn%y ) )
    call dlg_check ( nf90_put_var ( nc_id, Z_id, dataIn%z ) )

    call dlg_check ( nf90_close ( nc_id ) )

    write(*,*) 'DONE'

end subroutine write_data

end module write_output
