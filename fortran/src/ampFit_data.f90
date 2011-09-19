module ampFit_data

use constants
use ampFit_namelist
use netcdf
use dlg

implicit none

! Populate data structure variables

type(ampData), allocatable :: dataOriginal(:), dataShifted(:), dataHalfSphere(:)

real(kind=DBL), allocatable :: time(:), pseudo_sv_num(:), &
    plane_num(:), pos_eci(:,:), b_eci(:,:), &
    pseudo_sv_quality(:), data_splice(:)

contains

subroutine ampFit_read_data

    implicit none

    integer :: nSingular, dim_ids(2)
    integer :: nVec, nObs

    integer :: nc_id, nObs_id, nVec_id, time_id, pseudo_sv_num_id, &
        plane_num_id, pos_eci_id, b_eci_id, pseudo_sv_quality_id, &
        data_splice_id

    write(*,*) 'Reading ', trim ( deltab_fileName ), ' ...'

    call dlg_check ( nf90_open ( path = deltab_fileName, mode = nf90_nowrite, ncid = nc_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'time', time_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'pseudo_sv_num', pseudo_sv_num_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'plane_num', plane_num_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'pos_eci', pos_eci_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'b_eci', b_eci_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'pseudo_sv_quality', pseudo_sv_quality_id ) )
    call dlg_check ( nf90_inq_varId ( nc_id, 'data_splice', data_splice_id ) )

    call dlg_check ( nf90_inquire_variable ( nc_id, pos_eci_id, dimIds = dim_ids ) )
    call dlg_check ( nf90_inquire_dimension ( nc_id, dim_ids(2), len = nObs ) )
    call dlg_check ( nf90_inquire_dimension ( nc_id, dim_ids(1), len = nVec ) )

    write(*,*) '    Number of Pts: ', nObs
    write(*,*) '    Number of Vector Components: ', nVec

    allocate ( time(nObs), &
        pseudo_sv_num(nObs), &
        plane_num(nObs), &
        pos_eci(nVec,nObs), &
        b_eci(nVec,nObs), &
        pseudo_sv_quality(nObs), &
        data_splice(nObs) )

    call dlg_check ( nf90_get_var ( nc_id, time_id, time ) )
    call dlg_check ( nf90_get_var ( nc_id, pseudo_sv_num_id, pseudo_sv_num ) )
    call dlg_check ( nf90_get_var ( nc_id, plane_num_id, plane_num ) )
    call dlg_check ( nf90_get_var ( nc_id, pos_eci_id, pos_eci ) )
    call dlg_check ( nf90_get_var ( nc_id, b_eci_id, b_eci ) )
    call dlg_check ( nf90_get_var ( nc_id, pseudo_sv_quality_id, pseudo_sv_quality ) )
    call dlg_check ( nf90_get_var ( nc_id, data_splice_id, data_splice ) )

    call dlg_check ( nf90_close ( nc_id ) )

    write (*,*) 'DONE'

end subroutine ampFit_read_data


subroutine ampFit_fill_structures

    implicit none

    integer :: i, j, nSubSet
    integer, allocatable :: iiSubSet(:)

    write(*,*) 'Populating dataOriginal structure ...'

    nSubSet = count ( time >= sHr .and. time <= eHr )

    write(*,*) '    sHr: ', sHr
    write(*,*) '    eHr: ', eHr
    write(*,*) '    nSubSet: ', nSubSet

    allocate ( iiSubSet(nSubSet), dataOriginal(nSubSet), dataShifted(nSubSet) )

    iiSubSet = pack ( (/ (i, i=1, size(time)) /), mask = (time >= sHr .and. time <= eHr) )  

    dataOriginal%utc = time(iiSubSet)
    dataOriginal%iSat = pseudo_sv_num(iiSubSet)
    dataOriginal%iPln = plane_num(iiSubSet)
    dataOriginal%qual = pseudo_sv_quality(iiSubSet)
    dataOriginal%splice = data_splice(iiSubSet)

    dataOriginal%x = pos_ECI(1,iiSubSet)*1d-3
    dataOriginal%y = pos_ECI(2,iiSubSet)*1d-3
    dataOriginal%z = pos_ECI(3,iiSubSet)*1d-3

    dataOriginal%bx = b_ECI(1,iiSubSet)
    dataOriginal%by = b_ECI(2,iiSubSet)
    dataOriginal%bz = b_ECI(3,iiSubSet)


    ! Deallocate to full set of data, 
    ! i.e., keep only the subSet available 
    ! ------------------------------------

    deallocate ( &
            time, &
            pseudo_sv_num, &
            plane_num, &
            pseudo_sv_quality, &
            data_splice, &
            pos_ECI, &
            b_eCI )

    write(*,*) 'DONE'

end subroutine ampFit_fill_structures

subroutine create_dataHalfSphere ( dataIn )

    implicit none

    type(ampData), intent(in) :: dataIn(:)

    integer, allocatable :: iiSubSet(:)
    integer :: nHalfSphere, i
    logical, allocatable :: mask(:)

    write(*,*) 'Creating dataHalfSphere ...'

    allocate ( mask(size(dataIn) ) )

    mask    = dataIn%T * radToDeg <= 90

    nHalfSphere = count ( mask )

    write(*,*) '    nHalfSphere: ', nHalfSphere

    allocate ( iiSubSet(nHalfSphere), dataHalfSphere(nHalfSphere) )

    iiSubSet = pack ( (/ (i, i=1, size(dataIn)) /), mask = mask )  

    dataHalfSphere  = dataIn(iiSubSet)

    deallocate ( mask, iiSubSet )

end subroutine create_dataHalfSphere

end module ampFit_data

