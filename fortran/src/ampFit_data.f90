module ampFit_data
!
! subroutine ampFit_read_data
! subroutine ampFit_fill_structures
! subroutine ampFit_conv_GEI_GEO
! subroutine create_dataHalfSphere
!
! D.L. Green, C.L Waters
! Centre for Space Physics
! University of Newcastle
! Australia
! Nov 2011
!
! Ver : 201202
!
use constants
use netcdf
use ampFit_rotate

implicit none

real(kind=DBL), allocatable :: time(:), pseudo_sv_num(:), &
    plane_num(:), pos_eci(:,:), b_eci(:,:), &
    pseudo_sv_quality(:), data_splice(:)

contains

subroutine check_cdf( status )
  integer, intent(in) :: status

  if (status /= nf90_noerr) then
    print*, trim(nf90_strerror(status))
    err_stat = 1
  end if
  
end subroutine check_cdf
!
! ----------------------------------------------------------------------
!
subroutine ampFit_read_data (deltab_filename)
! Input:  deltab_filename
! Output: arrays of time,pseudo_sv_num,plane_num,pos_eci,b_eci, etc.

    implicit none

    character(len=100), intent(in) :: deltab_filename

    integer :: nSingular, dim_ids(2)
    integer :: nVec, nObs

    integer :: nc_id, nObs_id, nVec_id, time_id, pseudo_sv_num_id, &
        plane_num_id, pos_eci_id, b_eci_id, pseudo_sv_quality_id, &
        data_splice_id

    call check_cdf ( nf90_open ( path = deltab_fileName, mode = nf90_nowrite, ncid = nc_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'time', time_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'pseudo_sv_num', pseudo_sv_num_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'plane_num', plane_num_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'pos_eci', pos_eci_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'b_eci', b_eci_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'pseudo_sv_quality', pseudo_sv_quality_id ) )
    call check_cdf ( nf90_inq_varId ( nc_id, 'data_splice', data_splice_id ) )

    call check_cdf ( nf90_inquire_variable ( nc_id, pos_eci_id, dimIds = dim_ids ) )
    call check_cdf ( nf90_inquire_dimension ( nc_id, dim_ids(2), len = nObs ) )
    call check_cdf ( nf90_inquire_dimension ( nc_id, dim_ids(1), len = nVec ) )
    if (err_stat .eq. 1) return

!    write(*,*) '    Number of Pts: ', nObs
!    write(*,*) '    Number of Vector Components: ', nVec

    allocate ( time(nObs), &
        pseudo_sv_num(nObs), &
        plane_num(nObs), &
        pos_eci(nVec,nObs), &
        b_eci(nVec,nObs), &
        pseudo_sv_quality(nObs), &
        data_splice(nObs) )

    call check_cdf ( nf90_get_var ( nc_id, time_id, time ) )
    call check_cdf ( nf90_get_var ( nc_id, pseudo_sv_num_id, pseudo_sv_num ) )
    call check_cdf ( nf90_get_var ( nc_id, plane_num_id, plane_num ) )
    call check_cdf ( nf90_get_var ( nc_id, pos_eci_id, pos_eci ) )
    call check_cdf ( nf90_get_var ( nc_id, b_eci_id, b_eci ) )
    call check_cdf ( nf90_get_var ( nc_id, pseudo_sv_quality_id, pseudo_sv_quality ) )
    call check_cdf ( nf90_get_var ( nc_id, data_splice_id, data_splice ) )

    call check_cdf ( nf90_close ( nc_id ) )
    if (err_stat .eq. 1) return

end subroutine ampFit_read_data
!
! ----------------------------------------------------------------------------------
!
subroutine ampFit_fill_structures(deltab_filename, &
                  data_dir, base_name, &
                  sHr, eHr, sYr, sMn, sDy, extend, &
                  dataGEI_orig)

! Input: arrays of time,pseudo_sv_num,plane_num,pos_eci,b_eci, etc.
! Output: dataGEI_orig structure - see constants.f90

  use cnvtime

  implicit none

  character(len=100), intent(in) :: deltab_filename
  character(len=60), intent(in) :: data_dir
  character(len=20), intent(in)  :: base_name
  real(kind=DBL), intent(in) :: sHr, eHr
  integer, intent(in) :: sYr, sMn, sDy, extend
  type(ampData), intent(out), allocatable :: dataGEI_orig(:)

  character(len=100) :: exdb_name
  character(len=8) :: date_str
  character(len=4) :: yr_str, frmt
  integer :: i, j, nSubSet, sw, np, totnp, date_int, hr,mn,sc
  integer :: eYr, eMn, eDy
  integer, allocatable :: iiSubSet(:)

  integer(kind=DBL) :: ep_tme
  real(kind=DBL), allocatable :: t_time(:), e_time(:), &
                        t_pseudo_sv_num(:), e_pseudo_sv_num(:), &
                            t_plane_num(:), e_plane_num(:), &
                            t_pos_eci(:,:), e_pos_eci(:,:), &
                              t_b_eci(:,:), e_b_eci(:,:), &
                    t_pseudo_sv_quality(:), e_pseudo_sv_quality(:), &
                          t_data_splice(:), e_data_splice(:)

! Start of routine: read data file and fill structures for given time interval
!
!  read cdf datafile
  print*,'Extend = ',extend
  print*,'Filename = ',deltab_filename
  call ampFit_read_data (deltab_filename)
  np = size(time)
  if (extend.eq.0) then
    allocate ( t_time(np), t_pseudo_sv_num(np), &
               t_plane_num(np), t_pos_eci(3,np), &
               t_b_eci(3,np), t_pseudo_sv_quality(np), &
               t_data_splice(np) )
    t_time = time
    t_pseudo_sv_num = pseudo_sv_num
    t_plane_num = plane_num
    t_pos_eci = pos_eci
    t_b_eci = b_eci
    t_pseudo_sv_quality = pseudo_sv_quality
    t_data_splice = data_splice
  endif

  if (extend.eq.1) then
! store 1st data segment in e_ arrays
    allocate ( e_time(np), e_pseudo_sv_num(np), &
               e_plane_num(np), e_pos_eci(3,np), &
               e_b_eci(3,np), e_pseudo_sv_quality(np), &
               e_data_splice(np) )
    do i=1,np
      e_time(i) = time(i)
      e_pseudo_sv_num(i) = pseudo_sv_num(i)
      e_plane_num(i) = plane_num(i)
      do j = 1,3
        e_pos_eci(j,i) = pos_eci(j,i)
        e_b_eci(j,i) = b_eci(j,i)
      enddo
      e_pseudo_sv_quality(i) = pseudo_sv_quality(i)
      e_data_splice(i) = data_splice(i)
    enddo

! calc deltab_fname for next file
    ep_tme = cnv_YMDHMSToEpoch(sYr,sMn,sDy,1,0,0)
    ep_tme = ep_tme + 86400
    call cnv_EpochToYMDHMS(ep_tme,eYr,eMn,eDy,hr,mn,sc)
    frmt='i8'
    date_int = eYr*10000 + eMn*100 + eDy
    write(date_str,frmt) date_int   ! format using internal file
    frmt = '(i4)'
    write(yr_str,frmt) eYr
    exdb_name = trim(data_dir) // yr_str // '/' &
                // trim(date_str) // trim(base_name)
    
    call ampFit_read_data (exdb_name)
    totnp = np           ! num pnts of 1st day
    np = size(time)
    totnp = totnp + np
    allocate ( t_time(totnp), t_pseudo_sv_num(totnp), &
             t_plane_num(totnp), t_pos_eci(3,totnp), &
             t_b_eci(3,totnp), t_pseudo_sv_quality(totnp), &
             t_data_splice(totnp) )
    do i=1,totnp-np
      t_time(i) = e_time(i)
      t_pseudo_sv_num(i) = e_pseudo_sv_num(i)
      t_plane_num(i) = e_plane_num(i)
      do j=1,3
        t_pos_eci(j,i) = e_pos_eci(j,i)
        t_b_eci(j,i) = e_b_eci(j,i)
      enddo
      t_pseudo_sv_quality(i) = e_pseudo_sv_quality(i)
      t_data_splice(i) = e_data_splice(i)
    enddo

    do i = 1, np
      t_time(i+totnp-np) = time(i)
      t_pseudo_sv_num(i+totnp-np) = pseudo_sv_num(i)
      t_plane_num(i+totnp-np) = plane_num(i)
      do j=1,3
        t_pos_eci(j,i+totnp-np) = pos_eci(j,i)
        t_b_eci(j,i+totnp-np) = b_eci(j,i)
      enddo
      t_pseudo_sv_quality(i+totnp-np) = pseudo_sv_quality(i)
      t_data_splice(i+totnp-np) = data_splice(i)
    enddo 
  endif 

  nSubSet = count ( t_time >= sHr .and. t_time <= eHr )

  allocate ( iiSubSet(nSubSet), dataGEI_orig(nSubSet) )

  iiSubSet = pack ( (/ (i, i=1, size(t_time)) /), mask = (t_time >= sHr .and. t_time <= eHr) )  

  dataGEI_orig%utc = t_time(iiSubSet)
  dataGEI_orig%iSat = t_pseudo_sv_num(iiSubSet)
  dataGEI_orig%iPln = t_plane_num(iiSubSet)
  dataGEI_orig%qual = t_pseudo_sv_quality(iiSubSet)
  dataGEI_orig%splice = t_data_splice(iiSubSet)

  dataGEI_orig%x = t_pos_ECI(1,iiSubSet)*1d-3
  dataGEI_orig%y = t_pos_ECI(2,iiSubSet)*1d-3
  dataGEI_orig%z = t_pos_ECI(3,iiSubSet)*1d-3

  dataGEI_orig%bx = t_b_ECI(1,iiSubSet)
  dataGEI_orig%by = t_b_ECI(2,iiSubSet)
  dataGEI_orig%bz = t_b_ECI(3,iiSubSet)

  sw = -1
  call rtp_xyz_coord ( dataGEI_orig, sw)
  call xyz_to_rtp_vec ( dataGEI_orig )

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
          b_ECI )

  deallocate ( &
          t_time, &
          t_pseudo_sv_num, &
          t_plane_num, &
          t_pseudo_sv_quality, &
          t_data_splice, &
          t_pos_ECI, &
          t_b_ECI )

  if (extend.eq.1) then
    deallocate ( &
          e_time, &
          e_pseudo_sv_num, &
          e_plane_num, &
          e_pseudo_sv_quality, &
          e_data_splice, &
          e_pos_ECI, &
          e_b_ECI )
  endif

end subroutine ampFit_fill_structures
!
! ----------------------------------------------------------------------------------
!
subroutine ampFit_conv_GEI_GEO ( dataGEI_in, dataGEO_out )
!  assumes dataGEI structure has been populated

  implicit none

  type(ampData), intent(in) :: dataGEI_in(:)
  type(ampData), allocatable, intent(out) :: dataGEO_out(:)

  integer :: i, sw

! Initialise dataGEO_out. This also sets dB(r,tp) which are invariant here (CLW)
  allocate(dataGEO_out(size(dataGEI_in)))
  dataGEO_out = dataGEI_in

! GEI (x,y,z) -> GEO (x,y,z)
  sw = 1
  GEIxyz_to_GEOxyz: &
  do i=1,size(dataGEI_in) 
    call geigeo_08(dataGEI_in(i)%x,dataGEI_in(i)%y,dataGEI_in(i)%z, &
                dataGEO_out(i)%x, dataGEO_out(i)%y,dataGEO_out(i)%z, sw)

  end do GEIxyz_to_GEOxyz

! GEOxyz_to_GEOrtp
  sw = -1
  call rtp_xyz_coord ( dataGEO_out, sw)

! vec_rtp_to_GEOxyz
  call rtp_to_xyz_vec ( dataGEO_out )

end subroutine ampFit_conv_GEI_GEO
!
! -----------------------------------------------------------------------------
!
subroutine create_dataHalfSphere ( dataIn, dataHalf, south )

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampData), allocatable, intent(inout) :: dataHalf(:)
    integer, intent(in) :: south

    integer, allocatable :: iiSubSet(:)
    integer :: nHalfSphere, i
    logical, allocatable :: mask(:)

    allocate ( mask(size(dataIn) ) )

    if (south .eq. 1) then
      mask = dataIn%T * radToDeg >= 90
    else
      mask = dataIn%T * radToDeg <= 90
    end if

    nHalfSphere = count ( mask )
    if (nHalfSphere .lt. 10) then
      err_stat = 2
      err_msg='Insufficient data in create_dataHalfSphere'
      write(*,*) err_msg
      return
    endif

    allocate ( iiSubSet(nHalfSphere), dataHalf(nHalfSphere) )

    iiSubSet = pack ( (/ (i, i=1, size(dataIn)) /), mask = mask )  

    dataHalf  = dataIn(iiSubSet)

    deallocate ( mask, iiSubSet )

end subroutine create_dataHalfSphere

end module ampFit_data
!
! =============================================================================
