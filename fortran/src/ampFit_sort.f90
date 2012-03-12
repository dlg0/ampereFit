module ampFit_sort
! This module contains:
! 1. do_merge -> subroutine for merge_sort
! 2. merge_sort -> sort routine adapted from Wiki page on merge sort
!                  This routine has the sorted data plus the indexes
! 3. cross_p -> vector cross product
! 4. norm_vec -> normailise 3D vectors
! 5. gcirc_dist -> great circle distance (over a sphere) routine
! 6. my_where -> similar to IDL where function for 1D arrays
! 7. sort_struc -> sorts satellite track data points by great circle distance
! 8. tag_lon_strays -> finds stray satellites (off longitude track)
! 9. ampFit_ghost_aacgm -> adds ghost data near poles
!10. pop_ghost_struc -> used by ampFit_ghost_aacgm
!
! C. L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
! Aug 2011
!
! Last modified: 
!   Feb 2012: constrained ghosts with lon in GEO - CLW
!
! Ver 201202
!
  use constants
  use ISO_C_BINDING
  use f_aacgm

  contains

  subroutine do_merge(A,na,B,nb,C,nc,Ai,Bi,Ci)

    implicit none
  
    integer, intent(in) :: na, nb, nc
    real(kind=DBL), intent(in out) :: A(na), C(nc)
    real(kind=DBL), intent(in) :: B(nb)
    integer, intent(in) :: Bi(nb)
    integer, intent(in out) :: Ai(na), Ci(nc)
    integer :: i, j, k
    integer :: ii

    i=1; j=1; k=1;
    do while (i <= na .and. j <= nb)
      if (A(i) <= B(j)) then
        C(k)=A(i)
        Ci(k)=Ai(i)
        i=i+1
      else
        C(k)=B(j)
        Ci(k)=Bi(j)
        j=j+1
      end if
      k=k+1
    enddo
    do while (i <= na)
      C(k)=A(i)
      Ci(k)=Ai(i)
      k=k+1
      i=i+1
    enddo
    return
  end subroutine do_merge

!------------------------------------------------------------------
  
  recursive subroutine merge_sort(A,N,T,Ai,Ti)
! merge sort with recursive call
! calls do_merge subroutine
!   C.L. Waters:  22 Aug 2011
! The arrays A and Ai are replaced by the sorted results
! A : the input array of real numbers
! N : number of elements
! T : temp array for the sort
! Ai : array of indexes of A
! Ti : array of indexes for T
!
    implicit none

    integer, intent(in) :: N
    real(kind=DBL), dimension(N), intent(in out) :: A
    integer, dimension(N), intent(in out) :: Ai
    real(kind=DBL), dimension((N+1)/2), intent(out) :: T
    integer, dimension((N+1)/2), intent(out) :: Ti

    integer :: na, nb, val_i
    integer :: ii
    real(kind=DBL) :: val

    if (N < 2) return      ! only 1 value in in_arr

    if (N==2) then         ! 2 values in in_arr
      if (A(1) > A(2)) then
        val=A(1)
        val_i=Ai(1)
        A(1) = A(2)
        Ai(1)= Ai(2)
        A(2)=val
        Ai(2)=val_i
      endif
      return
    endif
    na=(N+1)/2
    nb=N-na

    call merge_sort(A,na,T,Ai,Ti)
    call merge_sort(A(na+1),nb,T,Ai(na+1),Ti)
    if (A(na) > A(na+1)) then
      T(1:na)=A(1:na)
      Ti(1:na)=Ai(1:na)

      call do_merge(T,na,A(na+1),nb,A,N,Ti,Ai(na+1),Ai)
    endif
    return

  end subroutine merge_sort

!--------------------------------------------------------------

  function cross_p(vec1, vec2)
! calc the cross product of two 3D vectors
    implicit none

    real(kind=DBL), dimension(3) :: cross_p
    real(kind=DBL), dimension(3), intent(in) :: vec1, vec2

    cross_p(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_p(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_p(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

  end function cross_p

!--------------------------------------------------------------

  function norm_vec(vecin)
! normalise a 3D vector: division by its magnitude
    implicit none

    real(kind=DBL), dimension(3) :: norm_vec
    real(kind=DBL), dimension(3), intent(in) :: vecin
    real(kind=DBL) :: vecmag

    vecmag = sqrt(vecin(1)**2 + vecin(2)**2 + vecin(3)**2)
    norm_vec = vecin/vecmag

  end function norm_vec

! -------------------------------------------------------------

  subroutine gcirc_dist(lat_s, lat_f, lon_s, lon_f, &
                        np, gcirc_res)
! Calc great circle distance over a sphere given:
! lat_s -> start latitude (deg)
! lon_s -> start longitude (deg)
! lat_f -> array of latitudes to calc distances
! lon_f -> array of corresponding longitudes
! np -> number of points in arrays: lat_f, lon_f
! gcirc_res -> results array of great circle distances (km) 
! uses rsat and rE in constants.f90
!
! C.L. Waters - Aug 2011
!
    implicit none

    integer, intent(in) :: np
    real(kind=DBL), intent(in) :: lat_s, lon_s
    real(kind=DBL), dimension(np), intent(in) :: lat_f, lon_f
    real(kind=DBL) :: phi_s, lmda_s
    real(kind=DBL), dimension(np) :: phi_f, lmda_f
    real(kind=DBL), dimension(np), intent(out) :: gcirc_res
    real(kind=DBL), dimension(np) :: d_lmda, tp, bt

    phi_s = lat_s*pi/180.0
    phi_f = lat_f*pi/180.0
    lmda_s= lon_s*pi/180.0
    lmda_f= lon_f*pi/180.0
    d_lmda= lmda_f - lmda_s

    tp = sqrt( ( cos(phi_f)*sin(d_lmda) )**2 + &
        ( cos(phi_s)*sin(phi_f) - sin(phi_s)*cos(phi_f)*cos(d_lmda))**2)
    bt = sin(phi_s)*sin(phi_f) + cos(phi_s)*cos(phi_f)*cos(d_lmda)
    gcirc_res = (rE + rsat)/1000.0 * atan2(tp,bt) 
  end subroutine gcirc_dist

!---------------------------------------------------------------

  subroutine my_where(logic_arr, np, idx_arr, ntrue)
! Subroutine to find 1D array indexes similar to IDL 'where' function
! inputs are 
!   logic_arr :   array of T or F for given logic statement
!   np :  size of logic_array
!
! output is index array of all the T and how many (ntrue)
!
! C.L. Waters :  Aug, 2011
!
    implicit none

    integer, intent(in) :: np
    integer :: i
    logical, dimension(np), intent(in) :: logic_arr
    integer, dimension(np), intent(inout) :: idx_arr
    integer, dimension(np) :: idx_tmp
    integer, intent(inout) :: ntrue

    do i=1,np
      idx_tmp(i) = i
    enddo
    ntrue = count (logic_arr)
    idx_arr = pack ( (/ (idx_tmp(i), i=1,np) /), mask = logic_arr)

  end subroutine my_where

! --------------------------------------------------------------

  subroutine sort_struc(struc_data, trk_order, south)
! Given the shiftedData structure, this routine:
! (i)   Cycles through each Track Number/hemisphere
! (ii)  Identify and label any strays (off longitude)
! (iii) Sort by great circle distance from an equator point
!  (iv) ensures the tracks are sequential in longitude
!
! Output is the modified shiftedData structure that has
! - Strays tagged in the typ property, where
!    data%typ = 0 -> normal data
!    data%typ = 1 -> stray data off longitude
!    data%typ = 2 -> ghost data points near the poles
!    data%typ = 3 -> additional ghosts where longitude separation is > 120 deg
!
! The resulting data structure (sorted): 
!  and track order by longitude, ready for ghost data calc
!
! C.L. Waters - Aug 2011
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Ver : 201202
!
    implicit none

    type(ampData), intent(inout) :: struc_data(:)
    integer, intent(inout) :: trk_order(6)
    integer, intent(in) :: south

    logical, dimension(1:size(struc_data)) :: tmplogic_arr
    integer, dimension(1:size(struc_data)) :: tmp_idx, sortII
    integer :: cnt, i, np, Tr_num, iiTrack
    integer :: iiTrack_lon
    integer, dimension(:), allocatable :: idx, idx_lon, idx_srt
    real(kind=DBL), dimension(:), allocatable :: gcDist
    real(kind=DBL) :: st_lat, st_lon
    real(kind=DBL), dimension(6) :: strt_lona
    real(kind=DBL), dimension(3) :: uvec, vvec, nvec, ncuvec
    integer, dimension(1) :: iiSt
! temp arrays for sort routine
    real(kind=DBL), dimension(:), allocatable :: T
    integer, dimension(:), allocatable :: Ti

    np = size(struc_data)
    cnt = 1

!print*,'Sorting each data track...'
    sort_track_loop: &
    do Tr_num = 0,5

! Get indexes where data has correct TrackNum and hemisphere
      tmplogic_arr = .false.          ! initialise logic array
      if (south .eq. 0) then
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg < 90.0
      else
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg > 90.0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack)
      allocate(idx(iiTrack))      ! get mem for these data
      idx = tmp_idx(1:iiTrack)    ! indexes for these vals

      if (iiTrack .lt. 20) then
        err_stat = 5
        err_msg = 'Insufficient data on track'
        write(*,*) err_msg
        print*,'Track Num:',Tr_num
        return
      endif

! Get indexes for (i) Track (ii) < 180 in lon 
      tmplogic_arr = .false.
      tmplogic_arr=struc_data(idx)%P*radToDeg <= 180.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_lon)
      allocate(idx_lon(iiTrack_lon))    ! get mem for these data
      idx_lon = tmp_idx(1:iiTrack_lon)  ! indexes for these vals
! These are the indexes of struc for correct hemis, Tr_num and lon <= 180

      if (south .eq. 0) then
        iiSt = maxloc(struc_data(idx(idx_lon))%T*radToDeg)  ! index of equator point
      else
        iiSt = minloc(struc_data(idx(idx_lon))%T*radToDeg)  ! index of equator point
      endif
      i = idx(idx_lon(iiSt(1)))
      st_lat = 90.0 - struc_data(i)%T*radToDeg  ! Lat of this point
      st_lon  = struc_data(i)%P*radToDeg          ! Lon of this point
! st_Lat,st_Lon is the coord of the great circle dist reference point

      strt_lona(Tr_num+1) = st_lon         ! store start longitude
      trk_order(Tr_num+1) = Tr_num
      allocate(gcDist(iiTrack))            ! iiTrack is the size of idx
      call gcirc_dist(st_lat, 90.0-struc_data(idx)%T*radToDeg, &
             st_lon, struc_data(idx)%P*radToDeg, iiTrack, gcDist)

! sort the distances, returns idx of sorted array
      allocate(T((iiTrack+1)/2))
      allocate(Ti((iiTrack+1)/2))
      allocate(idx_srt(iiTrack))
      idx_srt = idx              ! sort these indexes
      call merge_sort(gcDist, iiTrack, T, idx_srt, Ti)
      deallocate(Ti)
      deallocate(T)
      deallocate(gcDist)

      sortII(cnt:cnt+iiTrack-1) = idx_srt  ! save sorted indexes
      cnt = cnt + iiTrack
 
      deallocate(idx_lon)
      deallocate(idx)
      deallocate(idx_srt)

    enddo sort_track_loop
    struc_data = struc_data( sortII )   ! sort the structure

! sort the track order, returns trk_order
    allocate(T(6))
    allocate(Ti(6))
    call merge_sort(strt_lona, 6, T, trk_order, Ti)
    deallocate(Ti)
    deallocate(T)

  end subroutine sort_struc

! ----------------------------------------------------------------

  subroutine tag_lon_strays(struc_data, south)
! Given the shiftedData structure, and south=0 or 1, this routine:
! (i)   Cycles through each Track Number/hemis
! (ii)  Identify and label any strays (off longitude)
! Assumes the data for each track has been sorted by sort_struc
! - Tr0(Nth,Sth), Tr1(Nth,Sth) etc sequence sorted
! Output has the struc_data%typ flag modified
!
! C.L. Waters - Aug 2011
!
    implicit none
    type(ampData), intent(in out) :: struc_data(:)
    integer, intent(in) :: south

    logical, dimension(1:size(struc_data)) :: tmplogic_arr
    integer, dimension(1:size(struc_data)) :: tmp_idx
    integer :: i, np, Tr_num, iiTrack, iiSel_lat
    integer, dimension(1) :: iiMin_Loc
    integer, dimension(:), allocatable :: idx_n, idx_sel_lat
    real(kind=DBL), dimension(3) :: uvec, vvec, nvec, ncuvec
    real(kind=DBL), dimension(:), allocatable :: xarr, yarr, zarr, rarr, tarr
    real(kind=DBL), dimension(:), allocatable :: rVarr, LatVarr, LonVarr
    real(kind=DBL) :: dp, lon_diff

!    print*,'Tagging longitude strays...'
    np = size(struc_data)

    track_cyc_loop: &
    do Tr_num = 0,5

! select data from correct track and Nth hemis
      tmplogic_arr = .false.          ! initialise logic array
      if (south .eq. 1) then
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg >= 90.0
      else
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg <= 90.0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack)
      allocate(idx_n(iiTrack))      ! get mem for these data
      idx_n = tmp_idx(1:iiTrack)    ! indexes for these vals
      struc_data(idx_n)%typ = 0     ! initialise type as normal

! 1st point is on the equator - assumes data are sorted
      uvec(1) = struc_data(idx_n(1))%X
      uvec(2) = struc_data(idx_n(1))%Y
      uvec(3) = struc_data(idx_n(1))%Z
      uvec = norm_vec(uvec)

! Get data subset 1/2 about between pole and equator for 2nd point
      tmplogic_arr = .false.
      if (south .eq. 1) then
        tmplogic_arr=struc_data(idx_n)%T*radToDeg > 130.0 &
               .and. struc_data(idx_n)%T*radToDeg < 140.0
      else
        tmplogic_arr=struc_data(idx_n)%T*radToDeg > 40.0 &
               .and. struc_data(idx_n)%T*radToDeg < 50.0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiSel_lat)
      allocate(idx_sel_lat(iiSel_lat))    ! get mem for these data
      idx_sel_lat = tmp_idx(1:iiSel_lat)  ! indexes for these vals

      iiMin_loc = minloc ( &
                  abs( struc_data(idx_n(idx_sel_lat))%P*radToDeg &
                    -  struc_data(idx_n(1))%P*radToDeg ) )
! point 2 for great circle eqn
      vvec(1) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%X
      vvec(2) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%Y
      vvec(3) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%Z
      vvec = norm_vec(vvec)

      deallocate(idx_sel_lat)
! normal vector to uvec and vvec
      nvec = cross_p(uvec, vvec)
      nvec = norm_vec(nvec)

! normal vector cross uvec
      ncuvec = cross_p(nvec, uvec)
      
      allocate(rarr(iiTrack))
      rarr = sqrt(struc_data(idx_n)%X**2 + &
                  struc_data(idx_n)%Y**2 + &
                  struc_data(idx_n)%Z**2)
      allocate(tarr(iiTrack))
      tarr(1)=0.0d0
      do i=2,iiTrack
        dp = uvec(1)*struc_data(idx_n(i))%X + &
             uvec(2)*struc_data(idx_n(i))%Y + &
             uvec(3)*struc_data(idx_n(i))%Z
        tarr(i) = acos(dp/rarr(i))
      enddo

! get eqn of great circle curve in x,y,z
      allocate(xarr(iiTrack))
      allocate(yarr(iiTrack))
      allocate(zarr(iiTrack))
      xarr = rarr*cos(tarr)*uvec(1) + rarr*sin(tarr)*ncuvec(1)
      yarr = rarr*cos(tarr)*uvec(2) + rarr*sin(tarr)*ncuvec(2)
      zarr = rarr*cos(tarr)*uvec(3) + rarr*sin(tarr)*ncuvec(3)

! Conv from x,y,z to r,th,ph
      allocate(rVarr(iiTrack))
      allocate(LatVarr(iiTrack))
      allocate(LonVarr(iiTrack))

! Convert x,y,z to r,the,ph coords
      do i=1,iiTrack
        call sphcar_08 (rVarr(i), LatVarr(i), LonVarr(i), & 
                        xarr(i), yarr(i), zarr(i), -1)
      enddo
      LonVarr = LonVarr*180.0d0/pi ! conv to degrees
      deallocate(LatVarr)
      deallocate(rVarr)
      deallocate(zarr)
      deallocate(yarr)
      deallocate(xarr)
      deallocate(tarr)
      deallocate(rarr)

! ------ Check Lon of great circle eqn with the data for strays -----

! check for lon (from eqn) - lon (from data) is > 180 (across 360 deg)
! can't compare arrays with different indexes, so loop'
      do i=1,iiTrack
        lon_diff = abs(LonVarr(i) - struc_data(idx_n(i))%P*radToDeg)
        if (lon_diff > 180.0) then
! if we have lon diff > 180, need to subtract off 360
          if (((360.0 - lon_diff) > 15.0 ) .or. &
             ( lon_diff > 15.0 )) then
            struc_data(idx_n(i))%typ = 1
          endif
        else
          if ( lon_diff > 15.0 ) then
            struc_data(idx_n(i))%typ = 1
          endif
        endif
      enddo     ! search track data points

      deallocate(LonVarr)
      deallocate(idx_n)

    enddo track_cyc_loop
!    print*,'DONE'

  end subroutine tag_lon_strays
!
! ----------------------------------------------------------------------
!
  subroutine ampFit_ghosts_aacgm(struc_data, comb_data, trk_order, south)

! The southern hemisphere fit suffers from satellite intersection point i
!  offset in the Birkeland current symmetry location.  Solve this by adding
!  'ghost' data.
! Ghost the full 1/2 hemisphere tracks and use AACGM Lat coords to choose
!  which dB values to average
! Input data are assumed to be from ONE hemisphere (Nth or Sth)
! Ghosting is done for all points passed in
!
! Given:
!   shiftedData structure,
!   trk_order (integer array of 6 of the track order)
!
! Assumes the data for each track have been sorted by sort_struc
! - Tr0(Nth,Sth), Tr1(Nth,Sth) etc sequence sorted
!
! C. L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
! Oct 2011
!
! Last modified:
!  Feb 2012 - ghost using aacgm lat but GEO lon - CLW
!
! Ver : 201202
!
    implicit none

    type(ampData), intent(in) :: struc_data(:)
    type(ampData), allocatable, intent(out) :: comb_data(:)
    integer, dimension(6), intent(in) :: trk_order
    integer, intent(in) :: south

    type(ampData), allocatable :: gh_data(:)
    type(ampData) :: exgh_data, gh_struc, Trk1_p, Trk2_p

    logical, dimension(1:size(struc_data)) :: tmplogic_arr

    integer, dimension(1:size(struc_data)) :: tmp_idx
    integer, dimension(:), allocatable :: idx_Trk1, idx_Trk2

    integer :: np, Tr_num, iiTrk_1, iiTrk_2, nxt_trk, iiex_gh
    integer :: ii, jj, Nyq_Ok, cnt, iighost, ii_c, do_exgh
    integer :: tot_rec, tmp_open, s, flg, iij, Tr1_i, Tr2_i
    integer, dimension(1) :: iimn
    real(kind=DBL), dimension(:), allocatable :: gcDist

    real(kind=DBL), dimension(:), allocatable :: Tr1_GEO_coLat_deg, &
        Tr1_GEO_Lon_deg, Tr2_GEO_coLat_deg, Tr2_GEO_Lon_deg
    real(kind=DBL), dimension(:), allocatable :: Tr1_aacgm_coLat_deg, &
        Tr1_aacgm_Lon_deg, Tr2_aacgm_coLat_deg, Tr2_aacgm_Lon_deg
    real(kind=DBL), dimension(:), allocatable :: Tr1_coLat_deg, &
        Tr1_Lon_deg, Tr2_coLat_deg, Tr2_Lon_deg

    real(kind=DBL), target :: aacgm_clat_d, aacgm_lon_d, err
    real(kind=DBL) :: st_clat, st_lon, e_clat, e_lon, xs, ys, xe, ye
    real(kind=DBL) :: xv, yv, zv, coLat_wdth, hgt
    character(len=20) :: gh_fname

    integer :: i       ! for testing

    np = size(struc_data)
    iiex_gh = 0
    coLat_wdth = 2.0
    cnt = 1
    gh_fname='ex_ghosts.dat'
    tmp_open = 0
    hgt = rSat/1000.0d0    ! altitude in km for aacgm

!    print*,'Calc aacgm Ghost data...'

! Calc min number of ghost points (excludes extra ghosts)
    tmplogic_arr = .false.            ! initialise logic array
    tmplogic_arr = struc_data%typ==0  ! exclude lon strays
    call my_where(tmplogic_arr, np, tmp_idx, iighost)
    allocate(gh_data(iighost))        ! get mem for ghost_data struc

    track_pair_loop: &
    do Tr_num = 0,5
! mask for correct Tr_num and hemis
      tmplogic_arr = .false.        ! initialise logic array
      tmplogic_arr = struc_data%ipln==trk_order(Tr_num+1) .and. &
                     struc_data%typ==0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_1)
      allocate(idx_Trk1(iiTrk_1))   ! get mem for indexes for Trk_1
      idx_Trk1 = tmp_idx(1:iiTrk_1) ! indexes for these vals

      if (iiTrk_1 .lt. 20) then
        err_stat = 6
        err_msg = 'Insufficient data on Track in ghost routine'
        write(*,*) err_msg
        print*,'Track Num:',Tr_num
        return
      endif

      if (Tr_num .lt. 5) then       ! Track numbers 0..5
        nxt_trk=Tr_num+2            ! index numbers 1..6
      else
         nxt_trk=1                  ! fold back to 1st track in sorted sequence
      endif

      tmplogic_arr = .false.        ! initialise logic array
! Get data for adjacent satellite track -> Trk_2
      tmplogic_arr = struc_data%ipln==trk_order(nxt_trk) .and. &
                     struc_data%typ==0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_2)
      allocate(idx_Trk2(iiTrk_2))   ! get mem for indexes for Trk_2 data
      idx_Trk2 = tmp_idx(1:iiTrk_2) ! indexes for these vals

      if (iiTrk_2 .lt. 20) then
        err_stat = 6
        err_msg = 'Insufficient data on Track in ghost routine'
        write(*,*) err_msg
        print*,'Track Num:',Tr_num
        return
      endif

! Take these coords to AACGM
      allocate(Tr1_GEO_coLat_deg(iiTrk_1))
      Tr1_GEO_coLat_deg = struc_data(idx_Trk1)%T*radToDeg
      allocate(Tr1_GEO_Lon_deg(iiTrk_1))
      Tr1_GEO_Lon_deg = struc_data(idx_Trk1)%P*radToDeg

      allocate(Tr1_aacgm_coLat_deg(iiTrk_1))
      allocate(Tr1_aacgm_Lon_deg(iiTrk_1))

      flg = 0            ! 0 -> to aacgm, 1 -> geographic
      Tr1_aacgm_cnv: &
      do jj=1,iiTrk_1
        st_clat = Tr1_GEO_coLat_deg(jj)
        st_lon = Tr1_GEO_Lon_deg(jj)

        s = f_AACGMConvert(90.0-st_clat, st_lon, hgt, &
                            C_LOC(aacgm_cLat_d), &
                            C_LOC(aacgm_lon_d), &
                            C_LOC(err), flg)
        Tr1_aacgm_Lon_deg(jj) = aacgm_lon_d
        Tr1_aacgm_coLat_deg(jj) = 90.0-aacgm_clat_d
 
      enddo Tr1_aacgm_cnv     ! Track_1 th,ph conv to AACGM

      allocate(Tr2_GEO_coLat_deg(iiTrk_2))
      Tr2_GEO_coLat_deg = struc_data(idx_Trk2)%T*radToDeg
      allocate(Tr2_GEO_Lon_deg(iiTrk_2))
      Tr2_GEO_Lon_deg = struc_data(idx_Trk2)%P*radToDeg

      allocate(Tr2_aacgm_coLat_deg(iiTrk_2))
      allocate(Tr2_aacgm_Lon_deg(iiTrk_2))

      allocate(Tr1_coLat_deg(iiTrk_1))
      allocate(Tr1_Lon_deg(iiTrk_1))
      allocate(Tr2_coLat_deg(iiTrk_2))
      allocate(Tr2_Lon_deg(iiTrk_2))

      flg = 0            ! 0 -> to aacgm, 1 -> geographic
      Tr2_aacgm_cnv: &
      do jj=1,iiTrk_2
        st_clat = Tr2_GEO_coLat_deg(jj)
        st_lon = Tr2_GEO_Lon_deg(jj)

        s = f_AACGMConvert(90.0-st_clat, st_lon, hgt, &
                            C_LOC(aacgm_cLat_d), &
                            C_LOC(aacgm_lon_d), &
                            C_LOC(err), flg)
        Tr2_aacgm_Lon_deg(jj) = aacgm_lon_d
        Tr2_aacgm_coLat_deg(jj) = 90.0-aacgm_clat_d

      enddo Tr2_aacgm_cnv    ! Track_2 th,ph conv to AACGM

      allocate(gcDist(iiTrk_2))     ! mem for great circ dist values

! Loop thru all data point of Track_1
      point_ex_loop: &
      do jj=1,iiTrk_1
        st_clat = Tr1_GEO_coLat_deg(jj)
        st_lon = Tr1_GEO_Lon_deg(jj)
        call gcirc_dist(90.0-st_clat, 90.0-Tr2_GEO_coLat_deg, &
             st_lon, Tr2_GEO_Lon_deg, iiTrk_2, gcDist)

! iimn is the index of idx_Trk2 which has min dist
!  between Trk2 and the jj point of Trk1
        iimn = minloc(gcDist)        ! index of idx_Trk2 of min_dist

        if (south .eq. 1) st_clat=180.0-st_clat           ! Put into Nth hemis for XY calc
        xs = st_clat*cos(st_lon*pi/180.0)
        ys = st_clat*sin(st_lon*pi/180.0)
        e_clat = Tr2_GEO_coLat_deg(iimn(1))
        if (south .eq. 1) e_clat=180.0-e_clat
        e_lon = Tr2_GEO_Lon_deg(iimn(1))
        Nyq_Ok = 1

! Check if we need extra ghosts -> need to add 2 ghost points
        if ( abs(e_lon-st_lon) .ge. 120.0 .and. &
             abs(e_lon-st_lon) .lt. 180.0) then
          Nyq_Ok = 0
          iiex_gh = iiex_gh + 1
          if (tmp_open .eq. 0) then   ! 1st extra ghost
! open the extra ghost temp data file here
            open(unit=1,file=gh_fname,status='replace', &
                 action='write')
            tmp_open = 1   ! set file open switch
          endif
        endif

!  Calc the x,y coord of the Track_2 point
        xe = e_clat*cos(e_lon*pi/180.0)
        ye = e_clat*sin(e_lon*pi/180.0)

! Put ghost data into ghost data structure
        Tr1_i = jj
        Tr2_i = iimn(1)
        Trk1_p = struc_data(idx_Trk1(Tr1_i))
        Trk2_p = struc_data(idx_Trk2(Tr2_i))

! select GEO or AACGM coords depending on latitude
        If ( ((Trk1_p%T*radToDeg .lt. 70.0) .and. (Trk2_p%T*radToDeg .lt. 70.0)) .or. &
            ((Trk1_p%T*radToDeg .gt. 110.0) .and. (Trk2_p%T*radToDeg .gt. 110.0)) ) then
          Tr1_coLat_deg = Tr1_aacgm_coLat_deg
          Tr2_coLat_deg = Tr2_aacgm_coLat_deg
        else
          Tr1_coLat_deg = Tr1_GEO_coLat_deg
          Tr2_coLat_deg = Tr2_GEO_coLat_deg
        endif
! MOD CLW - 21 Feb 2012, work lon in GEO
        Tr1_Lon_deg = Tr1_GEO_Lon_deg
        Tr2_Lon_deg = Tr2_GEO_Lon_deg

        if (Nyq_Ok .eq. 0) then  ! if we need 2 extra ghosts
          do_exgh = 1
          call pop_ghost_struc(Trk1_p,Trk2_p, &
                   Tr1_coLat_deg, Tr1_Lon_deg, &
                   Tr2_coLat_deg, Tr2_Lon_deg, &
                   coLat_wdth, xs,ys,xe,ye,south,do_exgh, &
                   Tr1_i,Tr2_i,gh_struc,exgh_data)
          gh_data(cnt) = gh_struc

! Add in dB data
          gh_data(cnt)%bR = 2.0/3.0*struc_data(idx_Trk1(Tr1_i))%bR + &
                            1.0/3.0*struc_data(idx_Trk2(Tr2_i))%bR
          exgh_data%bR = 1.0/3.0*struc_data(idx_Trk1(Tr1_i))%bR + &
                         2.0/3.0*struc_data(idx_Trk2(Tr2_i))%bR

          gh_data(cnt)%bT = 2.0/3.0*struc_data(idx_Trk1(Tr1_i))%bT + &
                            1.0/3.0*struc_data(idx_Trk2(Tr2_i))%bT
          exgh_data%bT = 1.0/3.0*struc_data(idx_Trk1(Tr1_i))%bT + &
                         2.0/3.0*struc_data(idx_Trk2(Tr2_i))%bT

          gh_data(cnt)%bP = 2.0/3.0*struc_data(idx_Trk1(Tr1_i))%bP + &
                            1.0/3.0*struc_data(idx_Trk2(Tr2_i))%bP
          exgh_data%bP = 1.0/3.0*struc_data(idx_Trk1(Tr1_i))%bP + &
                         2.0/3.0*struc_data(idx_Trk2(Tr2_i))%bP

! Convert spherical dB vector to x,y,z
          call bspcar_08(gh_data(cnt)%T, &
                         gh_data(cnt)%P, &
                         gh_data(cnt)%bR, &
                         gh_data(cnt)%bT, &
                         gh_data(cnt)%bP, &
                         xv, yv, zv)
          gh_data(cnt)%bX = xv
          gh_data(cnt)%bY = yv
          gh_data(cnt)%bZ = zv
          call bspcar_08(exgh_data%T, &
                         exgh_data%P, &
                         exgh_data%bR, &
                         exgh_data%bT, &
                         exgh_data%bP, &
                         xv, yv, zv)
          exgh_data%bX = xv
          exgh_data%bY = yv
          exgh_data%bZ = zv

! write exgh_data to temp file
          if (exgh_data%T*radToDeg > 1.0 .and. &
              exgh_data%T*radToDeg < 179.0) then
            write(1,*) exgh_data
          else
            iiex_gh = iiex_gh - 1   ! don't count this exghost
          endif

        endif         ! if we need 2 extra ghosts

! if we do not need extras - i.e. ghost only 1 point in Longitude
        if (Nyq_Ok .eq. 1) then 
          do_exgh = 0
          call pop_ghost_struc(Trk1_p,Trk2_p, &
                   Tr1_coLat_deg, Tr1_Lon_deg, &
                   Tr2_coLat_deg, Tr2_Lon_deg, &
                   coLat_wdth, xs,ys,xe,ye,south,do_exgh, &
                   Tr1_i,Tr2_i,gh_struc,exgh_data)

          gh_data(cnt) = gh_struc

          gh_data(cnt)%bR = (struc_data(idx_Trk1(Tr1_i))%bR + &
! Add in dB data
                             struc_data(idx_Trk2(Tr2_i))%bR)/2.0d0

          gh_data(cnt)%bT = (struc_data(idx_Trk1(Tr1_i))%bT + &
                             struc_data(idx_Trk2(Tr2_i))%bT)/2.0d0

          gh_data(cnt)%bP = (struc_data(idx_Trk1(Tr1_i))%bP + &
                             struc_data(idx_Trk2(Tr2_i))%bP)/2.0d0

! Convert spherical dB vector to x,y,z
          call bspcar_08(gh_data(cnt)%T, &
                         gh_data(cnt)%P, &
                         gh_data(cnt)%bR, &
                         gh_data(cnt)%bT, &
                         gh_data(cnt)%bP, &
                         xv, yv, zv)
          gh_data(cnt)%bX = xv
          gh_data(cnt)%bY = yv
          gh_data(cnt)%bZ = zv

! ###### For TESTING ghost routine ######################
! if (Tr_num .eq. 0) then
!   write(13,*) jj, Tr1_i,Tr2_i, Trk1_p%T*radtodeg, Trk1_p%P*radtodeg, & 
!      gh_data(cnt)%T*radtodeg, gh_data(cnt)%P*radtodeg, &
!       Trk2_p%T*radtodeg, Trk2_p%P*radtodeg
!
!   write(14,*) jj, idx_Trk1(Tr1_i), idx_Trk2(Tr2_i), &
!            gh_data(cnt)%T*radtodeg, gh_data(cnt)%P*radtodeg, &
!                 struc_data(idx_Trk1(Tr1_i))%bT, &
!                 struc_data(idx_Trk2(Tr2_i))%bT, &
!                 gh_data(cnt)%bT
! endif
! #######################################################
        endif           ! If ghost 1 data point

        if (gh_data(cnt)%T*radToDeg > 1.0 .and. &
            gh_data(cnt)%T*radToDeg < 179.0) then
          cnt = cnt + 1
        endif

      enddo point_ex_loop

      deallocate(Tr2_Lon_deg)
      deallocate(Tr2_coLat_deg)
      deallocate(Tr1_Lon_deg)
      deallocate(Tr1_coLat_deg)
      deallocate(Tr2_aacgm_Lon_deg)
      deallocate(Tr2_aacgm_coLat_deg)
      deallocate(Tr2_GEO_Lon_deg)
      deallocate(Tr2_GEO_coLat_deg)
      deallocate(Tr1_aacgm_Lon_deg)
      deallocate(Tr1_aacgm_coLat_deg)
      deallocate(Tr1_GEO_Lon_deg)
      deallocate(Tr1_GEO_coLat_deg)
      deallocate(gcDist)
      deallocate(idx_Trk2)
      deallocate(idx_Trk1)

    enddo track_pair_loop

    if (tmp_open .eq. 1) then
      close(unit=1)
    endif

! allocate mem for output combined data structure
    tot_rec = np + cnt-1 + iiex_gh

    allocate(comb_data(tot_rec))
    comb_data(1:np) = struc_data(1:np)  ! load input data

    do jj=1,cnt-1                       ! load ghost data
     comb_data(np+jj) = gh_data(jj)
    enddo

    if (iiex_gh .gt. 0) then            ! load extra ghost data
      open(unit=1,file=gh_fname,status='old', &
           action='read')
      do jj=1,iiex_gh       ! loop for number of records in temp file
        read(1,*) exgh_data
        comb_data(np+cnt-1+jj) = exgh_data
      enddo
      close(unit=1)
    endif     ! if we have extra ghosts

!    print*,'DONE_aacgm_ghosts'

  end subroutine ampFit_ghosts_aacgm

! ---------------------------------------------------------------------------------------------------------

  subroutine pop_ghost_struc(Tr1_struc_p,Tr2_struc_p, &
                   Trk1_coLat_deg_arr,Trk1_Lon_deg_arr, &
                   Trk2_coLat_deg_arr,Trk2_Lon_deg_arr, &
                   coLat_wdth, xs,ys,xe,ye,south, c_exgh, &
                   Tr1_i,Tr2_i,gh_struc,exgh_struc)

! Populate the ghost data structure
!
! Tr1_struc_p : current data point on Track_1 (structure)
! Tr2_struc_p : closest point (by gcdist) on Track_2 (structure)
! Trk1_coLat_deg_arr : [np] of coLat_deg values. If Tr1_struc_p(coLat_deg) is within 20 deg
!                      of the equator then this array contains GEO coLat else AACGM coLats
! Trk1_Lon_deg_arr : [np] of Lon_deg values. If Tr1_struc_p(coLat_deg) is within 20 deg
!                      of the equator then this array contains GEO Lons else AACGM Lons
! coLat_wdth : p/m deg to search either side of gh_pivot location
! Xs,Xe,Ys,Ye : XY coords of Tr1 and Tr2 points (proj in Nth hemis)
! south : 0=Nth, 1=Sth hemis
! c_exgh : 0=no extra ghosts, 1=calc extra ghost point
! Tr1_i, Tr2_i : indexes of Tr1 and Tr2 that are closest to gh_pivot location
! gh_struc : output ghost data point (structure)
! exgh_struc : extra ghost data point (structure)
!
! Outputs:
!    Tr1_i, Tr2_i, gh_struc, exgh_struc (if c_exgh=1)
! 
! C.L. Waters
! Centre for Space Phyics
! University of Newcastle
! Australia
! Nov 2011
!
! Ver : 201202
!
    implicit none

    type(ampData), intent(in) :: Tr1_struc_p, Tr2_struc_p
    type(ampData), intent(out) :: gh_struc, exgh_struc

    real(kind=DBL), intent(in) :: Trk1_coLat_deg_arr(:), &
                                  Trk1_Lon_deg_arr(:), &
                                  Trk2_coLat_deg_arr(:), &
                                  Trk2_Lon_deg_arr(:)

    real(kind=DBL), intent(in) :: coLat_wdth,xs,ys,xe,ye

    integer, intent(in) :: south, c_exgh
    integer, intent(inout) :: Tr1_i, Tr2_i

    logical, dimension(:), allocatable :: wlogic_arr
    integer, dimension(:), allocatable :: idx_logic, c_idx, idx_srt

    real(kind=DBL), dimension(:), allocatable :: dff_coLat

    integer :: ii, jj, iij, flg, np, ii_c, s

    real(kind=DBL), target :: aacgm_lat_d, aacgm_lon_d, err
    real(kind=DBL) :: st_clat, st_lon
    real(kind=DBL) :: m1,m2,xg1,xg2,yg1,yg2,hgt
    real(kind=DBL) :: xv,yv,zv,mid_x,mid_y, dff_lon
    real(kind=DBL) :: gh_pivot_coLat_deg, gh_pivot_Lon_deg

! temp arrays for sort
    real(kind=DBL), dimension(:), allocatable :: T
    integer, dimension(:), allocatable :: Ti

    hgt = rSat/1000.0d0    ! altitude in km for aacgm
!  calc ghost structure components
    gh_struc%utc = (Tr1_struc_p%utc + Tr2_struc_p%utc)/2.0d0
    gh_struc%isat = -1
    gh_struc%iPln = Tr1_struc_p%iPln
    gh_struc%splice=Tr1_struc_p%splice
    exgh_struc = gh_struc

    if (c_exgh .eq. 0) then
      m1 = 0.5
      m2 = 0.5
      gh_struc%typ = 2
    else
      m1 = 2.0/3.0          ! mult on gh_struc
      m2 = 1.0/3.0
      gh_struc%typ = 3
      exgh_struc%typ = 3
    endif
    Xg1 = m1*Xs + m2*Xe
    Yg1 = m1*Ys + m2*Ye
! For extra ghosts
    Xg2=m2*Xs + m1*Xe
    Yg2=m2*Ys + m1*Ye

    gh_struc%qual = m1*Tr1_struc_p%qual + m2*Tr2_struc_p%qual
    gh_struc%R = m1*Tr1_struc_p%R + m2*Tr2_struc_p%R
    if (c_exgh .eq. 1) then
      exgh_struc%qual = m2*Tr1_struc_p%qual + m1*Tr2_struc_p%qual
      exgh_struc%R = m2*Tr1_struc_p%R + m1*Tr2_struc_p%R
    endif

    if (south .eq. 1) then   ! south hemis
      gh_struc%T = (180.0 - sqrt(xg1**2 + yg1**2))*degToRad
      if (c_exgh .eq. 1) exgh_struc%T = (180.0 - sqrt(xg2**2 + yg2**2))*degToRad
    else                     ! North hemisphere
      gh_struc%T = (sqrt(xg1**2 + yg1**2))*degToRad
      if (c_exgh .eq. 1) exgh_struc%T = (sqrt(xg2**2 + yg2**2))*degToRad
    endif

! Convert spherical coords to x,y,z
    gh_struc%P = atan2(yg1,xg1)
    call sphcar_08(gh_struc%R, &
                   gh_struc%T, &
                   gh_struc%P, &
                   xv, yv, zv, 1)
    gh_struc%X = xv
    gh_struc%Y = yv
    gh_struc%Z = zv

    if (c_exgh .eq. 1) then
      exgh_struc%P = atan2(yg2,xg2)
      call sphcar_08(exgh_struc%R, &
                     exgh_struc%T, &
                     exgh_struc%P, &
                   xv, yv, zv, 1)
      exgh_struc%X = xv
      exgh_struc%Y = yv
      exgh_struc%Z = zv
    endif

! Initialise pivot point at ghost in GEO - always at midpoint
    mid_x = (xs+xe)/2.0d0
    mid_y = (ys+ye)/2.0d0
    if (south .eq. 1) then 
      gh_pivot_coLat_deg = 180.0d0-sqrt(mid_x**2 + mid_y**2)
    else
      gh_pivot_coLat_deg = sqrt(mid_x**2 + mid_y**2)
    endif
    gh_pivot_Lon_deg = atan2(mid_y,mid_x)*180.0d0/pi

! Select dB on the basis of AACGM coLat if within valid AACGM Latitude
    if ( (gh_struc%T*radToDeg .lt. 70.0) .or. &
         (gh_struc%T*radToDeg .gt. 110.0) ) then
      flg=0
! (i) AACGM coords of gh point
      st_clat = gh_struc%T*radToDeg
      st_lon = gh_struc%P*radToDeg
      s = f_AACGMConvert(90.0-st_clat, st_lon, hgt, &
                         C_LOC(aacgm_lat_d), &
                         C_LOC(aacgm_lon_d), &
                         C_LOC(err), flg)
!      gh_pivot_lon_deg = aacgm_lon_d     ! work lon in GEO
      gh_pivot_coLat_deg = 90.0-aacgm_lat_d
    endif

! (ii) find idx for where Tr1 is within (+-)coLat_wdth
    np = size(Trk1_coLat_deg_arr)
    allocate(idx_logic(np))
    allocate(wlogic_arr(np))
    wlogic_arr = .false.        ! initialise logic array
    wlogic_arr = &
       (Trk1_coLat_deg_arr <= (gh_pivot_coLat_deg+coLat_wdth) .and. &
        Trk1_coLat_deg_arr >= (gh_pivot_coLat_deg-coLat_wdth))
    call my_where(wlogic_arr, np, idx_logic, ii_c)
    deallocate(wlogic_arr)
    allocate(c_idx(ii_c))         ! get mem for indexes for Trk_1
    c_idx = idx_logic(1:ii_c)     ! indexes for these vals
    deallocate(idx_logic)

! (iii) Search for closest location
    if (ii_c .eq. 1) Tr1_i = c_idx(1)
    if (ii_c .gt. 1) then
      allocate(dff_coLat(ii_c))
      do ii = 1,ii_c
        dff_coLat(ii) = abs(gh_pivot_coLat_deg-Trk1_coLat_deg_arr(c_idx(ii)) )
! ##### For testing ghost routine ###############################
! write(15,*) 'Trk_1 ',c_idx(ii),Trk1_coLat_deg_arr(c_idx(ii)), &
!              Trk1_lon_deg_arr(c_idx(ii))
! ###############################################################
      enddo
      allocate(idx_srt(ii_c))
      allocate(T((ii_c+1)/2))
      allocate(Ti((ii_c+1)/2))
      idx_srt = c_idx         ! sort these indexes
      call merge_sort(dff_coLat, ii_c, T, idx_srt, Ti)
      deallocate(Ti)
      deallocate(T)
      iij = 0
!  While loop to find min distance
      DO
        iij = iij+1
        dff_lon = abs(gh_pivot_Lon_deg-Trk1_lon_deg_arr(idx_srt(iij)))
        if (dff_lon .gt. 180.0) dff_lon=abs(360.0-dff_lon)
! ##### For testing ghost routine ############################
! write(15,*) iij, gh_pivot_coLat_deg, gh_pivot_lon_deg, &
!             Trk1_coLat_deg_arr(idx_srt(iij)), &
!             Trk1_lon_deg_arr(idx_srt(iij))
! ###########################################################
        if ( (dff_lon .lt. 30.0) .or. (iij .ge. ii_c) ) EXIT
      END DO
      Tr1_i = idx_srt(iij)
      deallocate(idx_srt)
      deallocate(dff_colat)
    endif
    deallocate(c_idx)
! ##### For testing ghost routine #####################
! write(15,*) 'Trk_1 idx = ',Tr1_struc_p%ipln, Tr1_i
! ####################################################

! (iv) Examine Track 2
! find idx for where Tr2 is within (+-)coLat_wdth
    np = size(Trk2_coLat_deg_arr)
    allocate(idx_logic(np))
    allocate(wlogic_arr(np))
    wlogic_arr = .false.        ! initialise logic array
    wlogic_arr = &
       (Trk2_coLat_deg_arr <= (gh_pivot_coLat_deg+coLat_wdth) .and. &
        Trk2_coLat_deg_arr >= (gh_pivot_coLat_deg-coLat_wdth))
    call my_where(wlogic_arr, np, idx_logic, ii_c)
    deallocate(wlogic_arr)
    allocate(c_idx(ii_c))       ! get mem for indexes for Trk_1
    c_idx = idx_logic(1:ii_c)   ! indexes for these vals
    deallocate(idx_logic)

! Search for closest location
    if (ii_c .eq. 1) Tr2_i = c_idx(1)
    if (ii_c .gt. 1) then
      allocate(dff_coLat(ii_c))
      do ii = 1,ii_c
        dff_coLat(ii) = abs(gh_pivot_coLat_deg-Trk2_coLat_deg_arr(c_idx(ii)) )  ! CHECK
! ##### For testing ghost routine #####################
! write(14,*) 'Trk_2 ',c_idx(ii),Trk2_coLat_deg_arr(c_idx(ii)), &
!              Trk2_lon_deg_arr(c_idx(ii))
! #####################################################
      enddo
      allocate(idx_srt(ii_c))
      allocate(T((ii_c+1)/2))
      allocate(Ti((ii_c+1)/2))
      idx_srt = c_idx         ! sort these indexes
      call merge_sort(dff_coLat, ii_c, T, idx_srt, Ti)
      deallocate(Ti)
      deallocate(T)
      iij = 0
      DO
        iij = iij+1
        dff_lon = abs(gh_pivot_Lon_deg-Trk2_lon_deg_arr(idx_srt(iij)))
        if (dff_lon .gt. 180.0) dff_lon=abs(360.0-dff_lon)
        if ( (dff_lon .lt. 30.0) .or. (iij .ge. ii_c) ) EXIT
      END DO
      Tr2_i = idx_srt(iij)
      deallocate(idx_srt)
      deallocate(dff_colat)
    endif
    deallocate(c_idx)

!write(14,*) 'Trk_2 idx = ',Tr2_i

  end subroutine pop_ghost_struc

end module ampFit_sort
!
! =====================================================================
