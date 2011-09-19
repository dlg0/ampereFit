module ampFit_sort
! This nodule contains:
! 1. do_merge -> subroutine for merge_sort
! 2. merge_sort -> sort routine adapted from Wiki page on merge sort
!                  This routine has the sorted data plus the indexes
! 3. cross_p -> vector cross product
! 4. norm_vec -> normailise 3D vectors
! 5. gcirc_dist -> great circle distance (over a sphere) routine
! 6. my_where -> similar to IDL where function for 1D arrays
! 7. sort_struc -> sorts satelliet track data points by great circle distance
! 8. tag_lon_strays -> finds stray satellites (off longitude track)
! 9. ampFit_ghost -> adds ghost data near poles
!
! CLW - Aug 2011
!
  use constants

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
!   CLW:  22 Aug 2011
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
! CLW - Aug 2011

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
! CLW:  Aug, 2011

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

  subroutine sort_struc(struc_data, trk_order)
! Given the shiftedData structure, this routine:
! (i)   Cycles through each Track Number/hemisphere
! (ii)  Identify and label any strays (off longitude)
! (iii) Sort by great circle distance from an equator point
!  (iv) ensures the tracks are sequential in longitude
!
! Each hemisphere is sorted separately.
! Output is the modified shiftedData structure that has
! - Strays tagged in the typ property, where
!    data%typ = 0 -> normal data
!    data%typ = 1 -> stray data off longitude
!    data%typ = 2 -> ghost data points near the poles
!    data%typ = 3 -> additional ghosts where longitude separation is > 120 deg
!
! The resulting data structure (sorted): 
! Tr0(Nth,Sth), Tr1(Nth,Sth) etc. sequence
!  and track order by longitude, ready for ghost data calc
!
! CLW - Aug 2011
!
    implicit none

    type(ampData), intent(inout) :: struc_data(:)
    integer, intent(inout) :: trk_order(6)

    logical, dimension(1:size(struc_data)) :: tmplogic_arr
    integer, dimension(1:size(struc_data)) :: tmp_idx, sortII
    integer :: cnt, i, np, Tr_num, iiTrack_n, iiTrack_s
    integer :: iiTrack_nlon
    integer, dimension(:), allocatable :: idx_n, idx_nlon, idx_nsrt
    integer, dimension(:), allocatable :: idx_s, idx_ssrt
    real(kind=DBL), dimension(:), allocatable :: gcDist_n, gcDist_s
    real(kind=DBL) :: st_lat, st_lon
    real(kind=DBL), dimension(6) :: strt_lona
    real(kind=DBL), dimension(3) :: uvec, vvec, nvec, ncuvec
    integer, dimension(1) :: iiSt_n
! temp arrays for sort routine
    real(kind=DBL), dimension(:), allocatable :: T
    integer, dimension(:), allocatable :: Ti

    np = size(struc_data)
    cnt = 1

print*,'Sorting each data track...'

    sort_track_loop: &
    do Tr_num = 0,5

! Get indexes where data has correct TrackNum and Bth hemisphere
      tmplogic_arr = .false.          ! initialise logic array
      tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg <= 90.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_n)
      allocate(idx_n(iiTrack_n))      ! get mem for these data
      idx_n = tmp_idx(1:iiTrack_n)    ! indexes for these vals

! Get indexes for (i) Track (ii) Nth hemis (iii) < 180 in lon 
      tmplogic_arr = .false.
      tmplogic_arr=struc_data(idx_n)%P*radToDeg <= 180.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_nlon)
      allocate(idx_nlon(iiTrack_nlon))    ! get mem for these data
      idx_nlon = tmp_idx(1:iiTrack_nlon)  ! indexes for these vals
! These are the indexes of struc for Nth, correct Tr_num and lon <= 180

      iiSt_n = maxloc(struc_data(idx_n(idx_nlon))%T*radToDeg)  ! index of equator point
      i = idx_n(idx_nlon(iiSt_n(1)))
      st_lat = 90.0 - struc_data(i)%T*radToDeg  ! Lat of this point
      st_lon  = struc_data(i)%P*radToDeg          ! Lon of this point
! st_Lat,st_Lon is the coord of the great circle dist reference point

      strt_lona(Tr_num+1) = st_lon         ! store start longitude
      trk_order(Tr_num+1) = Tr_num
      allocate(gcDist_n(iiTrack_n))        ! iiTrack_n is the size of idx_n
      call gcirc_dist(st_lat, 90.0-struc_data(idx_n)%T*radToDeg, &
             st_lon, struc_data(idx_n)%P*radToDeg, iiTrack_n, gcDist_n)

! sort the distances, returns idx of sorted array
      allocate(T((iiTrack_n+1)/2))
      allocate(Ti((iiTrack_n+1)/2))
      allocate(idx_nsrt(iiTrack_n))
      idx_nsrt = idx_n              ! sort these indexes
      call merge_sort(gcDist_n, iiTrack_n, T, idx_nsrt, Ti)
      deallocate(Ti)
      deallocate(T)
      deallocate(gcDist_n)

      sortII(cnt:cnt+iiTrack_n-1) = idx_nsrt  ! save sorted indexes
      cnt = cnt + iiTrack_n
 
      deallocate(idx_nlon)
      deallocate(idx_n)
      deallocate(idx_nsrt)

! ------  do Sth hemisphere ------

      tmplogic_arr = .false.          ! initialise logic array
      tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%T*radToDeg > 90.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_s)
      allocate(idx_s(iiTrack_s))      ! get mem for these data
      idx_s = tmp_idx(1:iiTrack_s)    ! indexes for these vals

      allocate(gcDist_s(iiTrack_s))           ! iiTrack_s is the size of idx_s
      call gcirc_dist(st_lat, 90.0-struc_data(idx_s)%T*radToDeg, &
             st_lon, struc_data(idx_s)%P*radToDeg, iiTrack_s, gcDist_s)

! sort the distances, returns idx of sorted array
      allocate(T((iiTrack_s+1)/2))
      allocate(Ti((iiTrack_s+1)/2))
      allocate(idx_ssrt(iiTrack_s))
      idx_ssrt = idx_s    ! sort these indexes
      call merge_sort(gcDist_s, iiTrack_s, T, idx_ssrt, Ti)
      deallocate(Ti)
      deallocate(T)
      deallocate(gcDist_s)

      sortII(cnt:cnt+iiTrack_s-1) = idx_ssrt  ! save sorted indexes
      cnt = cnt + iiTrack_s

      deallocate(idx_s)
      deallocate(idx_ssrt)

    enddo sort_track_loop
    struc_data = struc_data( sortII )   ! sort the structure
!do i=1,50
! print*,struc_data(i)%T*radToDeg,struc_data(i)%P*radToDeg
!enddo

! sort the track order, returns trk_order
    allocate(T(6))
    allocate(Ti(6))
    call merge_sort(strt_lona, 6, T, trk_order, Ti)
    deallocate(Ti)
    deallocate(T)

!print*,'Track_order = ',trk_order
print*,'DONE_sort'
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
! CLW - Aug 2011
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

    print*,'Tagging longitude strays...'
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
    print*,'DONE'

  end subroutine tag_lon_strays
!
! ----------------------------------------------------------------------
!
  subroutine ampFit_ghost(struc_data, comb_data, clat_lim, trk_order)
! Given:
!   shiftedData structure,
!   clat_lim, the pole limit for ghosting
!   trk_order (integer arry of 6 of the track order)
! This routine:
!  (i) Cycles through adjacent data track pairs
! (ii) Identify coLat limit and add data between longitudes
!(iii) Check if 2 ghosts are required. This is usually the case
!      in track the intesection area near the poles
!
! Assumes the data for each track have been sorted by sort_struc
! - Tr0(Nth,Sth), Tr1(Nth,Sth) etc sequence sorted
!
! CLW - Aug 2011
!
    implicit none

    type(ampData), intent(in) :: struc_data(:)
    type(ampData), allocatable, intent(out) :: comb_data(:)
    real(kind=DBL), intent(in) :: clat_lim
    integer, dimension(6), intent(in) :: trk_order

    type(ampData), allocatable :: gh_data(:)
    type(ampData) :: exgh_data
    logical, dimension(1:size(struc_data)) :: tmplogic_arr
    integer, dimension(1:size(struc_data)) :: tmp_idx
    integer, dimension(:), allocatable :: idx_Trk1, idx_Trk2
    integer :: np, Tr_num, iiTrk_1, iiTrk_2, nxt_trk, iiex_gh
    integer :: jj, Nyq_Ok, cnt, iighost, idx1, idx2
    integer :: tot_rec, tmp_open
    integer, dimension(1) :: iimn
    real(kind=DBL), dimension(:), allocatable :: gcDist
    real(kind=DBL) :: st_lat, st_lon, e_lat, e_lon, xe, ye
    real(kind=DBL) :: xs, ys, xg, yg, xg1, yg1, xg2, yg2
    real(kind=DBL) :: xv, yv, zv
    character(len=20) :: gh_fname

    integer :: i       ! for testing

    print*,'Adding ghost data for CoLat: ',clat_lim
    np = size(struc_data)
    iiex_gh = 0
    cnt = 1
    gh_fname='ex_ghosts.dat'
    tmp_open = 0

! Calc min number of ghost points (excludes extra ghosts)
    tmplogic_arr = .false.        ! initialise logic array
    if (clat_lim .lt. 90.0) then  ! Nth hemis data
      tmplogic_arr = struc_data%T*radToDeg <= clat_lim .and. &
                     struc_data%typ==0
    else                          ! Sth hemis data
      tmplogic_arr = struc_data%T*radToDeg >= clat_lim .and. &
                     struc_data%typ==0
    endif
    call my_where(tmplogic_arr, np, tmp_idx, iighost)
    allocate(gh_data(iighost)) ! get mem for ghost_data struc

    track_pair_loop: &
    do Tr_num = 0,5

! mask for correct Tr_num and hemis
      tmplogic_arr = .false.        ! initialise logic array
      if (clat_lim .lt. 90.0) then  ! Nth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(Tr_num+1) .and. &
                       struc_data%T*radToDeg <= clat_lim .and. &
                       struc_data%typ==0
      else                          ! Sth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(Tr_num+1) .and. &
                       struc_data%T*radToDeg >= clat_lim .and. &
                       struc_data%typ==0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_1)
      allocate(idx_Trk1(iiTrk_1))   ! get mem for indexes for Trk_1
      idx_Trk1 = tmp_idx(1:iiTrk_1) ! indexes for these vals

! if iiTrk1 = 0 then no data for ghosts - error

      if (Tr_num .lt. 5) then       ! Track number s 0..5
        nxt_trk=Tr_num+2            ! index numbers 1..6
      else
         nxt_trk=1                  ! fold back to 1st track in sorted sequence
      endif

      tmplogic_arr = .false.        ! initialise logic array
      if (clat_lim .lt. 90.0) then  ! Nth hemis data
! Get data for adjacent satellite track -> Trk_2
        tmplogic_arr = struc_data%ipln==trk_order(nxt_trk) .and. &
                       struc_data%T*radToDeg <= clat_lim .and. &
                       struc_data%typ==0
      else                          ! Sth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(nxt_trk) .and. &
                       struc_data%T*radToDeg >= clat_lim .and. &
                       struc_data%typ==0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_2)
      allocate(idx_Trk2(iiTrk_2))   ! get mem for indexes for Trk_2 data
      idx_Trk2 = tmp_idx(1:iiTrk_2) ! indexes for these vals
      allocate(gcDist(iiTrk_2))     ! mem for great circ dist values
! if iiTrk2 = 0 then no data for ghost - error

      point_ex_loop: &
      do jj=1,iiTrk_1                ! loop thru all points on Trk_1
        st_lat = struc_data(idx_Trk1(jj))%T*radToDeg
        st_lon = struc_data(idx_Trk1(jj))%P*radToDeg
        call gcirc_dist(90.0-st_lat, 90.0-struc_data(idx_Trk2)%T*radToDeg, &
             st_lon, struc_data(idx_Trk2)%P*radToDeg, iiTrk_2, gcDist)

! iimn is the index of idx_Trk2 which has min dist
!  between Trk2 and the jj point of Trk1
        iimn = minloc(gcDist)        ! index of idx_Trk2 of min_dist
        
        if (clat_lim .gt. 90.0) then ! Sth hemis data
          st_lat=180.0-st_lat
        endif
        xs = st_lat*cos(st_lon*pi/180.0)
        ys = st_lat*sin(st_lon*pi/180.0)
        e_lat=struc_data(idx_Trk2(iimn(1)))%T*radToDeg
        if (clat_lim .gt. 90.0) then
          e_lat=180.0-e_lat
        endif
        e_lon=struc_data(idx_Trk2(iimn(1)))%P*radToDeg
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

        xe = e_lat*cos(e_lon*pi/180.0)
        ye = e_lat*sin(e_lon*pi/180.0)
        if (Nyq_Ok .eq. 1) then   ! need 1 ghost point inserted
          xg = (xs + xe)/2.0
          yg = (ys + ye)/2.0
        else                      ! need 2 ghost points inserted
          xg1 = 2.0/3.0*xs + xe/3.0
          yg1 = 2.0/3.0*ys + ye/3.0
          xg2 = xs/3.0 + 2.0/3.0*xe
          yg2 = ys/3.0 + 2.0/3.0*ye
        endif

! Put ghost data into ghost data structure
        idx1 = idx_Trk1(jj)
        idx2 = idx_Trk2(iimn(1))
        gh_data(cnt)%utc = (struc_data(idx1)%utc+struc_data(idx2)%utc)/2.0
        gh_data(cnt)%isat = -1
        gh_data(cnt)%ipln = struc_data(idx1)%ipln
        gh_data(cnt)%splice = struc_data(idx1)%splice

        if (Nyq_Ok .eq. 0) then  ! if we need 2 extra ghosts
          exgh_data%utc = gh_data(cnt)%utc
          exgh_data%isat = -1
          exgh_Data%ipln = struc_data(idx1)%ipln

          gh_data(cnt)%qual = 2.0/3.0*struc_data(idx1)%qual + &
                                1.0/3.0*struc_data(idx2)%qual
          exgh_data%qual = 1.0/3.0*struc_data(idx1)%qual + &
                                2.0/3.0*struc_data(idx2)%qual

          exgh_data%splice = struc_data(idx1)%splice
          gh_data(cnt)%typ = 3
          exgh_data%typ = 3
          
          gh_data(cnt)%R = 2.0/3.0*struc_data(idx1)%R + &
                                1.0/3.0*struc_data(idx2)%R
          exgh_data%R = 1.0/3.0*struc_data(idx1)%R + &
                                2.0/3.0*struc_data(idx2)%R

          if (clat_lim .gt. 90.0) then   ! south hemis
            gh_data(cnt)%T = (180.0 - sqrt(xg1**2 + yg1**2))*degToRad
            exgh_data%T = (180.0 - sqrt(xg2**2 + yg2**2))*degToRad
          else                           ! North hemisphere
            gh_data(cnt)%T = (sqrt(xg1**2 + yg1**2))*degToRad
            exgh_data%T = (sqrt(xg2**2 + yg2**2))*degToRad
          endif

          gh_data(cnt)%T = gh_data(cnt)%T*radToDeg*pi/180.0
          exgh_data%T = exgh_data%T*radToDeg*pi/180.0

          gh_data(cnt)%P = atan2(yg1,xg1)
          exgh_data%P = atan2(yg2,xg2)
! Convert spherical coords to x,y,z
          call sphcar_08(gh_data(cnt)%R, &
                         gh_data(cnt)%T, &
                         gh_data(cnt)%P, &
                         xv, yv, zv, 1)
          gh_data(cnt)%X = xv
          gh_data(cnt)%Y = yv
          gh_data(cnt)%Z = zv
          call sphcar_08(exgh_data%R, &
                         exgh_data%T, &
                         exgh_data%P, &
                         xv, yv, zv, 1)
          exgh_data%X = xv
          exgh_data%Y = yv
          exgh_data%Z = zv

          gh_data(cnt)%bR = 2.0/3.0*struc_data(idx1)%bR + &
                                1.0/3.0*struc_data(idx2)%bR
          exgh_data%bR = 1.0/3.0*struc_data(idx1)%bR + &
                                2.0/3.0*struc_data(idx2)%bR

          gh_data(cnt)%bT = 2.0/3.0*struc_data(idx1)%bT + &
                                1.0/3.0*struc_data(idx2)%bT
          exgh_data%bT = 1.0/3.0*struc_data(idx1)%bT + &
                                2.0/3.0*struc_data(idx2)%bT

          gh_data(cnt)%bP = 2.0/3.0*struc_data(idx1)%bP + &
                                1.0/3.0*struc_data(idx2)%bP
          exgh_data%bP = 1.0/3.0*struc_data(idx1)%bP + &
                                2.0/3.0*struc_data(idx2)%bP
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
          gh_data(cnt)%qual = (struc_data(idx1)%qual + &
                               struc_data(idx2)%qual)/2.0
          gh_data(cnt)%typ = 2
          
          gh_data(cnt)%R = (struc_data(idx1)%R + &
                                   struc_data(idx2)%R)/2.0

          if (clat_lim .gt. 90.0) then   ! south hemis
            gh_data(cnt)%T = (180.0 - sqrt(xg**2 + yg**2))*degToRad
          else                           ! North hemisphere
            gh_data(cnt)%T = (sqrt(xg**2 + yg**2))*degToRad
          endif

          gh_data(cnt)%T = gh_data(cnt)%T*radToDeg*pi/180.0

          gh_data(cnt)%P = atan2(yg,xg)

          call sphcar_08(gh_data(cnt)%R, &
                         gh_data(cnt)%T, &
                         gh_data(cnt)%P, &
                         xv, yv, zv, 1)
          gh_data(cnt)%X = xv
          gh_data(cnt)%Y = yv
          gh_data(cnt)%Z = zv

          gh_data(cnt)%bR = struc_data(idx1)%bR + &
                                struc_data(idx2)%bR

          gh_data(cnt)%bT = (struc_data(idx1)%bT + &
                                     struc_data(idx2)%bT)/2.0

          gh_data(cnt)%bP = (struc_data(idx1)%bP + &
                                   struc_data(idx2)%bP)/2.0

          call bspcar_08(gh_data(cnt)%T, &
                         gh_data(cnt)%P, &
                         gh_data(cnt)%bR, &
                         gh_data(cnt)%bT, &
                         gh_data(cnt)%bP, &
                         xv, yv, zv)
          gh_data(cnt)%bX = xv
          gh_data(cnt)%bY = yv
          gh_data(cnt)%bZ = zv

        endif    ! one ghost point

        if (gh_data(cnt)%T*radToDeg > 1.0 .and. &
            gh_data(cnt)%T*radToDeg < 179.0) then
          cnt = cnt + 1
        endif

      enddo point_ex_loop

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

    print*,'DONE_ghosts'

  end subroutine ampFit_ghost

end module ampFit_sort
