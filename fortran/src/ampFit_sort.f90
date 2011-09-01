module ampFit_sort
! Merge Sort algorithm for AmpFit
! CLW - Aug 2011
! Adapted from the Wiki page on merge Sort - I added the index sort capability
!

  use constants
!  integer, parameter :: DBL = selected_real_kind ( p=13, r=200)

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
!print*,'In do_merge'
!print*,'A ', (A(ii),ii=1,na)
!print*,'Ai ', (Ai(ii),ii=1,na)
!print*,'B ', (B(ii),ii=1,nb)
!print*,'Bi ', (Bi(ii),ii=1,nb)
!print*,'C ', (C(ii),ii=1,nc)
!print*,'Ci ', (Ci(ii),ii=1,nc)
    do while (i <= na .and. j <= nb)
      if (A(i) <= B(j)) then
        C(k)=A(i)
        Ci(k)=Ai(i)
!print*,i,k,A(i),B(j)
        i=i+1
      else
        C(k)=B(j)
        Ci(k)=Bi(j)
!print*,j,k,B(j)
        j=j+1
      end if
      k=k+1
    enddo
    do while (i <= na)
      C(k)=A(i)
      Ci(k)=Ai(i)
!print*,i,k,A(i)
      k=k+1
      i=i+1
    enddo
!print*,'end do_merge***'
    return
  end subroutine do_merge

!------------------------------------------------------------------
  
  recursive subroutine merge_sort(A,N,T,Ai,Ti)
! merge sort with recursive call - CLW 22 Aug 2011
! The arrays A and Ai are replaced by t6he sorted results
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
!print*,'after do_merge'
!print*,'n,na,nb=',n,na,nb
!print*,'A ', (A(ii),ii=1,n)
!print*,'Ai ', (Ai(ii),ii=1,N)
!print*,'T ', (T(ii),ii=1,(n+1)/2)
!print*,'Ti ', (Ti(ii),ii=1,(n+1)/2)
    endif
    return

  end subroutine merge_sort

!--------------------------------------------------------------

  function cross_p(vec1, vec2)
    implicit none
    real(kind=DBL), dimension(3) :: cross_p
    real(kind=DBL), dimension(3), intent(in) :: vec1, vec2

    cross_p(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_p(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_p(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

  end function cross_p

!--------------------------------------------------------------

  function norm_vec(vecin)
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
!    real(kind=DBL) :: gcirc_dist
! calc great circle distance given coords as lat,lom
! CLW - Aug 2011

    implicit none

    integer, intent(in) :: np
    real(kind=DBL), intent(in) :: lat_s, lon_s
    real(kind=DBL), dimension(np), intent(in) :: lat_f, lon_f
    real(kind=DBL) :: phi_s, lmda_s
    real(kind=DBL), dimension(np) :: phi_f, lmda_f
    real(kind=DBL), dimension(np), intent(out) :: gcirc_res
    real(kind=DBL), dimension(np) :: d_lmda, tp, bt

!    integer :: i   ! for testing

    phi_s = lat_s*pi/180.0
    phi_f = lat_f*pi/180.0
    lmda_s= lon_s*pi/180.0
    lmda_f= lon_f*pi/180.0
    d_lmda= lmda_f - lmda_s

    tp = sqrt( ( cos(phi_f)*sin(d_lmda) )**2 + &
        ( cos(phi_s)*sin(phi_f) - sin(phi_s)*cos(phi_f)*cos(d_lmda))**2)
    bt = sin(phi_s)*sin(phi_f) + cos(phi_s)*cos(phi_f)*cos(d_lmda)
!print*,'In gcirc_dist'
    gcirc_res = (rE + rsat)/1000.0 * atan2(tp,bt) 
!print*,'gcirc_res=', (gcirc_res(i), i=1,20)
!stop
  end subroutine gcirc_dist

!---------------------------------------------------------------

  subroutine my_where(logic_arr, np, idx_arr, ntrue)
! inputs are 
!   logic_arr :   array of T or F
!   np :  size of logic_array
!
! output is index array of all the T and how many (ntrue)
! CLW Aug 24 2011
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
! (i)   Cycle through each Track Number/hemis
! (ii)  Identify and label any strays (off longitude)
! (iii) Sort by great circle distance from an equator point
! Each hemisphere is sorted separately.
! Output is the modified shiftedData structure that has
! - Strays tagged in the typ property
! - Tr0(Nth,Sth), Tr1(Nth,Sth) etc sequence sorted
!  and track order by longitude (for ghost calc)
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
! for sort routine
    real(kind=DBL), dimension(:), allocatable :: T
    integer, dimension(:), allocatable :: Ti

    np = size(struc_data)
    cnt = 1
print*,'Sorting data'

    sort_track_loop: &
    do Tr_num = 0,5
! mask for correct Tr_num and Nth hemis
!print*,'Track Num=',Tr_num

! Get indexes where data has correct TrackNum and Bth hemisphere
      tmplogic_arr = .false.          ! initialise logic array
      tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%GEI_coLat_deg <= 90.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_n)
      allocate(idx_n(iiTrack_n))      ! get mem for these data
      idx_n = tmp_idx(1:iiTrack_n)    ! indexes for these vals
!print*,'idx_n=', (idx_n(i), i=1,iiTrack_n)

! Get indexes for (i) Track (ii) Nth hemis (iii) < 180 in lon 
      tmplogic_arr = .false.
      tmplogic_arr=struc_data(idx_n)%GEI_lon_deg <= 180.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_nlon)
      allocate(idx_nlon(iiTrack_nlon))    ! get mem for these data
      idx_nlon = tmp_idx(1:iiTrack_nlon)  ! indexes for these vals
! These are the indexes of struc for Nth, correct Tr_num and lon <= 180
!print*,'idx_nlon=', (idx_nlon(i), i=1,iiTrack_nlon)

      iiSt_n = maxloc(struc_data(idx_n(idx_nlon))%GEI_coLat_deg)  ! index of equator point
!print*,'iiSt_n=',iiSt_n
      i = idx_n(idx_nlon(iiSt_n(1)))
      st_lat = 90.0 - struc_data(i)%GEI_coLat_deg  ! Lat of this point
      st_lon  = struc_data(i)%GEI_lon_deg          ! Lon of this point
! st_Lat,st_Lon is the coord of the great circle dist reference point
!print*,'st_lat, st_lon=',st_lat,st_lon

      strt_lona(Tr_num+1) = st_lon         ! store start longitude
      trk_order(Tr_num+1) = Tr_num
      allocate(gcDist_n(iiTrack_n))        ! iiTrack_n is the size of idx_n
      call gcirc_dist(st_lat, 90.0-struc_data(idx_n)%GEI_coLat_deg, &
             st_lon, struc_data(idx_n)%GEI_lon_deg, iiTrack_n, gcDist_n)
print*,'iiTrack_n = ',iiTrack_n
!print*,'gcirc_d=', (gcDist_n(i), i=1,iiTrack_n)
!stop
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
!print*,'cnt,iiTrack_n,np ',cnt,iiTrack_n,np
!print*,'cnt:cnt+iiTrack_n-1 ',cnt,cnt+iiTrack_n-1
      cnt = cnt + iiTrack_n
 
      deallocate(idx_nlon)
      deallocate(idx_n)
      deallocate(idx_nsrt)

! ------  do Sth hemisphere ------

      tmplogic_arr = .false.          ! initialise logic array
      tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%GEI_coLat_deg > 90.0
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack_s)
      allocate(idx_s(iiTrack_s))      ! get mem for these data
      idx_s = tmp_idx(1:iiTrack_s)    ! indexes for these vals
!print*,'idx_s=', (idx_s(i), i=1,iiTrack_s)

      allocate(gcDist_s(iiTrack_s))           ! iiTrack_s is the size of idx_s
      call gcirc_dist(st_lat, 90.0-struc_data(idx_s)%GEI_coLat_deg, &
             st_lon, struc_data(idx_s)%GEI_lon_deg, iiTrack_s, gcDist_s)
!print*,'gcirc_d=', (gcDist_s(i), i=1,iiTrack_s)
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
!print*,'cnt,iiTrack_s,np ',cnt,iiTrack_s,np
!print*,'cnt:cnt+iiTrack_s-1 ',cnt,cnt+iiTrack_s-1

      cnt = cnt + iiTrack_s

      deallocate(idx_s)
      deallocate(idx_ssrt)

    enddo sort_track_loop
    struc_data = struc_data( sortII )   ! sort the structure
!do i=1,50
! print*,struc_data(i)%GEI_coLat_deg,struc_data(i)%GEI_lon_deg
!enddo

! sort the track order, returns trk_order
    allocate(T(6))
    allocate(Ti(6))
!print*,'strt_lona = ',strt_lona
    call merge_sort(strt_lona, 6, T, trk_order, Ti)
    deallocate(Ti)
    deallocate(T)

!print*,'Track_order = ',trk_order
!print*,'DONE_sort'
!stop
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

    np = size(struc_data)

print*,'Tagging strays, south = ',south

    track_cyc_loop: &
    do Tr_num = 0,5

! select data from correct track and Nth hemis
      tmplogic_arr = .false.          ! initialise logic array
      if (south .eq. 1) then
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%GEI_coLat_deg >= 90.0
      else
        tmplogic_arr=struc_data%ipln==Tr_num .and. struc_data%GEI_coLat_deg <= 90.0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrack)
      allocate(idx_n(iiTrack))      ! get mem for these data
      idx_n = tmp_idx(1:iiTrack)    ! indexes for these vals
      struc_data(idx_n)%typ = 0     ! initialise type as normal
print*,'Track Num = ',Tr_num
print*,'iiTrack = ',iiTrack
! do i=1,iiTrack
! do i=1,20
!   print*,i, idx_n(i)
!   print*,struc_data(idx_n(i))%px, &
!          struc_data(idx_n(i))%py, &
!          struc_data(idx_n(i))%pz
!   print*,struc_data(idx_n(i))%GEI_coLat_deg, &
!          struc_data(idx_n(i))%GEI_lon_deg
! enddo
! 1st point is on the equator - assumes data are sorted
      uvec(1) = struc_data(idx_n(1))%px
      uvec(2) = struc_data(idx_n(1))%py
      uvec(3) = struc_data(idx_n(1))%pz
!print*,'uvec:',uvec
      uvec = norm_vec(uvec)

!stop
! Get data subset 1/2 about between pole and equator for 2nd point
      tmplogic_arr = .false.
      if (south .eq. 1) then
        tmplogic_arr=struc_data(idx_n)%GEI_coLat_deg > 130.0 &
               .and. struc_data(idx_n)%GEI_coLat_deg < 140.0
      else
        tmplogic_arr=struc_data(idx_n)%GEI_coLat_deg > 40.0 &
               .and. struc_data(idx_n)%GEI_coLat_deg < 50.0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiSel_lat)
      allocate(idx_sel_lat(iiSel_lat))    ! get mem for these data
      idx_sel_lat = tmp_idx(1:iiSel_lat)  ! indexes for these vals

      iiMin_loc = minloc ( &
                  abs( struc_data(idx_n(idx_sel_lat))%GEI_lon_deg &
                    -  struc_data(idx_n(1))%GEI_lon_deg ) )
! point 2 for great circle eqn
      vvec(1) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%px
      vvec(2) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%py
      vvec(3) = struc_data(idx_n(idx_sel_lat(iiMin_loc(1))))%pz
!print*,'vvec:',vvec
      vvec = norm_vec(vvec)

      deallocate(idx_sel_lat)
! normal vector to uvec and vvec
      nvec = cross_p(uvec, vvec)
!print*,'nvec:',nvec
      nvec = norm_vec(nvec)
! normal vector cross uvec
      ncuvec = cross_p(nvec, uvec)
      
      allocate(rarr(iiTrack))
      rarr = sqrt(struc_data(idx_n)%px**2 + &
                  struc_data(idx_n)%py**2 + &
                  struc_data(idx_n)%pz**2)
      allocate(tarr(iiTrack))
      tarr(1)=0.0d0
      do i=2,iiTrack
        dp = uvec(1)*struc_data(idx_n(i))%px + &
             uvec(2)*struc_data(idx_n(i))%py + &
             uvec(3)*struc_data(idx_n(i))%pz
        tarr(i) = acos(dp/rarr(i))
      enddo

!print*,'tarr=', (tarr(i), i=1,20)

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
!print*,'LonVarr=', (LonVarr(i), i=1,20)
      deallocate(LatVarr)  ! used to get LonVarr
      deallocate(rVarr)    ! used to get LonVarr
      deallocate(zarr)
      deallocate(yarr)
      deallocate(xarr)
      deallocate(tarr)
      deallocate(rarr)
!stop
! ------ Check Lon of great circle eqn with the data for strays -----

! check for lon (from eqn) - lon (from data) is > 180 (across 360 deg)
! can't compare arrays with different indexes, so loop'
      do i=1,iiTrack
        lon_diff = abs(LonVarr(i) - struc_data(idx_n(i))%GEI_Lon_deg)
        if (lon_diff > 180.0) then
! if we have lon diff > 180, need to subtract off 360
          if (((360.0 - lon_diff) > 15.0 ) .or. &
             ( lon_diff > 15.0 )) then
            struc_data(idx_n(i))%typ = 1
!print*,'lon_diff=',lon_diff
          endif
        else
          if ( lon_diff > 15.0 ) then
            struc_data(idx_n(i))%typ = 1
!print*,'lon_diff=',lon_diff
          endif
        endif
      enddo     ! search track data points

      deallocate(LonVarr)
      deallocate(idx_n)

    enddo track_cyc_loop

  end subroutine tag_lon_strays
!
! ----------------------------------------------------------------------
!
  subroutine ampFit_ghost(struc_data, comb_data, clat_lim, trk_order)
! Given:
!   shiftedData structure,
!   clat_lim, the pole limit for ghosting
!   trk_order (integer arry of 6 of the track order)
! this routine:
! (i)   Cycles through adjacent pair tracks
! (ii)  Identify  coLat limit and add data between lon
! Assumes the data for each track has been sorted by sort_struc
! - Tr0(Nth,Sth), Tr1(Nth,Sth) etc sequence sorted
!
! CLW - Aug 2011
!
    implicit none

    type(ampData), intent(in out) :: struc_data(:)
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

    np = size(struc_data)

print*,'In ampFit_ghost, np=',np
print*,'clat_lim=',clat_lim
    iiex_gh = 0
    cnt = 1
    gh_fname='ex_ghosts.dat'
    tmp_open = 0

! Calc min number of ghost points (excludes extra ghosts)
    tmplogic_arr = .false.        ! initialise logic array
    if (clat_lim .lt. 90.0) then  ! Nth hemis data
      tmplogic_arr = struc_data%GEI_coLat_deg <= clat_lim .and. &
                     struc_data%typ==0
    else                          ! Sth hemis data
      tmplogic_arr = struc_data%GEI_coLat_deg >= clat_lim .and. &
                     struc_data%typ==0
    endif
    call my_where(tmplogic_arr, np, tmp_idx, iighost)
    allocate(gh_data(iighost)) ! get mem for ghost_data struc

print*,'iighost=',iighost

    track_pair_loop: &
    do Tr_num = 0,5

! mask for correct Tr_num and hemis
      tmplogic_arr = .false.        ! initialise logic array
      if (clat_lim .lt. 90.0) then  ! Nth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(Tr_num+1) .and. &
                       struc_data%GEI_coLat_deg <= clat_lim .and. &
                       struc_data%typ==0
      else                          ! Sth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(Tr_num+1) .and. &
                       struc_data%GEI_coLat_deg >= clat_lim .and. &
                       struc_data%typ==0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_1)
      allocate(idx_Trk1(iiTrk_1))   ! get mem for indexes for Trk_1
      idx_Trk1 = tmp_idx(1:iiTrk_1) ! indexes for these vals

print*,'Track Num:',Tr_num
print*,'iiTrk_1=',iiTrk_1

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
                       struc_data%GEI_coLat_deg <= clat_lim .and. &
                       struc_data%typ==0
      else                          ! Sth hemis data
        tmplogic_arr = struc_data%ipln==trk_order(nxt_trk) .and. &
                       struc_data%GEI_coLat_deg >= clat_lim .and. &
                       struc_data%typ==0
      endif
      call my_where(tmplogic_arr, np, tmp_idx, iiTrk_2)
      allocate(idx_Trk2(iiTrk_2))   ! get mem for indexes for Trk_2 data
      idx_Trk2 = tmp_idx(1:iiTrk_2) ! indexes for these vals
print*,'iiTrk_2 = ',iiTrk_2

      allocate(gcDist(iiTrk_2))     ! mem for great circ dist values
! if iiTrk2 = 0 then no data for ghost - error

      point_ex_loop: &
      do jj=1,iiTrk_1                ! loop thru all points on Trk_1
        st_lat = struc_data(idx_Trk1(jj))%GEI_coLat_deg
        st_lon = struc_data(idx_Trk1(jj))%GEI_lon_deg
        call gcirc_dist(90.0-st_lat, 90.0-struc_data(idx_Trk2)%GEI_coLat_deg, &
             st_lon, struc_data(idx_Trk2)%GEI_lon_deg, iiTrk_2, gcDist)

! iimn is the index of idx_Trk2 which has min dist
!  between Trk2 and the jj point of Trk1
        iimn = minloc(gcDist)        ! index of idx_Trk2 of min_dist
!print*,'iimn = ',iimn
        
        if (clat_lim .gt. 90.0) then ! Sth hemis data
          st_lat=180.0-st_lat
        endif
        xs = st_lat*cos(st_lon*pi/180.0)
        ys = st_lat*sin(st_lon*pi/180.0)
        e_lat=struc_data(idx_Trk2(iimn(1)))%GEI_coLat_deg
        if (clat_lim .gt. 90.0) then
          e_lat=180.0-e_lat
        endif
        e_lon=struc_data(idx_Trk2(iimn(1)))%GEI_lon_deg
        Nyq_Ok = 1

! Check if we need extra ghosts -> need to add 2 ghost points
        if ( abs(e_lon-st_lon) .ge. 120.0 .and. &
             abs(e_lon-st_lon) .lt. 180.0) then
          Nyq_Ok = 0
          iiex_gh = iiex_gh + 1
          if (tmp_open .eq. 0) then   ! 1st extra ghost
! open the extra ghist temp data file here
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
          
          gh_data(cnt)%GEI_R_km = 2.0/3.0*struc_data(idx1)%GEI_R_km + &
                                1.0/3.0*struc_data(idx2)%GEI_R_km
          exgh_data%GEI_R_km = 1.0/3.0*struc_data(idx1)%GEI_R_km + &
                                2.0/3.0*struc_data(idx2)%GEI_R_km

          if (clat_lim .gt. 90.0) then   ! south hemis
            gh_data(cnt)%GEI_coLat_deg = 180.0 - sqrt(xg1**2 + yg1**2)
            exgh_data%GEI_coLat_deg = 180.0 - sqrt(xg2**2 + yg2**2)
          else                           ! North hemisphere
            gh_data(cnt)%GEI_coLat_deg = sqrt(xg1**2 + yg1**2)
            exgh_data%GEI_coLat_deg = sqrt(xg2**2 + yg2**2)
          endif

          gh_data(cnt)%GEI_coLat_rad = gh_data(cnt)%GEI_coLat_deg*pi/180.0
          exgh_data%GEI_coLat_rad = exgh_data%GEI_coLat_deg*pi/180.0

          gh_data(cnt)%GEI_lon_rad = atan2(yg1,xg1)
          exgh_data%GEI_lon_rad = atan2(yg2,xg2)
          gh_data(cnt)%GEI_lon_deg = atan2(yg1,xg1)*180.0/pi
          exgh_data%GEI_lon_deg = atan2(yg2,xg2)*180.0/pi
! Convert spherical coords to x,y,z
          call sphcar_08(gh_data(cnt)%GEI_R_km, &
                         gh_data(cnt)%GEI_coLat_rad, &
                         gh_data(cnt)%GEI_lon_rad, &
                         xv, yv, zv, 1)
          gh_data(cnt)%px = xv
          gh_data(cnt)%py = yv
          gh_data(cnt)%pz = zv
          call sphcar_08(exgh_data%GEI_R_km, &
                         exgh_data%GEI_coLat_rad, &
                         exgh_data%GEI_lon_rad, &
                         xv, yv, zv, 1)
          exgh_data%px = xv
          exgh_data%py = yv
          exgh_data%pz = zv

          gh_data(cnt)%br_GEI = 2.0/3.0*struc_data(idx1)%br_GEI + &
                                1.0/3.0*struc_data(idx2)%br_GEI
          exgh_data%br_GEI = 1.0/3.0*struc_data(idx1)%br_GEI + &
                                2.0/3.0*struc_data(idx2)%br_GEI

          gh_data(cnt)%btheta_GEI = 2.0/3.0*struc_data(idx1)%btheta_GEI + &
                                1.0/3.0*struc_data(idx2)%btheta_GEI
          exgh_data%btheta_GEI = 1.0/3.0*struc_data(idx1)%btheta_GEI + &
                                2.0/3.0*struc_data(idx2)%btheta_GEI

          gh_data(cnt)%bphi_GEI = 2.0/3.0*struc_data(idx1)%bphi_GEI + &
                                1.0/3.0*struc_data(idx2)%bphi_GEI
          exgh_data%bphi_GEI = 1.0/3.0*struc_data(idx1)%bphi_GEI + &
                                2.0/3.0*struc_data(idx2)%bphi_GEI
! Convert spherical dB vector to x,y,z
          call bspcar_08(gh_data(cnt)%GEI_coLat_rad, &
                         gh_data(cnt)%GEI_lon_rad, &
                         gh_data(cnt)%br_GEI, &
                         gh_data(cnt)%btheta_GEI, &
                         gh_data(cnt)%bphi_GEI, &
                         xv, yv, zv)
          gh_data(cnt)%dbx = xv
          gh_data(cnt)%dby = yv
          gh_data(cnt)%dbz = zv
          call bspcar_08(exgh_data%GEI_coLat_rad, &
                         exgh_data%GEI_lon_rad, &
                         exgh_data%br_GEI, &
                         exgh_data%btheta_GEI, &
                         exgh_data%bphi_GEI, &
                         xv, yv, zv)
          exgh_data%dbx = xv
          exgh_data%dby = yv
          exgh_data%dbz = zv
! write exgh_data to temp file
          if (exgh_data%GEI_coLat_deg > 1.0 .and. &
              exgh_data%GEI_coLat_deg < 179.0) then
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
          
          gh_data(cnt)%GEI_R_km = (struc_data(idx1)%GEI_R_km + &
                                   struc_data(idx2)%GEI_R_km)/2.0

          if (clat_lim .gt. 90.0) then   ! south hemis
            gh_data(cnt)%GEI_coLat_deg = 180.0 - sqrt(xg**2 + yg**2)
          else                           ! North hemisphere
            gh_data(cnt)%GEI_coLat_deg = sqrt(xg**2 + yg**2)
          endif

          gh_data(cnt)%GEI_coLat_rad = gh_data(cnt)%GEI_coLat_deg*pi/180.0

          gh_data(cnt)%GEI_lon_rad = atan2(yg,xg)
          gh_data(cnt)%GEI_lon_deg = atan2(yg,xg)*180.0/pi

          call sphcar_08(gh_data(cnt)%GEI_R_km, &
                         gh_data(cnt)%GEI_coLat_rad, &
                         gh_data(cnt)%GEI_lon_rad, &
                         xv, yv, zv, 1)
          gh_data(cnt)%px = xv
          gh_data(cnt)%py = yv
          gh_data(cnt)%pz = zv

          gh_data(cnt)%br_GEI = struc_data(idx1)%br_GEI + &
                                struc_data(idx2)%br_GEI

          gh_data(cnt)%btheta_GEI = (struc_data(idx1)%btheta_GEI + &
                                     struc_data(idx2)%btheta_GEI)/2.0

          gh_data(cnt)%bphi_GEI = (struc_data(idx1)%bphi_GEI + &
                                   struc_data(idx2)%bphi_GEI)/2.0

          call bspcar_08(gh_data(cnt)%GEI_coLat_rad, &
                         gh_data(cnt)%GEI_lon_rad, &
                         gh_data(cnt)%br_GEI, &
                         gh_data(cnt)%btheta_GEI, &
                         gh_data(cnt)%bphi_GEI, &
                         xv, yv, zv)
          gh_data(cnt)%dbx = xv
          gh_data(cnt)%dby = yv
          gh_data(cnt)%dbz = zv

        endif    ! one ghost point

        if (gh_data(cnt)%GEI_coLat_deg > 1.0 .and. &
            gh_data(cnt)%GEI_coLat_deg < 179.0) then
          cnt = cnt + 1
        endif

      enddo point_ex_loop

      deallocate(gcDist)
      deallocate(idx_Trk2)
      deallocate(idx_Trk1)
!print*,'gh_lon data=',(gh_data(i)%GEI_lon_deg, i=1,cnt-1)
!stop
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

  end subroutine ampFit_ghost

end module ampFit_sort
