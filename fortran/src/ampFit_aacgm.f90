module ampFit_aacgm
!
! create_aacgm_grid
! subroutine geosph_to_aacgmvec -> Converts db vectors from GEO to AACGM, for one point
!
! Ver : 201202
!
! C.L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! November, 2011
!
 use ISO_C_BINDING
 use f_aacgm
 use constants
 use ampFit_sort

 implicit none

contains

subroutine create_aacgm_grid (gridAACGM, gridGEO, &
           nLatGrid, nLonGrid, latStep, lonStep, south )

  use ampFit_rotate

  implicit none

  type(ampData), intent(inout) :: gridAACGM(:), gridGEO(:)
  real(kind=DBL), intent(in) :: latStep, lonStep
  integer, intent(in) :: nLatGrid, nLonGrid, south

  real(kind=DBL), target :: geoLat_deg, geoLon_deg, geoHgt_km
  real(kind=DBL) :: aacgmHgt_km, aacgmcolat_deg, aacgmlon_deg
  integer :: i, j, s, flg, idx

  aacgmHgt_km = rSat * 1e-3

  flg = 1 ! 0 -> to aacgm, 1 -> geographic
  do i=1,nLatGrid
    if (south .eq. 0) then
      aacgmcoLat_deg = i * latStep
    else
      aacgmcoLat_deg = 180.0 - (i * latStep)
    endif
    do j=1,nLonGrid
      aacgmLon_deg = j * lonStep

      s = f_AACGMConvert(90.0-aacgmcoLat_deg, &
              aacgmLon_deg, aacgmHgt_km,&
              C_LOC(geoLat_deg), &
              C_LOC(geoLon_deg), &
              C_LOC(geoHgt_km), flg)

      if(s/=0) then
        err_stat = 3
        err_msg = 'AACGM conv error in create_aacgm_grid'
        write(*,*) err_msg
        write(*,*) 'Error in f_AACGMConvert: aacgm->GEO'
        write(*,*) 'Inputs were: ', 90.0-aacgmcoLat_deg, &
            aacgmLon_deg, aacgmHgt_km, flg
         return
      endif

      idx = (i-1)*nLonGrid+j

      gridAACGM(idx)%bR = 0.0
      gridAACGM(idx)%bT = 0.0
      gridAACGM(idx)%bP = 0.0
 
      gridAACGM(idx)%T = aacgmcoLat_deg* degToRad
      gridAACGM(idx)%P = aacgmLon_deg * degToRad
      gridAACGM(idx)%R = aacgmHgt_km + rE*1e-3

      gridGEO(idx)%bR = 0.0
      gridGEO(idx)%bT = 0.0
      gridGEO(idx)%bP = 0.0

      gridGEO(idx)%T = (90.0-geoLat_deg)* degToRad
      gridGEO(idx)%P = geoLon_deg * degToRad
      gridGEO(idx)%R = geoHgt_km + rE*1e-3

    enddo
  enddo

  j = 1
  call rtp_xyz_coord( gridGEO, j)
  call rtp_xyz_coord( gridAACGM, j ) 

end subroutine create_aacgm_grid
!
! ----------------------------------------------------------------------
!
subroutine geosph_to_aacgmvec (dataGEO, &
             aacgm_clat_rad, aacgm_lon_rad, &
             dbthet_aacgm_th, dbphi_aacgm_th, &
             dbthet_aacgm_ph, dbphi_aacgm_ph)
!
! C.L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Ver : 201202
!
! geo_clat_rad : GEO coLat in radians
! geo_lon_rad : GEO longitude in radians
! geo_dB_thet : mag field in N-S
! geo_dB_phi : mag field in E-W
!
  implicit none

  type(ampData), intent(in) :: dataGEO(:)

  real(kind=DBL), intent(out) :: aacgm_clat_rad(:), &
                                 aacgm_lon_rad(:), &
                                dbthet_aacgm_th(:), &
                                dbphi_aacgm_th(:), &
                                dbthet_aacgm_ph(:), &
                                dbphi_aacgm_ph(:)

  real(kind=DBL) :: geo_r_km, geo_clat_rad, geo_lon_rad, &
                    geo_db_thet, geo_db_phi

  real(kind=DBL), dimension(3) :: geo_xyz, &
                                  geo_xyz_th, geo_xyz_ph, &
                                  geo_rtp_th, geo_rtp_ph, &
                         aacgm_xyz, aacgm_xyz_th, aacgm_xyz_ph 

  real(kind=DBL) :: geo_lat_deg, geo_lon_deg
  real(kind=DBL) :: aacgm_hgt
  real(kind=DBL), target :: mlat, mlon, mr
  real(kind=DBL) :: mclat_rad, mlon_rad
  real(kind=DBL) :: vx, vy, vz, br,bth,bph
  real(kind=DBL) :: rg_th, glat_th, cglat_th, glon_th

  real(kind=DBL), dimension(3) :: mxyz_r, mxyz_ruv, &
                                  mxyz_th, mxyz_thuv, &
                                  mxyz_ph, mxyz_phuv, &
                                  mxyz_ph_gth, mxyz_th_gph, &
                                  mvec_th, mvec_ph
  integer :: i, s, flg, np

  np = size(dataGEO)

! Rotate vectors from GEO to AACGM
  do i=1,np

    geo_r_km = dataGEO(i)%R
    geo_clat_rad = dataGEO(i)%T
    geo_lon_rad = dataGEO(i)%P
    geo_db_thet = dataGEO(i)%bT
    geo_db_phi = dataGEO(i)%bP

! GEO_spherical to GEO x,y,z
    call sphcar_08(geo_r_km, geo_clat_rad, geo_lon_rad, &
                   geo_xyz(1), geo_xyz(2), geo_xyz(3), 1)

! GEO_spherical to aacgm coords
    geo_lat_deg = 90.0 - geo_clat_rad*180.0/pi
    geo_lon_deg = geo_lon_rad*180.0/pi
    aacgm_hgt = rsat/1000.0
    flg = 0
  
    s = f_AACGMConvert(geo_lat_deg, geo_lon_deg, aacgm_hgt, &
                       C_LOC(mlat), C_LOC(mlon), C_LOC(mr), flg)

    if(s/=0) then
      err_stat = 4
      err_msg = 'AACGM conv error in geosph_to_aacgmvec'
      write(*,*) err_msg
      write(*,*) 'Error in f_AACGMConvert: GEO->AACGM'
      write(*,*) 'Inputs were: ', geo_lat_deg, &
          geo_lon_deg, aacgm_hgt
       return
    endif

    if (mlon .lt. 0.0) then
      mlon = mlon + 360.0
    endif
    aacgm_lon_rad(i) = mlon*pi/180.0

! AACGM_spherical to x,y,z
    aacgm_clat_rad(i) = (90.0 - mlat)*pi/180.0
    call sphcar_08(geo_r_km, aacgm_clat_rad(i), aacgm_lon_rad(i), &
                   aacgm_xyz(1), aacgm_xyz(2), aacgm_xyz(3), 1)

! - - - - - - - - - del_th in GEO - - - - - - - - - -
    br=0.0
    bth=1.0
    bph=0.0
    call bspcar_08(geo_clat_rad, geo_lon_rad, &
                   br, bth, bph, vx, vy, vz)
    geo_xyz_th(1) = geo_xyz(1) + vx
    geo_xyz_th(2) = geo_xyz(2) + vy
    geo_xyz_th(3) = geo_xyz(3) + vz

! Conv x,y,z to spherical
    call sphcar_08(geo_rtp_th(1), geo_rtp_th(2), geo_rtp_th(3), &
                   geo_xyz_th(1), geo_xyz_th(2), geo_xyz_th(3), -1)

    geo_rtp_th(2) = 90.0 - geo_rtp_th(2)*180.0/pi
    geo_rtp_th(3) = geo_rtp_th(3)*180.0/pi

! convert to AACGM (from GEO_th)
    flg = 0
    s = f_AACGMConvert(geo_rtp_th(2), geo_rtp_th(3), aacgm_hgt, &
                           C_LOC(mlat), C_LOC(mlon), C_LOC(mr), flg)

    if(s/=0) then
      err_stat = 4
      err_msg = 'AACGM conv error in geosph_to_aacgmvec'
      write(*,*) err_msg
      write(*,*) 'Error in f_AACGMConvert: GEO->AACGM'
      write(*,*) 'Inputs were: ', geo_rtp_th(2), &
          geo_rtp_th(3), aacgm_hgt
       return
    endif
    mclat_rad = (90.0 - mlat)*pi/180.0
    if (mlon .lt. 0.0) then
      mlon = mlon + 360.0
    endif
    mlon_rad = mlon*pi/180.0

! Get x,y,z coords of AACGM + dthet
    call sphcar_08(geo_r_km, mclat_rad, mlon_rad, &
                   aacgm_xyz_th(1), &
                   aacgm_xyz_th(2), &
                   aacgm_xyz_th(3), 1)

!  - - - - - - end of dth shift - - - - - - - - - - 

! - - - - -  - del_phi in GEO  - - - - - - - - - -
    br=0.0
    bth=0.0
    bph=1.0
    call bspcar_08(geo_clat_rad, geo_lon_rad, &
                   br, bth, bph, vx, vy, vz)
    geo_xyz_ph(1) = geo_xyz(1) + vx
    geo_xyz_ph(2) = geo_xyz(2) + vy
    geo_xyz_ph(3) = geo_xyz(3) + vz

! Conv x,y,z to spherical
    call sphcar_08(geo_rtp_ph(1), geo_rtp_ph(2), geo_rtp_ph(3), &
                   geo_xyz_ph(1), geo_xyz_ph(2), geo_xyz_ph(3), -1)

    geo_rtp_ph(2) = 90.0 - geo_rtp_ph(2)*180.0/pi
    geo_rtp_ph(3) = geo_rtp_ph(3)*180.0/pi

! convert to AACGM (from GEO_ph)
    flg = 0
    s = f_AACGMConvert(geo_rtp_ph(2), geo_rtp_ph(3), aacgm_hgt, &
          C_LOC(mlat), C_LOC(mlon), C_LOC(mr), flg)

    if(s/=0) then
      err_stat = 4
      err_msg = 'AACGM conv error in geosph_to_aacgmvec'
      write(*,*) err_msg
      write(*,*) 'Error in f_AACGMConvert: GEO->AACGM'
      write(*,*) 'Inputs were: ', geo_rtp_ph(2), &
          geo_rtp_ph(3), aacgm_hgt
       return
    endif

    mclat_rad = (90.0 - mlat)*pi/180.0
    if (mlon .lt. 0.0) then
      mlon = mlon + 360.0
    endif
    mlon_rad = mlon*pi/180.0

! Get x,y,z coords of AACGM + dthet
    call sphcar_08(geo_r_km, mclat_rad, mlon_rad, &
                   aacgm_xyz_ph(1), &
                   aacgm_xyz_ph(2), &
                   aacgm_xyz_ph(3), 1)

!  - - - - - - end of dph shift - - - - - - - - - - 

! calc AACGM radial unit vector
    mxyz_ruv = norm_vec(aacgm_xyz)

!   calc AACGM(x,y,z) unit vector for a GEO dth shift
    mxyz_th = aacgm_xyz_th - aacgm_xyz

    mxyz_thuv = norm_vec(mxyz_th)

! calc AACGM(x,y,z) unit vector for a GEO dph shift
    mxyz_ph = aacgm_xyz_ph - aacgm_xyz
    mxyz_phuv = norm_vec(mxyz_ph)

! calc cross products
    mxyz_ph_gth = cross_p(mxyz_ruv, mxyz_thuv)

! For a GEO dth shift -> AACGM dth
    call bcarsp_08(aacgm_xyz(1), aacgm_xyz(2), aacgm_xyz(3), &
                   mxyz_thuv(1), mxyz_thuv(2), mxyz_thuv(3), &
                   mvec_th(1), mvec_th(2), mvec_th(3))

! For a GEO dth shift -> AACGM dph
    call bcarsp_08(aacgm_xyz(1), aacgm_xyz(2), aacgm_xyz(3), &
                   mxyz_ph_gth(1), mxyz_ph_gth(2), mxyz_ph_gth(3), &
                   mvec_ph(1), mvec_ph(2), mvec_ph(3))
  
    dbthet_aacgm_th(i) = geo_db_thet*mvec_th(2) + geo_db_phi*mvec_ph(2)
    dbphi_aacgm_th(i)  = geo_db_thet*mvec_th(3) + geo_db_phi*mvec_ph(3)

! Now do the PHI components
! calc cross products
    mxyz_th_gph = cross_p(mxyz_phuv, mxyz_ruv)

! For a GEO dph shift -> AACGM dth
    call bcarsp_08(aacgm_xyz(1), aacgm_xyz(2), aacgm_xyz(3), &
                   mxyz_th_gph(1), mxyz_th_gph(2), mxyz_th_gph(3), &
                   mvec_th(1), mvec_th(2), mvec_th(3))

! For a GEO dph shift -> AACGM dph
    call bcarsp_08(aacgm_xyz(1), aacgm_xyz(2), aacgm_xyz(3), &
                   mxyz_phuv(1), mxyz_phuv(2), mxyz_phuv(3), &
                   mvec_ph(1), mvec_ph(2), mvec_ph(3))
  
    dbthet_aacgm_ph(i) = geo_db_thet*mvec_th(2) + geo_db_phi*mvec_ph(2)
    dbphi_aacgm_ph(i)  = geo_db_thet*mvec_th(3) + geo_db_phi*mvec_ph(3)

  enddo

end subroutine geosph_to_aacgmvec

end module ampFit_aacgm
!
! =============================================================================
