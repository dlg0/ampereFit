module ampFit_geivec_aacgmvec

 use ISO_C_BINDING
 use f_aacgm
 use constants
 use ampFit_sort

 implicit none

contains

subroutine geisph_to_aacgmvec (r_km, thet_rad, phi_rad, &
           db_thet, db_phi, &
           dbthet_aacgm_th, dbphi_aacgm_th, &
           dbthet_aacgm_ph, dbphi_aacgm_ph, &
           aacgm_lat_deg, aacgm_lon_deg )
! r_km : Re + 780 km 
! thet_rad : GEI coLat in radians
! phi_rad : GEI longitude in radians
! dB_thet : mag field in N-S
! dB_phi : mag field in E-W
!
! Colin Waters
! Sept 2011

  implicit none

  real(kind=DBL), intent(in) :: r_km, thet_rad, phi_rad, &
                                db_thet, db_phi
  real(kind=DBL), intent(out) :: dbthet_aacgm_th, dbphi_aacgm_th, &
                                 dbthet_aacgm_ph, dbphi_aacgm_ph, &
                                 aacgm_lat_deg, aacgm_lon_deg

  real(kind=DBL), dimension(3) :: gei_xyz, geo_xyz, &
                                  geo_xyz_th, geo_xyz_ph, &
                                  geo_rtp_th, geo_rtp_ph, &
                         aacgm_xyz, aacgm_xyz_th, aacgm_xyz_ph 

  real(kind=DBL) :: geo_r_km, geo_lat_deg, geo_lon_deg, geo_clat_rad, geo_lon_rad
  real(kind=DBL) :: aacgm_hgt, aacgm_colat_rad, aacgm_lon_rad
  real(kind=DBL), target :: mlat, mlon, mr
  real(kind=DBL) :: mclat_rad, mlon_rad
  real(kind=DBL) :: vx, vy, vz, br,bth,bph
  real(kind=DBL) :: rg_th, glat_th, cglat_th, glon_th

  real(kind=DBL), dimension(3) :: mxyz_r, mxyz_ruv, &
                                  mxyz_th, mxyz_thuv, &
                                  mxyz_ph, mxyz_phuv, &
                                  mxyz_ph_gth, mxyz_th_gph, &
                                  mvec_th, mvec_ph
  integer :: s, flg

! GEI_spherical to GEI x,y,z
  call sphcar_08(r_km, thet_rad, phi_rad, &
                 gei_xyz(1), gei_xyz(2), gei_xyz(3), 1)

! GEI_x,y,z to GEO_x,y,z
  call geigeo_08(gei_xyz(1), gei_xyz(2), gei_xyz(3), &
                 geo_xyz(1), geo_xyz(2), geo_xyz(3), 1)

! GEO_x,y,z to GEO_spherical
  call sphcar_08(geo_r_km, geo_clat_rad, geo_lon_rad, &
                 geo_xyz(1), geo_xyz(2), geo_xyz(3), -1)

! GEO_spherical to aacgm coords
  geo_lat_deg = 90.0 - geo_clat_rad*180.0/pi
  geo_lon_deg = geo_lon_rad*180.0/pi
  aacgm_hgt = rsat/1000.0
  flg = 0
  s = f_AACGMConvert(geo_lat_deg, geo_lon_deg, aacgm_hgt, &
        C_LOC(mlat), C_LOC(mlon), C_LOC(mr), flg)

  aacgm_lat_deg = mlat
  if (mlon .lt. 0.0) then
    mlon = mlon + 360.0
  endif
  aacgm_lon_deg = mlon
  aacgm_lon_rad = aacgm_lon_deg*pi/180.0

! AACGM_spherical to x,y,z
  aacgm_colat_rad = (90.0 - aacgm_lat_deg)*pi/180.0
  call sphcar_08(r_km, aacgm_colat_rad, aacgm_lon_rad, &
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

  mclat_rad = (90.0 - mlat)*pi/180.0
  if (mlon .lt. 0.0) then
    mlon = mlon + 360.0
  endif
  mlon_rad = mlon*pi/180.0

! Get x,y,z coords of AACGM + dthet
  call sphcar_08(r_km, mclat_rad, mlon_rad, &
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

  mclat_rad = (90.0 - mlat)*pi/180.0
  if (mlon .lt. 0.0) then
    mlon = mlon + 360.0
  endif
  mlon_rad = mlon*pi/180.0

! Get x,y,z coords of AACGM + dthet
  call sphcar_08(r_km, mclat_rad, mlon_rad, &
                 aacgm_xyz_ph(1), &
                 aacgm_xyz_ph(2), &
                 aacgm_xyz_ph(3), 1)

!  - - - - - - end of dph shift - - - - - - - - - - 

! calc AACGM radial unit vector
!  mxyz_r(1) = xm
!  mxyz_r(2) = ym
!  mxyz_r(3) = zm
  mxyz_ruv = norm_vec(aacgm_xyz)

! calc AACGM(x,y,z) unit vector for a GEO dth shift
  mxyz_th = aacgm_xyz_th - aacgm_xyz

!  mxyz_th(1) = xm_th - xm
!  mxyz_th(2) = ym_th - ym
!  mxyz_th(3) = zm_th - zm
  mxyz_thuv = norm_vec(mxyz_th)

! calc AACGM(x,y,z) unit vector for a GEO dph shift
  mxyz_ph = aacgm_xyz_ph - aacgm_xyz
!  mxyz_ph(1) = xm_ph - xm
!  mxyz_ph(2) = ym_ph - ym
!  mxyz_ph(3) = zm_ph - zm
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
  
  dbthet_aacgm_th = db_thet*mvec_th(2) + db_phi*mvec_ph(2)
  dbphi_aacgm_th  = db_thet*mvec_th(3) + db_phi*mvec_ph(3)

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
  
  dbthet_aacgm_ph = db_thet*mvec_th(2) + db_phi*mvec_ph(2)
  dbphi_aacgm_ph  = db_thet*mvec_th(3) + db_phi*mvec_ph(3)

  end subroutine geisph_to_aacgmvec

end module ampFit_geivec_aacgmvec
