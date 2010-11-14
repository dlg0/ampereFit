;
; Generate uniform AACGM grid
; DLG & CLW : Dec, 2009
;
pro aacgm_grid, calc_coLat, south, $
  aacgmGrid_R_km = aacgmGrid_R_km, $
	aacgmGrid_coLat_deg = aacgmGrid_coLat_deg, $
	aacgmGrid_lon_deg = aacgmGrid_lon_deg, $
  geiGrid_R_km = geiGrid_R_km, $
	geiGrid_coLat_rad = geiGrid_coLat_rad, $
	geiGrid_lon_rad = geiGrid_lon_rad, $
  xGrid_GEI=xGrid_GEI, $
  yGrid_GEI=yGrid_GEI, $
  zGrid_GEI=zGrid_GEI, $
	nLat = nLat, nLon = nLon, $
	year = year, yrSec = yrSec, $
	mltShift = mltShift, $
	epoch = epoch

	; Default keyword variables
	; -------------------------

	if (not keyword_set(nLat)) then nLat = calc_coLat
	if (not keyword_set(nLon)) then nLon = 24

	; Iridium altitude
	; ----------------

	Re_km = 6371.0
	R_km = 780.0 + Re_km

	; Generate uniform AACGM grid
	; ---------------------------

	if south eq 0 then $
		aacgm_coLat_deg = (fIndGen(nLat)+1)/nLat*calc_coLat $
	else $
	   aacgm_coLat_deg = 180.0-calc_coLat+(fIndGen(nLat)+1)/nLat*calc_coLat

	mlon0 = aacgm_mlong(year,yrSec,0.0)    ; Mag longitude for 0 MLT
	aacgm_lon_deg = (fIndGen(nLon))*(360.0-360.0/nLon)/(nLon-1)+mlon0
	aacgm_lon_deg = aacgm_lon_deg mod 360.0
	aacgmGrid_coLat_deg = rebin(aacgm_coLat_deg, nLat, nLon )
	aacgmGrid_lon_deg = transpose (rebin(aacgm_lon_deg, nLon, nLat))

	mltShift = aacgm_mlt(year, yrSec, 0.0 )  ; MLT at 0 MLon
	aacgmGrid_R_km = aacgmGrid_lon_deg * 0 + R_km

	; Convert AACGM grid to GEO
	; -------------------------

	aacgm_conv_coord, 90.0-aacgmGrid_coLat_deg, aacgmGrid_lon_deg, $
		aacgmGrid_R_km-Re_km, geogGrid_lat_deg, geogGrid_lon_deg, err, /to_geo

	; Convert to GEO XYZ locations (input theta,phi is radians)
	; ------------------------------------------------------

	geoPack_sphCar,aacgmGrid_R_km,(90.0-geogGrid_lat_deg)*!dtor, geogGrid_lon_deg*!dtor, $
		xGrid_GEO, yGrid_GEO, zGrid_GEO, /to_rect

	; Convert GEO XYZ to GEI XYZ
	; --------------------------

	geoPack_conv_coord, xGrid_GEO, yGrid_GEO, zGrid_GEO, $
		xGrid_GEI, yGrid_GEI, zGrid_GEI,/from_geo,/to_gei, $
		epoch = fltArr ( n_elements ( xGrid_GEO[*] ) ) + epoch

	; Convert GEI XYZ locations to GEI (r,thet,phi)
	; ------------------------------------------

	geoPack_sphCar, xGrid_GEI, yGrid_GEI, zGrid_GEI, $
		geiGrid_R_km, geiGrid_coLat_rad, geiGrid_lon_rad, /to_sphere
end
