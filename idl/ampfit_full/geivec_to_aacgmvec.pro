pro geiVec_to_aacgmVec, $
		gei_R_km, gei_coLat_rad, gei_lon_rad, $
		dBTheta_gei, dbPhi_gei, $
		dBTheta_aacgm_gth = dBTheta_aacgm_gth, $
		dBPhi_aacgm_gth = dBPhi_aacgm_gth, $
		dBTheta_aacgm_gph = dBTheta_aacgm_gph, $
		dBPhi_aacgm_gph = dBPhi_aacgm_gph, $
		aacgm_lat_deg = aacgm_lat_deg, $
		aacgm_lon_deg = aacgm_lon_deg

	@amp_fit_constants

	; Rotate the final dB vectors into AACGM
	; For this we need the GEO coords of the
	; grid and at this point we only have the
	; AACGM and GEI grid coords.
	; --------------------------------------

	; Get grid GEO coords
	; This should probably be done in the grid
	; generation rouinte. Will do in future.
	; ----------------------------------------

    geopack_sphcar, $
			gei_R_km, gei_coLat_rad, gei_lon_rad,$
           	gei_x, gei_y, gei_z, $
			/to_rect

    geoPack_conv_coord, $
			gei_x, gei_y, gei_z, $
    		geo_x, geo_y, geo_z, $
      		/from_gei, /to_geo

    geopack_sphcar, $
			geo_x, geo_y, geo_z, $
   			geo_rIrid_km, geo_coLat_rad, geo_lon_rad, $
			/to_sphere


    geo_lat_deg	= 90.0 - geo_coLat_rad * !radeg
    geo_lon_deg	= geo_lon_rad * !radeg
    aacgm_yr=2005    ; latest coeffs

    aacgm_conv_vec, $
		geo_lat_deg, geo_lon_deg, geo_rIrid_km-rE_km, $
		dBTheta_GEI, dBPhi_GEI, $
		aacgm_lat_deg, aacgm_lon_deg, $
		dBTheta_aacgm_gth, dBPhi_aacgm_gth, $
		dBTheta_aacgm_gph, dBphi_aacgm_gph, $
		err, /to_aacgm

    ;vec_geo2aacgm,path,aacgm_yr,geog_rIrid_km-rE_km,$
    ;                   gLat_a,gLon_a,dBTheta_GEI_grid,dBPhi_GEI_grid,$
    ;                   mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph


end
