

pro rotate_gei_to_aacgm, gei_R_km, gei_coLat_rad, gei_lon_rad, $
		aacgm_R_km, aacgm_coLat_rad, aacgm_lon_rad, $
		gei_dbTh, gei_dbPh, $
		aacgm_dbTh = aacgm_dbTh, $
		aacgm_dbPh = aacgm_dbPh, $
		year = year, $
		epoch = epoch

	if ( not keyword_set ( year ) ) then year = 2000

	;	OK, method 2 ... we need the r & theta unit vectors
	;	in terms of their x and y comonents for both GEI and
	;	AACGM at each location. Here we use polar coordinates 
	;	for the conversion, not sure how correct that is but
	;	that's how it's been done in the past.

	dR	= 0.01

	;	get the x/y coords

	GEI_x	= gei_colat_rad * !radeg * cos ( gei_lon_rad )
	GEI_y	= gei_colat_rad * !radeg * sin ( gei_lon_rad )

	;	get GEI r & theta unit vectors in terms of x and y components

		GEI_x_R	= ( gei_colat_rad * !radeg + dR ) * cos ( gei_lon_rad )
		GEI_y_R	= ( gei_colat_rad * !radeg + dR ) * sin ( gei_lon_rad )

		GEI_R_unit_x	= ( GEI_x_R - GEI_x ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
		GEI_R_unit_y	= ( GEI_y_R - GEI_y ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )

		;	rotate the r vector by 90 degrees to get the thetat vector

		GEI_th_unit_x	= -GEI_R_unit_y
		GEI_th_unit_y	= GEI_R_unit_x

	;	now get the same info for the aacgm r and theta unit vectors,
   	;	but there are more steps here since we need to convert the aacgm
	;	coords to gei	

		aacgm_load_coef, year<2000 ; once we have newer coeffs update this

		;	small step in aacgm r direction & convert to geo coords

 		aacgm_conv_coord, (90.0 - aacgm_coLat_rad * !radeg)+dR, aacgm_lon_rad * !radeg, $
			 gei_R_km-6357.0, GEO_lat_deg_R, GEO_lon_deg_R, err, /to_geo

	 	;	convert geo to gei
	 	
 		geoPack_sphCar, gei_R_km, ( 90.0 - GEO_lat_deg_R ) * !dtor, GEO_lon_deg_R * !dtor, $
		   GEO_x_R, GEO_y_R, GEO_z_R, /to_rect        
		geoPack_conv_coord, GEO_x_R, GEO_y_R, GEO_z_R, $  ; units of km
				GEI_x_R, GEI_y_R, GEI_z_R, $
				/from_geo, /to_gei, $
				epoch = fltArr(n_elements(GEO_x_R)) + epoch
		geoPack_sphCar, GEI_x_R, GEI_y_R, GEI_z_R, $
				GEI_R_km_R, GEI_colat_rad_R, GEI_lon_rad_R,$
			   	/to_spher       

		;	get the x/y coords of the small step in aacgm r

		GEI_x_R	= ( GEI_colat_rad_R * !radeg ) * cos ( GEI_lon_rad_R )
		GEI_y_R	= ( GEI_colat_rad_R * !radeg ) * sin ( GEI_lon_rad_R )

		;	get aacgm r and thetat unit vector x and y components

		AACGM_R_unit_x	= ( GEI_x_R - GEI_x ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
		AACGM_R_unit_y	= ( GEI_y_R - GEI_y ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
		AACGM_th_unit_x	= -AACGM_R_unit_y
		AACGM_th_unit_y	= AACGM_R_unit_x

	;	get the x/y components of db	

	dbX	= gei_dbTh * GEI_R_unit_x + gei_dbPh * GEI_th_unit_x
	dbY	= gei_dbTh * GEI_R_unit_y + gei_dbPh * GEI_th_unit_y

	;	get the aacgm r/theta components which are the th/ph components
	;	that we want

	aacgm_dbTh	= dbX * AACGM_R_unit_x + dbY * AACGM_R_unit_y
	aacgm_dbPh	= dbX * AACGM_th_unit_x + dbY * AACGM_th_unit_y

	;	*** now i do not think this is working properly for a couple of 
	;	reasons. first since the th/ph comonents appear to be switched and
	;	the th component appears negative. i am at the moment compensating for
	;	this by switching the order and sign of the plotting of these results
	;	in clw_amp_v2. also the mapping is almost, but not quite right i 
	;	think.

	stop
;	;	get GEI xyz coords
;
;	geoPack_sphCar, gei_R_km, gei_coLat_rad, gei_lon_rad, $
;		   gei_x, gei_y, gei_z, /to_rect        
;
;    geoPack_sphCar, aacgm_R_km, aacgm_coLat_rad, aacgm_lon_rad, $
;		   aacgm_x, aacgm_y, aacgm_z, /to_rect        
;
;	;	get GEI xyz vector components
;
;	geoPack_bSpCar, gei_coLat_rad, gei_lon_rad, $
;			gei_dbTh*0, gei_dbTh, gei_dbPh, $
;			gei_dbx, gei_dby, gei_dbz
;
;	;	create rotation matrix from GEI to AACGM 
;
;	dx	= 1.
;	dy	= 1.
;	dz	= 1.
;
;	geoPack_conv_coord, aacgm_x+dx, aacgm_y, aacgm_z, $  ; units of km
;			xGEO_x, yGEO_x, zGEO_x, $
;			/from_gei, /to_geo, $
;			epoch = fltArr(n_elements(gei_x)) + epoch
;	geoPack_conv_coord, aacgm_x, aacgm_y+dy, aacgm_z, $  ; units of km
;			xGEO_y, yGEO_y, zGEO_y, $
;			/from_gei, /to_geo, $
;			epoch = fltArr(n_elements(gei_x)) + epoch
;	geoPack_conv_coord, aacgm_x, aacgm_y, aacgm_z+dz, $  ; units of km
;			xGEO_z, yGEO_z, zGEO_z, $
;			/from_gei, /to_geo, $
;			epoch = fltArr(n_elements(gei_x)) + epoch
;	
;	geoPack_sphCar, xGEO_x, yGEO_x, zGEO_x, $
;		   gei_R_km_x, gei_coLat_rad_x, gei_lon_rad_x, /to_spher       
;	geoPack_sphCar, xGEO_y, yGEO_y, zGEO_y, $
;		   gei_R_km_y, gei_coLat_rad_y, gei_lon_rad_y, /to_spher       
;	geoPack_sphCar, xGEO_z, yGEO_z, zGEO_z, $
;		   gei_R_km_z, gei_coLat_rad_z, gei_lon_rad_z, /to_spher       
;
;	aacgm_load_coef, year<2000 ; once we have newer coeffs update this
;
; 	aacgm_conv_coord, 90.0 - gei_coLat_rad_x * !radeg, gei_lon_rad_x * !radeg, $
;		 gei_R_km_x-6357.0, aacgm_lat_deg_x, aacgm_lon_deg_x, err, /to_aacgm
; 	aacgm_conv_coord, 90.0 - gei_coLat_rad_y * !radeg, gei_lon_rad_y * !radeg, $
;		 gei_R_km_y-6357.0, aacgm_lat_deg_y, aacgm_lon_deg_y, err, /to_aacgm
;  	aacgm_conv_coord, 90.0 - gei_coLat_rad_z * !radeg, gei_lon_rad_z * !radeg, $
;		 gei_R_km_z-6357.0, aacgm_lat_deg_z, aacgm_lon_deg_z, err, /to_aacgm
; 
; 	aacgm_coLat_rad_x	= (90.0-aacgm_lat_deg_x) * !dtor
; 	aacgm_coLat_rad_y	= (90.0-aacgm_lat_deg_y) * !dtor
; 	aacgm_coLat_rad_z	= (90.0-aacgm_lat_deg_z) * !dtor
;
;	aacgm_lon_rad_x	= aacgm_lon_deg_x * !dtor
;	aacgm_lon_rad_y	= aacgm_lon_deg_y * !dtor
;	aacgm_lon_rad_z	= aacgm_lon_deg_z * !dtor
;
;	geoPack_sphCar, gei_R_km_x, aacgm_coLat_rad_x, aacgm_lon_rad_x, $
;		   aacgm_x_x, aacgm_y_x, aacgm_z_x, /to_rect        
;	geoPack_sphCar, gei_R_km_y, aacgm_coLat_rad_y, aacgm_lon_rad_y, $
;		   aacgm_x_y, aacgm_y_y, aacgm_z_y, /to_rect        
;	geoPack_sphCar, gei_R_km_z, aacgm_coLat_rad_z, aacgm_lon_rad_z, $
;		   aacgm_x_z, aacgm_y_z, aacgm_z_z, /to_rect        
;
;	rotMat	= fltArr ( 3, 3, n_elements ( aacgm_R_km ) )
;	rotMat[0,0,*]	= aacgm_x-aacgm_x_x
;	rotMat[1,0,*]	= aacgm_x-aacgm_x_y
;	rotMat[2,0,*]	= aacgm_x-aacgm_x_z
;
;	rotMat[0,1,*]	= aacgm_y-aacgm_y_x
;	rotMat[1,1,*]	= aacgm_y-aacgm_y_y
;	rotMat[2,1,*]	= aacgm_y-aacgm_y_z
;
;	rotMat[0,2,*]	= aacgm_z-aacgm_z_x
;	rotMat[1,2,*]	= aacgm_z-aacgm_z_y
;	rotMat[2,2,*]	= aacgm_z-aacgm_z_z
;
;	rotMat	= -rotMat
;stop
;	;	apply rotation
;
;	dbAACGM	= rotMat ## [ $
;			[ gei_dbx[*] ],$
;			[ gei_dby[*] ],$
;			[ gei_dbz[*] ] ]
;
;	geoPack_bCarSp, aacgm_x, aacgm_y, aacgm_z, $
;			dbAACGM[*,0], dbAACGM[*,1], dbAACGM[*,2], $
;			aacgm_dbR, aacgm_dbTh, aacgm_dbPh
;	;geoPack_bCarSp, gei_x, gei_y, gei_z, $
;	;		gei_dbx, gei_dby, gei_dbz, $
;	;		aacgm_dbR, aacgm_dbTh, aacgm_dbPh
;
;
end
