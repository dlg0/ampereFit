pro read_ampere_dlg, fileName, sHr, eHr, dataOut, year, t_arr, $
		yrSec = yrSec

	header	= 2
	ndat	= file_lines ( fileName ) - header

	data_struct	= { 	ipln:0,$
						sv:0,$
						sday:0d0,$
						rad:0d0,$
						lat:0d0,$
						lon:0d0,$
						dbx:0d0,$
						dby:0d0,$
						dbz:0d0,$
						px:0d0,$
						py:0d0,$
						pz:0d0,$
						vx:0d0,$
						vy:0d0,$
						vz:0d0}

	data	= replicate ( data_struct, ndat )
	openr, rUnit, fileName, /get_lun

	dateStr=''
	readf, rUnit, dateStr
	reads, strmid(dateStr,10,8), $
			year, month, day, $
			format='(i4,i2,i2)'

	skip_lun, rUnit, 1, /lines
	readf, rUnit, data
	free_lun, rUnit
;
; ***********************
    month=3
    day=20
; ************************ puts noon at top

;	Select out a time slice

	iiTime	= where ( data.sday gt min(data.sday)+60.0*(sHr-8)*60.0 $
					and data.sday lt min(data.sday)+60.0*(eHr-8)*60.0 $
					and data.pz gt 0 )

;	Convert to GEI coords

; Not sure why this is here - CW
;	geopack_sphcar, data.px, data.py, data.pz,$
;			posR, posTh, posPh, /to_sphere, /degree

	cdf_epoch, epoch0, year, month, day, /compute_epoch
	epoch	= data.sday * 1d3 + epoch0
	t_arr	= epoch[iiTime]
	avgHour	= fix ( (sHr + eHr) / 2.0 )
	avgMin	= ( (sHr + eHr) / 2.0 mod 1 ) * 60

	yrSec	= cnvTime ( year, month, day, $
			avgHour, avgMin, 0.0 )

	geoPack_reCalc, year, month, day, /date

	geoPack_conv_coord, (data.px)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $  ; units of km
			xGEI, yGEI, zGEI, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]

	geoPack_sphCar, xGEI, yGEI, zGEI, R, theta, phi, /to_sphere         ; km -> km, radians
	geoPack_sphCar, (data.px)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $
		   geog_R_km, geog_coLat_rad, geog_lon_rad, /to_sphere         

;	Create rotation matrix from xyz to xyzGEI

	dx	= 1.
	dy	= 1.
	dz	= 1.

	geoPack_conv_coord, (data.px+dx)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $
			xGEI_x, yGEI_x, zGEI_x, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]
	geoPack_conv_coord, (data.px)[iiTime], (data.py+dy)[iiTime], (data.pz)[iiTime], $
			xGEI_y, yGEI_y, zGEI_y, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]
	geoPack_conv_coord, (data.px)[iiTime], (data.py)[iiTime], (data.pz+dz)[iiTime], $
			xGEI_z, yGEI_z, zGEI_z, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]

	rotMat	= fltArr ( 3, 3, n_elements ( iiTime ) )
	rotMat[0,0,*]	= xGEI-xGEI_x
	rotMat[1,0,*]	= xGEI-xGEI_y
	rotMat[2,0,*]	= xGEI-xGEI_z

	rotMat[0,1,*]	= yGEI-yGEI_x
	rotMat[1,1,*]	= yGEI-yGEI_y
	rotMat[2,1,*]	= yGEI-yGEI_z

	rotMat[0,2,*]	= zGEI-zGEI_x
	rotMat[1,2,*]	= zGEI-zGEI_y
	rotMat[2,2,*]	= zGEI-zGEI_z

;	Rotate db vectors to GEI
	dbGEI	= rotMat ## [ $
			[ (data.dbx)[iiTime] ],$
			[ (data.dby)[iiTime] ],$
			[ (data.dbz)[iiTime] ] ]

;	Rotate from xyzGEI to sphericalGEI
	geoPack_bCarSp, xGEI, yGEI, zGEI, $
			dbGEI[*,0], dbGEI[*,1], dbGEI[*,2], $
			bR, bTheta, bPhi

;	Colin, i don't think what you have below is correct. The geoPack_conv_coord
;	routine takes coordinates, not vector components. If you look at the 2,2 
;	component of the rotation matrix above you will see it is -1 (or 1 depending 
;	on if i got it right or not) meaning that the GEO and GEI have the same z 
;	axis as you told me. However, the x and y directions in GEO and GEI will be
;	different as time goes by of course, the rotation matrix takes that in to account
;	but i think if you took the xyz GEO db components and converted them into
;	spherical and then did a rotation in longitude (about the z axis) for the correct
;	time at each point you should get the same answer as what the above give. However, 
;	i'm not sure how to do that exactly, so apart from the +/- sign on the 1 in the 
;	rotation matrix i THINK the above is on the right track. I did try a straight 
;	conversion from the GEO xyz db vectors to spherical without accounting for the time
;	but is does not look right.
; -------------------------------------------------
; Take dB from GEO->GEI	: CW
;	geoPack_conv_coord, (data.dbx)[iiTime], (data.dby)[iiTime], (data.dbz)[iiTime], $
;			dbxGEI, dbyGEI, dbzGEI, $
;			/from_geo, /to_gei, $
;			epoch = epoch[iiTime]
;
;;	Rotate from xyzGEI to sphericalGEI
;	geoPack_bCarSp, xGEI, yGEI, zGEI, $
;			dbxGEI, dbyGEI, dbzGEI, $
;			bR, bTheta, bPhi
; --------------------------------------------------
	dataOut	= { R : R, $
				theta : theta, $
				phi : phi, $
				dbR : bR, $
			   	dbTheta : bTheta, $
				dbPhi : bPhi, $
				geog_R_km : geog_R_km, $
			   	geog_coLat_rad : geog_coLat_rad, $
				geog_lon_rad : geog_lon_rad }
		
end
