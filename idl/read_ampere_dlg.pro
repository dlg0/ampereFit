pro read_ampere_dlg, $
	   	fileName, $
		sHr, eHr, $
		dataOut, $
		t_arr, $
		capSize, $
		yrSec = yrSec, $
		year = year, $
		month = month, $
		day = day, $
		avgYrSec = yrSecAvg, $
		avgEpoch = avgEpoch

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

;;
;; ***********************
;    month=3
;    day=20
;; ************************ puts noon at top - what?

;	Select out a time slice

	iiTime	= where ( data.sday gt min(data.sday)+60.0*(sHr-8)*60.0 $
					and data.sday lt min(data.sday)+60.0*(eHr-8)*60.0 $
					and data.pz gt 0, iiCnt )

	data	= data[iiTime]

;	Convert to GEI coords

; Not sure why this is here - CW
;	geopack_sphcar, data.px, data.py, data.pz,$
;			posR, posTh, posPh, /to_sphere, /degree

	cdf_epoch, epoch0, year, month, day, /compute_epoch
	epoch	= data.sday * 1d3 + epoch0
	avgEpoch	= mean ( epoch )
	t_arr	= epoch
	avgHour	= fix ( (sHr + eHr) / 2.0 )
	avgMin	= ( (sHr + eHr) / 2.0 mod 1 ) * 60

	yrSecAvg	= cnvTime ( year, month, day, $
			avgHour, avgMin, 0.0 )

	julTime	= julDay ( month, day, year, 0.0, 0.0, data.sday )
	calDat, julTime, month_, day_, year_, hour_, minute_, second_

	yrSec	= fltArr ( n_elements(year_) )
	for i=0,n_elements(year_)-1 do begin
		yrSec[i]	= cnvTime ( year_[i], month_[i], day_[i], $
			hour_[i], minute_[i], second_[i] )
	endfor

	geoPack_reCalc, year, month, day, /date

	geoPack_conv_coord, (data.px), (data.py), (data.pz), $  ; units of km
			xGEI, yGEI, zGEI, $
			/from_geo, /to_gei, $
			epoch = epoch

	geoPack_sphCar, xGEI, yGEI, zGEI, gei_R_km, gei_coLat_rad, gei_lon_rad, /to_sphere         ; km -> km, radians
	geoPack_sphCar, (data.px), (data.py), (data.pz), $
		   geog_R_km, geog_coLat_rad, geog_lon_rad, /to_sphere         


;	Create rotation matrix from xyzGEO to xyzGEI

	dx	= 1.
	dy	= 1.
	dz	= 1.

	geoPack_conv_coord, (data.px+dx), (data.py), (data.pz), $
			xGEI_x, yGEI_x, zGEI_x, $
			/from_geo, /to_gei, $
			epoch = epoch
	geoPack_conv_coord, (data.px), (data.py+dy), (data.pz), $
			xGEI_y, yGEI_y, zGEI_y, $
			/from_geo, /to_gei, $
			epoch = epoch
	geoPack_conv_coord, (data.px), (data.py), (data.pz+dz), $
			xGEI_z, yGEI_z, zGEI_z, $
			/from_geo, /to_gei, $
			epoch = epoch

	rotMat	= fltArr ( 3, 3, iiCnt )
	rotMat[0,0,*]	= xGEI-xGEI_x
	rotMat[1,0,*]	= xGEI-xGEI_y
	rotMat[2,0,*]	= xGEI-xGEI_z

	rotMat[0,1,*]	= yGEI-yGEI_x
	rotMat[1,1,*]	= yGEI-yGEI_y
	rotMat[2,1,*]	= yGEI-yGEI_z

	rotMat[0,2,*]	= zGEI-zGEI_x
	rotMat[1,2,*]	= zGEI-zGEI_y
	rotMat[2,2,*]	= zGEI-zGEI_z

	rotMat	= -rotMat

;	Rotate db vectors to GEI
	dbGEI	= rotMat ## [ $
			[ (data.dbx) ],$
			[ (data.dby) ],$
			[ (data.dbz) ] ]

;	Convert from xyzGEI to sphericalGEI
	geoPack_bCarSp, xGEI, yGEI, zGEI, $
			dbGEI[*,0], dbGEI[*,1], dbGEI[*,2], $
			bR, bTheta, bPhi

;	so this has now been tested and it is all consistent within a minus sign
;	with the method Colin has below

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


	aacgm_load_coef, year<2000 ; once we have newer coeffs update this
 	aacgm_conv_coord, 90.0 - geog_coLat_rad * !radeg, geog_lon_rad * !radeg, $
		 geog_R_km-6371.0, aacgm_lat_deg, aacgm_lon_deg, err, /to_aacgm
 	aacgm_coLat_deg	= 90.0-aacgm_lat_deg
	mlt	= aacgm_mlt ( fltArr(iiCnt) + year, fltArr(iiCnt) + yrSec, aacgm_lon_deg )

	iiNeg	= where ( aacgm_lon_deg lt 0, iiNegCnt )
	aacgm_lon_deg[iiNeg]	= aacgm_lon_deg[iiNeg] + 360

	;	select out those observations within the desired capSize

	iiCap	= where ( gei_coLat_rad lt capSize * !dtor, iiCapCnt )

	dataOut	= { gei_R_km : gei_R_km[iiCap], $
				gei_coLat_rad : gei_coLat_rad[iiCap], $
				gei_lon_rad : gei_lon_rad[iiCap], $
				dbR : bR[iiCap], $
			   	dbTheta : bTheta[iiCap], $
				dbPhi : bPhi[iiCap], $
				geog_R_km : geog_R_km[iiCap], $
			   	geog_coLat_rad : geog_coLat_rad[iiCap], $
				geog_lon_rad : geog_lon_rad[iiCap], $
			    aacgm_coLat_rad : aacgm_coLat_deg[iiCap] * !dtor, $
				aacgm_lon_rad : aacgm_lon_deg[iiCap] * !dtor, $
				mlt : mlt[iiCap] }
end
