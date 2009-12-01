pro read_ampere_sav, $
    savFileName = savFileName, $
    capSize = capSize, $
    dataOut = dataOut, $
    sHr = sHr, eHr = eHr, $
    year = year, $
    month = month, $
    day = day, $
    avgYrSec = yrSecAvg, $
    avgEpoch = avgEpoch, $
    south = south
    

    if not keyword_set ( savFileName ) then $
        savFileName = '~/code/ampereFit/data/Amp_invert.sav'
    if not keyword_set ( capSize ) then capSize = 40 * !dtor

    restore, savFileName

    ;   variables are 
    ;   
    ;   x_axis_frac_hour
    ;   xaxis_frac_hour UTC in fractional hours
    ;
    ;   plane_number_total
    ;   plane_number_total, is integer number 0-5
    ;   for the 6 planes that the sats fly in.
    ;
    ;   pos_ECI_total
    ;   pos_eci is ECI position in meters and has 
    ;   dimensions [3,x]  where x is number of points
    ;   and 0-2 in the first dimension is X,Y,Z

    ;   B_ECI
    ;   B_eci is the same dimension as pos_eci but has delta B in nT

    data_struct = { ipln: 0, $
                    utc: 0d0, $
                    px: 0d0, $
                    py: 0d0, $
                    pz: 0d0, $
                    dbx: 0d0, $
                    dby: 0d0, $
                    dbz: 0d0 }

    ;   here we do the spherical coord calculation for GEI twice
    ;   just to avoid the 10000 limit on the epoch vector size :-(

	geoPack_sphCar, pos_ECI_total[0,*]*1d-3, $
                    pos_ECI_total[1,*]*1d-3, $
                    pos_ECI_total[2,*]*1d-3, $
            gei_R_km_TMP, gei_coLat_rad_TMP, gei_lon_rad_TMP, $
             /to_sphere         ; km -> km, radians

    if keyword_set ( south ) then $
        gei_coLat_rad_Tmp   = !pi - gei_coLat_rad_Tmp

    iiTime  = where ( x_axis_frac_hour ge sHr $
                        and x_axis_frac_hour le eHr $
                        and gei_coLat_rad_TMP lt capSize * !dtor, iiTimeCnt )

    data    = replicate ( data_struct, iiTimeCnt )

    data.px = (pos_ECI_total[0,*])[iiTime]*1d-3
    data.py = (pos_ECI_total[1,*])[iiTime]*1d-3
    data.pz = (pos_ECI_total[2,*])[iiTime]*1d-3

    data.ipln   = plane_number_total[iiTime]
    data.utc    = x_axis_frac_hour[iiTime]

    data.dbx    = (B_ECI[0,*])[iiTime]
    data.dby    = (B_ECI[1,*])[iiTime]
    data.dbz    = (B_ECI[2,*])[iiTime]

;   get spherical GEI coords

	geoPack_sphCar, data.px, data.py, data.pz, $
            gei_R_km, gei_coLat_rad, gei_lon_rad, $
             /to_sphere         ; km -> km, radians

    if keyword_set ( south ) then $
        gei_coLat_rad   = !pi - gei_coLat_rad

;   get spherical GEI vector

    geoPack_bCarSp, data.px, data.py, data.pz, $
			data.dbx, data.dby, data.dbz, $
			bR_GEI, bTheta_GEI, bPhi_GEI

;   get the epoch and times for GEI to GEOG and AACGM conversion

	cdf_epoch, epoch0, year, month, day, /compute_epoch
	epoch	= data.utc * 3600d0 * 1d3 + epoch0
	avgEpoch	= mean ( epoch )
	t_arr	= epoch
	avgHour	= fix ( (sHr + eHr) / 2.0 )
	avgMin	= ( (sHr + eHr) / 2.0 mod 1 ) * 60

	yrSecAvg	= cnvTime ( year, month, day, $
			avgHour, avgMin, 0.0 )

	julTime	= julDay ( month, day, year, 0.0, 0.0, data.utc * 3600d0 )
	calDat, julTime, month_, day_, year_, hour_, minute_, second_

	yrSec	= fltArr ( n_elements(year_) )
	for i=0,n_elements(year_)-1 do begin
		yrSec[i]	= cnvTime ( year_[i], month_[i], day_[i], $
			hour_[i], minute_[i], second_[i] )
	endfor

;   convert to GEOG from GEI

	geoPack_reCalc, year, month, day, /date

    for i = 0, n_elements ( data.px ) / 10000 do begin  

        print, i
        thisIndex   = indGen(10000)+i*10000
        iiKeep  = where ( thisIndex lt n_elements(data.px) )
        thisIndex   = thisIndex[iiKeep]
	    geoPack_conv_coord, (data.px)[thisIndex], $
                            (data.py)[thisIndex], $
                            (data.pz)[thisIndex], $  ; units of km
			xGEO_tmp, yGEO_tmp, zGEO_tmp, $
			/from_gei, /to_geo, $
			epoch = epoch[thisIndex]

        if size ( xGEO, /dim ) ne 0 then begin

            xGEO    = [ xGEO, xGEO_tmp ]
            yGEO    = [ yGEO, yGEO_tmp ]
            zGEO    = [ zGEO, zGEO_tmp ]

        endif else begin

            xGEO    = xGEO_tmp
            yGEO    = yGEO_tmp
            zGEO    = zGEO_tmp
 
        endelse

    endfor

	geoPack_sphCar, xGEO, yGEO, zGEO, $
            geog_R_km, geog_coLat_rad, geog_lon_rad, $
             /to_sphere         ; km -> km, radians

    if keyword_set ( south ) then $
        geog_coLat_rad = !pi - geog_coLat_rad

;   get AACGM coords

   	aacgm_load_coef, year<2000 ; once we have newer coeffs update this
 	aacgm_conv_coord, 90.0 - geog_coLat_rad * !radeg, geog_lon_rad * !radeg, $
		 geog_R_km-6357.0, aacgm_lat_deg, aacgm_lon_deg, err, /to_aacgm
 	aacgm_coLat_deg	= 90.0-aacgm_lat_deg
	mlt	= aacgm_mlt ( fltArr(iiTimeCnt) + year, fltArr(iiTimeCnt) + yrSec, aacgm_lon_deg )

	iiNeg	= where ( aacgm_lon_deg lt 0, iiNegCnt )
	if iiNegCnt gt 0 then $
        aacgm_lon_deg[iiNeg] = aacgm_lon_deg[iiNeg] + 360 

;   select out data from capSize and create output data structure

	;iiCap	= where ( gei_coLat_rad lt capSize * !dtor, iiCapCnt )

	dataOut	= { gei_R_km : gei_R_km, $
				gei_coLat_rad : gei_coLat_rad, $
				gei_lon_rad : gei_lon_rad, $
				dbR : bR_GEI, $
			   	dbTheta : bTheta_GEI, $
				dbPhi : bPhi_GEI, $
				geog_R_km : geog_R_km, $
			   	geog_coLat_rad : geog_coLat_rad, $
				geog_lon_rad : geog_lon_rad, $
			    aacgm_coLat_rad : aacgm_coLat_deg * !dtor, $
				aacgm_lon_rad : aacgm_lon_deg * !dtor, $
				mlt : mlt } 

    stop
end
