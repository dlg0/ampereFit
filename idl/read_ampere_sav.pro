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
    south = south, $
    fillPole = fillPole, $
    fillPoleLine = fillPoleLine
    

    if not keyword_set ( savFileName ) then $
        savFileName = '~/code/ampereFit/data/20091023Amp_invert.sav'
    if not keyword_set ( capSize ) then capSize = 40 * !dtor

	dateStr=''
	reads, strmid(file_baseName ( savFileName ),0,8), $
			year, month, day, $
			format='(i4,i2,i2)'

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
                    isat: 0, $
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
    data.isat   = pseudoSVNum_total[iiTime]
    data.utc    = x_axis_frac_hour[iiTime]

    data.dbx    = (B_ECI[0,*])[iiTime]
    data.dby    = (B_ECI[1,*])[iiTime]
    data.dbz    = (B_ECI[2,*])[iiTime]

;   add some extra data near the pole ;-)

    if keyword_set ( fillPoleLine ) then begin

        ;   identify which two tracks surround the gap

        iiTenLat = where ( gei_coLat_rad_TMP[iiTime]*!radeg gt 6 $
                    and gei_coLat_rad_TMP[iiTime]*!radeg lt 8, iiTenCnt ) 

        tenLons = (gei_lon_rad_TMP[iiTime])[iiTenLat] * !radeg        
        iiSortTenLons   = sort ( tenLons )
        tenLons = tenLons[iiSortTenLons]

        derivLons   = fltArr ( n_elements ( tenLons ) )
        for i=0,n_elements ( tenLons ) - 2 do begin 
        
            derivLons[i]   = tenLons[i+1] - tenLons[i]
        
        endfor

        derivLons[n_elements(tenLons)-1]    = tenLons[0]+360 - tenLons[n_elements(tenLons)-1]

        iiMaxGap    = where ( derivLons eq max ( derivLons ) )

        track1  = (((plane_number_total[iiTime])[tenLons])[iiSortTenLons])[iiMaxGap]
        track2  = (((plane_number_total[iiTime])[tenLons])[iiSortTenLons])[iiMaxGap+1]

        if track1 eq track2 then begin

            print, 'Awww crap!, the detection of which tracks straddle the intersection'
            print, 'point hole has failed.'

        endif

        iiTrack1    = where ( plane_number_total[iiTime] eq track1 )
	    map_set, 90, 0, 0, $
	    	/ortho, $
   	    	/iso, $
        	limit = [ 90.0 - ( 50 + 2 ), 0, 90, 360 ],$
   	    	/noborder,$
   	    	/advance, title = 'GEI'
	    map_grid, label = 1, latDel = 10.0

	    ;plots, (gei_lon_rad_TMP[iiTime])[iiTrack1]*!radeg, $
        ;        90.0-(gei_coLat_rad_TMP[iiTime])[iiTrack1]*!radeg, psym = 4


        stop
    endif

    if keyword_set ( fillPole ) then begin

        nPts    = 5 
        spreadFac   = 200 
        newX    = rebin ( fIndGen ( nPts ) - nPts / 2, nPts, nPts ) * spreadFac
        newY    = transpose ( newX ) 

        iiKeep    = where ( newX ne 0 and newY ne 0, nExtra )
        if nExtra gt 0 then begin
            newX    = newX[iiKeep]
            newY    = newY[iiKeep]
        endif

        if keyword_set ( south ) then minVal = min ( data.pz ) $
            else minVal = max ( data.pz )
        newZ    = newX * 0 + minVal 

        minIdx  = where ( data.pz eq minVal )
        new_dbX = (data.dbx)[minIdx]
        new_dbY = (data.dby)[minIdx]
        new_dbZ = (data.dbz)[minIdx]

        data_    = replicate ( data_struct, iiTimeCnt + nExtra )

        data_.px = [ data.px, newX[*] ]
        data_.py = [ data.py, newY[*] ]
        data_.pz = [ data.pz, newZ[*] ]

        data_.dbx    = [ data.dbx, fltArr ( nExtra ) + new_dbX[0] ]
        data_.dby    = [ data.dby, fltArr ( nExtra ) + new_dbY[0] ]
        data_.dbz    = [ data.dbz, fltArr ( nExtra ) + new_dbZ[0] ]

        data_.isat   = [ data.isat, fltArr ( nExtra ) + 99 ]
        data_.utc    = [ data.utc, fltArr ( nExtra ) + data.utc[minIdx] ]

        data = data_

        iiTimeCnt += nExtra 
        
    endif

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

    if keyword_set ( south ) then $
        bTheta_GEI  = -bTheta_GEI

    ;for i = 0, n_elements ( data.px ) - 1 do begin

    ;    th  = gei_colat_rad[i]
    ;    ph  = gei_lon_rad[i]
    ;    r   = 1 
    ;    rotMat  = [ [ sin(th)*cos(ph), sin(th)*sin(ph), cos(th) ], $
    ;                [ r*cos(th)*cos(ph), r*cos(th)*sin(ph), -r*sin(th) ], $
    ;                [ -r*sin(ph), r*cos(ph), 0 ] ]

    ;    sphVec  = rotMat ## [   [(data.dbx)[i]], $
    ;                            [(data.dby)[i]], $
    ;                            [(data.dbz)[i]] ]

    ;    bR_GEI[i]   = sphVec[0]
    ;    bTheta_GEI[i]   = sphVec[1]
    ;    bPhi_GEI[i] = sphVec[2]
    ;endfor


	;plt_dat, gei_coLat_rad * !radeg, gei_lon_rad * !radeg, $
	;	n_elements ( gei_coLat_rad ), -bTheta_GEI, bPhi_GEI, [1,1], [1,2], $
	;	title = 'GEI - Raw Data', $
    ;    satNu = data.iSat, $
    ;    capSize = 50 

    ;plot_vec, gei_coLat_rad, gei_lon_rad, $
    ;    -bTheta_GEI, bPhi_GEI
stop	
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
				mlt : mlt, $
                iPln : data.iPln, $
                iSat : data.iSat, $
                utc : data.utc } 
stop    
end
