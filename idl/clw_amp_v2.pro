pro clw_amp_v2, $
		numerical = numerical, $
        plot_bFns = plot_bFns

	if strCmp ( !version.os, 'linux' ) or strCmp ( !version.os, 'darwin' ) then begin
		plotDev = 'X'
		path	= '~/code/ampereFit/idl/'
		pnmPath	= path + 'pnmSavs/pnmSav'
	endif else begin
		plotDev	= 'win'
		path	= 'd:\cwac\hi_res\'
		pnmPath	= path + 'pnmsavs\pnmSav'
	endelse

	capSize	= 55.0
	plotCapSize	= 50.0
	fileName = path + '20050515_a_RevB.dat'
	
	sHr = 23.4 
	eHr = 23.6
	
	;read_ampere_dlg, fileName, sHr, eHr, data, t_arr, capSize, $
	;   	 yrSec = yrSec, $
	;   	 year = year, $
	;   	 month = month, $
	;   	 day = day, $
	;	 avgYrSec = avgYrSec, $
	;   	 avgEpoch = avgEpoch

    savFileName = '~/code/ampereFit/data/20091023Amp_invert.sav'
	
    read_ampere_sav, $
        savFileName = savFileName, $
        capSize = capSize, $
        dataOut = data, $
        sHr = sHr, $
        eHr = eHr, $
        year = year, $
        month = month, $
        day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        south = 1, $
        fillPole = 0

	nLatGrid	= 200
	nLonGrid	= 10 

	aacgm_grid, geiGrid_coLat_rad = geiGrid_coLat_rad, $
	   	 geiGrid_lon_rad = geiGrid_lon_rad , $
	   	 geiGrid_R_km = geiGrid_R_km, $
	   	 aacgm_coLat_deg = aacgm_coLat_deg, $
	   	 aacgm_lon_deg = aacgm_lon_deg, $
	   	 aacgmGrid_R_km = aacgmGrid_R_km, $
	   	 nLat = nLatGrid, nLon = nLonGrid, $
	   	 yrSec = avgYrSec, year = year, $
	   	 mltShift = mltShift, $
	   	 aacgmGrid_coLat_deg = aacgmGrid_coLat_deg, $
	   	 aacgmGrid_lon_deg = aacgmGrid_lon_deg, $
	   	 epoch = avgEpoch
    
    mltShift    = mltShift[0]
	
	kMax    = 50 
	mMax    = 5

	schaBasisFunctions, kMax, mMax, capSize, data.gei_coLat_rad, data.gei_lon_rad, $
	     YkmBFns=YkmBFns, $
	   	 dYkmDthBFns=dYkmDthBFns, $
	   	 dYkmDphBFns=dYkmDphBFns, $
	     OUTNKVALUES=outNkValues, $
	   	 OUTKVALUES=outKValues, $
	   	 OUTMVALUES=outMValues, $
	   	 pnmPath = pnmPath;, $
         ;/evenSet

 	ampere_setupSHFns, 1.0, data.gei_coLat_rad, data.gei_lon_rad, $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns2, $
			bPhBFnArr = dYkmDPhBFns2, $
			YBFnArr = YkmBFns2, $
			mArr = outMValues2, $
			nArr = outNkValues2, $
			kArr = outKValues2;, $
            ;/bc2

    if keyword_set ( plot_bFns ) then begin
	    !p.multi = [0,3,4]
	    !p.charSize = 2.0
	    device, decomposed = 0
	    window, 0, xSize = 1200, ySize = 800
	    !p.background = 255

        m1  = 0
        k1  = 6
	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5, Y)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdTh', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdPh', $
	    		psym=4,color=0
	    
	    plot, data.gei_coLat_rad*!radeg, ykmbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
	    		title = 'numerical (m=2,k=5), Y', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1 ),*], $
	    		title = 'numerical (m=2,k=5), dYdTh', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1 ),*], $
	    		title = 'numerical (m=2,k=5), dYdPh', $
	    		psym=4,color=0

        m2  = 0 
        k2  = 8 
	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0

	    plot, data.gei_coLat_rad*!radeg, ykmbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0
	    
	    !p.multi = 0.0
	    !p.charSize = 1.0

    endif
	
	if keyword_set ( numerical ) then begin

        ;iiSort  = sort ( outNkValues2 )
		YkmBFns	= temporary ( YkmBFns2);[iiSort,*] )
		dYkmDthBFns	= temporary ( dYkmDthBFns2);[iiSort,*] )
		dYkmDphBFns	= temporary ( dYkmDphBFns2);[iiSort,*] )
		outNkValues	= temporary ( outNkValues2);[iiSort] )

	endif


;	Fit |dB| to dB.grad Ykm in the GEI system
;	-----------------------------------------	
	
	dBMag   = sqrt(data.dBTheta^2+data.dBPhi^2)

    bFuncs  = temporary ( $
                [ [ dYkmDThBfns, -dYkmDPhBfns ], $
                  [ dYkmDPhBfns,  dYkmDthBfns ] ] )

    print, 'bFuncs is', n_elements ( bFuncs[*] ) * 16 / ( 1024.0^2 ), 'MB' 

    coeffs_    = la_least_squares ( bFuncs, [data.dbTheta, data.dbPhi], $
                    status = stat, $
                    method = 3, $
                    /double, $
                    rank = rank )
    if stat ne 0 then begin
        print, 'ERROR: la_least_squares threw an error code: ', stat
    endif

    fit = bFuncs ## coeffs_
	
	fit_dBTheta     = fit[0:n_elements(dBMag)-1] 
	fit_dBPhi       = fit[n_elements(dBMag):*] 

	Re    = 6.371d6
	R     = Re + 780.0d3
	u0    = 4.0*!dpi*1.0d-7
	jPar  = (YkmBFns ## $
	        (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}

;	generate basis fns at regular grid
;	----------------------------------

	schaBasisFunctions, kMax, mMax, capSize, geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
         YkmBFns=YkmBFns_grid, $
		 dYkmDthBFns=dYkmDthBFns_grid, $
		 dYkmDphBFns=dYkmDphBFns_grid, $
		 pnmPath = pnmPath;,$
         ;/oddSet

	ampere_setupSHFns, 1.0, geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns_grid2, $
			bPhBFnArr = dYkmDPhBFns_grid2, $
			YBFnArr = YkmBFns_grid2, $
			mArr = outMValues_grid2, $
			nArr = outNkValues_grid2, $
			kArr = outKValues_grid2;, $
            ;/bc2

	if keyword_set ( numerical ) then begin

		YkmBFns_grid	= temporary ( YkmBFns_grid2 )
		dYkmDthBFns_grid	= temporary ( dYkmDthBFns_grid2 )
		dYkmDphBFns_grid	= temporary ( dYkmDphBFns_grid2 )

	endif

    bFuncs_grid  = temporary ( $
                    [ [ dYkmDThBfns_grid, -dYkmDPhBfns_grid ], $
                      [ dYkmDPhBfns_grid,  dYkmDthBfns_grid ] ] )


  	jParAACGM	= (YkmBFns_grid ## $
         (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
                         /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}
	jParAACGM	= reform ( jParAACGM, nLatGrid, nLonGrid )

    fit_grid    = bFuncs_grid ## coeffs_
   	dBTheta_GEI_grid     = fit_grid[0:nLatGrid*nLonGrid-1] 
 	dBPhi_GEI_grid       = fit_grid[nLatGrid*nLonGrid:*]
 
	rotate_gei_to_aacgm, geiGrid_R_km[*], geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
		aacgmGrid_R_km[*], aacgmGrid_coLat_deg[*]*!dtor, aacgmGrid_lon_deg[*]*!dtor, $
		dBTheta_GEI_grid, dBPhi_GEI_grid, $
		aacgm_dbTh = aacgm_dbTh, $
		aacgm_dbPh = aacgm_dbPh, $
		year = year, $
		epoch = avgEpoch, $
		/dot

;	plot stuff
;	----------

	set_plot, plotDev
	device,decomposed = 0
	!p.background = 255
	winNum	= 1

;	plot the db vectors
;	-------------------

	!p.multi = [0,2,2]

	window, winnum, xSize = 800, ySize = 800
	winNum++

	plt_dat, data.gei_coLat_rad * !radeg, data.gei_lon_rad * !radeg, $
		n_elements ( data.gei_coLat_rad ), -data.dbTheta, data.dbPhi, [1,1], [1,2], $
		title = 'GEI - Raw Data', $
        satNu = data.iSat
	plt_dat, data.gei_coLat_rad * !radeg, data.gei_lon_rad * !radeg, $
	         n_elements ( data.gei_coLat_rad ), -fit_dBTheta, fit_dBPhi, [1,1], [1,2], $
			 title = 'GEI - Fitted Data @ raw locs'
	plt_dat, geiGrid_coLat_rad[*] * !radeg, geiGrid_lon_rad[*] * !radeg, $
	         n_elements ( geiGrid_coLat_rad[*] ), -dBTheta_GEI_grid, dBPhi_GEI_grid, [1,1], [1,2], $
			 title = 'GEI - on regular aacgm grid'
	plt_dat, aacgmGrid_coLat_deg[*], aacgmGrid_lon_deg[*], $
	         n_elements ( geiGrid_coLat_rad[*] ), -aacgm_dbTh, aacgm_dbPh, [1,1], [1,2], $
			 title = 'AACGM-LON'

;	plot the FAC maps
;	-----------------

	loadct, 13, file = path + 'davect.tbl'
	window, winnum, xSize = 600, ySize = 600
	winnum++
	!p.backGround = 0
	!p.multi = [0,2,2]
	!p.font = 0

	map_set, 90, 0, 0, /ortho, /iso, $
	    limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	    /noborder, /advance, title = 'jPar GEI'
	jLevels = ( fIndGen ( 21 ) - 10 ) * 0.2;2.0
	colors  = bytScl ( jLevels, top = 253 ) + 1
	oldbck=!p.background
	!p.background=0
	jpar=reform(jpar)
    triangulate, data.gei_lon_rad*!radeg, 90-data.gei_colat_rad*!radeg, $
        tri, sphere = sphere
	contour, jPar, data.gei_lon_rad*!radeg, 90.0-data.gei_coLat_rad*!radeg, $
	   	 	c_labels=fltarr(n_elements(jLevels))+1,$
	       	/irreg, /over, levels = jLevels, $
	          	c_colors = colors, /fill
	map_grid, label = 1, latDel = 10.0 


   	map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar AACGM-MLT'
	jParTmp	= jParAACGM[*]
	lonTmp = ((aacgmGrid_lon_deg[*]/15.0+mltShift) mod 24 )*15.0 $
			+ randomu ( sysTime(/sec ), n_elements(jParTmp), /uni ) * 1e-5
	latTmp = 90.0-aacgmGrid_coLat_deg[*]
   	contour, jParTmp, $
			lonTmp, $
			latTmp, $
		 	c_labels=fltarr(n_elements(jLevels))+1,$
        	/over, levels = jLevels, $
           	c_colors = colors, /fill, /irreg
	map_grid, label = 1, $ 
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360], $
			latLab = 45, $
			latDel = 10.0


	map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar AACGM-LON'
    contour, transpose ( jParAACGM ), $
            aacgmGrid_lon_deg[0,*], $
            90.0 - aacgmGrid_coLat_deg[*,0], $
            /over, levels = jLevels, $
            c_colors = colors, /fill
	map_grid, label = 1, latDel = 10.0


	map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar AACGM-LON [GEI]'
	jParTmp	= jParAACGM[*]
	lonTmp = geiGrid_lon_rad[*]*!radeg
	latTmp = 90.0-geiGrid_coLat_rad[*]*!radeg
	contour, jParTmp, $
			lonTmp, $
			latTmp, $
		 	c_labels=fltarr(n_elements(jLevels))+1,$
        	/over, levels = jLevels, $
           	c_colors = colors, /fill, /irreg
	map_grid, label = 1, latDel = 10.0


;	plot the data locations
;	-----------------------

	!p.multi = [0,2,2]
	window, winnum, xSize = 600, ySize = 600
	winnum++
  
 	map_set, 90, 0, 0, $
		/ortho, $
   		/iso, $
    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
   		/noborder,$
   		/advance, title = 'AACGM-LON'
	map_grid, label = 1, latDel = 10.0

	plots, data.aacgm_lon_rad*!radeg, 90.0-data.aacgm_coLat_rad*!radeg, psym = 4

 	map_set, 90, 0, 0, $
		/ortho, $
   		/iso, $
    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
   		/noborder,$
   		/advance, title = 'AACGM-MLT'
	map_grid, label = 1, latDel = 10.0

	plots, data.mlt*15.0, 90.0-data.aacgm_coLat_rad*!radeg, psym = 4
	
   	map_set, 90, 0, 0, $
		/ortho, $
   		/iso, $
    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
   		/noborder,$
   		/advance, title = 'GEI'
	map_grid, label = 1, latDel = 10.0

	plots, data.gei_lon_rad*!radeg, 90.0-data.gei_coLat_rad*!radeg, psym = 4
	
	map_set, 90, 0, 0, $
		/ortho, $
   		/iso, $
    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
   		/noborder,$
   		/advance, title = 'GE0G'
	map_grid, label = 1, latDel = 10.0

	plots, data.geog_lon_rad*!radeg, 90.0-data.geog_coLat_rad*!radeg, psym = 4
	
	!p.multi = 0
 	!p.background = oldbck

 stop
end

