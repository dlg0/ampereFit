pro clw_amp_v2, $
		numerical = numerical

; @'d:\cwac\hi_res\davidg\jpar_ver2\schabasisfunctions.pro'

	if strCmp ( !version.os, 'linux' ) then begin
		plotDev = 'X'
		path	= '~/code/ampereFit/idl/'
		pnmPath	= path + 'pnmSavs/pnmSav'
	endif else begin
		plotDev	= 'win'
		path	= 'd:\cwac\hi_res\'
		pnmPath	= path + 'pnmsavs\pnmSav'
	endelse

	capSize	= 50.0
	plotCapSize	= 40.0
	;fileName = path + '20080105_a_RevA.dat'
	fileName = path + '20050515_a_RevB.dat'
	
	sHr = 11 
	eHr = 13 
	
	read_ampere_dlg, fileName, sHr, eHr, data, t_arr, capSize, $
	   	 yrSec = yrSec, $
	   	 year = year, $
	   	 month = month, $
	   	 day = day, $
		 avgYrSec = avgYrSec, $
	   	 avgEpoch = avgEpoch
	
	nLatGrid	= 20
	nLonGrid	= 24 

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
	
	kMax    = 30 
	mMax    = 5

	schaBasisFunctions, kMax, mMax, capSize, data.gei_coLat_rad, data.gei_lon_rad, $
	     YkmBFns=YkmBFns, $
	   	 dYkmDthBFns=dYkmDthBFns, $
	   	 dYkmDphBFns=dYkmDphBFns, $
	     OUTNKVALUES=outNkValues, $
	   	 OUTKVALUES=outKValues, $
	   	 OUTMVALUES=outMValues, $
	   	 pnmPath = pnmPath, /evenSet

 	ampere_setupSHFns, 1.0, data.gei_coLat_rad, data.gei_lon_rad, $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns2, $
			bPhBFnArr = dYkmDPhBFns2, $
			YBFnArr = YkmBFns2, $
			mArr = outMValues2, $
			nArr = outNkValues2, $
			kArr = outKValues2
stop
	!p.multi = [0,3,4]
	!p.charSize = 2.0
	device, decomposed = 0
	window, 0, xSize = 1200, ySize = 800
	!p.background = 255
	plot, data.gei_coLat_rad*!radeg, ykmbfns[(where(outmvalues eq 2))[5],*], $
			title = 'table lookup (m=2,k=5, Y)', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdthbfns[(where(outmvalues eq 2))[5],*], $
			title = 'table lookup (m=2,k=5), dYdTh', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdphbfns[(where(outmvalues eq 2))[5],*], $
			title = 'table lookup (m=2,k=5), dYdPh', $
			psym=4,color=0
	
	plot, data.gei_coLat_rad*!radeg, ykmbfns2[(where(outmvalues2 eq 2))[5],*], $
			title = 'numerical (m=2,k=5), Y', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[(where(outmvalues2 eq 2))[5],*], $
			title = 'numerical (m=2,k=5), dYdTh', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[(where(outmvalues2 eq 2))[5],*], $
			title = 'numerical (m=2,k=5), dYdPh', $
			psym=4,color=0

	plot, data.gei_coLat_rad*!radeg, ykmbfns[(where(outmvalues eq 4))[3],*], $
			title = 'table lookup (m=4,k=3)', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdthbfns[(where(outmvalues eq 4))[3],*], $
			title = 'table lookup (m=4,k=3)', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdphbfns[(where(outmvalues eq 4))[3],*], $
			title = 'table lookup (m=4,k=3)', $
			psym=4,color=0

	plot, data.gei_coLat_rad*!radeg, ykmbfns2[(where(outmvalues2 eq 4))[3],*], $
			title = 'numerical (m=4,k=3)', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[(where(outmvalues2 eq 4))[3],*], $
			title = 'numerical (m=4,k=3)', $
			psym=4,color=0
	plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[(where(outmvalues2 eq 4))[3],*], $
			title = 'numerical (m=4,k=3)', $
			psym=4,color=0
	
	!p.multi = 0.0
	!p.charSize = 1.0
	
	if keyword_set ( numerical ) then begin

		YkmBFns	= YkmBFns2
		dYkmDthBFns	= dYkmDthBFns2
		dYkmDphBFns	= dYkmDphBFns2
		outNkValues	= outNkValues2

	endif


;	Fit |dB| to dB.grad Ykm in the GEI system
;	-----------------------------------------	
	
	dBMag   = sqrt(data.dBTheta^2+data.dBPhi^2)
	kTh     = transpose(rebin(data.dBTheta/dBMag,$
	       n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))
	kPh     = transpose(rebin(data.dBPhi/dBMag,$
	        n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))
	
	bFuncs  = kTh*dYkmDphBFns + kPh*dYkmDthBFns
	alpha   = transpose(bFuncs) ## bFuncs
	beta_   = transpose( transpose(bFuncs) ## dBMag)
	la_svd, alpha, w, u, v, status=svdStatus
	kk      = where(w lt max ( w ) * 1d-5, kkcnt)
	if kkcnt gt 0 then w[kk]=0.0
	coeffs  = svsol(u, w, v, beta_)
	
	fit_dBTheta     = dYkmDPhBFns ## coeffs
	fit_dBPhi       = dYkmDThBFns ## coeffs
	fit_dBMag       = bFuncs ## coeffs
	
	Re=6.371d6
	R     = Re + 780.0d3
	u0    = 4.0*!dpi*1.0d-7
	jpar  = (YkmBFns ## $
	        (coeffs*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}

;	generate basis fns at regular grid
;	----------------------------------

	schaBasisFunctions, kMax, mMax, capSize, geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
         YkmBFns=YkmBFns_grid, $
		 dYkmDthBFns=dYkmDthBFns_grid, $
		 dYkmDphBFns=dYkmDphBFns_grid, $
		 pnmPath = pnmPath, /evenSet

	ampere_setupSHFns, 1.0, geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns_grid2, $
			bPhBFnArr = dYkmDPhBFns_grid2, $
			YBFnArr = YkmBFns_grid2, $
			mArr = outMValues_grid2, $
			nArr = outNkValues_grid2, $
			kArr = outKValues_grid2

	if keyword_set ( numerical ) then begin

		YkmBFns_grid	= YkmBFns_grid2
		dYkmDthBFns_grid	= dYkmDthBFns_grid2
		dYkmDphBFns_grid	= dYkmDphBFns_grid2

	endif

  	jParAACGM	= (YkmBFns_grid ## $
         (coeffs*(-outNkValues*(outNkValues+1.0))))$
                         /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}
	jParAACGM	= reform ( jParAACGM, nLatGrid, nLonGrid )

   	dBTheta_GEI_grid     = dYkmDPhBFns_grid ## coeffs
 	dBPhi_GEI_grid       = dYkmDThBFns_grid ## coeffs
 
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

	window, winnum, xSize = 600, ySize = 600
	winNum++

	plt_dat, data.gei_coLat_rad * !radeg, data.gei_lon_rad * !radeg, $
		n_elements ( data.gei_coLat_rad ), -data.dbTheta, data.dbPhi, [1,1], [1,2], $
		title = 'GEI - Raw Data'
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
	jLevels = ( fIndGen ( 21 ) - 10 ) * 0.5;2.0
	colors  = bytScl ( jLevels, top = 253 ) + 1
	oldbck=!p.background
	!p.background=0
	jpar=reform(jpar)
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
	jParTmp	= jParAACGM[*]
	lonTmp = aacgmGrid_lon_deg[*] $
			+ randomu ( sysTime(/sec ), n_elements(jParTmp), /uni ) * 1e-5
	latTmp = 90.0-aacgmGrid_coLat_deg[*]
	contour, jParTmp, $
			lonTmp, $
			latTmp, $
		 	c_labels=fltarr(n_elements(jLevels))+1,$
        	/over, levels = jLevels, $
           	c_colors = colors, /fill, /irreg
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
