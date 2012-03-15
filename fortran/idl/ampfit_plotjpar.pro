pro ampfit_plotjpar

	;filelist = ['output/ampData_gridFitCurlAACGM.nc','output/ampData_gridFitAACGM.nc'] 
	GrdFileList = file_search ( 'output/*grd.ncdf' )

	for f=0,n_elements(GrdFileList)-1 do begin

    	cdfid = ncdf_open ( GrdFileList[f], /nowrite ) 


			nCdf_varGet, cdfId, 'year', year
			nCdf_varGet, cdfId, 'month', month
			nCdf_varGet, cdfId, 'day', day
			nCdf_varGet, cdfId, 'hour', hour
			nCdf_varGet, cdfId, 'minut', minut
			nCdf_varGet, cdfId, 'sec', sec
			nCdf_varGet, cdfId, 'avgint', avgint

			;nCdf_varGet, cdfId, 'kmax', kmax
			;nCdf_varGet, cdfId, 'mmax', mmax

			;nCdf_varGet, cdfId, 'res_deg', res_deg
			;nCdf_varGet, cdfId, 'nLatGrid', nLatGrid
			;nCdf_varGet, cdfId, 'nLonGrid', nLonGrid

			nCdf_varGet, cdfId, 'cLat_deg', cLat_deg
			nCdf_varGet, cdfId, 'lon_hr', lon_hr
			nCdf_varGet, cdfId, 'dbth_th', dBth_th
			nCdf_varGet, cdfId, 'dbph_th', dBph_th
			nCdf_varGet, cdfId, 'dbth_ph', dBth_ph
			nCdf_varGet, cdfId, 'dbph_ph', dBph_ph

		    ;ncdf_varget, cdfid, 'bR', bR 
		    ;ncdf_varget, cdfid, 'bT', bT 
		    ;ncdf_varget, cdfid, 'bP', bP 

		    ;ncdf_varget, cdfid, 'bX', bX 
		    ;ncdf_varget, cdfid, 'bY', bY 
		    ;ncdf_varget, cdfid, 'bZ', bZ 

		    ;ncdf_varget, cdfid, 'T', T 
		    ;ncdf_varget, cdfid, 'P', P
		    ;ncdf_varget, cdfid, 'R', R

    	    ;ncdf_varget, cdfid, 'X', X
    	    ;ncdf_varget, cdfid, 'Y', Y
    	    ;ncdf_varget, cdfid, 'Z', Z

    	    ncdf_varget, cdfid, 'jPar', jPar 

		ncdf_close, cdfId


		; Get the X,Y coords
		rE = 6371.0
		rI = rE + 780.0
		scaleFac = 1/90.0*rI
		x = scaleFac * cLat_deg * cos ( lon_hr*15.0*!dtor)
		y = scaleFac * cLat_deg * sin ( lon_hr*15.0*!dtor)

		scale = 1.0;max(abs(jPar))*0.8
		print, 'jPar scale: ', scale, ' uAm^-2'
		nLevs = 6
		levels = (1+fIndGen(nLevs))/(nLevs)*scale
		colors = reverse(bytScl ( levels, top = 253 ) + 1)

		set_plot, 'ps'
		device, fileName = strJoin(['output/jPar_',file_basename(GrdFileList[f],'.nc'),'.eps']), $
				/encaps, xSize = 6, ySize = 6, /color, /inches
		loadct, 0

		meanr = 3000.0

    	coLatLines = [10,20,30,40,50,60,70,80,90]
    	lonLines = [0,90,180,270]
    	nPtsLat = 360
    	nPtsLon = 20
		print, 'Mean(r): ',meanr 
		title=file_basename(GrdFileList[f],'.nc')

		loadct, 1
		contour, jPar, x+randomn(1,n_elements(x))*1e-6, y, $
				/irreg, $
				levels = levels, $
				c_colors = colors, $
				position = [0.1,0.1,0.9,0.9], /norm, $
				xStyle = 1, yStyle = 1, $
				yRange = [-meanr, meanr], xRange=[-meanr,meanr], $
				c_labels = colors*0+1, $
				title = title, $
				c_charSize = 0.3, /fill
		loadct, 3
		contour, -jPar, x+randomn(1,n_elements(x))*1e-6, y, $
				/irreg, $
				levels = levels, $
				c_colors = colors, $
				c_labels = colors*0+1, $
				/over, $
				c_charSize = 0.3, /fill

		xyouts, 0.7, 0.95, 'Max: '+string(max(jPar)), /norm, charSize = 0.8
		xyouts, 0.7, 0.92, 'Min: '+string(min(jPar)), /norm, charSize = 0.8

    	for i=0,n_elements(coLatLines)-1 do begin

    	    lons    = fIndGen(nPtsLat)/(nPtsLat-1)*360
    	    lats    = fltArr(nPtsLat)+coLatLines[i]

    	    x   = meanr * sin ( lats * !dtor ) * cos ( lons * !dtor )
    	    y   = meanr * sin ( lats * !dtor ) * sin ( lons * !dtor )

    	    plots,  x, y, color = 0 

    	endfor

    	for i=0,n_elements(lonLines)-1 do begin

    	    lons    = fltArr(nPtsLon)+lonLines[i]
    	    lats    = fIndGen(nPtsLon)/(nPtsLon-1)*90

    	    x   = meanr * sin ( lats * !dtor ) * cos ( lons * !dtor )
    	    y   = meanr * sin ( lats * !dtor ) * sin ( lons * !dtor )

    	    plots, x, y, color = 0 

    	endfor

		p = plot(jPar[*])
		device, /close

    	ampfit_plotvecs, rI+cLat_deg*0, cLat_deg*!dtor, lon_hr*15*!dtor, $
    	    dBth_th, dBph_ph, title = file_basename(GrdFileList[f],'.nc'), $
			fileNameIn=strJoin(['output/',file_basename(GrdFileList[f],'.nc'),'.eps'])


	endfor

end
