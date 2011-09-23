pro ampfit_plotjpar

	filelist = ['output/ampData_gridFitCurlAACGM.nc','output/ampData_gridFitAACGM.nc'] 

	for f=0,n_elements(filelist)-1 do begin

    	cdfid = ncdf_open ( filelist[f], /nowrite ) 

		    ncdf_varget, cdfid, 'bR', bR 
		    ncdf_varget, cdfid, 'bT', bT 
		    ncdf_varget, cdfid, 'bP', bP 

		    ncdf_varget, cdfid, 'bX', bX 
		    ncdf_varget, cdfid, 'bY', bY 
		    ncdf_varget, cdfid, 'bZ', bZ 

		    ncdf_varget, cdfid, 'T', T 
		    ncdf_varget, cdfid, 'P', P
		    ncdf_varget, cdfid, 'R', R

    	    ncdf_varget, cdfid, 'X', X
    	    ncdf_varget, cdfid, 'Y', Y
    	    ncdf_varget, cdfid, 'Z', Z

    	    ncdf_varget, cdfid, 'jPar', jPar 

		ncdf_close, cdfid

		scale = 0.5;max(abs(jPar))*0.8
		print, 'jPar scale: ', scale, ' uAm^-2'
		nLevs = 11
		levels = (1+fIndGen(nLevs))/(nLevs)*scale
		colors = reverse(bytScl ( levels, top = 253 ) + 1)

		set_plot, 'ps'
		device, fileName = strJoin(['output/jPar_',file_basename(fileList[f],'.nc'),'.eps']), $
				/encaps, xSize = 6, ySize = 6, /color, /inches
		loadct, 0

    	coLatLines = [10,20,30,40,50,60,70,80,90]
    	lonLines = [0,90,180,270]
    	nPtsLat = 360
    	nPtsLon = 20
    	meanr   = mean ( r )
		print, 'Mean(r): ', mean(r)
		title=file_basename(fileList[f],'.nc')

		loadct, 1
		contour, jPar, x+randomn(1,n_elements(x))*1e-6, y, $
				/irreg, $
				levels = levels, $
				/fill, $
				c_colors = colors, $
				position = [0.1,0.1,0.9,0.9], /norm, $
				xStyle = 1, yStyle = 1, $
				yRange = [-meanr, meanr], xRange=[-meanr,meanr], $
				title = title
		loadct, 3
		contour, -jPar, x+randomn(1,n_elements(x))*1e-6, y, $
				/irreg, $
				levels = levels, $
				c_colors = colors, $
				/fill, $
				/over

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


		device, /close

	endfor

end
