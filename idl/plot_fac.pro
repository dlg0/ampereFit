pro plot_fac, jPar, pole, cap_coLat_deg, coLat_deg, lon_deg, $
		mlt_shift = mlt_shift, $
		title = title, $
		winNo = winNo, $
		south = south

	if not keyword_set(mlt_shift) then mlt_shift = 0
	if not keyword_set(winNo) then winNo = 0

	;loadct, 13, file = path + 'davect.tbl',/silent
	loadct, 13, file = 'davect.tbl',/silent

	window, winNo
   	map_set, pole, 0, 0, /ortho, /iso, $
     limit = [ pole - cap_coLat_deg, 0, pole, 360 ], xmargin=[1,1], ymargin=[1,10], $
	 /noborder, /advance, title = title

	jParTmp	= jPar[*]

	;jParTmp = jParTmp>mn_fac
	;jParTmp = jParTmp<mx_fac

	lonTmp = ((lon_deg[*]/15.0+mlt_shift) mod 24 )*15.0 $
			+ randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*1e-5-0.5e-5

	latTmp = 90.0-coLat_deg[*]


	facScale = 1.0

	; Contour positive and negative seperately.
	; This prevents zero from looking like one color.
	; Only a problem for filled contours.

	; Positive contours
	; -----------------

	nLevs = 10
	jLevels = (fIndGen(nLevs)+1)/nLevs * facScale
	colors  = (255-bytScl ( jLevels, top = 253 ) )

	loadct, 3, /silent

	iiPos = where ( jParTmp gt 0 )
	thisJPar = jParTmp[iiPos]
	thisLon = lonTmp[iiPos]
	thisLat = latTmp[iiPos]

   	contour, thisJPar, $
			thisLon, $
			thisLat, $
        	/over, levels = jLevels, $
           	c_colors = colors, /irreg;, /cell_fill

	; Negative contours
	; -----------------

	loadct, 1, /silent

	jParTmp = -jParTmp
	iiPos = where ( jParTmp gt 0 )
	thisJPar = jParTmp[iiPos]
	thisLon = lonTmp[iiPos]
	thisLat = latTmp[iiPos]

   	contour, thisJPar, $
			thisLon, $
			thisLat, $
        	/over, levels = jLevels, $
           	c_colors = colors, /irreg;, /cell_fill


	; Just grid both north and south hemispheres
	; ------------------------------------------

	loadct, 0, /silent

	sgn = 1
	map_grid, label = 1, $
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360], $
			lonlab=10, $
			latLab = 45, $
			lats = (fIndGen(16)+1)*10, $
			color = 0
stop
end
