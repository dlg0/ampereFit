pro ampfit_plotvecs, r, th, ph, bTh, bPh, title = title, fileNameIn=fileNameIn

    xVec    = fltArr ( n_elements(bTh) )
    yVec    = fltArr ( n_elements(bTh) )

    xPos     = fltArr ( n_elements(bTh) )
    yPos     = fltArr ( n_elements(bTh) )
    
    rot_    = fltArr(2,2)

   for i=0,n_elements(bTh)-1 do begin
       
        rot_[0,0]   = cos(ph[i])
        rot_[1,0]   = sin(ph[i])
        rot_[0,1]   = -sin(ph[i])
        rot_[1,1]   = cos(ph[i]) 

        xPos[i] = r[i]*sin(th[i]) * cos(ph[i])
        yPos[i] = r[i]*sin(th[i]) * sin(ph[i])

        rotVec  = rot_ # [ bTh[i],bPh[i] ] 

        xVec[i] = rotVec[0]
        yVec[i] = rotVec[1]

    endfor

    coLatLines = [10,20,30,40,50,60,70,80,90]
    lonLines = [0,90,180,270]
    nPtsLat = 360
    nPtsLon = 20
    meanr   = mean ( r )
	print, 'Mean(r): ', mean(r)

    ;myVec = vector ( xVec, yVec, $
    ;        xPos, yPos, $
    ;        auto_color = 1, $
    ;        rgb_table = 1, $
    ;        head_size = 0.3, $
	;		length = xVec, $
    ;        data_location = 0, $
    ;        aspect_ratio = 1, $
    ;        xRange = [-meanr,meanr], yRange = [-meanr,meanr], title = title, /noData )


	set_plot, 'ps'
	device, fileName = fileNameIn, /encaps, xSize = 6, ySize = 6, /color, /inches
	loadct, 0
	plot, xPos, yPos, xRange = [-meanr,meanr], yRange = [-meanr,meanr], title = title, $
			position = [0.1,0.1,0.9,0.9], /norm, color = 0, xStyle=1, yStyle=1, psym = 4, symsize=0.2

	loadct, 1
	sf = 4
	for i=0,n_elements(xVec)-1 do begin
			x0 = xPos[i]
			y0 = yPos[i]
			x1 = x0+xVec[i]*sf
			y1 = y0+yVec[i]*sf
			plots, [x0,x1], [y0,y1], color = (255-sqrt(xVec[i]^2+yVec[i]^2))>1
	endfor

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

end
