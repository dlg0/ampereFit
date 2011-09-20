pro dlg_plot_vecs, r, th, ph, bTh, bPh, title = title

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

    myVec = vector ( xVec, yVec, $
            xPos, yPos, $
            grid_units = 'degrees', $
            auto_color = 1, $
            rgb_table = 1, $
            length = 2, $
            head_size = 0.3, $
            data_location = 0, $
            aspect_ratio = 1, $
            xRange = [-meanr,meanr], yRange = [-meanr,meanr], title = title ) 

    for i=0,n_elements(coLatLines)-1 do begin

        lons    = fIndGen(nPtsLat)/(nPtsLat-1)*360
        lats    = fltArr(nPtsLat)+coLatLines[i]

        x   = meanr * sin ( lats * !dtor ) * cos ( lons * !dtor )
        y   = meanr * sin ( lats * !dtor ) * sin ( lons * !dtor )

        p = plot ( x, y, /over, color = 'black' )

    endfor

    for i=0,n_elements(lonLines)-1 do begin
        lons    = fltArr(nPtsLon)+lonLines[i]
        lats    = fIndGen(nPtsLon)/(nPtsLon-1)*90

        x   = meanr * sin ( lats * !dtor ) * cos ( lons * !dtor )
        y   = meanr * sin ( lats * !dtor ) * sin ( lons * !dtor )

        p = plot ( x, y, /over, color = 'black' )

    endfor

end
