pro plot_vec, coLat_rad, lon_rad, bTh=bTh, bPh=bPh, $
        capSize = capSize, $
        title = title, $
        iTool = iTool, $
        winNo = winNo, $
        ptsOnly  = ptsOnly

    if not keyword_set ( title ) then title = 'good morning dave ;-)'
    if not keyword_set ( capSize ) then capSize = 50
    if keyword_set ( winNo ) then window, winNo

    if keyword_set ( iTool ) then begin

        iMap, map_projection = 'Orthographic', $
          	limit = [ 90.0 - capSize, -180, 90, 180 ], $
            center_lat = 90, $
            center_lon = -90, $
            id = map
           
        sf =    sqrt(bTh^2+bPh^2) $
                    / sqrt(bTh^2+(bPh*sin(coLat_rad))^2)

        plotLon_rad    = lon_rad
        iiPi    = where ( plotLon_rad ge !pi )
        plotLon_rad[iiPi] -= 2*!Pi

        bMag    = sqrt ( bTh^2 + bPh^2 )
        color   = 255 - ( bytScl ( bMag, min = 0, max = 500, top = 253 ) + 1 )

        iVector, -bTh, -bPh * sin ( coLat_rad ), $
                plotLon_rad * !radeg, 90 - coLat_rad * !radeg, $
                /over, $
                data_location = 0, $
                rgb_table = 1, $
                vector_colors = color, $
                length_scale = 5
    endif else begin

        map_set, 90, 0, 0, $
    	    	/ortho, $
       	    	/iso, $
            	limit = [ 90.0 - ( capSize + 2 ), 0, 90, 360 ], $
       	    	/noborder, $
       	    	/advance, $
                title = title 
    
        map_grid, label = 1, latDel = 10.0


        ;   get the th and ph unit vectors in an XY system assuming the 
        ;   map projection is a polar coordinate system

        x   = coLat_rad * !radeg * cos ( lon_rad )
        y   = coLat_rad * !radeg * sin ( lon_rad )

        thStep  = 0.01; [deg]
        th_x    = ( coLat_rad * !radeg + thStep ) * cos ( lon_rad )
        th_y    = ( coLat_rad * !radeg + thStep ) * sin ( lon_rad )

        phStep  = 0.01; [deg]
        ph_x    = coLat_rad * !radeg * cos ( lon_rad + phStep * !dtor )
        ph_y    = coLat_rad * !radeg * sin ( lon_rad + phStep * !dtor )

        unitTh_x   = ( th_x - x ) / sqrt ( (th_x-x)^2 + (th_y-y)^2 )
        unitTh_y   = ( th_y - y ) / sqrt ( (th_x-x)^2 + (th_y-y)^2 )

        unitPh_x   = ( ph_x - x ) / sqrt ( (ph_x-x)^2 + (ph_y-y)^2 )
        unitPh_y   = ( ph_y - y ) / sqrt ( (ph_x-x)^2 + (ph_y-y)^2 )

        b_x = bTh * unitTh_x + bPh * unitPh_x
        b_y = bTh * unitTh_y + bPh * unitPh_y

        bMag_1  = sqrt ( bTh^2 + bPh^2 )
        bMag_2  = sqrt ( b_x^2 + b_y^2 )

        ;plot, bMag_1
        ;oPlot, bMag_2

        scale   = 50

        endX    = x + b_x / scale
        endY    = y + b_y / scale

        end_coLat_rad   = sqrt ( endX^2 + endY^2 ) * !dtor
        end_lon_rad   = atan ( endY, endX )  

        if keyword_set ( ptsOnly ) then begin

            plots, lon_rad*!radeg, 90-coLat_rad*!radeg, psym = 4

        endif else begin

            for i = 0, n_elements ( x ) - 1 do begin

                points  = [ [lon_rad[i]*!radeg,end_lon_rad[i]*!radeg],$ 
                            [90-coLat_rad[i]*!radeg, 90-end_coLat_rad[i]*!radeg] ]
                plots,  [lon_rad[i]*!radeg,end_lon_rad[i]*!radeg], $
                    [90-coLat_rad[i]*!radeg, 90-end_coLat_rad[i]*!radeg]
                plots,  [lon_rad[i]*!radeg,end_lon_rad[i]*!radeg], $
                    [90-coLat_rad[i]*!radeg, 90-end_coLat_rad[i]*!radeg]


            endfor

        endelse

    endelse    

end
