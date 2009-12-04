pro plot_vec, coLat_rad, lon_rad, bTh, bPh, $
        capSize = capSize, $
        title = title

    if not keyword_set ( title ) then title = 'good morning dave ;-)'
    if not keyword_set ( capSize ) then capSize = 50

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

    thStep  = 0.1; [deg]
    th_x    = ( coLat_rad * !radeg + thStep ) * cos ( lon_rad )
    th_y    = ( coLat_rad * !radeg + thStep ) * sin ( lon_rad )

    phStep  = 0.1; [deg]
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

    plot, bMag_1
    oPlot, bMag_2

    endX    = x + b_x 
    endY    = y + b_y 

    end_coLat_rad   = sqrt ( endX^2 + endY^2 ) * !dtor
    end_lon_rad   = atan ( endY, endX )  

    for i = 0, n_elements ( x ) - 1 do begin

        plots, [lon_rad[i],end_lon_rad[i]], [coLat_rad[i], end_coLat_rad[i]] 

    endfor

    stop
end
