pro new_pole_rotate, newPole_coLat_rad, newPole_lon_rad, $
    coLat_rad, lon_rad, bTh, bPh, $
    coLat_rad_out = coLat_rad_1, lon_rad_out = lon_rad_1, $
    bTh_out = bTheta_1, bPh_out = bPhi_1

    nPts    = n_elements ( coLat_rad )
    r       = 1

    ;   get cartesian pole location

    xPole   = r * sin ( newPole_coLat_rad ) * cos ( newPole_lon_rad )
    yPole   = r * sin ( newPole_coLat_rad ) * sin ( newPole_lon_rad )
    zPole   = r * cos ( newPole_coLat_rad )

    ;   get data cartesian locations

    x   = r * sin ( coLat_rad ) * cos ( lon_rad )
    y   = r * sin ( coLat_rad ) * sin ( lon_rad )
    z   = r * cos ( coLat_rad )

    ;   get cartesian vectors

    bx  = fltArr ( nPts )
    by  = fltArr ( nPts )
    bz  = fltArr ( nPts )
 
    rotMat  = fltArr ( 3, 3, nPts )
     
    for i = 0, n_elements ( coLat_rad ) - 1 do begin

        th  = colat_rad[i]
        ph  = lon_rad[i]
        rotMat  = [ [ sin(th)*cos(ph), sin(th)*sin(ph), cos(th) ], $
                    [ r*cos(th)*cos(ph), r*cos(th)*sin(ph), -r*sin(th) ], $
                    [ -r*sin(ph), r*cos(ph), 0 ] ]

        cartVec  = transpose ( rotMat ) ## [    [0], $
                                                       [bTh[i]], $
                                                       [bPh[i]] ]

        bx[i]   = cartVec[0]
        by[i]   = cartVec[1]
        bz[i]   = cartVec[2]

    endfor

    ;   we want to rotate the z axis to be parallel to the vector from the origin
    ;   to the point defined by the [xyz]Pole variables

    angle_y_to_z    = !pi/2-atan ( zPole, yPole ) 
    rot_y_to_z  = [ [ 1, 0, 0 ], $
                    [ 0, cos ( angle_y_to_z ), -sin ( angle_y_to_z) ], $
                    [ 0, sin ( angle_y_to_z ), cos ( angle_y_to_z ) ] ] 

    angle_z_to_x    = atan ( xPole, zPole ) 
    rot_z_to_x  = [ [ cos ( angle_z_to_x ), 0, sin ( angle_z_to_x ) ], $
                    [ 0, 1, 0 ], $
                    [ -sin ( angle_z_to_x ), 0, cos ( angle_z_to_x ) ] ] 

    angle_x_to_y    = atan ( yPole, xPole ) 
    rot_x_to_y  = [ [ cos ( angle_x_to_y ), -sin ( angle_x_to_y ), 0 ], $
                    [ sin ( angle_x_to_y ), cos ( angle_x_to_y ), 0 ], $
                    [ 0, 0, 1 ] ] 

    bx1 = fltArr ( nPts )
    by1 = fltArr ( nPts )
    bz1 = fltArr ( nPts )
 
    x1  = fltArr ( nPts )
    y1  = fltArr ( nPts )
    z1  = fltArr ( nPts )

    ;   rotate both the db vectors and the position vectors
    ;   of the coordinates

    totalRotMat = rot_z_to_x ## rot_x_to_y ## rot_y_to_z
 
    for i = 0, nPts - 1 do begin

        bRot1   = totalRotMat ## [   [ bx[i] ], $
                                    [ by[i] ], $
                                    [ bz[i] ] ]
        bx1[i]  = bRot1[0]
        by1[i]  = bRot1[1]
        bz1[i]  = bRot1[2]

        rot1   = totalRotMat ## [   [ x[i] ], $
                                    [ y[i] ], $
                                    [ z[i] ] ]
        x1[i]  = rot1[0]
        y1[i]  = rot1[1]
        z1[i]  = rot1[2]

    endfor

    ;   convert new x,y,z coords back to spherical

    r1   = sqrt ( x1^2 + y1^2 + z1^2 )
    coLat_rad_1 = aCos ( z1 / ( r1 ) )
    lon_rad_1   = atan ( y1, x1 )  
    iiNegLon    = where ( lon_rad_1 lt 0 ) 
    lon_rad_1 += 2*!pi
 
    ;   rotate vectors back to spherical using the new 
    ;   spherical coords
 
    bR_1 = fltArr ( nPts )
    bTheta_1 = fltArr ( nPts )
    bPhi_1 = fltArr ( nPts )
 
    for i = 0, nPts - 1 do begin

        th  = colat_rad_1[i]
        ph  = lon_rad_1[i]
        rotMat  = [ [ sin(th)*cos(ph), sin(th)*sin(ph), cos(th) ], $
                    [ r*cos(th)*cos(ph), r*cos(th)*sin(ph), -r*sin(th) ], $
                    [ -r*sin(ph), r*cos(ph), 0 ] ]


        sphVec  = rotMat ## [    [bx1[i]], $
                                        [by1[i]], $
                                        [bz1[i]] ]
        
        bR_1[i]   = sphVec[0]
        bTheta_1[i]   = sphVec[1]
        bPhi_1[i] = sphVec[2]
    
    endfor

   stop

end
