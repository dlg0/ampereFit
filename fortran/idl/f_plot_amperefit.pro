pro f_plot_amperefit

    iFName = '/Users/dg6/ampere/20100824Amp_invert.ncdf'
    oFName = '/Users/dg6/code/ampereFit/fortran/ampData.nc'

    cdfId = ncdf_open ( iFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'b_eci', in_b_eci

	ncdf_close, cdfId

    cdfId = ncdf_open ( oFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'data_bR_GEI', data_bR_GEI 
	    ncdf_varGet, cdfId, 'data_bT_GEI', data_bT_GEI 
	    ncdf_varGet, cdfId, 'data_bP_GEI', data_bP_GEI 

	    ncdf_varGet, cdfId, 'data_bx_GEI', data_bx_GEI 
	    ncdf_varGet, cdfId, 'data_by_GEI', data_by_GEI 
	    ncdf_varGet, cdfId, 'data_bz_GEI', data_bz_GEI 

	    ncdf_varGet, cdfId, 'fit_bT_GEI', fit_bT_GEI 
	    ncdf_varGet, cdfId, 'fit_bP_GEI', fit_bP_GEI 

	    ncdf_varGet, cdfId, 'data_GEI_coLat_rad', data_GEI_coLat_rad
	    ncdf_varGet, cdfId, 'data_GEI_lon_rad', data_GEI_lon_rad

        ncdf_varGet, cdfId, 'data_GEI_px', data_GEI_px
        ncdf_varGet, cdfId, 'data_GEI_py', data_GEI_py
        ncdf_varGet, cdfId, 'data_GEI_pz', data_GEI_pz

	ncdf_close, cdfId

    iiNorth = where ( data_GEI_pz gt 0 )
    myVec = vector ( data_bx_GEI[iiNorth], data_by_GEI[iiNorth], $
            data_GEI_px[iiNorth], data_GEI_py[iiNorth], $
            grid_units = 'degrees', $
            auto_color = 1, $
            rgb_table = 1, $
            length = 2, $
            head_size = 0.3, $
            data_location = 0, $
            /iso, $
            xRange = [-4e3,4e3], yRange = [-4e3,4e3] ) 

    myMap = map ( 'STEREOGRAPHIC', $
            center_latitude = 90, $
            center_longitude = 0, $
            latitude_min = 50 )

    myVec = vector ( data_bP_GEI, -data_bT_GEI*0, $
            data_GEI_lon_rad*!radeg-180, 90-data_GEI_coLat_rad*!radeg, $
            /overPlot, $
            grid_units = 'degrees', $
            auto_color = 1, $
            rgb_table = 1, $
            length = 2, $
            head_size = 0.3, $
            data_location = 0 ) 

    myMap = map ( 'STEREOGRAPHIC', $
            center_latitude = 90, $
            center_longitude = 0, $
            latitude_min = 50 )

    myVec = vector ( fit_bP_GEI, -fit_bT_GEI, $
            data_GEI_lon_rad*!radeg-180, 90-data_GEI_coLat_rad*!radeg, $
            /overPlot, $
            grid_units = 'degrees', $
            auto_color = 1, $
            rgb_table = 1, $
            length = 2, $
            head_size = 0.3, $
            data_location = 0 ) 


stop
	

end
