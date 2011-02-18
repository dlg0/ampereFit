pro f_plot_amperefit

    iFName = '/Users/dg6/ampere/20100824Amp_invert.ncdf'
    oFName = '/Users/dg6/code/ampereFit/fortran/jParIrid.nc'

    cdfId = ncdf_open ( iFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'b_eci', in_b_eci

	ncdf_close, cdfId

    cdfId = ncdf_open ( oFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'data_bR_GEI', data_bR_GEI 
	    ncdf_varGet, cdfId, 'data_bT_GEI', data_bT_GEI 
	    ncdf_varGet, cdfId, 'data_bP_GEI', data_bP_GEI 

	    ncdf_varGet, cdfId, 'fit_bT_GEI', fit_bT_GEI 
	    ncdf_varGet, cdfId, 'fit_bP_GEI', fit_bP_GEI 

	    ncdf_varGet, cdfId, 'data_GEI_coLat_rad', data_GEI_coLat_rad
	    ncdf_varGet, cdfId, 'data_GEI_lon_rad', data_GEI_lon_rad

	ncdf_close, cdfId


    myMap = map ( 'STEREOGRAPHIC', $
            center_latitude = 90, $
            center_longitude = 0, $
            latitude_min = 70 )

    myVec = vector ( data_bP_GEI, -data_bT_GEI, $
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
            latitude_min = 70 )

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
