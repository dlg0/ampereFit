pro f_plot_amperefit

    iFName = '/Users/dg6/ampere/20100824Amp_invert.ncdf'
    oFName = '/Users/dg6/code/ampereFit/fortran/ampData_original.nc'
    sFName = '/Users/dg6/code/ampereFit/fortran/ampData_shifted.nc'

    cdfId = ncdf_open ( iFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'b_eci', in_b_eci

	ncdf_close, cdfId

    cdfId = ncdf_open ( oFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'data_bR_GEI', orig_bR_GEI 
	    ncdf_varGet, cdfId, 'data_bT_GEI', orig_bT_GEI 
	    ncdf_varGet, cdfId, 'data_bP_GEI', orig_bP_GEI 

	    ncdf_varGet, cdfId, 'data_bx_GEI', orig_bx_GEI 
	    ncdf_varGet, cdfId, 'data_by_GEI', orig_by_GEI 
	    ncdf_varGet, cdfId, 'data_bz_GEI', orig_bz_GEI 

	    ncdf_varGet, cdfId, 'data_GEI_coLat_rad', orig_GEI_coLat_rad
	    ncdf_varGet, cdfId, 'data_GEI_lon_rad', orig_GEI_lon_rad
	    ncdf_varGet, cdfId, 'data_GEI_R_km', orig_GEI_R_km

        ncdf_varGet, cdfId, 'data_GEI_px', orig_GEI_px
        ncdf_varGet, cdfId, 'data_GEI_py', orig_GEI_py
        ncdf_varGet, cdfId, 'data_GEI_pz', orig_GEI_pz

	ncdf_close, cdfId

    cdfId = ncdf_open ( sFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'data_bR_GEI', shif_bR_GEI 
	    ncdf_varGet, cdfId, 'data_bT_GEI', shif_bT_GEI 
	    ncdf_varGet, cdfId, 'data_bP_GEI', shif_bP_GEI 

	    ncdf_varGet, cdfId, 'data_bx_GEI', shif_bx_GEI 
	    ncdf_varGet, cdfId, 'data_by_GEI', shif_by_GEI 
	    ncdf_varGet, cdfId, 'data_bz_GEI', shif_bz_GEI 

	    ncdf_varGet, cdfId, 'data_GEI_coLat_rad', shif_GEI_coLat_rad
	    ncdf_varGet, cdfId, 'data_GEI_lon_rad', shif_GEI_lon_rad
	    ncdf_varGet, cdfId, 'data_GEI_R_km', shif_GEI_R_km

        ncdf_varGet, cdfId, 'data_GEI_px', shif_GEI_px
        ncdf_varGet, cdfId, 'data_GEI_py', shif_GEI_py
        ncdf_varGet, cdfId, 'data_GEI_pz', shif_GEI_pz

	    ncdf_varGet, cdfId, 'fit_bT_GEI', shif_fit_bT_GEI 
	    ncdf_varGet, cdfId, 'fit_bP_GEI', shif_fit_bP_GEI 

	ncdf_close, cdfId


    iiNorth = where ( orig_GEI_coLat_rad*!radeg gt 0 and orig_GEI_pz gt 0 )
    dlg_plot_vecs, orig_GEI_R_km[iiNorth], orig_GEI_coLat_rad[iiNorth], orig_GEI_lon_rad[iiNorth], $
        orig_bT_GEI[iiNorth], orig_bP_GEI[iiNorth]

    iiNorth = where ( shif_GEI_pz gt 0 )
    dlg_plot_vecs, shif_GEI_R_km[iiNorth], shif_GEI_coLat_rad[iiNorth], shif_GEI_lon_rad[iiNorth], $
        shif_bT_GEI[iiNorth], shif_bP_GEI[iiNorth]

    dlg_plot_vecs, shif_GEI_R_km[iiNorth], shif_GEI_coLat_rad[iiNorth], shif_GEI_lon_rad[iiNorth], $
        shif_fit_bT_GEI[iiNorth], shif_fit_bP_GEI[iiNorth]

stop

end
