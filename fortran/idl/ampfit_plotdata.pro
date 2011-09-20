pro ampfit_plotdata

    iFName = 'data/20100824Amp_invert.ncdf'
    rawFName = 'output/ampData_original.nc'
    shiftedFName = 'output/ampData_shifted.nc'
	ghostedFName = 'output/ampData_ghosted.nc'
	ghostedFitFName = 'output/ampData_ghosted_fit.nc'

    cdfId = ncdf_open ( iFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'b_eci', in_b_eci

	ncdf_close, cdfId

    cdfId = ncdf_open ( rawFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_raw_GEI 
	    ncdf_varGet, cdfId, 'bT', bT_raw_GEI 
	    ncdf_varGet, cdfId, 'bP', bP_raw_GEI 

	    ncdf_varGet, cdfId, 'bX', bx_raw_GEI 
	    ncdf_varGet, cdfId, 'bY', by_raw_GEI 
	    ncdf_varGet, cdfId, 'bZ', bz_raw_GEI 

	    ncdf_varGet, cdfId, 'T', T_raw_GEI 
	    ncdf_varGet, cdfId, 'P', P_raw_GEI
	    ncdf_varGet, cdfId, 'R', R_raw_GEI

        ncdf_varGet, cdfId, 'X', X_raw_GEI
        ncdf_varGet, cdfId, 'Y', Y_raw_GEI
        ncdf_varGet, cdfId, 'Z', Z_raw_GEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( shiftedFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_shifted_GEI 
	    ncdf_varGet, cdfId, 'bT', bT_shifted_GEI 
	    ncdf_varGet, cdfId, 'bP', bP_shifted_GEI 

	    ncdf_varGet, cdfId, 'bX', bx_shifted_GEI 
	    ncdf_varGet, cdfId, 'bY', by_shifted_GEI 
	    ncdf_varGet, cdfId, 'bZ', bz_shifted_GEI 

	    ncdf_varGet, cdfId, 'T', T_shifted_GEI 
	    ncdf_varGet, cdfId, 'P', P_shifted_GEI
	    ncdf_varGet, cdfId, 'R', R_shifted_GEI

        ncdf_varGet, cdfId, 'X', X_shifted_GEI
        ncdf_varGet, cdfId, 'Y', Y_shifted_GEI
        ncdf_varGet, cdfId, 'Z', Z_shifted_GEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( ghostedFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_ghosted_GEI 
	    ncdf_varGet, cdfId, 'bT', bT_ghosted_GEI 
	    ncdf_varGet, cdfId, 'bP', bP_ghosted_GEI 

	    ncdf_varGet, cdfId, 'bX', bx_ghosted_GEI 
	    ncdf_varGet, cdfId, 'bY', by_ghosted_GEI 
	    ncdf_varGet, cdfId, 'bZ', bz_ghosted_GEI 

	    ncdf_varGet, cdfId, 'T', T_ghosted_GEI 
	    ncdf_varGet, cdfId, 'P', P_ghosted_GEI
	    ncdf_varGet, cdfId, 'R', R_ghosted_GEI

        ncdf_varGet, cdfId, 'X', X_ghosted_GEI
        ncdf_varGet, cdfId, 'Y', Y_ghosted_GEI
        ncdf_varGet, cdfId, 'Z', Z_ghosted_GEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( ghostedFitFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_ghostedFit_GEI 
	    ncdf_varGet, cdfId, 'bT', bT_ghostedFit_GEI 
	    ncdf_varGet, cdfId, 'bP', bP_ghostedFit_GEI 

	    ncdf_varGet, cdfId, 'bX', bx_ghostedFit_GEI 
	    ncdf_varGet, cdfId, 'bY', by_ghostedFit_GEI 
	    ncdf_varGet, cdfId, 'bZ', bz_ghostedFit_GEI 

	    ncdf_varGet, cdfId, 'T', T_ghostedFit_GEI 
	    ncdf_varGet, cdfId, 'P', P_ghostedFit_GEI
	    ncdf_varGet, cdfId, 'R', R_ghostedFit_GEI

        ncdf_varGet, cdfId, 'X', X_ghostedFit_GEI
        ncdf_varGet, cdfId, 'Y', Y_ghostedFit_GEI
        ncdf_varGet, cdfId, 'Z', Z_ghostedFit_GEI

	ncdf_close, cdfId

 
 
    iiNorth = where ( T_raw_GEI*!radeg gt 0 and Z_raw_GEI gt 0 )
    dlg_plot_vecs, R_raw_GEI[iiNorth], T_raw_GEI[iiNorth], P_raw_GEI[iiNorth], $
        bT_raw_GEI[iiNorth], bP_raw_GEI[iiNorth]

    iiNorth = where ( Z_shifted_GEI gt 0 )
    dlg_plot_vecs, R_shifted_GEI[iiNorth], T_shifted_GEI[iiNorth], P_shifted_GEI[iiNorth], $
        bT_shifted_GEI[iiNorth], bP_shifted_GEI[iiNorth]

    iiNorth = where ( Z_ghosted_GEI gt 0 )
    dlg_plot_vecs, R_ghosted_GEI[iiNorth], T_ghosted_GEI[iiNorth], P_ghosted_GEI[iiNorth], $
        bT_ghosted_GEI[iiNorth], bP_ghosted_GEI[iiNorth]

    iiNorth = where ( Z_ghostedFit_GEI gt 0 )
    dlg_plot_vecs, R_ghostedFit_GEI[iiNorth], T_ghostedFit_GEI[iiNorth], P_ghostedFit_GEI[iiNorth], $
        bT_ghostedFit_GEI[iiNorth], bP_ghostedFit_GEI[iiNorth]

stop

end
