pro ampfit_plotdata

    iFName = 'data/20100824Amp_invert.ncdf'
    rawFName = 'output/ampData_original.nc'
    shiftedFName = 'output/ampData_shifted.nc'
	ghostedFName = 'output/ampData_ghosted.nc'
	ghostedFitFName = 'output/ampData_ghosted_fit.nc'
	UnShiftedFitFName = 'output/ampData_unshifted_fit.nc'
	gridAACGMFName = 'output/ampData_gridAACGM.nc'
	gridGEIFName = 'output/ampData_gridGEI.nc'
	gridShiftedGEIFName = 'output/ampData_gridShiftedGEI.nc'
	gridUnShiftedGEIFName = 'output/ampData_gridUnShiftedGEI.nc'

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

        ncdf_varGet, cdfId, 'sv_number', sv_ghosted_GEI

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

        ncdf_varGet, cdfId, 'sv_number', sv_ghostedFit_GEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( UnShiftedFitFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_UnShiftedFit_GEI 
	    ncdf_varGet, cdfId, 'bT', bT_UnShiftedFit_GEI 
	    ncdf_varGet, cdfId, 'bP', bP_UnShiftedFit_GEI 

	    ncdf_varGet, cdfId, 'bX', bx_UnShiftedFit_GEI 
	    ncdf_varGet, cdfId, 'bY', by_UnShiftedFit_GEI 
	    ncdf_varGet, cdfId, 'bZ', bz_UnShiftedFit_GEI 

	    ncdf_varGet, cdfId, 'T', T_UnShiftedFit_GEI 
	    ncdf_varGet, cdfId, 'P', P_UnShiftedFit_GEI
	    ncdf_varGet, cdfId, 'R', R_UnShiftedFit_GEI

        ncdf_varGet, cdfId, 'X', X_UnShiftedFit_GEI
        ncdf_varGet, cdfId, 'Y', Y_UnShiftedFit_GEI
        ncdf_varGet, cdfId, 'Z', Z_UnShiftedFit_GEI

        ncdf_varGet, cdfId, 'sv_number', sv_UnShiftedFit_GEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( gridAACGMFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_gridAACGM 
	    ncdf_varGet, cdfId, 'bT', bT_gridAACGM 
	    ncdf_varGet, cdfId, 'bP', bP_gridAACGM 

	    ncdf_varGet, cdfId, 'bX', bx_gridAACGM 
	    ncdf_varGet, cdfId, 'bY', by_gridAACGM 
	    ncdf_varGet, cdfId, 'bZ', bz_gridAACGM 

	    ncdf_varGet, cdfId, 'T', T_gridAACGM 
	    ncdf_varGet, cdfId, 'P', P_gridAACGM
	    ncdf_varGet, cdfId, 'R', R_gridAACGM

        ncdf_varGet, cdfId, 'X', X_gridAACGM
        ncdf_varGet, cdfId, 'Y', Y_gridAACGM
        ncdf_varGet, cdfId, 'Z', Z_gridAACGM

	ncdf_close, cdfId


    cdfId = ncdf_open ( gridGEIFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_gridGEI 
	    ncdf_varGet, cdfId, 'bT', bT_gridGEI 
	    ncdf_varGet, cdfId, 'bP', bP_gridGEI 

	    ncdf_varGet, cdfId, 'bX', bx_gridGEI 
	    ncdf_varGet, cdfId, 'bY', by_gridGEI 
	    ncdf_varGet, cdfId, 'bZ', bz_gridGEI 

	    ncdf_varGet, cdfId, 'T', T_gridGEI 
	    ncdf_varGet, cdfId, 'P', P_gridGEI
	    ncdf_varGet, cdfId, 'R', R_gridGEI

        ncdf_varGet, cdfId, 'X', X_gridGEI
        ncdf_varGet, cdfId, 'Y', Y_gridGEI
        ncdf_varGet, cdfId, 'Z', Z_gridGEI

	ncdf_close, cdfId


    cdfId = ncdf_open ( gridShiftedGEIFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_gridShiftedGEI 
	    ncdf_varGet, cdfId, 'bT', bT_gridShiftedGEI 
	    ncdf_varGet, cdfId, 'bP', bP_gridShiftedGEI 

	    ncdf_varGet, cdfId, 'bX', bx_gridShiftedGEI 
	    ncdf_varGet, cdfId, 'bY', by_gridShiftedGEI 
	    ncdf_varGet, cdfId, 'bZ', bz_gridShiftedGEI 

	    ncdf_varGet, cdfId, 'T', T_gridShiftedGEI 
	    ncdf_varGet, cdfId, 'P', P_gridShiftedGEI
	    ncdf_varGet, cdfId, 'R', R_gridShiftedGEI

        ncdf_varGet, cdfId, 'X', X_gridShiftedGEI
        ncdf_varGet, cdfId, 'Y', Y_gridShiftedGEI
        ncdf_varGet, cdfId, 'Z', Z_gridShiftedGEI

	ncdf_close, cdfId

    cdfId = ncdf_open ( gridUnShiftedGEIFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'bR', bR_gridUnShiftedGEI 
	    ncdf_varGet, cdfId, 'bT', bT_gridUnShiftedGEI 
	    ncdf_varGet, cdfId, 'bP', bP_gridUnShiftedGEI 

	    ncdf_varGet, cdfId, 'bX', bx_gridUnShiftedGEI 
	    ncdf_varGet, cdfId, 'bY', by_gridUnShiftedGEI 
	    ncdf_varGet, cdfId, 'bZ', bz_gridUnShiftedGEI 

	    ncdf_varGet, cdfId, 'T', T_gridUnShiftedGEI 
	    ncdf_varGet, cdfId, 'P', P_gridUnShiftedGEI
	    ncdf_varGet, cdfId, 'R', R_gridUnShiftedGEI

        ncdf_varGet, cdfId, 'X', X_gridUnShiftedGEI
        ncdf_varGet, cdfId, 'Y', Y_gridUnShiftedGEI
        ncdf_varGet, cdfId, 'Z', Z_gridUnShiftedGEI

	ncdf_close, cdfId


    iiNorth = where ( T_raw_GEI*!radeg gt 0 and Z_raw_GEI gt 0 )
    ampfit_plotvecs, R_raw_GEI[iiNorth], T_raw_GEI[iiNorth], P_raw_GEI[iiNorth], $
        bT_raw_GEI[iiNorth], bP_raw_GEI[iiNorth], title = 'dataGEI', fileNameIn='output/rawData.eps'

    iiNorth = where ( Z_shifted_GEI gt 0 )
    ampfit_plotvecs, R_shifted_GEI[iiNorth], T_shifted_GEI[iiNorth], P_shifted_GEI[iiNorth], $
        bT_shifted_GEI[iiNorth], bP_shifted_GEI[iiNorth], title = 'dataShiftedGEI', fileNameIn='output/shiftedData.eps'

    iiNorth = where ( Z_ghosted_GEI gt 0 )
    ampfit_plotvecs, R_ghosted_GEI[iiNorth], T_ghosted_GEI[iiNorth], P_ghosted_GEI[iiNorth], $
        bT_ghosted_GEI[iiNorth], bP_ghosted_GEI[iiNorth], title = 'dataGhostedShiftedGEI', fileNameIn='output/ghostedData.eps'

    iiNorth = where ( Z_ghostedFit_GEI gt 0 )
    ampfit_plotvecs, R_ghostedFit_GEI[iiNorth], T_ghostedFit_GEI[iiNorth], P_ghostedFit_GEI[iiNorth], $
        bT_ghostedFit_GEI[iiNorth], bP_ghostedFit_GEI[iiNorth], title = 'dataGhostedShiftedGEI[fit]', fileNameIn='output/ghostedFit.eps'

    iiNorth = where ( Z_UnShiftedFit_GEI gt 0 )
    ampfit_plotvecs, R_UnShiftedFit_GEI[iiNorth], T_UnShiftedFit_GEI[iiNorth], P_UnShiftedFit_GEI[iiNorth], $
        bT_UnShiftedFit_GEI[iiNorth], bP_UnShiftedFit_GEI[iiNorth], title = 'dataGhostedUnShiftedGEI[fit]', fileNameIn='output/UnShiftedFit.eps'

    iiNorth = where ( Z_gridAACGM gt 0 )
    ampfit_plotvecs, R_gridAACGM[iiNorth], T_gridAACGM[iiNorth], P_gridAACGM[iiNorth], $
        bT_gridAACGM[iiNorth], bP_gridAACGM[iiNorth], title = 'gridAACGM', fileNameIn='output/gridAACGM.eps'

    iiNorth = where ( Z_gridGEI gt 0 )
    ampfit_plotvecs, R_gridGEI[iiNorth], T_gridGEI[iiNorth], P_gridGEI[iiNorth], $
        bT_gridGEI[iiNorth], bP_gridGEI[iiNorth], title = 'gridGEI', fileNameIn='output/gridGEI.eps'

    iiNorth = where ( Z_gridShiftedGEI gt 0 )
    ampfit_plotvecs, R_gridShiftedGEI[iiNorth], T_gridShiftedGEI[iiNorth], P_gridShiftedGEI[iiNorth], $
        bT_gridShiftedGEI[iiNorth], bP_gridShiftedGEI[iiNorth], title = 'gridShiftedGEI', fileNameIn='output/gridShiftedGEI.eps'

    iiNorth = where ( Z_gridUnShiftedGEI gt 0 )
    ampfit_plotvecs, R_gridUnShiftedGEI[iiNorth], T_gridUnShiftedGEI[iiNorth], P_gridUnShiftedGEI[iiNorth], $
        bT_gridUnShiftedGEI[iiNorth], bP_gridUnShiftedGEI[iiNorth], title = 'gridUnShiftedGEI', fileNameIn='output/gridUnShiftedGEI.eps'
end
