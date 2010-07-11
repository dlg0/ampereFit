; +
; AmpereFit
; ---------
;
; Library to take Iridium vector (horizontal only) dB 
; field and create a FAC and dB map on regular grid.
;
; C.L. Waters and D.L. Green
; Dec 2009
;

pro amp_fit, sHr, eHr, south, $
	plot_bFns = plot_bFns, $
	path = path, $
	aacgmpath = aacgmpath, $
	savFileName = SavFileName, $
	kmax = kmax, mmax = mmax, $
	thresh = thresh, $
	sigma = sigma, $
	nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
	mn_fac = mn_fac, mx_fac=mx_fac, $
	plt_png = plt_png, $
	plt_tracks = plt_tracks, $
	aacgm_cap_coLat_deg = aacgm_cap_coLat_deg, $
	shiftGEI = shiftGEI, $
	debug = debug


	; Check for reasonable inputs
	; ---------------------------

    if size(sHr,/type) eq 0 $
            or size(eHr,/type) eq 0 $
            or size(south,/type) eq 0 then begin

        print, 'INPUT ERROR: need sHr, eHr and south'
        print, '    e.g., amp_fit, 4, 5, 0'
        stop

    endif


	; Set relative path to data
	; (do NOT go back to absolute paths)
	; ----------------------------------

	if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN64') then begin
		if not keyword_set(path) then path	= 'data\'
	endif else begin                        
		if not keyword_set(path) then path = 'data/'
	endelse


	; Set default keyword values
	; --------------------------	

	@amp_fit_defaults


	; Read input data for full hemisphere
	; Cap is selected out later based on 
	; AACGM / GEOG pole offsets
	; -----------------------------------

	Print,'Reading AMPERE Data File....'

    read_ampere_sav, sHr, eHr, south, 90.0,$            
        savFileName = savFileName, $
        dataOrig = dataOriginal, $                            
		dataShif = dataShifted, $
        year = year, month = month, day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        rot_mat=rot_mat, debug = debug;, /noShift


	; Create uniform grids in AACGM and GEI
	; -------------------------------------

 	aacgm_grid, aacgm_cap_coLat_deg, south, $
	   	 aacgmGrid_coLat_deg = aacgmGrid_coLat_deg, $
	   	 aacgmGrid_lon_deg = aacgmGrid_lon_deg, $
	   	 aacgmGrid_R_km = aacgmGrid_R_km, $
	     geiGrid_coLat_rad = geiGrid_coLat_rad, $
	   	 geiGrid_lon_rad = geiGrid_lon_rad , $
	   	 geiGrid_R_km = geiGrid_R_km, $
	   	 nLat = nLatGrid, nLon = nLonGrid, $
	   	 year = year, yrSec = avgYrSec, $
	   	 mltShift = mltShift, $
	   	 epoch = avgEpoch


	; Shift GEI grid to the shifted GEI system
	; THIS SHOULD BE A SUBROUTINE
	; ----------------------------------------

    geopack_sphcar, geiGrid_R_km[*], geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
            geiGridX, geiGridY, geiGridZ, $
			/to_rect       

	geiGridX_shifted = fltArr ( size(geiGridX,/dim) )
	geiGridY_shifted = fltArr ( size(geiGridY,/dim) )
	geiGridZ_shifted = fltArr ( size(geiGridZ,/dim) )

    for ii=long(0),n_elements(geiGridX)-1 do begin
			
		xyzGrid = [geiGridX[ii],geiGridY[ii],geiGridZ[ii]]
		xyzGridShifted = matrix_multiply ( rot_mat, xyzGrid, /Atrans )

		geiGridX_shifted[ii] = xyzGridShifted[0]
		geiGridY_shifted[ii] = xyzGridShifted[1]
		geiGridZ_shifted[ii] = xyzGridShifted[2]

	endfor

    geopack_sphcar, geiGridX_shifted, geiGridY_shifted, geiGridZ_shifted,$
			geiGrid_rIrid_km_shifted, geiGrid_coLat_rad_shifted, geiGrid_lon_rad_shifted,$
			/to_sphere  


	; Choose shifted or unshifted data
	; --------------------------------

	if shiftGEI then begin
			data = dataShifted
	endif else begin
			data = dataOriginal
	endelse


	; Use the maximum coLat of the regular grid in the shifted GEI 
	; system as the capSize for the fit and hence also for the 
	; selection of how much data to use in the fit. So here we select
	; out data to cover the requested regular grid.
	; ---------------------------------------------------------------

	if south eq 1 then begin
		minCoLat_GEI_deg_shifted = min ( geiGrid_coLat_rad_shifted * !radeg ) - 2
		maxCoLat_GEI_deg_shifted = 180.0
	endif else begin
		minCoLat_GEI_deg_shifted = 0.0
		maxCoLat_GEI_deg_shifted = max ( geiGrid_coLat_rad_shifted * !radeg ) + 2
	endelse


	iiShiftedCap = where ( data.gei_coLat_rad * !radeg ge minCoLat_GEI_deg_shifted $
					and data.gei_coLat_rad*!radeg le maxCoLat_GEI_deg_shifted, $
					iiShiftedCapCnt )

	data = data[iiShiftedCap]


	;; Plot the AACGM grid in GEI, shifted and not
	;; -------------------------------------------

	;winNo = 0
	;window, winNo 
	;set_map, maxCoLat_GEI_deg_shifted, title = 'AACGM Grid (GEI)'
	;plots, geiGrid_lon_rad*!radeg, 90-geiGrid_coLat_rad*!radeg, psym = 4
	;winNo++
	;window, winNo
	;set_map, maxCoLat_GEI_deg_shifted, title = 'AACGM Grid Shifted (GEI)'
	;plots, geiGrid_lon_rad_shifted*!radeg, 90-geiGrid_coLat_rad_shifted*!radeg, psym = 4


	; Create basis functions at shifted GEI coords
	; --------------------------------------------

	Print,'Generating Basis Set at Input Data Locations....'

	; For Sth, 90->180 deg
 	ampere_setupSHFns, 1.0, data.gei_coLat_rad, data.gei_lon_rad, $   
			kMax, mMax, $
			minTheta = minCoLat_GEI_deg_shifted, $
			maxTheta = maxCoLat_GEI_deg_shifted, $
			bThBFnArr = dYkmDThBFns, $
			bPhBFnArr = dYkmDPhBFns, $
			YBFnArr = YkmBFns, $
			mArr = outMValues, $
			nArr = outNkValues, $
			kArr = outKValues, /bc2


	; Fit |dB| to dB.grad Ykm in the shifted GEI system
	; -------------------------------------------------

	Print,'Fitting dB Data to Basis Set....'

	dBMag   = sqrt(data.bTheta_GEI^2+data.bPhi_GEI^2)

	; Apply weighting [after Dec 10 discussion at APL]
	; Please include a description here of said discussion.
	; Since right now I do not see the point. Also, I think
	; the weighting should be applied to the LHS of the linear
	; system to be solved.

    if (sigma ne 1) then begin

     	w_idx=where(dBMag lt thresh)

     	If w_idx(0) gt -1 then begin

     	 	data[w_idx].bTheta_GEI=data[w_idx].dTheta_GEI/sigma         ; CHECK THIS
     	 	data[w_idx].bPhi_GEI=data[w_idx].dPhi_GEI/sigma

	 	endif
	endif

    bFuncs  = temporary ( $
                [[ dYkmDThBfns, -dYkmDPhBfns ], $
                 [ dYkmDPhBfns,  dYkmDthBfns ]])

    print, 'bFuncs is', n_elements(bFuncs[*])*16/(1024.0^2),' MB'


	; Fit to shifted dBs
	; ------------------

    coeffs_ = la_least_squares ( bFuncs, $
        [data.bTheta_GEI, data.bPhi_GEI], $
        status = stat, $
        method = 2, $
        /double, rank = rank )

	if stat ne 0 then begin

		print, 'ERROR: la_least_squares', status
		stop

	endif

	Print,'Calculating FAC at Input Data Locations....'

    fit = bFuncs ## coeffs_


	; Extract theta and phi component of fit from solution
	; ----------------------------------------------------

	fit_bTheta_GEI = fit[0:n_elements(dBMag)-1]  
	fit_bPhi_GEI   = fit[n_elements(dBMag):*] 


	; Calc rms error for these
	; ------------------------

	err_dbTheta	= total ( ( data.bTheta_GEI-fit_bTheta_GEI)^2 )
	err_dbTheta	= sqrt ( err_dbTheta/n_elements(fit_bTheta_GEI))
	err_dbPhi	= total ( ( data.bPhi_GEI-fit_bPhi_GEI)^2 )
	err_dbPhi	= sqrt ( err_dbPhi/n_elements(fit_bPhi_GEI))

	Print,'RMS Error : dbTheta = ',err_dbTheta
	Print,'RMS Error : dbPhi = ',err_dbPhi


	; Plot along track comparison of the raw and fitted data
	; ------------------------------------------------------

	if debug then begin
		@load_colors
		iPlot, view_grid = [2,3]
		for i=0,5 do begin
			iiTrack = where ( data.iPln eq i, iiTrackCnt )
			iPlot, (data.bPhi_GEI)[iiTrack], $
					sym_index = 4, $
					view_number = i+1, $
					/stretch_to_fit, $
					lineStyle = 6
			iPlot, fit_bPhi_GEI[iiTrack], $
					thick = 2, /over, $
					color = red
		endfor
	endif


	; Construct jPar in the GEI system
	; --------------------------------

	@constants
	jPar  = (YkmBFns ## $
	        (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}


	; Generate basis fns at regular grid 
	; ----------------------------------

	print,'Generating Basis Set over Uniform Grid....'

   	ampere_setupSHFns, 1.0, geiGrid_coLat_rad_shifted[*], geiGrid_lon_rad_shifted[*], $
			kMax, mMax, $
			minTheta = minCoLat_GEI_deg_shifted, $
			maxTheta = maxCoLat_GEI_deg_shifted, $
			bThBFnArr = dYkmDThBFns_grid, $
			bPhBFnArr = dYkmDPhBFns_grid, $
			YBFnArr = YkmBFns_grid, $
			mArr = outMValues_grid, $
			nArr = outNkValues_grid, $
			kArr = outKValues_grid, /bc2

	Print,'Calculating FAC over Uniform Grid....'
	
    bFuncs_grid  = temporary($
                    [[dYkmDThBfns_grid, -dYkmDPhBfns_grid ], $
                     [dYkmDPhBfns_grid,  dYkmDthBfns_grid ]])

  	jParAACGM	= (YkmBFns_grid ## $
         (coeffs_[n_elements(YkmBFns_grid[*,0]):*]*(-outNkValues_grid*(outNkValues_grid+1.0))))$
                         /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}

	jParAACGM	= reform(jParAACGM, nLatGrid, nLonGrid)


	; Plot jPar in various coord systems
	; ----------------------------------

	if debug then begin

		pole = 90
		capLimit = maxCoLat_GEI_deg_shifted
		if(south) then begin
			pole = -90
			capLimit = minCoLat_GEI_deg_shifted
		endif

		winNo = -1
		winNo++
		window, winNo, xSize = 600, ySize = 600
		!p.multi = [0,2,2]	
		!p.charSize = 1.0
	
		plot_fac, jPar, pole, capLimit, $
				data.gei_coLat_deg, data.gei_lon_deg, $
				title = 'jPar [GEI]', $
				south = south

		plot_fac, jParAACGM[*], pole, capLimit, $
				geiGRid_coLat_rad_shifted[*]*!radeg, geiGrid_lon_rad_shifted[*]*!radeg, $
				title = 'jPar on grid [GEI]', $
				south = south

		plot_fac, jParAACGM[*], pole, capLimit, $
				aacgmGrid_coLat_deg[*], aacgmGrid_lon_deg[*], $
				title = 'jPar on grid [AACGM]', $
				south = south

		plot_fac, jParAACGM[*], pole, capLimit, $
				aacgmGrid_coLat_deg[*], aacgmGrid_lon_deg[*]+mltShift*15, $
				title = 'jPar on grid [AACGM-MLT]', $
				south = south

	endif


	; Reconstruct the dB vector from the fit
	; --------------------------------------

    recon_dB_geiGrid	= bFuncs_grid ## coeffs_
   	dBTheta_geiGrid		= recon_dB_geiGrid[0:nLatGrid*nLonGrid-1]
 	dBPhi_geiGrid		= recon_dB_geiGrid[nLatGrid*nLonGrid:*]
	dBR_geiGrid			= dBPhi_geiGrid * 0


	; Un-Shift the dB Vectors
	; THIS SHOULD ALSO BE A SUBROUTINE
	; --------------------------------

	if shiftGEI then begin

		rev_rot_mat = transpose(rot_mat)       
		
		; conv dB to XYZ comp

    	geopack_bspcar,geiGrid_coLat_rad, geiGrid_lon_rad, $     
    	       dBR_geiGrid, dBTheta_geiGrid, dBPhi_geiGrid ,$
    	       vx_a,vy_a,vz_a

   ;	 AAARRRGGGHH!!!! USE USEFUL VARIABLE NAMES!!!

		for ii=Long(0),n_elements(dBTheta_geiGrid_sh)-1 do begin

    		dBx=vx_a[ii]*rev_rot_mat[0,0] + $
    		    vy_a[ii]*rev_rot_mat[1,0] + $
    		    vz_a[ii]*rev_rot_mat[2,0]
    		dBy=vx_a[ii]*rev_rot_mat[0,1] + $
    		    vy_a[ii]*rev_rot_mat[1,1] + $
    		    vz_a[ii]*rev_rot_mat[2,1]
    		dBz=vx_a[ii]*rev_rot_mat[0,2] + $
    		    vy_a[ii]*rev_rot_mat[1,2] + $
    		    vz_a[ii]*rev_rot_mat[2,2]

			; convert to spherical

    		geopack_bcarsp, geiGridX[ii], geiGridY[ii], geiGridZ[ii], $
    		      dbx, dBy, dBz, vbr,vbth,vbph

    		dBTheta_geiGrid[ii]	= vbth
    		dBPhi_geiGrid[ii]		= vbph

		endfor

	endif


	; Rotate fitted, raw and shifted (or not) dB vectors to AACGM
	; -----------------------------------------------------------

	geiVec_to_aacgmVec, $
		geiGrid_R_km, geiGrid_coLat_rad, geiGrid_lon_rad, $
		dBTheta_geiGrid, dbPhi_geiGrid, $
		dBTheta_aacgm_gth = dBTheta_aacgmGrid_gth, $
		dBPhi_aacgm_gph = dBPhi_aacgmGrid_gph

	geiVec_to_aacgmVec, $
		data.gei_R_km, data.gei_coLat_rad, data.gei_lon_rad, $
		data.bTheta_gei, data.bPhi_gei, $
		dBTheta_aacgm_gth = dBTheta_aacgmData_gth, $
		dBPhi_aacgm_gph = dBPhi_aacgmData_gph, $
		aacgm_lat_deg = aacgmData_lat_deg, $
		aacgm_lon_deg = aacgmData_lon_deg

	geiVec_to_aacgmVec, $
		dataOriginal.gei_R_km, dataOriginal.gei_coLat_rad, dataOriginal.gei_lon_rad, $
		dataOriginal.bTheta_gei, dataOriginal.bPhi_gei, $
		dBTheta_aacgm_gth = dBTheta_aacgmDataOriginal_gth, $
		dBPhi_aacgm_gph = dBPhi_aacgmDataOriginal_gph, $
		aacgm_lat_deg = aacgmDataOriginal_lat_deg, $
		aacgm_lon_deg = aacgmDataOriginal_lon_deg


	; Plot the dB vectors in AACGM-MLT
	; --------------------------------

	if debug then begin

		winNo ++
		window, winNo, xSize = 900, ySize = 300
		!p.multi = [0,3,1]
		!p.charSize = 2.0

		plt_dat, 90.0-aacgmDataOriginal_lat_deg[*], ((aacgmDataOriginal_lon_deg[*]/15.0+mltShift) mod 24)*15.0, $
			n_elements(aacgmDataOriginal_lat_deg[*] ), $
			-dBTheta_aacgmDataOriginal_gth, dBPhi_aacgmDataOriginal_gph, $
			south, title = 'Input dataOriginal: AACGM-MLT', $
			capSize = abs(aacgm_cap_coLat_deg)


		plt_dat, 90.0-aacgmData_lat_deg[*], ((aacgmData_lon_deg[*]/15.0+mltShift) mod 24)*15.0, $
		         n_elements(aacgmData_lat_deg[*] ), -dBTheta_aacgmData_gth, dBPhi_aacgmData_gph, south,$
				 title = 'Input data: AACGM-MLT', $
    	        capSize = abs(aacgm_cap_coLat_deg)

		; ***** Check if it should be 90.0-abs(m_lat_a)	*****
		plt_dat, aacgmGrid_coLat_deg[*], ((aacgmGrid_lon_deg[*]/15.0+mltShift) mod 24)*15.0, $
		        n_elements(aacgmGrid_coLat_deg[*] ), -dBTheta_aacgmGrid_gth, dBPhi_aacgmGrid_gph, south,$
			 	title = 'Fitted data: AACGM-MLT', $
 				capSize = abs(aacgm_cap_coLat_deg)

		!p.charSize = 1.0
		!p.multi = 0


		; Write PNG image
		; ---------------

		if plt_png eq 1 then begin

		    png_file=path+'amp_dB&FAC_'+yr_str+mnth_str+dy_str+' at '+hr_str+mn_str+sc_str+'_'+hem+'.png'
    		image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
    		write_PNG,png_file,image,r,g,b
    		Print,'PNG for FACs written to ', png_file

    	end

	endif

	;@big_nsf_plot

end

