; AMPERE dB to FAC code
; Given a set of Iridium dB data (in GEI coords), compute the associated dB vectors
; in AACGM and the FAC pattern using SCHA
;
; C.L. Waters and D.L. Green
; Dec 2009
;


; Form string expression from time
pro tmestr1,Hr,Mn,Sc,hr_str,mn_str,sc_str

 	hr_str=strtrim(string(Hr),2)
 	if Hr lt 10 then hr_str='0'+hr_str
 	mn_str=strtrim(string(Mn),2)
 	if Mn lt 10 then mn_str='0'+mn_str
 	sc_str=strtrim(string(fix(sc)),2)
 	if Sc lt 10 then sc_str='0'+sc_str

end

pro cross_p,v1,v2,vout

 	v1=reform(v1) & v2=reform(v2)
 	vout=[[v1(1)*v2(2)-v1(2)*v2(1)],$
 	      [v1(2)*v2(0)-v1(0)*v2(2)],$
 	      [v1(0)*v2(1)-v1(1)*v2(0)]]
 	vout=transpose(vout)

end

pro amp_fit, sHr, eHr, south, $
	plot_bFns = plot_bFns, $
	path = path, $
	aacgmpath = aacgmpath, $
	pnmPath = pnmpath, $
	savFileName = SavFileName, $
	kmax = kmax, mmax = mmax, $
	thresh = thresh, $
	sigma = sigma, $
	nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
	mn_fac = mn_fac, mx_fac=mx_fac, $
	plt_png = plt_png, $
	plt_tracks = plt_tracks, $
	aacgm_cap_coLat_deg = aacgm_cap_coLat_deg


	; Check for reasonable inputs
	; ---------------------------

    if size(sHr,/type) eq 0 $
            or size(eHr,/type) eq 0 $
            or size(south,/type) eq 0 then begin

        print, 'INPUT ERROR: need sHr, eHr and south'
        print, '    e.g., amp_fit, 4, 5, 0'
        stop

    endif


	; Set path to data files and pnmSav files
	; ---------------------------------------

	if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN64') then begin
		if not keyword_set(path) then path	= 'd:\cwac\hi_res\davidg\'
		if not keyword_set(pnmpath) then pnmPath = path + 'jpar_ver2\pnmsavs\pnmSav'
	endif else begin                        
		if not keyword_set(path) then path = '~/code/ampereFit/idl/'
		if not keyword_set(pnmpath) then pnmPath = path + 'pnmSavs/pnmSav'
	endelse


	; Set default keyword values
	; --------------------------	

	@amp_fit_defaults


	; Read input data
	; ---------------

	Print,'Reading AMPERE Data File....'

    read_ampere_sav, sHr, eHr, south, 90.0,$            
        savFileName = savFileName, $
        dataOrig = dataOriginal, $                            
		dataShif = dataShifted, $
        year = year, month = month, day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        rot_mat=rot_mat, /show, /noShift


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


	; Shift the regular gei grid into the shifted gei system
	; ------------------------------------------------------

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
			geiGrid_R_km_shifted, geiGrid_coLat_rad_shifted, geiGrid_lon_rad_shifted,$
			/to_sphere  


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

	iiShiftedCap = where ( dataShifted.gei_coLat_rad * !radeg ge minCoLat_GEI_deg_shifted $
					and dataShifted.gei_coLat_rad*!radeg le maxCoLat_GEI_deg_shifted, $
					iiShiftedCapCnt )

	data = dataShifted[iiShiftedCap]


	; Plot the AACGM grid in GEI, shifted and not
	; -------------------------------------------

	window, 0
	set_map, maxCoLat_GEI_deg_shifted, title = 'AACGM Grid (GEI)'
	plots, geiGrid_lon_rad*!radeg, 90-geiGrid_coLat_rad*!radeg, psym = 4
	window, 1
	set_map, maxCoLat_GEI_deg_shifted, title = 'AACGM Grid Shifted (GEI)'
	plots, geiGrid_lon_rad_shifted*!radeg, 90-geiGrid_coLat_rad_shifted*!radeg, psym = 4


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


	; Compare the lookup table and 'on-the-fly' basis fns
	; ---------------------------------------------------

    if keyword_set ( plot_bFns ) then begin            

	    !p.multi = [0,3,4]
	    !p.charSize = 2.0
	    device, decomposed = 0
	    window, 0, xSize = 1200, ySize = 800
	    !p.background = 255

        m1  = 0
        k1  = 6
	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
	    		title = 'table lookup (m=2,k=5, Y)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdTh', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdPh', $
	    		psym=4,color=0

        m2  = 0
        k2  = 8
	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
	    		title = 'table lookup (m=0,k=7)', $
	    		psym=4,color=0

	    !p.multi = 0.0
	    !p.charSize = 1.0

    endif


	; Fit |dB| to dB.grad Ykm in the shifted GEI system
	; -------------------------------------------------

	Print,'Fitting dB Data to Basis Set....'

	dBMag   = sqrt(data.bTheta_GEI^2+data.bPhi_GEI^2)

	; Apply weighting [after Dec 10 discussion at APL]
	; Please include a description here of said discussion.
	; Since right now I do not see the point. Also, I think
	; the weighting should be applied to the LHS of the linear
	; system to be solved.

    if (sigma ne 1.) then begin

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


	; Plot along track comparison of the raw and fitted data
	; ------------------------------------------------------

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


	@constants
	jPar  = (YkmBFns ## $
	        (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}

	pole = 90
	capLimit = maxCoLat_GEI_deg_shifted
	if(south) then begin
		pole = -90
		capLimit = minCoLat_GEI_deg_shifted
	endif
	plot_fac, jPar, pole, capLimit, $
			data.gei_coLat_deg, data.gei_lon_deg, $
			title = 'jPar GEI', $
			south = south
stop
	; jPar is an array[1,Num of dB values]

	; Calc rms error for these

	err_dbTheta=total( (data.bTheta_GEI-fit_bTheta_GEI)^2 )
	err_dbTheta=sqrt(err_dbTheta/n_elements(fit_bTheta_GEI))
	err_dbPhi=total( (data.bPhi_GEI-fit_bPhi_GEI)^2 )
	err_dbPhi=sqrt(err_dbPhi/n_elements(fit_bPhi_GEI))
	Print,'RMS Error : dbTheta = ',err_dbTheta
	Print,'RMS Error : dbPhi = ',err_dbPhi

	Print,'Generating Basis Set over Uniform Grid....'

	; generate basis fns at regular grid (shifted)
	; ----------------------------------

   	ampere_setupSHFns, 1.0, geiGrid_coLat_rad_shifted[*], geiGrid_lon_rad_shifted[*], $
			kMax, mMax, $
			minTheta = minTheta, $
			maxTheta = maxTheta, $
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

	; FAC array	of [nLatGrid, nLonGrid]

	jParAACGM	= reform(jParAACGM, nLatGrid, nLonGrid)

	; Still the shifted grid

    fit_grid_sh    = bFuncs_grid ## coeffs_
   	dBTheta_GEI_grid_sh = fit_grid_sh[0:nLatGrid*nLonGrid-1]
 	dBPhi_GEI_grid_sh   = fit_grid_sh[nLatGrid*nLonGrid:*]
	dBTheta_GEI_grid=dBTheta_GEI_grid_sh
	dBPhi_GEI_grid=dBPhi_GEI_grid_sh
	dBR_GEI_grid_sh=dBPhi_GEI_grid_sh*0.0

	Print,'Un-Shift the dB Vectors....'

	; Need to rotate the dBs back
	; undo rotation by transpose

	rev_rot_mat=transpose(rot_mat)       
	
	; conv dB to XYZ comp

    geopack_bspcar,geiGrid_coLat_rad_shifted, geiGrid_lon_rad_shifted, $     
           dBR_GEI_grid_sh, dBTheta_GEI_grid_sh, dBPhi_GEI_grid_sh ,$
           vx_a,vy_a,vz_a

	for ii=Long(0),n_elements(dBTheta_GEI_grid_sh)-1 do begin

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
    	dBTheta_GEI_grid[ii]=vbth
    	dBPhi_GEI_grid[ii]=vbph
	endfor

	; Take GEI -> GEOG -> AACGM

    Re_km = 6371.0
    R_km = 780.0 + Re_km

	; conv GEI r,thet,phi -> GEI XYZ

    geopack_sphcar,geiGrid_R_km, geiGrid_coLat_rad, geiGrid_lon_rad,$
           gei_x, gei_y, gei_z, /to_rect

	; conv GEI XYZ -> GEOG XYZ

    geoPack_conv_coord, gei_x, gei_y, gei_z, $
                        geo_x, geo_y, geo_z, $
                       /from_gei, /to_geo

    epoch = fltArr(n_elements(gei_x[*])) + Avgepoch

	; conv GEOG XYZ -> GEOG r,thet,phi

    geopack_sphcar, geo_x, geo_y, geo_z, $
                    geog_R_km, geog_coLat_rad, geog_lon_rad, /to_sphere

    gLat_a=90.0-geog_coLat_rad*!radeg
    gLon_a=geog_lon_rad*!radeg
    aacgm_yr=2005
    aacgm_conv_vec,gLat_a, gLon_a, geog_R_km-Re_km, $
                   dBTheta_GEI_grid, dBPhi_GEI_grid, $
                   mlat_a, mlon_a, $
                   mth_vec_gth, mph_vec_gth, mth_vec_gph, mph_vec_gph, err,/to_aacgm

    ;vec_geo2aacgm,path,aacgm_yr,geog_R_km-Re_km,$
    ;                   gLat_a,gLon_a,dBTheta_GEI_grid,dBPhi_GEI_grid,$
    ;                   mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph


	; Input data
	; Take GEI -> GEOG -> AACGM

    Re_km = 6371.0
    R_km = 780.0 + Re_km

	; conv GEI r,thet,phi -> GEI XYZ

    geopack_sphcar,data.gei_R_km, data.gei_coLat_rad, data.gei_lon_rad,$
           gei_x_in, gei_y_in, gei_z_in, /to_rect

	; conv GEI XYZ -> GEOG XYZ

    geoPack_conv_coord, gei_x_in, gei_y_in, gei_z_in, $
                        geo_x_in, geo_y_in, geo_z_in, $
                       /from_gei, /to_geo

    epoch = fltArr(n_elements(gei_x_in[*])) + Avgepoch

	; conv GEOG XYZ -> GEOG r,thet,phi

    geopack_sphcar, geo_x_in, geo_y_in, geo_z_in,$
                    geog_R_km_in, geog_coLat_rad_in, geog_lon_rad_in, /to_sphere

    gLat_a_in=90.0-geog_coLat_rad_in*!radeg
    gLon_a_in=geog_lon_rad_in*!radeg

    aacgm_conv_vec,gLat_a_in, gLon_a_in, geog_R_km_in-Re_km, $
                   data.bTheta_GEI, data.bPhi_GEI, $
                   mlat_a_in, mlon_a_in, $
                   mth_vec_gth_in, mph_vec_gth_in, mth_vec_gph_in, mph_vec_gph_in, err,/to_aacgm

    ;vec_geo2aacgm,path,aacgm_yr,geog_R_km_in-Re_km,$
    ;                   gLat_a_in,gLon_a_in,data.dBTheta,data.dBPhi,$
    ;                   mlat_a_in,mlon_a_in,$
    ;                   mth_vec_gth_in,mph_vec_gth_in,mth_vec_gph_in,mph_vec_gph_in



	; plot stuff
	; ----------


	!p.background = 255

    yr_str=strtrim(string(fix(year)),2)
    mnth_str=strtrim(string(fix(month)),2)
    If month lt 10 then mnth_str='0'+mnth_str
    dy_str=strtrim(string(fix(day)),2)
    If day lt 10 then dy_str='0'+dy_str
	date_str='dd/mm/yyyy= '+dy_str+'/'+mnth_str+'/'+dy_str

    StHr=fix(sHr)
    StMn=fix((SHr*3600.0-float(StHr)*3600.0)/60.0)
    StSc=sHr*3600.0-float(StHr)*3600.0-float(StMn)*60.0
	tmestr1,StHr,StMn,StSc,hr_str,mn_str,sc_str
	tme_str='hh:mm:ss= '+hr_str+':'+mn_str+':'+sc_str

	Print,'Plotting dB and FAC....'
	wait,0.001

	; plot the db vectors
	; -------------------
	; plot the FAC maps

    hem='north'
    if south then hem='south'
  	fac_div=mx_fac/10.
	pole=90
	sgn=1.

	if south then begin
	 	pole=-90
	 	aacgm_cap_coLat_deg=-aacgm_cap_coLat_deg
	 	sgn=-1.
	endif

	winNum = 0

	!p.backGround = 255
	device, decomposed = 0
	!p.font = 0
	!p.multi = [0,3,1,0]
	jLevels = (fIndGen(21)-10.)*fac_div
	colors  = bytScl ( jLevels, top = 253 ) + 1
	Loadct,0,/silent
	window, winnum, xSize = 1050, ySize = 400,title='dB and FAC for '+date_str+'  '+tme_str
	winnum++

	plt_dat, 90.0-mlat_a_in[*], ((mlon_a_in[*]/15.0+mltShift) mod 24)*15.0, $
	         n_elements(mlat_a_in[*] ), -mth_vec_gth_in, mph_vec_gph_in, [1,1], [1,2], south,$
			 title = 'Input data: AACGM-MLT', $
            capSize = abs(aacgm_cap_coLat_deg)

	; ***** Check if it should be 90.0-abs(m_lat_a)	*****
	plt_dat, 90.0-mlat_a[*], ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	         n_elements(mlat_a[*] ), -mth_vec_gth, mph_vec_gph, [1,1], [1,2], south,$
			 title = 'Fitted data: AACGM-MLT', $
            capSize = abs(aacgm_cap_coLat_deg)

	If plt_png eq 1 then begin
	    png_file=path+'amp_dB&FAC_'+yr_str+mnth_str+dy_str+' at '+hr_str+mn_str+sc_str+'_'+hem+'.png'
    	image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
    	Write_PNG,png_file,image,r,g,b
    	Print,'PNG for FACs written to ', png_file
    end

	;@big_nsf_plot

	Print,'Min:Max FAC = ',min(jParTmp),'  ',max(jParTmp)
	!p.multi = 0
end

