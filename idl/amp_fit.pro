; AMPERE dB to FAC code
; Given a set of Iridium dB data (in GEI coords), compute the associated dB vectors
; in AACGM and the FAC pattern using SCHA
;
; C.L. Waters and D.L. Green
; Dec 2009
;
Pro tmestr1,Hr,Mn,Sc,hr_str,mn_str,sc_str
; Form string expression from time
 hr_str=strtrim(string(Hr),2)
 If Hr lt 10 then hr_str='0'+hr_str
 mn_str=strtrim(string(Mn),2)
 If Mn lt 10 then mn_str='0'+mn_str
 sc_str=strtrim(string(fix(sc)),2)
 If Sc lt 10 then sc_str='0'+sc_str
end
;
; -----------------------------------------------------------------------------------
;
Pro cross_p,v1,v2,vout
 v1=reform(v1) & v2=reform(v2)
 vout=[[v1(1)*v2(2)-v1(2)*v2(1)],$
       [v1(2)*v2(0)-v1(0)*v2(2)],$
       [v1(0)*v2(1)-v1(1)*v2(0)]]
 vout=transpose(vout)
end
;
; -----------------------------------------------------------------------------------
;
pro amp_fit, sHr, eHr, south, $
        plot_bFns = plot_bFns, $
		path = path, $
		aacgmpath = aacgmpath, $
;		pnmPath = pnmpath, $
		savFileName = SavFileName, $
		kmax = kmax, mmax = mmax, $
		thresh = thresh, $
		sigma = sigma, $
		nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
		mn_fac = mn_fac, mx_fac=mx_fac, $
		plt_png = plt_png, $
		plt_tracks = plt_tracks, $
		plt_coLat = plt_coLat

	if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN62') then begin
		if not keyword_set(path) then path	= 'd:\cwac\hi_res\davidg\'
;		if not keyword_set(pnmpath) then pnmPath = path + 'jpar_ver2\pnmsavs\pnmSav'
	endif else begin                        ; other systems (linux, darwin etc.)
		if not keyword_set(path) then path = '~/code/ampereFit/idl/'
;		if not keyword_set(pnmpath) then pnmPath = path + 'pnmSavs/pnmSav'
	endelse
	if keyword_set(aacgmpath) then aacgm_set_path,aacgmpath else aacgm_set_path,path

; Input data file
	if not keyword_set(savFileName) then savFileName = path+'20091125Amp_invert.sav'

; Basis functions
	if not keyword_set(kmax) then kmax = 35
	if not keyword_set(mmax) then mmax = 5           ; 0 to 5

; Data weighting parameters
	if not keyword_set(thresh) then thresh=80.0      ; dBMag threshold, below this we apply 1/sigma: no weighting if thresh=0
	if not keyword_set(sigma) then sigma=1.0         ; sigma=2 would halve the amplitude of each component of the data

; Controls final solution grid
	if not keyword_set(nLatGrid) then nLatGrid	= 50 ; default fit grid
	if not keyword_set(nLonGrid) then nLonGrid	= 24

; Min and Max FAC for plot
	if not keyword_set(mn_fac) then mn_fac=0.07      ; do not plot abs(FAC) < mn_fac
	if not keyword_set(mx_fac) then mx_fac=0.50

; Switch to create PNG files of dB and FAC
	if not keyword_set(plt_png) then plt_png=0       ; default is to not output PNG files

; Switch to plot fitted dB by Iridium orbit track
	if not keyword_set(plt_tracks) then plt_tracks=0       ; default is to not to plot track data

; Plot range in Latitude
	if not keyword_set(plt_coLat) then plt_coLat = 50.0

	winNum	= 1

; Restore saveset of input data
; Calc Lat,Lon shift to place avg orbit intersection point at centre
	Print,'Reading AMPERE Data File....'
	wait,0.001
; Need to read data first to get shift -> read full hemisphere of data
    read_ampere_sav, sHr, eHr, south, 90.0,$            ; includes amp_shiftdata routine in here
        savFileName = savFileName, $
        dataOut = data_rd, $                            ; GEI are all shifted in structure data
        year = year, month = month, day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        rot_mat=rot_mat
; Output is in variable 'data'
; Calc uniform grid (output grid) in order to determine co_Lat limits of data to read
; 1. Calculate Uniform AACGM grid
; 2. Conv to geog -> then to GEI
 	aacgm_grid, plt_coLat, south, $
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
    mltShift    = mltShift[0]

; 3. Calc shifted GEI from uniform AACGM grid -> get max coLat
; Calc XYZ components of uniform unshifted GEI grid
    geopack_sphcar, geiGrid_R_km[*], geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
            gei_xGrid, gei_yGrid, gei_zGrid,/to_rect       ; XYZ coord of GEI grid

; create room for shifted [Lat,Lon]
	geiGrid_coLat_rad_sh=geiGrid_coLat_rad[*]
	geiGrid_lon_rad_sh=geiGrid_lon_rad[*]

; r and r_sh are the same as its a rotation about the same origin
    For ii=Long(0),n_elements(gei_xGrid)-1 do begin
; rotate coord locations
     xsh=gei_xGrid[ii]*rot_mat[0,0] + $
         gei_yGrid[ii]*rot_mat[1,0] + $
         gei_zGrid[ii]*rot_mat[2,0]
     ysh=gei_xGrid[ii]*rot_mat[0,1] + $
         gei_yGrid[ii]*rot_mat[1,1] + $
         gei_zGrid[ii]*rot_mat[2,1]
     zsh=gei_xGrid[ii]*rot_mat[0,2] + $
         gei_yGrid[ii]*rot_mat[1,2] + $
         gei_zGrid[ii]*rot_mat[2,2]
     geopack_sphcar,xsh,ysh,zsh,r_sh,th_sh,ph_sh,/to_sphere  ; r,thet,phi coord of GEI grid (shifted)
     geiGrid_coLat_rad_sh[ii]=th_sh
   	 geiGrid_lon_rad_sh[ii]=ph_sh
	end
; Get limits of shifted GEI grid
	mx_gei_coLat_sh=max(geiGrid_coLat_rad_sh*180.0/!dpi)

; 4. Prune data back to this co_Lat
    minTheta=0.0                   ; for Nth hemis
    maxTheta=round(mx_gei_coLat_sh+2.)
    if south eq 1 then begin
     minTheta=180.0-mx_gei_coLat_sh
     maxTheta=180.0
    end

    data_struct={ipln:0, isat:0, utc:0d0, $
              gei_R_km: 0d0, gei_coLat_rad: 0d0, gei_lon_rad: 0d0, $
              dbR: 0d0, dbTheta: 0d0, dbPhi: 0d0, $
              gei_R_km_sh : 0d0, gei_coLat_rad_sh: 0d0, gei_lon_rad_sh: 0d0, $
              dbR_sh: 0d0, dbTheta_sh: 0d0, dbPhi_sh: 0d0 }

	if south then idx_dat=where(data_rd.gei_coLat_rad_sh*180.0/!dpi ge minTheta) else $
	              idx_dat=where(data_rd.gei_coLat_rad_sh*180.0/!dpi le maxTheta)

    data=replicate(data_struct, n_elements(idx_dat))
	data.ipln             = data_rd.ipln[idx_dat]
	data.isat             = data_rd.isat[idx_dat]
	data.utc              = data_rd.utc[idx_dat]
	data.gei_R_km         = data_rd.gei_R_km[idx_dat]
	data.gei_coLat_rad    = data_rd.gei_coLat_rad[idx_dat]
	data.gei_lon_rad      = data_rd.gei_lon_rad[idx_dat]
	data.dbR              = data_rd.dbR[idx_dat]
	data.dbTheta          = data_rd.dbTheta[idx_dat]
	data.dbPhi            = data_rd.dbPhi[idx_dat]
	data.gei_R_km_sh      = data_rd.gei_R_km_sh[idx_dat]
	data.gei_coLat_rad_sh = data_rd.gei_coLat_rad_sh[idx_dat]
	data.gei_lon_rad_sh   = data_rd.gei_lon_rad_sh[idx_dat]
	data.dbR_sh           = data_rd.dbR_sh[idx_dat]
	data.dbTheta_sh       = data_rd.dbTheta_sh[idx_dat]
	data.dbPhi_sh         = data_rd.dbPhi_sh[idx_dat]

; Expand on input data locations using the shifted GEI coords
	Print,'Generating Basis Set at Input Data Locations....'
	wait,0.001
 	ampere_setupSHFns, 1.0, data.gei_coLat_rad_sh, data.gei_lon_rad_sh, $   ; For Sth, 90->180 deg
			kMax, mMax, $
			minTheta = minTheta, $
			maxTheta = maxTheta, $
			bThBFnArr = dYkmDThBFns, $
			bPhBFnArr = dYkmDPhBFns, $
			YBFnArr = YkmBFns, $
			mArr = outMValues, $
			nArr = outNkValues, $
			kArr = outKValues, /bc2

    if keyword_set ( plot_bFns ) then begin            ; set this keyword to plot basis set - used for diagnostics
	    !p.multi = [0,3,4]
	    !p.charSize = 2.0
	    device, decomposed = 0
	    window, 0, xSize = 1200, ySize = 800
	    !p.background = 255

        m1  = 0
        k1  = 6
;	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
;	    		title = 'table lookup (m=2,k=5, Y)', $
;	    		psym=4,color=0
;	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
;	    		title = 'table lookup (m=2,k=5), dYdTh', $
;	    		psym=4,color=0
;	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
;	    		title = 'table lookup (m=2,k=5), dYdPh', $
;	    		psym=4,color=0

        m2  = 0
        k2  = 8
;	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
;	    		title = 'table lookup (m=0,k=7)', $
;	    		psym=4,color=0
;	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
;	    		title = 'table lookup (m=0,k=7)', $
;	    		psym=4,color=0
;	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues eq m2 and outkValues eq k2),*], $
;	    		title = 'table lookup (m=0,k=7)', $
;	    		psym=4,color=0

	    !p.multi = 0.0
	    !p.charSize = 1.0

    endif

	Print,'Fitting dB Data to Basis Set....'
	wait,0.001
;	Fit |dB| to dB.grad Ykm in the shifted GEI system
;	-----------------------------------------
	dBMag   = sqrt(data.dBTheta_sh^2+data.dBPhi_sh^2)
    If (sigma ne 1.) then begin
; Apply weighting [after Dec 10 discussion at APL]
     w_idx=where(dBMag lt thresh)
     If w_idx(0) gt -1 then begin
      data[w_idx].dBTheta_sh=data[w_idx].dBTheta_sh/sigma         ; CHECK THIS
      data[w_idx].dBPhi_sh=data[w_idx].dBPhi_sh/sigma
     end
    end
;
    bFuncs  = temporary ( $
                [[ dYkmDThBfns, -dYkmDPhBfns ], $
                 [ dYkmDPhBfns,  dYkmDthBfns ]])

    print, 'bFuncs is', n_elements(bFuncs[*])*16/(1024.0^2),' MB'
; Fit to shifted dBs
    coeffs_ = la_least_squares ( bFuncs, [data.dbTheta_sh, data.dbPhi_sh], $
               status = stat, $
               method = 3, $
               /double, rank = rank )
    idx_ls=where(finite(coeffs_) eq 0, cnt_ls)         ; trap for NaN in solver
    if cnt_ls ne 0 then begin                          ; if NaN, use the slower method
     coeffs_ = la_least_squares ( bFuncs, [data.dbTheta_sh, data.dbPhi_sh], $
               status = stat, $
               method = 2, $
               /double, rank = rank )
    endif
    if stat ne 0 then print, 'ERROR: la_least_squares threw an error code: ', stat

	Print,'Calculating FAC at Input Data Locations....'
    fit = bFuncs ## coeffs_
; Shifted fit dBs
	fit_dBTheta_sh = fit[0:n_elements(dBMag)-1]  ; 1st half of fit
	fit_dBPhi_sh   = fit[n_elements(dBMag):*]    ; 2nd half of fit

	Re    = 6.371d6
	R     = Re + 780.0d3
	u0    = 4.0*!dpi*1.0d-7
	jPar  = (YkmBFns ## $
	        (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}

; jPar is an array[1,Num of dB values]
; plot dbTh and dbPh by Iridium orbit track
	if plt_tracks then begin
	 !p.multi=[0,1,6,0]
	 window,winnum,xsize=1100,ysize=800,title='db_Theta by Iridium Track'
	 For  tt=0,5 do begin
      idx=where(data.ipln eq tt)
      plot,data[idx].dbTheta_sh,background=255,color=0
      oplot,fit_dbTheta_sh[idx],color=0,psym=2
     end
     winnum++
	 window,winnum,xsize=1100,ysize=800,title='db_Phi by Iridium Track'
	 For  tt=0,5 do begin
      idx=where(data.ipln eq tt)
      plot,data[idx].dbPhi_sh,background=255,color=0
      oplot,fit_dbPhi_sh[idx],color=0,psym=2
     end
     winnum++
	end
; Calc rms error for these
	err_dbTheta=total( (data.dbTheta_sh-fit_dbTheta_sh)^2 )
	err_dbTheta=sqrt(err_dbTheta/n_elements(fit_dbTheta_sh))
	err_dbPhi=total( (data.dbPhi_sh-fit_dbPhi_sh)^2 )
	err_dbPhi=sqrt(err_dbPhi/n_elements(fit_dbPhi_sh))
	Print,'RMS Error : dbTheta = ',err_dbTheta
	Print,'RMS Error : dbPhi = ',err_dbPhi

	Print,'Generating Basis Set over Uniform Grid....'
	wait,0.001
;	generate basis fns at regular grid (shifted)
;	----------------------------------

	if south then idx_dat=where(geiGrid_coLat_rad_sh*180.0/!dpi ge minTheta) else $
	              idx_dat=where(geiGrid_coLat_rad_sh*180.0/!dpi le maxTheta)

   	ampere_setupSHFns, 1.0, geiGrid_coLat_rad_sh[*], geiGrid_lon_rad_sh[*], $
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
	wait,0.001
    bFuncs_grid  = temporary($
                    [[dYkmDThBfns_grid, -dYkmDPhBfns_grid ], $
                     [dYkmDPhBfns_grid,  dYkmDthBfns_grid ]])

  	jParAACGM	= (YkmBFns_grid ## $
         (coeffs_[n_elements(YkmBFns_grid[*,0]):*]*(-outNkValues_grid*(outNkValues_grid+1.0))))$
                         /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}

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
	wait,0.001
; Need to rotate the dBs back
	rev_rot_mat=transpose(rot_mat)       ; undo rotation by transpose
    geopack_bspcar,geiGrid_coLat_rad_sh, geiGrid_lon_rad_sh, $     ; conv dB to XYZ comp
           dBR_GEI_grid_sh, dBTheta_GEI_grid_sh, dBPhi_GEI_grid_sh ,$
           vx_a,vy_a,vz_a
	For ii=Long(0),n_elements(dBTheta_GEI_grid_sh)-1 do begin
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
     geopack_bcarsp, gei_xGrid[ii], gei_yGrid[ii], gei_zGrid[ii], $
          dbx, dBy, dBz, vbr,vbth,vbph
     dBTheta_GEI_grid[ii]=vbth
     dBPhi_GEI_grid[ii]=vbph
    end

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
;    vec_geo2aacgm,path,aacgm_yr,geog_R_km-Re_km,$
;                       gLat_a,gLon_a,dBTheta_GEI_grid,dBPhi_GEI_grid,$
;                       mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph

; ***** ***** ***** ***** ***** END Uniform Grid Calc ***** ***** ***** ***** *****

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
                   data.dBTheta, data.dBPhi, $
                   mlat_a_in, mlon_a_in, $
                   mth_vec_gth_in, mph_vec_gth_in, mth_vec_gph_in, mph_vec_gph_in, err,/to_aacgm

;    vec_geo2aacgm,path,aacgm_yr,geog_R_km_in-Re_km,$
;                       gLat_a_in,gLon_a_in,data.dBTheta,data.dBPhi,$
;                       mlat_a_in,mlon_a_in,$
;                       mth_vec_gth_in,mph_vec_gth_in,mth_vec_gph_in,mph_vec_gph_in

;	plot stuff
;	----------
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
;	plot the db vectors
;	-------------------
;	plot the FAC maps
    hem='north'
    if south then hem='south'
  	fac_div=mx_fac/10.
	pole=90
	sgn=1.
	if south then begin
	 pole=-90
	 plt_colat=-plt_colat
	 sgn=-1.
	end

	!p.backGround = 0
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
            capSize = abs(plt_coLat)

; ***** Check if it should be 90.0-abs(m_lat_a)	*****
	plt_dat, 90.0-mlat_a[*], ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	         n_elements(mlat_a[*] ), -mth_vec_gth, mph_vec_gph, [1,1], [1,2], south,$
			 title = 'Fitted data: AACGM-MLT', $
            capSize = abs(plt_coLat)

	loadct, 13, file = path + 'davect.tbl',/silent
   	map_set, pole, 0, 0, /ortho, /iso, $
     limit = [ pole - plt_coLat, 0, pole, 360 ], xmargin=[1,1], ymargin=[1,10], $
	 /noborder, /advance, title = 'Fitted data: jPar AACGM-MLT'
	jParTmp	= jParAACGM[*]
	fc_idx=where(abs(jParTmp) lt mn_fac)
	if fc_idx(0) gt -1 then jParTmp(fc_idx)=0.0      ; Check for min val, set to white
	fc_idx=where(jParTmp gt mx_fac)
	if fc_idx(0) gt -1 then jParTmp(fc_idx)=mx_fac
	fc_idx=where(jParTmp lt -mx_fac)
	if fc_idx(0) gt -1 then jParTmp(fc_idx)=-mx_fac
	lonTmp = ((aacgmGrid_lon_deg[*]/15.0+mltShift) mod 24 )*15.0 $
			+ randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*1e-5-0.5e-5
	latTmp = 90.0-aacgmGrid_coLat_deg[*]
   	contour, jParTmp, $
			lonTmp, $
			latTmp, $
		 	c_labels=fltarr(n_elements(jLevels))+1,$
        	/over, levels = jLevels, $
           	c_colors = colors, /fill, /irreg
	map_grid, label = 1, $
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360], $
			lonlab=pole-plt_colat, $
			latLab = 45, $
			latDel = sgn*10.0

	If plt_png eq 1 then begin
	    png_file=path+'amp_dB&FAC_'+yr_str+mnth_str+dy_str+' at '+hr_str+mn_str+sc_str+'_'+hem+'.png'
    	image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
    	Write_PNG,png_file,image,r,g,b
    	Print,'PNG for FACs written to ', png_file
    end

; Window 4, big GEO-UT Plot for NSF
;	window, winnum, xSize = 750, ySize = 700,title='GEO_UT FAC for '+date_str+'  '+tme_str
;	winnum++
;	!p.backGround = 0
;	!p.multi = 0
;	loadct, 13, file = path + 'davect.tbl'
;	If south eq 0 then begin                  ; Nth Hemisphere
;    map_set, 90, 0, ln_sh, /ortho, /iso, $
;      limit = [90.0-(plotCapSize + 2), 0, 90, 360 ], xmargin=[1,7], ymargin=[5,5], $
;	  /noborder, /advance, title = 'AMPERE: j_par GEO - UT'
;     jParTmp	= jParAACGM[*]
;     fc_idx=where(abs(jParTmp) lt mn_fac)
;     if fc_idx(0) gt -1 then jParTmp(fc_idx)=0.0
;     lonTmp = ((geog_lon_rad[*]*!radeg/15.0+0.0) mod 24 )*15.0 $
;	          + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*1.0e-5-0.5e-5
;	 latTmp = 90.0-geog_coLat_rad[*]*!radeg
;   	 contour, jParTmp, $
;			lonTmp, $
;			latTmp, $
;		 	c_labels=fltarr(n_elements(jLevels))+1,$
;        	/over, levels = jLevels, $
;           	c_colors = colors, /fill, /irreg
;    end
;
;	If south eq 1 then begin                  ; Sth Hemisphere
;     map_set, -90, 0, ln_sh, /ortho, /iso, $
;      limit = [-90, 0, -90+(plotCapSize + 2), 360 ], xmargin=[1,7], ymargin=[5,5], $
;	  /noborder, /advance, title = 'AMPERE: j_par GEO - UT'
;     jParTmp	= jParAACGM[*]
;     fc_idx=where(abs(jParTmp) lt mn_fac)
;     if fc_idx(0) gt -1 then jParTmp(fc_idx)=0.0
;     lonTmp = ((geog_lon_rad[*]*!radeg/15.0+0.0) mod 24 )*15.0 $
;	          + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*1.0e-5-0.5e-5
;	 latTmp = -90.0+geog_coLat_rad[*]*!radeg
;   	 contour, jParTmp, $
;			lonTmp, $
;			latTmp, $
;	          c_labels=fltarr(n_elements(jLevels))+1,$
;              /over, levels = jLevels, $
;           	  c_colors = colors, /fill, /irreg
;    end
;
;	map_grid, label = 1, $
;			lonNames	= ['0/24','6','12','18',''], $
;			lons = [0,90,180,270,360]+ln_sh, $
;			latLab = 45+ln_sh, $
;			latDel = 20.0
;	map_continents,color=60
;
; Color Bar
;    width = 0.04               ; changes the width of the color bar
;    XStrt = 0.88               ; changes the X location for the color bar, assuming vertical bar
;    YStrt = 0.30               ; changes the Y starting point for the color bar, assuming vertical bar
;    YEnd = 0.65                ; changes the height of the color bar, assuming vertical bar
;    cl=255
;    ttle = 'uAm!u-2!d'
;    cbar = Obj_New('Colorbar',$
;    	Position = [XStrt, YStrt, XStrt+width, YEnd],$
;	    Vertical = 1,$
;        Color=cl,$
;	    Range = [-mx_fac,mx_fac], format='(f4.1)',TickLen=-0.1,major=10);, $
	 ;   Title = Ttle)
;    cbar->Draw

;	If plt_png eq 1 then begin
;	    png_file=path+'GEO_UT_FAC_'+dy_str+mnth_str+yr_str+' at '+hr_str+'_'+mn_str+'_'+sc_str+'.png'
;    	image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
;    	Write_PNG,png_file,image,r,g,b
;    	Print,'PNG for GEO_UT FAC written to ', png_file
;    end


	Print,'Min:Max FAC = ',min(jParTmp),'  ',max(jParTmp)
	!p.multi = 0
end

