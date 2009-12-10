pro amp_fit, $
		numerical = numerical, $
        plot_bFns = plot_bFns

	if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN62') then begin
		plotDev	= 'win'
		path	= 'd:\cwac\hi_res\davidg\'
		pnmPath	= path + 'jpar_ver2\pnmsavs\pnmSav'
		aacgm_set_path,path
	endif else begin                        ; other systems (linux, darwin etc.)
		plotDev = 'X'
		path	= '~/code/ampereFit/idl/'
		pnmPath	= path + 'pnmSavs/pnmSav'
	endelse

	capSize	= 60.0
	plotCapSize	= 50.0

; Input data file
;	fileName = path + '20050515_a_RevB.dat'
;    savFileName = '~/code/ampereFit/data/20091023Amp_invert_test2.sav'
    savFileName = path+'20091124Amp_invert.sav'

; Time interval (UT)
	sHr = 12.0+1./6.
	eHr = 12.0+3./6.

; Hemisphere switch
    south=0               ; 0=North, 1=South

; Basis functions
	kMax    = 15
	mMax    = 5          ; 0 to 5

; Controls final solution grid
	nLatGrid	= 50
	nLonGrid	= 24

; Max FAC for plot
	mx_fac=1.0

; Select out relevant data from time and Hemisphere constraints
	;read_ampere_dlg, fileName, sHr, eHr, data, t_arr, capSize, $
	;   	 yrSec = yrSec, $
	;   	 year = year, $
	;   	 month = month, $
	;   	 day = day, $
	;	 avgYrSec = avgYrSec, $
	;   	 avgEpoch = avgEpoch

; Restore saveset of input data
; Calc Lat,Lon shift to place avg orbit intersection point at centre
; Rotate GEI coords and dBs to centred system
    read_ampere_sav, $
        savFileName = savFileName, $
        capSize = capSize, $
        dataOut = data, $         ; GEI are all shifted in structure data
        sHr = sHr, $
        eHr = eHr, $
        year = year, $
        month = month, $
        day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        south = south, $
        rot_mat=rot_mat, $
        fillPoleLine = 0

; Generate final output data grid
	aacgm_grid, $
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
; Expand on input data locations using the shifted GEI coords
	schaBasisFunctions, kMax, mMax, capSize, data.gei_coLat_rad, data.gei_lon_rad, $
	     YkmBFns=YkmBFns, $
	   	 dYkmDthBFns=dYkmDthBFns, $
	   	 dYkmDphBFns=dYkmDphBFns, $
	     OUTNKVALUES=outNkValues, $
	   	 OUTKVALUES=outKValues, $
	   	 OUTMVALUES=outMValues, $
	   	 pnmPath = pnmPath;, /evenSet
; Uses the shifted GEI locations
 	ampere_setupSHFns, 1.0, data.gei_coLat_rad, data.gei_lon_rad, $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns2, $
			bPhBFnArr = dYkmDPhBFns2, $
			YBFnArr = YkmBFns2, $
			mArr = outMValues2, $
			nArr = outNkValues2, $
			kArr = outKValues2;, $
            ;/bc2

    if keyword_set ( plot_bFns ) then begin
	    !p.multi = [0,3,4]
	    !p.charSize = 2.0
	    device, decomposed = 0
	    window, 0, xSize = 1200, ySize = 800
	    !p.background = 255

        m1  = 0
        k1  = 6
	    plot, data.gei_coLat_rad*!radeg, ykmbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5, Y)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdTh', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns[where(outmvalues eq m1 and outkValues eq k1),*], $
	    		title = 'table lookup (m=2,k=5), dYdPh', $
	    		psym=4,color=0

	    plot, data.gei_coLat_rad*!radeg, ykmbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1),*], $
	    		title = 'numerical (m=2,k=5), Y', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1 ),*], $
	    		title = 'numerical (m=2,k=5), dYdTh', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[where(outmvalues2 eq m1 and outkValues2 eq k1 ),*], $
	    		title = 'numerical (m=2,k=5), dYdPh', $
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

	    plot, data.gei_coLat_rad*!radeg, ykmbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdthbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0
	    plot, data.gei_coLat_rad*!radeg, dykmdphbfns2[where(outmvalues2 eq m2 and outkValues2 eq k2),*], $
	    		title = 'numerical (m=0,k=7)', $
	    		psym=4,color=0

	    !p.multi = 0.0
	    !p.charSize = 1.0

    endif

	if keyword_set(numerical) then begin
        ;iiSort  = sort ( outNkValues2 )
		YkmBFns	= temporary ( YkmBFns2);[iiSort,*] )
		dYkmDthBFns	= temporary ( dYkmDthBFns2);[iiSort,*] )
		dYkmDphBFns	= temporary ( dYkmDphBFns2);[iiSort,*] )
		outNkValues	= temporary ( outNkValues2);[iiSort] )
	endif

;	Fit |dB| to dB.grad Ykm in the shifted GEI system
;	-----------------------------------------
	dBMag   = sqrt(data.dBTheta^2+data.dBPhi^2)
    bFuncs  = temporary ( $
                [ [ dYkmDThBfns, -dYkmDPhBfns ], $
                  [ dYkmDPhBfns,  dYkmDthBfns ] ] )

    print, 'bFuncs is', n_elements(bFuncs[*])*16/(1024.0^2),' MB'
; Fit to shifted dBs
    coeffs_ = la_least_squares ( bFuncs, [data.dbTheta, data.dbPhi], $
               status = stat, $
               method = 3, $
               /double, rank = rank )
    if stat ne 0 then begin
        print, 'ERROR: la_least_squares threw an error code: ', stat
    endif

    fit = bFuncs ## coeffs_
; Shifted fit dBs
	fit_dBTheta = fit[0:n_elements(dBMag)-1]  ; 1st half of fit
	fit_dBPhi   = fit[n_elements(dBMag):*]    ; 2nd half of fit

	Re    = 6.371d6
	R     = Re + 780.0d3
	u0    = 4.0*!dpi*1.0d-7
	jPar  = (YkmBFns ## $
	        (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
	                        /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}
; jPar is an array[1,Num of dB values]

;	generate basis fns at regular grid (shifted)
;	----------------------------------
    sphcar, geiGrid_R_km[*], geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
            gei_xGrid, gei_yGrid, gei_zGrid,/to_rect       ; XYZ coord of GEI grid (unshifted)
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
     sphcar,r_sh,th_sh,ph_sh,xsh,ysh,zsh,/to_sphere  ; r,thet,phi coord of GEI grid (shifted)
     geiGrid_coLat_rad_sh[ii]=th_sh
   	 geiGrid_lon_rad_sh[ii]=ph_sh
	end

	schaBasisFunctions, kMax, mMax, capSize, geiGrid_coLat_rad_sh[*], geiGrid_lon_rad_sh[*], $
         YkmBFns=YkmBFns_grid, $
		 dYkmDthBFns=dYkmDthBFns_grid, $
		 dYkmDphBFns=dYkmDphBFns_grid, $
		 pnmPath = pnmPath;,$
         ;/oddSet

	ampere_setupSHFns, 1.0, geiGrid_coLat_rad_sh[*], geiGrid_lon_rad_sh[*], $
			kMax, mMax, $
			minTheta = 0.0, $
			maxTheta = capSize, $
			bThBFnArr = dYkmDThBFns_grid2, $
			bPhBFnArr = dYkmDPhBFns_grid2, $
			YBFnArr = YkmBFns_grid2, $
			mArr = outMValues_grid2, $
			nArr = outNkValues_grid2, $
			kArr = outKValues_grid2;, $
            ;/bc2

	if keyword_set(numerical) then begin
		YkmBFns_grid	= temporary ( YkmBFns_grid2 )
		dYkmDthBFns_grid	= temporary ( dYkmDthBFns_grid2 )
		dYkmDphBFns_grid	= temporary ( dYkmDphBFns_grid2 )
	endif

    bFuncs_grid  = temporary($
                    [[dYkmDThBfns_grid, -dYkmDPhBfns_grid ], $
                     [dYkmDPhBfns_grid,  dYkmDthBfns_grid ]])

  	jParAACGM	= (YkmBFns_grid ## $
         (coeffs_[n_elements(YkmBFns[*,0]):*]*(-outNkValues*(outNkValues+1.0))))$
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

; Need to rotate the dBs back
	rev_rot_mat=transpose(rot_mat)       ; undo rotation by transpose
    bspcar,geiGrid_coLat_rad_sh, geiGrid_lon_rad_sh, $     ; conv dB to XYZ comp
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
     bcarsp, gei_xGrid[ii], gei_yGrid[ii], gei_zGrid[ii], $
          dbx, dBy, dBz, vbr,vbth,vbph
     dBTheta_GEI_grid[ii]=vbth
     dBPhi_GEI_grid[ii]=vbph
    end

; Take GEI -> GEOG -> AACGM
    Re_km = 6357.0
    R_km = 110.0 + Re_km
; conv GEI r,thet,phi -> GEI XYZ
    sphcar,geiGrid_R_km, geiGrid_coLat_rad, geiGrid_lon_rad,$
           gei_x, gei_y, gei_z, /to_rect
; conv GEI XYZ -> GEOG XYZ
    geoPack_conv_coord, gei_x, gei_y, gei_z, $
                        geo_x, geo_y, geo_z, $
                       /from_gei, /to_geo

    epoch = fltArr(n_elements(gei_x[*])) + Avgepoch

; conv GEOG XYZ -> GEOG r,thet,phi
    sphcar,geog_R_km, geog_coLat_rad, geog_lon_rad, $
           geo_x, geo_y, geo_z, /to_sphere

    gLat_a=90.0-geog_coLat_rad*!radeg
    gLon_a=geog_lon_rad*!radeg
    aacgm_conv_vec, gLat_a, gLon_a, geog_R_km-Re_km, $
                    dBTheta_GEI_grid, dBPhi_GEI_grid, $
     mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph,err,/to_aacgm

;	rotate_gei_to_aacgm, geiGrid_R_km[*], geiGrid_coLat_rad[*], geiGrid_lon_rad[*], $
;		aacgmGrid_R_km[*], aacgmGrid_coLat_deg[*]*!dtor, aacgmGrid_lon_deg[*]*!dtor, $
;		dBTheta_GEI_grid, dBPhi_GEI_grid, $
;		aacgm_dbTh = aacgm_dbTh, $
;		aacgm_dbPh = aacgm_dbPh, $
;		year = year, $
;		epoch = avgEpoch

;	plot stuff
;	----------
	set_plot, plotDev
	device,decomposed = 0
	Loadct,0,/silent
	!p.background = 255
	winNum	= 1

;	plot the db vectors
;	-------------------
	!p.multi = [0,2,2]
	window, winnum, xSize=700, ySize=700
	winNum++
	plt_dat, data.gei_coLat_rad*!radeg, data.gei_lon_rad*!radeg, $
             n_elements(data.gei_coLat_rad), -data.dbTheta, data.dbPhi, [1,1], [1,2], $
             title = 'GEI - Shifted Input dB Data', $
             satNu = data.iSat, $
             capSize = plotCapSize
	plt_dat, data.gei_coLat_rad*!radeg, data.gei_lon_rad*!radeg, $
	         n_elements(data.gei_coLat_rad), -fit_dBTheta, fit_dBPhi, [1,1], [1,2], $
			 title = 'GEI - Fitted Input dB Data (shifted)', $
             capSize = plotCapSize
	plt_dat, geiGrid_coLat_rad[*]*!radeg, geiGrid_lon_rad[*]*!radeg, $
	         n_elements (geiGrid_coLat_rad[*] ), -dBTheta_GEI_grid, dBPhi_GEI_grid, [1,1], [1,2], $
			 title = 'GEI - on regular grid (uniform AACGM equiv)', $
             capSize = plotCapSize
; ***** Check if it should be 90.0-abs(m_lat_a)	*****
	plt_dat, 90.0-mlat_a[*], ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	         n_elements(mlat_a[*] ), -mth_vec_gth, mph_vec_gph, [1,1], [1,2], $
			 title = 'AACGM-MLT', $
            capSize = plotCapSize
;	plt_dat, aacgmGrid_coLat_deg[*], aacgmGrid_lon_deg[*], $
;	         n_elements ( geiGrid_coLat_rad[*] ), -aacgm_dbTh, aacgm_dbPh, [1,1], [1,2], $
;			 title = 'AACGM-LON', $
;            capSize = plotCapSize

;	plot the FAC maps
;	-----------------
	loadct, 13, file = path + 'davect.tbl'
	window, winnum, xSize = 650, ySize = 650
	winnum++
	!p.backGround = 0
	!p.multi = [0,2,2,0]
	!p.font = 0
	map_set, 90, 0, 0, /ortho, /iso, $
	    limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	    /noborder, /advance, title = 'jPar GEI (shifted)'
	jLevels = (fIndGen(21)-10.)*0.1  ;1.0
	colors  = bytScl ( jLevels, top = 253 ) + 1
	oldbck=!p.background
	!p.background=0
	jpar=reform(jpar)
    triangulate, data.gei_lon_rad*!radeg, 90-data.gei_colat_rad*!radeg, $
        tri, sphere = sphere
	contour, jPar, data.gei_lon_rad*!radeg, 90.0-data.gei_coLat_rad*!radeg, $
	   	 	c_labels=fltarr(n_elements(jLevels))+1,$
	       	/irreg, /over, levels = jLevels, $
	          	c_colors = colors, /fill
	map_grid, label = 1, $
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360], $
			latLab = 45, $
			latDel = 20.0


   	map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar AACGM-MLT'
	jParTmp	= jParAACGM[*]
	lonTmp = ((aacgmGrid_lon_deg[*]/15.0+mltShift) mod 24 )*15.0 $
			+ randomu ( sysTime(/sec ), n_elements(jParTmp), /uni ) * 1e-5
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
			latLab = 45, $
			latDel = 20.0

   	map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar AACGM'
	jParTmp	= jParAACGM[*]
	lonTmp = ((aacgmGrid_lon_deg[*]/15.0) mod 24 )*15.0 $
			+ randomu ( sysTime(/sec ), n_elements(jParTmp), /uni ) * 1e-5
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
			latLab = 45, $
			latDel = 20.0

	avg_hr = (sHr+eHr)/2.0
	ln_sh=(24.0-avg_hr)*15.0

	If south eq 0 then begin                  ; Nth Hemisphere
     map_set, 90, 0, ln_sh, /ortho, /iso, $
      limit = [90.0-(plotCapSize + 2), 0, 90, 360 ], $
	 /noborder, /advance, title = 'jPar GEOG - UT'
   	 contour, jParaacgm, geog_lon_rad*!radeg, 90.0-geog_coLat_rad*!radeg, $
		 	c_labels=fltarr(n_elements(jLevels))+1,$
        	/over, levels = jLevels, $
           	c_colors = colors, /fill, /irreg
    end

	If south eq 1 then begin                  ; Sth Hemisphere
     map_set, -90, 0, ln_sh, /ortho, /iso, $
      limit = [-90, 0, -90+(plotCapSize + 2), 360 ], $
	  /noborder, /advance, title = 'jPar GEOG - UT'
   	 contour, jParaacgm, geog_lon_rad*!radeg, -90.0+geog_coLat_rad*!radeg, $
	          c_labels=fltarr(n_elements(jLevels))+1,$
              /over, levels = jLevels, $
           	  c_colors = colors, /fill, /irreg
    end

	map_grid, label = 1, $
			lonNames	= ['0/24','6','12','18',''], $
			lons = [0,90,180,270,360]+ln_sh, $
			latLab = 45+ln_sh, $
			latDel = 20.0
	map_continents,color=60

;	map_set, 90, 0, 0, /ortho, /iso, $
;     limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ], $
;	 /noborder, /advance, title = 'jPar AACGM-LON [GEI]'
;	jParTmp	= jParAACGM[*]
;	lonTmp = geiGrid_lon_rad[*]*!radeg
;	latTmp = 90.0-geiGrid_coLat_rad[*]*!radeg
;	contour, jParTmp, $
;			lonTmp, $
;			latTmp, $
;		 	c_labels=fltarr(n_elements(jLevels))+1,$
;       	/over, levels = jLevels, $
;           	c_colors = colors, /fill, /irreg
;	map_grid, label = 1, latDel = 10.0


;	plot the data locations
;	-----------------------
;	!p.multi = [0,2,1]
;	window, winnum, xSize = 650, ySize = 450
;	winnum++
; 	map_set, 90, 0, 0, $
;		/ortho, $
;  		/iso, $
;    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
;  		/noborder,$
;  		/advance, title = 'AACGM-LON'
;	map_grid, label = 1, latDel = 10.0
;	plots, data.aacgm_lon_rad*!radeg, 90.0-data.aacgm_coLat_rad*!radeg, psym = 4
;
;   	map_set, 90, 0, 0, $
;		/ortho, /iso, $
;   	limit = [90.0-(plotCapSize + 2), 0, 90, 360 ],$
;   		/noborder,/advance, title = 'GEI (shifted)'
;	map_grid, label = 1, latDel = 10.0
;	plots, data.gei_lon_rad*!radeg, 90.0-data.gei_coLat_rad*!radeg, psym = 4
;
; 	map_set, 90, 0, 0, $
;		/ortho, /iso, $
;    	limit = [90.0-(plotCapSize + 2), 0, 90, 360 ],$
;   		/noborder,/advance, title = 'AACGM-MLT'
;	map_grid, label = 1, latDel = 10.0
;	plots, data.mlt*15.0, 90.0-data.aacgm_coLat_rad*!radeg, psym = 4
;
;	map_set, 90, 0, 0, $
;		/ortho, $
;   		/iso, $
;    	limit = [ 90.0 - ( plotCapSize + 2 ), 0, 90, 360 ],$
;   		/noborder,$
;   		/advance, title = 'GE0G'
;	map_grid, label = 1, latDel = 10.0
;	plots, data.geog_lon_rad*!radeg, 90.0-data.geog_coLat_rad*!radeg, psym = 4

	!p.multi = 0
 	!p.background = oldbck
; stop
 Print,'Finished'
end

