; +
; AmpereFit
; ---------
;
; Library to take Iridium vector (horizontal only) dB
; field and create a FAC and dB map on regular grid.
;
; C.L. Waters and D.L. Green
; colin.waters@newcastle.edu.au
;
; Dec 2009
;
; Full sphere version CLW - Sept 2010

pro amp_fit_full, sHr, eHr, $
	datFileName = datFileName, $
	kmax = kmax, mmax = mmax, $
	nLonGrid = nLonGrid, $                  ; deleted nLat - set to 1 deg
	mn_fac = mn_fac, mx_fac=mx_fac, $
	plt_png = plt_png, $
	debug = debug, $
	sepoch = sepoch, eepoch = eepoch, $			; start/end epoch times
	extend = extend, $                      ; get data across day if required
	SpacRes = SpacRes                       ; 1/2 wavelength res in deg (Latitude)

	@amp_fit_paths_full								; specify dir paths

; error trap routine - if error, then return to calling routine
  if (debug eq 0) then begin
    status=0
    catch,error_status
    if (error_status ne 0) then begin
      status=1
      Print,'Error in amp_fit_f.pro - returning'
      catch,/cancel
      return
    endif
  end
	; Check for reasonable inputs
	; ---------------------------

    if size(sHr,/type) eq 0 $
            or size(eHr,/type) eq 0 then begin

        print, 'INPUT ERROR: need sHr and eHr '
        print, '    e.g., amp_fit, 4, 5'
;        stop
    endif

	; Set default keyword values
	; --------------------------
  aacgm_set_path, aacgmpath

	@amp_fit_defaults

	; Read input data for full sphere
	Print,'Reading AMPERE Data File....'

  read_ampere_data_full, sHr, eHr, $
        datFileName = datFileName, $
        dataOriginal = dataOriginal, $
        dataShifted = dataShifted, $
        year = year, month = month, day = day, $
        avgYrSec = avgYrSec, $
        avgEpoch = avgEpoch, $
        rot_mat=rot_mat, $
        debug = debug, $
        extend = extend, $
        SpacRes = SpacRes, $
        status = status
; window,0 used in here if debug=1
; avgEpoch here is for 0 UT, avgYrSec carries the hr, mn info

; Trap for errors
  if (debug eq 0) then begin
   if (status ne 0) then return
  endif

; Create basis functions at shifted GEI coords of input data
; -----------------------------------------------------------

  Print,'Generating Basis Set at Input Data Locations....'

  ; For Sth, 90->180 deg

  ampere_setupSHFns_full, 1.0, dataShifted.gei_coLat_rad, dataShifted.gei_lon_rad, $
      kMax, mMax, $
      brBFnArr = brBFnArr, $
      bThBFnArr = dYkmDThBFns, $
      bPhBFnArr = dYkmDPhBFns, $
      YBFnArr = YBFnArr, $
      mArr = outMValues, $
      nArr = outNkValues, $
      kArr = outKValues

; Apply data qual weighting to BFs
  For nb=0,n_elements(outKValues)-1 do begin
   brBFnArr[nb,*]=brBFnArr[nb,*]/dataShifted[*].qual         ; see Num. Recipes p673
   dYkmDThBfns[nb,*]=dYkmDThBfns[nb,*]/dataShifted[*].qual
   dYkmDPhBfns[nb,*]=dYkmDPhBfns[nb,*]/dataShifted[*].qual
  end

  Print,'Fitting dB Data to Basis Set....'
  bFuncs  = temporary ( $
            [[ dYkmDThBfns, -dYkmDPhBfns ], $    ; thet comps of thet,Phi eqns
             [ dYkmDPhBfns,  dYkmDthBfns ]])     ; phi  comps of thet,phi eqns

; Free up memory for some arrays [IDL 8]
  dYkmDThBFns=!null
  dYkmDPhBFns=!null

  print, 'bFuncs is', n_elements(bFuncs[*])*16/(1024.0^2),' MB'

  ; Fit to shifted dBs
  ; ------------------
; Apply data qual weighting to dB data - see Num Recipes p673
  coeffs_ = la_least_squares ( bFuncs, $
        [dataShifted.bTheta_GEI/dataShifted.qual, dataShifted.bPhi_GEI/dataShifted.qual], $
        status = stat, $
        method = 2, $
        /double, rank = rank )

  if stat ne 0 then begin
    print, 'ERROR: la_least_squares', stat
;   stop
  endif

  fit = bFuncs ## coeffs_         ; get fitted data
  bfuncs=!null                    ; free memory taken by bfuncs

; Do dBR data fit - wait for a decision from APL
;  coeffs_r = la_least_squares ( brBFnArr, dataShifted.br_GEI/dataShifted.qual, $
;        status = stat, $
;        method = 2, $
;        /double, rank = rank )

;  fit_bR_GEI = brBFnArr ## coeffs_r
;  brBFnArr=!null                      ; free memory

  ; Extract theta and phi component of fit from solution
  ; ----------------------------------------------------

  fit_bTheta_GEI = fit[0:n_elements(dataShifted.bTheta_GEI)-1]
  fit_bPhi_GEI   = fit[n_elements(dataShifted.bTheta_GEI):*]
;
  ; Calc rms error for these
  ; ------------------------

;  err_dbR = total ( ( dataShifted.bR_GEI-fit_bR_GEI)^2 )
;  err_dbR = sqrt ( err_dbR/n_elements(fit_bR_GEI))
  err_dbTheta = total ( (( dataShifted.bTheta_GEI-fit_bTheta_GEI)/dataShifted.qual)^2 ) ; This is chiSq - Num Recipes p673
  err_dbTheta = sqrt ( err_dbTheta/n_elements(fit_bTheta_GEI))
  err_dbPhi = total ( (( dataShifted.bPhi_GEI-fit_bPhi_GEI)/dataShifted.qual)^2 )
  err_dbPhi = sqrt ( err_dbPhi/n_elements(fit_bPhi_GEI))

;  Print,'Min:Max input dbR = ',min(dataShifted.br_GEI),max(dataShifted.br_GEI)
;  Print,'Min:Max Fit dbR = ',min(fit_br_GEI),max(fit_br_GEI)
;  Print,'RMS Error : dbR = ',err_dbR
  Print,'Min:Max input dbTheta = ',min(dataShifted.bTheta_GEI),max(dataShifted.bTheta_GEI)
  Print,'Min:Max Fit dbTheta = ',min(fit_bTheta_GEI),max(fit_bTheta_GEI)
  Print,'RMS Error : dbTheta = ',err_dbTheta
  Print,'Min:Max input dbPhi = ',min(dataShifted.bPhi_GEI),max(dataShifted.bPhi_GEI)
  Print,'Min:Max Fit dbPhi = ',min(fit_bPhi_GEI),max(fit_bPhi_GEI)
  Print,'RMS Error : dbPhi = ',err_dbPhi

; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
  if debug then begin
  ; Plot along track comparison of the raw and fitted data
  ; ------------------------------------------------------
;   @load_colors
;   iPlot, view_grid = [2,3]
;   for i=0,5 do begin
;     iiTrack = where ( data.iPln eq i, iiTrackCnt )
;     iPlot, (data.bPhi_GEI)[iiTrack], $
;         sym_index = 4, $
;         view_number = i+1, $
;         /stretch_to_fit, $
;         lineStyle = 6
;     iPlot, fit_bPhi_GEI[iiTrack], $
;         thick = 2, /over, $
;         color = red
;   endfor
    !p.multi=[0,1,6,0]
    Loadct,0,/silent
;    window,1,xsize=1100,ysize=700,title='db_R by Iridium Track: data(solid), fit(dotted)'
;    For  tt=0,5 do begin
;      idx=where(dataShifted.ipln eq tt)
;      plot,dataShifted[idx].bR_GEI,background=255,color=200
;      oplot,fit_bR_GEI[idx],color=0,psym=3
;    end

    window,2,xsize=1100,ysize=700,title='db_Theta by Iridium Track: data(solid), fit(dotted)'
    For  tt=0,5 do begin
      idx=where(dataShifted.ipln eq tt)
      plot,dataShifted[idx].bTheta_GEI,background=255,color=200
      oplot,fit_bTheta_GEI[idx],color=0,psym=3
    end

    window,3,xsize=1100,ysize=700,title='db_Phi by Iridium Track: data(solid), fit(dotted)'
    For  tt=0,5 do begin
      idx=where(dataShifted.ipln eq tt)
      plot,dataShifted[idx].bPhi_GEI,background=255,color=200
      oplot,fit_bPhi_GEI[idx],color=0,psym=3
    end
  endif             ; end debug for track dB comparison
; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
  fit_bTheta_GEI = !null
  fit_bPhi_GEI = !null

  @amp_fit_constants

; ############ NORTH HEMISPHERE AACGM CAP ###############################
;
; Create uniform grids in AACGM and GEI
; -------------------------------------
; 1. Calculate Uniform AACGM grid in Nth hemisphere cap [AACGM not defined at equatorial lats]
; 2. Conv to geog -> then to GEI coord

  south=0
  cap_sz_deg=60.
  aacgm_grid, cap_sz_deg, south, $             ; get north hemisphere cap grid
       aacgmGrid_R_km = aacgmGrid_R_km_nth, $
	   	 aacgmGrid_coLat_deg = aacgmGrid_coLat_deg_nth, $
	   	 aacgmGrid_lon_deg = aacgmGrid_lon_deg_nth, $
       geiGrid_R_km = geiGrid_R_km_nth, $
	     geiGrid_coLat_rad = geiGrid_coLat_rad_nth, $
	   	 geiGrid_lon_rad = geiGrid_lon_rad_nth, $
       xGrid_GEI=xGrid_GEI_nth, $
       yGrid_GEI=yGrid_GEI_nth, $
       zGrid_GEI=zGrid_GEI_nth, $
	   	 nLat = cap_sz_deg, nLon = nLonGrid, $
	   	 year = year, yrSec = avgYrSec, $
	   	 mltShift = mltShift, $
	   	 epoch = avgEpoch

; 3. Calc shifted GEI from uniform AACGM grid

  ampere_shift_coord, xGrid_GEI_nth[*], yGrid_GEI_nth[*], zGrid_GEI_nth[*], $   ; input GEI coords
      geiGridX_shifted_nth, geiGridY_shifted_nth, geiGridZ_shifted_nth, $       ; output shifted GEI coords
      geiGrid_rIrid_km_shifted_nth, geiGrid_coLat_deg_shifted_nth, geiGrid_lon_deg_shifted_nth, $  ; spherical coords of shifted coords
      rot_mat          ; rotation matrix

  xGrid_GEI_nth=!null
  yGrid_GEI_nth=!null
  zGrid_GEI_nth = !null

  geiGrid_coLat_rad_shifted_nth = geiGrid_coLat_deg_shifted_nth[*]*!dtor
  geiGrid_lon_rad_shifted_nth = geiGrid_lon_deg_shifted_nth[*]*!dtor

	; Generate basis fns at regular grid
	; ----------------------------------

	print,'Generating Basis Set over Northern Uniform Grid....'

  ampere_setupSHFns_full, 1.0, geiGrid_coLat_rad_shifted_nth[*], geiGrid_lon_rad_shifted_nth[*], $
			kMax, mMax, $
      brBFnArr = brBFns_grid_nth, $
			bThBFnArr = dYkmDThBFns_grid_nth, $
			bPhBFnArr = dYkmDPhBFns_grid_nth, $
			YBFnArr = YkmBFns_grid_nth, $
			mArr = outMValues_grid_nth, $
			nArr = outNkValues_grid_nth, $
			kArr = outKValues_grid_nth

	Print,'Calculating FAC over Northern Uniform Grid....'

  bFuncs_grid_nth  = temporary($
                 [[dYkmDThBfns_grid_nth, -dYkmDPhBfns_grid_nth ], $
                  [dYkmDPhBfns_grid_nth,  dYkmDthBfns_grid_nth ]])

  jParAACGM_nth	= (YkmBFns_grid_nth ## $
         (coeffs_[n_elements(YkmBFns_grid_nth[*,0]):*]*(-outNkValues_grid_nth*(outNkValues_grid_nth+1.0)))) $
                         /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}

	jParAACGM_nth	= reform(jParAACGM_nth, cap_sz_deg, nLonGrid)

  print,'Min:Max jPar - North = ',min(jParAACGM_nth),max(jParAACGM_nth)

  ; Reconstruct the dB vector from the fit
  ; --------------------------------------

  recon_dB_geiGrid  = bFuncs_grid_nth ## coeffs_
  bFuncs_grid_nth=!null      ; release this memory

  dBTheta_geiGrid_sh_nth = recon_dB_geiGrid[0:cap_sz_deg*nLonGrid-1]
  dBPhi_geiGrid_sh_nth   = recon_dB_geiGrid[cap_sz_deg*nLonGrid:*]
  dBR_geiGrid_sh_nth   = dBPhi_geiGrid_sh_nth * 0.0d0   ; dBR -> 0 for now
;  dBR_geiGrid_sh_nth   = brBFns_grid_nth ## coeffs_r   ; dBR - actual code - APL decision
;  brBFns_grid_nth=!null                      ; free memory

  ; Un-Shift the dB Vectors
  ; --------------------------------

  rev_rot_mat = transpose(rot_mat)

;  conv fitted dB vectors from shifted spherical to shifted XYZ components
  geopack_bspcar,geiGrid_coLat_rad_shifted_nth[*], geiGrid_lon_rad_shifted_nth[*], $
           dBR_geiGrid_sh_nth[*], dBTheta_geiGrid_sh_nth[*], dBPhi_geiGrid_sh_nth[*], $
           dbx_geiGrid_sh_nth, dby_geiGrid_sh_nth, dbz_geiGrid_sh_nth

; unshift the fitted vectors
  ampere_shift_vec, geiGridX_shifted_nth[*], geiGridY_shifted_nth[*], geiGridZ_shifted_nth[*], $  ; input GEI locations
          dbx_geiGrid_sh_nth[*], dby_geiGrid_sh_nth[*], dbz_geiGrid_sh_nth[*], $                  ; input GEI vectors
          dbx_geiGrid, dby_geiGrid, dbz_geiGrid, $                    ; unshifted GEI vectors in XYZ
          dBR_geiGrid_nth, dBTheta_geiGrid_nth, dBPhi_geiGrid_nth, $  ; unshifted spherical dB's
          rev_rot_mat                                                 ; rotation matrix

  geiGridX_shifted_nth=!null
  geiGridY_shifted_nth=!null
  geiGridZ_shifted_nth=!null

  dbx_geiGrid_sh_nth = !null
  dby_geiGrid_sh_nth = !null
  dbz_geiGrid_sh_nth = !null

  dbx_geiGrid=!null
  dby_geiGrid=!null
  dbz_geiGrid=!null

  geiGrid_coLat_rad_shifted_nth = !null
  geiGrid_Lon_rad_shifted_nth = !null
  dBR_geiGrid_sh_nth = !null
  dBTheta_geiGrid_sh_nth = !null
  dBPhi_geiGrid_sh_nth = !null

; Conv uniform grid data to AACGM - North
  geiVec_to_aacgmVec, $
    geiGrid_R_km_nth[*], geiGrid_coLat_rad_nth[*], geiGrid_lon_rad_nth[*], $  ; input GEI coords
    dBTheta_geiGrid_nth[*], dbPhi_geiGrid_nth[*], $                           ; input GEI vecs
    dBTheta_aacgm_gth = dBTheta_aacgmGrid_gth_nth, $
    dBPhi_aacgm_gth = dBPhi_aacgmGrid_gth_nth, $
    dBTheta_aacgm_gph = dBTheta_aacgmGrid_gph_nth, $
    dBPhi_aacgm_gph = dBPhi_aacgmGrid_gph_nth

; Original input data to AACGM - for IRD file output
  geiVec_to_aacgmVec, $
    dataOriginal[*].gei_R_km, dataOriginal[*].gei_coLat_rad, dataOriginal[*].gei_lon_rad, $
    dataOriginal[*].bTheta_gei, dataOriginal[*].bPhi_gei, $
    dBTheta_aacgm_gth = dBTheta_aacgmDataOriginal_gth, $
    dBPhi_aacgm_gth = dBPhi_aacgmDataOriginal_gth, $
    dBTheta_aacgm_gph = dBTheta_aacgmDataOriginal_gph, $
    dBPhi_aacgm_gph = dBPhi_aacgmDataOriginal_gph, $
    aacgm_lat_deg = aacgmDataOriginal_lat_deg, $
    aacgm_lon_deg = aacgmDataOriginal_lon_deg

; Calc variables used in output files
  geopack_epoch,sepoch,yr0,mo0,day0,hh0,mm0,ss0,/breakdown
  geopack_epoch,eepoch,yr1,mo1,day1,hh1,mm1,ss1,/breakdown
  d_epoch=(eepoch-sepoch)/1000.
  yrstr=string(yr0,format='(i4.4)')
  daystr=string(mo0,format='(i2.2)')+string(day0,format='(i2.2)')
  t0str=string(hh0,format='(i2.2)')+string(mm0,format='(i2.2)')+string(ss0,format='(i2.2)')
  t1str=string(hh1,format='(i2.2)')+string(mm1,format='(i2.2)')+string(ss1,format='(i2.2)')
  amp_fits_subpath=amp_fits_path+yrstr+path_sep()+daystr+path_sep()
  exist=file_test(amp_fits_subpath,/directory)
  if not exist then file_mkdir,amp_fits_subpath

; Nth hemisphere data
  hemstr='north'
  resol_deg=SpacRes
  grdfname=amp_fits_subpath+'sphfit_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.grd'
  write_grd_file,grdfname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, cap_sz_deg, nLonGrid, $
    aacgmGrid_coLat_deg_nth[*], $
    ((aacgmGrid_Lon_deg_nth[*]/15.0+mltshift) mod 24), $
    -dBTheta_aacgmGrid_gth_nth[*], dBPhi_aacgmGrid_gth_nth[*], $
    -dBTheta_aacgmGrid_gph_nth[*], dBPhi_aacgmGrid_gph_nth[*],$
     jParAACGM_nth[*]

; Write ird data file
  irdfname=amp_fits_subpath+'sphfit_'+yrstr+daystr+'_'+t0str+'_'+t1str+'.ird'
  write_ird_file,irdfname, yrstr, daystr, t0str, d_epoch,$
    90.0-aacgmDataOriginal_lat_deg[*], $
    ((aacgmDataOriginal_lon_deg[*]/15.0+mltshift) mod 24), $
    -dBTheta_aacgmDataOriginal_gth[*], dBPhi_aacgmDataOriginal_gth[*], $
    -dBTheta_aacgmDataOriginal_gph[*], dBPhi_aacgmDataOriginal_gph[*], $
    dataOriginal[*].ipln, $
    dataOriginal[*].isat, $
    dataOriginal[*].qual, $
    dataOriginal[*].splice

; Write coeff file here
  coffname=amp_fits_subpath+'sphfit_'+yrstr+daystr+'_'+t0str+'_'+t1str+'.cof'
  write_cof_file,coffname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, outKvalues,outMvalues, coeffs_

; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
	if debug then begin
  ; Plot jPar in various coord systems
  ; ----------------------------------

		winNo = 3
		winNo++
		window, winNo, xSize = 680, ySize = 360, title='North - FAC'
		!p.multi = [0,2,1]
		!p.charSize = 1.0

    south=0
    pole = 90
    capLimit = 50.             ; coLat size

    plot_fac_cw, jParAACGM_nth[*], pole, capLimit, $
       geiGRid_coLat_rad_nth[*]*!radeg, geiGrid_lon_rad_nth[*]*!radeg, mn_fac, mx_fac, path, $
       title = 'North: jPar on GEI grid', $
       south = south

		plot_fac_cw, jParAACGM_nth[*], pole, capLimit, $
				aacgmGrid_coLat_deg_nth[*], aacgmGrid_lon_deg_nth[*]+mltShift*15, mn_fac, mx_fac, path, $
				title = 'North: jPar on grid [AACGM-MLT]', $
				south = south

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_polfacfull_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for North Polar FAC written to ', pngfname
    end

  ; Plot the dB vectors in AACGM-MLT - North
  ; --------------------------------
    winNo ++
    window, winNo, xSize = 680, ySize = 360, title='Fitted dB Vectors - North'
    !p.multi = [0,2,1]
    !p.charSize = 1.0
    loadct,0,/sil

    south=0
    plt_dat, 90.0-aacgmDataOriginal_lat_deg[*], ((aacgmDataOriginal_lon_deg[*]/15.0+mltShift) mod 24)*15.0, $
      n_elements(aacgmDataOriginal_lat_deg[*] ), $
      -dBTheta_aacgmDataOriginal_gth, dBPhi_aacgmDataOriginal_gph, $
      south, title = 'North Input dataOriginal: AACGM-MLT', $
      capSize = cap_sz_deg

    plt_dat, aacgmGrid_coLat_deg_nth[*], ((aacgmGrid_lon_deg_nth[*]/15.0+mltShift) mod 24)*15.0, $
            n_elements(aacgmGrid_coLat_deg_nth[*] ), $
           -dBTheta_aacgmGrid_gth_nth, dBPhi_aacgmGrid_gph_nth, south,$
        title = 'North Fitted data: AACGM-MLT', $
        capSize = cap_sz_deg

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_dBvec_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for North Vector dB written to ', pngfname
    endif
  end
; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
; ############ SOUTH HEMISPHERE AACGM CAP ###############################
;
; Create uniform grids in AACGM and GEI
; -------------------------------------
; 1. Calculate Uniform AACGM grid in Sth hemisphere cap [AACGM not defined at equatorial lats]
; 2. Conv to geog -> then to GEI coord

  south=1
  aacgm_grid, cap_sz_deg, south, $             ; get north hemisphere cap grid
       aacgmGrid_R_km = aacgmGrid_R_km_sth, $
       aacgmGrid_coLat_deg = aacgmGrid_coLat_deg_sth, $
       aacgmGrid_lon_deg = aacgmGrid_lon_deg_sth, $
       geiGrid_R_km = geiGrid_R_km_sth, $
       geiGrid_coLat_rad = geiGrid_coLat_rad_sth, $
       geiGrid_lon_rad = geiGrid_lon_rad_sth, $
       xGrid_GEI=xGrid_GEI_sth, $
       yGrid_GEI=yGrid_GEI_sth, $
       zGrid_GEI=zGrid_GEI_sth, $
       nLat = cap_sz_deg, nLon = nLonGrid, $
       year = year, yrSec = avgYrSec, $
       mltShift = mltShift, $
       epoch = avgEpoch

; 3. Calc shifted GEI from uniform AACGM grid

  ampere_shift_coord, xGrid_GEI_sth[*], yGrid_GEI_sth, zGrid_GEI_sth[*], $     ; input GEI coords
      geiGridX_shifted_sth, geiGridY_shifted_sth, geiGridZ_shifted_sth, $      ; output shifted GEI coords
      geiGrid_rIrid_km_shifted_sth, geiGrid_coLat_deg_shifted_sth, geiGrid_lon_deg_shifted_sth, $  ; spherical coords of shifted coords
      rot_mat          ; rotation matrix

  xGrid_GEI_sth=!null
  yGrid_GEI_sth=!null
  zGrid_GEI_sth=!null

  geiGrid_coLat_rad_shifted_sth = geiGrid_coLat_deg_shifted_sth[*]*!dtor
  geiGrid_lon_rad_shifted_sth = geiGrid_lon_deg_shifted_sth[*]*!dtor

  ; Generate basis fns at regular grid
  ; ----------------------------------

  print,'Generating Basis Set over Southern Uniform Grid....'

  ampere_setupSHFns_full, 1.0, geiGrid_coLat_rad_shifted_sth[*], geiGrid_lon_rad_shifted_sth[*], $
      kMax, mMax, $
      brBFnArr = brBFns_grid_sth, $
      bThBFnArr = dYkmDThBFns_grid_sth, $
      bPhBFnArr = dYkmDPhBFns_grid_sth, $
      YBFnArr = YkmBFns_grid_sth, $
      mArr = outMValues_grid_sth, $
      nArr = outNkValues_grid_sth, $
      kArr = outKValues_grid_sth

  Print,'Calculating FAC over Southern Uniform Grid....'

  bFuncs_grid_sth  = temporary($
                 [[dYkmDThBfns_grid_sth, -dYkmDPhBfns_grid_sth ], $
                  [dYkmDPhBfns_grid_sth,  dYkmDthBfns_grid_sth ]])

  jParAACGM_sth = (YkmBFns_grid_sth ## $
         (coeffs_[n_elements(YkmBFns_grid_sth[*,0]):*]*(-outNkValues_grid_sth*(outNkValues_grid_sth+1.0)))) $
                         /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}

  jParAACGM_sth = reform(jParAACGM_sth, cap_sz_deg, nLonGrid)

  print,'Min:Max jPar - South = ',min(jParAACGM_sth),max(jParAACGM_sth)

  ; Reconstruct the dB vector from the fit
  ; --------------------------------------

  recon_dB_geiGrid  = bFuncs_grid_sth ## coeffs_
  bFuncs_grid_sth=!null      ; release this memory

  dBTheta_geiGrid_sh_sth = recon_dB_geiGrid[0:cap_sz_deg*nLonGrid-1]
  dBPhi_geiGrid_sh_sth   = recon_dB_geiGrid[cap_sz_deg*nLonGrid:*]
  dBR_geiGrid_sh_sth   = dBPhi_geiGrid_sh_sth * 0.0d0   ; set to zero for now
;  dBR_geiGrid_sh_sth   = brBFns_grid_sth ## coeffs_r   ; dBR
;  brBFns_grid_sth=!null                      ; free some memory
  recon_dB_geiGrid = !null

  ; Un-Shift the dB Vectors
  ; --------------------------------

  rev_rot_mat = transpose(rot_mat)

;  conv fitted dB vectors from shifted spherical to shifted XYZ components
  geopack_bspcar,geiGrid_coLat_rad_shifted_sth, geiGrid_lon_rad_shifted_sth, $
           dBR_geiGrid_sh_sth, dBTheta_geiGrid_sh_sth, dBPhi_geiGrid_sh_sth, $
           dbx_geiGrid_sh_sth, dby_geiGrid_sh_sth, dbz_geiGrid_sh_sth

; unshift the fitted vectors
  ampere_shift_vec, geiGridX_shifted_sth[*], geiGridY_shifted_sth[*], geiGridZ_shifted_sth[*], $  ; input GEI locations
          dbx_geiGrid_sh_sth[*], dby_geiGrid_sh_sth[*], dbz_geiGrid_sh_sth[*], $                  ; input GEI vectors
          dbx_geiGrid, dby_geiGrid, dbz_geiGrid, $        ; unshifted GEI vectors in XYZ
          dBR_geiGrid_sth, dBTheta_geiGrid_sth, dBPhi_geiGrid_sth, $  ; unshifted spherical dB's
          rev_rot_mat                   ; rotation matrix

  dbx_geiGrid_sh_sth=!null
  dby_geiGrid_sh_sth=!null
  dbz_geiGrid_sh_sth=!null

  dbx_geiGrid=!null
  dby_geiGrid=!null
  dbz_geiGrid=!null

  geiGridX_shifted_sth=!null
  geiGridY_shifted_sth=!null
  geiGridZ_shifted_sth=!null

  geiGrid_coLat_rad_shifted_sth = !null
  geiGrid_Lon_rad_shifted_sth = !null
  dBR_geiGrid_sh_sth = !null
  dBTheta_geiGrid_sh_sth = !null
  dBPhi_geiGrid_sh_sth = !null

; Conv uniform grid data to AACGM - South
  geiVec_to_aacgmVec, $
    geiGrid_R_km_sth[*], geiGrid_coLat_rad_sth[*], geiGrid_lon_rad_sth[*], $   ; input GEI locations
    dBTheta_geiGrid_sth[*], dbPhi_geiGrid_sth[*], $                            ; input GEI vecs
    dBTheta_aacgm_gth = dBTheta_aacgmGrid_gth_sth, $
    dBPhi_aacgm_gth = dBPhi_aacgmGrid_gth_sth, $
    dBTheta_aacgm_gph = dBTheta_aacgmGrid_gph_sth, $
    dBPhi_aacgm_gph = dBPhi_aacgmGrid_gph_sth

; Sth hemisphere data
; Write GRD data file
  hemstr='south'
  grdfname=amp_fits_subpath+'sphfit_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.grd'
  write_grd_file,grdfname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, cap_sz_deg, nLonGrid, $
    aacgmGrid_coLat_deg_sth[*], $
    ((aacgmGrid_Lon_deg_sth[*]/15.0+mltshift) mod 24), $
    -dBTheta_aacgmGrid_gth_sth[*], dBPhi_aacgmGrid_gth_sth[*], $
    -dBTheta_aacgmGrid_gph_sth[*], dBPhi_aacgmGrid_gph_sth[*], $
     jParAACGM_sth[*]

; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
  if debug then begin
  ; Plot jPar in various coord systems
  ; ----------------------------------

    winNo++
    window, winNo, xSize = 680, ySize = 360, title='South - FAC'
    !p.multi = [0,2,1]
    !p.charSize = 1.0

    south=1
    pole = -90
    capLimit = 180.-capLimit             ; coLat

    plot_fac_cw, jParAACGM_sth[*], pole, capLimit, $
       geiGRid_coLat_rad_sth[*]*!radeg, geiGrid_lon_rad_sth[*]*!radeg, mn_fac, mx_fac, path, $
       title = 'South: jPar on GEI grid', $
       south = south

    plot_fac_cw, jParAACGM_sth[*], pole, capLimit, $
        aacgmGrid_coLat_deg_sth[*], aacgmGrid_lon_deg_sth[*]+mltShift*15, mn_fac, mx_fac, path, $
        title = 'South: jPar on grid [AACGM-MLT]', $
        south = south

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_polfacfull_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for South Polar FAC written to ', pngfname
    end

  ; Plot the dB vectors in AACGM-MLT - South
  ; --------------------------------
    winNo ++
    window, winNo, xSize = 680, ySize = 360, title='Fitted dB Vectors - South'
    !p.multi = [0,2,1]
    !p.charSize = 1.0
    loadct,0,/sil

;    cap_sz_deg=180.-cap_sz_deg          ; cap size
    south=1
    plt_dat, 90.0-aacgmDataOriginal_lat_deg[*], ((aacgmDataOriginal_lon_deg[*]/15.0+mltShift) mod 24)*15.0, $
      n_elements(aacgmDataOriginal_lat_deg[*] ), $
      -dBTheta_aacgmDataOriginal_gth, dBPhi_aacgmDataOriginal_gph, $
      south, title = 'South Input dataOriginal: AACGM-MLT', $
      capSize = cap_sz_deg

    plt_dat, aacgmGrid_coLat_deg_sth[*], ((aacgmGrid_lon_deg_sth[*]/15.0+mltShift) mod 24)*15.0, $
            n_elements(aacgmGrid_coLat_deg_sth[*] ), $
           -dBTheta_aacgmGrid_gth_sth, dBPhi_aacgmGrid_gph_sth, south,$
        title = 'South Fitted data: AACGM-MLT', $
        capSize = cap_sz_deg

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_dBvec_'+yrstr+daystr+'_'+t0str+'_'+t1str+'_'+hemstr+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for North Vector dB written to ', pngfname
    endif
    !p.charSize = 1.0
    !p.multi = 0

  end
; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
;
; ############## WHOLE SPHERE ###########################
; Whole sphere, GEI, GEO dB and FAC section
; Generate whole sphere grid

  Re_km = 6371.0d0
  R_km = 780.0d0 + Re_km
  nLatGrid=90                ; 2 deg in lat, change to 180 for 1 deg spacing
  nLonGrid=24
  geiuni_coLat_deg = (dindgen(nLatGrid)*2+1)
  geiuni_lon_deg = (dindgen(nLonGrid))*(360.0-360.0/nLonGrid)/(nLonGrid-1)
  geiuni_lon_deg = geiuni_lon_deg mod 360.0

  geiuni_coLat_deg = rebin(geiuni_coLat_deg, nLatGrid, nLonGrid )         ; 180 by 24
  geiuni_lon_deg = transpose (rebin(geiuni_lon_deg, nLonGrid, nLatGrid))  ;
  geiuni_R_km = geiuni_lon_deg * 0.0d0 + R_km

  ; Convert from GEI_sph to GEI_xyz locations
  ; ------------------------
;  geopack_epoch,AvgEpoch,yrA,moA,dayA,hhA,mmA,ssA,/breakdown
;  print,yrA,moA,dayA,hhA,mmA,ssA
;stop

;  geopack_recalc,yrA,moA,dayA,hhA,mmA,ssA,/DATE                ; for correct UT
  geoPack_sphCar, geiuni_R_km[*], geiuni_coLat_deg[*], geiuni_lon_deg[*], $
    xGrid_GEI, yGrid_GEI, zGrid_GEI, /to_rect, /degree

; Conv from GEI to GEO for output file
  geopack_conv_coord, xGrid_GEI, yGrid_GEI, zGrid_GEI, $
                      xGrid_GEO, yGrid_GEO, zGrid_GEO, $
                      /from_gei, /to_geo

  geoPack_sphCar, xGrid_GEO, yGrid_GEO, zGrid_GEO, $
                  geouni_R_km, geouni_coLat_deg, geouni_lon_deg, $
                  /to_sphere, /degree

; Get this grid in shifted coords
  ampere_shift_coord, xGrid_GEI[*], yGrid_GEI[*], zGrid_GEI[*], $     ; input GEI coords
      geiuniX_shifted, geiuniY_shifted, geiuniZ_shifted, $            ; output shifted GEI coords
      geiuni_rIrid_km_shifted, geiuni_coLat_deg_shifted, geiuni_lon_deg_shifted, $  ; spherical coords of shifted coords
      rot_mat          ; rotation matrix

  geiuni_coLat_rad_shifted = geiuni_coLat_deg_shifted[*]*!dtor
  geiuni_lon_rad_shifted = geiuni_lon_deg_shifted[*]*!dtor

  ; Generate basis fns at full sphere, regular GEI grid
  ; ----------------------------------

  print,'Generating Basis Set over Full Sphere Uniform Grid....'

  ampere_setupSHFns_full, 1.0, geiuni_coLat_rad_shifted[*], geiuni_lon_rad_shifted[*], $
      kMax, mMax, $
      brBFnArr = brBFns_geiuni, $
      bThBFnArr = dYkmDThBFns_geiuni, $
      bPhBFnArr = dYkmDPhBFns_geiuni, $
      YBFnArr = YkmBFns_geiuni, $
      mArr = outMValues_geiuni, $
      nArr = outNkValues_geiuni, $
      kArr = outKValues_geiuni

  bFuncs_geiuni  = temporary($
            [[dYkmDThBfns_geiuni, -dYkmDPhBfns_geiuni ], $
             [dYkmDPhBfns_geiuni,  dYkmDthBfns_geiuni ]])

  jPargeiuni = (YkmBFns_geiuni ## $
      (coeffs_[n_elements(YkmBFns_geiuni[*,0]):*]*(-outNkValues_geiuni*(outNkValues_geiuni+1.0))))$
                         /(u0*rIrid_m)*1.0d-9*1.0d6        ; uAm^{-2}

  jPargeiuni = reform(jPargeiuni, nLatGrid, nLonGrid)

; Reconstruct dB's over full sphere
  dB_geiuni  = bFuncs_geiuni ## coeffs_
  bFuncs_geiuni = !null

  dBTheta_geiuni = dB_geiuni[0:nLatGrid*nLonGrid-1]
  dBPhi_geiuni   = dB_geiuni[nLatGrid*nLonGrid:*]
  dBR_geiuni     = dBPhi_geiuni * 0.0d0        ; set to zero for now
;  dBR_geiuni     = brBFns_geiuni ## coeffs_r   ; dBR
;  brBFns_geiuni = !null                      ; free memory
  dB_geiuni = !null
  dB_mg=sqrt(dBTheta_geiuni[*]^2+dBPhi_geiuni[*]^2)

; Write to file
; ----------------------
  wspfname=amp_fits_subpath+'sphfit_'+yrstr+daystr+'_'+t0str+'_'+t1str+'.wsp'
  write_wsp_file,wspfname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, nLatGrid, nLonGrid, $
    geiuni_coLat_deg[*], ((geiuni_Lon_deg[*]/15.0) mod 24), ((geouni_Lon_deg[*]/15.0) mod 24), $
   -dBTheta_geiuni[*], dBPhi_geiuni[*], $
    jPargeiuni[*]

; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d
  if debug then begin

    winNo++
    window, winNo, xSize = 700, ySize = 550, title='GEO FAC'
    !p.multi = 0
    device,decomposed=0
    loadct, 13, file = path + 'davect.tbl', /silent
    lon_deg=geouni_lon_deg[*]
    idx=where(lon_deg gt 180.)
    lon_deg[idx]=lon_deg[idx]-360.
    plot_fac_moldw,jPargeiuni, 90.-geiuni_coLat_deg, lon_deg, $
        mn_fac, mx_fac, title = 'GEO FAC'

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_geofac_'+yrstr+daystr+'_'+t0str+'_'+t1str+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for GEO FAC written to ', pngfname
    end

    mn_dB=30.
    mx_dB=700.
    winNo++

    window, winNo, xSize = 700, ySize = 550, title='GEO |dB|'
    !p.multi = 0
    loadct,22,/silent
    tvlct,r,g,b,/get
    r[255]=0 & g[255]=0 & b[255]=0
    tvlct,r,g,b
    plot_fac_moldw,dB_mg, 90.-geiuni_coLat_deg, lon_deg, $
        mn_dB, mx_dB, title = '|dB|'

    if plt_png eq 1 then begin
      pngfname=amp_fits_subpath+'sphfit_geodbmag_'+yrstr+daystr+'_'+t0str+'_'+t1str+'.png'
      image = TVRD(0,0,!d.X_size,!d.Y_size,true=1)
      write_PNG,pngfname,image,r,g,b
      Print,'PNG for GEO |dB| written to ', pngfname
    end

	endif     ; if debug FAC maps
; d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d-d


  !p.multi = 0

	;@big_nsf_plot

end

