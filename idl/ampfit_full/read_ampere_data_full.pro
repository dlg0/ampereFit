; +
; Reads AMPERE data file (input data)
;
; David L Green & Colin L Waters
; Centre for Space Physics
; University of NEwcastle
; Dec 2009
;
; Comments:
; data selection -> full sphere version [Sept 2010 by CLW]
; modified may 2010 - sort data by orbit plane and sequentially along track (for diagnostics)
;                   - tweak ampdata shift routine [CLW]
;
; DLG 4-Jul-10
;	There are 2 data strucutres, the un-shifted and shifted.
;	Both are always created and then selected later.
;
; CLW 7 Sept 2010
; Added netcdf read capability
; After data read, set large arrays to !null to free up memory [CLW]
; Added data Ghost routine for longitude
;

pro read_ampere_data_full, sHr, eHr, $
	datFileName = datFileName, $
	dataOriginal = dataOriginal, $
  dataShifted = dataShifted, $
	year = year, month = month, day = day, $
	avgYrSec = yrSecAvg, $
	avgEpoch = avgEpoch, $
	rot_mat=rot_mat, $
	debug = debug, $
	extend = extend, $
	SpacRes = SpacRes, $
	status=status

	forward_function cnvTime

; Trap for errors
  if debug eq 0 then begin
    status=0
    catch,error_status
    if (error_status ne 0) then begin
		  status=1
      Print,'Error in read_ampere_data_full - returning'
      catch,/cancel
      return
    endif
  endif

; NetCDF read section
  read_ampere_ncdf,datFileName,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice

  dateStr=''
  reads, strmid(file_baseName ( datFileName ),0,8), year, month, day, format='(i4,i2,i2)'

  if keyword_set(extend) then begin
    geopack_epoch,epoch,year,month,day,/compute
    epoch=epoch+86400.0d3
    geopack_epoch,epoch,extyear,extmonth,extday,/breakdown
    extdatfilename=file_basename(datfilename)
    strput,extdatfilename,string(extyear,format='(i4.4)'),0
    strput,extdatfilename,string(extmonth,format='(i2.2)'),4
    strput,extdatfilename,string(extday,format='(i2.2)'),6
    extdatfilename=file_dirname(datfilename)+path_sep()+extdatfilename

    if file_test(extdatfilename) then begin
      t_x_axis_frac_hour=x_axis_frac_hour
      t_pseudosvnum_total=pseudosvnum_total
      t_plane_number_total=plane_number_total
      t_pos_eci_total=pos_eci_total
      t_b_eci=b_eci
      t_pseudo_sv_quality=pseudo_sv_quality
      t_data_splice=data_splice

      read_ampere_ncdf,extdatfilename,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice

      x_axis_frac_hour=[t_x_axis_frac_hour,x_axis_frac_hour+24]
      pseudosvnum_total=[t_pseudosvnum_total,pseudosvnum_total]
      plane_number_total=[t_plane_number_total,plane_number_total]
      pos_eci_total=[[t_pos_eci_total],[pos_eci_total]]
      b_eci=[[t_b_eci],[b_eci]]
      pseudo_sv_quality=[t_pseudo_sv_quality,pseudo_sv_quality]
      data_splice=[t_data_splice,data_splice]

    endif else begin
      status=1
      return
    endelse
  endif

	;
	; Sav file read section
	; ---------------------------

;	restore, savFileName

	; Variable in save file
	;
	;   x_axis_frac_hour : UTC in fractional hours
	; plane_number_total : integer number 0-5 for the 6 orbit planes.
	;      pos_ECI_total : ECI position in meters, dimensions [{X,Y,Z},n]
	;              B_ECI : the same dimension as pos_eci but has delta B in nT
	; ECI = GEI


	; Create data structure
	; ---------------------

 	data = { $
		utc: 0d0, $                ; UT time in dec hours
		isat: 0, $                 ; coded SV number (for Haje)
		iPln: 0, $                 ; orbit track number (0->5)
		qual: 0d0, $               ; data quality from Lars
		splice : 0, $              ; flag for spliced data - where missing data estimated
		typ : 0, $                 ; data type (for CLW)
		px: 0d0, py: 0d0, pz: 0d0, $
		dbx: 0d0, dby: 0d0, dbz: 0d0, $
		GEI_R_km : 0d0, GEI_coLat_deg: 0d0, GEI_lon_deg: 0d0, $
		GEI_coLat_rad: 0d0, GEI_lon_rad: 0d0, $
		br_GEI: 0d0, bTheta_GEI: 0d0, bPhi_GEI: 0d0 }

  ; Select subset of data, based on time interval
  ; ---------------------

;  nPts = n_elements( x_axis_frac_hour )

  @amp_fit_constants

  iiSubSet = where ( $
    x_axis_frac_hour ge sHr and $
    x_axis_frac_hour le eHr ,iiSubSetCnt )

  ; Fill data structure
  ; -------------------

  dataOriginal = temporary ( replicate ( data, iiSubSetCnt ) )

  dataOriginal.utc = x_axis_frac_hour[iiSubSet]
  dataOriginal.iSat = pseudoSVNum_total[iiSubSet]
  dataOriginal.iPln = plane_number_total[iiSubSet]
  dataOriginal.qual = pseudo_sv_quality[iiSubSet]
  dataOriginal.splice = data_splice[iiSubSet]

  ; GEI XYZ position in KM

  dataOriginal.px = (pos_ECI_total[0,*])[iiSubSet]*1d-3
  dataOriginal.py = (pos_ECI_total[1,*])[iiSubSet]*1d-3
  dataOriginal.pz = (pos_ECI_total[2,*])[iiSubSet]*1d-3

  ; GEI XYZ db vector

  dataOriginal.dbx = (B_ECI[0,*])[iiSubSet]
  dataOriginal.dby = (B_ECI[1,*])[iiSubSet]
  dataOriginal.dbz = (B_ECI[2,*])[iiSubSet]

  ; Get spherical coords of the GEI XYZ locations
  ; ---------------------------------------------

  geopack_sphcar, dataOriginal.px, dataOriginal.py, dataOriginal.pz, $
                GEI_R_km, GEI_coLat_deg, GEI_lon_deg, $
                /to_sphere, /degree

  ; coLat is 0->90 for Nth and 90->180 for Sth (deg)

  dataOriginal.GEI_R_km = GEI_R_km
  dataOriginal.GEI_coLat_deg = GEI_coLat_deg
  dataOriginal.GEI_lon_deg = GEI_lon_deg

  dataOriginal.GEI_coLat_rad = dataOriginal.GEI_coLat_deg * !dtor
  dataOriginal.GEI_lon_rad = dataOriginal.GEI_lon_deg * !dtor

  ; Rotate XYZ GEI db to spherical GEI
  ; ----------------------------------

  geopack_bcarsp, dataOriginal.px,  dataOriginal.py,  dataOriginal.pz, $
                  dataOriginal.dbx, dataOriginal.dby, dataOriginal.dbz, $
                  br_GEI, bTheta_GEI, bPhi_GEI

  dataOriginal.br_GEI = br_GEI
  dataOriginal.bTheta_GEI = bTheta_GEI
  dataOriginal.bPhi_GEI = bPhi_GEI

; Free up some memory
  x_axis_frac_hour = !null       ; goes to utc
  pseudosvnum_total = !null      ; goes to iSat
  plane_number_total = !null     ; goes to iPln
  pos_eci_total = !null          ; px, py, pz
  b_eci = !null                  ; dbx, dby, dbz
  pseudo_sv_quality = !null      ; data quality from Lars
  data_splice = !null            ; splice data flag

; Satellite track labels, examine longitudes - CLW Oct 2010
  For Tr=0,5 do begin
    idx=where(dataOriginal.iPln eq Tr)  ; index of data from Lars' Track assignment data
    np=n_elements(idx)
    bnsz=5.
    tr_hist=histogram(dataOriginal[idx].GEI_Lon_Deg,binsize=bnsz)   ; get Lons histogram (should be double peaked, extra track s gives 3 peaks
    amp_trk_lab,np,tr_hist,bnsz,Ln_trk
    if Ln_trk(0) gt -1 then begin
      if n_elements(Ln_trk) gt 2 Then Print,'WARNING: Orbit Track mis-label detected on Track ',Tr
    end else print,'WARNING : no peaks detected in Longitude Histogram Track identification routine.'
  end
; If Trk_Lons[Tr,2] ne -1 then we have track number mis-assigngment
; Also, abs(Trk_Lons[Tr,1] - Trk_Lons[Tr,0]) should be ~ 180 deg
;
; Shift the data to be more centered near the track intersection point.
; -------------------------------------------------------

    south=0
    ampere_get_rotmat_full, dataOriginal, rot_mat, south;,/diag   ; does not require data to be in any particular order
    ampere_shift_data, dataOriginal, dataShifted, rot_mat        ; shifts coords and vecs

;  endif else begin
;    dataShifted = dataOriginal          ; no rotation option
;    rot_mat = [[1.d0,0.d0,0.d0],$
;               [0.d0,1.d0,0.d0],$
;               [0.d0,0.d0,1.d0]]
;  endelse

 ;
 ; Sort data along each orbit track - to get ready for longitude ghost data
 ; --------------------------------------------------------------

  sortII = indGen ( iiSubSetCnt )
  cnt = 0
  stLon_a=fltarr(6)

	for trackNo = 0, 5 do begin
; get data from each track and Nth hemis 1/2
		iiThisTrack_nth = where( (dataShifted.iPln eq trackNo) and $
		                         (dataShifted.GEI_colat_deg le 90.), thisTrackCnt_nth )    ; select orbit plane and Nth hemis
;*********************************************************************************************
;; Check locations using this test plot command
;plot,dataShifted[iiThisTrack_nth].GEI_coLat_deg,dataShifted[iiThisTrack_nth].GEI_Lon_deg,psym=3,yrange=[-20,360],$
;     xtitle='GEI coLat',ytitle='GEI Longitude'
;print,'Num points available = ',n_elements(dataShifted[iiThisTrack_nth].GEI_coLat_deg)
; aa=10
; xyouts,dataShifted[iithistrack_nth[aa]].gei_colat_deg,dataShifted[iithistrack_nth[aa]].gei_lon_deg,strtrim(string(dataShifted[iithistrack_nth[aa]].isat),2)
; stop
;*********************************************************************************************
    lon_sel = where( (dataShifted.GEI_Lon_deg)[iiThisTrack_nth] le 180. )   ; ensure start in same 1/2 longitude space
    iiStartHere_nth = (where( dataShifted[iiThisTrack_nth[lon_sel]].GEI_coLat_deg eq max(dataShifted[iiThisTrack_nth[lon_sel]].GEI_coLat_deg) ))[0]
    strt_clat=dataShifted[iiThisTrack_nth[lon_sel[iiStartHere_nth]]].GEI_coLat_deg
    strt_lon =dataShifted[iiThisTrack_nth[lon_sel[iiStartHere_nth]]].GEI_Lon_deg
    stLon_a[TrackNo]=strt_lon
    gcDistance_nth=gc_dist( 90.-strt_clat, 90.-dataShifted[iiThisTrack_nth].GEI_coLat_deg, $
                            strt_lon, dataShifted[iiThisTrack_nth].GEI_Lon_deg)
    iiThisTrackSorted_nth = sort ( gcDistance_nth )
    sortII[cnt:cnt+thisTrackCnt_nth-1] = iiThisTrack_nth[iiThisTrackSorted_nth]
    cnt += thisTrackCnt_nth

; get data from each track and Sth hemis 1/2
    iiThisTrack_sth = where ( (dataShifted.iPln eq trackNo) and $
                 (dataShifted.GEI_colat_deg gt 90.), thisTrackCnt_sth )
; use same start point as nth - since very close to equator
    gcDistance_sth=gc_dist( 90.-strt_clat, (90.-dataShifted.GEI_coLat_deg)[iiThisTrack_sth], $
                            strt_lon, (dataShifted.GEI_Lon_deg)[iiThisTrack_sth])
    iiThisTrackSorted_sth = sort ( gcDistance_sth )
    sortII[cnt:cnt+thisTrackCnt_sth-1] = iiThisTrack_sth[iiThisTrackSorted_sth]
    cnt += thisTrackCnt_sth

	endfor

  dataShifted = dataShifted[sortII]
  o_np=n_elements(sortII)

; check for points too close to the pole
  vd_i=where( (dataShifted.GEI_coLat_deg gt 1.) and (dataShifted.GEI_coLat_deg lt 179.), n_np, complement=nvd_i)
  if (n_np lt o_np) then begin
    Print,'Number of data points close to the pole removed : ',o_np-n_np
    Print,'coLat of these : ',dataShifted[nvd_i].GEI_coLat_deg
    tmp_dat=dataShifted   ; make a copy
    dataShifted = tmp_dat[vd_i]
    tmp_dat = !null
  end
  dataOriginal = dataOriginal[sortII]
  trk_order=sort(stLon_a)              ; sorted by Lon

; Latitude Ghost routine - Lars to take care of
; Longitude ghost routine - CLW
; ----------------------------

; Find strays and label [see CLW for details]
; Data is sorted by track, latitude (start at equator in lon sector < 180 deg with adjacent track order known
; Added check for Lat gaps - CLW, computes great circle eqn for each track (using vectors)
  For Tr_ii = 0, 5 do begin
; Check Nth hemis
    ii_nth=where( (dataShifted.iPln eq Tr_ii) and (dataShifted.GEI_colat_deg le 90.), ii_nth_cnt)
    dataShifted[ii_nth].typ = 0              ; initialise data with 'normal' tag
    u = [dataShifted[ii_nth[0]].px, dataShifted[ii_nth[0]].py, dataShifted[ii_nth[0]].pz]
    u_mag = sqrt(u(0)^2 + u(1)^2 + u(2)^2)
    u=u/u_mag
; get second vector around midway, could double check with longitude histogram values from above
    sel_lat=where( (dataShifted[ii_nth].GEI_coLat_deg gt 40.) and $
                   (dataShifted[ii_nth].GEI_coLat_deg lt 50.) )  ; assumes strays are interleaved or poleward of this
                                                                 ; selects both longitude 1/2 spaces
    mn=min(abs(dataShifted[ii_nth[sel_lat]].GEI_Lon_deg-dataShifted[ii_nth[0]].GEI_Lon_deg),loc)   ; avoid strays and get correct longitude
    v = [dataShifted[ii_nth[sel_lat[loc]]].px, dataShifted[ii_nth[sel_lat[loc]]].py, dataShifted[ii_nth[sel_lat[loc]]].pz]
    v_mag = sqrt(v(0)^2 + v(1)^2 + v(2)^2)
    v=v/v_mag

    cross_p,u,v,n
    n_mag = sqrt(n(0)^2 + n(1)^2 + n(2)^2)
    n=n/n_mag

    cross_p,n,u,res

    r=sqrt(dataShifted[ii_nth].px^2 + dataShifted[ii_nth].py^2 + dataShifted[ii_nth].pz^2)
    t=dblarr(ii_nth_cnt)
    t(0)=0.0d0          ; IDL does a -NaN here - weird.
    for cc=1,ii_nth_cnt-1 do begin
      dp = u(0)*dataShifted[ii_nth[cc]].px + u(1)*dataShifted[ii_nth[cc]].py + u(2)*dataShifted[ii_nth[cc]].pz
      t(cc)=acos(dp/sqrt( (dataShifted[ii_nth[cc]].px)^2 + (dataShifted[ii_nth[cc]].py)^2 + (dataShifted[ii_nth[cc]].pz)^2 ))
    end
    xa=r*cos(t)*u(0) + r*sin(t)*res(0)
    ya=r*cos(t)*u(1) + r*sin(t)*res(1)
    za=r*cos(t)*u(2) + r*sin(t)*res(2)
    geopack_sphcar, xa, ya, za, Rv, Latv, lonv, /to_sphere, /degree
    bg_i=where( abs(lonv-dataShifted[ii_nth].GEI_Lon_deg) gt 180.0, complement=sm_i)        ; across 360 deg
    if bg_i(0) gt -1 then $
    i_bd=where( (360.-abs(lonv[bg_i]-dataShifted[ii_nth[bg_i]].GEI_Lon_deg) gt 15.0) or $
                (abs(lonv[sm_i]-dataShifted[ii_nth[sm_i]].GEI_Lon_deg) gt 15.) ,complement=i_gd ) else $
    i_bd=where( (abs(lonv-dataShifted[ii_nth].GEI_Lon_deg) gt 15.),complement=i_gd )
    if i_bd(0) gt -1 then dataShifted[ii_nth[0]+i_bd].typ=1      ; mark strays as type 1

; code to check sequential locations
    nclat_a=dataShifted[ii_nth[0]+i_gd].GEI_coLat_deg
;    nlon_a =dataShifted[ii_nth[0]+i_gd].GEI_Lon_deg
;    gcDis=gc_dist( 90.0d0-nclat_a(0), 90.0d0-nclat_a, nlon_a(0), nlon_a)
;    i_s=sort(gcDis)
;    print,i_s-shift(i_s,1)

    gap=where(abs(nclat_a-shift(nclat_a,1)) ge SpacRes)
    if gap(0) gt -1 then begin
      print,'WARNING : Data gaps in Nth on Track ',Tr_ii,' at coLats ',dataShifted[ii_nth[0]+gap].GEI_coLat_deg
;      nlon_a =dataShifted[ii_nth[0]+i_gd].GEI_Lon_deg
;      gcDis=gc_dist( 90.0d0-nclat_a(0), 90.0d0-nclat_a, nlon_a(0), nlon_a)
;      i_s=sort(gcDis)
;      lat_g=nclat_a[i_s]-shift(nclat_a[i_s],1)
;    print,i_s-shift(i_s,1)
;      stop
    end

; Check Sth hemis
    ii_sth=where( (dataShifted.iPln eq Tr_ii) and (dataShifted.GEI_colat_deg ge 90.), ii_sth_cnt)
    dataShifted[ii_sth].typ = 0
    u = [dataShifted[ii_sth[0]].px, dataShifted[ii_sth[0]].py, dataShifted[ii_sth[0]].pz]
    u_mag = sqrt(u(0)^2 + u(1)^2 + u(2)^2)
    u=u/u_mag

    sel_lat=where( (dataShifted[ii_sth].GEI_coLat_deg gt 130.) and $
                   (dataShifted[ii_sth].GEI_coLat_deg lt 140.) )  ; assumes strays are interleaved or poleward of this
;                                                                 ; selects both longitude 1/2 spaces
    mn=min(abs(dataShifted[ii_sth[sel_lat]].GEI_Lon_deg-dataShifted[ii_sth[0]].GEI_Lon_deg),loc)   ; avoid strays and get correct longitude
    v = [dataShifted[ii_sth[sel_lat[loc]]].px, dataShifted[ii_sth[sel_lat[loc]]].py, dataShifted[ii_sth[sel_lat[loc]]].pz]
    v_mag = sqrt(v(0)^2 + v(1)^2 + v(2)^2)
    v=v/v_mag

    cross_p,u,v,n
    n_mag = sqrt(n(0)^2 + n(1)^2 + n(2)^2)
    n=n/n_mag

    cross_p,n,u,res

    r=sqrt(dataShifted[ii_sth].px^2 + dataShifted[ii_sth].py^2 + dataShifted[ii_sth].pz^2)
    t=dblarr(ii_sth_cnt)
    t(0)=0.0d0
    for cc=1,ii_sth_cnt-1 do begin
      dp = u(0)*dataShifted[ii_sth[cc]].px + u(1)*dataShifted[ii_sth[cc]].py + u(2)*dataShifted[ii_sth[cc]].pz
      t(cc)=acos(dp/sqrt(dataShifted[ii_sth[cc]].px^2 + dataShifted[ii_sth[cc]].py^2 + dataShifted[ii_sth[cc]].pz^2))
    end
    xa=r*cos(t)*u(0) + r*sin(t)*res(0)
    ya=r*cos(t)*u(1) + r*sin(t)*res(1)
    za=r*cos(t)*u(2) + r*sin(t)*res(2)
    geopack_sphcar, xa, ya, za, Rv, Latv, Lonv, /to_sphere, /degree

    bg_i=where( abs(lonv-dataShifted[ii_sth].GEI_Lon_deg) gt 180.0, complement=sm_i)        ; across 360 deg
    if bg_i(0) gt -1 then $
    i_bd=where( (360.-abs(lonv[bg_i]-dataShifted[ii_sth[bg_i]].GEI_Lon_deg) gt 15.0) or $
                (abs(lonv[sm_i]-dataShifted[ii_sth[sm_i]].GEI_Lon_deg) gt 15.) ,complement=i_gd ) else $
    i_bd=where( (abs(lonv-dataShifted[ii_sth].GEI_Lon_deg) gt 15.),complement=i_gd )

    if i_bd(0) gt -1 then dataShifted[ii_sth[0]+i_bd].typ=1       ; mark as strays

; code to check sequential locations
    sclat_a=dataShifted[ii_sth[0]+i_gd].GEI_coLat_deg
;    slon_a =dataShifted[ii_sth[0]+i_gd].GEI_Lon_deg
;    gcDis=gc_dist( 90.0d0-sclat_a(0), 90.0d0-sclat_a, slon_a(0), slon_a)
;    i_s=sort(gcDis)
;    print,i_s-shift(i_s,1)

    gap=where(abs(sclat_a-shift(sclat_a,1)) ge SpacRes)
    if gap(0) gt -1 then print,print,'WARNING : Data gaps in Sth on Track ',Tr_ii,' at coLats ',dataShifted[ii_sth[0]+gap].GEI_coLat_deg
  end
  r = !null & t = !null
  xa = !null & ya = !null & za = !null
  Latv = !null & Lonv = !null
  ln_a = !null & lt_a = !null
  nclat_a = !null & sclat_a = !null
  i_bd=!null

;  Track order is now sequential and dataShifted structure is now sorted by:
;   (i) orbit plane
;  (ii) distance from equator point where
; (iii) Lon start point is in same semi-circle in longitude (all < 180 deg)

  @amp_fit_paths_full
  gh_clat_lim=25.               ; max coLat to put in ghosts
  ampere_ghost_calc, dataShifted, data, dataGhost_nth, gh_clat_lim, trk_order,path
  gh_clat_lim=180.-gh_clat_lim
  ampere_ghost_calc, dataShifted, data, dataGhost_sth, gh_clat_lim, trk_order,path

  dataShifted=[dataShifted,dataGhost_nth,dataGhost_sth]   ; extend data structure

	; Plot input data
	; --------------------------

	if (debug eq 1) then begin

		device, decomposed = 0
    !p.multi = 0

		window, 0, xSize = 680, ySize = 720, title='Original Data'
;    !p.multi = [0,2,2]
    !p.multi = [0,2,2]
		!p.charsize=1.0

    south=0
 		plt_dat,dataOriginal.GEI_coLat_deg, dataOriginal.GEI_lon_deg, $
			      n_elements(dataOriginal.GEI_coLat_deg), $
			      -dataOriginal.bTheta_GEI, dataOriginal.bPhi_GEI, $
      			south,title='North: Input DataOriginal', capSize=90

		plt_dat,dataShifted.GEI_coLat_deg, dataShifted.GEI_lon_deg, $
			      n_elements(dataShifted.GEI_coLat_deg), $
			      -dataShifted.bTheta_GEI, dataShifted.bPhi_GEI, $
			      south,title='North: Input DataShifted+Ghosts',capSize=90

    south=1
    plt_dat,dataOriginal.GEI_coLat_deg, dataOriginal.GEI_lon_deg, $
            n_elements(dataOriginal.GEI_coLat_deg), $
            -dataOriginal.bTheta_GEI, dataOriginal.bPhi_GEI, $
            south,title='South: Input DataOriginal', capSize=90

    plt_dat,dataShifted.GEI_coLat_deg, dataShifted.GEI_lon_deg, $
            n_elements(dataShifted.GEI_coLat_deg), $
            -dataShifted.bTheta_GEI, dataShifted.bPhi_GEI, $
            south,title='South: Input DataShifted+Ghosts',capSize=90
		!p.multi = 0

	endif

	; Get the average epoch and times for
	; GEI to GEOG and AACGM conversion
	; -----------------------------------

	print,'Year, Month, Day, Hour, Minute : ',year,' ',month,' ',day,' ',fix(shr),' ',round((shr mod 1)*60)
  avgHour = fix ( (sHr + eHr) / 2.0 )
  avgMin  = ( (sHr + eHr) / 2.0 mod 1 ) * 60
	cdf_epoch, epoch0, year, month, day, ((avgHour lt 24)?avgHour:avgHour-24), avgMin,/compute_epoch
	if (avgHour ge 24) then epoch0=epoch0+86400.0d3
	epoch = data.utc*3600d0*1d3 + epoch0
	avgEpoch = mean(epoch)
	yrSecAvg = cnvTime ( year, month, day, avgHour, avgMin, 0.0 )
	geoPack_reCalc, year, month, day, avgHour, avgMin,/date

end
