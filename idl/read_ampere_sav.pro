; Reads AMPERE data SAV file (input data)
;
; David L Green & Colin L Waters
; Centre for Space Physics
; University of NEwcastle
; Dec 2009
;
; Comments:
; data selection -> latitude range set to max_coLat
; modified may 2010 - sort data by orbit plane and sequentially alomg track (for diagnostics)
;                   - tweak ampdata shift routine
;
pro read_ampere_sav, sHr, eHr, south, rd_coLat, $
    savFileName = savFileName, $
    dataOut = dataOut, $
    year = year, month = month, day = day, $
    avgYrSec = yrSecAvg, $
    avgEpoch = avgEpoch, $
    rot_mat=rot_mat

 Forward_function cnvTime                ; IDL sometimes treats this as a variable. This stops it doing that

 if not keyword_set(max_coLat) then max_coLat=70.0
 dateStr=''
 reads, strmid(file_baseName ( savFileName ),0,8), year, month, day, format='(i4,i2,i2)'
 restore, savFileName
;   variables are
;
;   x_axis_frac_hour   : UTC in fractional hours
;   plane_number_total : integer number 0-5 for the 6 orbit planes.
;   pos_ECI_total      : ECI position in meters, dimensions [3,x] where x is number of points, 0-2 in the first dimension is X,Y,Z
;   B_ECI              : the same dimension as pos_eci but has delta B in nT

 data_struct={ipln: 0, $
              isat: 0, $
              utc: 0d0, $
              px: 0d0, py: 0d0, pz: 0d0, $
              dbx: 0d0, dby: 0d0, dbz: 0d0, $
              pr : 0d0, pth: 0d0, pph: 0d0, $
              dbr: 0d0, dbth: 0d0, dbph: 0d0 }
 ntot=n_elements(x_axis_frac_hour)
 rwpz_a=dblarr(ntot)
 rwpz_a=pos_eci_total(2,*)/1000.0   ; conv Z coord to km
 earth_rad=6371.0                   ; Earth radius in km
 sat_rad=780.0                      ; Iridium satellite orbit height in km

; Select data subset based on time interval and hemisphere : above 90-max_coLat deg lat.
 z0_lim= (earth_rad + sat_rad)*cos(rd_coLat*!pi/180.0)         ; z-comp of max_coLat

 If south then begin
  z0_lim=-z0_lim
  iiTime=where(x_axis_frac_hour ge sHr and x_axis_frac_hour le eHr $
                                      and rwpz_a lt z0_lim, iiTimeCnt)
 end else begin
  iiTime=where(x_axis_frac_hour ge sHr and x_axis_frac_hour le eHr $
                                      and rwpz_a gt z0_lim, iiTimeCnt)  ; Nth Hemisphere
 end
;stop
; Sort data along each orbit track - unnecessary for the fit, but very useful for debugging
 data=replicate(data_struct, iiTimeCnt)
 ipln_tmp= plane_number_total[iiTime]               ; orbit plane number
 isat_tmp= pseudoSVNum_total[iiTime]                ; satellite number
 utc_tmp = x_axis_frac_hour[iiTime]                 ; UT time

 px_tmp = (pos_ECI_total[0,*])[iiTime]*1d-3         ; Cartesian (geog) x posn, conv to km
 py_tmp = (pos_ECI_total[1,*])[iiTime]*1d-3
 pz_tmp = (pos_ECI_total[2,*])[iiTime]*1d-3

 dbx_tmp = (B_ECI[0,*])[iiTime]                     ; Cartesion (geog) bx vector
 dby_tmp = (B_ECI[1,*])[iiTime]
 dbz_tmp = (B_ECI[2,*])[iiTime]

 data.ipln= ipln_tmp
 data.isat= isat_tmp
 data.utc = utc_tmp

 data.px = px_tmp
 data.py = py_tmp
 data.pz = pz_tmp

 data.dbx = dbx_tmp
 data.dby = dby_tmp
 data.dbz = dbz_tmp

 st=0
 For tt=0,5 do begin
  tt_idx=where(ipln_tmp eq tt)              ; search selected data values (time, lat)

  ipln_tmp_pl= ipln_tmp[tt_idx]             ; orbit plane number, one orbit plane at a time
  isat_tmp_pl= isat_tmp[tt_idx]             ; satellite number
  utc_tmp_pl = utc_tmp[tt_idx]              ; UT time

  px_tmp_pl = px_tmp[tt_idx]                ; ECI location, one orbit plane
  py_tmp_pl = py_tmp[tt_idx]
  pz_tmp_pl = pz_tmp[tt_idx]

  dbx_tmp_pl = dbx_tmp[tt_idx]              ; dB, one orbit plane
  dby_tmp_pl = dby_tmp[tt_idx]
  dbz_tmp_pl = dbz_tmp[tt_idx]

  dist=sqrt(px_tmp_pl^2 + py_tmp_pl^2)      ; dist calc on x,y only
  mx=max(dist,mx_idx)                       ; find location furthest away
  rel_dist=sqrt( (px_tmp_pl-px_tmp_pl[mx_idx])^2 + (py_tmp_pl-py_tmp_pl[mx_idx])^2)   ; cal dist from this furthest point
  srt_idx=sort(rel_dist)                    ; array index sort pattern
  np_pln=n_elements(srt_idx)

  data[st:st+np_pln-1].ipln = ipln_tmp_pl[srt_idx]    ; put sorted data back into data structure
  data[st:st+np_pln-1].isat = isat_tmp_pl[srt_idx]
  data[st:st+np_pln-1].utc = utc_tmp_pl[srt_idx]

  data[st:st+np_pln-1].px = px_tmp_pl[srt_idx]
  data[st:st+np_pln-1].py = py_tmp_pl[srt_idx]
  data[st:st+np_pln-1].pz = pz_tmp_pl[srt_idx]

  data[st:st+np_pln-1].dbx = dbx_tmp_pl[srt_idx]
  data[st:st+np_pln-1].dby = dby_tmp_pl[srt_idx]
  data[st:st+np_pln-1].dbz = dbz_tmp_pl[srt_idx]
  st=st+np_pln
 end
; end of track sort section
;stop
;stop
; Get spherical coords of the GEI XYZ locations
; sphcar,rg_th,cglat,glon,data.px,data.py,data.pz,/to_sphere,/degree   ; Get GEI(r,thet,phi)
; np=n_elements(data.px)
; rg_th=dblarr(np) & cglat=dblarr(np) & glon=dblarr(np)
 geopack_sphcar,data.px,data.py,data.pz,rg_th,cglat,glon,/to_sphere,/degree   ; Get GEI(r,thet,phi)
 data.pr=rg_th
 data.pth=cglat                                   ; coLat: 0->90 for Nth and 90->180 for Sth (deg)
 data.pph=glon

; Put the XYZ GEI input db data into spherical components
  geopack_bcarsp, data.px, data.py, data.pz, data.dbx, data.dby, data.dbz, $
          dbr, dbth, dbph
 data.dbr = dbr
 data.dbth= dbth                                   ; coLat: 0->90 for Nth and 90->180 for Sth (deg)
 data.dbph= dbph

; For ii=Long(0),iiTimeCnt-Long(1) do begin               ; Nth has pz +ve, Sth has pz -ve
;  bcarsp, data[ii].px, data[ii].py, data[ii].pz, $
;          data[ii].dbx, data[ii].dby, data[ii].dbz, $
;          vbr,vbth,vbph
;  data[ii].dbr=vbr
;  data[ii].dbth=vbth
;  data[ii].dbph=vbph
; end

; plot the db vectors : for testing
; window, 1, xSize = 500, ySize = 500,title='Input Data'
; If south Then plt_dat,data.pth,data.pph,n_elements(data.pth),data.dbth,data.dbph,[1,1],[1,2],south,title='Input Data',capSize=max_coLat else $
;  plt_dat,data.pth,data.pph,n_elements(data.pth),-data.dbth, data.dbph,[1,1],[1,2],south,title='Input Data'
;
; ***   ***   ***   ***   ***   ***
; Shift the data
 amp_shiftdata,data,data_sh,rot_mat,south;,/diag
; ************** Bypass the data shift routine
; data_sh=data
; rot_mat=[[1.0,0.0,0.0],$
;          [0.0,1.0,0.0],$
;          [0.0,0.0,1.0]]
; *********************************************

; plot the db vectors : for testing
; window, 2, xSize = 500, ySize = 500,title='Input Data Shifted'
; If south Then plt_dat,data_sh.pth,data_sh.pph,n_elements(data_sh.pth),data_sh.dbth,data_sh.dbph,[1,1],[1,2],south,title='Input Data',capSize=max_coLat else $
;  plt_dat,data_sh.pth,data_sh.pph,n_elements(data_sh.pth),-data_sh.dbth, data_sh.dbph,[1,1],[1,2],south,title='Input Data: Shifted'
;
; input data
 gei_R_km=data.pr
 gei_coLat_rad=data.pth*!dpi/180.0
 gei_lon_rad=data.pph*!dpi/180.0

 br_GEI=data.dbr
 bTheta_GEI=data.dbth
 bPhi_GEI=data.dbph

; Shifted GEI data
 gei_R_km_sh=data_sh.pr
 gei_coLat_rad_sh=data_sh.pth*!dpi/180.0
 gei_lon_rad_sh=data_sh.pph*!dpi/180.0

 br_GEI_sh=data_sh.dbr
 bTheta_GEI_sh=data_sh.dbth
 bPhi_GEI_sh=data_sh.dbph

; ********* Check this ************
; If south then bTheta_GEI_sh = -bTheta_GEI_sh
;stop

;   get the epoch and times for GEI to GEOG and AACGM conversion
 Print,'Year, Month, Day : ',year,' ',month,' ',day
 cdf_epoch, epoch0, year, month, day, /compute_epoch
 epoch=data.utc*3600d0*1d3 + epoch0
 avgEpoch = mean(epoch)
 avgHour	= fix ( (sHr + eHr) / 2.0 )
 avgMin	= ( (sHr + eHr) / 2.0 mod 1 ) * 60
 yrSecAvg = cnvTime ( year, month, day, avgHour, avgMin, 0.0 )
 geoPack_reCalc, year, month, day, /date              ; initialise geoPack using data/time

; Reminder: Data in the output is
; Shifted GEI data
; gei_R_km_sh=data_sh.pr
; gei_coLat_rad_sh=data_sh.pth*!dpi/180.0
; gei_lon_rad_sh=data_sh.pph*!dpi/180.0
;
; Shifted dB data
; br_GEI_sh=data_sh.dbr
; bTheta_GEI_sh=data_sh.dbth
; bPhi_GEI_sh=data_sh.dbph
;
 dataOut = { gei_R_km : gei_R_km, $              ; Input data - unmodified
		     gei_coLat_rad : gei_coLat_rad, $
			 gei_lon_rad : gei_lon_rad, $
			 dbR : bR_GEI, $
			 dbTheta : bTheta_GEI, $
			 dbPhi : bPhi_GEI, $
			 gei_R_km_sh : gei_R_km_sh, $           ; Shifted input data
		     gei_coLat_rad_sh : gei_coLat_rad_sh, $
			 gei_lon_rad_sh : gei_lon_rad_sh, $
			 dbR_sh : bR_GEI_sh, $
			 dbTheta_sh : bTheta_GEI_sh, $
			 dbPhi_sh : bPhi_GEI_sh, $
             iPln : data.iPln, $
             iSat : data.iSat, $
             utc : data.utc }
end
