pro plot_fac_cw, jPar, pole, cap_coLat_deg, coLat_deg, lon_deg, mn_fac, mx_fac, path, $
		title = title, south = south

	if not keyword_set(noErase) then noErase = 0

	loadct, 13, file = path + 'davect.tbl',/silent

  pole=abs(pole)
  cap_sz_deg=cap_coLat_deg
  if south then begin
    cap_sz_deg=180.0-cap_sz_deg
	idx=where( (coLat_deg gt 92.) and (coLat_deg lt 178.) )
    if idx(0) gt -1 then begin
      clat_d=180.-coLat_deg[idx]         ; mirror to look thru Earth
      lon_d=lon_deg[idx]
      jParTmp=jPar[idx]
    end else begin
      clat_d=180.0-coLat_deg
      lon_d=Lon_Deg
      jParTmp=jPar
    endelse
  endif else begin
    idx=where( (coLat_deg lt 88.) and (coLat_deg gt 2.) )
    if idx(0) gt -1 then begin
      clat_d=coLat_deg[idx]
      lon_d=lon_deg[idx]
      jParTmp=jPar[idx]
    end else begin
      clat_d=coLat_deg
      lon_d=Lon_Deg
      jParTmp=jPar
    endelse
  endelse
  limit = [ pole - cap_sz_deg, 0, pole, 360 ] ; use Nth hemis for both Nth and Sth views

  map_set, pole, 0, 0, /ortho, /iso, $
     limit = limit, $
	 /noborder, title = title, /advance, $
	 yMargin = [0,3];, color = 0
  fc_idx=where(abs(jParTmp) lt mn_fac)
  if fc_idx(0) gt -1 then jParTmp(fc_idx)=0.0      ; Check for min val, set to white
  fc_idx=where(jParTmp gt mx_fac)
  if fc_idx(0) gt -1 then jParTmp(fc_idx)=mx_fac
  fc_idx=where(jParTmp lt -mx_fac)
  if fc_idx(0) gt -1 then jParTmp(fc_idx)=-mx_fac

  randFac = 1e-2
  lonTmp = ((lon_d[*]/15.0) mod 24 )*15.0 $
	  + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*randFac-randFac/2.0

  latTmp = 90.0-clat_d[*] $
	  + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*randFac-randFac/2.0

  nLevs = 10
  fac_div=mx_fac/10.0
  jLevels = (fIndGen(2*nLevs+1)-nLevs)*fac_div
  colors  = bytScl ( jLevels, top = 253 ) + 1

  contour, jParTmp, lonTmp, latTmp, /over,$
           levels = jLevels, $
;            c_colors = colors, /cell_fill, /irreg
            c_colors = colors, /fill, /irreg

	; Just grid both north and south hemispheres
	; ------------------------------------------

  nLats = fix((cap_sz_deg)/10)
  lats = 90-(fIndGen(nLats)+1)*10

  lnlbp=min(abs(lats),idx)
  lnpos=lats(idx)+3
  lt_lab=strtrim(string(lats,format='(i2)'),2)
  if south then lt_lab=strtrim(string(-lats,format='(i3)'),2)

	map_grid, label = 1, $
			lonNames	= ['0','6','12','18',''], $
			lons = [0,90,180,270,360], $
			lonlab=lnpos, $
			latLab = 45, $
			lats = lats, $
            latNames=lt_lab, $
			color = 255

end
