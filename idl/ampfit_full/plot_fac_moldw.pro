Pro plot_fac_moldw,jPar, Lat_deg, lon_deg, mn_fac, mx_fac, title = title
; Plot FAC on Mollweide projection - both Nth/Sth to show equatorial FAC
; CLW - Sept 2010
;
  pole=0.
  limit=[-90,-180,90,180]
  map_set, pole, 0, 0, /mollweide, limit = limit, /continents, /grid, color = 255, /iso

; order by small to large longitudes
  srt=sort(lon_deg)
  lonTmp=lon_deg[srt]
  LatTmp=Lat_deg[srt]
  jParTmp=jPar[srt]

  randFac = 1e-2
  lonTmp = ((lonTmp[*]/15.0) mod 24 )*15.0  $
      + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*randFac-randFac/2

  latTmp = LatTmp[*] $
      + randomu(sysTime(/sec ), n_elements(jParTmp), /uni)*randFac-randFac/2

  idx=where( abs(latTmp[*]) lt 90.0)   ; Check for correct Lat range
  if (idx(0) ne -1) then begin
    latTmp=latTmp[idx]
    lonTmp=lonTmp[idx]
    jParTmp=jParTmp[idx]
  end

  nLevs = 10
  fac_div=mx_fac/10.0
  jLevels = (fIndGen(2*nLevs+1)-nLevs)*fac_div
  colors  = bytScl ( jLevels, top = 253 ) + 1

;  transparency=50
;c= contour(final, lons*!radeg, 90.-lats*!radeg, $
;       /over, levels = jLevels, transparency=transparency, $
;            c_colors = colors, /fill, /irreg)

  contour, jParTmp, $        ; assumes data starts at long=-180, lat=-90
      lonTmp, $
      latTmp, $
       /over, levels = jLevels, $
            c_colors = colors, /irreg, /cell_fill
  map_continents

end