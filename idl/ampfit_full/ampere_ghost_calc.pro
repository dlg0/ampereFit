pro ampere_ghost_calc, data_in, data, ghdata_out, gh_clat_lim, trk_order, path
; Polar cap ghost calculation routine
;  Given that the Iridium orbit tracks do not necessarily cross at the poles at the same location
;  the data fit routine for high m-order does not have sufficient coverage to avoid a Nyquist prob
;  This routine 'ghosts' data points around the polar regions in longitude
;
; C. L. Waters
;  Space Physics Group
;  University of Newcastle
;  New South Wales, Australia
;  Sep 2010
; 
;  Inputs:
;        data_in         : Input Iridium data structure - shifted
;        data            : Defn of data structure
;        gh_clat_lim     : coLat of ghost Lat limit (< 90 for Nth, > 90 for Sth)
;        trk_order       ; orbit number sequence (array[6])
;        
;  Output:
;        ghdata_out      : ghost data with same structure as input data
;
; **********************************************************************************

dbg=0

  ex_gh_cnt=0
  
; cycle thru each track pair to add polar ghosts [CLW]
  if gh_clat_lim lt 90 then g_idx=where( (data_in.GEI_coLat_deg le gh_clat_lim) and $   ; longitude ghost for < gh_clat_lim deg coLat
                                         (data_in.typ eq 0) ) $
    else g_idx=where( (data_in.GEI_coLat_deg ge gh_clat_lim) and (data_in.typ eq 0) )
  cnt=0

  ghdata_out = replicate( data, n_elements(g_idx) )

  For TrN=0,5 do begin
    if (gh_clat_lim lt 90.) then begin                                    ; Nth hemis
      Track_1 = where( (data_in.iPln eq Trk_order[TrN]) and $             ; select track number
                       (data_in.GEI_colat_deg le gh_clat_lim) and $       ; select by coLat range
                       (data_in.typ eq 0), Tr_1_ii )                      ; exclude strays
    end else $                                                            ; Sth hemis
      Track_1 = where ( (data_in.iPln eq Trk_order[TrN]) and $
                        (data_in.GEI_colat_deg ge gh_clat_lim) and $
                        (data_in.typ eq 0), Tr_1_ii )

    if TrN lt 5 then NTrack=TrN+1 else NTrack=0                          ; cycle thru each Track - they are ordered sequentially

    if (gh_clat_lim lt 90.) then begin
      Track_2 = where ( (data_in.iPln eq Trk_order[NTrack]) and $
                        (data_in.GEI_colat_deg le gh_clat_lim) and $
                        (data_in.typ eq 0), Tr_2_ii ) 
    end else $
      Track_2 = where ( (data_in.iPln eq Trk_order[NTrack]) and $
                        (data_in.GEI_colat_deg ge gh_clat_lim) and $
                        (data_in.typ eq 0), Tr_2_ii )
;print,'Track number pairs:',Tr_srt[TrN],Tr_srt[NTrack]
;print,'start lat, lon : ',90.-data_in[Track_1[0]].GEI_coLat_deg, data_in[Track_1[0]].GEI_Lon_deg

if dbg eq 1 then begin
  device, decomposed = 0
  !p.multi = 0
  window, 0, xSize = 680, ySize = 720, title='Test Data'
  ncol    = 256
  col     = 0
  IrCol   = col
  cLatMin = 30.0
  nlatlab = cLatMin / 10
  mmg     = 7.0
  DMag    = 200.
  XSh     = 0.0
  YSh     = 0.0
  x       = dblArr(nlatlab,180)
  y       = dblArr(nlatlab,180)
  for i=0,nlatlab-1 do begin
    for j=0,179 do begin
      x(i,j)=float(i+1)*10.*cos(float(j)/90.*!pi)
      y(i,j)=float(i+1)*10.*sin(float(j)/90.*!pi)
      endfor
  endfor
  plot,x(nlatlab-1,*),y(nlatlab-1,*),$
    /isotropic, xmargin=[1,1], ymargin=[1,1], xstyle=5,ystyle=5, $
  background=255, color=0, linestyle=1, title = title
  for i=0,nlatlab-2 do oPlot,x(i,*),y(i,*), linestyle=1, color=0
  loadct, 1,/silent
end

    For jj=0,n_elements(Track_1)-1 do begin
      strt_lat=data_in[Track_1[jj]].GEI_coLat_deg                        ; cycle thru each point on 1st track
      strt_lon=data_in[Track_1[jj]].GEI_Lon_deg
      gc_d=gc_dist( 90.-strt_lat, 90.-data_in[Track_2].GEI_coLat_deg, $  ; search all of track 2 using great circle distance formula
                    strt_lon, data_in[Track_2].GEI_Lon_deg)
      gc_dmn=min(gc_d,mn_idx)                                            ; this is the closest point on the 2nd track

      ang=Strt_lon*!pi/180.                                              ; do dist calc in X,Y to get midpoint over weird Lons
      if gh_clat_lim gt 90. then strt_lat=180.-strt_lat
      Xs=Strt_lat*cos(ang)
      Ys=Strt_lat*sin(ang)
      end_lat=data_in[Track_2[mn_idx]].GEI_coLat_deg                     ; this is the coLat, Lon coord for the 2nd track point
      if gh_clat_lim gt 90. then end_lat=180.-end_lat
      end_lon=data_in[Track_2[mn_idx]].GEI_Lon_deg
      Nyq_Ok=1                                                           ; set Flag to true
      if ( abs(end_lon-strt_lon) ge 120. and $                           ; 2 times 6 orbit plane spacing (ghosts 1/2 the spacing)
           abs(end_lon-strt_lon) lt 180.) then begin
            print,'Nyquist violation in Longitude at CoLat = ',strt_lat, '  Lon gap = ',abs(end_lon-strt_lon),' fixing...'
            Nyq_Ok=0          ; toggle flag for Nyquist
            ex_gh_cnt++
            if ex_gh_cnt eq 1 then begin
              ex_gh_file=path+'more_ghosts.dat'
              openw,ug,ex_gh_file,/get_lun      ; open file in case its needed
              exghdata = replicate( data, 1 )         ; structure to push to file if extra ghosts are needed
            end
      end
; now calc the midpoint - in 'XY'
      ang=end_lon*!pi/180.
      Xe=end_lat*cos(ang)
      Ye=end_lat*sin(ang)
      If Nyq_Ok then begin
        Xg=(Xs+Xe)/2.
        Yg=(Ys+Ye)/2.
      end else begin
        Xg1=2.0d0/3.0d0*Xs + Xe/3.0d0
        Yg1=2.0d0/3.0d0*Ys + Ye/3.0d0
        Xg2=Xs/3.0d0 + 2.0d0/3.0d0*Xe
        Yg2=Ys/3.0d0 + 2.0d0/3.0d0*Ye
      end

if dbg eq 1 then begin
  xyouts,Xs,Ys,strtrim(string(data_in[Track_1[jj]].iPln),2),color=0
  If Nyq_Ok then xyouts,Xg,Yg,'G',color=0 else begin
    xyouts,Xg1,Yg1,'a',color=0
    xyouts,Xg2,Yg2,'b',color=0
  end
  xyouts,Xe,Ye,strtrim(string(data_in[Track_2[mn_idx]].iPln),2),color=0
end

; put data in Ghost data structure
      ghdata_out[cnt].utc=(data_in[Track_1[jj]].utc + data_in[Track_2[mn_idx]].utc)/2.0d0
      ghdata_out[cnt].isat=-1
      ghdata_out[cnt].iPln=data_in[Track_1[jj]].iPln
      ghdata_out[cnt].splice=data_in[Track_1[jj]].splice

      if Nyq_Ok eq 0 then begin         ; need extra ghosts
        exghdata.utc=ghdata_out[cnt].utc
        exghdata.isat=-1
        exghdata.iPln=data_in[Track_1[jj]].iPln

        ghdata_out[cnt].qual=2.0d0/3.0d0*data_in[Track_1[jj]].qual + 1.0d0/3.0d0*data_in[Track_2[mn_idx]].qual
        exghdata.qual=1.0d0/3.0d0*data_in[Track_1[jj]].qual + 2.0d0/3.0d0*data_in[Track_2[mn_idx]].qual

        exghdata.splice=data_in[Track_1[jj]].splice
        ghdata_out[cnt].typ=3
        exghdata.typ=3

        ghdata_out[cnt].GEI_R_km=2.d0/3.d0*data_in[Track_1[jj]].GEI_R_km + 1.d0/3.d0*data_in[Track_2[mn_idx]].GEI_R_km
        exghdata.GEI_R_km=1.d0/3.d0*data_in[Track_1[jj]].GEI_R_km + 2.d0/3.d0*data_in[Track_2[mn_idx]].GEI_R_km

        if gh_clat_lim gt 90. then begin      ; Sth hemis
          ghdata_out[cnt].GEI_coLat_deg=180.d0-sqrt(Xg1^2+Yg1^2)
          exghdata.GEI_coLat_deg=180.d0-sqrt(Xg2^2+Yg2^2)
        end else begin                        ; Nth hemis
          ghdata_out[cnt].GEI_coLat_deg=sqrt(Xg1^2+Yg1^2)
          exghdata.GEI_coLat_deg=sqrt(Xg2^2+Yg2^2)
        end
        ghdata_out[cnt].GEI_coLat_rad = ghdata_out[cnt].GEI_coLat_deg*!dtor
        exghdata.GEI_coLat_rad = exghdata.GEI_coLat_deg*!dtor

        ghdata_out[cnt].GEI_Lon_deg=atan(Yg1,Xg1)*180.d0/!dpi
        exghdata.GEI_Lon_deg=atan(Yg2,Xg2)*180.d0/!dpi
        ghdata_out[cnt].GEI_Lon_rad = ghdata_out[cnt].GEI_Lon_deg*!dtor
        exghdata.GEI_Lon_rad = exghdata.GEI_Lon_deg*!dtor

        geopack_sphcar, ghdata_out[cnt].GEI_R_km, ghdata_out[cnt].GEI_coLat_deg, ghdata_out[cnt].GEI_Lon_deg, $
              xv, yv, zv, /to_rect, /degree
        ghdata_out[cnt].px=xv & ghdata_out[cnt].py=yv & ghdata_out[cnt].pz=zv
        geopack_sphcar, exghdata.GEI_R_km, exghdata.GEI_coLat_deg, exghdata.GEI_Lon_deg, $
              xv, yv, zv, /to_rect, /degree
        exghdata.px=xv & exghdata.py=yv & exghdata.pz=zv

        ghdata_out[cnt].br_GEI=2.d0/3.d0*data_in[Track_1[jj]].br_GEI + 1.d0/3.d0*data_in[Track_2[mn_idx]].br_GEI
        exghdata.br_GEI=1.d0/3.d0*data_in[Track_1[jj]].br_GEI + 2.d0/3.d0*data_in[Track_2[mn_idx]].br_GEI

        ghdata_out[cnt].bTheta_GEI=2.d0/3.d0*data_in[Track_1[jj]].bTheta_GEI + 1.d0/3.d0*data_in[Track_2[mn_idx]].bTheta_GEI
        exghdata.bTheta_GEI=1.d0/3.d0*data_in[Track_1[jj]].bTheta_GEI + 2.d0/3.d0*data_in[Track_2[mn_idx]].bTheta_GEI

        ghdata_out[cnt].bPhi_GEI=2.d0/3.d0*data_in[Track_1[jj]].bPhi_GEI + 1.d0/3.d0*data_in[Track_2[mn_idx]].bPhi_GEI
        exghdata.bPhi_GEI=1.d0/3.d0*data_in[Track_1[jj]].bPhi_GEI + 2.d0/3.d0*data_in[Track_2[mn_idx]].bPhi_GEI

        geopack_bspcar, ghdata_out[cnt].GEI_coLat_deg, ghdata_out[cnt].GEI_Lon_deg, $
                        ghdata_out[cnt].br_GEI, ghdata_out[cnt].bTheta_GEI, ghdata_out[cnt].bPhi_GEI, $ 
                        bxv, byv, bzv, /degree
        ghdata_out[cnt].dbx=bxv & ghdata_out[cnt].dby=byv & ghdata_out[cnt].dbz=bzv
        geopack_bspcar, exghdata.GEI_coLat_deg, exghdata.GEI_Lon_deg, $
                        exghdata.br_GEI, exghdata.bTheta_GEI, exghdata.bPhi_GEI, $ 
                        bxv, byv, bzv, /degree
        exghdata.dbx=bxv & exghdata.dby=byv & exghdata.dbz=bzv

        writeu,ug,exghdata
      end

      if Nyq_Ok then begin   ; only need 1 ghost point
        ghdata_out[cnt].qual=(data_in[Track_1[jj]].qual + data_in[Track_2[mn_idx]].qual)/2.0d0
        ghdata_out[cnt].typ=2

        ghdata_out[cnt].GEI_R_km=(data_in[Track_1[jj]].GEI_R_km + data_in[Track_2[mn_idx]].GEI_R_km)/2.0d0
        if gh_clat_lim gt 90. then ghdata_out[cnt].GEI_coLat_deg=180.0d0-sqrt(Xg^2+Yg^2) else $   ; Sth
                                   ghdata_out[cnt].GEI_coLat_deg=sqrt(Xg^2+Yg^2)               ; Nth
        ghdata_out[cnt].GEI_Lon_deg=atan(Yg,Xg)*180.0d0/!dpi

        ghdata_out[cnt].GEI_coLat_rad = ghdata_out[cnt].GEI_coLat_deg*!dtor
        ghdata_out[cnt].GEI_Lon_rad = ghdata_out[cnt].GEI_Lon_deg*!dtor

        geopack_sphcar, ghdata_out[cnt].GEI_R_km, ghdata_out[cnt].GEI_coLat_deg, ghdata_out[cnt].GEI_Lon_deg, $
                        xv, yv, zv, /to_rect, /degree
        ghdata_out[cnt].px=xv & ghdata_out[cnt].py=yv & ghdata_out[cnt].pz=zv

        ghdata_out[cnt].br_GEI=(data_in[Track_1[jj]].br_GEI + data_in[Track_2[mn_idx]].br_GEI)/2.0d0
        ghdata_out[cnt].bTheta_GEI=(data_in[Track_1[jj]].bTheta_GEI + data_in[Track_2[mn_idx]].bTheta_GEI)/2.0d0
        ghdata_out[cnt].bPhi_GEI=(data_in[Track_1[jj]].bPhi_GEI + data_in[Track_2[mn_idx]].bPhi_GEI)/2.0d0

        geopack_bspcar, ghdata_out[cnt].GEI_coLat_deg, ghdata_out[cnt].GEI_Lon_deg, $
                        ghdata_out[cnt].br_GEI, ghdata_out[cnt].bTheta_GEI, ghdata_out[cnt].bPhi_GEI, $ 
                        bxv, byv, bzv, /degree
        ghdata_out[cnt].dbx=bxv & ghdata_out[cnt].dby=byv & ghdata_out[cnt].dbz=bzv
      end

      cnt++
    end           ; number of ghosts
    if dbg then stop

  endfor          ; number of tracks

  if ex_gh_cnt gt 0 then begin
    close,ug
    free_lun,ug
    mor_gh = replicate( data, ex_gh_cnt )
    exghdata = replicate( data, 1 )
    openr,ug,ex_gh_file,/get_lun      ; open file for read
    for n_ex=0,ex_gh_cnt-1 do begin
      readu,ug,exghdata
      mor_gh[n_ex]=exghdata
    end
    close,ug         ; release extra ghost file
    free_lun,ug
  end

  if ex_gh_cnt gt 0 then ghdata_out = [ ghdata_out[0:cnt-1], mor_gh] else ghdata_out = ghdata_out[0:cnt-1]
  gd_i=where( (ghdata_out.GEI_coLat_deg gt 1.) and (ghdata_out.GEI_coLat_deg lt 179.))
  tmp_dat=ghdata_out   ; make a copy
  ghdata_out = tmp_dat[gd_i]
  tmp_dat = !null

end  