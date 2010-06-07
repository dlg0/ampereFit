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


