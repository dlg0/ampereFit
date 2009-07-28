; Plot of dB vectors
; Sub_code to plot arrow vectors
;
; Main code is further down
;
Pro plt_dat,IrCLat,IrMLT,nirid,IrNCmp,IrECmp,xm,ym
 ncol=256
 col=0
 IrCol=col
 nlatlab=5
 CLatMin=50.
 mmg=7.0
 DMag=500.
 XSh=0.0
 YSh=0.0
 x=DblArr(nlatlab,180)
 y=DblArr(nlatlab,180)
 For i=0,nlatlab-1 do $
 Begin
  For j=0,179 do $
  Begin
   x(i,j)=float(i+1)*10.*cos(float(j)/90.*!pi)
   y(i,j)=float(i+1)*10.*sin(float(j)/90.*!pi)
  end
 end
 Plot,x(nlatlab-1,*),y(nlatlab-1,*),/isotropic,xmargin=[xm(0),xm(1)],ymargin=[ym(0),ym(1)],xstyle=5,ystyle=5, $
   background=255,color=0,linestyle=1
 For i=0,nlatlab-2 do OPlot,x(i,*),y(i,*),linestyle=1,color=0
 For ii=0,nirid-1 do $
 Begin
  SLon=IrMLT(ii)+180.                ; Shifted coords Lon,Lat
  SLat=IrCLat(ii)
  ang=(SLon-90.)*!pi/180.
  Xs=SLat*cos(ang)              ; Shifted coords, X,Y
  Ys=SLat*sin(ang)
; do N shift (N dir unit vector)
  Mg1=sqrt(Xs*Xs+Ys*Ys)
  xN=Xs/Mg1
  yN=Ys/Mg1
; do E shift (anticlk 90 deg rotation)
  xE=-yN
  yE=xN
  NVal=mmg*IrNCmp(ii)/DMag  ; shifted coords B Val in theta,phi coords
  EVal=mmg*IrECmp(ii)/DMag
  x1=Xs-NVal*xN+EVal*xE
  y1=Ys-NVal*yN+EVal*yE
  x0=Xs-YSh                         ; Convert to geomag but plot is rotated 90 Deg
  x1=x1-YSh
  y0=Ys+XSh
  y1=y1+XSh
  If (sqrt(x0*x0+y0*y0) le CLatMin) Then $
  Begin
   plots,[x0,x1],[y0,y1],color=IrCol,thick=1.5
; Now do arrow heads
   mgn=sqrt((x0-x1)^2+(y0-y1)^2)
   if (mgn gt 0.) then $
   begin
    xu=(x0-x1)/mgn
    yu=(y0-y1)/mgn
    x2=0.4*xu-0.4*yu
    y2=0.4*xu+0.4*yu
    plots,[x1,x1+x2],[y1,y1+y2],color=IrCol,thick=1.5
    x2=0.4*xu+0.4*yu
    y2=-0.4*xu+0.4*yu
    plots,[x1,x1+x2],[y1,y1+y2],color=IrCol,thick=1.5
   end
  end    ; If within Lat Range
 end     ; NIrid Loop
end      ; If Plot Irid
;
; ***********************************************************************************
;
pro clw_amp_v2
; @'d:\cwac\hi_res\davidg\jpar_ver2\schabasisfunctions.pro'
 CapSize=50.
 winnum=1
 fileName = 'd:\cwac\hi_res\20080105_a_RevA.dat'
 sHr = 19
 eHr = 21
 read_ampere_dlg, fileName, sHr, eHr, data, yr, t_arr

 set_plot,'win'
 device,decomposed=0
 window,winnum,xsize=600,ysize=550,title='dgIrid Data'
 winnum++
 LoadCT,0
 device, decomposed = 0
 !p.background=255
 plt_dat, data.theta * !radeg, data.phi * !radeg, $
	n_elements ( data.theta ), -data.dbTheta, data.dbPhi, [1,1], [1,1]
 kMax    = 16
 mMax    = 5
stop
 schaBasisFunctions, kMax, mMax, capSize, data.theta, data.phi, $
         YkmBFns=YkmBFns, dYkmDthBFns=dYkmDthBFns, dYkmDphBFns=dYkmDphBFns,$
         OUTNKVALUES=outNkValues, OUTKVALUES=outKValues, OUTMVALUES=outMValues;, /oddSet
;
;       Fit |dB| to dB.grad Ykm in the Iridium system
 dBMag   = sqrt(data.dBTheta^2+data.dBPhi^2)
 kTh     = transpose(rebin(data.dBTheta/dBMag,$
        n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))
 kPh     = transpose(rebin(data.dBPhi/dBMag,$
         n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))

 bFuncs  = kTh*dYkmDphBFns + kPh*dYkmDthBFns
 alpha   = transpose(bFuncs) ## bFuncs
 beta_   = transpose( transpose(bFuncs) ## dBMag)
 la_svd, alpha, w, u, v, status=svdStatus
 kk      = where(w lt max ( w ) * 1d-5, kkcnt)
 if kkcnt gt 0 then w[kk]=0.0
 coeffs  = svsol(u, w, v, beta_)

 fit_dBTheta     = dYkmDPhBFns ## coeffs
 fit_dBPhi       = dYkmDThBFns ## coeffs
 fit_dBMag       = bFuncs ## coeffs

 Re=6.371d6
 R     = Re + 780.0d3
 u0    = 4.0*!dpi*1.0d-7
 jpar  = (YkmBFns ## $
         (coeffs*(-outNkValues*(outNkValues+1.0))))$
                         /(u0*R)*1.0d-9*1.0d6        ; uAm^{-2}

 window, winnum
 winnum++
 plt_dat, data.theta * !radeg, data.phi * !radeg, $
          n_elements ( data.theta ), -fit_dBTheta, fit_dBPhi, [1,1], [1,1]

 window, winnum, xSize = 1600, ySize = 300
 winnum++
 !p.background=0
 plot, dBMag
 loadct, 12
 oPlot, fit_dBMag, color = 8*16-1

 r_arr=dblarr(n_elements(data.theta))
 r_arr(*)=r/1.0d3
 geoPack_sphCar, r_arr, data.theta, data.phi, xGEI, yGEI, zGEI, /to_rect   ; from GEI_SPH to GEI_xyz
 geoPack_conv_coord, xGEI, yGEI, zGEI, xGEO, yGEO, zGEO, /from_gei, /to_geo, epoch = t_arr  ; from GEI to GEO
 geoPack_sphCar, xGEO, yGEO, zGEO, Ra_g, theta_g, phi_g, /to_sphere  ; from GEO_xyz to GEO_sph
 yrr=2000
 aacgm_load_coef,yrr
 r_aacgm=Ra_g-Re/1.0d3
 aacgm_conv_coord,theta_g*!radeg,phi_g*!radeg,r_aacgm,mlat,mlon_a,err,/to_aacgm
 id=where(mlon_a lt 0.)
 mlon_a(id)=mlon_a(id)+360.

stop
 loadct, 13, file = 'd:\cwac\davidg\idl\save\davect.tbl'
 window, winnum,title='GEI FAC'
 winnum++
 map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( capSize + 2 ), 0, 90, 360 ], /noborder, /advance
 jLevels = ( fIndGen ( 21 ) - 10 ) * 0.1;2.0
 colors  = bytScl ( jLevels, top = 253 ) + 1
 oldbck=!p.background
 !p.background=0
 jpar=reform(jpar)
 contour, jPar, data.phi*!radeg+180., 90.0-data.theta*!radeg, c_labels=fltarr(n_elements(jLevels))+1,$
           /irreg, /over, levels = jLevels, $
           c_colors = colors, /fill
; tmpVar=drawgrid()
stop
 window, winnum,title='AACGM FAC'
 winnum++
 map_set, 90, 0, 0, /ortho, /iso, $
     limit = [ 90.0 - ( capSize + 2 ), 0, 90, 360 ], /noborder, /advance
 jLevels = ( fIndGen ( 21 ) - 10 ) * 0.1;2.0
 colors  = bytScl ( jLevels, top = 253 ) + 1
 contour, jPar, mlon_a, 90.0-mlat, c_labels=fltarr(n_elements(jLevels))+1,$
           /irreg, /over, levels = jLevels, $
           c_colors = colors, /fill
; tmpVar=drawgrid()
 !p.background=oldbck
 stop
end

