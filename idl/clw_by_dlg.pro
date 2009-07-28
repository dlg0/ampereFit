; ick
Pro plt_dat,IrCLat,IrMLT,nirid,IrNCmp,IrECmp,xm,ym
 ncol=256
 col=ncol-1
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
 Plot,x(nlatlab-1,*),y(nlatlab-1,*),/isotropic,xmargin=[xm(0),xm(1)],ymargin=[ym(0),ym(1)],xstyle=5,ystyle=5,$
   background=0,color=ncol-1,linestyle=1
 For i=0,nlatlab-2 do OPlot,x(i,*),y(i,*),linestyle=1,color=col
 For ii=0,nirid-1 do $
 Begin
  SLon=IrMLT(ii)+180.0                ; Shifted coords Lon,Lat
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


pro clw_by_dlg

	fileName = '20080105_a_RevA.dat'
	sHr = 14 
	eHr	= 16
	capSize	= 50

	read_ampere_dlg, fileName, sHr, eHr, capSize, data

	device, decomposed = 0
	window, 0
   	plt_dat, data.theta * !radeg, data.phi * !radeg, $
			n_elements ( data.theta ), -data.dbTheta, data.dbPhi, $
			[1,1], [1,1]

	kMax	= 16 
	mMax	= 5
	schaBasisFunctions, kMax, mMax, capSize, data.theta, data.phi, $
		YkmBFns=YkmBFns, dYkmDthBFns=dYkmDthBFns, dYkmDphBFns=dYkmDphBFns,$
		OUTNKVALUES=outNkValues, OUTKVALUES=outKValues, OUTMVALUES=outMValues;, /oddSet

;	Fit |dB| to dB.grad Ykm in the Iridium system

	dBMag	= sqrt(data.dBTheta^2+data.dBPhi^2)
	
	kTh	= transpose(rebin($
		data.dBTheta/dBMag,$
		n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))
	kPh	= transpose(rebin($
		data.dBPhi/dBMag,$
		n_elements(data.dBTheta),n_elements(dYkmDthBFns[*,0])))
	
	bFuncs	= kTh*dYkmDphBFns + $
		kPh*dYkmDthBFns

	alpha	= transpose(bFuncs) ## bFuncs 
	beta_	= transpose( transpose(bFuncs) ## dBMag)

	la_svd, alpha, w, u, v, status=svdStatus

	kk	= where(w lt max ( w ) * 1d-5, kkcnt)
	if kkcnt gt 0 then w[kk]=0
		
	coeffs	= svsol(u, w, v, beta_)

	fit_dBTheta	= dYkmDPhBFns ## coeffs
	fit_dBPhi	= dYkmDThBFns ## coeffs
	fit_dBMag	= bFuncs ## coeffs

	R	= 1
	u0	= 1
	jpar	= (YkmBFns ## $
		(coeffs*(-outNkValues*(outNkValues+1.0))))$
				/(u0*R)*1e-9*1e6	; uAm^{-2}

	window, 1
	plt_dat, data.theta * !radeg, data.phi * !radeg, $
			n_elements ( data.theta ), -fit_dBTheta, fit_dBPhi, $
			[1,1], [1,1]

	window, 2, xSize = 1600, ySize = 300
	plot, dBMag
	loadct, 12
	oPlot, fit_dBMag, color = 8*16-1

	window, 3
	loadct, 13, file = 'davect.tbl'

	map_set, 90, 0, 0, /ortho, /iso, $
		limit = [ 90.0 - ( capSize + 2 ), 0, 90, 360 ], /noborder, $
		/advance;

	jLevels	= ( fIndGen ( 21 ) - 10 ) * 2.0
	colors	= bytScl ( jLevels, top = 253 ) + 1
	contour, jPar, data.phi*!radeg, 90.0-data.theta*!radeg, $
			/irreg, $
			/over, $
			levels = jLevels, $
			c_colors = colors, $
			/fill
stop
end
