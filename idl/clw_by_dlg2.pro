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


pro clw_by_dlg2

 winnum=1
	fileName = '20080105_a_RevA.dat'
	sHr = 12
	eHr	= 15
	read_ampere_dlg2, fileName, sHr, eHr, data

 device,decomposed=0
 window,winnum,xsize=800,ysize=800,title='dgIrid Data'
 LoadCT,0
	device, decomposed = 0
   	plt_dat, data.theta * !radeg, data.phi * !radeg, $
			n_elements ( data.theta ), -data.dbTheta, data.dbPhi, $
			[1,1], [1,1]
stop
end
