; Plot of dB vectors
; Sub_code to plot arrow vectors
;
;	doing this without a map_set is terribly incovenient!
;	please fix to use a map_set so we can overlay things, 
;	use a real grid, etc ...
;


Pro plt_dat,IrCLat,IrMLT,nirid,IrNCmp,IrECmp,xm,ym, $
		title = title

ncol    = 256
col     = 0
IrCol   = col
nlatlab = 5
CLatMin = 50.
mmg     = 7.0
DMag    = 400.
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
    /isotropic,$
    xmargin=[xm(0),xm(1)],$
    ymargin=[ym(0),ym(1)],$
    xstyle=5,ystyle=5, $
    background=255,$
    color=0,$
    linestyle=1, $
    title = title

for i=0,nlatlab-2 do $
    oPlot,x(i,*),y(i,*),$
        linestyle=1,$
        color=0

loadct, 1

for ii=0,nirid-1 do begin

    SLon=IrMLT(ii);+180.                ; Shifted coords Lon,Lat
    SLat=IrCLat(ii)
    ang=(SLon-90.)*!pi/180.
    Xs=SLat*cos(ang)              ; Shifted coords, X,Y
    Ys=SLat*sin(ang)

;   do N shift (N dir unit vector)

    Mg1=sqrt(Xs*Xs+Ys*Ys)
    xN=Xs/Mg1
    yN=Ys/Mg1

;   do E shift (anticlk 90 deg rotation)

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

    if (sqrt(x0*x0+y0*y0) le CLatMin) then begin

        mgn=sqrt((x0-x1)^2+(y0-y1)^2)

        irCol   = 255 - ( bytScl ( mgn, top = 253, min = 0, max = 10 ) + 1 )
        plots,[x0,x1],[y0,y1],$
            color=IrCol,$
            thick=1

        ;; Now do arrow heads

        ;if (mgn gt 0.) then begin

        ;    xu=(x0-x1)/mgn
        ;    yu=(y0-y1)/mgn
        ;    x2=0.4*xu-0.4*yu
        ;    y2=0.4*xu+0.4*yu

        ;    plots,[x1,x1+x2],[y1,y1+y2],$
        ;        color=IrCol,$
        ;        thick=1

        ;    x2=0.4*xu+0.4*yu
        ;    y2=-0.4*xu+0.4*yu

        ;    plots,[x1,x1+x2],[y1,y1+y2],$
        ;        color=IrCol,$
        ;        thick=1

        ;endif

    endif   ; If within Lat Range

endfor     ; NIrid Loop

end     
