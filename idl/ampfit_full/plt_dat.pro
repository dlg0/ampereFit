; Plot of dB vectors
; Sub_code to plot arrow vectors
; C. L Waters
;
; full sphere version - plots one hemisphere
;
pro plt_dat,inCLat,IrMLT,nirid,IrNCmp,IrECmp, south, $
		title = title, $
        satNu = satNu, $
        arrowHeads = arrowHeads, $
        capSize = cLatMin

	xm = [1,1]
	ym = [1,1]

 IrCLat=inCLat                 ; get a copy of coLat data
 InNCmp=IrNCmp
 If south eq 1 then begin
  IrCLat=180.0-inCLat          ; mirror to Nth hemisphere perspective
  InNCmp=-IrNCmp               ; +ve is Nth
 end
 if keyword_set ( satNu ) then begin
    allowedColors   = [ 0,1,3,7,8]
    presentColorIndex    = 0
    satColors   = satNu * 0
    satMatch    = [[presentColorIndex],[satNu[0]]]
    for i = 0, n_elements ( satNu ) - 1 do begin

        iiMatch = where ( satMatch[*,1] eq satNu[i], iiMatchCnt )
        if iiMatchCnt eq 0 then begin

            ++presentColorIndex
            satMatch = [[satMatch[*,0], presentColorIndex mod 5], $
                        [satMatch[*,1], satNu[i] ] ]
            satColors[i]    = allowedColors[presentColorIndex mod 5]
        endif else begin
            satColors[i]    = allowedColors[satMatch[iiMatch,0]]
        endelse
    endfor
endif

ncol    = 256
col     = 0
IrCol   = col
if not keyword_set ( cLatMin ) then cLatMin = 50.0

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

loadct, 1,/silent
if keyword_set ( satNu ) then loadct, 12,/silent
db_Mag=fltarr(nirid)

for ii=0,nirid-1 do begin

    SLon=IrMLT(ii);+180.                ; Shifted coords Lon,Lat
    SLat=IrCLat(ii)
    ang=(SLon-90.0)*!pi/180.

    Xs=SLat*cos(ang)              ; Shifted coords, X,Y
    Ys=SLat*sin(ang)

;   do N shift (N dir unit vector)

    Mg1=sqrt(Xs*Xs+Ys*Ys)
    xN=Xs/Mg1
    yN=Ys/Mg1

;   do E shift (anticlk 90 deg rotation)

    xE=-yN
    yE=xN
    NVal=mmg*InNCmp(ii)/DMag  ; shifted coords B Val in theta,phi coords
    EVal=mmg*IrECmp(ii)/DMag
    dB_mag(ii)=sqrt(IrNCmp(ii)^2+IrECmp(ii)^2)
    x1=Xs-NVal*xN+EVal*xE
    y1=Ys-NVal*yN+EVal*yE
    x0=Xs-YSh                         ; Convert to geomag but plot is rotated 90 Deg
    x1=x1-YSh
    y0=Ys+XSh
    y1=y1+XSh

    if (sqrt(x0*x0+y0*y0) le CLatMin) then begin

        mgn=sqrt((x0-x1)^2+(y0-y1)^2)

        if keyword_set ( satNu ) then begin
            loadct, satColors[ii], /silent
            irCol   = 255 - ( bytScl (mgn, top = 200, min = 0, max = 10 ) + 1 )
        endif else begin
            irCol   = 255 - ( bytScl (mgn, top = 200, min = 0, max = 10 ) + 1 )
        endelse

        xyouts,x0,y0,'.',color=0     ; added to show data locations
        plots,[x0,x1],[y0,y1],$
            color=IrCol,$
            thick=1

        ; Now do arrow heads

        if (mgn gt 0. and keyword_set ( arrowHeads ) ) then begin

            xu=(x0-x1)/mgn
            yu=(y0-y1)/mgn
            x2=0.4*xu-0.4*yu
            y2=0.4*xu+0.4*yu

            plots,[x1,x1+x2],[y1,y1+y2],$
                color=IrCol,$
                thick=1

            x2=0.4*xu+0.4*yu
            y2=-0.4*xu+0.4*yu

            plots,[x1,x1+x2],[y1,y1+y2],$
                color=IrCol,$
                thick=1

        endif

    endif   ; If within Lat Range
endfor     ; NIrid Loop
;Print,'Min:Max DB_Mag = ',min(dB_mag),'  ',max(dB_mag)
end
