; Code to read and plot AMPERE DATA
; CLW June 2009
;
Pro doy_proc,yr,mnth,dy,dofyr
 d2date=0             ; 1=convert day_of_year to day/month 0=convert day/moth/year to day_of_year
 dMnth=IntArr(12)
 dMnth(0)=31       ; Jan
 dMnth(1)=28       ; Feb
 dMnth(2)=31       ; Mar
 dMnth(3)=30       ; Apr
 dMnth(4)=31       ; May
 dMnth(5)=30       ; Jun
 dMnth(6)=31       ; Jul
 dMnth(7)=31       ; Aug
 dMnth(8)=30       ; Sep
 dMnth(9)=31       ; Oct
 dMnth(10)=30      ; Nov
 dMnth(11)=31      ; Dec
 If (Yr mod 400) eq 0 Then dMnth(1)=29
 If (d2date eq 0) Then $
 Begin
  dofyr=0
  For ii=0,Mnth-2 do dofyr=dofyr+dMnth(ii)
  dofyr=dofyr+Dy
 end
 If (d2date eq 1) Then $
 Begin
  tmp=dofyr & kk=0
  While tmp gt 0 do $
  Begin
   tmp=tmp-dMnth(kk)
   kk=kk+1
  end
  tmp=0
  Mnth=kk
  For ii=0,kk-2 do tmp=tmp+dMnth(ii)
  dy=dofyr-tmp
;  Print,'Day  Month  Year : ',Dy,' ',sMnth(Mnth-1), Yr
;  Print,'Day of Year = ',dofyr
  end
 end
; ************************************************************************************
;
Pro sun,IYEAR,IDAY,IHOUR,MN,ISEC,GST,SLONG,SRASN,SDEC
;  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
;  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
;
;-------  INPUT PARAMETERS:
;  IYR,IDAY,IHOUR,MN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
;    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
;
;-------  OUTPUT PARAMETERS:
;  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
;  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
;  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
;  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
;
 rad=180.0d0/!dpi
 IF (IYEAR LT 1901 OR IYEAR GT 2099) then RETURN
 FDAY=(IHOUR*3600.0d0+MN*60.0d0+ISEC)/86400.D0
 DJ=365.0d0*(IYEAR-1900.0d0)+(IYEAR-1901.0d0)/4.0d0+IDAY-0.5D0+FDAY
 T=DJ/36525.
; VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
 VL=(279.696678+0.9856473354*DJ) mod 360.D0
; GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
 GST=(279.690983+0.9856473354*DJ+360.*FDAY+180.0) mod 360.D0
 GST=GST/rad
; G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
 G=(358.475845+0.985600267*DJ) mod 360.D0
 G=G/rad
 SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
 IF (SLONG GT 2.0d0*!dpi) then SLONG=SLONG-2.0d0*!dpi
 IF (SLONG LT 0.0d0) then SLONG=SLONG+2.0d0*!dpi
 OBLIQ=(23.45229-0.0130125*T)/RAD
 SOB=SIN(OBLIQ)
 SLP=SLONG-9.924E-5
 SIND=SOB*SIN(SLP)
 COSD=SQRT(1.-SIND*SIND)
 SC=SIND/COSD
 SDEC=ATAN(SC)
 SRASN=3.141592654-ATAN(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
end

; **************************************************************************
;
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
;
; ********************************************************************
;
;pro clw_ampere,fname,data,hemis
pro clw_ampere
; fname : Input data file name
; data  : returned structure
; hemis : 0=Nth, 1=Sth
;
; Read Ampere data file. Data coords are Geographic
;
 fname='20080105_a_RevA.dat'
 hemis=0
 st_hr=8.0
 st_mn=0.0
 en_hr=11.0
 en_mn=0.0
;
 header=2                          ; 2 lines of data info
 ndat=file_lines(fname)-header     ; number of data points in the file
;
; Data structure matches data file contents
; OrbitPlane, Satellite, t_sec, radius(km), Latitude(deg), Longitude(deg), dB_xyz(nT)
; position_xyz(km), velocity_xyz(km/s)
;
 data0_struct={ipln:0,sv:0,sday:0.0d0,rad:0.0d0,lat:0.0d0,lon:0.0d0,dbx:0.0d0,dby:0.0d0,dbz:0.0d0,$
   px:0.0d0,py:0.0d0,pz:0.0d0,vx:0.0d0,vy:0.0d0,vz:0.0d0}
 data0=replicate(data0_struct,ndat)

 openr,runit,fname,/get_lun
 dstr=''
 readf,runit,dstr
 reads,strmid(dstr,10,8),year,month,day,format='(i4,i2,i2)'  ; Get year, month, day for time calc later
 readf,runit,dstr
 readf,runit,data0        ; read whole data file
 close,runit
 free_lun,runit

; Select on hemis and time
 st_sc=st_hr*3600.0+st_mn*60.0
 en_sc=en_hr*3600.0+en_mn*60.0
; Get Lat/Lon from X,Y,Z - Lat/Lon in file is suspect
 r_sph=sqrt(data0.px^2+data0.py^2+data0.pz^2)
 th_sph=acos(data0.pz/r_sph)
;  ph_sph=atan(data0.py,data0.px)
 th_deg=th_sph*180./!pi              ; CoLat in degrees

 if hemis eq 0 then np_i=where((th_deg lt 90.0) and (data0.sday ge st_sc) and (data0.sday le en_sc))  ; Nth hemis
 if hemis eq 1 then np_i=where((th_deg gt 90.0) and (data0.sday ge st_sc) and (data0.sday le en_sc))  ; Sth hemis
 np_sel=n_elements(np_i)

; Make a subset of data0 structure
 data_struct={ipln:0,sv:0,sday:0.0d0,rad:0.0d0,lat:0.0d0,lon:0.0d0,dbx:0.0d0,dby:0.0d0,dbz:0.0d0,$
   px:0.0d0,py:0.0d0,pz:0.0d0,vx:0.0d0,vy:0.0d0,vz:0.0d0}
 data=replicate(data_struct,np_sel)
 struct_assign,data0(np_i),data

; Cast data into GEI
 x_gei=dblarr(np_sel)
 y_gei=dblarr(np_sel)
 z_gei=dblarr(np_sel)
 bx_gei=dblarr(np_sel)
 by_gei=dblarr(np_sel)
 bz_gei=dblarr(np_sel)
 doy_proc,year,month,day,dofyr    ; get DOY
 For ii=0,np_sel-1 do begin
  ut=data(ii).sday/3600.0         ; get time off data file
  doy=dofyr
  if ut gt 24. then begin
   ut=ut-24.
   doy++
  end
  ut_hr=fix(ut)
  ut_mn=fix((ut-ut_hr)*60.0)
  ut_sc=fix(data(ii).sday-ut_hr*3600.0-ut_mn*60.0)
  sun,year,doy,ut_hr,ut_mn,ut_sc,GST,SLONG,SRASN,SDEC       ; figure out where the Sun is
  tr_geo2gei=[[ cos(gst),-sin(gst), 0.0],$                  ; rotate into GEI : Checked with Fortran geopack_2005 - Ok
              [ sin(gst),cos(gst), 0.0],$
              [ 0.0,       0.0,    1.0]]
  xyz_gei=[data(ii).px, data(ii).py, data(ii).pz]
  p_gei=reform(tr_geo2gei##xyz_gei)               ; remove dimensions of 1
  x_gei(ii)=p_gei(0)
  y_gei(ii)=p_gei(1)
  z_gei(ii)=p_gei(2)

; This is suspect  - need to check - I think its wrong ************
  b_xyz_gei=[data(ii).dbx,data(ii).dby,data(ii).dbz]
  b_gei=reform(tr_geo2gei##b_xyz_gei)               ; remove dimensions of 1
  bx_gei(ii)=b_gei(0)
  by_gei(ii)=b_gei(1)
  bz_gei(ii)=b_gei(2)
 end

; Get Lat/Lon from X,Y,Z - Lat/Lon in file is suspect
 r_sph=sqrt(x_gei^2+y_gei^2+z_gei^2)
 th_sph=acos(z_gei/r_sph)
 ph_sph=atan(y_gei,x_gei)
 th_deg=th_sph*180./!pi              ; CoLat in degrees
 data.lat=th_sph
 data.lon=ph_sph

; Transformation matrix from X,Y,Z to spherical
; r = sin_th*cos_ph x + sin_th*sin_ph y + cos_th z
;th = cos_th*cos_ph x + cos_th*sin_ph y - sin_th z
;ph =       -sin_ph x +        cos_ph y +      0 z
; IDL does Column, Row
; Fix this later....
 IrNCmp=dblarr(np_sel)
 IrECmp=dblarr(np_sel)
 For ii=0,np_sel-1 do begin
;  tr_a=[[sin(data(ii).lat)*cos(data(ii).lon),sin(data(ii).lat)*sin(data(ii).lon), cos(data(ii).lat)],$
;        [cos(data(ii).lat)*cos(data(ii).lon),cos(data(ii).lat)*sin(data(ii).lon),-sin(data(ii).lat)],$
;        [               -sin(data(ii).lon),                cos(data(ii).lon), 0.0]]
; b_xyz=[data(ii).dbx,data(ii).dby,data(ii).dbz]

; Try routine From geopack
 RHO2=x_gei(ii)^2+y_gei(ii)^2
 R=SQRT(RHO2+z_gei(ii)^2)
 RHO=SQRT(RHO2)
 IF (RHO NE 0.0) THEN begin
  CPHI=x_gei(ii)/RHO
  SPHI=y_gei(ii)/RHO
 endif else begin
  CPHI=1.
  SPHI=0.
 ENDelse
 CT=z_gei(ii)/R
 ST=RHO/R
; BR=(x_gei(ii)*bx_gei(ii) + y_gei(ii)*by_gei(ii) + z_gei(ii)*bz_gei(ii))/R
 IrNCmp(ii)=(bx_gei(ii)*CPHI + by_gei(ii)*SPHI)*CT-bz_gei(ii)*ST
 IrECmp(ii)= by_gei(ii)*CPHI - bx_gei(ii)*SPHI

; b_sph=reform(tr_a##b_xyz)               ; remove dimensions of 1
; IrNCmp(ii)=b_sph(1)
; IrECmp(ii)=b_sph(2)
;stop
 end

 set_plot,'X'
 device,decomposed=0
 window,0,xsize=800,ysize=800,title='Irid Data'
 LoadCT,0
 xm=[1,1]
 ym=[1,1]
; plt_dat,IrCLat,IrMLT,nirid,IrNCmp,IrECmp
 plt_dat,data.lat*180./!pi,data.lon*180./!pi,np_sel,IrNCmp,IrECmp,xm,ym

 Print,'Finished'
end
