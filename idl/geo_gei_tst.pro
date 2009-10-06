; Code to test GEO to GEI vector translation
; CLW Aug 2009
;
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
 VL=(279.696678+0.9856473354*DJ) mod 360.D0
 GST=(279.690983+0.9856473354*DJ+360.*FDAY+180.0) mod 360.D0
 GST=GST/rad
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
pro geo_gei_tst
 winnum=0

; select a date and UT time
 year=2005
 month=1
 day=1
 doy=1

 ut_hr=12
 ut_mn=0
 ut_sc=0
 np_sel=100        ; Num points to generate
;
; Generate some locations in GEO coords
 data_struc={dbx:0.0d0,dby:0.0d0,dbz:0.0d0,px:0.0d0,py:0.0d0,pz:0.0d0}
 data=replicate(data_struc,np_sel)
 data.px=30.*(randomu(s1,np_sel)-0.5)
 data.py=30.*(randomu(s2,np_sel)-0.5)
 data.pz=30.*(randomu(s3,np_sel)-0.5)

; Generate some dB values in GEO coords
 data.dbx=50.*(randomu(s4,np_sel)-0.5)
 data.dby=50.*(randomu(s5,np_sel)-0.5)
 data.dbz=50.*(randomu(s6,np_sel)-0.5)

; Cast data into GEI
 x_gei=dblarr(np_sel)
 y_gei=dblarr(np_sel)
 z_gei=dblarr(np_sel)
 bx_gei=dblarr(np_sel)
 by_gei=dblarr(np_sel)
 bz_gei=dblarr(np_sel)
 sun,year,doy,ut_hr,ut_mn,ut_sc,GST,SLONG,SRASN,SDEC       ; figure out where the Sun is

; GEO to GEI transformation matrix (see Kivelson and Russell book)
 tr_geo2gei=[[ cos(gst),-sin(gst), 0.0],$                  ; rotate into GEI : Checked with Fortran geopack_2005 - Ok
             [ sin(gst),cos(gst), 0.0],$
             [ 0.0,       0.0,    1.0]]

 For ii=0,np_sel-1 do begin
; Transform the GEO locations - should be fine
  xyz_gei=[data(ii).px, data(ii).py, data(ii).pz]
  p_gei=reform(tr_geo2gei##xyz_gei)               ; remove dimensions of 1
  x_gei(ii)=p_gei(0)
  y_gei(ii)=p_gei(1)
  z_gei(ii)=p_gei(2)

; Transform the dB vectors, should work since the z axes of GEo and GEI are common
  b_xyz_gei=[data(ii).dbx,data(ii).dby,data(ii).dbz]
  b_gei=reform(tr_geo2gei##b_xyz_gei)               ; remove dimensions of 1
  bx_gei(ii)=b_gei(0)
  by_gei(ii)=b_gei(1)
  bz_gei(ii)=b_gei(2)
 end

; DLG version
; Create rotation matrix from xyz to xyzGEI
 dx = 1.
 dy = 1.
 dz = 1.
 iiCnt=np_sel
 sday=dblarr(np_sel)
 sday(0:np_sel-1)=doy
 cdf_epoch, epoch0, year, month, day, /compute_epoch
 epoch = sday * 1d3 + epoch0

 geoPack_reCalc, year, month, day, /date
 geoPack_conv_coord, (data.px), (data.py), (data.pz), $ ; units of km
   xGEI, yGEI, zGEI, /from_geo, /to_gei, epoch = epoch

geoPack_conv_coord, (data.px+dx), (data.py), (data.pz), $
xGEI_x, yGEI_x, zGEI_x, $
/from_geo, /to_gei, $
epoch = epoch
geoPack_conv_coord, (data.px), (data.py+dy), (data.pz), $
xGEI_y, yGEI_y, zGEI_y, $
/from_geo, /to_gei, $
epoch = epoch
geoPack_conv_coord, (data.px), (data.py), (data.pz+dz), $
xGEI_z, yGEI_z, zGEI_z, $
/from_geo, /to_gei, $
epoch = epoch

rotMat = fltArr ( 3, 3, iiCnt )
rotMat[0,0,*] = xGEI-xGEI_x
rotMat[1,0,*] = xGEI-xGEI_y
rotMat[2,0,*] = xGEI-xGEI_z

rotMat[0,1,*] = yGEI-yGEI_x
rotMat[1,1,*] = yGEI-yGEI_y
rotMat[2,1,*] = yGEI-yGEI_z

rotMat[0,2,*] = zGEI-zGEI_x
rotMat[1,2,*] = zGEI-zGEI_y
rotMat[2,2,*] = zGEI-zGEI_z

rotMat = -rotMat

; Rotate db vectors to GEI
dbGEI = rotMat ## [ $
[ (data.dbx) ],$
[ (data.dby) ],$
[ (data.dbz) ] ]

 set_plot,'win'
 device,decomposed=0
 window,winnum,xsize=900,ysize=900,title='Irid Data'
 LoadCT,0
 !p.multi=[0,2,3,0]
 !p.charsize=1.8
 plot,x_gei,xgei,title='X',xtitle='CLW',ytitle='DLG'
 plot,y_gei,ygei,title='Y',xtitle='CLW',ytitle='DLG'
 plot,z_gei,zgei,title='Z',xtitle='CLW',ytitle='DLG'

 plot,bx_gei,dbGEI(*,0),title='dbx',xtitle='CLW',ytitle='DLG'
 plot,by_gei,dbGEI(*,1),title='dby',xtitle='CLW',ytitle='DLG'
 plot,bz_gei,dbGEI(*,2),title='dbz',xtitle='CLW',ytitle='DLG'

 Print,'Finished'
end
