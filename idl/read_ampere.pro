pro read_ampere,fname,data

  @daves_paths

	device, decomposed = 0
	loadct, 0
	!p.background = 255

; read ampere file
  header=2
  ndat=file_lines(fname)-header
  data0_struct={ipln:0,sv:0,sday:0.0d0,rad:0.0d0,lat:0.0d0,lon:0.0d0,dbx:0.0d0,dby:0.0d0,dbz:0.0d0,$
    px:0.0d0,py:0.0d0,pz:0.0d0,vx:0.0d0,vy:0.0d0,vz:0.0d0}
  data0=replicate(data0_struct,ndat)

  openr,runit,fname,/get_lun

  dstr=''
  readf,runit,dstr
  reads,strmid(dstr,10,8),year,month,day,format='(i4,i2,i2)'
  readf,runit,dstr

  readf,runit,data0

  close,runit
  free_lun,runit

; extend data structure  
  data_struct={ipln:0,sv:0,sday:0.0d0,rad:0.0d0,lat:0.0d0,lon:0.0d0,dbx:0.0d0,dby:0.0d0,dbz:0.0d0,$
    px:0.0d0,py:0.0d0,pz:0.0d0,vx:0.0d0,vy:0.0d0,vz:0.0d0,epoch:0.0d0,alat:0.0d0,alon:0.0d0,hrmlt:0.0d0,$
    xcrd:0.0d0,ycrd:0.0d0,ba:0.0d0,nxa:0.0d0,nya:0.0d0,bc:0.0d0,nxc:0.0d0,nyc:0.0d0}
  data=replicate(data_struct,n_elements(data0))
  struct_assign,data0,data

;; eliminate bad point
;  idx=where((data.lat ge -90) and (data.lat le 90) and (data.lon ge -180) and (data.lon le 360))
;  data=data[idx]
  
; change time to epoch
  cdf_epoch,epoch0,year,month,day,/compute_epoch
  data.epoch=data.sday*1.0d3+epoch0
  
; compute aacgm latitude and mlt
  ;aacgm_set_path,aacgm_path,/quiet
  aacgm_load_coef,(round(year/5.)*5)<2000;,/quiet
  geopack_sphcar,data.px,data.py,data.pz,zrad,zclt,zlon,/to_sphere,/degree
  aacgm_conv_coord,90.0-zclt,zlon,zrad-6371.2,alat,alon,error,/to_aacgm
  cdf_epoch,epoch0,year,1,1,/compute_epoch
  t0=(data.epoch-epoch0)*1.0d-3
  hrmlt=aacgm_mlt(fltarr(n_elements(t0))+year,t0,alon)
  data.alat=alat
  data.alon=alon
  data.hrmlt=hrmlt

; compute iridium coordinates
  lon=hrmlt*15
  deg2rad=(lon-180)*!pi/180
  xcrd=(90.0-abs(alat))*cos(deg2rad)
  ycrd=(90.0-abs(alat))*sin(deg2rad)
  data.xcrd=xcrd
  data.ycrd=ycrd
  

; compute along/cross-track magnetic field
; compute sv velocity unit vector
  vmag=sqrt(data.vx^2+data.vy^2+data.vz^2)
  uvx=data.vx/vmag
  uvy=data.vy/vmag
  uvz=data.vz/vmag
;compute nadir unit vector
  pmag=sqrt(data.px^2+data.py^2+data.pz^2)
  upx=-data.px/pmag
  upy=-data.py/pmag
  upz=-data.pz/pmag  
; compute cross-track unit vector and renormalize
  ucx=upy*uvz-upz*uvy                                                                                                                                                                                 
  ucy=upz*uvx-upx*uvz
  ucz=upx*uvy-upy*uvx
  cmag=sqrt(ucx^2+ucy^2+ucz^2)
  ucx=ucx/cmag
  ucy=ucy/cmag
  ucz=ucz/cmag
; compute along-track unit vector and renormalize 
  uax=ucy*upz-ucz*upy                                                                                                                                                                                 
  uay=ucz*upx-ucx*upz
  uaz=ucx*upy-ucy*upx
  amag=sqrt(uax^2+uay^2+uaz^2)
  uax=uax/amag
  uay=uay/amag
  uaz=uaz/amag
; compute magnetic field vector in along/cross-track and radial direction
  ba=data.dbx*uax+data.dby*uay+data.dbz*uaz
  bc=data.dbx*ucx+data.dby*ucy+data.dbz*ucz
  br=data.dbx*upx+data.dby*upy+data.dbz*upz                                                                              
  data.ba=ba
  data.bc=bc


; compute along/cross-track unit vectors in iridium coordinates
  px0=data.px
  py0=data.py
  pz0=data.pz
  geopack_sphcar,px0,py0,pz0,pr0,pt0,pp0,/to_sphere,/degree
  aacgm_conv_coord,90.0-pt0,pp0,pr0-6371.2,mlat0,mlon0,error0,/to_aacgm
  cdf_epoch,epoch0,year,1,1,/compute_epoch
  t0=(data.epoch-epoch0)*1.0d-3
  mlt0=aacgm_mlt(fltarr(n_elements(t0))+year,t0,mlon0)
  mltlon0=mlt0*15
  deg2rad=(mltlon0-180)*!dtor
  xcrd0=(90.0-abs(mlat0))*cos(deg2rad)
  ycrd0=(90.0-abs(mlat0))*sin(deg2rad)

  px1=data.px+uax
  py1=data.py+uay
  pz1=data.pz+uaz
  geopack_sphcar,px1,py1,pz1,pr1,pt1,pp1,/to_sphere,/degree
  aacgm_conv_coord,90.0-pt1,pp1,pr1-6371.2,mlat1,mlon1,error1,/to_aacgm
  cdf_epoch,epoch0,year,1,1,/compute_epoch
  t0=(data.epoch-epoch0)*1.0d-3
  mlt1=aacgm_mlt(fltarr(n_elements(t0))+year,t0,mlon1)
  mltlon1=mlt1*15
  deg2rad=(mltlon1-180)*!dtor
  xcrd1=(90.0-abs(mlat1))*cos(deg2rad)
  ycrd1=(90.0-abs(mlat1))*sin(deg2rad)

  px2=data.px+ucx
  py2=data.py+ucy
  pz2=data.pz+ucz
  geopack_sphcar,px2,py2,pz2,pr2,pt2,pp2,/to_sphere,/degree
  aacgm_conv_coord,90.0-pt2,pp2,pr2-6371.2,mlat2,mlon2,error2,/to_aacgm
  cdf_epoch,epoch0,year,1,1,/compute_epoch
  t0=(data.epoch-epoch0)*1.0d-3
  mlt2=aacgm_mlt(fltarr(n_elements(t0))+year,t0,mlon2)
  mltlon2=mlt2*15
  deg2rad=(mltlon2-180)*!dtor
  xcrd2=(90.0-abs(mlat2))*cos(deg2rad)
  ycrd2=(90.0-abs(mlat2))*sin(deg2rad)

  nxa=xcrd1-xcrd0
  nya=ycrd1-ycrd0
  nmag=sqrt(nxa^2+nya^2)
  nxa=nxa/nmag
  nya=nya/nmag  

  nxc=xcrd2-xcrd0
  nyc=ycrd2-ycrd0
  nmag=sqrt(nxc^2+nyc^2)
  nxc=nxc/nmag
  nyc=nyc/nmag  

  data.nxa=nxa
  data.nya=nya
  data.nxc=nxc
  data.nyc=nyc

return
end
