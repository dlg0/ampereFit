function gc_dist,lat_s,lat_f,lon_s,lon_f
; Great circle distance calc
;
; C. L. Waters
; Centre for Space Physics
; University of Newcastle
; Sept 2010
;
; Input : lat_s : start latitude (deg) - float
;         lat_f : array of end latitudes (deg)
;         lon_s : start longitude (deg) - float
;         lon_f : array of end longitudes )deg)
;
  re=6371.d0 + 780.d0
  phi_s=lat_s*!dpi/180.0
  phi_f=lat_f*!dpi/180.0
  lmda_s=lon_s*!dpi/180.0
  lmda_f=lon_f*!dpi/180.0
  d_lmda=lmda_f-lmda_s
  
  tp=sqrt((cos(phi_f)*sin(d_lmda))^2 + (cos(phi_s)*sin(phi_f)-sin(phi_s)*cos(phi_f)*cos(d_lmda))^2)
  bt=sin(phi_s)*sin(phi_f) + cos(phi_s)*cos(phi_f)*cos(d_lmda)
  ang=atan(tp,bt)
  dis=re*ang
  return,dis
;  Print,'Dist = ',dis
end