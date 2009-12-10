pro gei_to_aacgm, gei_R_km, gei_coLat_rad, gei_lon_rad, $
 aacgm_R_km, aacgm_coLat_rad, aacgm_lon_rad, $
 year = year, epoch = epoch

 if (not keyword_set(year)) then year = 2000
 aacgm_load_coef, year<2000 ; once we have newer coeffs update this
 Re_km = 6357.0
 R_km = 110.0 + Re_km

; conv GEI r,thet,phi -> GEI XYZ
 geoPack_sphCar, gei_R_km, gei_coLat_rad, gei_lon_rad, $
                 gei_x, gei_y, gei_z, /to_rect

; conv GEI XYZ -> GEOG XYZ
 geoPack_conv_coord, gei_x, gei_y, gei_z, $
                     geo_x, geo_y, geo_z, $
                     /from_gei, /to_geo

 epoch = fltArr(n_elements(gei_x[*])) + epoch

; conv GEOG XYZ -> GEOG r,thet,phi
 geoPack_sphCar, geo_x, geo_y, geo_z, $
                 geog_R_km, geog_coLat_rad, geog_lon_rad, /to_sphere

; conv GEOG r,thet,phi -> AACGM [lat,lon]
 aacgm_conv_coord, 90.0-geog_coLat_rad*!radeg, geog_lon_rad*!radeg, $
                   geog_R_km - Re_km, $
                   aacgm_lat_deg, aacgm_lon_deg, err, /to_aacgm
 aacgm_coLat_deg = 90.0-aacgm_lat_deg
 aacgm_coLat_rad = aacgm_coLat_deg * !dtor
 aacgm_lon_rad = aacgm_lon_deg * !dtor
end


; there are two (actually three) versions of the rotation
; code in here, but only the dot product version (/dot_method)
; is working right now :-(

pro rotate_gei_to_aacgm, gei_R_km, gei_coLat_rad, gei_lon_rad, $
aacgm_R_km, aacgm_coLat_rad, aacgm_lon_rad, $
gei_dbTh, gei_dbPh, $
aacgm_dbTh = aacgm_dbTh, $
aacgm_dbPh = aacgm_dbPh, $
year = year, $
epoch = epoch, $
dot_method = dot_method, $
rot_method = rot_method


if ( (not keyword_set ( dot_method )) and (not keyword_set ( rot_method )) ) then rot_method = 1
if ( not keyword_set ( year ) ) then year = 2000

aacgm_load_coef, year<2000 ; once we have newer coeffs update this
Re_km = 6357.0
R_km = 780.0 + Re_km

if keyword_set ( dot_method ) then begin

; OK, method 2 ... we need the r & theta unit vectors
; in terms of their x and y comonents for both GEI and
; AACGM at each location. Here we use polar coordinates
; for the conversion, not sure how correct that is but
; that's how it's been done in the past.

dR = 0.01

; get the x/y coords

GEI_x = gei_colat_rad * !radeg * cos ( gei_lon_rad )
GEI_y = gei_colat_rad * !radeg * sin ( gei_lon_rad )

; get GEI r & theta unit vectors in terms of x and y components

GEI_x_R = ( gei_colat_rad * !radeg + dR ) * cos ( gei_lon_rad )
GEI_y_R = ( gei_colat_rad * !radeg + dR ) * sin ( gei_lon_rad )

GEI_R_unit_x = ( GEI_x_R - GEI_x ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
GEI_R_unit_y = ( GEI_y_R - GEI_y ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )

; rotate the r vector by 90 degrees to get the thetat vector

GEI_th_unit_x = -GEI_R_unit_y
GEI_th_unit_y = GEI_R_unit_x

; now get the same info for the aacgm r and theta unit vectors,
    ; but there are more steps here since we need to convert the aacgm
; coords to gei. and yes, i know using R and th for the polar coords
; is somewhat stupid in that it is confusion when we also have th/ph
; coords for the rest of the stuff

; small step in aacgm r direction & convert to geo coords

  aacgm_conv_coord, (90.0 - aacgm_coLat_rad * !radeg)+dR, aacgm_lon_rad * !radeg, $
gei_R_km-6357.0, GEO_lat_deg_R, GEO_lon_deg_R, err, /to_geo

; convert geo to gei

  geoPack_sphCar, gei_R_km, ( 90.0 - GEO_lat_deg_R ) * !dtor, GEO_lon_deg_R * !dtor, $
GEO_x_R, GEO_y_R, GEO_z_R, /to_rect
geoPack_conv_coord, GEO_x_R, GEO_y_R, GEO_z_R, $ ; units of km
GEI_x_R, GEI_y_R, GEI_z_R, $
/from_geo, /to_gei, $
epoch = fltArr(n_elements(GEO_x_R)) + epoch
geoPack_sphCar, GEI_x_R, GEI_y_R, GEI_z_R, $
GEI_R_km_R, GEI_colat_rad_R, GEI_lon_rad_R,$
/to_spher

; get the x/y coords of the small step in aacgm r

GEI_x_R = ( GEI_colat_rad_R * !radeg ) * cos ( GEI_lon_rad_R )
GEI_y_R = ( GEI_colat_rad_R * !radeg ) * sin ( GEI_lon_rad_R )

; get aacgm r and thetat unit vector x and y components

AACGM_R_unit_x = ( GEI_x_R - GEI_x ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
AACGM_R_unit_y = ( GEI_y_R - GEI_y ) / sqrt ( (GEI_x_R-GEI_x)^2+ (GEI_y_R-GEI_y)^2 )
AACGM_th_unit_x = -AACGM_R_unit_y
AACGM_th_unit_y = AACGM_R_unit_x

; get the x/y components of db

dbX = gei_dbTh * GEI_R_unit_x + gei_dbPh * GEI_th_unit_x
dbY = gei_dbTh * GEI_R_unit_y + gei_dbPh * GEI_th_unit_y

; get the aacgm r/theta components which are the th/ph components
; that we want

aacgm_dbTh = -1.0 * ( dbX * AACGM_R_unit_x + dbY * AACGM_R_unit_y )
aacgm_dbPh = -1.0 * ( dbX * AACGM_th_unit_x + dbY * AACGM_th_unit_y )

endif

if keyword_set ( rot_method ) then begin

; get GEI xyz coords

geoPack_sphCar, gei_R_km, gei_coLat_rad, gei_lon_rad, $
gei_x, gei_y, gei_z, /to_rect

; check GEI to AACGM conversion

geoPack_conv_coord, gei_x, gei_y, gei_z, $
geo_x, geo_y, geo_z, $
/from_gei, /to_geo;, $
epoch = fltArr ( n_elements ( gei_x[*] ) ) + epoch
geoPack_sphCar, geo_x, geo_y, geo_z, $
geog_R_km, geog_coLat_rad, geog_lon_rad, /to_sphere
  aacgm_conv_coord, 90.0 - geog_coLat_rad * !radeg, geog_lon_rad * !radeg, $
geog_R_km - Re_km, $
aacgm_lat_deg, aacgm_lon_deg, err, /to_aacgm
aacgm_coLat_deg = 90 - aacgm_lat_deg
aacgm_coLat_rad1 = aacgm_coLat_deg * !dtor
aacgm_lon_rad1 = aacgm_lon_deg * !dtor

map_set, 90, 0, 0, /ortho, /iso, $
limit = [ 90.0 - ( 40 + 2 ), 0, 90, 360 ], $
/noborder, /advance, title = 'testing aacgm coords'
plots, aacgm_lon_rad*!radeg, 90-aacgm_coLat_rad*!radeg, $
psym = 4, color = 0
plots, aacgm_lon_rad1*!radeg, 90-aacgm_coLat_rad1*!radeg, $
psym = 4, color = 0
map_grid, label = 1, latDel = 10.0, color = 0

stop

; get AACGM xyz coords

     geoPack_sphCar, aacgm_R_km, aacgm_coLat_rad, aacgm_lon_rad, $
aacgm_x, aacgm_y, aacgm_z, /to_rect

; get GEI xyz vector components

geoPack_bSpCar, gei_coLat_rad, gei_lon_rad, $
gei_dbTh*0, gei_dbTh, gei_dbPh, $
gei_dbx, gei_dby, gei_dbz

; try a direct GEI to AACGM rotation matrix

dTh = 0.1*!dtor
dPh = 0.1*!dtor

gei_to_aacgm, gei_R_km, gei_coLat_rad, gei_lon_rad+dPh, $
aacgm_R_km, aacgm_coLat_rad_ph, aacgm_lon_rad_ph, $
year = year, epoch = epoch
gei_to_aacgm, gei_R_km, gei_coLat_rad+dTh, gei_lon_rad, $
aacgm_R_km, aacgm_coLat_rad_th, aacgm_lon_rad_th, $
year = year, epoch = epoch

iiNeg = where ( aacgm_lon_rad_ph lt 0, iiNegCnt )
if iiNegCnt gt 0 then aacgm_lon_rad_ph[iiNeg] = aacgm_lon_rad_ph[iiNeg] + 2*!pi
iiNeg = where ( aacgm_lon_rad_th lt 0, iiNegCnt )
if iiNegCnt gt 0 then aacgm_lon_rad_th[iiNeg] = aacgm_lon_rad_th[iiNeg] + 2*!pi

rotMat = fltArr ( 2, 2, n_elements ( aacgm_R_km ) )

for i=0, n_elements ( aacgm_R_km[*] )-1 do begin

rotMat[0,0,i] = ( aacgm_coLat_rad[i]-aacgm_coLat_rad_th[i] ) / dTh
rotMat[0,1,i] = ( aacgm_lon_rad[i]-aacgm_lon_rad_th[i] ) / dTh
rotMat[1,0,i] = ( aacgm_coLat_rad[i]-aacgm_coLat_rad_ph[i] ) / dPh
rotMat[1,1,i] = ( aacgm_lon_rad[i]-aacgm_lon_rad_ph[i] ) / dPh

if (aacgm_lon_rad[i]-aacgm_lon_rad_th[i]) lt -1.0*!pi then $
rotMat[0,1,i] = (( aacgm_lon_rad[i]-aacgm_lon_rad_th[i] )+2*!pi) / dTh
if (aacgm_lon_rad[i]-aacgm_lon_rad_th[i]) gt !pi then $
rotMat[0,1,i] = (( aacgm_lon_rad[i]-aacgm_lon_rad_th[i] )+2*!pi) / dTh

if aacgm_lon_rad[i]-aacgm_lon_rad_ph[i] ge 2*!pi then $
rotMat[1,1,i] = (( aacgm_lon_rad[i]-aacgm_lon_rad_ph[i] )-2*!pi) / dPh
;if (aacgm_lon_rad[i]-aacgm_lon_rad_ph[i]) lt 0 then $
; rotMat[1,1,i] = (( aacgm_lon_rad[i]-aacgm_lon_rad_ph[i] )+2*!pi) / dPh

endfor

dbAACGM = rotMat ## [ $
[ gei_dbTh[*] ],$
[ gei_dbPh[*] ] ]
stop
; create rotation matrix from GEI to AACGM

dx = 1.
dy = 1.
dz = 1.

geoPack_conv_coord, aacgm_x+dx, aacgm_y, aacgm_z, $ ; units of km
xGEO_x, yGEO_x, zGEO_x, $
/from_gei, /to_geo, $
epoch = fltArr(n_elements(gei_x)) + epoch
geoPack_conv_coord, aacgm_x, aacgm_y+dy, aacgm_z, $ ; units of km
xGEO_y, yGEO_y, zGEO_y, $
/from_gei, /to_geo, $
epoch = fltArr(n_elements(gei_x)) + epoch
geoPack_conv_coord, aacgm_x, aacgm_y, aacgm_z+dz, $ ; units of km
xGEO_z, yGEO_z, zGEO_z, $
/from_gei, /to_geo, $
epoch = fltArr(n_elements(gei_x)) + epoch

geoPack_sphCar, xGEO_x, yGEO_x, zGEO_x, $
gei_R_km_x, gei_coLat_rad_x, gei_lon_rad_x, /to_spher
geoPack_sphCar, xGEO_y, yGEO_y, zGEO_y, $
gei_R_km_y, gei_coLat_rad_y, gei_lon_rad_y, /to_spher
geoPack_sphCar, xGEO_z, yGEO_z, zGEO_z, $
gei_R_km_z, gei_coLat_rad_z, gei_lon_rad_z, /to_spher

  aacgm_conv_coord, 90.0 - gei_coLat_rad_x * !radeg, gei_lon_rad_x * !radeg, $
gei_R_km_x-Re_km, aacgm_lat_deg_x, aacgm_lon_deg_x, err, /to_aacgm
  aacgm_conv_coord, 90.0 - gei_coLat_rad_y * !radeg, gei_lon_rad_y * !radeg, $
gei_R_km_y-Re_km, aacgm_lat_deg_y, aacgm_lon_deg_y, err, /to_aacgm
   aacgm_conv_coord, 90.0 - gei_coLat_rad_z * !radeg, gei_lon_rad_z * !radeg, $
gei_R_km_z-Re_km, aacgm_lat_deg_z, aacgm_lon_deg_z, err, /to_aacgm

  aacgm_coLat_rad_x = (90.0-aacgm_lat_deg_x) * !dtor
  aacgm_coLat_rad_y = (90.0-aacgm_lat_deg_y) * !dtor
  aacgm_coLat_rad_z = (90.0-aacgm_lat_deg_z) * !dtor

aacgm_lon_rad_x = aacgm_lon_deg_x * !dtor
aacgm_lon_rad_y = aacgm_lon_deg_y * !dtor
aacgm_lon_rad_z = aacgm_lon_deg_z * !dtor

geoPack_sphCar, gei_R_km_x, aacgm_coLat_rad_x, aacgm_lon_rad_x, $
aacgm_x_x, aacgm_y_x, aacgm_z_x, /to_rect
geoPack_sphCar, gei_R_km_y, aacgm_coLat_rad_y, aacgm_lon_rad_y, $
aacgm_x_y, aacgm_y_y, aacgm_z_y, /to_rect
geoPack_sphCar, gei_R_km_z, aacgm_coLat_rad_z, aacgm_lon_rad_z, $
aacgm_x_z, aacgm_y_z, aacgm_z_z, /to_rect

rotMat = fltArr ( 3, 3, n_elements ( aacgm_R_km ) )
rotMat[0,0,*] = aacgm_x-aacgm_x_x
rotMat[1,0,*] = aacgm_x-aacgm_x_y
rotMat[2,0,*] = aacgm_x-aacgm_x_z

rotMat[0,1,*] = aacgm_y-aacgm_y_x
rotMat[1,1,*] = aacgm_y-aacgm_y_y
rotMat[2,1,*] = aacgm_y-aacgm_y_z

rotMat[0,2,*] = aacgm_z-aacgm_z_x
rotMat[1,2,*] = aacgm_z-aacgm_z_y
rotMat[2,2,*] = aacgm_z-aacgm_z_z

rotMat = -rotMat
stop
; apply rotation

dbAACGM = rotMat ## [ $
[ gei_dbx[*] ],$
[ gei_dby[*] ],$
[ gei_dbz[*] ] ]

geoPack_bCarSp, aacgm_x, aacgm_y, aacgm_z, $
dbAACGM[*,0], dbAACGM[*,1], dbAACGM[*,2], $
aacgm_dbR, aacgm_dbTh, aacgm_dbPh
;geoPack_bCarSp, gei_x, gei_y, gei_z, $
; gei_dbx, gei_dby, gei_dbz, $
; aacgm_dbR, aacgm_dbTh, aacgm_dbPh

endif


end