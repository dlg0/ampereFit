pro ampere_shift_data, data_in, data_out, r_mat
; Given coords and vectors, process through the rotation matrix
; input : data_in structure - modifies the following tags
;  px, py, pz
;  dbx, dby, dbz
;  GEI_R_km, GEI_coLat_deg, GEI_Lon_Deg
;  GEI_coLat_rad, GEI_lon_rad
;  br_GEI, bTheta_GEI, bPhi_GEI
;  
; C.L. Waters
; JHU/APL
; Sept 2010
;
  data_out = data_in        ; define output data structure

; Rotate input GEI (x,y,z) location data to shifted coords
  For ii=Long(0),n_elements(data_in)-1 do begin

; rotate coord locations
    xyz_dat = [data_in[ii].px, data_in[ii].py, data_in[ii].pz]
    xyz_shift = matrix_multiply ( r_mat, xyz_dat, /Atrans )

    data_out[ii].px = xyz_shift[0]
    data_out[ii].py = xyz_shift[1]
    data_out[ii].pz = xyz_shift[2]

  end

; Get shifted spherical coord locations
  geopack_sphcar, data_out.px, data_out.py, data_out.pz, $
                  rg_th, cglat, glon, /to_sphere, /degree   ; find shifted GEI(r,thet,phi)
  data_out.GEI_R_km=rg_th
  data_out.GEI_coLat_deg=cglat
  data_out.GEI_lon_deg=glon

  data_out.GEI_coLat_rad = data_out.GEI_coLat_deg * !dtor
  data_out.GEI_lon_rad = data_out.GEI_lon_deg * !dtor

; Get spherical dB components, calc at shifted (x,y,z) coords
; rotate vectors
  For ii=Long(0),n_elements(data_in)-1 do begin
    vec_dat = [data_in[ii].dbx, data_in[ii].dby, data_in[ii].dbz]
    vec_shift = matrix_multiply ( r_mat, vec_dat, /Atrans )

    data_out[ii].dbx = vec_shift[0]
    data_out[ii].dby = vec_shift[1]
    data_out[ii].dbz = vec_shift[2]
    geopack_bcarsp, data_out[ii].px,  data_out[ii].py,  data_out[ii].pz, $
                    data_out[ii].dbx, data_out[ii].dby, data_out[ii].dbz,$
                    vbr, vbth, vbph
    data_out[ii].br_GEI=vbr
    data_out[ii].bTheta_GEI=vbth
    data_out[ii].bPhi_GEI=vbph
  end
end
