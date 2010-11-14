pro ampere_shift_coord, pi_x, pi_y, pi_z, $   ; input GEI coords
                        po_x, po_y, po_z, $   ; shifted GEI coords
                        r_gei, th_gei_deg, ph_gei_deg, $  ; spherical coords of shifted coords
                        r_mat          ; rotation matrix
; Given coords (or vectors), process through the rotation matrix
; 
; C.L. Waters and D.L. Green
; JHU/APL
; Sept 2010
;
  po_x = pi_x
  po_y = pi_y
  po_z = pi_z

; Rotate input GEI (x,y,z) location data to shifted coords
  For ii=Long(0),n_elements(pi_x)-1 do begin

; rotate coord locations
    xyz_dat = [pi_x[ii], pi_y[ii], pi_z[ii]]
    xyz_shift = matrix_multiply ( r_mat, xyz_dat, /Atrans )

    po_x[ii] = xyz_shift[0]
    po_y[ii] = xyz_shift[1]
    po_z[ii] = xyz_shift[2]

  end
; Get shifted spherical coord locations
  geopack_sphcar, po_x, po_y, po_z, $
                  r_gei, th_gei_deg, ph_gei_deg, /to_sphere, /degree   ; find shifted GEI(r,thet,phi)

end
