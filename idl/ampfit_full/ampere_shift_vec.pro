pro ampere_shift_vec, pi_x, pi_y, pi_z, $     ; input GEI locations
                      vi_x, vi_y, vi_z, $     ; input GEI vectors
                      vs_x, vs_y, vs_z, $     ; shifted GEI vectors
                      db_r, db_th, db_ph, $   ; spherical components of shifted vecs
                      r_mat                   ; rotation matrix
;
; C.L. Waters and D.L. Green
; JHU/APL
; Sept 2010
;
  vs_x = vi_x
  vs_y = vi_y
  vs_z = vi_z

  db_r  = vi_x
  db_th = vi_x
  db_ph = vi_x

; Rotate input vecs
  For ii=Long(0),n_elements(vi_x)-1 do begin
    vec_dat = [vi_x[ii], vi_y[ii], vi_z[ii]]
    vec_shift = matrix_multiply ( r_mat, vec_dat, /Atrans )

    vs_x[ii] = vec_shift[0]
    vs_y[ii] = vec_shift[1]
    vs_z[ii] = vec_shift[2]

; Get spherical dB components, calc at shifted (x,y,z) coords
    geopack_bcarsp, pi_x[ii],  pi_y[ii],  pi_z[ii], $
                  vs_x[ii], vs_y[ii], vs_z[ii],$
                  b_r, b_th, b_ph
    db_r[ii] = b_r
    db_th[ii] = b_th
    db_ph[ii] = b_ph
  end

end
