pro write_grd_file,out_fname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, nLatGrid, nLonGrid, $
    coLat_arr, Lon_arr, $
    dbThet_th_arr, dbPh_th_arr, $
    dbThet_ph_arr, dbPh_ph_arr, $
    jPar_arr
; Writes a GRD file to disk
; C.L. waters
; JHU/APL
; Sept 2010
;
  openw,wunit,out_fname,/get_lun
  printf,wunit,yrstr, ' ',strmid(daystr,0,2),' ',strmid(daystr,2,2)
  printf,wunit,strmid(t0str,0,2),' ',strmid(t0str,2,2),' ',strmid(t0str,4,2)
  printf,wunit,strtrim(string(long(d_epoch)),2)
  printf,wunit,strtrim(string(long(KMax)),2),' ',strtrim(string(long(Mmax)),2),' ',resol_deg
  printf,wunit,strtrim(string(long(nLatGrid)),2),' ',strtrim(string(long(nLonGrid)),2)
  
  for i=0,n_elements(coLat_arr[*])-1 do begin    ; theta is +ve Nth, phi is +ve east
    printf,wunit,coLat_arr[i], Lon_arr[i], $
      dBThet_th_arr[i], dBPh_th_arr[i], $
      dBThet_ph_arr[i], dBPh_ph_arr[i], $
      jPar_arr[i], $
      format='(f6.2,1x,f8.5,1x,5(f8.2,1x))'
  endfor            
  close,wunit
  free_lun,wunit
end
