pro write_wsp_file, out_fname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, nLatGrid, nLonGrid, $
    GEI_coLat_arr, GEI_Lon_arr, GEO_Lon_arr, $
    dbThet_arr, dbPh_arr, $
    jPar_arr

; Writes a WSP (whole sphere) file to disk
; C.L. Waters
; Univ of Newcastle
; Australia
; Oct 2010
;
  openw,wunit,out_fname,/get_lun
  printf,wunit,yrstr, ' ',strmid(daystr,0,2),' ',strmid(daystr,2,2)
  printf,wunit,strmid(t0str,0,2),' ',strmid(t0str,2,2),' ',strmid(t0str,4,2)
  printf,wunit,strtrim(string(long(d_epoch)),2)
  printf,wunit,strtrim(string(long(KMax)),2),' ',strtrim(string(long(Mmax)),2),' ',resol_deg
  printf,wunit,strtrim(string(long(nLatGrid)),2),' ',strtrim(string(long(nLonGrid)),2)
  
  for i=0,n_elements(GEI_coLat_arr[*])-1 do begin    ; theta is +ve Nth, phi is +ve east
    printf,wunit,GEI_coLat_arr[i], GEI_Lon_arr[i], GEO_Lon_arr[i], $
                 dBThet_arr[i], dBPh_arr[i], $
                 jPar_arr[i], $
           format='(f6.2,1x,2(f8.5,1x),3(f8.2,1x))'
  endfor            
  close,wunit
  free_lun,wunit
end
