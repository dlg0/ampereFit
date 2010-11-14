pro write_ird_file,out_fname, yrstr,daystr,t0str,d_epoch, $
   coLat_deg_arr, Lon_arr, $
   dbThet_th_arr, dbPh_th_arr, dbThet_ph_arr, dbPh_ph_arr, $
   orbPln_arr, iSat_arr, qual_arr, splice_arr
; Writes a IRD file to disk
; C.L. waters
; JHU/APL
; Sept 2010
;
  openw,wunit,out_fname,/get_lun
  printf,wunit,yrstr, ' ',strmid(daystr,0,2),' ',strmid(daystr,2,2)
  printf,wunit,strmid(t0str,0,2),' ',strmid(t0str,2,2),' ',strmid(t0str,4,2)
  printf,wunit,strtrim(string(long(d_epoch)),2)
  
  for i=0,n_elements(coLat_deg_arr)-1 do begin  
    printf,wunit,coLat_deg_arr[i], Lon_arr[i], $
      dBThet_th_arr[i],dBPh_th_arr[i], $
      dBThet_ph_arr[i],dBPh_ph_arr[i], $
      orbPln_arr[i], iSat_arr[i], qual_arr[i], splice_arr[i], $
      format='(f6.2,1x,f8.5,1x,4(f8.2,1x),2(i3,1x),f5.2,1x,i2)'
  endfor   
  close,wunit
  free_lun,wunit
end