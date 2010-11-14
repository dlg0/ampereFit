pro write_cof_file,out_fname, yrstr, daystr, t0str, d_epoch, $
    Kmax, Mmax, resol_deg, kArr,mArr, coeff_a
;
; Writes the solution coefficients to a COF file
; C.L. waters
; JHU/APL
; Sept 2010
;
  nB=n_elements(kArr)
  openw,wunit,out_fname,/get_lun
  printf,wunit,yrstr, ' ',strmid(daystr,0,2),' ',strmid(daystr,2,2)
  printf,wunit,strmid(t0str,0,2),' ',strmid(t0str,2,2),' ',strmid(t0str,4,2)
  printf,wunit,strtrim(string(long(d_epoch)),2)
  printf,wunit,strtrim(string(long(KMax)),2),' ',strtrim(string(long(Mmax)),2),' ',resol_deg
  printf,wunit,strtrim(string(long(nB)),2)
  
  for i=0,nB-1 do begin
    printf,wunit,kArr[i],mArr[i],coeff_a[i], coeff_a[nB+i], $   ; k,m,thet,phi
      format='(f8.4,1x,i4,1x,2(f8.4,1x))'      ; k is float - spherical cap code has k - non-integer
  endfor            
  close,wunit
  free_lun,wunit
end
