Pro amp_trk_lab,nump,in_hist,bn_sz,Ln_a
; Code to find maxima in histogram function - used to check Lons and Track labels
; C. Waters
; October 2010
; Inputs:
;    nump : number of points (sum of histogram)
;    in_hist : histogram array
; Outputs:
;         Ln_a : array of longs that have 'most' data
;
  thresh=nump/30.
  ftn_a=[in_hist,in_hist[0]]         ; extend array to wrap 0->360 deg
  dv=ftn_a-shift(ftn_a,1)            ; histogram deriv

  n_pk=0                             ; initialise num peaks(lons) to zero
  pkLoc_a=intarr(n_elements(ftn_a))  ; ensure we ahve plenty of space for peak counting
  For ii=1,n_elements(ftn_a)-1 do begin
    if ( (dv(ii-1) ge 0.0) and (dv(ii) lt 0.0) ) then begin   ; look for +ve->ive transition to detect a max
      pkLoc_a(n_pk)=ii-1    ; store array index (prev)
      n_pk++                ; add peak number
    end
  end
  Ln_idx=where( ftn_a[pkLoc_a[0:n_pk-1]] gt thresh)  ; histogram array index of peaks (lons)
  Ln_a=pkLoc_a[Ln_idx]*bn_sz                         ; should be size 2 for normal operation

end





