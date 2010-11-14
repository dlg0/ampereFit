pro ampere_plot

; @'d:\cwac\hi_res\davidg\jpar_ver2\schabasisfunctions.pro'

 if strCmp (!version.os, 'linux' or !version.os, 'darwin' ) then begin    ; Linux machine for Dave
  plotDev = 'X'
  path = '~/code/ampereFit/idl/'
  pnmPath = path + 'pnmSavs/pnmSav'
 endif else begin                                ; Windows machine for Colin
  plotDev = 'win'
  path = 'd:\cwac\hi_res\davidg\'
 endelse

; InF='20091022Amp_invert.sav'
 InF='~/code/ampereFit/data/20091023Amp_invert_test2.sav'
 plotCapSize = 40.0
 fileName = InF

 sHr = 10.1
 eHr = 10.3

 restore,filename     ; B_ECI[3,np], PLANE_NUMBER_TOTAL[np], POS_ECI_TOTAL[3,np], X_AXIS_FRAC_HOUR[np], PSEUDO_SVNUM_TOTAL[np]

 ndat = n_elements(plane_number_total)
 data_struct = { ipln:0,$
 sv:0,$
 stme:0d0,$
 dbx:0d0,$
 dby:0d0,$
 dbz:0d0,$
 px:0d0,$
 py:0d0,$
 pz:0d0}
 data = replicate ( data_struct, ndat )

 data.stme=x_axis_frac_hour
 data.dbx=reform(B_ECI(0,*))
 data.dby=reform(B_ECI(1,*))
 data.dbz=reform(B_ECI(2,*))
 data.px=reform(POS_ECI_TOTAL(0,*))/1000.           ; conv to km
 data.py=reform(POS_ECI_TOTAL(1,*))/1000.
 data.pz=reform(POS_ECI_TOTAL(2,*))/1000.

; Select out a time slice
 iiTime = where (data.stme gt sHr $
          and data.stme lt EHr $
          and data.pz gt 0, iiCnt)          ; set for Nth Hemis with data.pz
 data = data[iiTime]

 geopack_sphcar,data.px,data.py,data.pz,rg_th,cglat_th,glon_th,/to_sphere,/degree   ; find GEO(r,thet,phi)
; Put the x,y,z GEO input db data into spherical components
 geoPack_bcarsp, data.px, data.py, data.pz, data.dbx, data.dby, data.dbz, bR, bTheta, bPhi

; plot stuff
; ----------

 set_plot, plotDev
 device,decomposed = 0
 !p.background = 255
 winNum = 1

; plot the db vectors
; -------------------

 window, winnum, xSize = 600, ySize = 600
 winNum++

 plt_dat, cglat_th, glon_th, $
  n_elements ( cglat_th), -bTheta, bPhi, [1,1], [1,2], title = 'Raw Data: SHr='+strtrim(string(SHr),2)

 stop
end
