; Widget wrapper for AMPERE data fitting code - AMP_FIT - full sphere
; Sets parameters for FAC calculation from AMPERE data and then calls AMP_FIT_F
;
; C Waters
; Space Physics
; University of Newcastle
; New South Wales, Australia
; Sept, 2010
;
; Modifications:
;
; -------------------------------------------------------------------------------------
;
; Widget event routines

Pro GetAmpPrms_full_event,ev
 Common WidgBlk_full,Base1,Base2,HrSld,MnSld,ScSld,DneBut
 Common WFileBlk_full,path,datFileName,amp_data_path
 Common WParBlk_full,StHr,StMn,StSc,GetSc,kmax,SpacRes,SpTxt,Mmax,plt_png,debug

 Widget_Control,ev.id,Get_UValue=UVal
 Case UVal of
 'gfi' : Begin
          Ttle='Select Input AMPERE File'
;          datFileName=Dialog_PickFile(Title=Ttle,path=amp_data_path,filter='*.sav')
          datFileName=Dialog_PickFile(Title=Ttle,path=amp_data_path,filter='*.ncdf')
          Widget_Control,HrSld,Set_Value=StHr
          Widget_Control,MnSld,Set_Value=StMn
          Widget_Control,ScSld,Set_Value=StSc
          Widget_Control,HrSld,Sensitive=1
          Widget_Control,MnSld,Sensitive=1
          Widget_Control,ScSld,Sensitive=1
          Widget_Control,DneBut,Sensitive=1
         end
 'shr' : Widget_Control,ev.id,Get_Value=StHr
 'smn' : Widget_Control,ev.id,Get_Value=StMn
 'ssc' : Widget_Control,ev.id,Get_Value=StSc
 'asc' : Widget_Control,ev.id,Get_Value=GetSc
 'lto' : begin
           Widget_Control,ev.id,Get_Value=kmax
; not necessary - just do a 180/order - CW ( here from cap code)
;           coLat = (fIndGen(180)+1.)*!dpi/180.d0
;           legendreFns = ampere_Legendre_full ( coLat, kmax, 1, $
;              nVals = nVals, $
;;              dLegendreFns = dLegendreFns, $
;              kVals = kVals)
;; Use an FFT to get the frequency/wavelength of the function, not necessary here but useful for cap fit code
;            amp   = abs(fft(legendreFns[*,kmax-1]))               ; Find the major frequency in power spectrum
;            mx    = max(amp,idx)
;            fr    = float(idx)
;            spacRes = 90.0/fr      ; Resolution
            SpacRes=180.0/kmax
            Widget_Control,SpTxt,Set_value='1/2 Wavelength res = '+strtrim(string(SpacRes,format='(f4.1)'),2)+' deg'
          end
 ; -------------------------------------------------------------------
 'lno' : Widget_Control,ev.id,Get_Value=mmax
 'pfc' : Widget_Control,ev.id,Get_Value=plt_png
 'dbg' : Widget_Control,ev.id,Get_Value=debug
 'dne' : Begin
          Widget_Control,ev.top,/DESTROY
          PrmFile=path+'amp_fit_full.par'
          OpenW,u1,PrmFile,/Get_Lun
          PrintF,u1,StHr,StMn,StSc,GetSc
          PrintF,u1,kmax,mmax,plt_png,debug
          Free_Lun,u1
         end
 end
end
;
Pro GetAmpPrms_full
 Common WidgBlk_full,Base1,Base2,HrSld,MnSld,ScSld,DneBut
 Common WFileBlk_full,path,datFileName,amp_data_path
 Common WParBlk_full,StHr,StMn,StSc,GetSc,kmax,SpacRes,SpTxt,Mmax,plt_png,debug

 WXPos=1 & WYPos=10 & WXSz=200  ; Widget Placement
 Base1 = WIDGET_BASE(Title='AMPERE FAC Menu',XOFFSET=WXPos,YOFFSET=WYPos,/COLUMN,XSize=WXSz)
 ADFBut = Widget_Button(base1,value='Select a New File',UValue='gfi')
 Widget_Control, ADFBut,Sensitive=1

 PrmFile=path+'amp_fit_full.par'
 OpenR,u1,PrmFile,/Get_Lun
 StHr=8 & StMn=0 & StSc=0 & GetSc=3600
 ReadF,u1,StHr,StMn,StSc,GetSc
 south=0 & kmax=35 & mmax=5 & mxclat=50. & sigma=2. & plt_png=1
 ReadF,u1,kmax,mmax,plt_png,debug
 Free_Lun,u1

 HrSld = Widget_SLider(base1,value=StHr,UValue='shr',Minimum=0,Maximum=23,Title='Start Hour')
 Widget_Control,HrSld,Sensitive=0
 MnSld = Widget_Slider(base1,value=StMn,UValue='smn',Minimum=0,Maximum=60,Title='Start Minute')
 Widget_Control,MnSld,Sensitive=0
 ScSld = Widget_SLider(base1,value=StSc,UValue='ssc',Minimum=0,Maximum=60,Title='Start Second')
 Widget_Control,ScSld,Sensitive=0
 AvSld = Widget_SLider(base1,value=GetSc,UValue='asc',Minimum=15,Maximum=7200,Title='Averaging Time (seconds)')
 DneBut = Widget_Button(base1,value='DO IT!',UValue='dne')
; Widget_Control,DneBut,Sensitive=0
 Widget_Control,DneBut,Sensitive=1
 WIDGET_CONTROL, base1, /REALIZE

 Base2 = Widget_Base(Title='AMPERE FAC Menu',Group_Leader=Base1,XOffset=WXPos+WXSz+8,YOffset=WYPos,/Column,XSize=WXSz)
 LatSld= Widget_SLider(base2,value=kmax,UValue='lto',Minimum=2,Maximum=65,Title='Latitude Order')
 SpacRes=180.0/kmax
 SpTxt = Widget_Text(base2,value='1/2 Wavelength res = '+strtrim(string(SpacRes,format='(f3.1)'),2)+' deg')
 LonSld= Widget_SLider(base2,value=mmax,UValue='lno',Minimum=2,Maximum=6,Title='Longitude Order')
;
 z1Arr=StrArr(2)
 For i=0,1 do z1Arr(i)=StrTrim(String(i),2)
 PFac = CW_BGroup(base2,z1Arr,UValue='pfc',Label_top='PNG FAC Output',Set_Value=plt_png,/row,/Exclusive)
;
 dbgch=StrArr(2)
 For i=0,1 do dbgch(i)=StrTrim(String(i),2)
 dbSB = CW_BGroup(base2,dbgch,UValue='dbg',Label_top='Debug',Set_Value=debug,/row,/Exclusive)
;
 WIDGET_CONTROL, base2, /REALIZE
 XMANAGER, 'GetAmpPrms_full',base1,/NO_BLOCK
 XMANAGER, 'GetAmpPrms_full',base2
end
;
; ------------------------------------------------------------------------------------------------------------
;
; Main driver code starts here
Pro amp_fit_wij_full
 Common WFileBlk_full,path,datFileName,amp_data_path
 Common WParBlk_full,StHr,StMn,StSc,GetSc,kmax,SpacRes,SpTxt,Mmax,plt_png,debug

; set default values
; set path, aacgmpath, amp_data_path, amp_fits_path (for output)
  @amp_fit_paths_full

  if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN64') then begin
	plotDev	= 'win'
;	path	= 'd:\cwac\hi_res\davidg\'
  endif else begin                        ; other systems (linux, darwin etc.)
	plotDev = 'X'
;	path = '~/code/ampereFit/idl/'
  endelse
  set_plot,plotdev
  loadct,0
  device,decomposed=0
  mn_fac=0.05              ; do not plot abs(FAC) lt than this
  mx_fac=1.0

  GetAmpPrms_full

  nLonGrid	= 24

  datestr=strmid(file_basename(datFileName),0,8)
  syr=fix(strmid(datestr,0,4))
  smo=fix(strmid(datestr,4,2))
  sday=fix(strmid(datestr,6,2))
  geopack_epoch,sepoch,syr,smo,sday,StHr,StMn,StSc,/compute
  eepoch=sepoch+GetSc*1.d3
  geopack_epoch,eepoch,eyr,emo,eday,enHr,enMn,enSc,/breakdown

  sHr = float(StHr)+float(StMn)/60.0+StSc/3600.0
  eHr = (eday-sday)*24.0d0+float(enHr)+float(enMn)/60.0+enSc/3600.0
  if ((eday-sday) eq 0) then extend=0 else extend=1

  amp_fit_full, sHr, eHr, $
		datFileName = datFileName, $
		kmax = kmax, mmax = mmax, $
		nLonGrid = nLonGrid, $
		mn_fac = mn_fac, mx_fac=mx_fac, $
		plt_png = plt_png, $
		debug=debug, $
		sepoch = sepoch, eepoch = eepoch, $
        extend = extend, $
        SpacRes = SpacRes

  print,'Finished'
 end
