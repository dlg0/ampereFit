; Widget wrapper for AMPERE data fitting code - AMP_FIT
; Sets parameters for FAC calculation from AMPERE data and then calls AMP_FIT
;
; C Waters
; January, 2010
;
; Modifications:
;
;
; -------------------------------------------------------------------------------------
;
; Widget event routines
Pro GetAmpPrms_event,ev
 Common WidgBlk,Base1,Base2,HrBut,MnSld,ScSld,DneBut
 Common WValBlk1,path,SavFileName,StHr,StMn,StSc,GetSc,TmSep,NumPlts,thresh,mxclat,sigma,kmax,mmax,$
  plt_png,HemSp,south
 Widget_Control,ev.id,Get_UValue=UVal
 Case UVal of
 'gfi' : Begin
          Ttle='Select Input AMPERE File'
          SavFileName=Dialog_PickFile(Title=Ttle,path=path,filter='*.sav')
          Widget_Control,HrBut,Set_Value='Start Hour = '+StrTrim(String(StHr),2)+' UT'
          Widget_Control,MnSld,Set_Value=StMn
          Widget_Control,ScSld,Set_Value=StSc
          Widget_Control,HrBut,Sensitive=1
          Widget_Control,MnSld,Sensitive=1
          Widget_Control,ScSld,Sensitive=1
          Widget_Control,DneBut,Sensitive=1
         end
 'shr' : Begin
          StHr=StHr+1
          If StHr gt 23 Then StHr=0
          Widget_Control,HrBut,Set_Value='Start Hour = '+StrTrim(String(StHr),2)+' UT'
         end
 'smn' : Widget_Control,ev.id,Get_Value=StMn
 'ssc' : Widget_Control,ev.id,Get_Value=StSc
 'asc' : Widget_Control,ev.id,Get_Value=GetSc
 'ins' : Widget_Control,ev.id,Get_Value=TmSep
 'npl' : Widget_Control,ev.id,Get_Value=NumPlts
 'sgs' : Widget_Control,ev.id,Get_Value=sigma
 'bth' : Widget_Control,ev.id,Get_Value=thresh
 'clt' : Widget_Control,ev.id,Get_Value=mxclat
 'lto' : Widget_Control,ev.id,Get_Value=kmax
 'lno' : Widget_Control,ev.id,Get_Value=mmax
 'pfc' : Widget_Control,ev.id,Get_Value=plt_png
 'hms' : Begin
          Widget_Control,ev.id,Get_Value=south
          If south eq 0 then HemSp='N'
          If south eq 1 then HemSp='S'
         end
 'dne' : Begin
          Widget_Control,ev.top,/DESTROY
          PrmFile=path+'amp_fit.par'
          OpenW,u1,PrmFile,/Get_Lun
          PrintF,u1,StHr,StMn,StSc,GetSc,TmSep,NumPlts
          PrintF,u1,south,kmax,mmax,thresh,mxclat,sigma,plt_png
          Free_Lun,u1
         end
 end
end
;
Pro GetAmpPrms
 Common WidgBlk,Base1,Base2,HrBut,MnSld,ScSld,DneBut
 Common WValBlk1,path,SavFileName,StHr,StMn,StSc,GetSc,TmSep,NumPlts,thresh,mxclat,sigma,kmax,mmax,$
  plt_png,HemSp,south

 WXPos=1 & WYPos=10 & WXSz=200  ; Widget Placement
 Base1 = WIDGET_BASE(Title='AMPERE FAC Menu',XOFFSET=WXPos,YOFFSET=WYPos,/COLUMN,XSize=WXSz)
 ADFBut = Widget_Button(base1,value='Select a New File',UValue='gfi')
 Widget_Control, ADFBut,Sensitive=1

 PrmFile=path+'amp_fit.par'
 OpenR,u1,PrmFile,/Get_Lun
 StHr=8 & StMn=0 & StSc=0 & GetSc=3600 & TmSep=900 & NumPlts=1
 ReadF,u1,StHr,StMn,StSc,GetSc,TmSep,NumPlts
 south=0 & kmax=35 & mmax=5 & thresh=80 & mxclat=50. & sigma=2. & plt_png=1
 ReadF,u1,south,kmax,mmax,thresh,mxclat,sigma,plt_png
 Free_Lun,u1

 HrBut = Widget_Button(base1,value='Start Hour = '+StrTrim(String(StHr),2)+' UT',UValue='shr')
 Widget_Control,HrBut,Sensitive=0
 MnSld = Widget_Slider(base1,value=StMn,UValue='smn',Minimum=0,Maximum=60,Title='Start Minute')
 Widget_Control,MnSld,Sensitive=0
 ScSld = Widget_SLider(base1,value=StSc,UValue='ssc',Minimum=0,Maximum=60,Title='Start Second')
 Widget_Control,ScSld,Sensitive=0
 AvSld = Widget_SLider(base1,value=GetSc,UValue='asc',Minimum=15,Maximum=7200,Title='Averaging Time (seconds)')
 IncSld= Widget_SLider(base1,value=TmSep,UValue='ins',Minimum=10,Maximum=3600,Title='Seconds to Shift')
 NmPSld= Widget_SLider(base1,value=NumPlts,UValue='npl',Minimum=1,Maximum=100,Title='Number of Plots to do')
 DneBut = Widget_Button(base1,value='GO',UValue='dne')
; Widget_Control,DneBut,Sensitive=0
 Widget_Control,DneBut,Sensitive=1
 WIDGET_CONTROL, base1, /REALIZE

 Base2 = Widget_Base(Title='AMPERE FAC Menu',Group_Leader=Base1,XOffset=WXPos+WXSz+8,YOffset=WYPos,/Column,XSize=WXSz)
 SgSld = Widget_SLider(base2,value=sigma,UValue='sgs',Minimum=1,Maximum=10,Title='Data Sigma')
 bthSld= Widget_SLider(base2,value=thresh,UValue='bth',Minimum=1,Maximum=200,Title='dB Threshold [nT]')
 pltSld= Widget_SLider(base2,value=mxclat,UValue='clt',Minimum=20,Maximum=90,Title='Plot CoLat Range [Deg]')
 LatSld= Widget_SLider(base2,value=kmax,UValue='lto',Minimum=2,Maximum=65,Title='Latitude Order')
 LonSld= Widget_SLider(base2,value=mmax,UValue='lno',Minimum=2,Maximum=6,Title='Longitude Order')
 z1Arr=StrArr(2)
 For i=0,1 do z1Arr(i)=StrTrim(String(i),2)
 PFac = CW_BGroup(base2,z1Arr,UValue='pfc',Label_top='PNG FAC Output',Set_Value=plt_png,/row,/Exclusive)
 hArr=StrArr(2)
 hArr(0)='N'
 hArr(1)='S'
 HmSB = CW_BGroup(base2,hArr,UValue='hms',Label_top='Hemisphere',Set_Value=south,/row,/Exclusive)
 HemSp=hArr(south)
 WIDGET_CONTROL, base2, /REALIZE
 XMANAGER, 'GetAmpPrms',base1,/NO_BLOCK
 XMANAGER, 'GetAmpPrms',base2
end
;
; ------------------------------------------------------------------------------------------------------------
;
; Main driver code starts here
Pro amp_fit_wij
 Common WValBlk1,path,SavFileName,StHr,StMn,StSc,GetSc,TmSep,NumPlts,thresh,mxclat,sigma,kmax,mmax,$
  plt_png,HemSp,south

; set default values
  if strCmp (!version.os, 'Win32') or strCmp (!version.os, 'WIN62') then begin
	plotDev	= 'win'
	path	= 'd:\cwac\hi_res\davidg\'
  endif else begin                        ; other systems (linux, darwin etc.)
	plotDev = 'X'
	path = '~/code/ampereFit/idl/'
;	pnmPath = path + 'pnmSavs/pnmSav'
  endelse
  set_plot,plotdev
  loadct,0
  device,decomposed=0
  aacgmpath=path                          ; location of AACGM coeff file (*.asc)
  mn_fac=0.05              ; do not plot abs(FAC) lt than this
  mx_fac=1.0

  GetAmpPrms
  nLatGrid	= mxclat    ; fit grid defaults
  plt_coLat=mxclat
  nLonGrid	= 24
  plot_bFns = 0            ; plot switch for basis set diagnostics
  sHr = float(StHr)+float(StMn)/60.0+StSc/3600.0
  eHr = sHr+GetSc/3600.0
  plt_tracks = 0           ; plot dbTh and dbPh data by Iridium orbit track

  amp_fit, sHr, eHr, south, $
        plot_bFns = plot_bFns, $
		path = path, $
		aacgmpath = aacgmpath, $
;		pnmPath = pnmpath, $
		savFileName = SavFileName, $
		kmax = kmax, mmax = mmax, $
		thresh = thresh, $
		sigma = sigma, $
		nLatGrid = nLatGrid, nLonGrid = nLonGrid, $
		mn_fac = mn_fac, mx_fac=mx_fac, $
		plt_png = plt_png, $
		plt_tracks = plt_tracks, $
		plt_coLat = plt_coLat
  print,'Finished'
 end
