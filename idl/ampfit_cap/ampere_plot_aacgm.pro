; Plot utitliy for AMPERE Data
; Presents data in GEI (straight from IDL SAV file) and also in AACGM
;
; C.L. Waters
; Centre for Space Physics
; University of Newcastle, Australia
; Dec 2009
;
; Notes:
;  Data structure is
;          data_rw_struct = {ipln:0, sv:0, utc:0d0,$
;          px:0d0, py:0d0, pz:0d0,$
;          dbx:0d0, dby:0d0, dbz:0d0,$
;          pr:0d0, pth:0d0, pph:0d0,$
;          dbr:0d0, dbth:0d0, dbph:0d0}
;
Pro tmestr,Hr,Mn,Sc,t_str
; Form string expression from time
 hr_str=strtrim(string(Hr),2)
 If Hr lt 10 then hr_str='0'+hr_str
 mn_str=strtrim(string(Mn),2)
 If Mn lt 10 then mn_str='0'+mn_str
 sc_str=strtrim(string(fix(sc)),2)
 If Sc lt 10 then sc_str='0'+sc_str
 t_str=hr_str+':'+mn_str+':'+sc_str
end
;
; -----------------------------------------------------------------------------
;
Pro AmpPar_mag_event,ev
 Common WidgBlk,base_in,FNme,DFNm,NxtBut,PrvBut,HrSld,MnSld,ScSld,DneBut
 Common WValBlk,InF,HemSp,StHr,StMn,StSc,GetSc,TmSep,path
 Common DataBlk,rw_data, year, month, day
 Widget_Control,ev.id,Get_UValue=UVal
 Case UVal of
 'gfi' : Begin
          Ttle='Select Input Data File'
          InF=Dialog_PickFile(Title=Ttle,path=path,filter='*.sav')
          restore,InF     ; B_ECI[3,np], PLANE_NUMBER_TOTAL[np], POS_ECI_TOTAL[3,np], X_AXIS_FRAC_HOUR[np], PSEUDO_SVNUM_TOTAL[np]
          reads, strmid(file_baseName (InF),0,8), year, month, day, format='(i4,i2,i2)'
          ndat = n_elements(plane_number_total)
          data_rw_struct = {ipln:0, sv:0, utc:0d0,$
          px:0d0, py:0d0, pz:0d0,$
          dbx:0d0, dby:0d0, dbz:0d0,$
          pr:0d0, pth:0d0, pph:0d0,$
          dbr:0d0, dbth:0d0, dbph:0d0}
          rw_data = replicate(data_rw_struct,ndat)
          rw_data.ipln=PLANE_NUMBER_TOTAL
          rw_data.utc=x_axis_frac_hour
          rw_data.px=reform(POS_ECI_TOTAL(0,*))/1000.           ; conv to km. These are GEI, X,Y,Z coords
          rw_data.py=reform(POS_ECI_TOTAL(1,*))/1000.
          rw_data.pz=reform(POS_ECI_TOTAL(2,*))/1000.
          rw_data.dbx=reform(B_ECI(0,*))
          rw_data.dby=reform(B_ECI(1,*))
          rw_data.dbz=reform(B_ECI(2,*))
          st_tme=float(StHr)+float(StMn)/60.0+StSc/3600.
          en_tme=st_tme+float(StHr)+float(TmSep)/3600.
          tmestr,StHr,StMn,StSc,st_str
; Select out data for this time interval
          if HemSp eq 'N' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz gt 0., iiCnt) ; set for Nth Hemis with data.pz
          if HemSp eq 'S' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz lt 0., iiCnt) ; set for Sth Hemis with data.pz
          data_in = rw_data[iiTime]
          np=n_elements(data_in)
; Get spherical coords of XYZ locations
          sphcar,rg_th,cglat,glon,data_in.px,data_in.py,data_in.pz,/to_sphere,/degree   ; find GEO(r,thet,phi)
          data_in[*].pr=rg_th
          data_in[*].pth=cglat
          data_in[*].pph=glon
; Put the x,y,z GEI input db data into spherical components
          For ii=Long(0),Long(np)-Long(1) do begin
           bcarsp, data_in[ii].px, data_in[ii].py, data_in[ii].pz, data_in[ii].dbx, data_in[ii].dby, data_in[ii].dbz, vbr,vbth,vbph
            data_in[ii].dbr=vbr
            data_in[ii].dbth=vbth
            data_in[ii].dbph=vbph
          end
          If HemSp eq 'S' then begin
           idx=where(cglat gt 90.)             ; Find Sth hemisphere points
           If (idx(0) gt -1) then cglat[idx]=180.0-cglat[idx] ; deal with Sth hemis
          end
; plot the db vectors
          window, 1, xSize = 500, ySize = 500,title='Input Data'
          If HemSp eq 'S' Then plt_dat,cglat,glon,np, data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
          If HemSp eq 'N' Then plt_dat,cglat,glon,np,-data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
; Take GEI -> GEOG -> AACGM
          Re_km = 6371.0
          R_km = 110.0 + Re_km
; conv GEI r,thet,phi -> GEI XYZ
          sphcar,data_in.pr, data_in.pth, data_in.pph,$
           gei_x, gei_y, gei_z, /degree, /to_rect
; set up GEOPACK
          geoPack_reCalc, year, month, day, /date
          avgHour	= fix ( (st_tme + en_tme) / 2.0 )
          avgMin	= ( (st_tme + en_tme) / 2.0 mod 1 ) * 60
          yrSecAvg	= cnvTime(year, month, day, avgHour, avgMin, 0.0)
; conv GEI XYZ -> GEOG XYZ
          geoPack_conv_coord, gei_x, gei_y, gei_z, $
                              geo_x, geo_y, geo_z, $
                              /from_gei, /to_geo

          cdf_epoch, epoch0, year, month, day, /compute_epoch
          epoch=data_in.utc*3600d0*1d3 + epoch0
          avgEpoch = mean(epoch)
          epoch = fltArr(n_elements(gei_x[*])) + AvgEpoch
          mltShift = aacgm_mlt(year, yrSecAvg, 0.0 )
; conv GEOG XYZ -> GEOG r,thet,phi
          sphcar,geog_R_km, geog_coLat_deg, geog_lon_deg, $
                 geo_x, geo_y, geo_z, /degree ,/to_sphere

          gLat_a=90.0-geog_coLat_deg         ; convert to Latitude
          aacgm_conv_vec, gLat_a, geog_lon_deg, geog_R_km-Re_km, $
                    data_in.dbth, data_in.dbph, $
     mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph,err,/to_aacgm

; Plot AACGM data
          window, 3, xSize = 500, ySize = 500,title='AACGM-MLT'
          If HemSp eq 'S' Then 	$
           plt_dat, 90.0-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	                np, mth_vec_gth, mph_vec_gph, [1,1], [1,2], title =st_str
          If HemSp eq 'N' Then $
           plt_dat,90.-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
                    np,-mth_vec_gth, mph_vec_gph, [1,1], [1,2], title=st_str

; Turn on widget controls
          Widget_Control,ScSld,Set_Value=StSc
          Widget_Control,HrSld,Sensitive=1
          Widget_Control,MnSld,Sensitive=1
          Widget_Control,ScSld,Sensitive=1
          Widget_Control,NxtBut,Sensitive=1
          Widget_Control,PrvBut,Sensitive=1
         end
 'hms' : Begin
          Widget_Control,ev.id,Get_Value=Hm
          If Hm eq 0 then HemSp='N'
          If Hm eq 1 then HemSp='S'
         end
 'shr' : Widget_Control,ev.id,Get_Value=StHr
 'smn' : Widget_Control,ev.id,Get_Value=StMn
 'ssc' : Widget_Control,ev.id,Get_Value=StSc
 'asc' : Widget_Control,ev.id,Get_Value=GetSc
 'ins' : Widget_Control,ev.id,Get_Value=TmSep
 'nxt' : Begin
          st_tme_sec=float(StHr)*3600.0+float(StMn)*60.0+StSc
          st_tme_sec=st_tme_sec+TmSep                             ; Advance the time
          en_tme_sec=st_tme_sec+GetSc
          If en_tme_sec gt 86400.0 Then begin
           st_tme_sec=86400.0-GetSc
           en_tme_sec=86400.0
          end
          st_tme=st_tme_sec/3600.0
          en_tme=en_tme_sec/3600.0
          if HemSp eq 'N' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz gt 0., iiCnt) ; set for Nth Hemis with data.pz
          if HemSp eq 'S' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz lt 0., iiCnt) ; set for Sth Hemis with data.pz
          data_in = rw_data[iiTime]
          np=n_elements(data_in)
          StHr=fix(st_tme_sec/3600.0)
          StMn=fix((st_tme_sec-float(StHr)*3600.0)/60.0)
          StSc=st_tme_sec-float(StHr)*3600.0-float(StMn)*60.0
          Widget_Control,HrSld,Set_Value=StHr
          Widget_Control,MnSld,Set_Value=StMn
          Widget_Control,ScSld,Set_Value=StSc
          tmestr,StHr,StMn,StSc,st_str
; Get spherical coords of XYZ locations
          sphcar,rg_th,cglat,glon,data_in.px,data_in.py,data_in.pz,/to_sphere,/degree   ; find GEO(r,thet,phi)
          data_in[*].pr=rg_th
          data_in[*].pth=cglat
          data_in[*].pph=glon
; Put the x,y,z GEO input db data into spherical components
          For ii=Long(0),Long(np)-Long(1) do begin
           bcarsp, data_in[ii].px, data_in[ii].py, data_in[ii].pz, data_in[ii].dbx, data_in[ii].dby, data_in[ii].dbz, vbr,vbth,vbph
            data_in[ii].dbr=vbr
            data_in[ii].dbth=vbth
            data_in[ii].dbph=vbph
          end
          If HemSp eq 'S' then begin
           idx=where(cglat gt 90.)             ; Find Sth hemisphere points
           If (idx(0) gt -1) then cglat[idx]=180.0-cglat[idx] ; deal with Sth hemis
          end
; plot the db vectors
          window, 1, xSize = 500, ySize = 500,title='Input Data'
          If HemSp eq 'S' Then plt_dat,cglat,glon,np, data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
          If HemSp eq 'N' Then plt_dat,cglat,glon,np,-data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
; Take GEI -> GEOG -> AACGM
          Re_km = 6371.0
          R_km = 110.0 + Re_km
; conv GEI r,thet,phi -> GEI XYZ
          sphcar,data_in.pr, data_in.pth, data_in.pph,$
           gei_x, gei_y, gei_z, /degree, /to_rect
; set up GEOPACK
          geoPack_reCalc, year, month, day, /date
          avgHour	= fix ( (st_tme + en_tme) / 2.0 )
          avgMin	= ( (st_tme + en_tme) / 2.0 mod 1 ) * 60
          yrSecAvg	= cnvTime(year, month, day, avgHour, avgMin, 0.0)
; conv GEI XYZ -> GEOG XYZ
          geoPack_conv_coord, gei_x, gei_y, gei_z, $
                              geo_x, geo_y, geo_z, $
                              /from_gei, /to_geo

          cdf_epoch, epoch0, year, month, day, /compute_epoch
          epoch=data_in.utc*3600d0*1d3 + epoch0
          avgEpoch = mean(epoch)
          epoch = fltArr(n_elements(gei_x[*])) + AvgEpoch
          mltShift = aacgm_mlt(year, yrSecAvg, 0.0 )
; conv GEOG XYZ -> GEOG r,thet,phi
          sphcar,geog_R_km, geog_coLat_deg, geog_lon_deg, $
                 geo_x, geo_y, geo_z, /degree ,/to_sphere

          gLat_a=90.0-geog_coLat_deg         ; convert to Latitude
          aacgm_conv_vec, gLat_a, geog_lon_deg, geog_R_km-Re_km, $
                    data_in.dbth, data_in.dbph, $
     mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph,err,/to_aacgm

; Plot AACGM data
          window, 3, xSize = 500, ySize = 500,title='AACGM-MLT'
          If HemSp eq 'S' Then 	$
           plt_dat, 90.0-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	                np, mth_vec_gth, mph_vec_gph, [1,1], [1,2], title =st_str
          If HemSp eq 'N' Then $
           plt_dat,90.-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
                    np,-mth_vec_gth, mph_vec_gph, [1,1], [1,2], title=st_str
         end
 'prv' : Begin
          st_tme_sec=float(StHr)*3600.0+float(StMn)*60.0+StSc
          st_tme_sec=st_tme_sec-TmSep
          en_tme_sec=st_tme_sec+GetSc
          If st_tme_sec lt 0.0 Then begin
           st_tme_sec=long(0)
           en_tme_sec=GetSc
          end
          st_tme=st_tme_sec/3600.0
          en_tme=en_tme_sec/3600.0
          if HemSp eq 'N' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz gt 0., iiCnt) ; set for Nth Hemis with data.pz
          if HemSp eq 'S' then iiTime=where(rw_data.utc gt st_tme and rw_data.utc lt en_tme and rw_data.pz lt 0., iiCnt) ; set for Sth Hemis with data.pz
          data_in = rw_data[iiTime]
          np=n_elements(data_in)
          StHr=fix(st_tme_sec/3600.0)
          StMn=fix((st_tme_sec-float(StHr)*3600.0)/60.0)
          StSc=st_tme_sec-float(StHr)*3600.0-float(StMn)*60.0
          Widget_Control,HrSld,Set_Value=StHr
          Widget_Control,MnSld,Set_Value=StMn
          Widget_Control,ScSld,Set_Value=StSc
          tmestr,StHr,StMn,StSc,st_str
; Get spherical coords of XYZ locations
          sphcar,rg_th,cglat,glon,data_in.px,data_in.py,data_in.pz,/to_sphere,/degree   ; find GEO(r,thet,phi)
          data_in[*].pr=rg_th
          data_in[*].pth=cglat
          data_in[*].pph=glon
; Put the x,y,z GEO input db data into spherical components
          For ii=Long(0),Long(np)-Long(1) do begin
           bcarsp, data_in[ii].px, data_in[ii].py, data_in[ii].pz, data_in[ii].dbx, data_in[ii].dby, data_in[ii].dbz, vbr,vbth,vbph
            data_in[ii].dbr=vbr
            data_in[ii].dbth=vbth
            data_in[ii].dbph=vbph
          end
          If HemSp eq 'S' then begin
           idx=where(cglat gt 90.)             ; Find Sth hemisphere points
           If (idx(0) gt -1) then cglat[idx]=180.0-cglat[idx] ; deal with Sth hemis
          end
; plot the db vectors
          window, 1, xSize = 500, ySize = 500,title='Input Data'
          If HemSp eq 'S' Then plt_dat,cglat,glon,np, data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
          If HemSp eq 'N' Then plt_dat,cglat,glon,np,-data_in.dbth, data_in.dbph,[1,1],[1,2],title=st_str
; Take GEI -> GEOG -> AACGM
          Re_km = 6371.0
          R_km = 110.0 + Re_km
; conv GEI r,thet,phi -> GEI XYZ
          sphcar,data_in.pr, data_in.pth, data_in.pph,$
           gei_x, gei_y, gei_z, /degree, /to_rect
; set up GEOPACK
          geoPack_reCalc, year, month, day, /date
          avgHour	= fix ( (st_tme + en_tme) / 2.0 )
          avgMin	= ( (st_tme + en_tme) / 2.0 mod 1 ) * 60
          yrSecAvg	= cnvTime(year, month, day, avgHour, avgMin, 0.0)
; conv GEI XYZ -> GEOG XYZ
          geoPack_conv_coord, gei_x, gei_y, gei_z, $
                              geo_x, geo_y, geo_z, $
                              /from_gei, /to_geo

          cdf_epoch, epoch0, year, month, day, /compute_epoch
          epoch=data_in.utc*3600d0*1d3 + epoch0
          avgEpoch = mean(epoch)
          epoch = fltArr(n_elements(gei_x[*])) + AvgEpoch
          mltShift = aacgm_mlt(year, yrSecAvg, 0.0 )
; conv GEOG XYZ -> GEOG r,thet,phi
          sphcar,geog_R_km, geog_coLat_deg, geog_lon_deg, $
                 geo_x, geo_y, geo_z, /degree ,/to_sphere

          gLat_a=90.0-geog_coLat_deg         ; convert to Latitude
          aacgm_conv_vec, gLat_a, geog_lon_deg, geog_R_km-Re_km, $
                    data_in.dbth, data_in.dbph, $
     mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph,err,/to_aacgm

; Plot AACGM data
          window, 3, xSize = 500, ySize = 500,title='AACGM-MLT'
          If HemSp eq 'S' Then 	$
           plt_dat, 90.0-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
	                np, mth_vec_gth, mph_vec_gph, [1,1], [1,2], title =st_str
          If HemSp eq 'N' Then $
           plt_dat,90.-abs(mlat_a), ((mlon_a[*]/15.0+mltShift) mod 24)*15.0, $
                    np,-mth_vec_gth, mph_vec_gph, [1,1], [1,2], title=st_str
         end
 'hms' : Begin
          Widget_Control,ev.id,Get_Value=Hm
          If Hm eq 0 then HemSp='N'
          If Hm eq 1 then HemSp='S'
         end
 'dne' : Begin
          Widget_Control,ev.top,/DESTROY
;          PrmFile='wirid.par'
;          OpenW,u1,PrmFile,/Get_Lun
;          PrintF,u1,GetSc,TmSep,NumPlts
;          PrintF,u1,sig,bsatsm,thrshld,L,MLim
;          If (HemSp eq 'N') Then Hm=0
;          If (HemSp eq 'S') Then Hm=1
;          PrintF,u1,PltFAC,InsOpt,StrTrack,Hm
;          Free_Lun,u1
         end
 end
end
;
; --------------------------------------------------------------------------------------------------
;
Pro AmpPar_mag
 Common WidgBlk,base_in,FNme,DFNm,NxtBut,PrvBut,HrSld,MnSld,ScSld,DneBut
 Common WValBlk,InF,HemSp,StHr,StMn,StSc,GetSc,TmSep,path

 base_in= Widget_Base(/COLUMN,XSize=250,title='AMPERE dB Plot')
 FNme  = Widget_Button(base_in,value='Select Input IDL Save File',UValue='gfi')
 DFNm  = Widget_Text(base_in,value='File: ',UValue='dfn')
 hArr=StrArr(2)
 Hm=1
 hArr(0)='N'
 hArr(1)='S'
 HmSB = CW_BGroup(base_in,hArr,UValue='hms',Label_top='Hemisphere',Set_Value=Hm,/row,/Exclusive)
 HemSp=hArr(Hm)
 HrSld = Widget_Slider(base_in,value=StHr,UValue='shr',Minimum=0,Maximum=23,Title='Start Hour')
 Widget_Control,HrSld,Sensitive=0
 MnSld = Widget_Slider(base_in,value=StMn,UValue='smn',Minimum=0,Maximum=59,Title='Start Minute')
 Widget_Control,MnSld,Sensitive=0
 ScSld = Widget_SLider(base_in,value=StSc,UValue='ssc',Minimum=0,Maximum=59,Title='Start Second')
 Widget_Control,ScSld,Sensitive=0
 AvSld = Widget_SLider(base_in,value=GetSc,UValue='asc',Minimum=15,Maximum=7200,Title='Display Time (seconds)')
 IncSld= Widget_SLider(base_in,value=TmSep,UValue='ins',Minimum=10,Maximum=3600,Title='Seconds to Shift')
 NxtBut=Widget_Button(base_in,value='Next Plot',UValue='nxt')
 Widget_Control,NxtBut,Sensitive=0
 PrvBut=Widget_Button(base_in,value='Previous Plot',UValue='prv')
 Widget_Control,PrvBut,Sensitive=0
 DneBut = Widget_Button(base_in,value='Exit Program',UValue='dne')
 WIDGET_CONTROL, base_in, /REALIZE
 XMANAGER, 'AmpPar_mag',base_in
end
;
; ------------------------------------------------------------------------------------------------------
;
pro ampere_plot_aacgm
 Common WValBlk,InF,HemSp,StHr,StMn,StSc,GetSc,TmSep,path

 if strCmp (!version.os, 'Win32' ) then begin  ; win macine for Colin
  plotDev = 'win'
  path = 'd:\cwac\hi_res\davidg\'
 endif else begin                            ; Linux machine for Dave
  plotDev = 'X'
  path = '~/code/ampereFit/idl/'
 endelse
; set plot system
 set_plot, !d.name
 device,decomposed = 0
 !p.background = 255
; Initialise some values
 StHr=0 & StMn=20 & StSc=0
 GetSc=600 & TmSep=600
 plotCapSize = 40.0
 AmpPar_mag
 !p.background = 0
 Print,'Finished'
end
