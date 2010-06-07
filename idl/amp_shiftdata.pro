; Procedure to centre AMPERE data orbit plane tracks
; Input : data_in (data structure)
;         'south' -> if set, data is for Sth hemisphere
; output : data_out (data structure)
;          rot_mat : rotation matrix for GEI XYZ to centre orbit intersections
;          south=south -> data is Sth Hemis
;          diag=diag -> print diagnostic messages
;
; Data structure:
;          data_struct = {ipln:0, sv:0, stme:0d0,$
;          px:0d0, py:0d0, pz:0d0,$
;          dbx:0d0, dby:0d0, dbz:0d0,$
;          pr:0d0, pth:0d0, pph:0d0,$
;          br_GEI:0d0, bTheta_GEI:0d0, bPhi_GEI:0d0,$
;          xcrd:0d0, ycrd:0d0}
;
; C.L. Waters
; Centre for Space Pysics
; University of Newcastle, Australia
; Dec 2009
;
; Comments:
;  Will need to deal with any cross hemisphere shift - e.g. full sphere fits
;
; Modifications:
;   Added Vec_rot and rot_ang to output list so that the rotation can be reversed [Dec, 2009]
;
;
Pro amp_shiftdata,data_in,data_out,rot_mat,south,diag=diag
; Estimate the average track intersection location
; Step 1: Calc x,y coords (funky cylindrical) of locations
 npl=6                     ; 6 orbit planes
 mnpnts=20                 ; Min number of points in a track to proceed

; Make Co_Lat with 0 at centre (for track fit)
 clat_a=data_in.GEI_coLat_deg               ; co_Lat in degrees

 idx=where(data_in.GEI_coLat_deg gt 90.)    ; Sth hemisphere data
 if idx(0) gt -1 then clat_a[idx]=180.0-data_in.GEI_coLat_deg
 xcrd_a=clat_a*cos(data_in.GEI_lon_deg*!dpi/180.0)
 ycrd_a=clat_a*sin(data_in.GEI_lon_deg*!dpi/180.0)
;
 p_cnt_pl=intarr(npl)                    ; space to store num points per orbit plane
 For mp=0,npl-1 do $                     ; Get Merid Planes in order, ipln runs from 0->5
 Begin
  idx1=where(data_in.ipln eq mp,cnt1)
  p_cnt_pl[mp]=cnt1
 end
 If keyword_set(diag) then Print,'Num. Points read from Input File=',n_elements(data_in)
 If keyword_set(diag) then Print,'Num. Points per Orbit track:',p_cnt_pl

; Check for sufficient data on all orbit tracks
 If (min(p_cnt_pl) le mnpnts) then $
 Begin
  status=1                                  ; Return if insufficient data on any track
;  return
  Print,'Insufficient data on Orbit Track'
  stop
 endif

; Calc the coeffs of the parabolic eqn of each orbit track (in x,y)
 pcoef=dblarr(NPl,3)
 t_cnt=0
 For mp=0,NPl-1 do begin
; WARNING : SVDFIT crashes if we get a straight track along the 0600-1800 line: CLW
  idx=where(data_in.ipln eq mp)                     ; select out each orbit track
  res_a=SVDFIT(xcrd_a[idx],ycrd_a[idx],3,status=status)  ; parabola fit to each track location data
;print,'t_cnt=',t_cnt
  pcoef(t_cnt,*)=res_a                              ; store parabola eqn coefficients
  t_cnt++
  if status ne 0 then begin
   print,'WARNING: amp_shiftdata - Track intersection solver, excluding track'
   t_cnt=t_cnt-1
  end
 end
 n_com=total(indgen(t_cnt))
 XInt=dblarr(n_com)                  ; There are n_com orbit track intersection combinations
 YInt=XInt
 IntCnt=0
; Loop through every combination of the 6 planes to find intersections
 For aa=0,t_cnt-1 do $
 Begin
  For bb=1,t_cnt-1 do $
  Begin
   If (bb gt aa) then $
   Begin
    adiff=pcoef(aa,2)-pcoef(bb,2)
    bdiff=pcoef(aa,1)-pcoef(bb,1)
    cdiff=pcoef(aa,0)-pcoef(bb,0)
    sqtrm=bdiff^2-4.*adiff*cdiff
    if (sqtrm ge 0.) then $           ; trap for -ve in sqrt of quadratic eqn
    Begin
     x1=0.5*(-bdiff+sqrt(bdiff^2-4.*adiff*cdiff))/adiff
     x2=0.5*(-bdiff-sqrt(bdiff^2-4.*adiff*cdiff))/adiff
     y1=pcoef(aa,2)*x1*x1+pcoef(aa,1)*x1+pcoef(aa,0)
     y2=pcoef(aa,2)*x2*x2+pcoef(aa,1)*x2+pcoef(aa,0)
     RDst1=sqrt(x1*x1+y1*y1)
     RDst2=sqrt(x2*x2+y2*y2)
     If (RDst2 gt RDst1) then $   ; Find the correct quadratic solution
     Begin
      XInt(IntCnt)=x1
      YInt(IntCnt)=y1
     end
     If (RDst1 gt RDst2) then $
     Begin
      XInt(IntCnt)=x2
      YInt(IntCnt)=y2
     end
     IntCnt++
    end    ; +ve quadratic sqrt
   end     ; if bb gt aa -> MPlane combinations
  end      ; bb MPlane Loop
 end       ; aa MPlane Loop

; Calc equally weighted, average location of track intersections
 XAvg=mean(XInt)
 YAvg=mean(YInt)
; ********** Turn off shift
; XAvg=0.0 & YAvg=0.0
; ************************

 If keyword_set(diag) then Print,'Intersect X Values: ',XInt
 If keyword_set(diag) then Print,'Intersect Y Values: ',YInt
 If keyword_set(diag) then begin
  Print,'PCoeffs for each Track [x^2, x, c] : '
  For aa=0,npl-1 do print,transpose(pcoef(aa,*))
 end
; If keyword_set(diag) then Print,'Coordinate Shift is [x,y]: ',XAvg,YAvg
 Print,'In amp_shiftdata: Coordinate Shift is [x,y]: ',XAvg,YAvg
 CRLat=sqrt(XAvg*XAvg+YAvg*YAvg)                   ; Find Lat,Lon corresponding to XAvg,YAvg
 RLon=180.0+atan(YAvg,XAvg)*180.0/!dpi
 If south then CRLat=180.0-CRLat            ; Report CRLat in correct Hemisphere
 If keyword_set(diag) then Print,'Shift in Lat,Lon [deg]=',CRLat,RLon

; Get rotation axis (perp to plane defined by z and shifted point and the origin)
 r0=mean(sqrt(data_in.px^2+data_in.py^2+data_in.pz^2))      ; always +ve
; sphcar,abs(r0),CRLat,RLon,x_sh,y_sh,z_sh,/to_rect,/degree  ; XYZ coord of new pole in old system (in correct Hemis)
 geopack_sphcar,abs(r0),CRLat,RLon,x_sh,y_sh,z_sh,/to_rect,/degree  ; XYZ coord of new pole in old system (in correct Hemis)
 sh_mag=sqrt(x_sh^2+y_sh^2+z_sh^2)                 ; magnitude of shifted pole vector
 If south then r0=-r0                 ; point r vec along Z axis in correct direction for Hemis
 cross_p,[0.0,0.0,r0],[x_sh,y_sh,z_sh],v_rot       ; rotation is about this axis, v_rot should always have a zero z component
 c_arg=(r0*z_sh)/(sh_mag*abs(r0))
 rot_ang=acos(c_arg)                               ; rotation angle, from dot product
 If keyword_set(diag) then Print,'Rotation off Pole : ',rot_ang*180.0/!dpi
 v_rot_mag=sqrt(v_rot(0)^2+v_rot(1)^2+v_rot(2)^2)
 urx=v_rot(0)/v_rot_mag                            ; unit rotation axis vector (perp to plane containing Z and shifted origin)
 ury=v_rot(1)/v_rot_mag
 urz=v_rot(2)/v_rot_mag
 p_a=[[urx*urx, urx*ury, urx*urz],$                ; matrices for Rodrigues rotation formula; R=P+(I-P)cos(thet)+Qsin(thet)
      [urx*ury, ury*ury, ury*urz],$
      [urx*urz, ury*urz, urz*urz]]
 i_a=[[1.0, 0.0, 0.0],$                            ; identity matrix
      [0.0, 1.0, 0.0],$
      [0.0, 0.0, 1.0]]
 q_a=[[0.0, -urz, ury],$
      [urz, 0.0, -urx],$
      [-ury, urx, 0.0]]
 rot_mat = p_a + (i_a-p_a)*cos(rot_ang) + q_a*sin(rot_ang)  ; IDL orders by [column,row]

 data_out=data_in
 np=n_elements(data_out)

 For ii=Long(0),Long(np)-1 do begin
; rotate coord locations
  data_out[ii].px=data_in[ii].px*rot_mat[0,0] + $
                  data_in[ii].py*rot_mat[1,0] + $
                  data_in[ii].pz*rot_mat[2,0]
  data_out[ii].py=data_in[ii].px*rot_mat[0,1] + $
                  data_in[ii].py*rot_mat[1,1] + $
                  data_in[ii].pz*rot_mat[2,1]
  data_out[ii].pz=data_in[ii].px*rot_mat[0,2] + $
                  data_in[ii].py*rot_mat[1,2] + $
                  data_in[ii].pz*rot_mat[2,2]
; rotate field components
  data_out[ii].dbx=data_in[ii].dbx*rot_mat[0,0] + $
                   data_in[ii].dby*rot_mat[1,0] + $
                   data_in[ii].dbz*rot_mat[2,0]
  data_out[ii].dby=data_in[ii].dbx*rot_mat[0,1] + $
                   data_in[ii].dby*rot_mat[1,1] + $
                   data_in[ii].dbz*rot_mat[2,1]
  data_out[ii].dbz=data_in[ii].dbx*rot_mat[0,2] + $
                   data_in[ii].dby*rot_mat[1,2] + $
                   data_in[ii].dbz*rot_mat[2,2]
; convert to spherical ready for fitting
  geopack_bcarsp, data_out[ii].px, data_out[ii].py, data_out[ii].pz, $
          data_out[ii].dbx, data_out[ii].dby, data_out[ii].dbz, vbr,vbth,vbph
  data_out[ii].br_GEI=vbr
  data_out[ii].bTheta_GEI=vbth
  data_out[ii].bPhi_GEI=vbph
 end
 geopack_sphcar,data_out.px,data_out.py,data_out.pz,rg_th,cglat,glon,/to_sphere,/degree   ; find shifted GEI(r,thet,phi)
 data_out.GEI_R_km=rg_th
 data_out.GEI_coLat_deg=cglat
 data_out.GEI_lon_deg=glon
end
