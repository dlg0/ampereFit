; Procedure to centre AMPERE data orbit plane tracks and calc xyz rotation matrix
; Full sphere version (CLW - Sept 2010)
;
; Input : data_in (data structure)
;         'south' -> if set, calc rotation using Sth hemisphere orbit plane intersection
; output : data_out (data structure)
;          rot_mat : rotation matrix for GEI XYZ to centre orbit intersections
;          diag=diag -> print diagnostic messages
;
; Data structure:
;  data = { ipln: 0, isat: 0, utc: 0d0, $
;           px: 0d0, py: 0d0, pz: 0d0, $
;           dbx: 0d0, dby: 0d0, dbz: 0d0, $
;        GEI_R_km : 0d0, GEI_coLat_deg: 0d0, GEI_lon_deg: 0d0, $
;        GEI_coLat_rad: 0d0, GEI_lon_rad: 0d0, $
;           br_GEI: 0d0, bTheta_GEI: 0d0, bPhi_GEI: 0d0 }
;
; C.L. Waters
; Centre for Space Pysics
; University of Newcastle, Australia
; Dec 2009
;
; Comments:
;
; Modifications:
;  Added Vec_rot and rot_ang to output list so that the rotation can be reversed [Dec, 2009]
;  Deleted coord and vec roation section - put into AMPERE_SHIFT_FULL.PRO
;
Pro ampere_get_rotmat_full, data_in, rot_mat, south, diag=diag
; If south = 1 then rotate to south hemisphere intersection else use Nth hemisphere
;
; Error check
  status=0
  catch,error_status
  if (error_status ne 0) then begin
    status=1
    Print,'Error in ampere_get_rotmat_full; returning...'
    catch,/cancel
    return
  endif
; 
; Estimate the average track intersection location
; Step 1: Calc x,y coords (funky cylindrical) of Iridium data locations
  npl=6                     ; 6 orbit planes
  mnpnts=20                 ; Min number of points in a track to proceed

  if south eq 0 then us_idx=where(data_in.GEI_coLat_deg lt 90.)
  if south eq 1 then us_idx=where(data_in.GEI_coLat_deg ge 90.)

; Get a copy of co_lat data
  clat_a=data_in[us_idx].GEI_coLat_deg               ; co_Lat in degrees [0->180]
  trnum_a=data_in[us_idx].iPln
  
; Calc x,y coords from co_lat long location data
  xcrd_a=clat_a*cos(data_in[us_idx].GEI_lon_rad)
  ycrd_a=clat_a*sin(data_in[us_idx].GEI_lon_rad)

; Count number of data points in each Iridium orbit track
  p_cnt_pl=intarr(npl)                    ; space to store num points per orbit plane
  For mp=0,npl-1 do begin                 ; Get Merid Planes in order, ipln runs from 0->5
    idx1=where(trnum_a eq mp,cnt1)
    p_cnt_pl[mp]=cnt1                   ; number of data points per plane
  end
  If keyword_set(diag) then Print,'Num. Points read from Input File=',n_elements(data_in)
  If keyword_set(diag) then Print,'Num. Points per Orbit track:',p_cnt_pl

; Check for sufficient data on all orbit tracks
  If (min(p_cnt_pl) le mnpnts) then begin
    status=1                                  ; Return if insufficient data on any track
    Print,'Insufficient data on Orbit Track'
    return
;  stop
  endif

; Each orbit track is assumed to be a parabola curve in x,y
; Calc the coeffs of the parabolic eqn of each orbit track (in x,y)
  pcoef=dblarr(NPl,3)
  t_cnt=0

  For mp=0,NPl-1 do begin
; WARNING : SVDFIT crashes if we get a straight track along the 0600-1800 line: CLW
    idx=where(trnum_a eq mp)                          ; select out each orbit track
    res_a=SVDFIT(xcrd_a[idx],ycrd_a[idx],3,status=status)  ; parabola fit to track location data
    pcoef(t_cnt,*)=res_a                                   ; store parabola eqn coefficients
    t_cnt++
    if status ne 0 then begin
      print,'Error: ampere_get_rotmat_full - Track intersection solver, excluding track'
      t_cnt=t_cnt-1
    end
  end
;
; Now calculate the intersection coords of these orbit track parabolas  
  n_com=total(indgen(t_cnt))
  XInt=dblarr(n_com)                  ; There are n_com orbit track intersection combinations
  YInt=XInt & IntCnt=0
; Loop through every combination of the 6 planes to find intersections
  For aa=0,t_cnt-1 do begin
    For bb=1,t_cnt-1 do begin
      If (bb gt aa) then begin
        adiff=pcoef(aa,2)-pcoef(bb,2)
        bdiff=pcoef(aa,1)-pcoef(bb,1)
        cdiff=pcoef(aa,0)-pcoef(bb,0)
        sqtrm=bdiff^2-4.*adiff*cdiff
        If (sqtrm ge 0.) then begin           ; trap for -ve in sqrt of quadratic eqn
          x1=0.5*(-bdiff+sqrt(bdiff^2-4.*adiff*cdiff))/adiff
          x2=0.5*(-bdiff-sqrt(bdiff^2-4.*adiff*cdiff))/adiff
          y1=pcoef(aa,2)*x1*x1+pcoef(aa,1)*x1+pcoef(aa,0)
          y2=pcoef(aa,2)*x2*x2+pcoef(aa,1)*x2+pcoef(aa,0)
          RDst1=sqrt(x1*x1+y1*y1)
          RDst2=sqrt(x2*x2+y2*y2)
          If (RDst2 gt RDst1) then begin   ; Find the correct quadratic solution
            XInt(IntCnt)=x1
            YInt(IntCnt)=y1
          end
          If (RDst1 gt RDst2) then begin
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

; Release some memory
  us_idx=!null
  clat_a=!null
  trnum_a=!null
  xcrd_a=!null
  ycrd_a=!null

  If keyword_set(diag) then begin
    Print,'In ampere_get_rotmat_full: Coordinate Shift is [x,y]: ',XAvg,YAvg
    Print,'Intersect X Values: ',XInt
    Print,'Intersect Y Values: ',YInt
    Print,'PCoeffs for each Track [x^2, x, c] : '
    For aa=0,npl-1 do print,transpose(pcoef(aa,*))
  end
;
  CRLat=sqrt(XAvg*XAvg+YAvg*YAvg)             ; Find Lat,Lon corresponding to XAvg,YAvg
  RLon=180.0+atan(YAvg,XAvg)*180.0/!dpi
  If south then CRLat=180.0-CRLat             ; Report CRLat in correct Hemisphere
  If keyword_set(diag) then Print,'Shift in Lat,Lon [deg]=',CRLat,RLon

; Given xsh, ysh - need to calc rotation matrix that shifts coords so xsh.ysh is the origin
; Get rotation axis (perp to plane defined by z and shifted point and the origin)
  r0=mean(sqrt(data_in.px^2+data_in.py^2+data_in.pz^2))              ; mag of radial vector, always +ve
; sphcar,abs(r0),CRLat,RLon,x_sh,y_sh,z_sh,/to_rect,/degree          ; XYZ coord of new pole in old system (in correct Hemis)
  geopack_sphcar,abs(r0),CRLat,RLon,x_sh,y_sh,z_sh,/to_rect,/degree  ; XYZ coord of new pole in old system (in correct Hemis)
  sh_mag=sqrt(x_sh^2+y_sh^2+z_sh^2)                                  ; magnitude of shifted pole vector
  If south then r0=-r0                                               ; point r vec along Z axis in correct direction for Hemis
  cross_p,[0.0,0.0,r0],[x_sh,y_sh,z_sh],v_rot                        ; rotation is about this axis, v_rot should always have a zero z component
  c_arg=(r0*z_sh)/(sh_mag*abs(r0))
  rot_ang=acos(c_arg)                                                ; rotation angle, calc from dot product
  If keyword_set(diag) then Print,'Rotation off input GEI Pole : ',rot_ang*180.0/!dpi,' deg'
  v_rot_mag=sqrt(v_rot(0)^2+v_rot(1)^2+v_rot(2)^2)                   ; magnitude of rotation vector
  urx=v_rot(0)/v_rot_mag                                             ; unit rotation axis vector (perp to plane containing Z and shifted origin)
  ury=v_rot(1)/v_rot_mag
  urz=v_rot(2)/v_rot_mag
; Matrices for rotation
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
end
