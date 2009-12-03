; Code to convert Geographic vector (velocity, dB etc.) to AACGM by moving in dth then dphi GEO
; The conversion process uses XYZ to ensure correct length scaling. Output are the contra and covariant vectors
;
; NB: The AACGM Lat,Long coords are non-orthogonal
;
; Inputs:
; aacgm_path : dir for AACGM coeffs
; gLat_a,gLon_a [dblarr(n_Lat*n_Long)] : locations of vectors in GEO (theta,phi) in degrees (lat=90 is the north pole)
; gVec_th,gVec_ph [dblarr(n_Lat*n_Long)] : GEO vector data in component form.
;
; Outputs:
; mlat_a=dblarr(n_Lat*n_Long)    AACGM location information
; mlon_a=dblarr(n_Lat*n_Long)
; mth_vec_gth=dblarr(n_Lat*n_Long)    Vec for dth shift in GEO
; mph_vec_gth=dblarr(n_Lat*n_Long)
; mth_vec_gph=dblarr(n_Lat*n_Long)    Vec for dph shift in GEO
; mph_vec_gph=dblarr(n_Lat*n_Long)
;
; C.L. Waters,
; Centre for Space Physics
; University of Newcastle
; New South Wales, Australia
; Nov 2009
;
; Modifications:
;
; -----------------------------------------------------------------------------------
;      SUBROUTINE SPHCAR (R,THETA,PHI,X,Y,Z,J)
;   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
;
;  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 (WHEN CONVERTING
;        FROM CARTESIAN TO SPHERICAL COORDS, I.E., FOR J<0)
;   AUTHOR:  N. A. TSYGANENKO
; Converted to IDL: C.L. Waters
;
Pro sphcar,r_a,thet_a,phi_a,x_a,y_a,z_a,$
           to_rect=to_rect,to_sphere=to_sphere,degree=degree
; Conv from XYZ to Spherical
 If keyword_set(to_sphere) then begin
  np=n_elements(x_a)
  phi_a=dblarr(np)
  thet_a=dblarr(np)
  sq=x_a^2+y_a^2
  r_a=SQRT(SQ+z_a^2)         ; have radial array
  For mm=0,np-1 do begin
   If sq(mm) eq 0. then begin
    phi_a(mm)=0.0
    If z_a(mm) LT 0.0 then thet_a(mm)=!dpi else thet_a(mm)=0.0
   end else begin
    sq(mm)=SQRT(sq(mm))
    phi_a(mm)=atan(y_a(mm),x_a(mm))
    thet_a(mm)=atan(sq(mm),z_a(mm))
    If phi_a(mm) lt 0.0 then phi_a(mm)=phi_a(mm)+2.0*!dpi
   end
  end       ; mm loop
 end        ; if keyword to_sphere
 If keyword_set(to_rect) then begin
  np=n_elements(r_a)
  x_a=dblarr(np)
  y_a=dblarr(np)
  z_a=dblarr(np)
; Check for deg or radians in theta,phi
  If keyword_set(degree) then begin
   thet_a=thet_a*!dpi/180.0
   phi_a=phi_a*!dpi/180.0
  end
  sq=r_a*sin(thet_a)
  x_a=sq*cos(phi_a)
  y_a=sq*sin(phi_a)
  z_a=r_a*cos(thet_a)
 end
 If keyword_set(degree) then begin
  thet_a=thet_a*180.0/!dpi
  phi_a=phi_a*180.0/!dpi
 end
END
;
; --------------------------------------------------------------------------
;
;      SUBROUTINE BSPCAR (THETA,PHI,BR,BTHETA,BPHI,BX,BY,BZ)
;
;   CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
;-----INPUT:   THETA,PHI - SPHERICAL ANGLES OF THE POINT
;              BR,BTHETA,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
;-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
;
;   WRITTEN BY:  N. A. TSYGANENKO
;   Conv to IDL: C.L. Waters
;
Pro bspcar,thet_a,phi_a,vr_a,vth_a,vph_a,vx_a,vy_a,vz_a,degree=degree
 If keyword_set(degree) then begin
  thet_a=thet_a*!dpi/180.0
  phi_a=phi_a*!dpi/180.0
 end
 s=sin(thet_a)
 c=cos(thet_a)
 sf=sin(phi_a)
 cf=cos(phi_a)
 BE=vr_a*S + vth_a*C
 vx_a=BE*CF-vph_a*SF
 vy_a=BE*SF+vph_a*CF
 vz_a=vr_a*C-vth_a*S
 If keyword_set(degree) then begin
  thet_a=thet_a*180.0/!dpi
  phi_a=phi_a*180.0/!dpi
 end
END
;
; -----------------------------------------------------------------------
;
;      SUBROUTINE BCARSP (X,Y,Z,BX,BY,BZ,BR,BTHETA,BPHI)
;CALCULATES SPHERICAL FIELD COMPONENTS FROM THOSE IN CARTESIAN SYSTEM
;
;-----INPUT:   X,Y,Z  - CARTESIAN COMPONENTS OF THE POSITION VECTOR
;              BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD VECTOR
;-----OUTPUT:  BR,BTHETA,BPHI - SPHERICAL COMPONENTS OF THE FIELD VECTOR
;
;  AT THE POLES (THETA=0 OR THETA=PI), ASSUME PHI=0,
;        AND HENCE BTHETA=BX, BPHI=BY
;   AUTHOR:   N. A. TSYGANENKO
;  Conv to IDL : C.L. Waters
;
Pro bcarsp,x,y,z,vx,vy,vz,vr,vth,vph
 RHO2=x^2+y^2
 R=SQRT(RHO2+z^2)
 RHO=SQRT(RHO2)
 If RHO ne 0. then begin
  CPHI=X/RHO
  SPHI=Y/RHO
 end else begin
  CPHI=1.
  SPHI=0.
 end
 CT=Z/R
 ST=RHO/R
 vr=(X*vx+Y*vy+Z*vz)/R
 vth=(vx*CPHI+vy*SPHI)*CT-vz*ST
 vph=vy*CPHI-vx*SPHI
END
;
;---------------------------------------------------------------------
;
Pro cross_p,v1,v2,vout
 v1=reform(v1) & v2=reform(v2)
 vout=[[v1(1)*v2(2)-v1(2)*v2(1)],$
       [v1(2)*v2(0)-v1(0)*v2(2)],$
       [v1(0)*v2(1)-v1(1)*v2(0)]]
 vout=transpose(vout)
end
;
; --------------------------------------------------------------------
;
Pro vec_geo2aacgm,aacgm_coef_path,aacgm_yr,aacgm_hgt,$
                       gLat_a,gLon_a,gVec_th,gVec_ph,$
                       mlat_a,mlon_a,mth_vec_gth,mph_vec_gth,mth_vec_gph,mph_vec_gph
 dgth=1.0 & dgph=1.0               ; shift in thet and phi for unit vecs
 NTOT=n_elements(gLat_a)
; d2r=!dpi/180.0
 mlat_a=dblarr(NTOT)
 mlon_a=dblarr(NTOT)
 mth_vec_gth=dblarr(NTOT)
 mph_vec_gth=dblarr(NTOT)
 mth_vec_gph=dblarr(NTOT)
 mph_vec_gph=dblarr(NTOT)

; Load AACGM coeff data
 aacgm_set_path,aacgm_coef_path    ; set aacgm coeff path
 aacgm_load_coef,aacgm_yr          ; load the coeffs
; Get XYZ coords of GEO coords
 cglat_a=(90.0-glat_a)             ; Co_Latitude in degrees
; Conv GEO(Co_Lat,Long) to GEO(X,Y,Z)
 sphcar,aacgm_hgt,cglat_a,glon_a,xg_a,yg_a,zg_a,/to_rect,/degree
; geopack_sphcar,aacgm_hgt,cglat_a,glon_a,xgg_a,ygg_a,zgg_a,/to_rect,/degree
; Conv GEO(Lat,Long) to AACGM(Lat,Long)
 aacgm_conv_coord,gLat_a,gLon_a,aacgm_hgt,mlat_a,mlon_a,err,/to_aacgm,order=10
 ln_fxi=where(mlon_a lt 0.0)
 if ln_fxi(0) gt -1 then mlon_a(ln_fxi)=mlon_a(ln_fxi)+360.0
; Conv AACGM(CoLat,Long) to XYZ
 cmlat_a=(90.0-mlat_a)
 sphcar,aacgm_hgt,cmlat_a,mlon_a,xm_a,ym_a,zm_a,/to_rect,/degree

 vr=dblarr(NTOT)
 vth=dblarr(NTOT)
 vph=dblarr(NTOT)
;  del_th in GEO
 vr(*)=0.0 & vth(*)=dgth & vph(*)=0.0              ; Construct a (0,theta,0) vector
 bspcar,cglat_a,glon_a,vr,vth,vph,vx,vy,vz,/degree ; Get this vector in XYZ components
; mg=sqrt(vx^2+vy^2+vz^2)                           ; Make this a unit vector
; vx=vx/mg & vy=vy/mg & vz=vz/mg
 xg_th=xg_a+vx & yg_th=yg_a+vy & zg_th=zg_a+vz     ; Get GEO(X,Y,Z) for GEO(0,dth,0) (unit vec) shift
 sphcar,rg_th,cglat_th,glon_th,xg_th,yg_th,zg_th,/to_sphere,/degree
 glat_th=90.0-cglat_th
;  Convert to AACGM mlat,mlon coords (from geo+dth)
 aacgm_conv_coord,glat_th,glon_th,aacgm_hgt,mlat_th,mlon_th,err,/to_aacgm,order=10
 ln_fxi=where(mlon_th lt 0.0)
 if ln_fxi(0) gt -1 then mlon_th(ln_fxi)=mlon_th(ln_fxi)+360.0
 cmlat_th=(90.0-mlat_th)  ; get AACGM colatitude
;  Get XYZ coords of AACGM+dth(geo) coords
 sphcar,aacgm_hgt,cmlat_th,mlon_th,xm_th,ym_th,zm_th,/to_rect,/degree

;  del_ph in GEO
 vr(*)=0.0 & vth(*)=0.0 & vph(*)=dgph              ; Construct a (0,0,phi) vector
 bspcar,cglat_a,glon_a,vr,vth,vph,vx,vy,vz,/degree ; Get this vector in XYZ coords
; mg=sqrt(vx^2+vy^2+vz^2)                           ; Make this a unit vector
; vx=vx/mg & vy=vy/mg & vz=vz/mg
 xg_ph=xg_a+vx & yg_ph=yg_a+vy & zg_ph=zg_a+vz     ; Get GEO(X,Y,Z) for GEO(0,0,dph) (unit vector) shift
 sphcar,rg_ph,cglat_ph,glon_ph,xg_ph,yg_ph,zg_ph,/to_sphere,/degree
 glat_ph=90.0-cglat_ph
;  Convert to AACGM mlat,mlon coords (from geo+dph)
 aacgm_conv_coord,glat_ph,glon_ph,aacgm_hgt,mlat_ph,mlon_ph,err,/to_aacgm,order=10
 ln_fxi=where(mlon_ph lt 0.0)
 if ln_fxi(0) gt -1 then mlon_ph(ln_fxi)=mlon_ph(ln_fxi)+360.0
 cmlat_ph=(90.0-mlat_ph)  ; get AACGM colatitude
;  Get XYZ coords of AACGM+dph(geo) coords
 sphcar,aacgm_hgt,cmlat_ph,mlon_ph,xm_ph,ym_ph,zm_ph,/to_rect,/degree

; Get aacgm radial unit vector
 mxyz_r=[[xm_a],[ym_a],[zm_a]]
 mxyz_mg=sqrt(mxyz_r(*,0)^2+mxyz_r(*,1)^2+mxyz_r(*,2)^2)
 mxyz_ruv=mxyz_r
 mxyz_ruv(*,0)=mxyz_r(*,0)/mxyz_mg    ; this forces IDL to create a (*,3) array
 mxyz_ruv(*,1)=mxyz_r(*,1)/mxyz_mg
 mxyz_ruv(*,2)=mxyz_r(*,2)/mxyz_mg
; Get AACGM(X,Y,Z) unit vec for a GEO dth shift
 mxyz_th=[[xm_th-xm_a],[ym_th-ym_a],[zm_th-zm_a]]
 mxyz_mg=sqrt(mxyz_th(*,0)^2+mxyz_th(*,1)^2+mxyz_th(*,2)^2)
 mxyz_thuv=mxyz_th
 mxyz_thuv(*,0)=mxyz_th(*,0)/mxyz_mg    ; this forces IDL to create a (*,3) array
 mxyz_thuv(*,1)=mxyz_th(*,1)/mxyz_mg
 mxyz_thuv(*,2)=mxyz_th(*,2)/mxyz_mg
; Get AACGM(X,Y,Z) unit vec for a GEO dph shift
 mxyz_ph=[[xm_ph-xm_a],[ym_ph-ym_a],[zm_ph-zm_a]]
 mxyz_mg=sqrt(mxyz_ph(*,0)^2+mxyz_ph(*,1)^2+mxyz_ph(*,2)^2)
 mxyz_phuv=mxyz_ph
 mxyz_phuv(*,0)=mxyz_ph(*,0)/mxyz_mg    ; this forces IDL to create a (*,3) array
 mxyz_phuv(*,1)=mxyz_ph(*,1)/mxyz_mg
 mxyz_phuv(*,2)=mxyz_ph(*,2)/mxyz_mg
 mxyz_ph_gth=dblarr(NTOT,3)  ; initialise cross prod result arrays
 mxyz_th_gph=dblarr(NTOT,3)

; Loop over number of points to do cross products (get orthogonal directions)
 For kk=0,NTOT-1 do begin
  cross_p,mxyz_ruv(kk,*),mxyz_thuv(kk,*),res  ; get dph unit vector for r x dth
  mxyz_ph_gth(kk,*)=res
; For a GEO dth shift -> AACGM_dth
  bcarsp,xm_a(kk),ym_a(kk),zm_a(kk),mxyz_thuv(kk,0),mxyz_thuv(kk,1),mxyz_thuv(kk,2),mrv_th,mthv_th,mphv_th
; For a GEO dth shift -> AACGM_dph
  bcarsp,xm_a(kk),ym_a(kk),zm_a(kk),mxyz_ph_gth(kk,0),mxyz_ph_gth(kk,1),mxyz_ph_gth(kk,2),mrv_ph,mthv_ph,mphv_ph
  mth_vec_gth(kk)=gVec_th(kk)*mthv_th+gVec_ph(kk)*mthv_ph
  mph_vec_gth(kk)=gVec_th(kk)*mphv_th+gVec_ph(kk)*mphv_ph

  cross_p,mxyz_phuv(kk,*),mxyz_ruv(kk,*),res  ; get dth unit vector for dph x r
  mxyz_th_gph(kk,*)=res
; For a GEO dph shift -> AACGM_dth
  bcarsp,xm_a(kk),ym_a(kk),zm_a(kk),mxyz_th_gph(kk,0),mxyz_th_gph(kk,1),mxyz_th_gph(kk,2),mrv_th,mthv_th,mphv_th
; For a GEO dph shift -> AACGM_dph
  bcarsp,xm_a(kk),ym_a(kk),zm_a(kk),mxyz_phuv(kk,0),mxyz_phuv(kk,1),mxyz_phuv(kk,2),mrv_ph,mthv_ph,mphv_ph
  mth_vec_gph(kk)=gVec_th(kk)*mthv_th+gVec_ph(kk)*mthv_ph
  mph_vec_gph(kk)=gVec_th(kk)*mphv_th+gVec_ph(kk)*mphv_ph
 end
end