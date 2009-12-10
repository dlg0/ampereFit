; -----------------------------------------------------------------------------------
;      SUBROUTINE SPHCAR (R,THETA,PHI,X,Y,Z)
;   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
;
;  NOTE: AT THE POLES (X=0 AND Y=0) WE ASSUME PHI=0 (WHEN CONVERTING
;        FROM CARTESIAN TO SPHERICAL COORDS, I.E., FOR J<0)
;   AUTHOR:  N. A. TSYGANENKO
; Converted to IDL: C.L. Waters
;
Pro sphcar,r_a,thet_a,phi_a,x_a,y_a,z_a,to_rect=to_rect,to_sphere=to_sphere,degree=degree
; Conv from XYZ to Spherical
 If keyword_set(to_sphere) then begin
  np=Long(n_elements(x_a))
  phi_a=dblarr(np)
  thet_a=dblarr(np)
  sq=x_a^2+y_a^2
  r_a=SQRT(SQ+z_a^2)         ; have radial array
  For mm=Long(0),Long(np)-Long(1) do begin
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
  np=Long(n_elements(r_a))
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
