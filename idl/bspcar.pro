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
