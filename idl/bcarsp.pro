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
;
;  Conv to IDL : C.L. Waters
;
Pro bcarsp,x,y,z,vx,vy,vz,vr,vth,vph
 RHO2=x^2+y^2
 R=SQRT(RHO2+z^2)
 RHO=SQRT(RHO2)
 If (RHO ne 0.0) then begin   ; check for location at the poles
  CPHI=X/RHO
  SPHI=Y/RHO
 end else begin
  CPHI=1.0
  SPHI=0.0
 end
 CT=Z/R         ; works with abs(Z)/R
 ST=RHO/R
 vr=(X*vx+Y*vy+Z*vz)/R
 vth=(vx*CPHI+vy*SPHI)*CT-vz*ST
 vph=vy*CPHI-vx*SPHI
END
