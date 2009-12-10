Pro cross_p,v1,v2,vout
; Calc vector cross product
; C.L. Waters, Dec 2009
;
 v1=reform(v1) & v2=reform(v2)
 vout=[[v1(1)*v2(2)-v1(2)*v2(1)],$
       [v1(2)*v2(0)-v1(0)*v2(2)],$
       [v1(0)*v2(1)-v1(1)*v2(0)]]
end
