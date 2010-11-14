;	ampere_Legendre_f.pro
; Full sphere version (CLW Sept 2010)
; ;
;	theta:		co-lat, radians (ARRAY,IN)
;
;	K : 		max integer index of basis function.
;
;	m : 		order of function (SCALAR,IN)
;
;	dLegendreFns:	latitudinal derivative of the legendre functions
;		returned by the function. Will be an array of size [n_elements(theta_) , K] (OPTIONAL, ARRAY, OUT)
;
;	kVals:		the integer index k values for the functions that are returned (OPTIONAL, ARRAY[K], OUT)
;
;	plotFns: 	set to an integer specifying how many
;		of the legendre functions to plot (SCALAR,IN)
;

function ampere_Legendre_full, theta, K, m, $    ; this is a function of 'm'
	NVALS = nVals, $
	DLEGENDREFNS = dLegendreFns, $
	KVALS	= kVals

  kVals = fltarr(K)
  legendreFns = dblArr( n_elements(theta), K )
  mm=abs(m)
  x_arr=cos(theta)

  eig_cnt=0
  For kk=0,K-1 do begin
    Schmdt=1.d0
    L=kk+1
    If (mm le L) Then begin
      If (m ne 0) Then Schmdt=sqrt(2.d0*factorial(L-mm)/(factorial(L+mm)))
;      If (m ne 0) Then Schmdt=sqrt( (2.0*L+1) * factorial(L-mm)/(factorial(L+mm)))
      lgpf=(-1.)^mm*Schmdt*LEGENDRE(x_arr,L,mm,/double)
      legendreFns[*,eig_cnt]=lgpf
      kVals[eig_cnt] = kk+1
      eig_cnt++
    end
  end
  idx=where(kVals)       ; find non-zero elements
  kVals=kVals[idx]
  nVals=kVals

; Create dP/d_theta of the legendre functions - only take 'K' of the eigenvectors

  dLegendreFns  = dblArr( n_elements(theta), K )
; Loop in K (or L)
  eig_cnt=0
  for kk = 0, K - 1 do begin
    L=kk+1
;stop
    if (mm le L) then begin
      if (m eq 0) then begin
        mlt=sqrt(0.5d0*L*(L+1))
        Schmdt=sqrt(2.d0*factorial(L-1)/(factorial(L+1)))    ;do for M=1
        trm=Schmdt*LEGENDRE(x_arr,L,1,/double)
        dLegendreFns[*,eig_cnt] = mlt*trm
      end else begin
;        If (mm eq 1) Then gm=2.0d0 else gm=1.0d0  ; from old code (does not agree with numerical soln - CLW)
        gm=1.0d0
; Do 1st term; P(L,M-1)
        Schmdt=sqrt(2.0d0*factorial(L-mm+1)/(factorial(L+mm-1)))
        lgp=(-1.)^(mm-1.)*Schmdt*LEGENDRE(x_arr,L,mm-1,/double)
        Ltrm=0.5*sqrt(gm*(L+mm)*(L-mm+1.))*lgp   ;Lhs
; do RH term P(L,M+1)
        If (L eq mm) Then lgp=0.0d0 else begin
          Schmdt=sqrt(2.0d0*factorial(L-mm-1)/(factorial(L+mm+1)))
          lgp=(-1.)^(mm+1)*Schmdt*LEGENDRE(x_arr,L,mm+1,/double)
        end
        Rtrm=-0.5*sqrt((l+mm+1)*(l-mm))*lgp  ;Rhs
        dLegendreFns[*,eig_cnt] = Ltrm+Rtrm
      end
      eig_cnt++
    end         ; if mm le L
  end

	return, legendreFns

end
