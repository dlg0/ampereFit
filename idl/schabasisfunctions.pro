pro schaBasisFunctions, kMax, mMax, capSize, inCoLats, inLons, $
	EVENSET=evenSet, ODDSET=oddSet,$
	YkmBFns=YkmBFns, dYkmDthBFns=dYkmDthBFns, dYkmDphBFns=dYkmDphBFns,$
	OUTNKVALUES=outNkValues, OUTKVALUES=outKValues, OUTMVALUES=outMValues, $
	PNMPATH = pnmPath

	;	Restore the pnm for specific capSize

	restore, pnmPath+strcompress(string(fix(capSize),format='(i2.2)'),/rem)+'.sav'
	coLat	= coLat[1:*]*!dtor	;	There is a -Nan @ theta=0 for the derivative so discard
														;	this point for use in the interpolation

	;	Select appropriate functions

	if keyword_set(evenSet) then begin	;	Function == 0 @ capSize
		iiKeep = where((fix(kValues)-fix(abs(mValues))) MOD 2 EQ 0 and kValues le kMax and abs(mValues) le mMax)
		outNkValues	= nkValues[iiKeep]
		outMValues	= mValues[iiKeep]
		outKValues	= kValues[iiKeep]
		Pnm	= Pnm[iiKeep,1:*]
		dPnm	= dPnm[iiKeep,1:*]
	endif else if keyword_set(oddSet) then begin	;	Derivative == 0 @ capSize
		iiKeep = where((fix(kValues)-fix(abs(mValues))) MOD 2 EQ 1 and kValues le kMax and abs(mValues) le mMax)
		outNkValues	= nkValues[iiKeep]
		outMValues	= mValues[iiKeep]
		outKValues	= kValues[iiKeep]
		Pnm	= Pnm[iiKeep,1:*]
		dPnm	= dPnm[iiKeep,1:*]
	endif else begin
		iiKeep = where(kValues le kMax and abs(mValues) le mMax)
		outNkValues	= nkValues[iiKeep]
		outMValues	= mValues[iiKeep]
		outKValues	= kValues[iiKeep]
		Pnm	= Pnm[iiKeep,1:*]
		dPnm	= dPnm[iiKeep,1:*]
	endelse

	;	Interpolate to desired inCoLats and create Ykm values

	PnmAtInCoLats	= fltarr(n_elements(Pnm[*,0]),n_elements(inCoLats[*]))
	dPnmAtInCoLats	= fltarr(n_elements(Pnm[*,0]),n_elements(inCoLats[*]))
	YkmBFns	= fltarr(n_elements(Pnm[*,0]),n_elements(inCoLats[*]))
	dYkmDthBFns	= fltarr(n_elements(Pnm[*,0]),n_elements(inCoLats[*]))
	dYkmDphBFns	= fltarr(n_elements(Pnm[*,0]),n_elements(inCoLats[*]))

	for jj=0,n_elements(Pnm[*,0])-1 do begin
		PnmAtInCoLats[jj,*]	= interpol((Pnm[jj,*])[*],coLat[*],inCoLats[*])
		dPnmAtInCoLats[jj,*]	= interpol((dPnm[jj,*])[*],coLat[*],inCoLats[*])
		Knm = sqrt(2.)*sqrt(factorial(outNkValues[jj]-abs(outMValues[jj])) $
			/factorial(outNkValues[jj]+abs(outMValues[jj])))
		if outMValues[jj] eq 0 then Knm=1
		if outMValues[jj] ge 0 then begin
			YkmBFns[jj,*] = Knm*cos(abs(outMValues[jj])*inLons[*])*PnmAtInCoLats[jj,*]
			dYkmDphBFns[jj,*]	= Knm*abs(outMValues[jj])/sin(inCoLats[*]) $
				*sin(abs(outMValues[jj])*inLons[*])*PnmAtInColats[jj,*]
			dYkmDthBFns[jj,*]	= Knm*cos(abs(outMValues[jj])*inLons[*]) $
				*dPnmAtInCoLats[jj,*]
		endif else begin
			YkmBFns[jj,*] = Knm*sin(abs(outMValues[jj])*inLons[*])*PnmAtInCoLats[jj,*]
			dYkmDphBFns[jj,*]	= -Knm*abs(outMValues[jj])/sin(inCoLats[*]) $
				*cos(abs(outMValues[jj])*inLons[*])*PnmAtInColats[jj,*]
			dYkmDthBFns[jj,*]	= Knm*sin(abs(outMValues[jj])*inLons[*]) $
				*dPnmAtInCoLats[jj,*]
		endelse

	endfor

	;	Create Ynm from Pnm

	;Knm = sqrt(2.)*sqrt(factorial(nkValues-abs(outMValues))/factorial(nkValues+abs(outMValues)))

;	iiNeg	= where(outMValues lt 0)
;	iiPos	= iiNeg EQ 0
;	YkmBFns	= fltarr(n_elements(PnmAtInCoLats[*,0]),n_elements(PnmAtInCoLats[0,*]))
;	YkmBFns[iiPos,*]	= (cos(abs(outMValues[iiPos])*inLons))[*,intarr(n_elements(PnmAtInCoLats[0,*]))]*PnmAtInCoLats[iiPos,*]
;	YkmBFns[iiNeg,*]	= (sin(abs(outMValues[iiPos])*inLons))[*,intarr(n_elements(PnmAtInCoLats[0,*]))]*PnmAtInCoLats[iiPos,*]

;	loadct, 12
;	device, decomposed=0
;	!p.multi=[0,1,2]
;
;	for jj=0,n_elements(Pnm[*,0])-1 do begin
;		plot,	coLat*!radeg, Pnm[jj,*]
;		oplot, [0,100],[0,0]
;		oplot, [capSize,capSize],[-1e7,1e7]
;		oplot, inCoLats*!radeg, PnmAtInCoLats[jj,*], psym=4, color=9*16-1
;		plot,	coLat*!radeg, dPnm[jj,*]
;		oplot, [0,100],[0,0]
;		oplot, [capSize,capSize],[-1e7,1e7]
;		oplot, inCoLats*!radeg, dPnmAtInCoLats[jj,*], psym=4, color=9*16-1
;
;
;	endfor
;	!p.multi=0
end
