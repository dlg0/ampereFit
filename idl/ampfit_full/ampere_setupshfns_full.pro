pro ampere_setupSHFns_full, r, colat, lon, $
  maxK, maxM, $
	brBFnArr = brBFnArr, $
	bThBFnArr = bThBFnArr, $
	bPhBFnArr = bPhBFnArr, $
	YBFnArr = YBFnArr, $
	mArr = mArr, $
	nArr = nArr, $
	kArr = kArr

	rE	= 1d0;
	nPts	= n_elements ( coLat )      ; number of dB data points
	nBFns = 0                         ; recoded to only have single basis set [CLW sept 2010]
  For i=0,maxM do nBFns=nBFns + (2*i+1)
  n_left=maxK-maxM
  nBFns=nBFns + n_left*(2*maxM+1)-1   ; Num of basis functions given maxK and maxM excluding zero term

	mArr	= intArr ( nBFns )
	nArr	= dblArr ( nBFns )
	kArr	= intArr ( nBFns )

	YBFnArr	= dblArr ( nBFns, nPts )
	brBFnArr	= dblArr ( nBFns, nPts )
	bThBFnArr	= dblArr ( nBFns, nPts )
	bPhBFnArr	= dblArr ( nBFns, nPts )

	bFnCnt	= 0;
	for m = 0, maxM do begin        ; loop in 'm'

; returns basis function values at requested coLats - interpolate is done in ampere_realLegendre
		legendreFns	= ampere_Legendre_full ( coLat, maxK, m, $
      nVals = nVals, $
      dLegendreFns = dLegendreFns, $
      kVals = kVals)

		for i = 0, n_elements ( kVals ) - 1 do begin
			rDep	= r^nVals[i] * $
				( 1d0 + nVals[i] / ( nVals[i] + 1d0 ) * ( rE / r )^( 2d0 * nVals[i] + 1d0 ) );

;			dYdr_rDep	= nVals[i] * r^( nVals[i] - 1d0 ) * $
;				( 1d0 - ( rE / r )^( 2d0 * nVals[i] + 1d0 ) );
      dYdr_rDep = 1.0d0

; calc for +m basis set
			YBFnArr[bFnCnt,*]	= rDep * legendreFns[*,i] * cos ( m * lon )
			brBFnArr[bFnCnt,*]	= dYdr_rDep * legendreFns[*,i] * cos ( m * lon )
			bThBFnArr[bFnCnt,*]	= 1d0 / r * rDep * dLegendreFns[*,i] * cos ( m * lon )
			bPhBFnArr[bFnCnt,*]	= -m * rDep / ( r * sin ( coLat ) ) * legendreFns[*,i] * sin ( m * lon )

			mArr[bFnCnt]	= m
			nArr[bFnCnt]	= nVals[i]
			kArr[bFnCnt]	= kVals[i]
			++bFnCnt

; calc for -m basis set
      If (m gt 0 ) then begin
        YBFnArr[bFnCnt,*]	= rDep * legendreFns[*,i] * sin ( m * lon )
        brBFnArr[bFnCnt,*]	= dYdr_rDep * legendreFns[*,i] * sin ( m * lon )
        bThBFnArr[bFnCnt,*]	= 1d0/r * rDep * dLegendreFns[*,i] * sin ( m * lon )
        bPhBFnArr[bFnCnt,*]	= m*rDep/(r*sin(coLat))*legendreFns[*,i]*cos( m * lon)
        mArr[bFnCnt] = -m
        nArr[bFnCnt]	= nVals[i]
        kArr[bFnCnt]	= kVals[i]
        ++bFnCnt
      end
		endfor
	endfor        ; m loop

end
