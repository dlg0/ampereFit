pro ampere_setupSHFns, r, coLat, lon, maxK, maxM, $
	brBFnArr = brBFnArr, $
	bThBFnArr = bThBFnArr, $
	bPhBFnArr = bPhBFnArr, $
	YBFnArr = YBFnArr, $
	MINTHETA = minTheta, $
	MAXTHETA = maxTheta, $
	DLAT = dLat, $
	DLON = dLon, $
	mArr = mArr, $
	nArr = nArr, $
	kArr = kArr, $
    bc1 = bc1, $
    bc2 = bc2

	;r	= 1d0;
	rE	= 1d0;

	;dLat	= abs ( coLat[0] - coLat[1] );
	;dLon	= abs ( lon[0] - lon[1] );

;------
;	Test variables
;	maxM = 2;
;	maxK	= 10;
;	r	= 1.02d0;
;;	coLat	= (rebin ( ( fIndgen ( 181 ) ) * !dtor, 181, 12 ))[*];
;;	lon	= (transpose(rebin(fIndGen ( 12 ) * 2*!pi/12, 12, 181 )))[*];
;	coLat	= fIndgen ( 181 ) * !dtor ;+ 20 * !dtor;
;	lon	= coLat * 0 + !pi / 4;
;	dLon	= 2d0 * !pi / 12.0;
;	dLat	= 1d0 * !dtor;
;	minTheta	= 20d0;
;	maxTheta	= 70d0;
;------

	nPts	= n_elements ( coLat );

	nBFns = 0;
	for k = 1, maxK do begin
		nBFns	= nBFns + (k<maxM) * 2;
	endfor
	nBFns	= nBFns + maxK;

	nBFns	= nBFns * 2;

	mArr	= intArr ( nBFns );
	nArr	= dblArr ( nBFns );
	kArr	= intArr ( nBFns );

	YBFnArr	= dblArr ( nBFns, nPts );
	brBFnArr	= dblArr ( nBFns, nPts );
	bThBFnArr	= dblArr ( nBFns, nPts );
	bPhBFnArr	= dblArr ( nBFns, nPts );

	bFnCnt	= 0;
	for m = -maxM, maxM do begin
		print, 'M: ',m;
	print, 'Calling realLegendre ...';
		legendreFnsBC1	= ampere_realLegendre ( coLat * !radeg, maxK, abs(m), $
			lowEndBC = 0, plotFns = 0, nVals = nValsBC1, $
			dLegendreFns = dLegendreFnsBC1, $
			kVals = kValsBC1, highEndBC = 0 , $
			minTheta = minTheta, maxTheta = maxTheta );

		legendreFnsBC2	= ampere_realLegendre ( coLat * !radeg, maxK, abs(m), $
			lowEndBC = 1, plotFns = 0, nVals = nValsBC2, $
			dLegendreFns = dLegendreFnsBC2, $
			kVals = kValsBC2, highEndBC = 0 , $
			minTheta = minTheta, maxTheta = maxTheta );

        ;kValsBC1    = kValsBC1 - 1

        if keyword_set ( bc1 ) then begin

		    legendreFns	= legendreFnsBC1;
		    dLegendreFns	= dLegendreFnsBC1;

		    nVals	= nValsBC1;
		    kVals	= kValsBC1;

        endif else if keyword_set ( bc2 ) then begin

		    legendreFns	= legendreFnsBC2;
		    dLegendreFns	= dLegendreFnsBC2;

		    nVals	= nValsBC2;
		    kVals	= kValsBC2;

        endif else begin

            print, 'both sets'
		    legendreFns	= [ [ legendreFnsBC2 ], [ legendreFnsBC1 ] ];
		    dLegendreFns	= [ [ dLegendreFnsBC2 ], [ dLegendreFnsBC1 ] ];

		    nVals	= [ nValsBC2, nValsBC1 ];
		    kVals	= [	kValsBC2, kValsBC1 ];

        endelse

		print, '+++++++++++++++++++++++'
		print, n_elements ( kVals )
		print, maxK - abs ( m)

		for i = 0, n_elements ( kVals ) - 1 do begin

			rDep	= r^nVals[i] * $
				( 1d0 + nVals[i] / ( nVals[i] + 1d0 ) * ( rE / r )^( 2d0 * nVals[i] + 1d0 ) );
			dYdr_rDep	= nVals[i] * r^( nVals[i] - 1d0 ) * $
				( 1d0 - ( rE / r )^( 2d0 * nVals[i] + 1d0 ) );

			;print, rDep, dYdr_rDep;

			if m ge 0 then begin

				YBFnArr[bFnCnt,*]	= rDep * legendreFns[*,i] * cos ( m * lon );
				brBFnArr[bFnCnt,*]	= dYdr_rDep * legendreFns[*,i] * cos ( m * lon );
				bThBFnArr[bFnCnt,*]	= 1d0 / r * rDep * dLegendreFns[*,i] * cos ( m * lon );
				bPhBFnArr[bFnCnt,*]	= -m * rDep / ( r * sin ( coLat ) ) * $
					legendreFns[*,i] * sin ( m * lon );

				;orthoNormFactor	=  sqrt ( 1d0 / ( total ( YBFnArr[bFnCnt,*]^2 * $
				;	dLat * dLon * sin ( coLat ) ) ) );

				;YBFnArr[bFnCnt,*]	= orthoNormFactor * YBFnArr[bFnCnt,*];
				;brBFnArr[bFnCnt,*]	= orthoNormFactor * brBFnArr[bFnCnt,*];
				;bThBFnArr[bFnCnt,*]	= orthoNormFactor * bThBFnArr[bFnCnt,*];
				;bPhBFnArr[bFnCnt,*]	= orthoNormFactor * bPhBFnArr[bFnCnt,*];

				;print, 'Inner Prod: ', sqrt ( 1d0 / ( total ( legendreFns[*,i]^2 * $
				;	dLat * dLon * sin ( coLat ) ) ) );
				;print, sqrt ( ( 2d0 * nVals[i] + 1d0 ) / ( 4d0 * !dpi ) * $
				;	gamma ( nVals[i] - m + 1d0 ) / gamma ( nVals[i] + m + 1d0 ) );

				mArr[bFnCnt]	= m;
				nArr[bFnCnt]	= nVals[i];
				kArr[bFnCnt]	= kVals[i];

				++bFnCnt;

			endif else begin

		;		YBFnArr[bFnCnt,*]	= rDep * legendreFns[*,i] * cos ( m * lon );
		;		brBFnArr[bFnCnt,*]	= dYdr_rDep * legendreFns[*,i] * cos ( m * lon );
		;		bThBFnArr[bFnCnt,*]	= 1d0 / r * rDep * dLegendreFns[*,i] * cos ( m * lon );
		;		bPhBFnArr[bFnCnt,*]	= -m * rDep / ( r * sin ( coLat ) ) * $
		;			legendreFns[*,i] * sin ( m * lon );
		;
		;		mArr[bFnCnt]	= m;
		;		nArr[bFnCnt]	= nVals[i];
		;		kArr[bFnCnt]	= kVals[i];
		;		++bFnCnt;

				YBFnArr[bFnCnt,*]	= rDep * legendreFns[*,i] * sin ( abs(m) * lon );
				brBFnArr[bFnCnt,*]	= dYdr_rDep * legendreFns[*,i] * sin ( abs(m) * lon );
				bThBFnArr[bFnCnt,*]	= 1d0 / r * rDep * dLegendreFns[*,i] * sin ( abs(m) * lon );
				bPhBFnArr[bFnCnt,*]	= abs(m) * rDep / ( r * sin ( coLat ) ) * $
					legendreFns[*,i] * cos ( abs(m) * lon );

				;orthoNormFactor	=  sqrt ( 1d0 / ( total ( legendreFns[*,i]^2 * $
				;	dLat * dLon * sin ( coLat ) ) ) );

				;YBFnArr[bFnCnt,*]	= orthoNormFactor * YBFnArr[bFnCnt,*];
				;brBFnArr[bFnCnt,*]	= orthoNormFactor * brBFnArr[bFnCnt,*];
				;bThBFnArr[bFnCnt,*]	= orthoNormFactor * bThBFnArr[bFnCnt,*];
				;bPhBFnArr[bFnCnt,*]	= orthoNormFactor * bPhBFnArr[bFnCnt,*];

				mArr[bFnCnt]	= m;
				nArr[bFnCnt]	= nVals[i];
				kArr[bFnCnt]	= kVals[i];
				++bFnCnt;

			endelse

		endfor

	endfor

;	Remove zero functions if only one set chosen

	YBFnArr	= YBFnArr[0:bFnCnt-1,*];
	brBFnArr	= brBFnArr[0:bFnCnt-1,*];
	bThBFnArr	= bThBFnArr[0:bFnCnt-1,*];
	bPhBFnArr	= -bPhBFnArr[0:bFnCnt-1,*];
	nArr	= nArr[0:bFnCnt-1];
	mArr	= mArr[0:bFnCnt-1];
	kArr	= kArr[0:bFnCnt-1];

end
