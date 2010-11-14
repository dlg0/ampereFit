;	ampere_realLegendre.pro
;
;	Function to calculate real legendre associated polynomials
;	of real degree (n) and integer order (m). Method is the
;	direct solution of Laplace's equation in spherical coordinates
;	for appropriate boundary conditions. i.e., a boundary value
;	problem.
;
;	Copyright 2008, 2009 Murray D. Sciffer and David L. Green
;   Value_added: Colin Waters -> 2009, May 2010 (added normalisation)
;
;	This program is free software: you can redistribute it and/or modify
;	it under the terms of the GNU General Public License as published by
;	the Free Software Foundation, either version 3 of the License, or
;	(at your option) any later version.
;	This program is distributed in the hope that it will be useful,
;	but WITHOUT ANY WARRANTY; without even the implied warranty of
;	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;	GNU General Public License for more details.
;
;	You should have received a copy of the GNU General Public License
;	along with this program.  If not, see <http://www.gnu.org/licenses/>.
;
;	Variables:
;	---------
;
;	returnValue: a number of legendre functions. The array will be
;		of size n_elements(theta_) x approx. K/2 depending on the
;		selected boundary conditions.
;
;	theta_:		co-lat, degrees (ARRAY,IN)
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
;	lowEndBC:	0 for P=0 and 1 for dPdTh = 0 at poleward boundary. There is no
;		switch for both as yet. Instead, just call the function
;		twice and sort according to the returned kVals.
;		(OPTIONAL, SCALAR, IN)
;
;	highEndBC:	0 for P=0 and 1 for dPdTh = 0 at equatorward boundary. This is the same description as the lowEndBC
;		(OPTIONAL, SCALAR, IN)
;
;	minTheta:	set to zero for a cap
;		(OPTIONAL, SCALAR, IN)
;
;	maxTheta:	max co_Lat
;
;	plotFns: 	set to an integer specifying how many
;		of the legendre functions to plot (SCALAR,IN)
;
;	Example call:
;	------------
;
;	P = ampere_reallegendre (fIndGen(30),10,1,plot=5,lowEndBC=1,kVals=k,dLegendreFns=dP)
;
;	or some high order functions (n>150)
;
;	P = ampere_reallegendre (fIndGen(30),50,0,plot=20,lowEndBC=1,kVals=k,dLegendreFns=dP)

function ampere_realLegendre, theta_, K, m, $    ; this is a function of 'm'
	lowEndBC, highEndBC, $                       ; do not make these keywords (0 does not trigger it)
	PLOTFNS = plotFns, $
	NVALS = eigenvalues_, $
	DLEGENDREFNS = dLegendreFns, $
	KVALS	= kVals, $
	MINTHETA = minTheta, $
	MAXTHETA = maxTheta;

    high_end_bc_tmp=highEndBC                    ; CLW added
;	if m eq 0 then highEndBC = 1                 ; pole dY/dthet=0 for m=0

	if not keyword_set ( minTheta ) then minTheta = 0.0
	if not keyword_set ( maxTheta ) then maxTheta = max(theta_)

	nTh	= K * 4                                           ; finer grid compared with the data
	dTh	= (maxTheta - minTheta)/( nTh - 1 )*!dtor
	theta = fIndGen(nTh) * dTh + minTheta * !dtor
	coeffMatrix	= dblArr( nTh-2, nTh-2 )                  ; used for finite diff calc
	Lp_a=dblarr(nTh-2)                                    ; used for eigenvalues

	for j = 0, nTh - 3 do begin
        if j gt 0 and j lt nTh-3 then begin               ; populate finite diff matrix
			coeffMatrix[j-1,j] =  1d0/(dTh^2) - cos(theta[j+1])/sin(theta[j+1])/(2d0*dTh)
			coeffMatrix[j,j]   = -2d0/(dTh^2) - m^2/(sin(theta[j+1])^2)
			coeffMatrix[j+1,j] =  1d0/(dTh^2) + cos(theta[j+1])/sin(theta[j+1])/(2d0*dTh)
        endif

		if j eq 0 then begin                               ; poleward BC  [CLW changed]
			if highEndBC eq 1 then p0 = 1d0/dTh^2 - cos(theta[1])/(2d0*sin(theta[1])*dTh) else p0=0.d0
			coeffMatrix[0,0] = -2d0/dTh^2 - m^2/sin(theta[1])^2 + 4d0/3d0*p0
			coeffMatrix[1,0] =  1d0/dTh^2 + cos(theta[1])/(2d0*sin(theta[1])*dTh) - 1d0/3d0*p0
		endif

		if j eq nTh - 3 then begin                         ; equatorward end
			if lowEndBC eq 1 then pn = 1d0/dTh^2 + cos(theta[j+1])/sin(theta[j+1])/(2d0*dTh) else pn=0.0d0
			coeffMatrix[j-1,j] =  1d0/dTh^2 - cos(theta[j+1])/(2d0*sin(theta[j+1])*dTh) - 1d0/3d0*pn
            coeffMatrix[j,j]   = -2d0/dTh^2 - m^2/sin(theta[j+1])^2 + 4d0/3d0*pn
		endif
	endfor
                                                           ; get eigenvalues and eigenvectors
	eigenValues	= la_eigenproblem(coeffMatrix, eigenVectors = eigenVectors, /double)
	iiOrder	= sort(eigenvalues)
	eigenValues	= eigenValues[iiOrder]
	eigenVectors= eigenVectors[*,iiOrder]

   ;   add end points
    eigenVectorsAll = fltArr(nTh, nTh)
    eigenVectorsAll[1:nTh-2,1:nTh-2] = eigenVectors
    eigenVectors = temporary ( eigenVectorsAll )

    eigenValuesAll  = fltArr(nTh)
    eigenValuesAll[1:nTh-2] = eigenValues
    eigenValues = temporary ( eigenValuesAll )

	;	Here eigenValues (a) are actually a=-n(n+1) so we need
	;	to find the roots of f(n)=n^2+n+a to get back the values
	;	of n we need.

	eigenValues_	= dblArr ( n_elements ( eigenValues ) )
    eigenValues_    =  ( -1d0 + sqrt ( 1d0 - 4d0 * eigenValues ) ) / 2d0
	eigenValues_keep	= where ( eigenValues_ gt m>1 )           ; exclude eigenval=0, keep +ve ones only
	eigenValues_	= eigenValues_[eigenValues_keep]
	eigenVectors	= eigenVectors[*,eigenValues_keep]

;	Reconstruct the end points if m=0 or lowEndBC = 1

	nVectors	= n_elements ( eigenVectors[0,*] );
	if lowEndBC then begin
    	for i = 0, nVectors - 1 do $
         eigenVectors[nTh-1,i] = 4d0/3d0*eigenVectors[nTh-2,i] - 1d0/3d0*eigenVectors[nTh-3,i]
	endif else eigenVectors[nTh-1,0:nVectors-1]=0.0d0

	if highEndBC eq 1 then begin
		for i = 0, nVectors - 1 do $
         eigenVectors[0,i] = 4d0/3d0*eigenVectors[1,i] - 1d0/3d0*eigenVectors[2,i]
	endif else eigenVectors[0,0:nVectors-1]=0.0d0

;; Normalise Legendre functions (Trapezoidal rule)
;	ortho_n = dblarr(nVectors,nVectors)
;	For ll = 0,nVectors-1 do begin
;		For ll1 = 0,nVectors-1 do begin
;		 End1 = eigenVectors(0,ll) * eigenVectors(0,ll1) * sin(theta(0))
;		 End2 = eigenVectors(nTh-1,ll) * eigenVectors(nTh-1,ll1) * sin(theta(nTh-1))
;		 ortho_n(ll1,ll) = ( 2.0 * (reform(eigenVectors(*,ll) * sqrt(sin(theta)))) ## $
;		                   Transpose(Reform((eigenVectors(*,ll1) * sqrt(sin(theta)))))-End1-End2)*(dTh/2.0d0)
;	   end
; 	end

;; Normalise for azimuth as well for spherical harmonic (with !dpi in azimuth)
;	For ll=0,nVectors-1 do eigenVectors(*,ll) = eigenVectors(*,ll)/sqrt(!dpi*(ortho_n(ll,ll))) ; Basis functions are now an orthonormal basis set

;	Create dY/d_theta of the legendre functions

	dEigenVectors = dblArr(n_elements(theta), K)
	for i = 0, K - 1 do begin
		dEigenVectors[*,i] = deriv(theta, eigenVectors[*,i])
	endfor

;	Interpolate to requested theta_ from equi-spaced theta.

	legendreFns	= dblArr ( n_elements ( theta_ ), K )
	dLegendreFns	= dblArr ( n_elements ( theta_ ), K )

	for i = 0, K - 1 do begin
		legendreFns[*,i]	= interpolate ( eigenVectors[*,i], $
			( theta_ - minTheta ) / ( maxTheta - minTheta ) * ( nTh - 1 ) )
		dLegendreFns[*,i]	= interpolate ( dEigenVectors[*,i], $
			( theta_ - minTheta ) / ( maxTheta - minTheta ) * ( nTh - 1 ) )
	endfor

	if keyword_set ( plotFns ) then begin
		if plotFns gt K then plotFns = K
		window, 0, ySize = 700
		device, decomposed = 0
		!p.background = 255
		loadct, 12
		!p.multi = [0,1,2]

		for i = 0, plotFns-1 do begin
			if i eq 0 then $
				plot, theta * !radeg, eigenVectors[*,i],$
				color = (i+1) * 16 - 1, thick = 1, $
				yRange = [-0.2,0.2], yStyle = 1, $
				xRange = [ minTheta, maxTheta ], xStyle = 1
			if i ne 0 then $
				oPlot, theta * !radeg, eigenVectors[*,i],$
				color = (i+1) * 16 - 1, thick = 1
			oPlot, theta_, legendreFns[*,i],$
				color = (i+1) * 16 - 1, thick = 2.0, $
                psym = 4
		endfor

		for i = 0, plotFns-1 do begin
			if i eq 0 then $
				plot, theta * !radeg, dEigenVectors[*,i],$
				color = (i+1) * 16 - 1, thick = 1, $
				yRange = [-10,10], yStyle = 1, $
				xRange = [ minTheta, maxTheta ], xStyle = 1;
			if i ne 0 then $
				oPlot, theta * !radeg, dEigenVectors[*,i],$
				color = (i+1) * 16 - 1, thick = 1;
			oPlot, theta_, dlegendreFns[*,i],$
				color = (i+1) * 16 - 1, thick = 2.0, $
                psym = 4
		endfor

		!p.multi = 0;

	endif

	if m eq 0 then begin
		if lowEndBC then $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + 2 $
		else $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + 1
	endif else begin
		if lowEndBC then $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + abs ( m ) $
		else $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + abs ( m ) + 1
	endelse

	iiGoodK	= where ( kVals le K )
	legendreFns	= legendreFns[*,iiGoodK];
	dLegendreFns	= dLegendreFns[*,iiGoodK];
	eigenValues_	= eigenValues_[iiGoodK];
	kVals	= kVals[iiGoodK];
    highEndBC=high_end_bc_tmp              ; set back to original value

	return, legendreFns

end
