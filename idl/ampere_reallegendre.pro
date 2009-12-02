;	ampere_realLegendre.pro
;
;	Function to calculate real legendre associated polynomials
;	of real degree (n) and integer order (m). Method is the 
;	direct solution of Laplace's equation in spherical coordinates
;	for appropriate boundary conditions. i.e., a boundary value
;	problem.
;
;	Copyright 2008, 2009 Murray D. Sciffer and David L. Green 
;
;	This program is free software: you can redistribute it and/or modify
;	it under the terms of the GNU General Public License as published by
;	the Free Software Foundation, either version 3 of the License, or
;	(at your option) any later version.
;	
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
;		This is the k index in Haines, i.e., when considering
;		a cap it will be even for lowEndBC=0 and odd for lowEndBC=1.
;		For a latitudinal slice ( i.e., min(theta_) ne 0 ) i'm not
;		sure of its meaning ... (SCALAR,IN)
;
;	m : 		order of function (SCALAR,IN)
;
;	dLegendreFns:	latitudinal derivative of the legendre functions
;		returned by the function. Will be an array of size
;		n_elements(theta_)
;		(OPTIONAL, ARRAY, OUT)
;
;	kVals:		the integer index k values for the functions that
;		are returned.
;		(OPTIONAL, ARRAY, OUT)
;
;	nVals:		the degree of the selected Legendre functions.
;		This will be an array of length	n_elements(kVals)	
;		(OPTIONAL, ARRAY, OUT)
;
;	lowEndBC:	0 for P=0, 1 for dPdTh = 0. There is no
;		switch for both as yet. Instead, just call the function
;		twice and sort according to the returned kVals. 
;		(OPTIONAL, SCALAR, IN)
;
;	highEndBC:	0 for P=0, 1 for dPdTh = 0. This is the same as
;		the lowEndBC expect for the case of a latitudinal section, 
;		NOT a cap, i.e., min(theta_) gt 0. When using a cap I would
;		not set this. 
;		(OPTIONAL, SCALAR, IN)
;
;	minTheta:	in general the limits of the cap/section are taken
;		as 0 to the limit of the theta_ array. However, to set other
;		limits where the boundary conditions apply set minTheta
;		and/or maxTheta to those alternate values.
;		(OPTIONAL, SCALAR, IN)
;
;	maxTheta:	see minTheta.
;
;	plotFns: 	set to an integer specifying how many
;		of the legendre functions to plot (SCALAR,IN)
;	
;	Example call:
;	------------
;
;	P = ampere_reallegendre (fIndGen(30),10,1,plot=5,lowEndBC=1,nVals=n,kVals=k,dLegendreFns=dP)
;	
;	or some high order functions (n>150)
;
;	P = ampere_reallegendre (fIndGen(30),50,0,plot=20,lowEndBC=1,nVals=n,kVals=k,dLegendreFns=dP)

function ampere_realLegendre, theta_, K, m, $
	LOWENDBC = lowEndBC, $
	HIGHENDBC	= highEndBC, $
	PLOTFNS = plotFns, $
	NVALS = eigenValues_, $
	DLEGENDREFNS = dLegendreFns, $
	KVALS	= kVals, $
	MINTHETA = minTheta, $
	MAXTHETA = maxTheta;

	if NOT keyword_set ( lowEndBC ) then lowEndBC = 0;
	if NOT keyword_set ( highEndBC ) then begin
		if m gt 0 then highEndBC = 0;
		if m eq 0 then highEndBC = 1;
	endif
	if not keyword_set ( minTheta ) then minTheta = 0.0
	if not keyword_set ( maxTheta ) then maxTheta = max ( theta_ )

	nTh	= K * 4;
	dTh	= ( maxTheta - minTheta ) / ( nTh - 1 ) * !dtor;
	theta	= fIndGen ( nTh ) * dTh + minTheta * !dtor;

	coeffMatrix	= dblArr ( nTh-2, nTh-2 );

	for j = 0, nTh - 3 do begin

        if j gt 0 and j lt nTh-3 then begin
	
			term1	= 1d0 / ( dTh^2 ) $
				- cos ( theta[j+1] ) / sin ( theta[j+1] ) / ( 2d0 * dTh );
			
			term2	= 1d0 / ( dTh^2 ) $
				+ cos ( theta[j+1] ) / sin ( theta[j+1] ) / ( 2d0 * dTh );

			term3	= - 2d0 / ( dTh^2 ) $
				- m^2 / ( sin ( theta[j+1] )^2 );

			coeffMatrix[j,j]	= term3;
			coeffMatrix[j-1,j]	= term1;
			coeffMatrix[j+1,j]	= term2;

        endif

		if j eq 0 then begin
			
			if highEndBC eq 1 or m eq 0 then begin

                    p0  = 1 / dTh^2 - cos ( theta[j+1] ) / sin ( theta[j+1] ) / ( 2 * dTh )
					term0	= -2d0 / dTh^2 $
						+ 4.0 / 3.0 * p0 $
						- m^2 / sin ( theta[j+1] )^2;
					term1	= 1d0 / dTh^2 $
						+ cos ( theta[j+1] ) / ( 2d0 * sin ( theta[j+1] ) * dTh ) $
                        - 1.0 / 3.0 * p0;

					coeffMatrix[j,j]	= term0;
					coeffMatrix[j+1,j]	= term1;

			endif else begin

                    p0  = 0
					term0	= -2d0 / dTh^2 $
						+ 4.0 / 3.0 * p0 $
						- m^2 / sin ( theta[j+1] )^2;
					term1	= 1d0 / dTh^2 $
						+ cos ( theta[j+1] ) / ( 2d0 * sin ( theta[j+1] ) * dTh ) $
                        - 1.0 / 3.0 * p0;

					coeffMatrix[j,j]	= term0;
					coeffMatrix[j+1,j]	= term1;

			endelse

		endif

		if j eq nTh - 3 then begin

			if lowEndBC then begin

                pn  = 1 / dTh^2 + cos ( theta[j+1] ) / sin ( theta[j+1] ) / ( 2 * dTh )
                term0   = -2d0 / dTh^2 $
                    + 4.0 / 3.0 * pn $
                    - m^2 / sin ( theta[j+1] )^2

				term1	= 1d0 / dTh^2 $
					-	cos ( theta[j+1] ) / ( 2d0 * sin ( theta[j+1] ) * dTh ) $
                    - 1.0 / 3.0 * pn
	
				coeffMatrix[j,j]	= term0;
				coeffMatrix[j-1,j]	= term1;

			endif else begin

                pn  = 0
	            term0   = -2d0 / dTh^2 $
                    + 4.0 / 3.0 * pn $
                    - m^2 / sin ( theta[j+1] )^2

				term1	= 1d0 / dTh^2 $
					-	cos ( theta[j+1] ) / ( 2d0 * sin ( theta[j+1] ) * dTh ) $
                    - 1.0 / 3.0 * pn
	
				coeffMatrix[j,j]	= term0;
				coeffMatrix[j-1,j]	= term1;
		
			endelse

		endif

	endfor

	eigenValues	= la_eigenproblem ( coeffMatrix, $
		eigenVectors = eigenVectors, /double );
	iiOrder	= sort ( eigenValues );

	eigenValues	= eigenValues[iiOrder];
	eigenVectors	= eigenVectors[*,iiOrder];

    ;   add end points

    eigenVectorsAll    = fltArr ( nTh, nTh )
    eigenVectorsAll[1:nTh-2,1:nTh-2]    = eigenVectors
    eigenVectors    = temporary ( eigenVectorsAll )

    eigenValuesAll  = fltArr ( nTh )
    eigenValuesAll[1:nTh-2] = eigenValues
    eigenValues = temporary ( eigenValuesAll )

	;	Here eigenValues (a) are actually a=-n(n+1) so we need
	;	to find the roots of f(n)=n^2+n+a to get back the values
	;	of n we need. Here I use Newton's method. DLG Update: this 
	;	can easily be replaced by using the quadratic equation
	;	formula! ID10T!

	eigenValues_	= dblArr ( n_elements ( eigenValues ) );

    eigenValues_    =  ( -1 + sqrt ( 1 - 4 * eigenValues ) ) / 2
	eigenValues_keep	= where ( eigenValues_ gt m>1 );
	eigenValues_	= eigenValues_[eigenValues_keep];
	eigenVectors	= eigenVectors[*,eigenValues_keep];

	;	Reconstruct the end points if m=0 or lowEndBC = 1

	nVectors	= n_elements ( eigenVectors[0,*] );
	if lowEndBC then begin
		
		for i = 0, nVectors - 1 do begin

			eigenVectors[nTh-1,i]	= 4d0 / 3d0 * eigenVectors[nTh-2,i] $
				- 1d0 / 3d0 * eigenVectors[nTh-3,i];

		endfor

	endif

	if highEndBC eq 1 then begin

		for i = 0, nVectors - 1 do begin

			eigenVectors[0,i]	= 4d0 / 3d0 * eigenVectors[1,i] $
				- 1d0 / 3d0 * eigenVectors[2,i];

		endfor

	endif

	if highEndBC eq 2 AND m eq 0 then begin

		for i = 0, nVectors - 1 do begin

			eigenVectors[0,i]	= 4d0 / 3d0 * eigenVectors[1,i] $
				- 1d0 / 3d0 * eigenVectors[2,i];

		endfor

	endif


	;	Create the first derivative of the legendre functions

	dEigenVectors	= dblArr ( n_elements ( theta ), K );
	for i = 0, K - 1 do begin

		dEigenVectors[*,i]	= deriv ( theta, eigenVectors[*,i] );

	endfor

	;	Interpolate to requested theta_ from equi-spaced theta.

	legendreFns	= dblArr ( n_elements ( theta_ ), K );
	dLegendreFns	= dblArr ( n_elements ( theta_ ), K );

	for i = 0, K - 1 do begin
		
		legendreFns[*,i]	= interpolate ( eigenVectors[*,i], $
			( theta_ - minTheta ) / ( maxTheta - minTheta ) * ( nTh - 1 ) );
		dLegendreFns[*,i]	= interpolate ( dEigenVectors[*,i], $
			( theta_ - minTheta ) / ( maxTheta - minTheta ) * ( nTh - 1 ) );

	endfor

	if keyword_set ( plotFns ) then begin

		if plotFns gt K then plotFns = K

		window, 0, ySize = 700;
		device, decomposed = 0;
		!p.background = 255;
		loadct, 12;
		!p.multi = [0,1,2];

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
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + 1;
		
	endif else begin

		if lowEndBC then $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + abs ( m ) $
		else $
			kVals	= indGen ( n_elements ( eigenValues_ ) ) * 2 + abs ( m ) + 1;
	
	endelse
	
	iiGoodK	= where ( kVals le K );

	legendreFns	= legendreFns[*,iiGoodK];
	dLegendreFns	= dLegendreFns[*,iiGoodK];
	eigenValues_	= eigenValues_[iiGoodK];
	kVals	= kVals[iiGoodK];

	return, legendreFns;

end
