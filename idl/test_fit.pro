pro test_fit, fileName, sHr, eHr

;	Read data

	fileName	= '20080105_a_RevA.dat'
	header	= 2
	ndat	= file_lines ( fileName ) - header

	data_struct	= { 	ipln:0,$
						sv:0,$
						sday:0d0,$
						rad:0d0,$
						lat:0d0,$
						lon:0d0,$
						dbx:0d0,$
						dby:0d0,$
						dbz:0d0,$
						px:0d0,$
						py:0d0,$
						pz:0d0,$
						vx:0d0,$
						vy:0d0,$
						vz:0d0}

	data	= replicate ( data_struct, ndat )
	openr, rUnit, fileName, /get_lun

	dateStr=''
	readf, rUnit, dateStr
	reads, strmid(dateStr,10,8), $
			year, month, day, $
			format='(i4,i2,i2)'

	skip_lun, rUnit, 1, /lines
	readf, rUnit, data
	free_lun, rUnit

;	Select out a time slice

	sHr = 13 
	eHr	= 16 

	iiTime	= where ( data.sday gt min(data.sday)+60.0*(sHr-8)*60.0 $
					and data.sday lt min(data.sday)+60.0*(eHr-8)*60.0 $
					and data.pz gt 0 )

;	Convert to GEI coords

	geopack_sphcar, data.px, data.py, data.pz,$
			posR, posTh, posPh, /to_sphere, /degree

	cdf_epoch, epoch0, year, month, day, /compute_epoch
	epoch	= data.sday * 1d3 + epoch0
 
	geoPack_reCalc, year, month, day, /date
	
	geoPack_conv_coord, (data.px)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $
			xGEI, yGEI, zGEI, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]

	geoPack_sphCar, xGEI, yGEI, zGEI, $
			R, theta, phi, /to_sphere

;	Create rotation matrix from xyz to xyzGEI

	dx	= 1
	dy	= 1
	dz	= 1

	geoPack_conv_coord, (data.px+dx)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $
			xGEI_x, yGEI_x, zGEI_x, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]
	geoPack_conv_coord, (data.px)[iiTime], (data.py+dy)[iiTime], (data.pz)[iiTime], $
			xGEI_y, yGEI_y, zGEI_y, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]
	geoPack_conv_coord, (data.px)[iiTime], (data.py)[iiTime], (data.pz+dz)[iiTime], $
			xGEI_z, yGEI_z, zGEI_z, $
			/from_geo, /to_gei, $
			epoch = epoch[iiTime]

	rotMat	= fltArr ( 3, 3, n_elements ( iiTime ) )

	rotMat[0,0,*]	= xGEI-xGEI_x
	rotMat[1,0,*]	= xGEI-xGEI_y
	rotMat[2,0,*]	= xGEI-xGEI_z

	rotMat[0,1,*]	= yGEI-yGEI_x
	rotMat[1,1,*]	= yGEI-yGEI_y
	rotMat[2,1,*]	= yGEI-yGEI_z

	rotMat[0,2,*]	= zGEI-zGEI_x
	rotMat[1,2,*]	= zGEI-zGEI_y
	rotMat[2,2,*]	= zGEI-zGEI_z

;	Rotate db vectors to GEI

	dbGEI	= rotMat ## [ $
			[ (data.dbx)[iiTime] ],$
			[ (data.dby)[iiTime] ],$
			[ (data.dbz)[iiTime] ] ] 

;	Rotate from xyzGEI to sphericalGEI

	geoPack_bCarSp, xGEI, yGEI, zGEI, $
			dbGEI[*,0], dbGEI[*,1], dbGEI[*,2], $
			bR, bTheta, bPhi

	loadct, 12
	device, decomposed = 0
	!p.background = 255
	!p.multi = [0,2,2]
	plot, (data.px)[iiTime], (data.py)[iiTime], $
			psym = 4, $
			color = 0

	plot, xGEI, yGEI, $
			psym = 4, $
			color = 0


;	Use density to find the intersection point

	nGrid	= 20
	gridRange	= 20000.0
	gridMin	= -10000.0
	xGrid	= fIndGen ( nGrid ) / ( nGrid - 1 ) * gridRange + gridMin
	gridSize	= abs ( xGrid[0] - xGrid[1] )
	density3D	= fltArr ( nGrid, nGrid, nGrid )

	xIndex	= (xGEI-gridMin)/gridRange*nGrid
	yIndex	= (yGEI-gridMin)/gridRange*nGrid
	zIndex	= (zGEI-gridMin)/gridRange*nGrid

	for i = 0, n_elements ( iiTime ) - 1 do begin
		++density3D[xIndex[i],yIndex[i],zIndex[i]]	
	endfor

	iiMax	= where ( density3D eq max ( density3D ) )
	iiArr	= array_indices ( density3D, iiMax )
	intPt	= xGrid[iiArr] + gridSize / 2

	plots, [intPt[0],intPt[0]], [intPt[1],intPt[1]], $
			psym = 2, $
		   	color = 8 * 16 - 1
	!p.multi = 0
	plot_b_vectors, bPhi, bTheta, $
			theta, phi,$
		   	40, 0, 2, 0, /no_cont, color1=255, thick1=1

	stop
	iPlot, (data.px)[iiTime], (data.py)[iiTime], (data.pz)[iiTime], $
			psym = 4

end
