; +
; Reads AMPERE data SAV file (input data)
;
; David L Green & Colin L Waters
; Centre for Space Physics
; University of NEwcastle
; Dec 2009
;
; Comments:
; data selection -> latitude range set to max_coLat
; modified may 2010 - sort data by orbit plane and sequentially alomg track (for diagnostics)
;                   - tweak ampdata shift routine
;
; DLG 4-Jul-10
;	There are 2 data strucutres, the un-shifted and shifted.
;	Both are always created and then selected from later.
;


pro read_ampere_sav, sHr, eHr, south, cap_coLat_deg, $
	savFileName = savFileName, $
	dataOriginal = data, $
	dataShifted = dataShifted, $
	year = year, month = month, day = day, $
	avgYrSec = yrSecAvg, $
	avgEpoch = avgEpoch, $
	rot_mat=rot_mat, $
	debug = debug, $
	noShift = noShift

	; IDL sometimes treats this as a variable. 
	; This stops it doing that. 
	; This should not be needed and probably uncovers
   	; a problem with your IDL setup.	

	;forward_function cnvTime                

	if not keyword_set( max_coLat ) then max_coLat = 70.0

	dateStr=''
	reads, strmid(file_baseName ( savFileName ),0,8), year, month, day, format='(i4,i2,i2)'


	; Read in data from save file
	; ---------------------------

	restore, savFileName

	; Variable in save file
	; 
	;   x_axis_frac_hour : UTC in fractional hours
	; plane_number_total : integer number 0-5 for the 6 orbit planes.
	;      pos_ECI_total : ECI position in meters, dimensions [{X,Y,Z},n]  
	;              B_ECI : the same dimension as pos_eci but has delta B in nT


	; Create data structure
	; ---------------------

 	data = { $
		ipln: 0, $
		isat: 0, $
		utc: 0d0, $
		px: 0d0, py: 0d0, pz: 0d0, $
		dbx: 0d0, dby: 0d0, dbz: 0d0, $
		GEI_R_km : 0d0, GEI_coLat_deg: 0d0, GEI_lon_deg: 0d0, $
		GEI_coLat_rad: 0d0, GEI_lon_rad: 0d0, $
		br_GEI: 0d0, bTheta_GEI: 0d0, bPhi_GEI: 0d0 }


	; Select subset of data
	; ---------------------

	nPts = n_elements ( x_axis_frac_hour )

	@constants

	min_z = rIrid_km * cos ( cap_coLat_deg * !dtor )

	if south then begin

  		min_z = -min_z
  		iiSubSet = where ( $
			x_axis_frac_hour ge sHr and $
			x_axis_frac_hour le eHr and $
			pos_eci_total[2,*]*1d-3 lt min_z, iiSubSetCnt)

	endif else begin

 		iiSubSet = where ( $
			x_axis_frac_hour ge sHr and $
			x_axis_frac_hour le eHr and $
			pos_eci_total[2,*]*1d-3 gt min_z, iiSubSetCnt )  

	endelse


	; Fill data structure
	; -------------------

	data = temporary ( $
		replicate ( data, iiSubSetCnt ) )

	data.iPln = plane_number_total[iiSubSet] 
	data.iSat = pseudoSVNum_total[iiSubSet]
	data.utc = x_axis_frac_hour[iiSubSet]

	; GEO XYZ position in KM

	data.px = (pos_ECI_total[0,*])[iiSubSet]*1d-3         
	data.py = (pos_ECI_total[1,*])[iiSubSet]*1d-3
	data.pz = (pos_ECI_total[2,*])[iiSubSet]*1d-3
	
	; GEO XYZ db vector

	data.dbx = (B_ECI[0,*])[iiSubSet]
	data.dby = (B_ECI[1,*])[iiSubSet]
	data.dbz = (B_ECI[2,*])[iiSubSet]


	; Sort data along each orbit track
	; --------------------------------

	sortII = indGen ( iiSubSetCnt )
	cnt = 0

	for trackNo = 0, 5 do begin

		iiThisTrack = where ( data.iPln eq trackNo, thisTrackCnt )
		polarDistance = sqrt ( (data.px)[iiThisTrack]^2 + (data.py)[iiThisTrack]^2 )
		iiStartHere = (where ( polarDistance eq max(polarDistance) ))[0]
		ptDistance = sqrt ( $
				(((data.px)[iiThisTrack])[iiStartHere] $
				- (data.px)[iiThisTrack])^2 $
			+ (((data.py)[iiThisTrack])[iiStartHere] $
				- (data.py)[iiThisTrack])^2 )
		iiThisTrackSorted = sort ( ptDistance )
		sortII[cnt:cnt+thisTrackCnt-1] = iiThisTrack[iiThisTrackSorted]
		cnt += thisTrackCnt

	endfor

	data = data[sortII]


	; Get spherical coords of the GEI XYZ locations
	; ---------------------------------------------

	geopack_sphcar, data.px,data.py,data.pz, $
			GEI_R_km, GEI_coLat_deg, GEI_lon_deg, $
			/to_sphere, /degree   

	; coLat is 0->90 for Nth and 90->180 for Sth (deg)

	data.GEI_R_km = GEI_R_km
	data.GEI_coLat_deg = GEI_coLat_deg
	data.GEI_lon_deg = GEI_lon_deg


	; Rotate XYZ GEI db to spherical GEI 
	; ----------------------------------

	geopack_bcarsp, data.px, data.py, data.pz, $
		data.dbx, data.dby, data.dbz, $
		br_GEI, bTheta_GEI, bPhi_GEI

	data.br_GEI = br_GEI
	data.bTheta_GEI = bTheta_GEI
	data.bPhi_GEI = bPhi_GEI


	; Shift the data to be more centered near the 
	; track intersection point. 
	;
	; You are kidding right! This routine quite unpleasant :(
	; -------------------------------------------------------

	if not keyword_set ( noShift ) then begin

		amp_shiftData, data, dataShifted, rot_mat, south;,/diag

	endif else begin

		dataShifted = data
		rot_mat = [[1.0,0.0,0.0],$
		         [0.0,1.0,0.0],$
		         [0.0,0.0,1.0]]
	endelse


	; Plot data and shifted data
	; --------------------------

	if keyword_set ( debug ) then begin

		device, decomposed = 0
		window, 9, xSize = 600, ySize = 300
		!p.multi = [0,2,1]

 		plt_dat,data.GEI_coLat_deg,data.GEI_lon_deg,$
			n_elements(data.GEI_coLat_deg),-data.bTheta_GEI, data.bPhi_GEI,$
			south,title='Input DataOriginal', capSize=90

		plt_dat,dataShifted.GEI_coLat_deg,dataShifted.GEI_lon_deg,$
			n_elements(dataShifted.GEI_coLat_deg),dataShifted.bTheta_GEI,dataShifted.bPhi_GEI,$
			south,title='Input DataShifted',capSize=90

		!p.multi = 0

	endif


	data.GEI_coLat_rad = data.GEI_coLat_deg * !dtor
	data.GEI_lon_rad = data.GEI_lon_deg * !dtor

	dataShifted.GEI_coLat_rad = dataShifted.GEI_coLat_deg * !dtor
	dataShifted.GEI_lon_rad = dataShifted.GEI_lon_deg * !dtor


	; ********* Check this ************
	; If south then bTheta_GEI_sh = -bTheta_GEI_sh


	; Get the average epoch and times for 
	; GEI to GEOG and AACGM conversion
	; -----------------------------------

	print,'Year, Month, Day : ',year,' ',month,' ',day
	cdf_epoch, epoch0, year, month, day, /compute_epoch
	epoch = data.utc*3600d0*1d3 + epoch0
	avgEpoch = mean(epoch)
	avgHour	= fix ( (sHr + eHr) / 2.0 )
	avgMin	= ( (sHr + eHr) / 2.0 mod 1 ) * 60
	yrSecAvg = cnvTime ( year, month, day, avgHour, avgMin, 0.0 )
	geoPack_reCalc, year, month, day, /date           


end
