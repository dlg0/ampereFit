pro ampfit_plotdata

    iFName = 'data/20100824Amp_invert.ncdf'

	fileList = file_search ( 'output/ampData_*nc' )

    cdfId = ncdf_open ( iFName, /noWrite ) 

	    ncdf_varGet, cdfId, 'b_eci', in_b_eci

	ncdf_close, cdfId

	for i=0,n_elements(fileList)-1 do begin

    	cdfId = ncdf_open ( fileList[i], /noWrite ) 

		    ncdf_varGet, cdfId, 'bR', bR 
		    ncdf_varGet, cdfId, 'bT', bT 
		    ncdf_varGet, cdfId, 'bP', bP 

		    ncdf_varGet, cdfId, 'bX', bx 
		    ncdf_varGet, cdfId, 'bY', by 
		    ncdf_varGet, cdfId, 'bZ', bz 

		    ncdf_varGet, cdfId, 'T', T 
		    ncdf_varGet, cdfId, 'P', P
		    ncdf_varGet, cdfId, 'R', R

    	    ncdf_varGet, cdfId, 'X', X
    	    ncdf_varGet, cdfId, 'Y', Y
    	    ncdf_varGet, cdfId, 'Z', Z

		ncdf_close, cdfId

    	iiNorth = where ( T*!radeg gt 0 and Z gt 0 )
    	ampfit_plotvecs, R[iiNorth], T[iiNorth], P[iiNorth], $
    	    bT[iiNorth], bP[iiNorth], title = file_basename(fileList[i],'.nc'), $
			fileNameIn=strJoin(['output/',file_basename(fileList[i],'.nc'),'.eps'])

		print, fileList[i]
		print, 'Max dB_P: ', max(abs(bP[iiNorth]))

	endfor

end
