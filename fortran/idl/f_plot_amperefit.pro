pro f_plot_amperefit

    iFName = '/Users/dg6/ampere/20100824Amp_invert.ncdf'
    oFName = '/Users/dg6/code/ampereFit/fortran/jParIrid.nc'

    cdfId = ncdf_open ( iFName, /noWrite ) 
	ncdf_varGet, cdfId, 'b_eci', in_b_eci
	
    cdfId = ncdf_open ( oFName, /noWrite ) 
	ncdf_varGet, cdfId, 'b_eci', out_b_eci

stop
	

end
