
	; Set keyword default values
	; --------------------------

	; Input data file
	if not keyword_set(datFileName) then datFileName = path+'20100602Amp_invert.ncdf'
	; Basis functions
	if not keyword_set(kmax) then kmax = 30
	; 0 to 5
	if not keyword_set(mmax) then mmax = 4
;
	if not keyword_set(nLonGrid) then nLonGrid	= 24
	; Min and Max FAC for plot
	; do not plot abs(FAC) < mn_fac
	if not keyword_set(mn_fac) then mn_fac=0.07
	if not keyword_set(mx_fac) then mx_fac=0.50
	; Switch to create PNG files of dB and FAC
	; default is to not output PNG files
	if not keyword_set(plt_png) then plt_png=0
	; Shift coord system so pole of fitting GEI system matches data intersction point.
	if not keyword_set(debug) then debug = 0
	if not keyword_set(extend) then extend = 0
