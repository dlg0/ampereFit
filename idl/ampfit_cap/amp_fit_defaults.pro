
	; Set keyword default values
	; --------------------------

	if keyword_set(aacgmpath) then aacgm_set_path,aacgmpath else aacgm_set_path,path
	; Input data file
	if not keyword_set(savFileName) then savFileName = path+'20100215Amp_invert.sav'
	; Basis functions
	if not keyword_set(kmax) then kmax = 25
	; 0 to 5
	if not keyword_set(mmax) then mmax = 5           
	; Data weighting parameters
	; dBMag threshold, below this we apply 1/sigma: no weighting if thresh=0
	if not keyword_set(thresh) then thresh=80.0
	; sigma=2 would halve the amplitude of each component of the data
	if not keyword_set(sigma) then sigma=1         
	; Controls final solution grid
	; default fit grid
	if not keyword_set(nLatGrid) then nLatGrid	= 50 
	if not keyword_set(nLonGrid) then nLonGrid	= 24
	; Min and Max FAC for plot
	; do not plot abs(FAC) < mn_fac
	if not keyword_set(mn_fac) then mn_fac=0.07      
	if not keyword_set(mx_fac) then mx_fac=0.50
	; Switch to create PNG files of dB and FAC
	; default is to not output PNG files
	if not keyword_set(plt_png) then plt_png=0       
	; Switch to plot fitted dB by Iridium orbit track
	; default is to not to plot track data
	if not keyword_set(plt_tracks) then plt_tracks=0       
	; Plot range in Latitude
	if not keyword_set(aacgm_cap_coLat_deg) then aacgm_cap_coLat_deg = 50.0
	; Shift coord system so pole of fitting GEI system matches data intersction point.
	if not keyword_set(shiftGEI) then shiftGEI = 0
	if not keyword_set(debug) then debug = 0



