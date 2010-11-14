pro read_ampere_ncdf,ncdfname,x_axis_frac_hour,pseudosvnum_total,plane_number_total,pos_eci_total,b_eci,pseudo_sv_quality,data_splice
;
; Haje Korth
; JHU/APL
; June 2010
;
; error handling.
	catch, theerror
	if theerror ne 0 then begin
		catch, /cancel
		obj_destroy, fileobj
		return
	endif

; initial definition of objects.
	fileobj = obj_new()

; open the source file in read-only mode.
	fileobj = obj_new('ncdf_file', ncdfname)
	if obj_valid(fileobj) eq 0 then message, 'Invalid file object returned from ncdf_file init.'

; load variables
	x_axis_frac_hour = fileobj -> getvardata('time')
	pseudosvnum_total = fileobj -> getvardata('pseudo_sv_num')
	plane_number_total = fileobj -> getvardata('plane_num')
	pos_eci_total = fileobj -> getvardata('pos_eci')
	b_eci = fileobj -> getvardata('b_eci')
	pseudo_sv_quality = fileobj -> getvardata('pseudo_sv_quality')
	data_splice = fileobj -> getvardata('data_splice')

; destroy the file object.
	obj_destroy, fileobj
	return
end