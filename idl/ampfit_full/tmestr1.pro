; Form string expression from time
pro tmestr1,Hr,Mn,Sc,hr_str,mn_str,sc_str

 	hr_str=strtrim(string(Hr),2)
 	if Hr lt 10 then hr_str='0'+hr_str
 	mn_str=strtrim(string(Mn),2)
 	if Mn lt 10 then mn_str='0'+mn_str
 	sc_str=strtrim(string(fix(sc)),2)
 	if Sc lt 10 then sc_str='0'+sc_str

end


