pro set_map, capSize, title = title

	map_set, 90, 0, 0, $
		/ortho, $
		/iso, $
		limit = [ 90.0 - ( capSize + 2 ), 0, 90, 360 ], $
		/noborder, $
		/advance, $
		title = title, $
		clip = 1
	
	map_grid, label = 1, latDel = 10.0

end
