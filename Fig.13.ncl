load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	wks  = gsn_open_wks("eps","Fig.12")              ; send graphics to PNG file
    
	res                           = True              ; plot mods desired
  
  res@mpCenterLonF = 180
  res@mpMaxLatF    = 50.
  res@mpMinLatF    = -20
  res@mpMaxLonF    = 150
  res@mpMinLonF    = 10
  
  res@tmYRMode             = "Automatic"
  res@tmXTOn               = False
  res@tmYROn               = False
  
  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02
  
  plot = gsn_csm_map_ce(wks,res)                    ; draw global map
  
end