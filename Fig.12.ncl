load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	;===read data
	
	file_EXP  = addfile("~/model/LBM/linear.v1.nc","r")
	
	hgt = file_EXP->z(:,11,:,:)
	printVarSummary(hgt)
	hgt_anom = dim_avg_n_Wrap(hgt(16:64,:,:),0)
	printVarSummary(hgt_anom)
	
	file_CTL  = addfile("~/model/LBM/frc.nc","r")
	
	frc = file_CTL->t(0,6,:,:)
	printVarSummary(frc)
	
	frc_rate = frc*86400
	copy_VarMeta(frc,frc_rate)
  printVarSummary(frc_rate)
	;======
	print("...read in data files...")
	
  ;===
  f_hgt  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc", "r")
  
	TIME    = f_hgt->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
  ymdStrt_Mar = 19790301                         ; start yyyymmdd
  ymdLast_Mar = 20190301                          
  ymdStrt_Apr = 19790401                         ; start yyyymmdd
  ymdLast_Apr = 20190401                         ; last  
  ymdStrt_May = 19790501                         ; start yyyymmdd
  ymdLast_May = 20190501                         ; last  
  
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20190601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20190701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20190801                         ; last  
  
  iStrt_Mar   = ind(YMD.eq.ymdStrt_Mar)              ; index start
  iLast_Mar   = ind(YMD.eq.ymdLast_Mar)              ; index last 
  iStrt_Apr   = ind(YMD.eq.ymdStrt_Apr)              ; index start
  iLast_Apr   = ind(YMD.eq.ymdLast_Apr)              ; index last 
  iStrt_May   = ind(YMD.eq.ymdStrt_May)              ; index start
  iLast_May   = ind(YMD.eq.ymdLast_May)              ; index last 
  
  iStrt_Jun   = ind(YMD.eq.ymdStrt_Jun)              ; index start
  iLast_Jun   = ind(YMD.eq.ymdLast_Jun)              ; index last 
  iStrt_Jul   = ind(YMD.eq.ymdStrt_Jul)              ; index start
  iLast_Jul   = ind(YMD.eq.ymdLast_Jul)              ; index last 
  iStrt_Aug   = ind(YMD.eq.ymdStrt_Aug)              ; index start
  iLast_Aug   = ind(YMD.eq.ymdLast_Aug)              ; index last 
  
  hgt_Mar   = short2flt(f_hgt->hgt(iStrt_Mar:iLast_Mar:12,{200},::-1,:)) 
  hgt_Apr   = short2flt(f_hgt->hgt(iStrt_Apr:iLast_Apr:12,{200},::-1,:)) 
  hgt_May   = short2flt(f_hgt->hgt(iStrt_May:iLast_May:12,{200},::-1,:)) 
  
  hgt_Jun   = f_hgt->hgt(iStrt_Jun:iLast_Jun:12,{200},::-1,:) 
  hgt_Jul   = f_hgt->hgt(iStrt_Jul:iLast_Jul:12,{200},::-1,:) 
  hgt_Aug   = f_hgt->hgt(iStrt_Aug:iLast_Aug:12,{200},::-1,:) 
  
  hgt_JJA   = (hgt_Jun+hgt_Jul+hgt_Aug)/3.
  copy_VarMeta(hgt_Mar,hgt_JJA)
  
  hgt_JJA_clim = dim_avg_n_Wrap(hgt_JJA,0)
	;===
	
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.12")
  gsn_define_colormap(wks,"MPL_bwr")
  
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@gsnSpreadColors = True
  
  res@gsnContourZeroLineThicknessF = 2 ;
  res@gsnContourNegLineDashPattern = 2 ; 
  
  res@cnLinesOn        = True
  res@cnMonoLineColor  = True
  res@cnLineColor      = 0
  res@cnInfoLabelOn    = False
  res@cnLineLabelsOn   = False
  res@cnFillOn         = True
  res@cnLineThicknessF = 1.3
  res@cnLineLabelFont  = "helvetica"
  
  res@tiMainString       = " "
  
  res@cnLineLabelFontHeightF   = 0.01
  res@cnLineLabelPlacementMode = "Computed"
  res@cnLineLabelInterval      = 1
  res@cnLineLabelFormat        = "@*+^sg"
  
  res@tmYRMode             = "Automatic"
  res@tmXTOn               = False
  res@tmYROn               = False
  
  res@mpShapeMode          = "FreeAspect"
  res@vpWidthF             = 0.7
  res@vpHeightF            = 0.35
  
  res@mpFillOn                    = True
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  res@mpCenterLonF = 180
  res@mpMaxLatF    = 60.
  res@mpMinLatF    = -20
  res@mpMaxLonF    = 180
  res@mpMinLonF    = 10
  
  res@tmXBLabelFontHeightF = 0.021
  res@tmYLLabelFontHeightF = 0.021
  
  res@gsnLeftStringFontHeightF   = 0.024
  res@gsnRightStringFontHeightF  = 0.024
  res@gsnCenterStringFontHeightF = 0.021
  
  res@gsnCenterStringOrthogonalPosF = 0.015
  res@gsnRightStringOrthogonalPosF  = 0.015
  res@gsnLeftStringOrthogonalPosF   = 0.015
  
  res@tmXBOn               = True
  res@tmYLOn               = True
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(b)"
  res@gsnRightString              = "16-65 days mean (Z200)"
  
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/-4,-3,-2,-1,1,2,3,4,5,6/)
  
  plot = new(2,graphic)
  
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.16
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.02
  
  plot(1)  = gsn_csm_contour_map_ce(wks,dim_avg_n_Wrap(hgt(16:64,:,:),0),res)
  
  delete(res@cnLineColor)
  res@cnLinesOn        = True
  res@cnMonoLineColor  = True
  res@cnLineColor      = "red"
  res@cnInfoLabelOn    = False
  res@cnLineLabelsOn   = False
  res@cnFillOn         = True
  res@cnLineThicknessF = 1.3
  res@cnLineLabelFont  = "helvetica"
  
  res@mpFillOn                    = False
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  delete(res@cnLevels)
  res@cnLevels             = (/0.1,0.4,0.7,1.0,1.3/)
  res@cnFillColors         = (/65,71,81,90,100,113/)
  
  res@gsnSpreadColors      = False
  
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(a)"
  res@gsnRightString              = "LBM Forcing"
  plot(0)  = gsn_csm_contour_map_ce(wks,frc_rate,res)
  
  delete(res)
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@gsnAddCyclic = False
  
  res@gsnSpreadColors = False
  
  res@gsnContourZeroLineThicknessF = 2 ;
  res@gsnContourNegLineDashPattern = 2 ; 
  
  res@cnLinesOn        = True
  res@cnInfoLabelOn    = False
  res@cnLineLabelsOn   = False
  res@cnFillOn         = False
  res@cnLineThicknessF = 3.8
  res@cnLineLabelFont  = "helvetica"
  
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/12500/) ;
  
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  
  plot_cn = gsn_csm_contour(wks,hgt_JJA_clim,res)
  
  overlay(plot(1),plot_cn)
  
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarParallelPosF   = -0.02 ;¸ºÖµ×óÒÆ
  resp@gsnPanelRight           = 0.97
  resp@gsnPanelBottom          = 0.1
  resp@pmLabelBarOrthogonalPosF = 0.005
  resp@pmLabelBarWidthF         = 0.7
  resp@pmLabelBarHeightF        = 0.07
  resp@lbLabelFontHeightF       = 0.013
  resp@lbLabelAutoStride        = True
  
  resp@lbBoxLineColor           = 0
  resp@lbBoxEndCapStyle         = "TriangleBothEnds"
  resp@lbLabelAutoStride        = True
  
  resp@gsnPanelYF               = (/0.95,0.51/);-0.05
  resp@gsnPanelXF               = (/0.1,0.49,0.1,0.49,0.1,0.49,0.1,0.49/);-0.05
  
  gsn_panel(wks,plot,(/2,1/),resp) 
  
end