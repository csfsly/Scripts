load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
begin
	
  f_hgt  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc", "r")
  
  ;===
	TIME    = f_hgt->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
	;===
	            
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20190601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20190701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20190801                         ; last  
  
  iStrt_Jun   = ind(YMD.eq.ymdStrt_Jun)              ; index start
  iLast_Jun   = ind(YMD.eq.ymdLast_Jun)              ; index last 
  iStrt_Jul   = ind(YMD.eq.ymdStrt_Jul)              ; index start
  iLast_Jul   = ind(YMD.eq.ymdLast_Jul)              ; index last 
  iStrt_Aug   = ind(YMD.eq.ymdStrt_Aug)              ; index start
  iLast_Aug   = ind(YMD.eq.ymdLast_Aug)              ; index last 
  
  ;===
  hgt_Jun   = short2flt(f_hgt->hgt(iStrt_Jun:iLast_Jun:12,{200},::-1,:))  ; 1949.12 - 2011.12 63定
  hgt_Jul   = short2flt(f_hgt->hgt(iStrt_Jul:iLast_Jul:12,{200},::-1,:))  ; 1950.01 - 2012.01 63定
  hgt_Aug   = short2flt(f_hgt->hgt(iStrt_Aug:iLast_Aug:12,{200},::-1,:))  ; 1950.02 - 2012.02 63定
  ;===
  hgt_JJA   = (hgt_Jun+hgt_Jul+hgt_Aug)/3.
  
  copy_VarMeta(hgt_Jun,hgt_JJA)
  
  printVarSummary(hgt_JJA)
  
  hgt_JJA_clim = dim_avg_n_Wrap(hgt_JJA,0)
  
  ;===
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  hgt_JJA_IA   = wgt_runave_Wrap(hgt_JJA(lat|:,lon|:,time|:),wgt,1)
  ;===   
  
  ;===
  I_NWSE = dim_avg_n_Wrap(dim_avg_n_Wrap(hgt_JJA_IA({27.5:35},{50:80},:),0),0)-dim_avg_n_Wrap(dim_avg_n_Wrap(hgt_JJA_IA({20:27.5},{85:115},:),0),0)
  I_NWSE = dim_standardize(I_NWSE,1)               
  ;===
  
  ;=========================hgt EOF========================
  rad    = 0.01745329
  clat   = hgt_JJA_IA&lat           
  clat   = sqrt( cos(rad*clat) )
  
  ;===
  whgt   = hgt_JJA_IA              ; copy meta data
  ;===
  neof   = 3        ; number of EOFs
  optEOF = False       
  optETS = False
  ;===
  NH_whgt = whgt({15:40},{20:135},:)  ;======= 
  
  EOF_hgt = NH_whgt(lat|:,lon|:,time|:)
  
  optEOF      = True
  optEOF@jopt = 1
  optETS      = True
  optETS@jopt = 1

  eof_sp_hgt = eofunc_Wrap(EOF_hgt, neof, optEOF)      
  eof_ts_hgt = eofunc_ts_Wrap(EOF_hgt, eof_sp_hgt, optETS)
  
  print(eof_sp_hgt@pcvar)
  
  EOF1_PC     = dim_standardize(eof_ts_hgt(0,:),1)*(-1.)
	EOF2_PC     = dim_standardize(eof_ts_hgt(1,:),1)*(1.)
	EOF3_PC     = dim_standardize(eof_ts_hgt(2,:),1)
	
	;===
	EOF1_Mode   = eof_sp_hgt(0,:,:) * dim_stddev(eof_ts_hgt(0,:))*(-10.)
	EOF2_Mode   = eof_sp_hgt(1,:,:) * dim_stddev(eof_ts_hgt(1,:))*(10.)
  EOF3_Mode   = eof_sp_hgt(2,:,:) * dim_stddev(eof_ts_hgt(2,:))      
  
  copy_VarMeta(eof_sp_hgt(0,:,:),EOF1_Mode)
  copy_VarMeta(eof_sp_hgt(1,:,:),EOF2_Mode)
  copy_VarMeta(eof_sp_hgt(2,:,:),EOF3_Mode)
  
  ;===
  print(escorc(EOF1_PC,I_NWSE))
  print(escorc(EOF2_PC,I_NWSE))
  ;===
  reg_coef_hgt_PC1  = regCoef(EOF1_PC(:) ,hgt_JJA_IA(lat|:,lon|:,time|:)) 
  reg_coef_hgt_PC2  = regCoef(EOF2_PC(:) ,hgt_JJA_IA(lat|:,lon|:,time|:)) 
  
  hgt_JJA_test_PC1      = new((/73,144/),"float")
  hgt_JJA_test_PC2      = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  hgt_JJA_test_PC1(j,i)   = reg_coef_hgt_PC1@tval(k)
  hgt_JJA_test_PC2(j,i)   = reg_coef_hgt_PC2@tval(k)
  k=k+1
  end do
  end do
  ;===
  reg_coef_hgt_PC1!0   = "lat"
  reg_coef_hgt_PC1!1   = "lon"
  reg_coef_hgt_PC1&lat = hgt_Jun&lat
  reg_coef_hgt_PC1&lon = hgt_Jun&lon
  
  copy_VarMeta(reg_coef_hgt_PC1,reg_coef_hgt_PC2)
  copy_VarMeta(reg_coef_hgt_PC1,hgt_JJA_test_PC1)
  copy_VarMeta(reg_coef_hgt_PC1,hgt_JJA_test_PC2)
  
  ;===
  print("start ERA-interim")
  f_hgt_ERA  = addfile("/disk_raid6/data/ERA-interim/pressure/monthly/hgt.mon.mean.1.0.nc", "r")
  
  delete(TIME)
  delete(YMD)
  
  ;===
	TIME    = f_hgt_ERA->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
	;===
	            
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20190601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20190701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20190801                         ; last  
  
  iStrt_Jun   = ind(YMD.eq.ymdStrt_Jun)              ; index start
  iLast_Jun   = ind(YMD.eq.ymdLast_Jun)              ; index last 
  iStrt_Jul   = ind(YMD.eq.ymdStrt_Jul)              ; index start
  iLast_Jul   = ind(YMD.eq.ymdLast_Jul)              ; index last 
  iStrt_Aug   = ind(YMD.eq.ymdStrt_Aug)              ; index start
  iLast_Aug   = ind(YMD.eq.ymdLast_Aug)              ; index last 
  
  ;===
  hgt_ERA_Jun   = short2flt(f_hgt_ERA->z(iStrt_Jun:iLast_Jun:12,{200},::-1,:))  ; 1949.12 - 2011.12 63定
  hgt_ERA_Jul   = short2flt(f_hgt_ERA->z(iStrt_Jul:iLast_Jul:12,{200},::-1,:))  ; 1950.01 - 2012.01 63定
  hgt_ERA_Aug   = short2flt(f_hgt_ERA->z(iStrt_Aug:iLast_Aug:12,{200},::-1,:))  ; 1950.02 - 2012.02 63定
  ;===
  hgt_ERA_JJA   = (hgt_ERA_Jun+hgt_ERA_Jul+hgt_ERA_Aug)/3.
  
  copy_VarMeta(hgt_ERA_Jun,hgt_ERA_JJA)
  
  printVarSummary(hgt_ERA_JJA)
  
  hgt_ERA_JJA_clim = dim_avg_n_Wrap(hgt_ERA_JJA,0)
  
  ;===
  hgt_ERA_JJA_IA   = wgt_runave_Wrap(hgt_ERA_JJA(latitude|:,longitude|:,time|:),wgt,1)
  
  ;=========================hgt_ERA EOF========================
  delete(clat)
  rad    = 0.01745329
  clat   = hgt_ERA_JJA_IA&latitude           
  clat   = sqrt( cos(rad*clat) )
  
  ;===
  whgt_ERA   = hgt_ERA_JJA_IA              ; copy meta data
  ;===
  neof   = 3        ; number of EOFs
  optEOF = False       
  optETS = False
  ;===
  NH_whgt_ERA = whgt_ERA({15:40},{20:135},:)  ;======= 
  
  EOF_hgt_ERA = NH_whgt_ERA(latitude|:,longitude|:,time|:)
  
  optEOF      = True
  optEOF@jopt = 1
  optETS      = True
  optETS@jopt = 1

  eof_sp_hgt_ERA = eofunc_Wrap(EOF_hgt_ERA, neof, optEOF)      
  eof_ts_hgt_ERA = eofunc_ts_Wrap(EOF_hgt_ERA, eof_sp_hgt_ERA, optETS)
  
  print(eof_sp_hgt_ERA@pcvar)
  
  EOF1_ERA_PC     = dim_standardize(eof_ts_hgt_ERA(0,:),1)*(-1.)
	EOF2_ERA_PC     = dim_standardize(eof_ts_hgt_ERA(1,:),1)*(1.)
	EOF3_ERA_PC     = dim_standardize(eof_ts_hgt_ERA(2,:),1)
	
	;===
	EOF1_ERA_Mode   = eof_sp_hgt_ERA(0,:,:) * dim_stddev(eof_ts_hgt_ERA(0,:))*(-10.)
	EOF2_ERA_Mode   = eof_sp_hgt_ERA(1,:,:) * dim_stddev(eof_ts_hgt_ERA(1,:))*(10.)
  EOF3_ERA_Mode   = eof_sp_hgt_ERA(2,:,:) * dim_stddev(eof_ts_hgt_ERA(2,:))      
  
  copy_VarMeta(eof_sp_hgt_ERA(0,:,:),EOF1_ERA_Mode)
  copy_VarMeta(eof_sp_hgt_ERA(1,:,:),EOF2_ERA_Mode)
  copy_VarMeta(eof_sp_hgt_ERA(2,:,:),EOF3_ERA_Mode)
  
  ;===
  colors =   (/(/255,255,255/),(/0,0,0/),(/0,255,0/),(/0,0,120/),(/0,0,160/),(/0,0,210/),(/0,0,255/),(/70,70,255/),(/150,150,255/),(/220,220,255/)\
           ,(/255,255,255/),(/255,220,220/),(/255,150,150/),(/255,70,70/),(/255,0,0/),(/210,0,0/),(/160,0,0/),(/120,0,0/),(/110,110,110/) /)*1.0
           
  colors = colors / 255.0
  ;===
  wks  = gsn_open_wks("eps","Fig.1")
  gsn_define_colormap(wks,"MPL_bwr");
  
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@gsnAddCyclic = False
  
  res@gsnSpreadColors = True
  
  res@gsnContourZeroLineThicknessF = 2 ;
  res@gsnContourNegLineDashPattern = 2 ; 
  
  res@cnLinesOn        = True
  res@cnInfoLabelOn    = False
  res@cnLineLabelsOn   = False
  res@cnFillOn         = True
  res@cnLineThicknessF = 0.3
  res@cnLineLabelFont  = "helvetica"
  
  res@tiMainString       = " "
  res@tiYAxisString      = " "
  
  res@cnLineLabelFontHeightF   = 0.01
  res@cnLineLabelPlacementMode = "Computed"
  res@cnLineLabelInterval      = 1
  res@cnLineLabelFormat        = "@*+^sg"
  
  res@tmYRMode             = "Automatic"
  res@tmXTOn               = False
  res@tmYROn               = False
  
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = (/10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10/)
  
  res@mpFillOn                    = True
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  res@mpCenterLonF = 0
  res@mpMaxLatF    = 40.
  res@mpMinLatF    = 15
  res@mpMaxLonF    = 135
  res@mpMinLonF    = 20
  
  res@mpFillOn                    = True
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  res@mpGeophysicalLineColor      = "grey77" 
  res@mpGeophysicalLineThicknessF = 1
  
  res@gsnRightStringFontHeightF   = 0.03
  res@gsnLeftStringFontHeightF    = 0.03
  res@gsnCenterStringFontHeightF  = 0.03
  
  res@tmXBLabelFontHeightF = 0.028
  res@tmYLLabelFontHeightF = 0.028
  
  res@mpShapeMode  = "FreeAspect"
  res@vpWidthF     = 0.7
  res@vpHeightF    = 0.35
  
  res@gsnLeftStringOrthogonalPosF  = 0.03
  res@gsnRightStringOrthogonalPosF = 0.03
  res@gsnCenterStringOrthogonalPosF = 0.05
  
  plot = new(8,graphic)
  
  res@tmXBOn               = False
  
  res@lbLabelBarOn                = False
  res@gsnLeftString               = "(a)"
  res@gsnRightString              = "EOF1 (53%)"
  res@gsnCenterString             = "NCEP"
  plot(0) = gsn_csm_contour_map_ce(wks,EOF1_Mode,res)
  
  res@gsnLeftString               = "(c)"
  res@gsnRightString              = "EOF2 (24%)"
  plot(2) = gsn_csm_contour_map_ce(wks,EOF2_Mode,res)
  
  res@gsnCenterString             = "ERA"
  res@gsnLeftString               = "(e)"
  res@gsnRightString              = "EOF1 (53%)"
  plot(4) = gsn_csm_contour_map_ce(wks,EOF1_ERA_Mode,res)
  
  res@lbLabelBarOn             = False
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  
  res@tmXBOn               = True
  res@gsnLeftString               = "(g)"
  res@gsnRightString              = "EOF2 (25%)"
  plot(6) = gsn_csm_contour_map_ce(wks,EOF2_ERA_Mode,res)
  
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
  
  plot_cn = new(4,graphic)
  plot_cn(0) = gsn_csm_contour(wks,hgt_JJA_clim,res)
  plot_cn(1) = gsn_csm_contour(wks,hgt_JJA_clim,res)
  res@cnLevels             = (/122500/) ;
  
  plot_cn(2) = gsn_csm_contour(wks,hgt_ERA_JJA_clim,res)
  plot_cn(3) = gsn_csm_contour(wks,hgt_ERA_JJA_clim,res)
  
  overlay(plot(0),plot_cn(0))
  overlay(plot(2),plot_cn(1))
  overlay(plot(4),plot_cn(2))
  overlay(plot(6),plot_cn(3))
  
  res_t = True
  
  res_t@gsnFrame         = False
  res_t@gsnDraw          = False
  
  res_t@cnLinesOn        = False
  res_t@cnLineLabelsOn   = False
  res_t@cnFillOn         = False
  res_t@cnInfoLabelOn    = False
  
  res_t@tmXBLabelFontHeightF = 0.03
  res_t@tmYLLabelFontHeightF = 0.03
  
  res_t@cnFillOn             = True             ; fill contour intervals
  res_t@cnLevelSelectionMode = "ExplicitLevels"
  res_t@cnLevels             = (/-2.110,2.110/) ;
  
  res_t@lbLabelBarOn         = False
  
  res_t@cnMonoFillPattern    = False            ; want multiple patterns 
  res_t@cnFillPatterns       = (/17,-1,17/) ; the patterns
  res_t@cnMonoFillScale      = False            ; want different densities 
  res_t@cnFillScales         = (/1.,1,1./) ; the densities
  
  res_t@cnFillColors         = (/1,3,1/)
  
  res_t@cnFillDotSizeF       = 0.002            ; increase dot size of contours
  
  plot_test     = new(2,graphic)
  
  plot_test(0)  = gsn_csm_contour(wks,hgt_JJA_test_PC1   ,res_t)
  plot_test(1)  = gsn_csm_contour(wks,hgt_JJA_test_PC2    ,res_t)
  
  res_ts   = True
  
  res_ts@gsnDraw               = False
  res_ts@gsnFrame              = False
  
  res_ts@gsnXYBarChartBarWidth = 0.8
  
  res_ts@gsnYRefLine           = 0.              ; reference line   
  res_ts@gsnXYBarChart         = True            ; create bar chart 
  res_ts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  res_ts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
  
  res_ts@trYMinF   =  -3        
  res_ts@trYMaxF   =  3.3        
  
  res_ts@trXMinF   =  1978        
  res_ts@trXMaxF   =  2020       
  
  res_ts@vpWidthF     = 0.7
  res_ts@vpHeightF    = 0.35
  
  res_ts@tmYRMode             = "Automatic"
  res_ts@tmXTOn               = False
  res_ts@tmYROn               = False
  
  res_ts@tmXBLabelFontHeightF = 0.03
  res_ts@tmYLLabelFontHeightF = 0.03
  
  res_ts@tmYLMode          = "Explicit"
  res_ts@tmYLValues        = (/"-3","-2","-1","0","1","2","3"/)
  res_ts@tmYLLabels        = (/"-3","-2","-1","0","1","2","3"/)
  
  ;===
  year = fspan(1979,2019,41)
  
  ;===standarize
	EOF1_PC     = dim_standardize(EOF1_PC,1)
	EOF2_PC     = dim_standardize(EOF2_PC,1) 
	;===
  
  res_ts@gsnRightStringFontHeightF   = 0.03
  res_ts@gsnLeftStringFontHeightF    = 0.03
  res_ts@gsnCenterStringFontHeightF  = 0.03
  
  res_ts@gsnLeftStringOrthogonalPosF   = 0.03
  res_ts@gsnRightStringOrthogonalPosF  = 0.03
  res_ts@gsnCenterStringOrthogonalPosF = 0.05
  
  res_ts@tmXBMajorOutwardLengthF = 0.022	; set tickmarks 
  res_ts@tmXBMinorOutwardLengthF = 0.007
  res_ts@tmXBMajorLengthF        = 0.022
  res_ts@tmXBMinorLengthF        = 0.007
  
  res_ts@tmYLMajorOutwardLengthF = 0.017	; set tickmarks 
  res_ts@tmYLMinorOutwardLengthF = 0.007
  res_ts@tmYLMajorLengthF        = 0.017
  res_ts@tmYLMinorLengthF        = 0.007
  
  res_ts@tiYAxisString      = ""
  res_ts@tiYAxisFontHeightF = 0.023
  res_ts@tiYAxisOffsetXF    = 0.011
  
  res_ts@tmXBOn               = False
  
  res_ts@gsnLeftString               = "(b) "
  res_ts@gsnRightString              = "PC1"
  res_ts@gsnCenterString             = "NCEP"
  plot(1) = gsn_csm_xy (wks,year,EOF1_PC,res_ts)
  
  res_ts@tiYAxisString      = ""
  res_ts@gsnLeftString               = "(d) "
  res_ts@gsnRightString              = "PC2"
  plot(3) = gsn_csm_xy (wks,year,EOF2_PC,res_ts)
  
  res_ts@tiYAxisString      = ""
  res_ts@gsnLeftString               = "(f) "
  res_ts@gsnRightString              = "PC1"
  res_ts@gsnCenterString             = "ERA"
  plot(5) = gsn_csm_xy (wks,year,EOF1_ERA_PC,res_ts)
  
  res_ts@tmXBOn               = True
  
  res_ts@tiYAxisString      = ""
  res_ts@gsnLeftString               = "(h) "
  res_ts@gsnRightString              = "PC2"
  plot(7) = gsn_csm_xy (wks,year,EOF2_ERA_PC,res_ts)
  ;===
  ;===
  resp = True
  
  resp@gsnPanelLabelBar         = True
  resp@pmLabelBarOrthogonalPosF = 0.005
  resp@pmLabelBarParallelPosF   = -0.21
  
  resp@pmLabelBarWidthF         = 0.3
  resp@pmLabelBarHeightF        = 0.06
  resp@lbLabelFontHeightF       = 0.015
  
  resp@gsnPanelBottom          = 0.1
  
  resp@lbBoxLineColor           = 0
  resp@lbLabelAutoStride        = True
  resp@lbBoxLineThicknessF       = 3
  resp@gsnPanelRight           = 0.97
  
  resp@gsnPanelYF               = (/0.95,0.95,0.73,0.73,0.51,0.51,0.29,0.29/);-0.05
  resp@gsnPanelXF               = (/0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5/);-0.05
  
  gsn_panel(wks,plot,(/4,2/),resp) 
  
end