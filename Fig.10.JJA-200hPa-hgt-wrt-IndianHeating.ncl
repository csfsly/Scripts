load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	;===
	IndianHeating   = fbindirread ("IndianHeating-JJA-1979-2017.dat" ,0,39,"float")
	
	IndianHeating   = dim_standardize(IndianHeating,1)
		
	print("...read in data files...")
	
  f_sst  =  addfile("/disk_raid6/data/noaa.ersst.v5/sst.mnmean.nc", "r")
	
	TIME    = f_sst->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
	ymdStrt_Dec = 19781201                         ; start yyyymmdd
  ymdLast_Dec = 20161201                          
  ymdStrt_Jan = 19790101                         ; start yyyymmdd
  ymdLast_Jan = 20170101                         ; last  
  ymdStrt_Feb = 19790201                         ; start yyyymmdd
  ymdLast_Feb = 20170201                         ; last  
  
  ymdStrt_Mar = 19790301                         ; start yyyymmdd
  ymdLast_Mar = 20170301                          
  ymdStrt_Apr = 19790401                         ; start yyyymmdd
  ymdLast_Apr = 20170401                         ; last  
  ymdStrt_May = 19790501                         ; start yyyymmdd
  ymdLast_May = 20170501                         ; last  
  
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20170601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20170701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20170801                         ; last  
  
  
  iStrt_Dec   = ind(YMD.eq.ymdStrt_Dec)              ; index start
  iLast_Dec   = ind(YMD.eq.ymdLast_Dec)              ; index last 
  iStrt_Jan   = ind(YMD.eq.ymdStrt_Jan)              ; index start
  iLast_Jan   = ind(YMD.eq.ymdLast_Jan)              ; index last 
  iStrt_Feb   = ind(YMD.eq.ymdStrt_Feb)              ; index start
  iLast_Feb   = ind(YMD.eq.ymdLast_Feb)              ; index last 
  
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
  
  sst_Dec   = short2flt(f_sst->sst(iStrt_Dec:iLast_Dec:12,::-1,:)) 
  sst_Jan   = short2flt(f_sst->sst(iStrt_Jan:iLast_Jan:12,::-1,:)) 
  sst_Feb   = short2flt(f_sst->sst(iStrt_Feb:iLast_Feb:12,::-1,:)) 
  
  sst_Mar   = short2flt(f_sst->sst(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  sst_Apr   = short2flt(f_sst->sst(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  sst_May   = short2flt(f_sst->sst(iStrt_May:iLast_May:12,::-1,:)) 
  
  sst_Jun   = short2flt(f_sst->sst(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  sst_Jul   = short2flt(f_sst->sst(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  sst_Aug   = short2flt(f_sst->sst(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	sst_DJF   = (sst_Dec+sst_Jan+sst_Feb)/3.
  ;sst_DJF   = (sst_May)/3.
  
  copy_VarMeta(sst_Mar,sst_DJF)
	;===
	;===高通滤波
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  sst_DJF_IA    = wgt_runave_Wrap(sst_DJF(lat|:,lon|:,time|:),wgt,1)
  
  ;===NINO34 index
	Nino34_DJF   = dim_avg_n_Wrap(dim_avg_n_Wrap(sst_DJF_IA({-5:5},{190:240},:),0),0)	
	Nino34_DJF   = dim_standardize(Nino34_DJF,1)
  
  ;===remove ENSO influence
	rc_IndianHeating      = regCoef(Nino34_DJF,IndianHeating)
	IndianHeating_Indep   = IndianHeating  - rc_IndianHeating * Nino34_DJF  ;#########
	;===
  
	;======
	print("...read in data files...")
	
  ;===
  f_hgt  = addfile("/disk_raid6/data/NCEP-2.5/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc", "r")
  
  delete(TIME)
  delete(YMD)
  
	TIME    = f_hgt->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
  ymdStrt_Mar = 19790301                         ; start yyyymmdd
  ymdLast_Mar = 20170301                          
  ymdStrt_Apr = 19790401                         ; start yyyymmdd
  ymdLast_Apr = 20170401                         ; last  
  ymdStrt_May = 19790501                         ; start yyyymmdd
  ymdLast_May = 20170501                         ; last  
  
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20170601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20170701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20170801                         ; last  
  
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
	;===高通滤波
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  hgt_JJA_IA    = wgt_runave_Wrap(hgt_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regressionIndianHeating_Indep
  reg_hgt_JJA = regCoef(IndianHeating,hgt_JJA_IA(lat|:,lon|:,time|:))
  reg_hgt_JJA_Indep = regCoef(IndianHeating_Indep,hgt_JJA_IA(lat|:,lon|:,time|:))
  
  hgt_JJA_test = new((/73,144/),"float")
  hgt_JJA_test_Indep = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  hgt_JJA_test(j,i) = reg_hgt_JJA@tval(k)
  hgt_JJA_test_Indep(j,i) = reg_hgt_JJA_Indep@tval(k)
  k = k +1 
  end do
  end do
  
  reg_hgt_JJA!0   = "lat"
  reg_hgt_JJA!1   = "lon"
  reg_hgt_JJA&lat = hgt_JJA_IA&lat
  reg_hgt_JJA&lon = hgt_JJA_IA&lon
  
  copy_VarMeta(reg_hgt_JJA,reg_hgt_JJA_Indep)
  
  copy_VarMeta(reg_hgt_JJA,hgt_JJA_test)
  copy_VarMeta(reg_hgt_JJA,hgt_JJA_test_Indep)
  ;===
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.10.JJA-200hPa-hgt-wrt-IndianHeating")
  gsn_define_colormap(wks,"BlueDarkRed18")
  
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@gsnSpreadColors = True
  
  res@gsnContourZeroLineThicknessF = 2 ;
  res@gsnContourNegLineDashPattern = 2 ; 
  
  res@cnLinesOn        = False
  res@cnMonoLineColor  = True
  res@cnLineColor      = 0
  res@cnInfoLabelOn    = False
  res@cnLineLabelsOn   = False
  res@cnFillOn         = True
  res@cnLineThicknessF = 1.3
  res@cnLineLabelFont  = "helvetica"
  
  res@tiMainString       = " "
  ;res@tiYAxisstring      = " "
  
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
  res@mpMaxLatF    = 50.
  res@mpMinLatF    = -20
  res@mpMaxLonF    = 150
  res@mpMinLonF    = 10
  
  res@tmXBLabelFontHeightF = 0.023
  res@tmYLLabelFontHeightF = 0.023
  
  res@gsnLeftStringFontHeightF   = 0.023
  res@gsnRightStringFontHeightF  = 0.023
  res@gsnCenterStringFontHeightF = 0.023
  
  res@gsnCenterStringOrthogonalPosF = 0.015
  res@gsnRightStringOrthogonalPosF  = 0.015
  res@gsnLeftStringOrthogonalPosF   = 0.015
  
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  ;res@lbBoxEndCapStyle         = "TriangleBothEnds"
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.16
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  
  ;res@tmXBOn               = False
  ;res@tmYLOn               = True
  
  ;res@tmXBMode          = "Explicit"
  ;res@tmXBValues        = (/"-180","-120","-60","0","60","120","180","240","300","360"/)
  ;res@tmXBLabels        = (/"180","120W","60W","0","60E","120E","180","120W","60W","0"/)
  
  res@tmXBOn               = True
  res@tmYLOn               = True
  res@gsnCenterString             = "HGT(200hPa) wrt South Asian Heating"
  res@gsnLeftString               = ""
  res@gsnRightString              = ""
  
  res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels             = (/-18,-14,-10,-6,-2,2,6,10,14,18/)
  res@cnLevels             = (/-13,-10,-7,-5,-2,2,5,7,10,13/)
  plot  = gsn_csm_contour_map_ce(wks,reg_hgt_JJA,res)
  
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
  
  overlay(plot,plot_cn)
  ;===
  
  res_t = True
  
  res_t@gsnFrame         = False
  res_t@gsnDraw          = False
  
  res_t@cnLinesOn        = False
  res_t@cnLineLabelsOn   = False
  res_t@cnFillOn         = False
  res_t@cnInfoLabelOn    = False
  
  res_t@cnFillOn             = True             ; fill contour intervals
  res_t@cnLevelSelectionMode = "ExplicitLevels"
  res_t@cnLevels             = (/-1.76,1.76/) ;
  
  res_t@lbLabelBarOn         = False
  
  res_t@cnMonoFillPattern    = False            ; want multiple patterns 
  res_t@cnFillPatterns       = (/17,-1,17/) ; the patterns
  res_t@cnMonoFillScale      = False            ; want different densities 
  res_t@cnFillScales         = (/1,1,1/) ; the densities
  res_t@cnFillColors         = (/1,3,1/)
  res_t@cnFillDotSizeF       = 0.002            ; increase dot size of contours
    
  plot_test   = gsn_csm_contour(wks,hgt_JJA_test    ,res_t)
  
  ;===
  overlay(plot,plot_test)
  
  draw(plot)
  frame(wks)
  
end