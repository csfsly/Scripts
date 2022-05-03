load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	;===
	SMJI = fbindirread ("~/SJI_AM_Indep-1979-2019.dat" ,0,41,"float")
  SMJI = dim_standardize(SMJI,1)
	
	;======
	print("...read in data files...")
	;===
  f_slp  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface/slp.mon.mean.nc", "r")
  
	TIME    = f_slp->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
	ymdStrt_Dec = 19781201                         ; start yyyymmdd
  ymdLast_Dec = 20181201                          
  ymdStrt_Jan = 19790101                         ; start yyyymmdd
  ymdLast_Jan = 20190101                         ; last  
  ymdStrt_Feb = 19790201                         ; start yyyymmdd
  ymdLast_Feb = 20190201                         ; last  
  
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
  
  slp_Dec   = f_slp->slp(iStrt_Dec:iLast_Dec:12,::-1,:) 
  slp_Jan   = f_slp->slp(iStrt_Jan:iLast_Jan:12,::-1,:) 
  slp_Feb   = f_slp->slp(iStrt_Feb:iLast_Feb:12,::-1,:) 
  
  slp_Mar   = short2flt(f_slp->slp(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  slp_Apr   = short2flt(f_slp->slp(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  slp_May   = short2flt(f_slp->slp(iStrt_May:iLast_May:12,::-1,:)) 
  
  slp_Jun   = f_slp->slp(iStrt_Jun:iLast_Jun:12,::-1,:) 
  slp_Jul   = f_slp->slp(iStrt_Jul:iLast_Jul:12,::-1,:) 
  slp_Aug   = f_slp->slp(iStrt_Aug:iLast_Aug:12,::-1,:) 
  
  slp_AM   = (slp_Apr+slp_May)/2.
  slp_JF   = (slp_Dec+slp_Jan)/2.
  slp_MA   = (slp_Mar+slp_Feb)/2.
  
  copy_VarMeta(slp_Mar,slp_AM)
  copy_VarMeta(slp_Mar,slp_JF)
  copy_VarMeta(slp_Mar,slp_MA)
  
  slp_AM_clim = dim_avg_n_Wrap(slp_AM,0)
  slp_JF_clim = dim_avg_n_Wrap(slp_JF,0)
  slp_MA_clim = dim_avg_n_Wrap(slp_MA,0)
  
	;===
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  slp_AM_IA    = wgt_runave_Wrap(slp_AM(lat|:,lon|:,time|:),wgt,1)
  slp_JF_IA    = wgt_runave_Wrap(slp_JF(lat|:,lon|:,time|:),wgt,1)
  slp_MA_IA    = wgt_runave_Wrap(slp_MA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_slp_AM = regCoef(SMJI,slp_AM_IA(lat|:,lon|:,time|:))
  reg_slp_JF = regCoef(SMJI,slp_JF_IA(lat|:,lon|:,time|:))
  reg_slp_MA = regCoef(SMJI,slp_MA_IA(lat|:,lon|:,time|:))
  
  slp_AM_test = new((/73,144/),"float")
  slp_JF_test = new((/73,144/),"float")
  slp_MA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  slp_AM_test(j,i) = reg_slp_AM@tval(k)
  slp_JF_test(j,i) = reg_slp_JF@tval(k)
  slp_MA_test(j,i) = reg_slp_MA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_slp_AM!0   = "lat"
  reg_slp_AM!1   = "lon"
  reg_slp_AM&lat = slp_JF_IA&lat
  reg_slp_AM&lon = slp_JF_IA&lon
  
  copy_VarMeta(reg_slp_AM,slp_AM_test)
  copy_VarMeta(reg_slp_AM,slp_JF_test)
  copy_VarMeta(reg_slp_AM,slp_MA_test)
  
  copy_VarMeta(reg_slp_AM,reg_slp_JF)
  copy_VarMeta(reg_slp_AM,reg_slp_MA)
  ;===
  ;======
	print("...read in data files...")
	delete(TIME)
	delete(YMD)
	
  f_uwnd  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/uwnd.mon.mean.nc", "r")
  f_vwnd  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc", "r")
  
	TIME    = f_uwnd->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
	ymdStrt_Dec = 19781201                         ; start yyyymmdd
  ymdLast_Dec = 20181201                          
  ymdStrt_Jan = 19790101                         ; start yyyymmdd
  ymdLast_Jan = 20190101                         ; last  
  ymdStrt_Feb = 19790201                         ; start yyyymmdd
  ymdLast_Feb = 20190201                         ; last  
  
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
  
  uwnd_Dec   = short2flt(f_uwnd->uwnd(iStrt_Dec:iLast_Dec:12,{925},::-1,:)) 
  uwnd_Jan   = short2flt(f_uwnd->uwnd(iStrt_Jan:iLast_Jan:12,{925},::-1,:)) 
  uwnd_Feb   = short2flt(f_uwnd->uwnd(iStrt_Feb:iLast_Feb:12,{925},::-1,:)) 
  uwnd_Mar   = short2flt(f_uwnd->uwnd(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  uwnd_Apr   = short2flt(f_uwnd->uwnd(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  uwnd_May   = short2flt(f_uwnd->uwnd(iStrt_May:iLast_May:12,{925},::-1,:)) 
  uwnd_Jun   = short2flt(f_uwnd->uwnd(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  uwnd_Jul   = short2flt(f_uwnd->uwnd(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  uwnd_Aug   = short2flt(f_uwnd->uwnd(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	uwnd_AM   = (uwnd_Apr+uwnd_May)/2.
  uwnd_JF   = (uwnd_Jan+uwnd_Dec)/2.
  uwnd_MA   = (uwnd_Mar+uwnd_Feb)/2.
  
  copy_VarMeta(uwnd_Mar,uwnd_AM)
  copy_VarMeta(uwnd_Mar,uwnd_JF)
  copy_VarMeta(uwnd_Mar,uwnd_MA)
	;===
  uwnd_AM_IA    = wgt_runave_Wrap(uwnd_AM(lat|:,lon|:,time|:),wgt,1)
  uwnd_JF_IA    = wgt_runave_Wrap(uwnd_JF(lat|:,lon|:,time|:),wgt,1)
  uwnd_MA_IA    = wgt_runave_Wrap(uwnd_MA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_uwnd_AM = regCoef(SMJI,uwnd_AM_IA(lat|:,lon|:,time|:))
  reg_uwnd_JF = regCoef(SMJI,uwnd_JF_IA(lat|:,lon|:,time|:))
  reg_uwnd_MA = regCoef(SMJI,uwnd_MA_IA(lat|:,lon|:,time|:))
  
  uwnd_AM_test = new((/73,144/),"float")
  uwnd_JF_test = new((/73,144/),"float")
  uwnd_MA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  uwnd_AM_test(j,i) = reg_uwnd_AM@tval(k)
  uwnd_JF_test(j,i) = reg_uwnd_JF@tval(k)
  uwnd_MA_test(j,i) = reg_uwnd_MA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_uwnd_AM!0   = "lat"
  reg_uwnd_AM!1   = "lon"
  reg_uwnd_AM&lat = uwnd_JF_IA&lat
  reg_uwnd_AM&lon = uwnd_JF_IA&lon
  
  copy_VarMeta(reg_uwnd_AM,uwnd_AM_test)
  copy_VarMeta(reg_uwnd_AM,uwnd_JF_test)
  copy_VarMeta(reg_uwnd_AM,uwnd_MA_test)
  
  copy_VarMeta(reg_uwnd_AM,reg_uwnd_JF)
  copy_VarMeta(reg_uwnd_AM,reg_uwnd_MA)
  ;===
  vwnd_Dec   = short2flt(f_vwnd->vwnd(iStrt_Dec:iLast_Dec:12,{925},::-1,:)) 
  vwnd_Jan   = short2flt(f_vwnd->vwnd(iStrt_Jan:iLast_Jan:12,{925},::-1,:)) 
  vwnd_Feb   = short2flt(f_vwnd->vwnd(iStrt_Feb:iLast_Feb:12,{925},::-1,:)) 
  vwnd_Mar   = short2flt(f_vwnd->vwnd(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  vwnd_Apr   = short2flt(f_vwnd->vwnd(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  vwnd_May   = short2flt(f_vwnd->vwnd(iStrt_May:iLast_May:12,{925},::-1,:)) 
  vwnd_Jun   = short2flt(f_vwnd->vwnd(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  vwnd_Jul   = short2flt(f_vwnd->vwnd(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  vwnd_Aug   = short2flt(f_vwnd->vwnd(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	vwnd_AM   = (vwnd_Apr+vwnd_May)/2.
  vwnd_JF   = (vwnd_Jan+vwnd_Dec)/2.
  vwnd_MA   = (vwnd_Mar+vwnd_Feb)/2.
  
  copy_VarMeta(vwnd_Mar,vwnd_AM)
  copy_VarMeta(vwnd_Mar,vwnd_JF)
  copy_VarMeta(vwnd_Mar,vwnd_MA)
	;===
  vwnd_AM_IA    = wgt_runave_Wrap(vwnd_AM(lat|:,lon|:,time|:),wgt,1)
  vwnd_JF_IA    = wgt_runave_Wrap(vwnd_JF(lat|:,lon|:,time|:),wgt,1)
  vwnd_MA_IA    = wgt_runave_Wrap(vwnd_MA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_vwnd_AM = regCoef(SMJI,vwnd_AM_IA(lat|:,lon|:,time|:))
  reg_vwnd_JF = regCoef(SMJI,vwnd_JF_IA(lat|:,lon|:,time|:))
  reg_vwnd_MA = regCoef(SMJI,vwnd_MA_IA(lat|:,lon|:,time|:))
  
  vwnd_AM_test = new((/73,144/),"float")
  vwnd_JF_test = new((/73,144/),"float")
  vwnd_MA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  vwnd_AM_test(j,i) = reg_vwnd_AM@tval(k)
  vwnd_JF_test(j,i) = reg_vwnd_JF@tval(k)
  vwnd_MA_test(j,i) = reg_vwnd_MA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_vwnd_AM!0   = "lat"
  reg_vwnd_AM!1   = "lon"
  reg_vwnd_AM&lat = vwnd_JF_IA&lat
  reg_vwnd_AM&lon = vwnd_JF_IA&lon
  
  copy_VarMeta(reg_vwnd_AM,vwnd_AM_test)
  copy_VarMeta(reg_vwnd_AM,vwnd_JF_test)
  copy_VarMeta(reg_vwnd_AM,vwnd_MA_test)
  
  copy_VarMeta(reg_vwnd_AM,reg_vwnd_JF)
  copy_VarMeta(reg_vwnd_AM,reg_vwnd_MA)
  ;===
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.14")
  gsn_define_colormap(wks,"MPL_bwr")
  
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
  res@vpWidthF             = 0.72
  res@vpHeightF            = 0.35
  
  res@mpFillOn                    = True
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  res@mpCenterLonF = 180
  res@mpMaxLatF    = 30.
  res@mpMinLatF    = -70
  res@mpMaxLonF    = 120
  res@mpMinLonF    = 20
  
  res@tmXBLabelFontHeightF = 0.026
  res@tmYLLabelFontHeightF = 0.026
  
  res@gsnLeftStringFontHeightF   = 0.028
  res@gsnRightStringFontHeightF  = 0.028
  res@gsnCenterStringFontHeightF = 0.028
  
  res@gsnCenterStringOrthogonalPosF = 0.025
  res@gsnRightStringOrthogonalPosF  = 0.015
  res@gsnLeftStringOrthogonalPosF   = 0.015
  
  res@lbLabelBarOn             = False
  
  plot = new(6,graphic)
  
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = (/-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9/)
  
  res@tmXBOn               = False
  res@gsnCenterString             = " "
  res@gsnLeftString               = "(a) "
  res@gsnRightString              = "D(-1)J(0)"
  plot(0)  = gsn_csm_contour_map_ce(wks,reg_slp_JF,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(b) "
  res@gsnRightString              = "FM(0)"
  plot(2)  = gsn_csm_contour_map_ce(wks,reg_slp_MA,res)
  
  res@gsnCenterString             = ""
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 5
  res@pmLabelBarOrthogonalPosF = 0.2
  res@pmLabelBarWidthF         = 0.6
  res@pmLabelBarHeightF        = 0.07
  res@lbLabelFontHeightF       = 0.025
  res@gsnLeftString               = "(c) "
  res@gsnRightString              = "AM(0)"
  res@tmXBOn               = True
  plot(4)  = gsn_csm_contour_map_ce(wks,reg_slp_AM,res)
  
  delete(res)
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
  res@cnLineThicknessF = 1.8
  res@cnLineLabelFont  = "helvetica"
  
  res@cnLineLabelFontHeightF   = 0.01
  res@cnLineLabelPlacementMode = "Computed"
  res@cnLineLabelInterval      = 1
  res@cnLineLabelFormat        = "@*+^sg"
  
  res@mpFillOn                    = False
  res@mpDataBaseVersion           = "Lowres"        ;Lowres, Mediumres
  res@mpDataSetName               = "Earth..4"
  res@mpFillAreaSpecifiers        = (/"Water","Land"/)
  res@mpSpecifiedFillColors       = (/0,0/)
  res@mpFillDrawOrder             = "PreDraw"
  
  res@mpShapeMode          = "FreeAspect"
  res@vpWidthF             = 0.72
  res@vpHeightF            = 0.35
  
  res@mpCenterLonF = 180
  res@mpMaxLatF    = 30.
  res@mpMinLatF    = -70
  res@mpMaxLonF    = 120
  res@mpMinLonF    = 20
  
  res@lbLabelBarOn                = False
  res@gsnCenterString             = "slp(850hPa)"
  res@gsnLeftString               = "(d) "
  res@gsnRightString              = "D(-1)J(0)"
  plot(1)  = gsn_csm_contour_map_ce(wks,slp_JF_clim,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(e) "
  res@gsnRightString              = "FM(0)"
  plot(3)  = gsn_csm_contour_map_ce(wks,slp_MA_clim,res)
  res@gsnLeftString               = "(f) "
  res@gsnRightString              = "MA(0)"
  res@gsnCenterString             = ""
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.2
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  res@tmXBOn                   = True
  plot(5)  = gsn_csm_contour_map_ce(wks,slp_AM_clim,res)
  
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
  res_t@cnLevels             = (/-2,2/) ;
  
  res_t@lbLabelBarOn         = False
  
  res_t@cnMonoFillPattern    = False            ; want multiple patterns 
  res_t@cnFillPatterns       = (/17,-1,17/) ; the patterns
  res_t@cnMonoFillScale      = False            ; want different densities 
  res_t@cnFillScales         = (/1,1,1/) ; the densities
  res_t@cnFillColors         = (/1,3,1/)
  res_t@cnFillDotSizeF       = 0.002            ; increase dot size of contours
  
  plot_test     = new(3,graphic)
  
  plot_test(0)   = gsn_csm_contour(wks,slp_JF_test    ,res_t)
  plot_test(1)   = gsn_csm_contour(wks,slp_MA_test    ,res_t)
  plot_test(2)   = gsn_csm_contour(wks,slp_AM_test    ,res_t)
  
  overlay(plot(0),plot_test(0))
  overlay(plot(2),plot_test(1))
  overlay(plot(4),plot_test(2))
  
  ;===
  res_vc = True

  res_vc@gsnDraw  = False
  res_vc@gsnFrame = False
  
  res_vc@tmYRMode             = "Automatic"
  res_vc@tmXTOn               = False
  res_vc@tmYROn               = False
  
  res_vc@vcRefAnnoString2        = ""
  res_vc@vcRefAnnoFont           = "helvetica"
  res_vc@vcRefAnnoFontHeightF    = 0.018
  res_vc@vcRefAnnoSide           = "Bottom"
  res_vc@vcRefAnnoOrthogonalPosF = -1.
  res_vc@vcRefAnnoParallelPosF   = 1.
  res_vc@vcMinFracLengthF        = 0.2
  res_vc@vcRefLengthF            = 0.055
  
  reg_uwnd_JF = mask(reg_uwnd_JF,reg_uwnd_JF.lt.0.2 .and. reg_uwnd_JF.gt.-0.2.and.reg_vwnd_JF.lt.0.1 .and. reg_vwnd_JF.gt.-0.1,False)
  reg_uwnd_MA = mask(reg_uwnd_MA,reg_uwnd_MA.lt.0.2 .and. reg_uwnd_MA.gt.-0.2.and.reg_vwnd_MA.lt.0.1 .and. reg_vwnd_MA.gt.-0.1,False)
  reg_uwnd_AM = mask(reg_uwnd_AM,reg_uwnd_AM.lt.0.2 .and. reg_uwnd_AM.gt.-0.2.and.reg_vwnd_AM.lt.0.1 .and. reg_vwnd_AM.gt.-0.1,False)
                                                                                                                                   
  reg_vwnd_JF  = mask(reg_vwnd_JF,reg_uwnd_JF.lt.0.2 .and. reg_uwnd_JF.gt.-0.2.and.reg_vwnd_JF.lt.0.1 .and. reg_vwnd_JF.gt.-0.1,False)
  reg_vwnd_MA  = mask(reg_vwnd_MA,reg_uwnd_MA.lt.0.2 .and. reg_uwnd_MA.gt.-0.2.and.reg_vwnd_MA.lt.0.1 .and. reg_vwnd_MA.gt.-0.1,False)
  reg_vwnd_AM  = mask(reg_vwnd_AM,reg_uwnd_AM.lt.0.2 .and. reg_uwnd_AM.gt.-0.2.and.reg_vwnd_AM.lt.0.1 .and. reg_vwnd_AM.gt.-0.1,False)
  
  plot_vc = new(4,graphic)
    
  res_vc@vcRefMagnitudeF         = 0.5
  
  res_vc@gsnRightStringFontHeightF   = 0.025
  res_vc@gsnLeftStringFontHeightF    = 0.025
  
  res_vc@tmXBLabelFontHeightF = 0.025
  res_vc@tmYLLabelFontHeightF = 0.025
  
  res_vc@gsnLeftString               = ""
  res_vc@gsnRightString              = " "
  plot_vc(0)  = gsn_csm_vector(wks,reg_uwnd_JF(::2,::2),reg_vwnd_JF(::2,::2),res_vc)
  
  plot_vc(1)  = gsn_csm_vector(wks,reg_uwnd_MA(::2,::2),reg_vwnd_MA(::2,::2),res_vc)
  res_vc@vcRefAnnoOrthogonalPosF = -1.4
  plot_vc(2)  = gsn_csm_vector(wks,reg_uwnd_AM(::2,::2),reg_vwnd_AM(::2,::2),res_vc)
  
  plot_vc(3)  = gsn_csm_vector(wks,reg_uwnd_AM(::2,::2),reg_vwnd_AM(::2,::2),res_vc)
  
  overlay(plot(0),plot_vc(0))
  overlay(plot(2),plot_vc(1))
  overlay(plot(4),plot_vc(2))
  
  rest                   = True                      ; polyline mods desired
  rest@gsLineColor       = "red"                    ; color of lines
  rest@gsLineThicknessF  = 2                       ; thickness of lines
  rest@gsLineDashPattern = 0
  
  ;===
  ypts_A = (/ 5.0,  5.0,  -5.0,  -5.0, 5.0/)      ;latitude
  xpts_A = (/ 40.0,  50.0,  50.0,   40.0,  40.0/)    ;longitude
  
  dum_A = new(4,graphic)
  dum_B = new(4,graphic)
  dum_C = new(4,graphic)
  dum_D = new(4,graphic)
  
  do i = 0, 3
   dum_C(i)=gsn_add_polyline(wks,plot(4),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
  end do
  
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarParallelPosF   = -0.02 
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
  
  resp@gsnPanelYF               = (/0.95,0.95,0.7,0.7,0.45,0.45/)
  resp@gsnPanelXF               = (/0.1,0.54,0.1,0.54,0.1,0.54/)
  
  gsn_panel(wks,(/plot(0),plot(2),plot(4)/),(/3,1/),resp)    
  
end