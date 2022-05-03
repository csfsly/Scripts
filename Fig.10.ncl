load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	;===
	PC2_NCEP   = fbindirread ("~/hgt200-PC2-NCEP1-1979-2019.dat" ,0,41,"float")
	INWSE_NCEP = fbindirread ("~/Index-NWSE-NCEP1-1979-2019.dat" ,0,41,"float")
	
	PC2_NCEP   = dim_standardize(PC2_NCEP,1)
	INWSE_NCEP = dim_standardize(INWSE_NCEP,1)
	
	SMJI = fbindirread ("~/SJI_AM_Indep-1979-2019.dat" ,0,41,"float")
  SMJI = dim_standardize(SMJI,1)
	
	print(escorc(INWSE_NCEP,SMJI))
	
	;======
	print("...read in data files...")
	
  f_precip  =  addfile("~/UDEL/precip.mon.total.v501.nc", "r")
	
	TIME    = f_precip->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
  ymdStrt_Jun = 19790601                         ; start yyyymmdd
  ymdLast_Jun = 20170601                          
  ymdStrt_Jul = 19790701                         ; start yyyymmdd
  ymdLast_Jul = 20170701                         ; last  
  ymdStrt_Aug = 19790801                         ; start yyyymmdd
  ymdLast_Aug = 20170801                         ; last  
  
  iStrt_Jun   = ind(YMD.eq.ymdStrt_Jun)              ; index start
  iLast_Jun   = ind(YMD.eq.ymdLast_Jun)              ; index last 
  iStrt_Jul   = ind(YMD.eq.ymdStrt_Jul)              ; index start
  iLast_Jul   = ind(YMD.eq.ymdLast_Jul)              ; index last 
  iStrt_Aug   = ind(YMD.eq.ymdStrt_Aug)              ; index start
  iLast_Aug   = ind(YMD.eq.ymdLast_Aug)              ; index last 
  
  precip_Jun   = short2flt(f_precip->precip(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  precip_Jul   = short2flt(f_precip->precip(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  precip_Aug   = short2flt(f_precip->precip(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
  precip_JJA  = (precip_Jun+precip_Jul+precip_Aug)/3.
  
  copy_VarMeta(precip_Jun,precip_JJA)
	;===
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  precip_JJA_IA   = wgt_runave_Wrap(precip_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===define precipitation index
  IHeat    = dim_avg_n_Wrap(dim_avg_n_Wrap(precip_JJA_IA({10:34},{75:84},:),0),0)
  IHeat    = dim_standardize(IHeat,1)
  
  ;===regression
  reg_precip_JJA = regCoef(SMJI(1979-1979:2017-1979),precip_JJA_IA(lat|:,lon|:,time|:))
  
  precip_JJA_test = new((/360,720/),"float")
  
  k = 0 
  do j=0,359
  do i=0,719
  precip_JJA_test(j,i) = reg_precip_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_precip_JJA!0   = "lat"
  reg_precip_JJA!1   = "lon"
  reg_precip_JJA&lat = precip_JJA_IA&lat
  reg_precip_JJA&lon = precip_JJA_IA&lon
  
  copy_VarMeta(reg_precip_JJA,precip_JJA_test)
  
  ;===
  f_omega  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/omega.mon.mean.nc", "r")
  
  delete(TIME)
  delete(YMD)
  TIME    = f_omega->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
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
  
  omega_Jun   = short2flt(f_omega->omega(iStrt_Jun:iLast_Jun:12,{500},::-1,:)) 
  omega_Jul   = short2flt(f_omega->omega(iStrt_Jul:iLast_Jul:12,{500},::-1,:)) 
  omega_Aug   = short2flt(f_omega->omega(iStrt_Aug:iLast_Aug:12,{500},::-1,:)) 
  
  omega_JJA   = (omega_Jun+omega_Jul+omega_Aug)/3. *(-1000)
  
  copy_VarMeta(omega_Jun,omega_JJA)
	;===
  omega_JJA_IA    = wgt_runave_Wrap(omega_JJA(lat|:,lon|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_omega_JJA = regCoef(SMJI,omega_JJA_IA(lat|:,lon|:,time|:))
  
  omega_JJA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  omega_JJA_test(j,i) = reg_omega_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_omega_JJA!0   = "lat"
  reg_omega_JJA!1   = "lon"
  reg_omega_JJA&lat = omega_JJA_IA&lat
  reg_omega_JJA&lon = omega_JJA_IA&lon
  
  copy_VarMeta(reg_omega_JJA,omega_JJA_test)
  ;===
  
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.10")
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
  res@mpMaxLatF    = 40.
  res@mpMinLatF    = 0
  res@mpMaxLonF    = 100
  res@mpMinLonF    = 60
  
  res@tmXBLabelFontHeightF = 0.023
  res@tmYLLabelFontHeightF = 0.023
  
  res@gsnLeftStringFontHeightF   = 0.025
  res@gsnRightStringFontHeightF  = 0.025
  res@gsnCenterStringFontHeightF = 0.025
  
  res@gsnCenterStringOrthogonalPosF = 0.1
  res@gsnRightStringOrthogonalPosF  = 0.015
  res@gsnLeftStringOrthogonalPosF   = 0.015
  
  res@lbLabelBarOn             = True
  
  plot = new(4,graphic)
  
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = (/-2.2,-1.8,-1.4,-1,-0.6,-0.2,0.2,0.6,1,1.4,1.8,2.2/)
  
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.16
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  
  res@tmXBOn               = True
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(a) "
  res@gsnRightString              = "Precip"
  plot(0)  = gsn_csm_contour_map_ce(wks,reg_precip_JJA,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(b) "
  res@gsnRightString              = "Omega(500hPa)"
  
  delete(res@cnLevels)
  res@cnLevels                    = (/-5,-4,-3,-2,-1,1,2,3,4,5/)
  
  plot(1)  = gsn_csm_contour_map_ce(wks,reg_omega_JJA,res)
  
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
  
  plot_test     = new(2,graphic)
  
  plot_test(0)   = gsn_csm_contour(wks,precip_JJA_test    ,res_t)
  plot_test(1)   = gsn_csm_contour(wks,omega_JJA_test     ,res_t)
  
  ;===
  do i=0,1
  overlay(plot(i),plot_test(i))
  end do
  
  rest                   = True                      ; polyline mods desired
  rest@gsLineColor       = "black"                    ; color of lines
  rest@gsLineThicknessF  = 2                       ; thickness of lines
  rest@gsLineDashPattern = 0
  
  ;===
  ypts_A = (/ 34.0,  34.0,  10.0,  10.0, 34.0/)      ;latitude
  xpts_A = (/ 75.0,  84.0,  84.0,   75.0,  75.0/)    ;longitude
  
  dum_A = new(4,graphic)
  
  do i = 0, 3
   dum_A(i)=gsn_add_polyline(wks,plot(0),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
  end do
  
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarParallelPosF   = -0.02 ;¸ºÖµ×óÒÆ
  resp@gsnPanelRight           = 0.97
  resp@gsnPanelBottom          = 0.07
  resp@pmLabelBarOrthogonalPosF = 0.005
  resp@pmLabelBarWidthF         = 0.7
  resp@pmLabelBarHeightF        = 0.07
  resp@lbLabelFontHeightF       = 0.013
  resp@lbLabelAutoStride        = True
  
  resp@lbBoxLineColor           = 0
  resp@lbBoxEndCapStyle         = "TriangleBothEnds"
  resp@lbLabelAutoStride        = True
  
  resp@gsnPanelYF               = (/0.95,0.5/);-0.05
  resp@gsnPanelXF               = (/0.1,0.49,0.1,0.49,0.1,0.49,0.1,0.49/);-0.05
  
  gsn_panel(wks,plot,(/2,1/),resp) 
  
end