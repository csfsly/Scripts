load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
begin
	
	;===
  ihp   = 1
  sigma = 1.0
  nWgt  = 5
  fca   = 1./5.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  ;===
  path = "/home/csf/2021/Sjet-SAH/select-1/"
  SJI_AM  = fbindirread(path+"SJI_AM-1979-2019.dat"    ,0,41,"float")
  SJI_MAM = fbindirread(path+"SJI_MAM-1979-2019.dat"   ,0,41,"float")
  SAH_JJA = fbindirread(path+"hgt200-PC2-NCEP1-1979-2019.dat"   ,0,41,"float")
  
	SJI_AM  = dim_standardize_n_Wrap(SJI_AM,  1, 0)
	SJI_MAM = dim_standardize_n_Wrap(SJI_MAM, 1, 0)
	SAH_JJA = dim_standardize_n_Wrap(SAH_JJA, 1, 0)
  
  ;===
  path     = "/home/csf/2021/Sjet-SAH/select-1/reg_model/output/"
	pre_SAH  = asciiread(path+"cross_y.SMJ.AM.txt" ,41,"float") ;28:13
  
	SJSAH_index      = new((/3,41/),"float") ;1951-2015
  SJSAH_index(0,:) = dim_standardize_n_Wrap(SJI_MAM, 1, 0)
  SJSAH_index(1,:) = dim_standardize_n_Wrap(SAH_JJA, 1, 0)
  SJSAH_index(2,:) = dim_standardize_n_Wrap(pre_SAH, 1, 0)
  
  print("EPM")
  print(escorc(pre_SAH,SAH_JJA))
    
  ;===
  wks = gsn_open_wks("eps" ,"Fig.4-SMJ.SAH")
  
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@trXMaxF  = 2015
  res@trXMinF  = 1981
  
  res@vpWidthF  = 0.7
  res@vpHeightF = 0.35
  
  res@xyLineThicknessF  = 3
  
  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .39                 ; move units right
  res@pmLegendOrthogonalPosF = -1.22                ; move units down
  res@pmLegendWidthF         = 0.06                ; Change width and
  res@pmLegendHeightF        = 0.08                ; height of legend.
  res@lgLabelFontHeightF     = .019                 ; change font height
  res@xyExplicitLegendLabels = (/"SMJ: AM","SAH: JJA","SAH: JJA (hindcast)"/)         ; explicit labels
  res@lgPerimOn              = False 
  ;===
  
  res@xyLineColors           = (/"red","blue","black"/)
  res@xyLineThicknesses      = (/2.5,3.5,4.5/)
  res@xyDashPatterns         = (/0,0,0/)
  
  res@xyMarkLineMode    = "Lines"
  res@xyMarkers         = (/16,16,16/)     
  res@xyMonoMarkerSize  = True
  res@xyMarkerSizeF     = 0.006
  res@xyMarkerColors    = (/"red","blue","green"/)
  
  res@tiMainString       = " "
  res@tiYAxisString      = " "
  
  plot = new(2,graphic)
  
  res@tmXBLabelFontHeightF = 0.021
  res@tmYLLabelFontHeightF = 0.021
  
  res@tmYLMode   = "Explicit"
  res@tmYLValues = (/"-3","-2","-1","0","1","2","3"/)
  res@tmYLLabels = (/"-3","-2","-1","0","1","2","3"/)
  
  ;res@tmYRMode             = "Automatic"
  ;res@tmXTOn               = False
  ;res@tmYROn               = False
  
  res@tmYLMinorLengthF  = 0.005
  res@tmYLMajorLengthF  = 0.01
  res@tmXBMinorLengthF  = 0.005
  res@tmXBMajorLengthF  = 0.01
  
  year  = fspan(1979.,2019.,41)
	
	res@gsnLeftStringOrthogonalPosF   = 0.03
  res@gsnRightStringOrthogonalPosF  = 0.03
  
	res@gsnCenterStringFontHeightF   = 0.023
  res@gsnLeftStringFontHeightF     = 0.023
  res@gsnCenterString              = "Empirical Prediton Model"
  res@gsnLeftString                = ""
  
  res@trYMaxF  = 3
  res@trYMinF  = -3
  
  plot(0)       = gsn_csm_xy(wks,year(:),SJSAH_index,res)
  res@gsnRightString              = "Empirical Prediton model"
  res@gsnLeftString               = "(b)"
  
  delete(res@pmLegendDisplayMode)
  
  plot(1)       = gsn_csm_xy(wks,year(:),SJSAH_index(0,:),res)
  
  delete(res@xyLineThicknesses)
  
  delete(res@xyLineColors     )
  delete(res@xyDashPatterns   )
  
  delete(res@pmLegendDisplayMode)
  
  res@xyLineColors       = "black"
  res@xyLineThicknesses  =  1
  
  res@xyMarkLineMode    = "Lines"
  
  delete(res@gsnRightStringFontHeightF )  
  delete(res@gsnLeftStringFontHeightF  )  
  delete(res@gsnRightString            )  
  delete(res@gsnLeftString             )  
  
  res@xyDashPatterns     =  2
  plot_t1                = gsn_csm_xy(wks,year(:),SJSAH_index(0,:)*0.,res)
  plot_t2                = gsn_csm_xy(wks,year(:),SJSAH_index(0,:)*0.,res)
  res@xyDashPatterns     =  0
  
  res@xyLineColors       = "blue"
  res@xyLineThicknesses  =  3
  
  plot_train             = gsn_csm_xy(wks,year(0:27) ,pre_SAH(0:27),res)
  res@xyDashPatterns     =  3
  plot_hindcast          = gsn_csm_xy(wks,year(28:40),pre_SAH(28:40),res)
  
  overlay(plot(0),plot_t1)
  overlay(plot(1),plot_t2)
  overlay(plot(1),plot_train)
  overlay(plot(1),plot_hindcast)
  
  tres                  = True                    ; text mods desired
  tres@txFontHeightF    = 0.021                 ; change font height
  tres@txFont           = "helvetica"
  tres@txFontThicknessF = 1.0
  tres@txFontColor      = "blue"
  
  gsn_text_ndc(wks,"r = 0.57",0.68,0.5, tres)
  
  
  tres@txFontHeightF    = 0.021                  ; change font height
  tres@txFont           = "helvetica"
  tres@txFontThicknessF = 1.0
  tres@txFontColor      = "blue"
  
  gsn_text_ndc(wks,"r = 0.52",0.5,0.10, tres)
  gsn_text_ndc(wks,"hindcast",0.5,0.14, tres)
  
  gsn_text_ndc(wks,"r = 0.57",0.8,0.10, tres)
  gsn_text_ndc(wks,"forecast",0.8,0.14, tres)
  
  rest                   = True                      ; polyline mods desired
  rest@gsLineColor       = "black"                     ; color of lines
  rest@gsLineThicknessF  = 2                       ; thickness of lines
  rest@gsLineDashPattern = 12
  
  ypts_E = (/ 3,  -3  /)
  xpts_E = (/ 2007.5,  2007.5/)
   
  dum_E =gsn_add_polyline(wks,plot(1),xpts_E(0:1),ypts_E(0:1),rest) 
  
  draw(plot(0))
  frame(wks)
  
  ;gsn_panel(wks,plot,(/2,1/),False)
  
end
