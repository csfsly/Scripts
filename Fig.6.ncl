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
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  ;===
	print("...read in data files...")
	
  f_lhtfl   = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/lhtfl.sfc.mon.mean.nc","r"); Monthly Mean of Latent Heat Net Flux
  
	TIME    = f_lhtfl->time                       ; days since ...
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
  
  lhtfl_Mar   = short2flt(f_lhtfl->lhtfl(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  lhtfl_Apr   = short2flt(f_lhtfl->lhtfl(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  lhtfl_May   = short2flt(f_lhtfl->lhtfl(iStrt_May:iLast_May:12,::-1,:)) 
  
  lhtfl_Jun   = short2flt(f_lhtfl->lhtfl(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  lhtfl_Jul   = short2flt(f_lhtfl->lhtfl(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  lhtfl_Aug   = short2flt(f_lhtfl->lhtfl(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	lhtfl_AM   = (lhtfl_Apr+lhtfl_May)/2. *(-1.)     /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  lhtfl_MJ   = (lhtfl_May+lhtfl_Jun)/2. *(-1.)     /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  lhtfl_JJ   = (lhtfl_Jun+lhtfl_Jul)/2. *(-1.)     /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  lhtfl_JJA  = (lhtfl_Jun+lhtfl_Jul+lhtfl_Aug)/3. *(-1.) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month;++++++
  
  copy_VarMeta(lhtfl_Mar,lhtfl_AM)
  copy_VarMeta(lhtfl_Mar,lhtfl_MJ)
  copy_VarMeta(lhtfl_Mar,lhtfl_JJ)
  copy_VarMeta(lhtfl_Mar,lhtfl_JJA)
  
  ;===
  lhtfl_AM   = smth9_Wrap(lhtfl_AM, 0.50, 0.25, False) ; light local smoothing
  lhtfl_MJ   = smth9_Wrap(lhtfl_MJ, 0.50, 0.25, False) ; light local smoothing
  lhtfl_JJ   = smth9_Wrap(lhtfl_JJ, 0.50, 0.25, False) ; light local smoothing
  lhtfl_JJA  = smth9_Wrap(lhtfl_JJA, 0.50, 0.25, False) ; light local smoothing
	;===
  lhtfl_AM_IA    = wgt_runave_Wrap(lhtfl_AM(lat|:,lon|:,time|:),wgt,1)
  lhtfl_MJ_IA    = wgt_runave_Wrap(lhtfl_MJ(lat|:,lon|:,time|:),wgt,1)
  lhtfl_JJ_IA    = wgt_runave_Wrap(lhtfl_JJ(lat|:,lon|:,time|:),wgt,1)
  lhtfl_JJA_IA    = wgt_runave_Wrap(lhtfl_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===
  reg_lhtfl_AM = regCoef(SMJI,lhtfl_AM_IA(lat|:,lon|:,time|:))
  reg_lhtfl_MJ = regCoef(SMJI,lhtfl_MJ_IA(lat|:,lon|:,time|:))
  reg_lhtfl_JJ = regCoef(SMJI,lhtfl_JJ_IA(lat|:,lon|:,time|:))
  reg_lhtfl_JJA = regCoef(SMJI,lhtfl_JJA_IA(lat|:,lon|:,time|:))
  
  lhtfl_AM_test = new((/94,192/),"float")
  lhtfl_MJ_test = new((/94,192/),"float")
  lhtfl_JJ_test = new((/94,192/),"float")
  lhtfl_JJA_test = new((/94,192/),"float")
  
  k = 0 
  do j=0,93
  do i=0,191
  lhtfl_AM_test(j,i) = reg_lhtfl_AM@tval(k)
  lhtfl_MJ_test(j,i) = reg_lhtfl_MJ@tval(k)
  lhtfl_JJ_test(j,i) = reg_lhtfl_JJ@tval(k)
  lhtfl_JJA_test(j,i) = reg_lhtfl_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_lhtfl_AM!0   = "lat"
  reg_lhtfl_AM!1   = "lon"
  reg_lhtfl_AM&lat = lhtfl_AM_IA&lat
  reg_lhtfl_AM&lon = lhtfl_AM_IA&lon
  
  copy_VarMeta(reg_lhtfl_AM,reg_lhtfl_MJ)
  copy_VarMeta(reg_lhtfl_AM,reg_lhtfl_JJ)
  copy_VarMeta(reg_lhtfl_AM,reg_lhtfl_JJA)
                    
  copy_VarMeta(reg_lhtfl_AM,lhtfl_AM_test)
  copy_VarMeta(reg_lhtfl_AM,lhtfl_MJ_test)
  copy_VarMeta(reg_lhtfl_AM,lhtfl_JJ_test)
  copy_VarMeta(reg_lhtfl_AM,lhtfl_JJA_test)
  ;===
  f_shtfl   = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/shtfl.sfc.mon.mean.nc","r"); Monthly Mean of Latent Heat Net Flux
  
  shtfl_Mar   = short2flt(f_shtfl->shtfl(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  shtfl_Apr   = short2flt(f_shtfl->shtfl(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  shtfl_May   = short2flt(f_shtfl->shtfl(iStrt_May:iLast_May:12,::-1,:)) 
  
  shtfl_Jun   = short2flt(f_shtfl->shtfl(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  shtfl_Jul   = short2flt(f_shtfl->shtfl(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  shtfl_Aug   = short2flt(f_shtfl->shtfl(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	shtfl_AM   = (shtfl_Apr+shtfl_May)/2. *(-1.) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  shtfl_MJ   = (shtfl_May+shtfl_Jun)/2. *(-1.) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  shtfl_JJ   = (shtfl_Jun+shtfl_Jul)/2. *(-1.) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  shtfl_JJA  = (shtfl_Jun+shtfl_Jul+shtfl_Aug)/3. *(-1.) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month;++++++ 
  
  copy_VarMeta(shtfl_Mar,shtfl_AM)
  copy_VarMeta(shtfl_Mar,shtfl_MJ)
  copy_VarMeta(shtfl_Mar,shtfl_JJ)
  copy_VarMeta(shtfl_Mar,shtfl_JJA)
  
  ;===nine point smoothing
  shtfl_AM   = smth9_Wrap(shtfl_AM, 0.50, 0.25, False) ; light local smoothing
  shtfl_MJ   = smth9_Wrap(shtfl_MJ, 0.50, 0.25, False) ; light local smoothing
  shtfl_JJ   = smth9_Wrap(shtfl_JJ, 0.50, 0.25, False) ; light local smoothing
  shtfl_JJA  = smth9_Wrap(shtfl_JJA, 0.50, 0.25, False) ; light local smoothing
	;===
  shtfl_AM_IA    = wgt_runave_Wrap(shtfl_AM(lat|:,lon|:,time|:),wgt,1)
  shtfl_MJ_IA    = wgt_runave_Wrap(shtfl_MJ(lat|:,lon|:,time|:),wgt,1)
  shtfl_JJ_IA    = wgt_runave_Wrap(shtfl_JJ(lat|:,lon|:,time|:),wgt,1)
  shtfl_JJA_IA    = wgt_runave_Wrap(shtfl_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_shtfl_AM = regCoef(SMJI,shtfl_AM_IA(lat|:,lon|:,time|:))
  reg_shtfl_MJ = regCoef(SMJI,shtfl_MJ_IA(lat|:,lon|:,time|:))
  reg_shtfl_JJ = regCoef(SMJI,shtfl_JJ_IA(lat|:,lon|:,time|:))
  reg_shtfl_JJA = regCoef(SMJI,shtfl_JJA_IA(lat|:,lon|:,time|:))
  
  shtfl_AM_test = new((/94,192/),"float")
  shtfl_MJ_test = new((/94,192/),"float")
  shtfl_JJ_test = new((/94,192/),"float")
  shtfl_JJA_test = new((/94,192/),"float")
  
  k = 0 
  do j=0,93
  do i=0,191
  shtfl_AM_test(j,i) = reg_shtfl_AM@tval(k)
  shtfl_MJ_test(j,i) = reg_shtfl_MJ@tval(k)
  shtfl_JJ_test(j,i) = reg_shtfl_JJ@tval(k)
  shtfl_JJA_test(j,i) = reg_shtfl_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_shtfl_AM!0   = "lat"
  reg_shtfl_AM!1   = "lon"
  reg_shtfl_AM&lat = shtfl_AM_IA&lat
  reg_shtfl_AM&lon = shtfl_AM_IA&lon
  
  copy_VarMeta(reg_shtfl_AM,reg_shtfl_MJ)
  copy_VarMeta(reg_shtfl_AM,reg_shtfl_JJ)
  copy_VarMeta(reg_shtfl_AM,reg_shtfl_JJA)
                    
  copy_VarMeta(reg_shtfl_AM,shtfl_AM_test)
  copy_VarMeta(reg_shtfl_AM,shtfl_MJ_test)
  copy_VarMeta(reg_shtfl_AM,shtfl_JJ_test)
  copy_VarMeta(reg_shtfl_AM,shtfl_JJA_test)
  ;===
  f_dswrf  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/dswrf.sfc.mon.mean.nc","r"); Monthly Mean Downward Solar Radiation Flux at surface
  f_uswrf  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/uswrf.sfc.mon.mean.nc","r"); Monthly Mean Upward Solar Radiation Flux at surface  
  ;===
  dswrf_Mar   = short2flt(f_dswrf->dswrf(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  dswrf_Apr   = short2flt(f_dswrf->dswrf(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  dswrf_May   = short2flt(f_dswrf->dswrf(iStrt_May:iLast_May:12,::-1,:)) 
  
  dswrf_Jun   = short2flt(f_dswrf->dswrf(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  dswrf_Jul   = short2flt(f_dswrf->dswrf(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  dswrf_Aug   = short2flt(f_dswrf->dswrf(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	dswrf_AM   = (dswrf_Apr+dswrf_May)/2. *(-1.)
  dswrf_MJ   = (dswrf_May+dswrf_Jun)/2. *(-1.)
  dswrf_JJ   = (dswrf_Jun+dswrf_Jul)/2. *(-1.)
  dswrf_JJA  = (dswrf_Jun+dswrf_Jul+dswrf_Aug)/3. *(-1.);++++++
  
  copy_VarMeta(dswrf_Mar,dswrf_AM)
  copy_VarMeta(dswrf_Mar,dswrf_MJ)
  copy_VarMeta(dswrf_Mar,dswrf_JJ)
  copy_VarMeta(dswrf_Mar,dswrf_JJA)
  
  uswrf_Mar   = short2flt(f_uswrf->uswrf(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  uswrf_Apr   = short2flt(f_uswrf->uswrf(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  uswrf_May   = short2flt(f_uswrf->uswrf(iStrt_May:iLast_May:12,::-1,:)) 
  
  uswrf_Jun   = short2flt(f_uswrf->uswrf(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  uswrf_Jul   = short2flt(f_uswrf->uswrf(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  uswrf_Aug   = short2flt(f_uswrf->uswrf(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	uswrf_AM   = (uswrf_Apr+uswrf_May)/2. *(-1.)
  uswrf_MJ   = (uswrf_May+uswrf_Jun)/2. *(-1.)
  uswrf_JJ   = (uswrf_Jun+uswrf_Jul)/2. *(-1.)
  uswrf_JJA  = (uswrf_Jun+uswrf_Jul+uswrf_Aug)/3. *(-1.);++++++
  
  copy_VarMeta(uswrf_Mar,uswrf_AM)
  copy_VarMeta(uswrf_Mar,uswrf_MJ)
  copy_VarMeta(uswrf_Mar,uswrf_JJ)
  copy_VarMeta(uswrf_Mar,uswrf_JJA)
  
  ;===net shortwave radiation
  Net_swr_AM  = (uswrf_AM - dswrf_AM)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_swr_MJ  = (uswrf_MJ - dswrf_MJ)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_swr_JJ  = (uswrf_JJ - dswrf_JJ)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_swr_JJA = (uswrf_JJA - dswrf_JJA) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  ;===
  copy_VarMeta(uswrf_Mar,Net_swr_AM )
  copy_VarMeta(uswrf_Mar,Net_swr_MJ )
  copy_VarMeta(uswrf_Mar,Net_swr_JJ )
  copy_VarMeta(uswrf_Mar,Net_swr_JJA)
  
  ;===nine point smoothing
  Net_swr_AM   = smth9_Wrap(Net_swr_AM, 0.50, 0.25, False) ; light local smoothing
  Net_swr_MJ   = smth9_Wrap(Net_swr_MJ, 0.50, 0.25, False) ; light local smoothing
  Net_swr_JJ   = smth9_Wrap(Net_swr_JJ, 0.50, 0.25, False) ; light local smoothing
  Net_swr_JJA  = smth9_Wrap(Net_swr_JJA, 0.50, 0.25, False) ; light local smoothing
  ;===
  Net_swr_AM_IA    = wgt_runave_Wrap(Net_swr_AM(lat|:,lon|:,time|:),wgt,1)
  Net_swr_MJ_IA    = wgt_runave_Wrap(Net_swr_MJ(lat|:,lon|:,time|:),wgt,1)
  Net_swr_JJ_IA    = wgt_runave_Wrap(Net_swr_JJ(lat|:,lon|:,time|:),wgt,1)
  Net_swr_JJA_IA    = wgt_runave_Wrap(Net_swr_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_Net_swr_AM = regCoef(SMJI,Net_swr_AM_IA(lat|:,lon|:,time|:))
  reg_Net_swr_MJ = regCoef(SMJI,Net_swr_MJ_IA(lat|:,lon|:,time|:))
  reg_Net_swr_JJ = regCoef(SMJI,Net_swr_JJ_IA(lat|:,lon|:,time|:))
  reg_Net_swr_JJA = regCoef(SMJI,Net_swr_JJA_IA(lat|:,lon|:,time|:))
  
  Net_swr_AM_test = new((/94,192/),"float")
  Net_swr_MJ_test = new((/94,192/),"float")
  Net_swr_JJ_test = new((/94,192/),"float")
  Net_swr_JJA_test = new((/94,192/),"float")
  
  k = 0 
  do j=0,93
  do i=0,191
  Net_swr_AM_test(j,i) = reg_Net_swr_AM@tval(k)
  Net_swr_MJ_test(j,i) = reg_Net_swr_MJ@tval(k)
  Net_swr_JJ_test(j,i) = reg_Net_swr_JJ@tval(k)
  Net_swr_JJA_test(j,i) = reg_Net_swr_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_Net_swr_AM!0   = "lat"
  reg_Net_swr_AM!1   = "lon"
  reg_Net_swr_AM&lat = Net_swr_AM_IA&lat
  reg_Net_swr_AM&lon = Net_swr_AM_IA&lon
  
  copy_VarMeta(reg_Net_swr_AM,reg_Net_swr_MJ)
  copy_VarMeta(reg_Net_swr_AM,reg_Net_swr_JJ)
  copy_VarMeta(reg_Net_swr_AM,reg_Net_swr_JJA)
                    
  copy_VarMeta(reg_Net_swr_AM,Net_swr_AM_test)
  copy_VarMeta(reg_Net_swr_AM,Net_swr_MJ_test)
  copy_VarMeta(reg_Net_swr_AM,Net_swr_JJ_test)
  copy_VarMeta(reg_Net_swr_AM,Net_swr_JJA_test)
  ;===
  f_dlwrf  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/dlwrf.sfc.mon.mean.nc","r"); Monthly Mean Downward Solar Radiation Flux at surface
  f_ulwrf  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/surface_gauss/ulwrf.sfc.mon.mean.nc","r"); Monthly Mean Upward Solar Radiation Flux at surface  
  ;===
  dlwrf_Mar   = short2flt(f_dlwrf->dlwrf(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  dlwrf_Apr   = short2flt(f_dlwrf->dlwrf(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  dlwrf_May   = short2flt(f_dlwrf->dlwrf(iStrt_May:iLast_May:12,::-1,:)) 
  
  dlwrf_Jun   = short2flt(f_dlwrf->dlwrf(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  dlwrf_Jul   = short2flt(f_dlwrf->dlwrf(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  dlwrf_Aug   = short2flt(f_dlwrf->dlwrf(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	dlwrf_AM   = (dlwrf_Apr+dlwrf_May)/2. *(-1.)
  dlwrf_MJ   = (dlwrf_May+dlwrf_Jun)/2. *(-1.)
  dlwrf_JJ   = (dlwrf_Jun+dlwrf_Jul)/2. *(-1.)
  dlwrf_JJA  = (dlwrf_Jun+dlwrf_Jul+dlwrf_Aug)/3. *(-1.);++++++
  
  copy_VarMeta(dlwrf_Mar,dlwrf_AM)
  copy_VarMeta(dlwrf_Mar,dlwrf_MJ)
  copy_VarMeta(dlwrf_Mar,dlwrf_JJ)
  copy_VarMeta(dlwrf_Mar,dlwrf_JJA)
  
  ulwrf_Mar   = short2flt(f_ulwrf->ulwrf(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  ulwrf_Apr   = short2flt(f_ulwrf->ulwrf(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  ulwrf_May   = short2flt(f_ulwrf->ulwrf(iStrt_May:iLast_May:12,::-1,:)) 
  
  ulwrf_Jun   = short2flt(f_ulwrf->ulwrf(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  ulwrf_Jul   = short2flt(f_ulwrf->ulwrf(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  ulwrf_Aug   = short2flt(f_ulwrf->ulwrf(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	ulwrf_AM   = (ulwrf_Apr+ulwrf_May)/2. *(-1.)
  ulwrf_MJ   = (ulwrf_May+ulwrf_Jun)/2. *(-1.)
  ulwrf_JJ   = (ulwrf_Jun+ulwrf_Jul)/2. *(-1.)
  ulwrf_JJA  = (ulwrf_Jun+ulwrf_Jul+ulwrf_Aug)/3. *(-1.);++++++
  
  copy_VarMeta(ulwrf_Mar,ulwrf_AM)
  copy_VarMeta(ulwrf_Mar,ulwrf_MJ)
  copy_VarMeta(ulwrf_Mar,ulwrf_JJ)
  copy_VarMeta(ulwrf_Mar,ulwrf_JJA)
  
  ;===net longwave radiation
  Net_lwr_AM  = (ulwrf_AM - dlwrf_AM)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_lwr_MJ  = (ulwrf_MJ - dlwrf_MJ)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_lwr_JJ  = (ulwrf_JJ - dlwrf_JJ)   /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  Net_lwr_JJA = (ulwrf_JJA - dlwrf_JJA) /30./3996./1035. * (30.*24.*60.*60.) ;Change unit to degC/month
  ;===
  copy_VarMeta(ulwrf_Mar,Net_lwr_AM )
  copy_VarMeta(ulwrf_Mar,Net_lwr_MJ )
  copy_VarMeta(ulwrf_Mar,Net_lwr_JJ )
  copy_VarMeta(ulwrf_Mar,Net_lwr_JJA)
  ;===nine point smoothing
  Net_lwr_AM   = smth9_Wrap(Net_lwr_AM, 0.50, 0.25, False) ; light local smoothing
  Net_lwr_MJ   = smth9_Wrap(Net_lwr_MJ, 0.50, 0.25, False) ; light local smoothing
  Net_lwr_JJ   = smth9_Wrap(Net_lwr_JJ, 0.50, 0.25, False) ; light local smoothing
  Net_lwr_JJA  = smth9_Wrap(Net_lwr_JJA, 0.50, 0.25, False) ; light local smoothing
  ;===
  Net_lwr_AM_IA    = wgt_runave_Wrap(Net_lwr_AM(lat|:,lon|:,time|:),wgt,1)
  Net_lwr_MJ_IA    = wgt_runave_Wrap(Net_lwr_MJ(lat|:,lon|:,time|:),wgt,1)
  Net_lwr_JJ_IA    = wgt_runave_Wrap(Net_lwr_JJ(lat|:,lon|:,time|:),wgt,1)
  Net_lwr_JJA_IA    = wgt_runave_Wrap(Net_lwr_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_Net_lwr_AM = regCoef(SMJI,Net_lwr_AM_IA(lat|:,lon|:,time|:))
  reg_Net_lwr_MJ = regCoef(SMJI,Net_lwr_MJ_IA(lat|:,lon|:,time|:))
  reg_Net_lwr_JJ = regCoef(SMJI,Net_lwr_JJ_IA(lat|:,lon|:,time|:))
  reg_Net_lwr_JJA = regCoef(SMJI,Net_lwr_JJA_IA(lat|:,lon|:,time|:))
  
  Net_lwr_AM_test = new((/94,192/),"float")
  Net_lwr_MJ_test = new((/94,192/),"float")
  Net_lwr_JJ_test = new((/94,192/),"float")
  Net_lwr_JJA_test = new((/94,192/),"float")
  
  k = 0 
  do j=0,93
  do i=0,191
  Net_lwr_AM_test(j,i) = reg_Net_lwr_AM@tval(k)
  Net_lwr_MJ_test(j,i) = reg_Net_lwr_MJ@tval(k)
  Net_lwr_JJ_test(j,i) = reg_Net_lwr_JJ@tval(k)
  Net_lwr_JJA_test(j,i) = reg_Net_lwr_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_Net_lwr_AM!0   = "lat"
  reg_Net_lwr_AM!1   = "lon"
  reg_Net_lwr_AM&lat = Net_lwr_AM_IA&lat
  reg_Net_lwr_AM&lon = Net_lwr_AM_IA&lon
  
  copy_VarMeta(reg_Net_lwr_AM,reg_Net_lwr_MJ)
  copy_VarMeta(reg_Net_lwr_AM,reg_Net_lwr_JJ)
  copy_VarMeta(reg_Net_lwr_AM,reg_Net_lwr_JJA)
                    
  copy_VarMeta(reg_Net_lwr_AM,Net_lwr_AM_test)
  copy_VarMeta(reg_Net_lwr_AM,Net_lwr_MJ_test)
  copy_VarMeta(reg_Net_lwr_AM,Net_lwr_JJ_test)
  copy_VarMeta(reg_Net_lwr_AM,Net_lwr_JJA_test)
  ;===
  ;===net heat flux
  NHF_AM  = (lhtfl_AM+shtfl_AM+Net_swr_AM+Net_lwr_AM)
  NHF_MJ  = (lhtfl_MJ+shtfl_MJ+Net_swr_MJ+Net_lwr_MJ)
  NHF_JJ  = (lhtfl_JJ+shtfl_JJ+Net_swr_JJ+Net_lwr_JJ)
  NHF_JJA = (lhtfl_JJA+shtfl_JJA+Net_swr_JJA+Net_lwr_JJA)
  
  copy_VarMeta(uswrf_Mar,NHF_AM)
  copy_VarMeta(uswrf_Mar,NHF_MJ)
  copy_VarMeta(uswrf_Mar,NHF_JJ)
  copy_VarMeta(uswrf_Mar,NHF_JJA)
  
  ;===
  NHF_AM_IA    = wgt_runave_Wrap(NHF_AM(lat|:,lon|:,time|:),wgt,1)
  NHF_MJ_IA    = wgt_runave_Wrap(NHF_MJ(lat|:,lon|:,time|:),wgt,1)
  NHF_JJ_IA    = wgt_runave_Wrap(NHF_JJ(lat|:,lon|:,time|:),wgt,1)
  NHF_JJA_IA    = wgt_runave_Wrap(NHF_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_NHF_AM  = regCoef(SMJI,NHF_AM_IA(lat|:,lon|:,time|:))
  reg_NHF_MJ  = regCoef(SMJI,NHF_MJ_IA(lat|:,lon|:,time|:))
  reg_NHF_JJ  = regCoef(SMJI,NHF_JJ_IA(lat|:,lon|:,time|:))
  reg_NHF_JJA = regCoef(SMJI,NHF_JJA_IA(lat|:,lon|:,time|:))
  
  NHF_AM_test  = new((/94,192/),"float")
  NHF_MJ_test  = new((/94,192/),"float")
  NHF_JJ_test  = new((/94,192/),"float")
  NHF_JJA_test = new((/94,192/),"float")
  
  k = 0 
  do j=0,93
  do i=0,191
  NHF_AM_test(j,i) = reg_NHF_AM@tval(k)
  NHF_MJ_test(j,i) = reg_NHF_MJ@tval(k)
  NHF_JJ_test(j,i) = reg_NHF_JJ@tval(k)
  NHF_JJA_test(j,i) = reg_NHF_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_NHF_AM!0   = "lat"
  reg_NHF_AM!1   = "lon"
  reg_NHF_AM&lat = NHF_AM_IA&lat
  reg_NHF_AM&lon = NHF_AM_IA&lon
  
  copy_VarMeta(reg_NHF_AM,reg_NHF_MJ)
  copy_VarMeta(reg_NHF_AM,reg_NHF_JJ)
  copy_VarMeta(reg_NHF_AM,reg_NHF_JJA)
                    
  copy_VarMeta(reg_NHF_AM,NHF_AM_test)
  copy_VarMeta(reg_NHF_AM,NHF_MJ_test)
  copy_VarMeta(reg_NHF_AM,NHF_JJ_test)
  copy_VarMeta(reg_NHF_AM,NHF_JJA_test)
  ;===
  fmsk   =  addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")    
  lsmask = landsea_mask(fmsk->LSMASK,reg_NHF_AM&lat,reg_NHF_AM&lon)
   
  reg_NHF_AM    = mask(reg_NHF_AM ,lsmask.eq.1,False) 
  reg_NHF_MJ    = mask(reg_NHF_MJ ,lsmask.eq.1,False) 
  reg_NHF_JJ    = mask(reg_NHF_JJ ,lsmask.eq.1,False) 
  reg_NHF_JJA   = mask(reg_NHF_JJA,lsmask.eq.1,False) 
  NHF_AM_test        = mask(NHF_AM_test ,lsmask.eq.1,False) 
  NHF_MJ_test        = mask(NHF_MJ_test ,lsmask.eq.1,False) 
  NHF_JJ_test        = mask(NHF_JJ_test ,lsmask.eq.1,False) 
  NHF_JJA_test       = mask(NHF_JJA_test,lsmask.eq.1,False) 
  
  reg_lhtfl_AM    = mask(reg_lhtfl_AM ,lsmask.eq.1,False) 
  reg_lhtfl_MJ    = mask(reg_lhtfl_MJ ,lsmask.eq.1,False) 
  reg_lhtfl_JJ    = mask(reg_lhtfl_JJ ,lsmask.eq.1,False) 
  reg_lhtfl_JJA   = mask(reg_lhtfl_JJA,lsmask.eq.1,False) 
  lhtfl_AM_test        = mask(lhtfl_AM_test ,lsmask.eq.1,False) 
  lhtfl_MJ_test        = mask(lhtfl_MJ_test ,lsmask.eq.1,False) 
  lhtfl_JJ_test        = mask(lhtfl_JJ_test ,lsmask.eq.1,False) 
  lhtfl_JJA_test       = mask(lhtfl_JJA_test,lsmask.eq.1,False) 
  
  reg_shtfl_AM    = mask(reg_shtfl_AM ,lsmask.eq.1,False) 
  reg_shtfl_MJ    = mask(reg_shtfl_MJ ,lsmask.eq.1,False) 
  reg_shtfl_JJ    = mask(reg_shtfl_JJ ,lsmask.eq.1,False) 
  reg_shtfl_JJA   = mask(reg_shtfl_JJA,lsmask.eq.1,False) 
  shtfl_AM_test        = mask(shtfl_AM_test ,lsmask.eq.1,False) 
  shtfl_MJ_test        = mask(shtfl_MJ_test ,lsmask.eq.1,False) 
  shtfl_JJ_test        = mask(shtfl_JJ_test ,lsmask.eq.1,False) 
  shtfl_JJA_test       = mask(shtfl_JJA_test,lsmask.eq.1,False) 
  
  reg_Net_swr_AM    = mask(reg_Net_swr_AM ,lsmask.eq.1,False) 
  reg_Net_swr_MJ    = mask(reg_Net_swr_MJ ,lsmask.eq.1,False) 
  reg_Net_swr_JJ    = mask(reg_Net_swr_JJ ,lsmask.eq.1,False) 
  reg_Net_swr_JJA   = mask(reg_Net_swr_JJA,lsmask.eq.1,False) 
  Net_swr_AM_test        = mask(Net_swr_AM_test ,lsmask.eq.1,False) 
  Net_swr_MJ_test        = mask(Net_swr_MJ_test ,lsmask.eq.1,False) 
  Net_swr_JJ_test        = mask(Net_swr_JJ_test ,lsmask.eq.1,False) 
  Net_swr_JJA_test       = mask(Net_swr_JJA_test,lsmask.eq.1,False) 
  
  reg_Net_lwr_AM    = mask(reg_Net_lwr_AM ,lsmask.eq.1,False) 
  reg_Net_lwr_MJ    = mask(reg_Net_lwr_MJ ,lsmask.eq.1,False) 
  reg_Net_lwr_JJ    = mask(reg_Net_lwr_JJ ,lsmask.eq.1,False) 
  reg_Net_lwr_JJA   = mask(reg_Net_lwr_JJA,lsmask.eq.1,False) 
  Net_lwr_AM_test        = mask(Net_lwr_AM_test ,lsmask.eq.1,False) 
  Net_lwr_MJ_test        = mask(Net_lwr_MJ_test ,lsmask.eq.1,False) 
  Net_lwr_JJ_test        = mask(Net_lwr_JJ_test ,lsmask.eq.1,False) 
  Net_lwr_JJA_test       = mask(Net_lwr_JJA_test,lsmask.eq.1,False) 
  ;===wind speed
  f_uwnd  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/uwnd.mon.mean.nc", "r")
  f_vwnd  = addfile("~/NCEP-2.5/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc", "r")
  
  uwnd_Mar   = short2flt(f_uwnd->uwnd(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  uwnd_Apr   = short2flt(f_uwnd->uwnd(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  uwnd_May   = short2flt(f_uwnd->uwnd(iStrt_May:iLast_May:12,{925},::-1,:)) 
  uwnd_Jun   = short2flt(f_uwnd->uwnd(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  uwnd_Jul   = short2flt(f_uwnd->uwnd(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  uwnd_Aug   = short2flt(f_uwnd->uwnd(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	uwnd_AM   = (uwnd_Apr+uwnd_May)/2.
  uwnd_MJ   = (uwnd_May+uwnd_Jun)/2.
  uwnd_JJ   = (uwnd_Jun+uwnd_Jul)/2.
  uwnd_JJA   = (uwnd_Jun+uwnd_Jul+uwnd_Aug)/3.
  
  vwnd_Mar   = short2flt(f_vwnd->vwnd(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  vwnd_Apr   = short2flt(f_vwnd->vwnd(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  vwnd_May   = short2flt(f_vwnd->vwnd(iStrt_May:iLast_May:12,{925},::-1,:)) 
  vwnd_Jun   = short2flt(f_vwnd->vwnd(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  vwnd_Jul   = short2flt(f_vwnd->vwnd(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  vwnd_Aug   = short2flt(f_vwnd->vwnd(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	vwnd_AM   = (vwnd_Apr+vwnd_May)/2.
  vwnd_MJ   = (vwnd_May+vwnd_Jun)/2.
  vwnd_JJA   = (vwnd_Jun+vwnd_Jul+vwnd_Aug)/3.
  vwnd_JJ   = (vwnd_Jun+vwnd_Jul)/2.
  
  wspd_AM = sqrt(uwnd_AM^2+vwnd_AM^2)
  wspd_MJ = sqrt(uwnd_MJ^2+vwnd_MJ^2)
  wspd_JJ = sqrt(uwnd_JJ^2+vwnd_JJ^2)
  wspd_JJA = sqrt(uwnd_JJA^2+vwnd_JJA^2)
  
  copy_VarMeta(uwnd_Mar,wspd_AM) ;###
  copy_VarMeta(uwnd_Mar,wspd_MJ) ;###
  copy_VarMeta(uwnd_Mar,wspd_JJ) ;###
  copy_VarMeta(uwnd_Mar,wspd_JJA) ;###
	;===
  wspd_AM_IA     = wgt_runave_Wrap(wspd_AM(lat|:,lon|:,time|:),wgt,1)
  wspd_MJ_IA     = wgt_runave_Wrap(wspd_MJ(lat|:,lon|:,time|:),wgt,1)
  wspd_JJ_IA     = wgt_runave_Wrap(wspd_JJ(lat|:,lon|:,time|:),wgt,1)
  wspd_JJA_IA    = wgt_runave_Wrap(wspd_JJA(lat|:,lon|:,time|:),wgt,1)  
  ;===regression
  reg_wspd_AM = regCoef(SMJI,wspd_AM_IA(lat|:,lon|:,time|:))
  reg_wspd_MJ = regCoef(SMJI,wspd_MJ_IA(lat|:,lon|:,time|:))
  reg_wspd_JJ = regCoef(SMJI,wspd_JJ_IA(lat|:,lon|:,time|:))
  reg_wspd_JJA = regCoef(SMJI,wspd_JJA_IA(lat|:,lon|:,time|:))
  
  wspd_AM_test = new((/73,144/),"float")
  wspd_MJ_test = new((/73,144/),"float")
  wspd_JJ_test = new((/73,144/),"float")
  wspd_JJA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  wspd_AM_test(j,i) = reg_wspd_AM@tval(k)
  wspd_MJ_test(j,i) = reg_wspd_MJ@tval(k)
  wspd_JJ_test(j,i) = reg_wspd_JJ@tval(k)
  wspd_JJA_test(j,i) = reg_wspd_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_wspd_AM!0   = "lat"
  reg_wspd_AM!1   = "lon"
  reg_wspd_AM&lat = wspd_AM_IA&lat
  reg_wspd_AM&lon = wspd_AM_IA&lon
  
  copy_VarMeta(reg_wspd_AM,reg_wspd_MJ)
  copy_VarMeta(reg_wspd_AM,reg_wspd_JJ)
  copy_VarMeta(reg_wspd_AM,reg_wspd_JJA)
                    
  copy_VarMeta(reg_wspd_AM,wspd_AM_test)
  copy_VarMeta(reg_wspd_AM,wspd_MJ_test)
  copy_VarMeta(reg_wspd_AM,wspd_JJ_test)
  copy_VarMeta(reg_wspd_AM,wspd_JJA_test)
  ;===
  delete(fmsk)
  delete(lsmask)
  
  fmsk   =  addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")    
  lsmask = landsea_mask(fmsk->LSMASK,reg_wspd_AM&lat,reg_wspd_AM&lon)
  
  reg_wspd_AM    = mask(reg_wspd_AM ,lsmask.eq.1,False) 
  reg_wspd_MJ    = mask(reg_wspd_MJ ,lsmask.eq.1,False) 
  reg_wspd_JJ    = mask(reg_wspd_JJ ,lsmask.eq.1,False) 
  reg_wspd_JJA   = mask(reg_wspd_JJA,lsmask.eq.1,False) 
  wspd_AM_test        = mask(wspd_AM_test ,lsmask.eq.1,False) 
  wspd_MJ_test        = mask(wspd_MJ_test ,lsmask.eq.1,False) 
  wspd_JJ_test        = mask(wspd_JJ_test ,lsmask.eq.1,False) 
  wspd_JJA_test       = mask(wspd_JJA_test,lsmask.eq.1,False)
  
  ;===sst tendency
  f_sst  =  addfile("/disk_raid6/data/noaa.ersst.v5/sst.mnmean.nc", "r")
	
	delete(TIME)
	delete(YMD)
	
	TIME    = f_sst->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
  ymdStrt_Nov = 19781101                         ; start yyyymmdd
  ymdLast_Nov = 20181101                          
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
  ymdStrt_Sep = 19790901                         ; start yyyymmdd
  ymdLast_Sep = 20190901                         ; last  
  
  iStrt_Nov   = ind(YMD.eq.ymdStrt_Nov)              ; index start
  iLast_Nov   = ind(YMD.eq.ymdLast_Nov)              ; index last 
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
  iStrt_Sep   = ind(YMD.eq.ymdStrt_Sep)              ; index start
  iLast_Sep   = ind(YMD.eq.ymdLast_Sep)              ; index last 
  
  sst_Nov   = short2flt(f_sst->sst(iStrt_Nov:iLast_Nov:12,::-1,:)) 
  sst_Dec   = short2flt(f_sst->sst(iStrt_Dec:iLast_Dec:12,::-1,:)) 
  sst_Jan   = short2flt(f_sst->sst(iStrt_Jan:iLast_Jan:12,::-1,:)) 
  sst_Feb   = short2flt(f_sst->sst(iStrt_Feb:iLast_Feb:12,::-1,:)) 
  
  sst_Mar   = short2flt(f_sst->sst(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  sst_Apr   = short2flt(f_sst->sst(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  sst_May   = short2flt(f_sst->sst(iStrt_May:iLast_May:12,::-1,:)) 
  
  sst_Jun   = short2flt(f_sst->sst(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  sst_Jul   = short2flt(f_sst->sst(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  sst_Aug   = short2flt(f_sst->sst(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  sst_Sep   = short2flt(f_sst->sst(iStrt_Sep:iLast_Sep:12,::-1,:)) 
  
  ;===sst tendency
  sst_AM   = ( (sst_May-sst_Mar)/2. + (sst_Jun-sst_Apr)/2.)/2. ;sst tendency
  sst_MJ   = ( (sst_Jun-sst_Apr)/2. + (sst_Jul-sst_May)/2.)/2. ;sst tendency
  sst_JJ   = ( (sst_Jul-sst_May)/2. + (sst_Aug-sst_Jun)/2.)/2. ;sst tendency
  sst_JJA  = ( (sst_Jul-sst_May)/2. + (sst_Aug-sst_Jun)/2. + (sst_Sep-sst_Jul)/2.)/3. ;sst tendency
  
  copy_VarMeta(sst_Mar,sst_AM)
  copy_VarMeta(sst_Mar,sst_MJ)
  copy_VarMeta(sst_Mar,sst_JJ)
  copy_VarMeta(sst_Mar,sst_JJA)
  
  sst_AM_IA    = wgt_runave_Wrap(sst_AM(lat|:,lon|:,time|:),wgt,1)
  sst_MJ_IA    = wgt_runave_Wrap(sst_MJ(lat|:,lon|:,time|:),wgt,1)
  sst_JJ_IA    = wgt_runave_Wrap(sst_JJ(lat|:,lon|:,time|:),wgt,1)
  sst_JJA_IA   = wgt_runave_Wrap(sst_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_sst_AM = regCoef(SMJI,sst_AM_IA(lat|:,lon|:,time|:))
  reg_sst_MJ = regCoef(SMJI,sst_MJ_IA(lat|:,lon|:,time|:))
  reg_sst_JJ = regCoef(SMJI,sst_JJ_IA(lat|:,lon|:,time|:))
  reg_sst_JJA = regCoef(SMJI,sst_JJA_IA(lat|:,lon|:,time|:))
  
  sst_AM_test = new((/89,180/),"float")
  sst_MJ_test = new((/89,180/),"float")
  sst_JJ_test = new((/89,180/),"float")
  sst_JJA_test = new((/89,180/),"float")
  
  k = 0 
  do j=0,88
  do i=0,179
  sst_AM_test(j,i) = reg_sst_AM@tval(k)
  sst_MJ_test(j,i) = reg_sst_MJ@tval(k)
  sst_JJ_test(j,i) = reg_sst_JJ@tval(k)
  sst_JJA_test(j,i) = reg_sst_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_sst_AM!0   = "lat"
  reg_sst_AM!1   = "lon"
  reg_sst_AM&lat = sst_AM_IA&lat
  reg_sst_AM&lon = sst_AM_IA&lon
  
  copy_VarMeta(reg_sst_AM,reg_sst_MJ)
  copy_VarMeta(reg_sst_AM,reg_sst_JJ)
  copy_VarMeta(reg_sst_AM,reg_sst_JJA)
                    
  copy_VarMeta(reg_sst_AM,sst_AM_test)
  copy_VarMeta(reg_sst_AM,sst_MJ_test)
  copy_VarMeta(reg_sst_AM,sst_JJ_test)
  copy_VarMeta(reg_sst_AM,sst_JJA_test)
  ;===
	;===
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.6")
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
  res@mpMaxLatF    = 30.
  res@mpMinLatF    = -20
  res@mpMaxLonF    = 120
  res@mpMinLonF    = 20
  
  res@tmXBLabelFontHeightF = 0.028
  res@tmYLLabelFontHeightF = 0.028
  
  res@gsnLeftStringFontHeightF   = 0.03
  res@gsnRightStringFontHeightF  = 0.03
  res@gsnCenterStringFontHeightF = 0.03
  
  res@gsnCenterStringOrthogonalPosF = 0.025
  res@gsnRightStringOrthogonalPosF  = 0.015
  res@gsnLeftStringOrthogonalPosF   = 0.015
  
  res@lbLabelBarOn             = False
  
  plot = new(6,graphic)
  
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = (/-0.14,-0.11,-0.08,-0.05,-0.02,0.02,0.05,0.08,0.11,0.14/)
  
  res@tmXBOn               = False
  res@gsnCenterString             = "SST tendency"
  res@gsnLeftString               = "(a) "
  res@gsnRightString              = "AM"
  plot(0)  = gsn_csm_contour_map_ce(wks,reg_sst_AM,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(b) "
  res@gsnRightString              = "MJ"
  plot(2)  = gsn_csm_contour_map_ce(wks,reg_sst_MJ,res)
  
  res@gsnCenterString             = ""
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.2
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  res@gsnLeftString               = "(c) "
  res@gsnRightString              = "JJA"
  res@tmXBOn               = True
  plot(4)  = gsn_csm_contour_map_ce(wks,reg_sst_JJA,res)
  
  res@tmYLOn               = False
  res@tmXBOn               = False

  res@lbLabelBarOn                = False
  res@gsnCenterString             = "Net Heat Flux"
  res@gsnLeftString               = "(d) "
  res@gsnRightString              = "AM"
  plot(1)  = gsn_csm_contour_map_ce(wks,reg_NHF_AM,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(e) "
  res@gsnRightString              = "MJ"
  plot(3)  = gsn_csm_contour_map_ce(wks,reg_NHF_MJ,res)
  res@gsnLeftString               = "(f) "
  res@gsnRightString              = "JJA"
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
  plot(5)  = gsn_csm_contour_map_ce(wks,reg_NHF_JJA,res)
  
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
  
  plot_test     = new(6,graphic)
  
  plot_test(0)   = gsn_csm_contour(wks,sst_AM_test    ,res_t)
  plot_test(2)   = gsn_csm_contour(wks,sst_MJ_test    ,res_t)
  plot_test(4)   = gsn_csm_contour(wks,sst_JJA_test    ,res_t)
  
  plot_test(1)   = gsn_csm_contour(wks,NHF_AM_test    ,res_t)
  plot_test(3)   = gsn_csm_contour(wks,NHF_MJ_test    ,res_t)
  plot_test(5)   = gsn_csm_contour(wks,NHF_JJA_test    ,res_t)
  
  ;===
  do i=0,5
  overlay(plot(i),plot_test(i))
  end do
  ;===
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
  
  resp@gsnPanelYF               = (/0.95,0.95,0.7,0.7,0.45,0.45/);-0.05
  resp@gsnPanelXF               = (/0.1,0.54,0.1,0.54,0.1,0.54/);-0.05
  
  gsn_panel(wks,plot(0:5),(/3,2/),resp)
  
end