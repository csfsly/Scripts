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
	
  f_sst  =  addfile("/disk_raid6/data/noaa.ersst.v5/sst.mnmean.nc", "r")
	
	TIME    = f_sst->time                       ; days since ...
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
  
  sst_Dec   = short2flt(f_sst->sst(iStrt_Dec:iLast_Dec:12,::-1,:)) 
  sst_Jan   = short2flt(f_sst->sst(iStrt_Jan:iLast_Jan:12,::-1,:)) 
  sst_Feb   = short2flt(f_sst->sst(iStrt_Feb:iLast_Feb:12,::-1,:)) 
  
  sst_Mar   = short2flt(f_sst->sst(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  sst_Apr   = short2flt(f_sst->sst(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  sst_May   = short2flt(f_sst->sst(iStrt_May:iLast_May:12,::-1,:)) 
  
  sst_Jun   = short2flt(f_sst->sst(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  sst_Jul   = short2flt(f_sst->sst(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  sst_Aug   = short2flt(f_sst->sst(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	sst_AM   = (sst_Apr+sst_May)/2.
  sst_MJ   = (sst_May+sst_Jun)/2.
  sst_JJ   = (sst_Jun+sst_Jul)/2.
  sst_JJA  = (sst_Jun+sst_Jul+sst_Aug)/3.
    
  copy_VarMeta(sst_Mar,sst_AM)
  copy_VarMeta(sst_Mar,sst_MJ)
  copy_VarMeta(sst_Mar,sst_JJ)
  copy_VarMeta(sst_Mar,sst_JJA)
	;===
	;===高通滤波
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  sst_AM_IA    = wgt_runave_Wrap(sst_AM(lat|:,lon|:,time|:),wgt,1)
  sst_MJ_IA    = wgt_runave_Wrap(sst_MJ(lat|:,lon|:,time|:),wgt,1)
  sst_JJ_IA    = wgt_runave_Wrap(sst_JJ(lat|:,lon|:,time|:),wgt,1)
  sst_JJA_IA   = wgt_runave_Wrap(sst_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===define WIOSST index
  WIOSSTI    = dim_avg_n_Wrap(dim_avg_n_Wrap(sst_JJA_IA({-10:20},{40:90},:),0),0)
  WIOSSTI    = dim_standardize(WIOSSTI,1)
  
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
  f_precip  =  addfile("~/rainfall-grid/gpcp/precip.mon.mean.nc", "r")
	
	delete(TIME)
	delete(YMD)
	
	TIME    = f_precip->time                       ; days since ...
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
  
  precip_Mar   = short2flt(f_precip->precip(iStrt_Mar:iLast_Mar:12,::-1,:)) 
  precip_Apr   = short2flt(f_precip->precip(iStrt_Apr:iLast_Apr:12,::-1,:)) 
  precip_May   = short2flt(f_precip->precip(iStrt_May:iLast_May:12,::-1,:)) 
  
  precip_Jun   = short2flt(f_precip->precip(iStrt_Jun:iLast_Jun:12,::-1,:)) 
  precip_Jul   = short2flt(f_precip->precip(iStrt_Jul:iLast_Jul:12,::-1,:)) 
  precip_Aug   = short2flt(f_precip->precip(iStrt_Aug:iLast_Aug:12,::-1,:)) 
  
	precip_AM   = (precip_Apr+precip_May)/2.
  precip_MJ   = (precip_May+precip_Jun)/2.
  precip_JJ   = (precip_Jun+precip_Jul)/2.
  precip_JJA  = (precip_Jun+precip_Jul+precip_Aug)/3.
    
  copy_VarMeta(precip_Mar,precip_AM)
  copy_VarMeta(precip_Mar,precip_MJ)
  copy_VarMeta(precip_Mar,precip_JJ)
  copy_VarMeta(precip_Mar,precip_JJA)
  
  precip_AM_IA    = wgt_runave_Wrap(precip_AM(lat|:,lon|:,time|:),wgt,1)
  precip_MJ_IA    = wgt_runave_Wrap(precip_MJ(lat|:,lon|:,time|:),wgt,1)
  precip_JJ_IA    = wgt_runave_Wrap(precip_JJ(lat|:,lon|:,time|:),wgt,1)
  precip_JJA_IA   = wgt_runave_Wrap(precip_JJA(lat|:,lon|:,time|:),wgt,1)
  
  ;===regression
  reg_precip_AM  = regCoef(SMJI,precip_AM_IA(lat|:,lon|:,time|:))
  reg_precip_MJ  = regCoef(SMJI,precip_MJ_IA(lat|:,lon|:,time|:))
  reg_precip_JJ  = regCoef(SMJI,precip_JJ_IA(lat|:,lon|:,time|:))
  reg_precip_JJA = regCoef(SMJI,precip_JJA_IA(lat|:,lon|:,time|:))
  
  precip_AM_test  = new((/72,144/),"float")
  precip_MJ_test  = new((/72,144/),"float")
  precip_JJ_test  = new((/72,144/),"float")
  precip_JJA_test = new((/72,144/),"float")
  
  k = 0 
  do j=0,71
  do i=0,143
  precip_AM_test(j,i) = reg_precip_AM@tval(k)
  precip_MJ_test(j,i) = reg_precip_MJ@tval(k)
  precip_JJ_test(j,i) = reg_precip_JJ@tval(k)
  precip_JJA_test(j,i) = reg_precip_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_precip_AM!0   = "lat"
  reg_precip_AM!1   = "lon"
  reg_precip_AM&lat = precip_AM_IA&lat
  reg_precip_AM&lon = precip_AM_IA&lon
  
  copy_VarMeta(reg_precip_AM,reg_precip_MJ)
  copy_VarMeta(reg_precip_AM,reg_precip_JJ)
  copy_VarMeta(reg_precip_AM,reg_precip_JJA)
                    
  copy_VarMeta(reg_precip_AM,precip_AM_test)
  copy_VarMeta(reg_precip_AM,precip_MJ_test)
  copy_VarMeta(reg_precip_AM,precip_JJ_test)
  copy_VarMeta(reg_precip_AM,precip_JJA_test)
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
  
  uwnd_Mar   = short2flt(f_uwnd->uwnd(iStrt_Mar:iLast_Mar:12,{850},::-1,:)) 
  uwnd_Apr   = short2flt(f_uwnd->uwnd(iStrt_Apr:iLast_Apr:12,{850},::-1,:)) 
  uwnd_May   = short2flt(f_uwnd->uwnd(iStrt_May:iLast_May:12,{850},::-1,:)) 
  uwnd_Jun   = short2flt(f_uwnd->uwnd(iStrt_Jun:iLast_Jun:12,{850},::-1,:)) 
  uwnd_Jul   = short2flt(f_uwnd->uwnd(iStrt_Jul:iLast_Jul:12,{850},::-1,:)) 
  uwnd_Aug   = short2flt(f_uwnd->uwnd(iStrt_Aug:iLast_Aug:12,{850},::-1,:)) 
  
	uwnd_AM   = (uwnd_Apr+uwnd_May)/2.
  uwnd_MJ   = (uwnd_May+uwnd_Jun)/2.
  uwnd_JJ   = (uwnd_Jun+uwnd_Jul)/2.
  uwnd_JJA   = (uwnd_Jun+uwnd_Jul+uwnd_Aug)/3.
  
  copy_VarMeta(uwnd_Mar,uwnd_AM)
  copy_VarMeta(uwnd_Mar,uwnd_MJ)
  copy_VarMeta(uwnd_Mar,uwnd_JJA)
  copy_VarMeta(uwnd_Mar,uwnd_JJ)
	;===
  uwnd_AM_IA    = wgt_runave_Wrap(uwnd_AM(lat|:,lon|:,time|:),wgt,1)
  uwnd_MJ_IA    = wgt_runave_Wrap(uwnd_MJ(lat|:,lon|:,time|:),wgt,1)
  uwnd_JJA_IA    = wgt_runave_Wrap(uwnd_JJA(lat|:,lon|:,time|:),wgt,1)
  uwnd_JJ_IA    = wgt_runave_Wrap(uwnd_JJ(lat|:,lon|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_uwnd_MJ  = regCoef(SMJI,uwnd_MJ_IA(lat|:,lon|:,time|:))
  reg_uwnd_AM  = regCoef(SMJI,uwnd_AM_IA(lat|:,lon|:,time|:))
  reg_uwnd_JJ  = regCoef(SMJI,uwnd_JJ_IA(lat|:,lon|:,time|:))
  reg_uwnd_JJA = regCoef(SMJI,uwnd_JJA_IA(lat|:,lon|:,time|:))
  
  uwnd_AM_test = new((/73,144/),"float")
  uwnd_MJ_test = new((/73,144/),"float")
  uwnd_JJA_test = new((/73,144/),"float")
  uwnd_JJ_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  uwnd_AM_test(j,i) = reg_uwnd_AM@tval(k)
  uwnd_MJ_test(j,i) = reg_uwnd_MJ@tval(k)
  uwnd_JJA_test(j,i) = reg_uwnd_JJA@tval(k)
  uwnd_JJ_test(j,i) = reg_uwnd_JJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_uwnd_AM!0   = "lat"
  reg_uwnd_AM!1   = "lon"
  reg_uwnd_AM&lat = uwnd_AM_IA&lat
  reg_uwnd_AM&lon = uwnd_AM_IA&lon
  
  copy_VarMeta(reg_uwnd_AM,reg_uwnd_MJ)
  copy_VarMeta(reg_uwnd_AM,reg_uwnd_JJA)
  copy_VarMeta(reg_uwnd_AM,reg_uwnd_JJ)
  
  copy_VarMeta(reg_uwnd_AM,uwnd_AM_test)
  copy_VarMeta(reg_uwnd_AM,uwnd_MJ_test)
  copy_VarMeta(reg_uwnd_AM,uwnd_JJA_test)
  copy_VarMeta(reg_uwnd_AM,uwnd_JJ_test)
  ;===
  
  vwnd_Mar   = short2flt(f_vwnd->vwnd(iStrt_Mar:iLast_Mar:12,{850},::-1,:)) 
  vwnd_Apr   = short2flt(f_vwnd->vwnd(iStrt_Apr:iLast_Apr:12,{850},::-1,:)) 
  vwnd_May   = short2flt(f_vwnd->vwnd(iStrt_May:iLast_May:12,{850},::-1,:)) 
  vwnd_Jun   = short2flt(f_vwnd->vwnd(iStrt_Jun:iLast_Jun:12,{850},::-1,:)) 
  vwnd_Jul   = short2flt(f_vwnd->vwnd(iStrt_Jul:iLast_Jul:12,{850},::-1,:)) 
  vwnd_Aug   = short2flt(f_vwnd->vwnd(iStrt_Aug:iLast_Aug:12,{850},::-1,:)) 
  
	vwnd_AM   = (vwnd_Apr+vwnd_May)/2.
  vwnd_MJ   = (vwnd_May+vwnd_Jun)/2.
  vwnd_JJA   = (vwnd_Jun+vwnd_Jul+vwnd_Aug)/3.
  vwnd_JJ   = (vwnd_Jun+vwnd_Jul)/2.
  
  copy_VarMeta(vwnd_Mar,vwnd_AM)
  copy_VarMeta(vwnd_Mar,vwnd_MJ)
  copy_VarMeta(vwnd_Mar,vwnd_JJA)
  copy_VarMeta(vwnd_Mar,vwnd_JJ)
	;===
  vwnd_AM_IA    = wgt_runave_Wrap(vwnd_AM(lat|:,lon|:,time|:),wgt,1)
  vwnd_MJ_IA    = wgt_runave_Wrap(vwnd_MJ(lat|:,lon|:,time|:),wgt,1)
  vwnd_JJA_IA    = wgt_runave_Wrap(vwnd_JJA(lat|:,lon|:,time|:),wgt,1)
  vwnd_JJ_IA    = wgt_runave_Wrap(vwnd_JJ(lat|:,lon|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_vwnd_AM  = regCoef(SMJI,vwnd_AM_IA(lat|:,lon|:,time|:))
  reg_vwnd_MJ  = regCoef(SMJI,vwnd_MJ_IA(lat|:,lon|:,time|:))
  reg_vwnd_JJA = regCoef(SMJI,vwnd_JJA_IA(lat|:,lon|:,time|:))
  reg_vwnd_JJ  = regCoef(SMJI,vwnd_JJ_IA(lat|:,lon|:,time|:))
  
  vwnd_AM_test = new((/73,144/),"float")
  vwnd_MJ_test = new((/73,144/),"float")
  vwnd_JJA_test = new((/73,144/),"float")
  vwnd_JJ_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  vwnd_AM_test(j,i) = reg_vwnd_AM@tval(k)
  vwnd_MJ_test(j,i) = reg_vwnd_MJ@tval(k)
  vwnd_JJA_test(j,i) = reg_vwnd_JJA@tval(k)
  vwnd_JJ_test(j,i) = reg_vwnd_JJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_vwnd_AM!0   = "lat"
  reg_vwnd_AM!1   = "lon"
  reg_vwnd_AM&lat = vwnd_AM_IA&lat
  reg_vwnd_AM&lon = vwnd_AM_IA&lon
  
  copy_VarMeta(reg_vwnd_AM,reg_vwnd_MJ)
  copy_VarMeta(reg_vwnd_AM,reg_vwnd_JJA)
  copy_VarMeta(reg_vwnd_AM,reg_vwnd_JJ)
  
  copy_VarMeta(reg_vwnd_AM,vwnd_AM_test)
  copy_VarMeta(reg_vwnd_AM,vwnd_MJ_test)
  copy_VarMeta(reg_vwnd_AM,vwnd_JJA_test)
  copy_VarMeta(reg_vwnd_AM,vwnd_JJ_test)
  ;===
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.5")
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
  
  plot = new(8,graphic)
  
  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = (/-0.17,-0.14,-0.11,-0.08,-0.05,-0.02,0.02,0.05,0.08,0.11,0.14,0.17/)
  
  res@tmXBOn               = False
  res@gsnCenterString             = "SST&winds(850hPa)"
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
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.2
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  res@gsnLeftString               = "(d) "
  res@gsnRightString              = "JJA"
  res@tmXBOn               = True
  plot(6)  = gsn_csm_contour_map_ce(wks,reg_sst_JJA,res)
  
  delete(res@cnLevels)
  res@cnLevels                    = (/-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.1,0.2,0.3,0.4,0.5,0.6/)
  
  res@tmYLOn               = False
  res@tmXBOn               = False

  res@lbLabelBarOn                = False
  res@gsnCenterString             = "Precip(GPCP)"
  res@gsnLeftString               = "(d) "
  res@gsnRightString              = "AM"
  plot(1)  = gsn_csm_contour_map_ce(wks,reg_precip_AM,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(e) "
  res@gsnRightString              = "MJ"
  plot(3)  = gsn_csm_contour_map_ce(wks,reg_precip_MJ,res)
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
  plot(5)  = gsn_csm_contour_map_ce(wks,reg_precip_JJA,res)
  
  res@lbLabelBarOn             = True
  res@lbLabelAutoStride        = True
  res@lbBoxLineColor           = 0
  res@lbBoxLineThicknessF      = 3
  res@pmLabelBarOrthogonalPosF = 0.2
  res@pmLabelBarWidthF         = 0.65
  res@pmLabelBarHeightF        = 0.09
  res@lbLabelFontHeightF       = 0.025
  res@gsnLeftString            = "(h) "
  res@gsnRightString           = "JJA"
  res@tmXBOn                   = True
  plot(7)  = gsn_csm_contour_map_ce(wks,reg_precip_JJA,res)
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
  
  plot_test     = new(8,graphic)
  
  plot_test(0)   = gsn_csm_contour(wks,sst_AM_test    ,res_t)
  plot_test(2)   = gsn_csm_contour(wks,sst_MJ_test    ,res_t)
  plot_test(4)   = gsn_csm_contour(wks,sst_JJA_test    ,res_t)
  plot_test(6)   = gsn_csm_contour(wks,sst_JJA_test    ,res_t)
  
  plot_test(1)   = gsn_csm_contour(wks,precip_AM_test    ,res_t)
  plot_test(3)   = gsn_csm_contour(wks,precip_MJ_test    ,res_t)
  plot_test(5)   = gsn_csm_contour(wks,precip_JJA_test    ,res_t)
  plot_test(7)   = gsn_csm_contour(wks,precip_JJA_test    ,res_t)
  
  ;===
  do i=0,7
  overlay(plot(i),plot_test(i))
  end do
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
  
  reg_uwnd_AM = mask(reg_uwnd_AM,reg_uwnd_AM.lt.0.1 .and. reg_uwnd_AM.gt.-0.1.and.reg_vwnd_AM.lt.0.1 .and. reg_vwnd_AM.gt.-0.1,False)
  reg_uwnd_MJ = mask(reg_uwnd_MJ,reg_uwnd_MJ.lt.0.1 .and. reg_uwnd_MJ.gt.-0.1.and.reg_vwnd_MJ.lt.0.1 .and. reg_vwnd_MJ.gt.-0.1,False)
  reg_uwnd_JJA = mask(reg_uwnd_JJA,reg_uwnd_JJA.lt.0.1 .and. reg_uwnd_JJA.gt.-0.1.and.reg_vwnd_JJA.lt.0.1 .and. reg_vwnd_JJA.gt.-0.1,False)
  reg_uwnd_JJ = mask(reg_uwnd_JJ,reg_uwnd_JJ.lt.0.1 .and. reg_uwnd_JJ.gt.-0.1.and.reg_vwnd_JJ.lt.0.1 .and. reg_vwnd_JJ.gt.-0.1,False)
                                                                                                                                    
  reg_vwnd_AM  = mask(reg_vwnd_AM,reg_uwnd_AM.lt.0.1 .and. reg_uwnd_AM.gt.-0.1.and.reg_vwnd_AM.lt.0.1 .and. reg_vwnd_AM.gt.-0.1,False)
  reg_vwnd_MJ  = mask(reg_vwnd_MJ,reg_uwnd_MJ.lt.0.1 .and. reg_uwnd_MJ.gt.-0.1.and.reg_vwnd_MJ.lt.0.1 .and. reg_vwnd_MJ.gt.-0.1,False)
  reg_vwnd_JJA = mask(reg_vwnd_JJA,reg_uwnd_JJA.lt.0.1 .and. reg_uwnd_JJA.gt.-0.1.and.reg_vwnd_JJA.lt.0.1 .and. reg_vwnd_JJA.gt.-0.1,False)
  reg_vwnd_JJ  = mask(reg_vwnd_JJ,reg_uwnd_JJ.lt.0.1 .and. reg_uwnd_JJ.gt.-0.1.and.reg_vwnd_JJ.lt.0.1 .and. reg_vwnd_JJ.gt.-0.1,False)
  
  plot_vc = new(4,graphic)
    
  res_vc@vcRefMagnitudeF         = 0.5
  
  res_vc@gsnRightStringFontHeightF   = 0.025
  res_vc@gsnLeftStringFontHeightF    = 0.025
  
  res_vc@tmXBLabelFontHeightF = 0.025
  res_vc@tmYLLabelFontHeightF = 0.025
  
  res_vc@gsnLeftString               = ""
  res_vc@gsnRightString              = " "
  plot_vc(0)  = gsn_csm_vector(wks,reg_uwnd_AM(::2,::2),reg_vwnd_AM(::2,::2),res_vc)
  
  plot_vc(1)  = gsn_csm_vector(wks,reg_uwnd_MJ(::2,::2),reg_vwnd_MJ(::2,::2),res_vc)
  res_vc@vcRefAnnoOrthogonalPosF = -1.455
  plot_vc(2)  = gsn_csm_vector(wks,reg_uwnd_JJA(::2,::2),reg_vwnd_JJA(::2,::2),res_vc)
  
  plot_vc(3)  = gsn_csm_vector(wks,reg_uwnd_JJA(::2,::2),reg_vwnd_JJA(::2,::2),res_vc)
  
  overlay(plot(0),plot_vc(0))
  overlay(plot(2),plot_vc(1))
  overlay(plot(4),plot_vc(2))
  overlay(plot(6),plot_vc(3))
  
  rest                   = True                      ; polyline mods desired
  rest@gsLineColor       = "red"                    ; color of lines
  rest@gsLineThicknessF  = 2                       ; thickness of lines
  rest@gsLineDashPattern = 0
  
  ;===
  ypts_A = (/ 20.0,  20.0,  -10.0,  -10.0, 20.0/)      ;latitude
  xpts_A = (/ 40.0,  90.0,  90.0,   40.0,  40.0/)    ;longitude
  
  dum_A = new(4,graphic)
  
  do i = 0, 3
   dum_A(i)=gsn_add_polyline(wks,plot(4),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
  end do
  
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarParallelPosF   = -0.02 ;负值左移
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