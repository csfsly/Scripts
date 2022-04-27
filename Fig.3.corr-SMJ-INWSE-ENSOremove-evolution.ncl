load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
  
begin
	
	;===
	PC2_NCEP   = fbindirread ("hgt200-PC2-NCEP1-1979-2019.dat" ,0,41,"float")
	INWSE_NCEP = fbindirread ("Index-NWSE-NCEP1-1979-2019.dat" ,0,41,"float")
	
	PC2_NCEP   = dim_standardize(PC2_NCEP,1)
	INWSE_NCEP = dim_standardize(INWSE_NCEP,1)
	
	print(escorc(INWSE_NCEP,PC2_NCEP))
	
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
	rc_PC2_NCEP      = regCoef(Nino34_DJF,PC2_NCEP)
	PC2_NCEP_Indep   = PC2_NCEP  - rc_PC2_NCEP  *Nino34_DJF  ;#########
  
	;===
  
	;======
	print("...read in data files...")
	
  f_vwnd  = addfile("/disk_raid6/data/NCEP-2.5/ncep.reanalysis.derived/pressure/vwnd.mon.mean.nc", "r")
  
  delete(TIME)
  delete(YMD)
  
	TIME    = f_vwnd->time                       ; days since ...
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
  
  vwnd_Dec   = f_vwnd->vwnd(iStrt_Dec:iLast_Dec:12,{925},::-1,:) 
  vwnd_Jan   = f_vwnd->vwnd(iStrt_Jan:iLast_Jan:12,{925},::-1,:) 
  vwnd_Feb   = f_vwnd->vwnd(iStrt_Feb:iLast_Feb:12,{925},::-1,:) 
  
  vwnd_Mar   = short2flt(f_vwnd->vwnd(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  vwnd_Apr   = short2flt(f_vwnd->vwnd(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  vwnd_May   = short2flt(f_vwnd->vwnd(iStrt_May:iLast_May:12,{925},::-1,:)) 
  
  vwnd_Jun   = f_vwnd->vwnd(iStrt_Jun:iLast_Jun:12,{925},::-1,:) 
  vwnd_Jul   = f_vwnd->vwnd(iStrt_Jul:iLast_Jul:12,{925},::-1,:) 
  vwnd_Aug   = f_vwnd->vwnd(iStrt_Aug:iLast_Aug:12,{925},::-1,:) 
  
	vwnd_AM   = (vwnd_Apr+vwnd_May)/2.
  vwnd_MA   = (vwnd_Mar+vwnd_Apr)/2.
  vwnd_MJ   = (vwnd_May+vwnd_Jun)/2.
  vwnd_JJ   = (vwnd_Jun+vwnd_Jul)/2.
  vwnd_MAM  = (vwnd_Mar+vwnd_Apr+vwnd_May)/3.
  
  copy_VarMeta(vwnd_Mar,vwnd_AM)
  copy_VarMeta(vwnd_Mar,vwnd_MA)
  copy_VarMeta(vwnd_Mar,vwnd_MJ)
  copy_VarMeta(vwnd_Mar,vwnd_JJ)
  copy_VarMeta(vwnd_Mar,vwnd_MAM)
	;===
	;===高通滤波
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  vwnd_AM_IA    = wgt_runave_Wrap(vwnd_AM(lat|:,lon|:,time|:),wgt,1)
  vwnd_MA_IA    = wgt_runave_Wrap(vwnd_MA(lat|:,lon|:,time|:),wgt,1)
  vwnd_MJ_IA    = wgt_runave_Wrap(vwnd_MJ(lat|:,lon|:,time|:),wgt,1)
  vwnd_JJ_IA    = wgt_runave_Wrap(vwnd_JJ(lat|:,lon|:,time|:),wgt,1)
  vwnd_MAM_IA   = wgt_runave_Wrap(vwnd_MAM(lat|:,lon|:,time|:),wgt,1)
  
  vwnd_Mar_IA    = wgt_runave_Wrap(vwnd_Mar(lat|:,lon|:,time|:),wgt,1)
  vwnd_Apr_IA    = wgt_runave_Wrap(vwnd_Apr(lat|:,lon|:,time|:),wgt,1)
  vwnd_May_IA    = wgt_runave_Wrap(vwnd_May(lat|:,lon|:,time|:),wgt,1)
  vwnd_Jun_IA    = wgt_runave_Wrap(vwnd_Jun(lat|:,lon|:,time|:),wgt,1)
  vwnd_Jul_IA    = wgt_runave_Wrap(vwnd_Jul(lat|:,lon|:,time|:),wgt,1)
  vwnd_Aug_IA    = wgt_runave_Wrap(vwnd_Aug(lat|:,lon|:,time|:),wgt,1)
  
  
  SJI_AM = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_AM_IA({-5:5},{40:50},:),0),0)
  SJI_MA = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_MA_IA({-5:5},{40:50},:),0),0)
  SJI_MJ = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_MJ_IA({-5:5},{40:50},:),0),0)
  SJI_JJ = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_JJ_IA({-5:5},{40:50},:),0),0)
  SJI_MAM = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_MAM_IA({-5:5},{40:50},:),0),0)
  
  SJI_Mar = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_Mar_IA({-5:5},{40:50},:),0),0)
  SJI_Apr = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_Apr_IA({-5:5},{40:50},:),0),0)
  SJI_May = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_May_IA({-5:5},{40:50},:),0),0)
  SJI_Jun = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_Jun_IA({-5:5},{40:50},:),0),0)
  SJI_Jul = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_Jul_IA({-5:5},{40:50},:),0),0)
  SJI_Aug = dim_avg_n_Wrap(dim_avg_n_Wrap(vwnd_Aug_IA({-5:5},{40:50},:),0),0)
  
  rc_SJI_Mar      = regCoef(Nino34_DJF,SJI_Mar)
  rc_SJI_Apr      = regCoef(Nino34_DJF,SJI_Apr)
  rc_SJI_May      = regCoef(Nino34_DJF,SJI_May)
  rc_SJI_Jun      = regCoef(Nino34_DJF,SJI_Jun)
  rc_SJI_Jul      = regCoef(Nino34_DJF,SJI_Jul)
  rc_SJI_Aug      = regCoef(Nino34_DJF,SJI_Aug)
  
	SJI_Mar_Indep   = SJI_Mar     - rc_SJI_Mar  *Nino34_DJF  ;#########
	SJI_Apr_Indep   = SJI_Apr     - rc_SJI_Apr  *Nino34_DJF  ;#########
	SJI_May_Indep   = SJI_May     - rc_SJI_May  *Nino34_DJF  ;#########
	SJI_Jun_Indep   = SJI_Jun     - rc_SJI_Jun  *Nino34_DJF  ;#########
	SJI_Jul_Indep   = SJI_Jul     - rc_SJI_Jul  *Nino34_DJF  ;#########
	SJI_Aug_Indep   = SJI_Aug     - rc_SJI_Aug  *Nino34_DJF  ;#########
	
	rc_SJI_AM      = regCoef(Nino34_DJF,SJI_AM)
	SJI_AM_Indep   = SJI_AM     - rc_SJI_AM  *Nino34_DJF  ;#########
	
	rc_SJI_MAM     = regCoef(Nino34_DJF,SJI_MAM)
	SJI_MAM_Indep  = SJI_MAM     - rc_SJI_MAM  *Nino34_DJF  ;#########
	
  ;===
  print("two monthly")
  print(escorc(SJI_MA,PC2_NCEP))
  print(escorc(SJI_AM,PC2_NCEP))
  print(escorc(SJI_MJ,PC2_NCEP))
  print(escorc(SJI_JJ,PC2_NCEP))
  
  print("monthly")
  
  print(escorc(SJI_Mar,PC2_NCEP))
  print(escorc(SJI_Apr,PC2_NCEP))
  print(escorc(SJI_May,PC2_NCEP))
  print(escorc(SJI_Jun,PC2_NCEP))
  print(escorc(SJI_Jul,PC2_NCEP))
  print(escorc(SJI_Aug,PC2_NCEP))
  
  print("two monthly Indep")
  print(escorc(SJI_MA,PC2_NCEP_Indep))
  print(escorc(SJI_AM,PC2_NCEP_Indep))
  print(escorc(SJI_MJ,PC2_NCEP_Indep))
  print(escorc(SJI_JJ,PC2_NCEP_Indep))
  
  print("monthly Indep")
  
  print(escorc(SJI_Mar_Indep,PC2_NCEP))
  print(escorc(SJI_Apr_Indep,PC2_NCEP))
  print(escorc(SJI_May_Indep,PC2_NCEP))
  print(escorc(SJI_Jun_Indep,PC2_NCEP))
  print(escorc(SJI_Jul_Indep,PC2_NCEP))
  print(escorc(SJI_Aug_Indep,PC2_NCEP))
  
  print("Indep")
  print(escorc(SJI_May_Indep,PC2_NCEP))
  print(escorc(SJI_AM_Indep,PC2_NCEP))
  print(escorc(SJI_MAM_Indep,PC2_NCEP))
  
  ;===
  system("rm -rf SJI_AM_Indep-1979-2019.dat")
  system("rm -rf SJI_MAM_Indep-1979-2019.dat")
  
  fbindirwrite("SJI_AM_Indep-1979-2019.dat"  ,SJI_AM_Indep)
  fbindirwrite("SJI_MAM_Indep-1979-2019.dat"  ,SJI_MAM_Indep)
  
  system("rm -rf SJI_AM-1979-2019.dat")
  system("rm -rf SJI_MAM-1979-2019.dat")
  
  fbindirwrite("SJI_AM-1979-2019.dat"  ,SJI_AM)
  fbindirwrite("SJI_MAM-1979-2019.dat"  ,SJI_MAM)
  ;===
  
  corr = new(6,"float")
  corr_Indep = new(6,"float")
  
  corr(0) = escorc(SJI_Mar,PC2_NCEP)
  corr(1) = escorc(SJI_Apr,PC2_NCEP)
  corr(2) = escorc(SJI_May,PC2_NCEP)
  corr(3) = escorc(SJI_Jun,PC2_NCEP)
  corr(4) = escorc(SJI_Jul,PC2_NCEP)
  corr(5) = escorc(SJI_Aug,PC2_NCEP)
  
  corr_Indep(0) = escorc(SJI_Mar_Indep,PC2_NCEP)
  corr_Indep(1) = escorc(SJI_Apr_Indep,PC2_NCEP)
  corr_Indep(2) = escorc(SJI_May_Indep,PC2_NCEP)
  corr_Indep(3) = escorc(SJI_Jun_Indep,PC2_NCEP)
  corr_Indep(4) = escorc(SJI_Jul_Indep,PC2_NCEP)
  corr_Indep(5) = escorc(SJI_Aug_Indep,PC2_NCEP)
  
  ;===
  f_hgt  = addfile("/disk_raid6/data/NCEP-2.5/ncep.reanalysis.derived/pressure/hgt.mon.mean.nc", "r")
  
	TIME    = f_hgt->time                       ; days since ...
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
  
  hgt_Dec   = f_hgt->hgt(iStrt_Dec:iLast_Dec:12,{200},::-1,:) 
  hgt_Jan   = f_hgt->hgt(iStrt_Jan:iLast_Jan:12,{200},::-1,:) 
  hgt_Feb   = f_hgt->hgt(iStrt_Feb:iLast_Feb:12,{200},::-1,:) 
  
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
  
  ;===regression
  reg_hgt_JJA = regCoef(SJI_AM_Indep,hgt_JJA_IA(lat|:,lon|:,time|:))
  
  hgt_JJA_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  hgt_JJA_test(j,i) = reg_hgt_JJA@tval(k)
  k = k +1 
  end do
  end do
  
  reg_hgt_JJA!0   = "lat"
  reg_hgt_JJA!1   = "lon"
  reg_hgt_JJA&lat = hgt_JJA_IA&lat
  reg_hgt_JJA&lon = hgt_JJA_IA&lon
  
  copy_VarMeta(reg_hgt_JJA,hgt_JJA_test)
  ;===
  ;===plot
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.3.corr-SMJ-INWSE-ENSOremove-evlution")
  gsn_define_colormap(wks,"BlueDarkRed18")
  
  
  res = True
  
  res@gsnDraw  = False
  res@gsnFrame = False
  
  res@trXMaxF  = 6.5
  res@trXMinF  = 0.5
  
  res@vpWidthF  = 0.7
  res@vpHeightF = 0.35
  
  res@xyLineColors           = (/"darkslategray1","blue","blue","blue","blue"/)
  res@xyLineThicknesses      = (/3,1.,1.,1.,1./)
  res@xyDashPatterns         = (/0,0,3,4,5/)
  
  res@xyMarkLineMode    = "Lines"
  res@xyMarkers         = (/16,16/)     
  res@xyMonoMarkerSize  = True
  res@xyMarkerSizeF     = 0.006
  res@xyMarkerColors    = (/"red","blue"/)
  
  res@tiMainString       = " "
  res@tiYAxisString      = " "
  
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  
  res@tmYRMode             = "Automatic"
  res@tmXTOn               = False
  res@tmYROn               = False
  
  res@tmYLMinorLengthF  = 0.005
  res@tmYLMajorLengthF  = 0.01
  res@tmXBMinorLengthF  = 0.005
  res@tmXBMajorLengthF  = 0.01
  
  year  = fspan(1.,6.,6)
  
	res@gsnRightStringFontHeightF   = 0.025
  res@gsnLeftStringFontHeightF    = 0.025
  res@gsnRightString              = ""
  res@gsnLeftString               = ""
  
  res@gsnRightString              = ""
  res@gsnLeftString               = ""
      
  res@tmXBMajorOutwardLengthF = 0.01	; set tickmarks 
  res@tmXBMinorOutwardLengthF = 0.005
  res@tmXBMajorLengthF        = 0.01
  res@tmXBMinorLengthF        = 0.005
  
  res@tmYLMajorOutwardLengthF = 0.01	; set tickmarks 
  res@tmYLMinorOutwardLengthF = 0.005
  res@tmYLMajorLengthF        = 0.01
  res@tmYLMinorLengthF        = 0.005
  
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/"1","2","3","4","5","6"/)
  res@tmXBLabels = (/"Mar","Apr","May","Jun","Jul","Aug"/)
  
  res@trYMaxF  = 0.61
  res@trYMinF  = 0.0
  
  res@gsnRightString               = ""
  res@gsnLeftString                = ""
  res@gsnLeftStringOrthogonalPosF   = 0.001
  res@gsnRightStringOrthogonalPosF  = 0.001
  
  res@gsnYRefLine           = 0.              ; reference line   
  res@gsnXYBarChart         = True            ; create bar chart 
  res@gsnAboveYRefLineColor = "darkslategray1"           ; above ref line fill red
  res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
  res@gsnXYBarChartBarWidth = 0.4
  
  res@tmYLMode          = "Explicit"
  res@tmYLValues        = (/"0","0.2","0.4","0.6"/)
  res@tmYLLabels        = (/"0","0.2","0.4","0.6"/)
  
  plot     = new(3,graphic)
  
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(a) "
  res@gsnRightString              = "corr(SMJI&SAH)"
  plot(0)       = gsn_csm_xy(wks,year(:),corr,res)
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(b) "
  res@gsnRightString              = "corr(SMJI&SAH) ENSO remove"
  plot(1)       = gsn_csm_xy(wks,year(:),corr_Indep,res)
   
  res_t = True
  
  res_t@gsnDraw  = False
  res_t@gsnFrame = False
  
  res_t@xyLineColors           = (/"red","black"/)
  res_t@xyLineThicknesses      = (/2,1.,1.,1.,1./)
  res_t@xyDashPatterns         = (/12,0,3,4,5/)
  res_t@gsnXYBarChart         = False
  
  plot_test     = new(3,graphic)
  year1  = fspan(0.,7.,8)
  plot_test(0)        = gsn_csm_xy(wks,year1,year1*0+0.31,res_t)
  plot_test(1)        = gsn_csm_xy(wks,year1,year1*0+0.31,res_t)
  
  overlay(plot(0),plot_test(0))
  overlay(plot(1),plot_test(1))
  
  delete(res)
  delete(res_t)
  
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
  
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  
  res@gsnLeftStringFontHeightF   = 0.025
  res@gsnRightStringFontHeightF  = 0.025
  res@gsnCenterStringFontHeightF = 0.025
  
  res@gsnCenterStringOrthogonalPosF = 0.1
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
  res@gsnCenterString             = ""
  res@gsnLeftString               = "(c) "
  res@gsnRightString              = "HGT(200hPa) wrt SMJI"
  
  res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels             = (/-18,-14,-10,-6,-2,2,6,10,14,18/)
  res@cnLevels             = (/-19,-16,-13,-10,-7,-5,-2,2,5,7,10,13,16,19/)
  plot(2)  = gsn_csm_contour_map_ce(wks,reg_hgt_JJA,res)
  
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
  
  plot_cn = new(2,graphic)
  plot_cn(0) = gsn_csm_contour(wks,hgt_JJA_clim,res)
  plot_cn(1) = gsn_csm_contour(wks,hgt_JJA_clim,res)
  
  overlay(plot(2),plot_cn(1))
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
    
  plot_test(2)   = gsn_csm_contour(wks,hgt_JJA_test    ,res_t)
  
  ;===
  overlay(plot(2),plot_test(2))
  
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarOrthogonalPosF = 0.01
  ;resp@pmLabelBarParallelPosF   = 0.02
  
  resp@pmLabelBarWidthF         = 0.85
  resp@pmLabelBarHeightF        = 0.08
  resp@lbLabelFontHeightF       = 0.015
  
  resp@lbBoxLineColor           = 0
  ;resp@lbBoxEndCapStyle         = "TriangleBothEnds"
  resp@lbLabelAutoStride        = True
  res@lbBoxLineThicknessF       = 3
  resp@gsnPanelRight           = 0.97
  
  resp@gsnPanelYF               = (/0.95,0.65,0.35/);-0.05
  
  gsn_panel(wks,plot,(/3,1/),resp)
  
end