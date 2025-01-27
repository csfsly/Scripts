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
	
	print(escorc(INWSE_NCEP,PC2_NCEP))
	;======
	print("...read in data files...")
	
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
  
	uwnd_AMJ   = (uwnd_Apr+uwnd_May+uwnd_Jun)/3.
  uwnd_MAM   = (uwnd_Mar+uwnd_Apr+uwnd_May)/3.
  uwnd_MJJ   = (uwnd_May+uwnd_Jun+uwnd_Jul)/3.
  uwnd_JJA   = (uwnd_Jun+uwnd_Jul+uwnd_Aug)/3.
  
  copy_VarMeta(uwnd_Mar,uwnd_AMJ)
  copy_VarMeta(uwnd_Mar,uwnd_MAM)
  copy_VarMeta(uwnd_Mar,uwnd_JJA)
  copy_VarMeta(uwnd_Mar,uwnd_MJJ)
	;===
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  uwnd_AMJ_IA    = wgt_runave_Wrap(uwnd_AMJ(lat|:,lon|:,time|:),wgt,1)
  uwnd_MAM_IA    = wgt_runave_Wrap(uwnd_MAM(lat|:,lon|:,time|:),wgt,1)
  uwnd_JJA_IA    = wgt_runave_Wrap(uwnd_JJA(lat|:,lon|:,time|:),wgt,1)
  uwnd_MJJ_IA    = wgt_runave_Wrap(uwnd_MJJ(lat|:,lon|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_uwnd_MAM = regCoef(PC2_NCEP,uwnd_MAM_IA(lat|:,lon|:,time|:))
  reg_uwnd_AMJ = regCoef(PC2_NCEP,uwnd_AMJ_IA(lat|:,lon|:,time|:))
  reg_uwnd_MJJ = regCoef(PC2_NCEP,uwnd_MJJ_IA(lat|:,lon|:,time|:))
  reg_uwnd_JJA = regCoef(PC2_NCEP,uwnd_JJA_IA(lat|:,lon|:,time|:))
  
  uwnd_AMJ_test = new((/73,144/),"float")
  uwnd_MAM_test = new((/73,144/),"float")
  uwnd_JJA_test = new((/73,144/),"float")
  uwnd_MJJ_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  uwnd_AMJ_test(j,i) = reg_uwnd_AMJ@tval(k)
  uwnd_MAM_test(j,i) = reg_uwnd_MAM@tval(k)
  uwnd_JJA_test(j,i) = reg_uwnd_JJA@tval(k)
  uwnd_MJJ_test(j,i) = reg_uwnd_MJJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_uwnd_AMJ!0   = "lat"
  reg_uwnd_AMJ!1   = "lon"
  reg_uwnd_AMJ&lat = uwnd_AMJ_IA&lat
  reg_uwnd_AMJ&lon = uwnd_AMJ_IA&lon
  
  copy_VarMeta(reg_uwnd_AMJ,reg_uwnd_MAM)
  copy_VarMeta(reg_uwnd_AMJ,reg_uwnd_JJA)
  copy_VarMeta(reg_uwnd_AMJ,reg_uwnd_MJJ)
  
  copy_VarMeta(reg_uwnd_AMJ,uwnd_AMJ_test)
  copy_VarMeta(reg_uwnd_AMJ,uwnd_MAM_test)
  copy_VarMeta(reg_uwnd_AMJ,uwnd_JJA_test)
  copy_VarMeta(reg_uwnd_AMJ,uwnd_MJJ_test)
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
  
	vwnd_AMJ   = (vwnd_Apr+vwnd_May+vwnd_Jun)/3.
  vwnd_MAM   = (vwnd_Mar+vwnd_Apr+vwnd_May)/3.
  vwnd_JJA   = (vwnd_Jun+vwnd_Jul+vwnd_Aug)/3.
  vwnd_MJJ   = (vwnd_May+vwnd_Jun+vwnd_Jul)/3.
  
  copy_VarMeta(vwnd_Mar,vwnd_AMJ)
  copy_VarMeta(vwnd_Mar,vwnd_MAM)
  copy_VarMeta(vwnd_Mar,vwnd_JJA)
  copy_VarMeta(vwnd_Mar,vwnd_MJJ)
	;===
	;===��ͨ�˲�
  ihp   = 1
  sigma = 1.0
  nWgt  = 9
  fca   = 1./9.
  wgt   = filwgts_lanczos(nWgt,ihp,fca,-999.,sigma)
  
  vwnd_AMJ_IA    = wgt_runave_Wrap(vwnd_AMJ(lat|:,lon|:,time|:),wgt,1)
  vwnd_MAM_IA    = wgt_runave_Wrap(vwnd_MAM(lat|:,lon|:,time|:),wgt,1)
  vwnd_JJA_IA    = wgt_runave_Wrap(vwnd_JJA(lat|:,lon|:,time|:),wgt,1)
  vwnd_MJJ_IA    = wgt_runave_Wrap(vwnd_MJJ(lat|:,lon|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_vwnd_AMJ = regCoef(PC2_NCEP,vwnd_AMJ_IA(lat|:,lon|:,time|:))
  reg_vwnd_MAM = regCoef(PC2_NCEP,vwnd_MAM_IA(lat|:,lon|:,time|:))
  reg_vwnd_JJA = regCoef(PC2_NCEP,vwnd_JJA_IA(lat|:,lon|:,time|:))
  reg_vwnd_MJJ = regCoef(PC2_NCEP,vwnd_MJJ_IA(lat|:,lon|:,time|:))
  
  vwnd_AMJ_test = new((/73,144/),"float")
  vwnd_MAM_test = new((/73,144/),"float")
  vwnd_JJA_test = new((/73,144/),"float")
  vwnd_MJJ_test = new((/73,144/),"float")
  
  k = 0 
  do j=0,72
  do i=0,143
  vwnd_AMJ_test(j,i) = reg_vwnd_AMJ@tval(k)
  vwnd_MAM_test(j,i) = reg_vwnd_MAM@tval(k)
  vwnd_JJA_test(j,i) = reg_vwnd_JJA@tval(k)
  vwnd_MJJ_test(j,i) = reg_vwnd_MJJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_vwnd_AMJ!0   = "lat"
  reg_vwnd_AMJ!1   = "lon"
  reg_vwnd_AMJ&lat = vwnd_AMJ_IA&lat
  reg_vwnd_AMJ&lon = vwnd_AMJ_IA&lon
  
  copy_VarMeta(reg_vwnd_AMJ,reg_vwnd_MAM)
  copy_VarMeta(reg_vwnd_AMJ,reg_vwnd_JJA)
  copy_VarMeta(reg_vwnd_AMJ,reg_vwnd_MJJ)
  
  copy_VarMeta(reg_vwnd_AMJ,vwnd_AMJ_test)
  copy_VarMeta(reg_vwnd_AMJ,vwnd_MAM_test)
  copy_VarMeta(reg_vwnd_AMJ,vwnd_JJA_test)
  copy_VarMeta(reg_vwnd_AMJ,vwnd_MJJ_test)
  ;===
  common_test_AMJ_1d  = new(10512,"float")
  common_test_AMJ     = new((/73,144/),"float")
  
  uwnd_test_AMJ_1d = ndtooned(uwnd_AMJ_test)
  vwnd_test_AMJ_1d = ndtooned(vwnd_AMJ_test)
  
  ;===
  
  do i=0,10511
  ;===
	if abs(uwnd_test_AMJ_1d(i)>1.8).or.abs(vwnd_test_AMJ_1d(i)>1.8) then
		if abs(uwnd_test_AMJ_1d(i)).gt.abs(vwnd_test_AMJ_1d(i)) then
		common_test_AMJ_1d(i)= uwnd_test_AMJ_1d(i)
	  end if
	  
	  if abs(uwnd_test_AMJ_1d(i)).lt.abs(vwnd_test_AMJ_1d(i)) then
		common_test_AMJ_1d(i)= vwnd_test_AMJ_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_AMJ = onedtond(common_test_AMJ_1d,(/73,144/))
  
  copy_VarMeta(reg_vwnd_AMJ,common_test_AMJ)
  ;===
  common_test_MJJ_1d  = new(10512,"float")
  common_test_MJJ     = new((/73,144/),"float")
  
  uwnd_test_MJJ_1d = ndtooned(uwnd_MJJ_test)
  vwnd_test_MJJ_1d = ndtooned(vwnd_MJJ_test)
  
  ;===
  
  do i=0,10511
  ;===
	if abs(uwnd_test_MJJ_1d(i)>1.8).or.abs(vwnd_test_MJJ_1d(i)>1.8) then
		if abs(uwnd_test_MJJ_1d(i)).gt.abs(vwnd_test_MJJ_1d(i)) then
		common_test_MJJ_1d(i)= uwnd_test_MJJ_1d(i)
	  end if
	  
	  if abs(uwnd_test_MJJ_1d(i)).lt.abs(vwnd_test_MJJ_1d(i)) then
		common_test_MJJ_1d(i)= vwnd_test_MJJ_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_MJJ = onedtond(common_test_MJJ_1d,(/73,144/))
  
  copy_VarMeta(reg_vwnd_MJJ,common_test_MJJ)
  ;===
  common_test_MAM_1d  = new(10512,"float")
  common_test_MAM     = new((/73,144/),"float")
  
  uwnd_test_MAM_1d = ndtooned(uwnd_MAM_test)
  vwnd_test_MAM_1d = ndtooned(vwnd_MAM_test)
  
  ;===
  
  do i=0,10511
  ;===
	if abs(uwnd_test_MAM_1d(i)>1.8).or.abs(vwnd_test_MAM_1d(i)>1.8) then
		if abs(uwnd_test_MAM_1d(i)).gt.abs(vwnd_test_MAM_1d(i)) then
		common_test_MAM_1d(i)= uwnd_test_MAM_1d(i)
	  end if
	  
	  if abs(uwnd_test_MAM_1d(i)).lt.abs(vwnd_test_MAM_1d(i)) then
		common_test_MAM_1d(i)= vwnd_test_MAM_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_MAM = onedtond(common_test_MAM_1d,(/73,144/))
  
  copy_VarMeta(reg_vwnd_MAM,common_test_MAM)
  ;===
  common_test_JJA_1d  = new(10512,"float")
  common_test_JJA     = new((/73,144/),"float")
  
  uwnd_test_JJA_1d = ndtooned(uwnd_JJA_test)
  vwnd_test_JJA_1d = ndtooned(vwnd_JJA_test)
  
  ;===
  
  do i=0,10511
  ;===
	if abs(uwnd_test_JJA_1d(i)>1.8).or.abs(vwnd_test_JJA_1d(i)>1.8) then
		if abs(uwnd_test_JJA_1d(i)).gt.abs(vwnd_test_JJA_1d(i)) then
		common_test_JJA_1d(i)= uwnd_test_JJA_1d(i)
	  end if
	  
	  if abs(uwnd_test_JJA_1d(i)).lt.abs(vwnd_test_JJA_1d(i)) then
		common_test_JJA_1d(i)= vwnd_test_JJA_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_JJA = onedtond(common_test_JJA_1d,(/73,144/))
  
  copy_VarMeta(reg_vwnd_JJA,common_test_JJA)
  ;===
  ;======
	print("...read in data files...")
	
  f_uwnd_ERA  = addfile("/disk_raid6/data/ERA-interim/pressure/monthly/uwnd.mon.mean.1.0.nc", "r")
  f_vwnd_ERA  = addfile("/disk_raid6/data/ERA-interim/pressure/monthly/vwnd.mon.mean.1.0.nc", "r")
  
  delete(TIME)
  delete(YMD)
  
	TIME    = f_uwnd_ERA->time                       ; days since ...
  YMD     = ut_calendar(TIME, -2)                ; entire (time,6)
  
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
  
  uwnd_ERA_Jan   = short2flt(f_uwnd_ERA->u(iStrt_Jan:iLast_Jan:12,{925},::-1,:)) 
  uwnd_ERA_Feb   = short2flt(f_uwnd_ERA->u(iStrt_Feb:iLast_Feb:12,{925},::-1,:)) 
  
  uwnd_ERA_Mar   = short2flt(f_uwnd_ERA->u(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  uwnd_ERA_Apr   = short2flt(f_uwnd_ERA->u(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  uwnd_ERA_May   = short2flt(f_uwnd_ERA->u(iStrt_May:iLast_May:12,{925},::-1,:)) 
  
  uwnd_ERA_Jun   = short2flt(f_uwnd_ERA->u(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  uwnd_ERA_Jul   = short2flt(f_uwnd_ERA->u(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  uwnd_ERA_Aug   = short2flt(f_uwnd_ERA->u(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	uwnd_ERA_AMJ   = (uwnd_ERA_Apr+uwnd_ERA_May+uwnd_ERA_Jun)/3.
  uwnd_ERA_MAM   = (uwnd_ERA_Mar+uwnd_ERA_Apr+uwnd_ERA_May)/3.
  uwnd_ERA_MJJ   = (uwnd_ERA_May+uwnd_ERA_Jun+uwnd_ERA_Jul)/3.
  uwnd_ERA_JJA   = (uwnd_ERA_Jun+uwnd_ERA_Jul+uwnd_ERA_Aug)/3.
  
  copy_VarMeta(uwnd_ERA_Mar,uwnd_ERA_AMJ)
  copy_VarMeta(uwnd_ERA_Mar,uwnd_ERA_MAM)
  copy_VarMeta(uwnd_ERA_Mar,uwnd_ERA_JJA)
  copy_VarMeta(uwnd_ERA_Mar,uwnd_ERA_MJJ)
	;===
  uwnd_ERA_AMJ_IA    = wgt_runave_Wrap(uwnd_ERA_AMJ(latitude|:,longitude|:,time|:),wgt,1)
  uwnd_ERA_MAM_IA    = wgt_runave_Wrap(uwnd_ERA_MAM(latitude|:,longitude|:,time|:),wgt,1)
  uwnd_ERA_JJA_IA    = wgt_runave_Wrap(uwnd_ERA_JJA(latitude|:,longitude|:,time|:),wgt,1)
  uwnd_ERA_MJJ_IA    = wgt_runave_Wrap(uwnd_ERA_MJJ(latitude|:,longitude|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_uwnd_ERA_MAM = regCoef(PC2_NCEP,uwnd_ERA_MAM_IA(latitude|:,longitude|:,time|:))
  reg_uwnd_ERA_AMJ = regCoef(PC2_NCEP,uwnd_ERA_AMJ_IA(latitude|:,longitude|:,time|:))
  reg_uwnd_ERA_MJJ = regCoef(PC2_NCEP,uwnd_ERA_MJJ_IA(latitude|:,longitude|:,time|:))
  reg_uwnd_ERA_JJA = regCoef(PC2_NCEP,uwnd_ERA_JJA_IA(latitude|:,longitude|:,time|:))
  
  uwnd_ERA_AMJ_test = new((/181,360/),"float")
  uwnd_ERA_MAM_test = new((/181,360/),"float")
  uwnd_ERA_JJA_test = new((/181,360/),"float")
  uwnd_ERA_MJJ_test = new((/181,360/),"float")
  
  k = 0 
  do j=0,180
  do i=0,359
  uwnd_ERA_AMJ_test(j,i) = reg_uwnd_ERA_AMJ@tval(k)
  uwnd_ERA_MAM_test(j,i) = reg_uwnd_ERA_MAM@tval(k)
  uwnd_ERA_JJA_test(j,i) = reg_uwnd_ERA_JJA@tval(k)
  uwnd_ERA_MJJ_test(j,i) = reg_uwnd_ERA_MJJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_uwnd_ERA_AMJ!0   = "lat"
  reg_uwnd_ERA_AMJ!1   = "lon"
  reg_uwnd_ERA_AMJ&lat = uwnd_ERA_AMJ_IA&latitude
  reg_uwnd_ERA_AMJ&lon = uwnd_ERA_AMJ_IA&longitude
  
  copy_VarMeta(reg_uwnd_ERA_AMJ,reg_uwnd_ERA_MAM)
  copy_VarMeta(reg_uwnd_ERA_AMJ,reg_uwnd_ERA_JJA)
  copy_VarMeta(reg_uwnd_ERA_AMJ,reg_uwnd_ERA_MJJ)
  
  copy_VarMeta(reg_uwnd_ERA_AMJ,uwnd_ERA_AMJ_test)
  copy_VarMeta(reg_uwnd_ERA_AMJ,uwnd_ERA_MAM_test)
  copy_VarMeta(reg_uwnd_ERA_AMJ,uwnd_ERA_JJA_test)
  copy_VarMeta(reg_uwnd_ERA_AMJ,uwnd_ERA_MJJ_test)
  ;===
  vwnd_ERA_Jan   = short2flt(f_vwnd_ERA->v(iStrt_Jan:iLast_Jan:12,{925},::-1,:)) 
  vwnd_ERA_Feb   = short2flt(f_vwnd_ERA->v(iStrt_Feb:iLast_Feb:12,{925},::-1,:)) 
  
  vwnd_ERA_Mar   = short2flt(f_vwnd_ERA->v(iStrt_Mar:iLast_Mar:12,{925},::-1,:)) 
  vwnd_ERA_Apr   = short2flt(f_vwnd_ERA->v(iStrt_Apr:iLast_Apr:12,{925},::-1,:)) 
  vwnd_ERA_May   = short2flt(f_vwnd_ERA->v(iStrt_May:iLast_May:12,{925},::-1,:)) 
  
  vwnd_ERA_Jun   = short2flt(f_vwnd_ERA->v(iStrt_Jun:iLast_Jun:12,{925},::-1,:)) 
  vwnd_ERA_Jul   = short2flt(f_vwnd_ERA->v(iStrt_Jul:iLast_Jul:12,{925},::-1,:)) 
  vwnd_ERA_Aug   = short2flt(f_vwnd_ERA->v(iStrt_Aug:iLast_Aug:12,{925},::-1,:)) 
  
	vwnd_ERA_AMJ   = (vwnd_ERA_Apr+vwnd_ERA_May+vwnd_ERA_Jun)/3.
  vwnd_ERA_MAM   = (vwnd_ERA_Mar+vwnd_ERA_Apr+vwnd_ERA_May)/3.
  vwnd_ERA_JJA   = (vwnd_ERA_Jun+vwnd_ERA_Jul+vwnd_ERA_Aug)/3.
  vwnd_ERA_MJJ   = (vwnd_ERA_May+vwnd_ERA_Jun+vwnd_ERA_Jul)/3.
  
  copy_VarMeta(vwnd_ERA_Mar,vwnd_ERA_AMJ)
  copy_VarMeta(vwnd_ERA_Mar,vwnd_ERA_MAM)
  copy_VarMeta(vwnd_ERA_Mar,vwnd_ERA_JJA)
  copy_VarMeta(vwnd_ERA_Mar,vwnd_ERA_MJJ)
	;===
  vwnd_ERA_AMJ_IA    = wgt_runave_Wrap(vwnd_ERA_AMJ(latitude|:,longitude|:,time|:),wgt,1)
  vwnd_ERA_MAM_IA    = wgt_runave_Wrap(vwnd_ERA_MAM(latitude|:,longitude|:,time|:),wgt,1)
  vwnd_ERA_JJA_IA    = wgt_runave_Wrap(vwnd_ERA_JJA(latitude|:,longitude|:,time|:),wgt,1)
  vwnd_ERA_MJJ_IA    = wgt_runave_Wrap(vwnd_ERA_MJJ(latitude|:,longitude|:,time|:),wgt,1)
  ;===
  
  ;===regression
  reg_vwnd_ERA_AMJ = regCoef(PC2_NCEP,vwnd_ERA_AMJ_IA(latitude|:,longitude|:,time|:))
  reg_vwnd_ERA_MAM = regCoef(PC2_NCEP,vwnd_ERA_MAM_IA(latitude|:,longitude|:,time|:))
  reg_vwnd_ERA_JJA = regCoef(PC2_NCEP,vwnd_ERA_JJA_IA(latitude|:,longitude|:,time|:))
  reg_vwnd_ERA_MJJ = regCoef(PC2_NCEP,vwnd_ERA_MJJ_IA(latitude|:,longitude|:,time|:))
  
  vwnd_ERA_AMJ_test = new((/181,360/),"float")
  vwnd_ERA_MAM_test = new((/181,360/),"float")
  vwnd_ERA_JJA_test = new((/181,360/),"float")
  vwnd_ERA_MJJ_test = new((/181,360/),"float")
  
  k = 0 
  do j=0,180
  do i=0,359
  vwnd_ERA_AMJ_test(j,i) = reg_vwnd_ERA_AMJ@tval(k)
  vwnd_ERA_MAM_test(j,i) = reg_vwnd_ERA_MAM@tval(k)
  vwnd_ERA_JJA_test(j,i) = reg_vwnd_ERA_JJA@tval(k)
  vwnd_ERA_MJJ_test(j,i) = reg_vwnd_ERA_MJJ@tval(k)
  k = k +1 
  end do
  end do
  
  reg_vwnd_ERA_AMJ!0   = "lat"
  reg_vwnd_ERA_AMJ!1   = "lon"
  reg_vwnd_ERA_AMJ&lat = vwnd_ERA_AMJ_IA&latitude
  reg_vwnd_ERA_AMJ&lon = vwnd_ERA_AMJ_IA&longitude
  
  copy_VarMeta(reg_vwnd_ERA_AMJ,reg_vwnd_ERA_MAM)
  copy_VarMeta(reg_vwnd_ERA_AMJ,reg_vwnd_ERA_JJA)
  copy_VarMeta(reg_vwnd_ERA_AMJ,reg_vwnd_ERA_MJJ)
  
  copy_VarMeta(reg_vwnd_ERA_AMJ,vwnd_ERA_AMJ_test)
  copy_VarMeta(reg_vwnd_ERA_AMJ,vwnd_ERA_MAM_test)
  copy_VarMeta(reg_vwnd_ERA_AMJ,vwnd_ERA_JJA_test)
  copy_VarMeta(reg_vwnd_ERA_AMJ,vwnd_ERA_MJJ_test)
  ;===
  common_test_ERA_AMJ_1d  = new(65160,"float")
  common_test_ERA_AMJ     = new((/181,360/),"float")
  
  uwnd_ERA_test_AMJ_1d = ndtooned(uwnd_ERA_AMJ_test)
  vwnd_ERA_test_AMJ_1d = ndtooned(vwnd_ERA_AMJ_test)
  
  ;===
  
  do i=0,65159
  ;===
	if abs(uwnd_ERA_test_AMJ_1d(i)>1.8).or.abs(vwnd_ERA_test_AMJ_1d(i)>1.8) then
		if abs(uwnd_ERA_test_AMJ_1d(i)).gt.abs(vwnd_ERA_test_AMJ_1d(i)) then
		common_test_ERA_AMJ_1d(i)= uwnd_ERA_test_AMJ_1d(i)
	  end if
	  
	  if abs(uwnd_ERA_test_AMJ_1d(i)).lt.abs(vwnd_ERA_test_AMJ_1d(i)) then
		common_test_ERA_AMJ_1d(i)= vwnd_ERA_test_AMJ_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_ERA_AMJ = onedtond(common_test_ERA_AMJ_1d,(/181,360/))
  
  copy_VarMeta(reg_vwnd_ERA_AMJ,common_test_ERA_AMJ)
  ;===
  common_test_ERA_MJJ_1d  = new(65160,"float")
  common_test_ERA_MJJ     = new((/181,360/),"float")
  
  uwnd_ERA_test_MJJ_1d = ndtooned(uwnd_ERA_MJJ_test)
  vwnd_ERA_test_MJJ_1d = ndtooned(vwnd_ERA_MJJ_test)
  
  ;===
  
  do i=0,65159
  ;===
	if abs(uwnd_ERA_test_MJJ_1d(i)>1.8).or.abs(vwnd_ERA_test_MJJ_1d(i)>1.8) then
		if abs(uwnd_ERA_test_MJJ_1d(i)).gt.abs(vwnd_ERA_test_MJJ_1d(i)) then
		common_test_ERA_MJJ_1d(i)= uwnd_ERA_test_MJJ_1d(i)
	  end if
	  
	  if abs(uwnd_ERA_test_MJJ_1d(i)).lt.abs(vwnd_ERA_test_MJJ_1d(i)) then
		common_test_ERA_MJJ_1d(i)= vwnd_ERA_test_MJJ_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_ERA_MJJ = onedtond(common_test_ERA_MJJ_1d,(/181,360/))
  
  copy_VarMeta(reg_vwnd_ERA_MJJ,common_test_ERA_MJJ)
  ;===
  common_test_ERA_MAM_1d  = new(65160,"float")
  common_test_ERA_MAM     = new((/181,360/),"float")
  
  uwnd_ERA_test_MAM_1d = ndtooned(uwnd_ERA_MAM_test)
  vwnd_ERA_test_MAM_1d = ndtooned(vwnd_ERA_MAM_test)
  
  ;===
  
  do i=0,65159
  ;===
	if abs(uwnd_ERA_test_MAM_1d(i)>1.8).or.abs(vwnd_ERA_test_MAM_1d(i)>1.8) then
		if abs(uwnd_ERA_test_MAM_1d(i)).gt.abs(vwnd_ERA_test_MAM_1d(i)) then
		common_test_ERA_MAM_1d(i)= uwnd_ERA_test_MAM_1d(i)
	  end if
	  
	  if abs(uwnd_ERA_test_MAM_1d(i)).lt.abs(vwnd_ERA_test_MAM_1d(i)) then
		common_test_ERA_MAM_1d(i)= vwnd_ERA_test_MAM_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_ERA_MAM = onedtond(common_test_ERA_MAM_1d,(/181,360/))
  
  copy_VarMeta(reg_vwnd_ERA_MAM,common_test_ERA_MAM)
  ;===
  common_test_ERA_JJA_1d  = new(65160,"float")
  common_test_ERA_JJA     = new((/181,360/),"float")
  
  uwnd_ERA_test_JJA_1d = ndtooned(uwnd_ERA_JJA_test)
  vwnd_ERA_test_JJA_1d = ndtooned(vwnd_ERA_JJA_test)
  
  ;===
  do i=0,65159
  ;===
	if abs(uwnd_ERA_test_JJA_1d(i)>1.8).or.abs(vwnd_ERA_test_JJA_1d(i)>1.8) then
		if abs(uwnd_ERA_test_JJA_1d(i)).gt.abs(vwnd_ERA_test_JJA_1d(i)) then
		common_test_ERA_JJA_1d(i)= uwnd_ERA_test_JJA_1d(i)
	  end if
	  
	  if abs(uwnd_ERA_test_JJA_1d(i)).lt.abs(vwnd_ERA_test_JJA_1d(i)) then
		common_test_ERA_JJA_1d(i)= vwnd_ERA_test_JJA_1d(i)
	  end if
	end if
	;===
  end do
  
  common_test_ERA_JJA = onedtond(common_test_ERA_JJA_1d,(/181,360/))
  
  copy_VarMeta(reg_vwnd_ERA_JJA,common_test_ERA_JJA)
  ;===
  colors =   (/(/255,255,255/),(/0,0,0/),(/101,101,255/),(/255,255,255/),\
           (/255,101,101/) , (/255,0,0/)/)*1.0
  colors = colors / 255.0
  
  ;===
  wks  = gsn_open_wks("eps","Fig.2")
  gsn_define_colormap(wks,"BlueDarkRed18")
  
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
  
  res_vc@gsnLeftStringOrthogonalPosF   = 0.03
  res_vc@gsnRightStringOrthogonalPosF  = 0.03
  res_vc@gsnCenterStringOrthogonalPosF = 0.05
  
  res_vc@mpCenterLonF = 180
  res_vc@mpMaxLatF    = 35.
  res_vc@mpMinLatF    = -20
  res_vc@mpMaxLonF    = 125
  res_vc@mpMinLonF    = 30
  
  reg_uwnd_AMJ = mask(reg_uwnd_AMJ,reg_uwnd_AMJ.lt.0.1 .and. reg_uwnd_AMJ.gt.-0.1.and.reg_vwnd_AMJ.lt.0.1 .and. reg_vwnd_AMJ.gt.-0.1,False)
  reg_uwnd_MAM = mask(reg_uwnd_MAM,reg_uwnd_MAM.lt.0.1 .and. reg_uwnd_MAM.gt.-0.1.and.reg_vwnd_MAM.lt.0.1 .and. reg_vwnd_MAM.gt.-0.1,False)
  reg_uwnd_JJA = mask(reg_uwnd_JJA,reg_uwnd_JJA.lt.0.1 .and. reg_uwnd_JJA.gt.-0.1.and.reg_vwnd_JJA.lt.0.1 .and. reg_vwnd_JJA.gt.-0.1,False)
  reg_uwnd_MJJ = mask(reg_uwnd_MJJ,reg_uwnd_MJJ.lt.0.1 .and. reg_uwnd_MJJ.gt.-0.1.and.reg_vwnd_MJJ.lt.0.1 .and. reg_vwnd_MJJ.gt.-0.1,False)
                                                                                                                                    
  reg_vwnd_AMJ = mask(reg_vwnd_AMJ,reg_uwnd_AMJ.lt.0.1 .and. reg_uwnd_AMJ.gt.-0.1.and.reg_vwnd_AMJ.lt.0.1 .and. reg_vwnd_AMJ.gt.-0.1,False)
  reg_vwnd_MAM = mask(reg_vwnd_MAM,reg_uwnd_MAM.lt.0.1 .and. reg_uwnd_MAM.gt.-0.1.and.reg_vwnd_MAM.lt.0.1 .and. reg_vwnd_MAM.gt.-0.1,False)
  reg_vwnd_JJA = mask(reg_vwnd_JJA,reg_uwnd_JJA.lt.0.1 .and. reg_uwnd_JJA.gt.-0.1.and.reg_vwnd_JJA.lt.0.1 .and. reg_vwnd_JJA.gt.-0.1,False)
  reg_vwnd_MJJ = mask(reg_vwnd_MJJ,reg_uwnd_MJJ.lt.0.1 .and. reg_uwnd_MJJ.gt.-0.1.and.reg_vwnd_MJJ.lt.0.1 .and. reg_vwnd_MJJ.gt.-0.1,False)
  
  reg_uwnd_ERA_AMJ = mask(reg_uwnd_ERA_AMJ,reg_uwnd_ERA_AMJ.lt.0.1 .and. reg_uwnd_ERA_AMJ.gt.-0.1.and.reg_vwnd_ERA_AMJ.lt.0.1 .and. reg_vwnd_ERA_AMJ.gt.-0.1,False)
  reg_uwnd_ERA_MAM = mask(reg_uwnd_ERA_MAM,reg_uwnd_ERA_MAM.lt.0.1 .and. reg_uwnd_ERA_MAM.gt.-0.1.and.reg_vwnd_ERA_MAM.lt.0.1 .and. reg_vwnd_ERA_MAM.gt.-0.1,False)
  reg_uwnd_ERA_JJA = mask(reg_uwnd_ERA_JJA,reg_uwnd_ERA_JJA.lt.0.1 .and. reg_uwnd_ERA_JJA.gt.-0.1.and.reg_vwnd_ERA_JJA.lt.0.1 .and. reg_vwnd_ERA_JJA.gt.-0.1,False)
  reg_uwnd_ERA_MJJ = mask(reg_uwnd_ERA_MJJ,reg_uwnd_ERA_MJJ.lt.0.1 .and. reg_uwnd_ERA_MJJ.gt.-0.1.and.reg_vwnd_ERA_MJJ.lt.0.1 .and. reg_vwnd_ERA_MJJ.gt.-0.1,False)
                                                                                                                                    
  reg_vwnd_ERA_AMJ = mask(reg_vwnd_ERA_AMJ,reg_uwnd_ERA_AMJ.lt.0.1 .and. reg_uwnd_ERA_AMJ.gt.-0.1.and.reg_vwnd_ERA_AMJ.lt.0.1 .and. reg_vwnd_ERA_AMJ.gt.-0.1,False)
  reg_vwnd_ERA_MAM = mask(reg_vwnd_ERA_MAM,reg_uwnd_ERA_MAM.lt.0.1 .and. reg_uwnd_ERA_MAM.gt.-0.1.and.reg_vwnd_ERA_MAM.lt.0.1 .and. reg_vwnd_ERA_MAM.gt.-0.1,False)
  reg_vwnd_ERA_JJA = mask(reg_vwnd_ERA_JJA,reg_uwnd_ERA_JJA.lt.0.1 .and. reg_uwnd_ERA_JJA.gt.-0.1.and.reg_vwnd_ERA_JJA.lt.0.1 .and. reg_vwnd_ERA_JJA.gt.-0.1,False)
  reg_vwnd_ERA_MJJ = mask(reg_vwnd_ERA_MJJ,reg_uwnd_ERA_MJJ.lt.0.1 .and. reg_uwnd_ERA_MJJ.gt.-0.1.and.reg_vwnd_ERA_MJJ.lt.0.1 .and. reg_vwnd_ERA_MJJ.gt.-0.1,False)
  
  plot_vc = new(8,graphic)
    
  res_vc@vcRefMagnitudeF             = 0.5
  
  res_vc@gsnRightStringFontHeightF   = 0.028
  res_vc@gsnLeftStringFontHeightF    = 0.028
  res_vc@gsnCenterStringFontHeightF  = 0.028
  
  res_vc@tmXBLabelFontHeightF = 0.028
  res_vc@tmYLLabelFontHeightF = 0.028
  
  res_vc@tmXBOn               = False
  res_vc@tmYLOn               = True
  res_vc@vcLineArrowColor        = "black"
  res_vc@vcMonoLineArrowColor    = True
  res_vc@gsnLeftString               = "(a)"
  res_vc@gsnRightString              = "MAM "
  res_vc@gsnCenterString             = "NCEP "
  plot_vc(0)  = gsn_csm_vector_map_ce(wks,reg_uwnd_MAM(::2,::2),reg_vwnd_MAM(::2,::2),res_vc)
  res_vc@gsnLeftString               = "(c)"
  res_vc@gsnRightString              = "AMJ "
  plot_vc(2)  = gsn_csm_vector_map_ce(wks,reg_uwnd_AMJ(::2,::2),reg_vwnd_AMJ(::2,::2),res_vc)
  res_vc@gsnLeftString               = "(e)"
  res_vc@gsnRightString              = "MJJ "
  plot_vc(4)  = gsn_csm_vector_map_ce(wks,reg_uwnd_MJJ(::2,::2),reg_vwnd_MJJ(::2,::2),res_vc)
  res_vc@gsnLeftString               = "(g)"
  res_vc@gsnRightString              = "JJA "
  res_vc@tmXBOn               = True
  plot_vc(6)  = gsn_csm_vector_map_ce(wks,reg_uwnd_JJA(::2,::2),reg_vwnd_JJA(::2,::2),res_vc)
  
  res_vc@gsnLeftString               = "(b)"
  res_vc@gsnRightString              = "MAM "
  res_vc@gsnCenterString             = "ERA "
  res_vc@tmXBOn               = False
  res_vc@tmYLOn               = False
  plot_vc(1)  = gsn_csm_vector_map_ce(wks,reg_uwnd_ERA_MAM(::4,::4),reg_vwnd_ERA_MAM(::4,::4),res_vc)
  res_vc@gsnLeftString               = "(d)"
  res_vc@gsnRightString              = "AMJ "
  plot_vc(3)  = gsn_csm_vector_map_ce(wks,reg_uwnd_ERA_AMJ(::4,::4),reg_vwnd_ERA_AMJ(::4,::4),res_vc)
  res_vc@gsnLeftString               = "(f)"
  res_vc@gsnRightString              = "MJJ "
  plot_vc(5)  = gsn_csm_vector_map_ce(wks,reg_uwnd_ERA_MJJ(::4,::4),reg_vwnd_ERA_MJJ(::4,::4),res_vc)
  res_vc@gsnLeftString               = "(h)"
  res_vc@gsnRightString              = "JJA "
  res_vc@tmXBOn               = True
  plot_vc(7)  = gsn_csm_vector_map_ce(wks,reg_uwnd_ERA_JJA(::4,::4),reg_vwnd_ERA_JJA(::4,::4),res_vc)
  
  plot_test     = new(8,graphic)
  
  res_vc@vcLineArrowColor        = "blue"
  res_vc@vcMonoLineArrowColor    = True
  
  ;===vector test
  reg_uwnd_AMJ = mask(reg_uwnd_AMJ,common_test_AMJ.lt.2. .and. common_test_AMJ.gt.-2.,False)
  reg_uwnd_MAM = mask(reg_uwnd_MAM,common_test_MAM.lt.2. .and. common_test_MAM.gt.-2.,False)
  reg_uwnd_JJA = mask(reg_uwnd_JJA,common_test_JJA.lt.2. .and. common_test_JJA.gt.-2.,False)
  reg_uwnd_MJJ = mask(reg_uwnd_MJJ,common_test_MJJ.lt.2. .and. common_test_MJJ.gt.-2.,False)
                                                      
  reg_vwnd_AMJ = mask(reg_vwnd_AMJ,common_test_AMJ.lt.2. .and. common_test_AMJ.gt.-2.,False)
  reg_vwnd_MAM = mask(reg_vwnd_MAM,common_test_MAM.lt.2. .and. common_test_MAM.gt.-2.,False)
  reg_vwnd_JJA = mask(reg_vwnd_JJA,common_test_JJA.lt.2. .and. common_test_JJA.gt.-2.,False)
  reg_vwnd_MJJ = mask(reg_vwnd_MJJ,common_test_MJJ.lt.2. .and. common_test_MJJ.gt.-2.,False)
  
  reg_uwnd_ERA_AMJ = mask(reg_uwnd_ERA_AMJ,common_test_ERA_AMJ.lt.2. .and. common_test_ERA_AMJ.gt.-2.,False)
  reg_uwnd_ERA_MAM = mask(reg_uwnd_ERA_MAM,common_test_ERA_MAM.lt.2. .and. common_test_ERA_MAM.gt.-2.,False)
  reg_uwnd_ERA_JJA = mask(reg_uwnd_ERA_JJA,common_test_ERA_JJA.lt.2. .and. common_test_ERA_JJA.gt.-2.,False)
  reg_uwnd_ERA_MJJ = mask(reg_uwnd_ERA_MJJ,common_test_ERA_MJJ.lt.2. .and. common_test_ERA_MJJ.gt.-2.,False)
                                                      
  reg_vwnd_ERA_AMJ = mask(reg_vwnd_ERA_AMJ,common_test_ERA_AMJ.lt.2. .and. common_test_ERA_AMJ.gt.-2.,False)
  reg_vwnd_ERA_MAM = mask(reg_vwnd_ERA_MAM,common_test_ERA_MAM.lt.2. .and. common_test_ERA_MAM.gt.-2.,False)
  reg_vwnd_ERA_JJA = mask(reg_vwnd_ERA_JJA,common_test_ERA_JJA.lt.2. .and. common_test_ERA_JJA.gt.-2.,False)
  reg_vwnd_ERA_MJJ = mask(reg_vwnd_ERA_MJJ,common_test_ERA_MJJ.lt.2. .and. common_test_ERA_MJJ.gt.-2.,False)
  
  ;===
  delete(res_vc@mpCenterLonF)
  delete(res_vc@mpMaxLatF   )
  delete(res_vc@mpMinLatF   )
  delete(res_vc@mpMaxLonF   )
  delete(res_vc@mpMinLonF   )
  
  res_vc@gsnLeftString               = ""
  res_vc@gsnRightString              = ""
  res_vc@gsnCenterString             = ""
  
  res_vc@vcRefAnnoOn = "False"
  
  plot_test(0)  = gsn_csm_vector(wks,reg_uwnd_MAM(::2,::2),reg_vwnd_MAM(::2,::2),res_vc)
  plot_test(2)  = gsn_csm_vector(wks,reg_uwnd_AMJ(::2,::2),reg_vwnd_AMJ(::2,::2),res_vc)
  plot_test(4)  = gsn_csm_vector(wks,reg_uwnd_MJJ(::2,::2),reg_vwnd_MJJ(::2,::2),res_vc)
  plot_test(6)  = gsn_csm_vector(wks,reg_uwnd_JJA(::2,::2),reg_vwnd_JJA(::2,::2),res_vc)
  
  plot_test(1)  = gsn_csm_vector(wks,reg_uwnd_ERA_MAM(::4,::4),reg_vwnd_ERA_MAM(::4,::4),res_vc)
  plot_test(3)  = gsn_csm_vector(wks,reg_uwnd_ERA_AMJ(::4,::4),reg_vwnd_ERA_AMJ(::4,::4),res_vc)
  plot_test(5)  = gsn_csm_vector(wks,reg_uwnd_ERA_MJJ(::4,::4),reg_vwnd_ERA_MJJ(::4,::4),res_vc)
  plot_test(7)  = gsn_csm_vector(wks,reg_uwnd_ERA_JJA(::4,::4),reg_vwnd_ERA_JJA(::4,::4),res_vc)
  
  do i=0,7
  overlay(plot_vc(i),plot_test(i))
  end do  
  ;====================================
  rest                   = True                      ; polyline mods desired
  rest@gsLineColor       = "red"                    ; color of lines
  rest@gsLineThicknessF  = 2                       ; thickness of lines
  rest@gsLineDashPattern = 0
  
  ;===
  ypts_A = (/ 5.0,  5.0,  -5.0,  -5.0, 5.0/)      ;latitude
  xpts_A = (/ 40.0,  50.0,  50.0,   40.0,  40.0/)    ;longitude
  
  ypts_B = (/ 70.0,  70.0,  40.0,  40.0, 70.0/)      ;latitude
  xpts_B = (/ 110.0,  140.0,  140.0,   110.0,  110.0/)    ;longitude
  
  dum_A1 = new(4,graphic)
  dum_B1 = new(4,graphic)
  dum_C1 = new(4,graphic)
  dum_D1 = new(4,graphic)
  
  dum_A2 = new(4,graphic)
  dum_B2 = new(4,graphic)
  dum_C2 = new(4,graphic)
  dum_D2 = new(4,graphic)
  
  do i = 0, 3
   dum_A1(i)=gsn_add_polyline(wks,plot_vc(0),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_B1(i)=gsn_add_polyline(wks,plot_vc(1),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_C1(i)=gsn_add_polyline(wks,plot_vc(2),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_D1(i)=gsn_add_polyline(wks,plot_vc(3),xpts_A(i:i+1),ypts_A(i:i+1),rest)
   
   dum_A2(i)=gsn_add_polyline(wks,plot_vc(4),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_B2(i)=gsn_add_polyline(wks,plot_vc(5),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_C2(i)=gsn_add_polyline(wks,plot_vc(6),xpts_A(i:i+1),ypts_A(i:i+1),rest) 
   dum_D2(i)=gsn_add_polyline(wks,plot_vc(7),xpts_A(i:i+1),ypts_A(i:i+1),rest)
  end do
  
  ;===
  ;===
  resp = True
  
  resp@gsnPanelLabelBar         = False
  resp@pmLabelBarOrthogonalPosF = 0.01
  
  resp@pmLabelBarWidthF         = 0.55
  resp@pmLabelBarHeightF        = 0.08
  resp@lbLabelFontHeightF       = 0.015
  
  resp@gsnPanelBottom          = 0.06
  
  resp@lbBoxLineColor           = 0
  resp@lbLabelAutoStride        = True
  resp@lbBoxLineThicknessF      = 3
  
  resp@gsnPanelYF               = (/0.95,0.95,0.72,0.72,0.49,0.49,0.26,0.26/);-0.05
  resp@gsnPanelXF               = (/0.13,0.5,0.13,0.5,0.13,0.5,0.13,0.5/);-0.05
  
  gsn_panel(wks,plot_vc,(/4,2/),resp)
  
end