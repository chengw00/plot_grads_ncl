#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

set DOMAIN = $1
set MEMBER = $2
set NCARG_ROOT = $3

if ( $1 == "" || $2 == "" || $3 == ""  ) then
 echo "DOMAIN NOT DEFINED"
 echo "or MEMBER not defined"
 echo "or NCARG_ROOT not defined"
 exit
endif

#setenv NCARG_ROOT /raid/cycles/GEXCEL/chengw/bin/ncl_ncarg-5.1.0
#setenv NCARG_LIB  /raid/cycles/GEXCEL/chengw/bin/ncl_ncarg-5.1.0/lib

set filelist_wrf = ` ls -1 wrfout_${DOMAIN}* `

foreach ff3 ( $filelist_wrf )
 set year  = ` echo $ff3 | cut -c12-15 `
 set month = ` echo $ff3 | cut -c17-18 `
 set day   = ` echo $ff3 | cut -c20-21 `
 set hour  = ` echo $ff3 | cut -c23-24 `
 #if ( $hour == 00 || $hour == 03 || $hour == 06 || $hour == 09 || $hour == 12 || $hour == 15 || $hour == 18 || $hour == 21 ) then
  set time_index = ` echo $ff3 | cut -c12-24 ` 
  echo $ff3 $year $month $day $hour $time_index
  set file_prefix = ` echo $ff3 | cut -c1-30 `

 # ln -s $ff3 ${ff3}.nc
  rm -f all*dat 
  rm -f terr_mapfc_latlon_wrf_realtime.gdat 
  rm -f wrfout_2d_realtime.gdat 
  rm -f wrfout_3d_realtime.gdat

cat <<EOF>test.ncl
;********************************************************
; WRF: latitude-z cross section.
;********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
;load "/home/chengw/bin/ncl4.2.0/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "/home/chengw/bin/ncl4.2.0/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "/usr/local/ncarg-5.0.0/lib/ncarg/nclscripts/csm/contributed.ncl"   
;load "./roux/StationModel.ncl"
begin
  ; ==================
  setfileoption("bin","WriteByteOrder","LittleEndian")
  file_in = "${file_prefix}.nc"
  file_ter = "terr_mapfc_latlon_wrf_realtime.gdat"
  file_2d = "wrfout_2d_realtime.gdat"
  file_3d = "wrfout_3d_realtime.gdat"
   
  fall      = addfile (file_in, "r")  
  ;========================

  tlat1     = fall@TRUELAT1
  tlat2     = fall@TRUELAT2
  clat      = fall@MOAD_CEN_LAT
  clon      = fall@STAND_LON
  dx        = fall@DX

  it =0 

  ; ====== set 1 ================
  ;xlat = wrf_user_getvar(fall,"lat",it) 
  ;xlon = wrf_user_getvar(fall,"lon",it) 
  ;hgt  = wrf_user_getvar(fall,"ter",it)    ; terrain

  xlat_3d = fall->XLAT
  xlon_3d = fall->XLONG
  hgt_3d  = fall->HGT
  pblh    = fall->PBLH
  mapf    = fall->MAPFAC_M

  xlat = xlat_3d(0,:,:)
  xlon = xlon_3d(0,:,:)
  hgt  = hgt_3d(0,:,:)

  ; ====== set 2 ==========
  qv2    = fall->Q2
  rainc  = fall->RAINC
  rainnc = fall->RAINNC
  th2    = fall->TH2
  psfc   = fall->PSFC

  hfx    = fall->HFX
  lh     = fall->LH

  psfc   = psfc/100.    ; convert to hPa

  slp = wrf_user_getvar(fall,"slp",it)    ; slp (Pa)
  u10 = wrf_user_getvar(fall,"U10",it)    ; u at 10 m, mass point
  v10 = wrf_user_getvar(fall,"V10",it)    ; v at 10 m, mass point

  ; ====== set 3 ========
  ua = wrf_user_getvar(fall,"ua",it)
  va = wrf_user_getvar(fall,"va",it)
  wa = wrf_user_getvar(fall,"wa",it)
  tc = wrf_user_getvar(fall,"tc",it)      ; 3D tc
  pressure = wrf_user_getvar(fall,"pressure",it) ; (hPa)
  height = wrf_user_getvar(fall,"z",it)

  qv = fall->QVAPOR
  qcloud = fall->QCLOUD
  qrain  = fall->QRAIN
  qice   = fall->QICE
  qsnow  = fall->QSNOW
  qgraup = fall->QGRAUP

  qhydromet = qcloud + qrain + qice + qsnow + qgraup
  qhydromet_liq = qcloud + qrain
  
  pb     = fall->PB

  dsizes_tc = dimsizes(tc)
  n_z = dsizes_tc(0)
  n_y = dsizes_tc(1) 
  n_x = dsizes_tc(2)

  print(n_x)
  print(n_y)
  print(n_z)


  ; ====== define parameters for GrADS control file ======
  xlat_sw = xlat(0,0)
  xlon_sw = xlon(0,0)

  xlat_se = xlat(0,n_x-1)
  xlon_se = xlon(0,n_x-1)

  xlat_ne = xlat(n_y-1,n_x-1)
  xlon_ne = xlon(n_y-1,n_x-1)

  xlat_nw = xlat(n_y-1,0)
  xlon_nw = xlon(n_y-1,0)

  print(xlat_sw)
  print(xlon_sw)

  xlon_sw_grads = 0.5*(xlon_sw+xlon_nw)
  dlon_grads    = 0.5*(xlon_ne+xlon_se-xlon_nw-xlon_sw)/(n_x - 1)

  print(xlon_sw_grads)
  print(dlon_grads)

  xlat_sw_grads = 0.5*(xlat_sw+xlat_se)
  dlat_grads    = 0.5*(xlat_nw+xlat_ne-xlat_sw-xlat_se)/(n_y - 1)

  print(xlat_sw_grads)
  print(dlat_grads)

  ; ======== write out fields in GrADS format ===========

  fbindirwrite(file_ter,hgt)
  fbindirwrite(file_ter,mapf(0,:,:))
  fbindirwrite(file_ter,xlat)
  fbindirwrite(file_ter,xlon)

  fbindirwrite(file_2d,qv2(0,:,:))
  fbindirwrite(file_2d,rainc(0,:,:))
  fbindirwrite(file_2d,rainnc(0,:,:))
  fbindirwrite(file_2d,slp)
  fbindirwrite(file_2d,th2(0,:,:))
  fbindirwrite(file_2d,u10)
  fbindirwrite(file_2d,v10)
  fbindirwrite(file_2d,psfc(0,:,:))
  fbindirwrite(file_2d,pblh(0,:,:))  
  fbindirwrite(file_2d,hfx(0,:,:))
  fbindirwrite(file_2d,lh(0,:,:))


  fbindirwrite(file_3d,ua)
  fbindirwrite(file_3d,va)
  fbindirwrite(file_3d,wa)
  fbindirwrite(file_3d,tc)
  fbindirwrite(file_3d,pressure)
  fbindirwrite(file_3d,height)
  fbindirwrite(file_3d,qv(0,:,:,:))
  fbindirwrite(file_3d,qhydromet(0,:,:,:))
  fbindirwrite(file_3d,qhydromet_liq(0,:,:,:))
  fbindirwrite(file_3d,pb(0,:,:,:))

  end
EOF

ncl test.ncl
diag_wrf_v2.exe 
mv all.dat press_level_${time_index}.gdat
mv all_zetap.dat zetap_level_${time_index}.gdat
mv all_2d.dat store_2d_${time_index}_${DOMAIN}.gdat
mv all_3d.dat store_3d_${time_index}_${DOMAIN}.gdat
mv fort.201 pbl_data_${time_index}_${MEMBER}_${DOMAIN}.dat

 #endif
#plot_tcp_ncl.csh ${DOMAIN} ${time_index} 18 $NCARG_ROOT

#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT 1000
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  975
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  950
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  925
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  900
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  850
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  700 
#plot_rh_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  500

#plot_rh_2m_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT 2m

##plot_terrain_with_rh_2m_90contour.csh $DOMAIN $time_index 1 $NCARG_ROOT 2m
#plot_terrain_with_rh_2m_70contour_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT 2m

#plot_terrain_with_tcp_contour_ncl.csh ${DOMAIN} ${time_index} 18 $NCARG_ROOT 2m

#plot_tmp_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT 1000
#plot_tmp_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  850
#plot_tmp_plev_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT  700

#plot_t2m_ncl.csh $DOMAIN $time_index 1 $NCARG_ROOT 2m

end

#rm -f *.gdat
