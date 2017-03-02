#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

#setenv NCARG_ROOT /raid/cycles/GEXCEL/chengw/bin/ncl_ncarg-5.1.0
#setenv NCARG_LIB  /raid/cycles/GEXCEL/chengw/bin/ncl_ncarg-5.1.0/lib

set filelist_wrf_d03 = ` ls -1 wrfout_d03* `

foreach ff3 ( $filelist_wrf_d03 )
 set year  = ` echo $ff3 | cut -c12-15 `
 set month = ` echo $ff3 | cut -c17-18 `
 set day   = ` echo $ff3 | cut -c20-21 `
 set hour  = ` echo $ff3 | cut -c23-24 `
 if ( $hour == 00 || $hour == 03 || $hour == 06 || $hour == 09 || $hour == 12 || $hour == 15 || $hour == 18 || $hour == 21 ) then
  set time_index = ` echo $ff3 | cut -c12-30 ` 
  echo $ff3 $year $month $day $hour $time_index
  set file_prefix = ` echo $ff3 | cut -c1-30 `

 # ln -s $ff3 ${ff3}.nc
  rm -f test3d.gdat

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
  file_out = "test3d.gdat"
   
  fall      = addfile (file_in, "r")  
  ;========================

  tlat1     = fall@TRUELAT1
  tlat2     = fall@TRUELAT2
  clat      = fall@MOAD_CEN_LAT
  clon      = fall@STAND_LON
  dx        = fall@DX

  it =0 

  ;xlat = wrf_user_getvar(fall,"lat",it) 
  ;xlon = wrf_user_getvar(fall,"lon",it) 
  ;hgt  = wrf_user_getvar(fall,"ter",it)    ; terrain

  xlat_3d = fall->XLAT
  xlon_3d = fall->XLONG
  hgt_3d  = fall->HGT
  pblh    = fall->PBLH

  xlat = xlat_3d(0,:,:)
  xlon = xlon_3d(0,:,:)
  hgt  = hgt_3d(0,:,:)

  slp = wrf_user_getvar(fall,"slp",it)    ; slp
   wrf_smooth_2d( slp, 3 )            ; smooth slp

  tc2 = wrf_user_getvar(fall,"T2",it)     ; T2 in Kelvin
       tc2 = tc2-273.16                  ; T2 in C
       tc2@units = "C"

  u10 = wrf_user_getvar(fall,"U10",it)    ; u at 10 m, mass point
  v10 = wrf_user_getvar(fall,"V10",it)    ; v at 10 m, mass point

  tc = wrf_user_getvar(fall,"tc",it)      ; 3D tc
  ;td = wrf_user_getvar(fall,"td",it)      ; 3D td
  uvmet = wrf_user_getvar(fall,"uvmet",it)

  dsizes_tc = dimsizes(tc)
  n_z = dsizes_tc(0)
  n_y = dsizes_tc(1) 
  n_x = dsizes_tc(2)

  print(n_x)
  print(n_y)
  print(n_z)

  printVarSummary(uvmet)


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

  fbindirwrite(file_out,xlat)
  fbindirwrite(file_out,xlon)
  fbindirwrite(file_out,hgt)

  fbindirwrite(file_out,tc2)
  fbindirwrite(file_out,u10)
  fbindirwrite(file_out,v10)
  fbindirwrite(file_out,slp)
  fbindirwrite(file_out,pblh(0,:,:))  

  fbindirwrite(file_out,tc(0,:,:))
  fbindirwrite(file_out,uvmet(0,0,:,:))
  fbindirwrite(file_out,uvmet(1,0,:,:))

  end
EOF

ncl test.ncl
mv test3d.gdat eta_level_$year$month$day${hour}.gdat

 endif

end
