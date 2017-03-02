#!/bin/csh

set NCARG_ROOT = $1
set NCARG_LIB  = $2
set DOMAIN =     $3

setenv NCARG_ROOT $NCARG_ROOT
setenv NCARG_LIB  $NCARG_LIB

#set DOMAIN = ( d01 d02 d03 )
#set DOMAIN = ( d03 )

foreach dd ( $DOMAIN )
 rm -r -f $dd
 mkdir $dd
 set filelist = `  ls -1 wrfout*${dd}*:00:00*  `
 foreach ff1 ( $filelist )
  set time_index = ` echo $ff1 | cut -c12-24 `
  set ff         = ` echo $ff1 | cut -c1-30 `
  set hour       = ` echo $ff1 | cut -c23-24 `

  if ( $hour == 00 || $hour == 03 || $hour == 06 || $hour == 09 || $hour == 12 || \
       $hour == 15 || $hour == 18 || $hour == 21 ) then 
   echo $dd $ff $time_index

cat <<EOF> test.ncl
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
  file_in = "${ff}.nc"
  file_out =     "${dd}/pcpac_sigz0000_${time_index}.gdat"
  file_tmp_out = "${dd}/tempc_plev0850_${time_index}.gdat"

  fall      = addfile (file_in, "r")
  system("rm -f ${dd}/pcpac_sigz0000_${time_index}.gdat  ")
  system("rm -f ${dd}/tempc_plev0850_${time_index}.gdat  ")

  plev = 85000   ; pressure level to be interpolated (Pa)
  ;========================

  tlat1     = fall@TRUELAT1
  tlat2     = fall@TRUELAT2
  clat      = fall@MOAD_CEN_LAT
  clon      = fall@STAND_LON
  dx        = fall@DX

  print(dx)

  hgt  =    fall->HGT
  xlat =    fall->XLAT
  xlon =    fall->XLONG

  theta_p = fall->T
  theta_add_300 = theta_p(0,:,:,:) + 300
  ;theta_3d = wrf_user_getvar(fall,"theta",0)      ; 3D t

  rainc     = fall->RAINC
  rainnc    = fall->RAINNC
  accum_tot = rainc+rainnc
  
  ;======== calculate height =========== 
  ph        = fall->PH
  phb       = fall->PHB

  z = (ph + phb)/9.81

  ; ======= get pressure (Pa) ========
  p         = fall->P
  pb        = fall->PB
  
  press = p+pb
  ; ======= get dimensions =========
  znu = fall->ZNU
  znw = fall->ZNW

  n_znu = dimsizes(znu)
  n_znw = dimsizes(znw)
  n_z_half = n_znu(1)
  n_z_full = n_znw(1)

  dsizes_theta_p = dimsizes(theta_p)
  n_z = dsizes_theta_p(1)
  n_y = dsizes_theta_p(2)
  n_x = dsizes_theta_p(3)

  ; ====== define parameters for GrADS control file ======
  xlat_sw = xlat(0,0,0)
  xlon_sw = xlon(0,0,0)

  xlat_se = xlat(0,0,n_x-1)
  xlon_se = xlon(0,0,n_x-1)

  xlat_ne = xlat(0,n_y-1,n_x-1)
  xlon_ne = xlon(0,n_y-1,n_x-1)

  xlat_nw = xlat(0,n_y-1,0)
  xlon_nw = xlon(0,n_y-1,0)

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

  ; ======== create new arrays ======
  z_half_eta = z(0,:,:,:)
  z_half_eta@_FillValue = -999999
  printVarSummary(z_half_eta)

  print(n_z_full)

  do k=0,n_z_full-2
   z_half_eta(k,:,:)=0.5*(z(0,k,:,:)+z(0,k+1,:,:))
  end do
  z_half_eta(n_z_full-1,:,:)=-999999
  printVarSummary(z_half_eta)


  ; ======= interpolate height to pressure levels ====
  theta_press   = new( (/n_y, n_x/), float)
  theta_press@_FillValue = -999999

  theta_press = wrf_user_intrp3d(theta_add_300,press(0,:,:,:),"h",plev,0.,False)
  tmpc_press = theta_press*(plev/100000.)^0.287

  tmpc_press = tmpc_press - 273.15

  ; ============================================================
  fbindirwrite(file_out,accum_tot(0,:,:))

  fbindirwrite(file_tmp_out,tmpc_press)
  ;fbindirwrite(file_tmp_out,theta_3d(1,:,:))

end
EOF

${NCARG_ROOT}/bin/ncl test.ncl

  endif

 end
end

