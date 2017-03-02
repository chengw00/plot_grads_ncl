#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts


set filelist = ` ls -1 wrfout_d03_2009-04-17_00:00:00.nc `
set lat_in = 39.78 
set lon_in = -104.86

foreach ff ( $filelist )
 set time_index  = ` echo $ff | cut -c12-24 `
 set file_prefix = ` echo $ff | cut -c1-30 `
 
 echo $ff $file_prefix
 
 rm -f sounding.gdat

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
  file_out = "sounding.gdat"

  fall      = addfile (file_in, "r")
  ;========================
  it = 0
  
  tc = wrf_user_getvar(fall,"tc",it)      ; 3D tc
  pressure = wrf_user_getvar(fall,"pressure",it) ; (hPa)
  height = wrf_user_getvar(fall,"z",it)
  uvmet  = wrf_user_getvar(fall,"uvmet",it)
  td     = wrf_user_getvar(fall,"td",it)

  umet   = uvmet(0,:,:,:)
  vmet   = uvmet(1,:,:,:)

  loc = wrf_user_latlon_to_ij (fall, ${lat_in}, ${lon_in})
  print(loc)


  ; ======== output profile ======
  fbindirwrite(file_out,pressure(:,loc(0),loc(1)))
  fbindirwrite(file_out,height(:,loc(0),loc(1)))
  fbindirwrite(file_out,tc(:,loc(0),loc(1)))
  fbindirwrite(file_out,umet(:,loc(0),loc(1)))
  fbindirwrite(file_out,vmet(:,loc(0),loc(1)))
  fbindirwrite(file_out,td(:,loc(0),loc(1)))


end
EOF
ncl test.ncl
rm -f all.dat
rewrite_sounding.exe
mv all.dat sounding.gdat

cat <<EOF>part1.ctl
dset sounding.gdat
options little_endian
undef -888888
xdef  1 linear -126.31445 0.14846323
ydef  1 linear   30.19919 0.1101682
EOF

cat <<EOF>part2.ctl
tdef 1 linear 12z12dec2008 1mn
vars 6
p            36 99   pressure (hPa)
z            36 99   hgt (m)
tempc        36 99   temp (deg C)
windspeed    36 99   wind speed (m/s)
windmeteodir 36 99   wind dir (deg)
dewptc       36 99   dewp (deg C)
endvars
EOF

cat part1.ctl fort.99 part2.ctl > sounding.ctl

cat <<EOF>comm.scr
'open sounding.ctl '
*'set lev 900 100'
'run plotskew_aug2002.gs'
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif sounding_${time_index}.gif

end
