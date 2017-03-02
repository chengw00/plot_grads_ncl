#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

set DOMAIN = $1
set MEMBER = $2
if ( $1 == "" || $2 == ""  ) then
 echo "DOMAIN NOT DEFINED"
 echo "or MEMBER not defined"
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
 if ( $DOMAIN == d01 || $DOMAIN == d02 || $DOMAIN == d03 ) then  # replace the above with this line
  set time_index = ` echo $ff3 | cut -c12-24 ` 
  echo $ff3 $year $month $day $hour $time_index
  set file_prefix = ` echo $ff3 | cut -c1-30 `

if ( $DOMAIN == d01 ) then

cat <<EOF>press_level_${time_index}.ctl
dset press_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  279  194 lcc 9.226906 65.54395  1.  1.  10.0  40.0 100.0 27000. 27000.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  26 levels 1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50 30 20 10
tdef 1 linear 00z04apr2003 1hr
vars 34
HGTsfc     1 99  terrain (m)
ptrop      1 99  dynamic tropopause pressure (hPa)
utrop      1 99  zonal wind on dynamic tropopause (m/s)
vtrop      1 99  meridional wind on dynamic tropopause (m/s)
slpmm5     1 99  sea-level pressure (hPa)
slpwrf     1 99  sea-level pressure (hPa)
psfc       1 99  surface pressure (hPa)
cape       1 99  CAPE (J/kg)
cin        1 99  CIN (J/kg)
acc_conv   1 99  accumulated subgrid-scale precip (mm)
acc_grid   1 99  accumulated grid-scale precip (mm)
tc2        1 99  temperature at 2 m AGL (deg C)
ept2       1 99  equiv. potential temp at 2 m AGL (K)
rh2        1 99  relative hum at 2 m AGL
u10        1 99  zonal wind at 2 m AGL (m/s)
v10        1 99  meridional wind at 2 m AGL (m/s)
ppw        1 99  precipitable water (mm)
hyd_path   1 99  total hydrometeor path (mm)
hyd_path_l 1 99 liquid hydrometeor path (mm)
cltopt     1 99  cloud top temperature (deg C)
pblh       1 99  PBL height (m)
flux_lh    1 99  sfc latent heat flux (W/m**2)
flux_sh    1 99  sfc sensible heat flux (W/m**2)
u         26 99  zonal wind at pressure levels (m/s)
v         26 99  meridional wind at pressure levels (m/s)
w         26 99  vertical velocity (m/s)
z         26 99  geopotential height (m)
tempc     26 99  temperature (deg C)
rh        26 99  relative humidity (%)
qhydro    26 99  total hydrometeor mixing ratio (kg/kg)
qhydro_l  26 99  liquid hydrometeor mixing ratio (kg/kg)
pv        26 99  potential vorticity (PVU)
relvort   26 99  relative vorticity (10^-5 s^-1)
absvort   26 99  absolute vorticity (10^-5 s^-1)
endvars
EOF

cat <<EOF>zetap_level_${time_index}.ctl
dset zetap_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  159  169 lcc  18.95885 109.7372  1.  1.  30.0  60.0 127.0 24300. 24300.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  5 levels 1 2 3 4 5
tdef 1 linear 00z04apr2003 1hr
vars 9
z          5 99 height (m)
press      5 99 pressure (hPa)
tempc      5 99 temperature (deg C)
rh         5 99 RH (%)
u          5 99 zonal wind (m/s)
v          5 99 meridional wind (m/s)
w          5 99 vertical velocity (m/s)
qh_tot     5 99 total hydrometeor mr (kg/kg)
qh_l       5 99 liquid hydrometeor mr (kg/kg)
endvars
EOF



set NBARB = 30

else if ( $DOMAIN == d02 ) then

cat <<EOF>press_level_${time_index}.ctl
dset press_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  651  483 lcc 15.52342 73.7645  1.  1.  10.0  40.0 100.0 9000. 9000.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  26 levels 1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50 30 20 10
tdef 1 linear 00z04apr2003 1hr
vars 34
HGTsfc     1 99  terrain (m) 
ptrop      1 99  dynamic tropopause pressure (hPa)
utrop      1 99  zonal wind on dynamic tropopause (m/s)
vtrop      1 99  meridional wind on dynamic tropopause (m/s)
slpmm5     1 99  sea-level pressure (hPa)
slpwrf     1 99  sea-level pressure (hPa)
psfc       1 99  surface pressure (hPa)
cape       1 99  CAPE (J/kg)
cin        1 99  CIN (J/kg)
acc_conv   1 99  accumulated subgrid-scale precip (mm)
acc_grid   1 99  accumulated grid-scale precip (mm)
tc2        1 99  temperature at 2 m AGL (deg C)
ept2       1 99  equiv. potential temp at 2 m AGL (K) 
rh2        1 99  relative hum at 2 m AGL 
u10        1 99  zonal wind at 2 m AGL (m/s)
v10        1 99  meridional wind at 2 m AGL (m/s)
ppw        1 99  precipitable water (mm)
hyd_path   1 99  total hydrometeor path (mm)
hyd_path_l 1 99 liquid hydrometeor path (mm)
cltopt     1 99  cloud top temperature (deg C)
pblh       1 99  PBL height (m) 
flux_lh    1 99  sfc latent heat flux (W/m**2)
flux_sh    1 99  sfc sensible heat flux (W/m**2)
u         26 99  zonal wind at pressure levels (m/s)
v         26 99  meridional wind at pressure levels (m/s)
w         26 99  vertical velocity (m/s)
z         26 99  geopotential height (m)
tempc     26 99  temperature (deg C)
rh        26 99  relative humidity (%)
qhydro    26 99  total hydrometeor mixing ratio (kg/kg)
qhydro_l  26 99  liquid hydrometeor mixing ratio (kg/kg)
pv        26 99  potential vorticity (PVU)
relvort   26 99  relative vorticity (10^-5 s^-1)
absvort   26 99  absolute vorticity (10^-5 s^-1)
endvars
EOF

cat <<EOF>zetap_level_${time_index}.ctl
dset zetap_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  651  483 lcc 15.52342 73.7645  1.  1.  10.0  40.0 100.0 9000. 9000.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  5 levels 1 2 3 4 5
tdef 1 linear 00z04apr2003 1hr
vars 9
z          5 99 height (m)
press      5 99 pressure (hPa)
tempc      5 99 temperature (deg C)
rh         5 99 RH (%)
u          5 99 zonal wind (m/s)
v          5 99 meridional wind (m/s)
w          5 99 vertical velocity (m/s)
qh_tot     5 99 total hydrometeor mr (kg/kg)
qh_l       5 99 liquid hydrometeor mr (kg/kg)
endvars
EOF

set NBARB = 30

else if ( $DOMAIN == d03 ) then

cat <<EOF>press_level_${time_index}.ctl
dset press_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  513 486 lcc  29.85936 117.885 1.  1.  30.00000  60.00000 127.0 2700. 2700.
xdef 513 linear   116.9581  0.0309388
ydef 486 linear    30.00962 0.02463833
zdef  26 levels 1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100 70 50 30 20 10
tdef 1 linear 00z04apr2003 1hr
vars 34
HGTsfc     1 99  terrain (m)
ptrop      1 99  dynamic tropopause pressure (hPa)
utrop      1 99  zonal wind on dynamic tropopause (m/s)
vtrop      1 99  meridional wind on dynamic tropopause (m/s)
slpmm5     1 99  sea-level pressure (hPa)
slpwrf     1 99  sea-level pressure (hPa)
psfc       1 99  surface pressure (hPa)
cape       1 99  CAPE (J/kg)
cin        1 99  CIN (J/kg)
acc_conv   1 99  accumulated subgrid-scale precip (mm)
acc_grid   1 99  accumulated grid-scale precip (mm)
tc2        1 99  temperature at 2 m AGL (deg C)
ept2       1 99  equiv. potential temp at 2 m AGL (K) 
rh2        1 99  relative hum at 2 m AGL 
u10        1 99  zonal wind at 2 m AGL (m/s)
v10        1 99  meridional wind at 2 m AGL (m/s)
ppw        1 99  precipitable water (mm)
hyd_path   1 99  total hydrometeor path (mm)
hyd_path_l 1 99 liquid hydrometeor path (mm)
cltopt     1 99  cloud top temperature (deg C)
pblh       1 99  PBL height (m) 
flux_lh    1 99  sfc latent heat flux (W/m**2)
flux_sh    1 99  sfc sensible heat flux (W/m**2)
u         26 99  zonal wind at pressure levels (m/s)
v         26 99  meridional wind at pressure levels (m/s)
w         26 99  vertical velocity (m/s)
z         26 99  geopotential height (m)
tempc     26 99  temperature (deg C)
rh        26 99  relative humidity (%)
qhydro    26 99  total hydrometeor mixing ratio (kg/kg)
qhydro_l  26 99  liquid hydrometeor mixing ratio (kg/kg)
pv        26 99  potential vorticity (PVU)
relvort   26 99  relative vorticity (10^-5 s^-1)
absvort   26 99  absolute vorticity (10^-5 s^-1)
endvars
EOF

cat <<EOF>zetap_level_${time_index}.ctl
dset zetap_level_${time_index}.gdat
options little_endian
undef -999999
title  OUTPUT FROM WRF V3.0 MODEL
pdef  513 486 lcc  29.85936 117.885 1.  1.  30.00000  60.00000 127.0 2700. 2700.
xdef 513 linear   116.9581  0.0309388
ydef 486 linear    30.00962 0.02463833
zdef  5 levels 1 2 3 4 5
tdef 1 linear 00z04apr2003 1hr
vars 9
z          5 99 height (m)
press      5 99 pressure (hPa)
tempc      5 99 temperature (deg C)
rh         5 99 RH (%)
u          5 99 zonal wind (m/s)
v          5 99 meridional wind (m/s)
w          5 99 vertical velocity (m/s)
qh_tot     5 99 total hydrometeor mr (kg/kg)
qh_l       5 99 liquid hydrometeor mr (kg/kg)
endvars
EOF

set NBARB = 15

else
 echo "DOMAIN entered not valid"
 echo "Must exit"
 exit

endif

#mv all.dat press_level_${time_index}.gdat

# ==== 250-hPa ===
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_hgt250.exec $NBARB '
'draw title $MEMBER : $DOMAIN 250-hPa height (m): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif hgt250_${time_index}_${DOMAIN}.gif



# ======= plot 500-hPa abs vort ======
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_absvort500.exec $NBARB '
'draw title $MEMBER : $DOMAIN 500-hPa abs vort (10-5 s-1): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif absvort500_${time_index}_${DOMAIN}.gif
# ======
# RH

cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_rh.exec 975 $NBARB '
'draw title $MEMBER : $DOMAIN 975-hPa RH: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif rh975_${time_index}_${DOMAIN}.gif
   # ===============

cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_rh.exec 950 $NBARB '
'draw title $MEMBER : $DOMAIN 975-hPa RH: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif rh975_${time_index}_${DOMAIN}.gif
   # ===============


cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_rh.exec 900 $NBARB '
'draw title $MEMBER : $DOMAIN 900-hPa RH: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif rh900_${time_index}_${DOMAIN}.gif

   # ===============
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_rh.exec 850 $NBARB '
'draw title $MEMBER : $DOMAIN 850-hPa RH: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif rh850_${time_index}_${DOMAIN}.gif

   # ===============


cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_rh.exec 700 $NBARB '
'draw title $MEMBER : $DOMAIN 700-hPa RH: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif rh700_${time_index}_${DOMAIN}.gif

# ======== plot 850-hPa height =====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_hgt850.exec $NBARB '
'draw title $MEMBER : $DOMAIN 850-hPa height (m): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif hgt850_${time_index}_${DOMAIN}.gif

# ===== plot 850-hPa tempc ====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_tmp850.exec $NBARB '
'draw title $MEMBER : $DOMAIN 850-hPa tempc (C): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif tmp850_${time_index}_${DOMAIN}.gif

# ======== plot 1000-hPa height =====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_hgt1000.exec $NBARB '
'draw title $MEMBER : $DOMAIN 1000-hPa height (m): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif hgt1000_${time_index}_${DOMAIN}.gif


# ===== plot 2-m tempc ====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_tmp2m.exec $NBARB '
'draw title $MEMBER : $DOMAIN 2-m tempc (C): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif tmp2m_${time_index}_${DOMAIN}.gif

# ===== plot SLP ====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
*'exec plot_slp.exec $NBARB '
'exec plot_slp_magen.exec $NBARB '
'draw title $MEMBER : $DOMAIN SLP (hPa): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif slp_${time_index}_${DOMAIN}.gif


# ======= plot cloud top temperature =======
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_cloud_top_tempc.exec '
'draw title $MEMBER : $DOMAIN cloud top tempc (deg C): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif cloud_top_tempc_${time_index}_${DOMAIN}.gif

# ====== plot CAPE =======
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_cape.exec '
'draw title $MEMBER : $DOMAIN CAPE (J/kg): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif cape_${time_index}_${DOMAIN}.gif

# ===== plot CIN ==========
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_cin.exec '
'draw title $MEMBER : $DOMAIN CIN (J/kg): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif cin_${time_index}_${DOMAIN}.gif

# ====== plot PBL height =====
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_pblh.exec '
'draw title $MEMBER : $DOMAIN PBL height (m): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif pblh_${time_index}_${DOMAIN}.gif

# ====== plot total condensate path ======
cat <<EOF>comm.scr
'open press_level_${time_index}.ctl'
'exec plot_total_cond_wp.exec '
'draw title $MEMBER : $DOMAIN total cond path (mm): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif tot_cond_path_${time_index}_${DOMAIN}.gif

# ====== plot zetap 
cat <<EOF>comm.scr
'open zetap_level_${time_index}.ctl'
'exec plot_qh_zetap.exec 3 '
'draw title $MEMBER : $DOMAIN qh_tot (g/kg) zetap 3: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif qh_tot_sigz0003_${time_index}_${DOMAIN}.gif

# ====== plot zetap
cat <<EOF>comm.scr
'open zetap_level_${time_index}.ctl'
'exec plot_qh_zetap.exec 5 '
'draw title $MEMBER : $DOMAIN qh_tot (g/kg) zetap 5: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif qh_tot_sigz0005_${time_index}_${DOMAIN}.gif


 endif

end

#rm -f *.gdat
