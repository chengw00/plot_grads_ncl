#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads 

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

set filelist = ` ls -1 eta_level_*.gdat `
# eta_level_2009041703.gdat

set ncount = 0

foreach ff ( $filelist )

@ ncount = $ncount + 1
if ( $ncount > $#filelist ) exit

set time_index  = ` echo $filelist[$ncount] | cut -c11-20 `
set file_prefix = ` echo $filelist[$ncount] | cut -c1-20 `
 
echo $filelist[$ncount] $time_index $file_prefix

cat <<EOF> ${file_prefix}.ctl
dset $filelist[$ncount]
options little_endian
undef 1.e30
title  OUTPUT FROM WRF V3.0 MODEL
pdef  135 255 lcc  28.24144 32.33243 1.  1.  30.00000  60.00000 35.00000   3333.300   3333.300
xdef 135 linear   32.17249  0.03579462
ydef 255 linear   28.25152  0.03020343
zdef    1 levels 1
tdef    1 linear 00Z18APR2009      60MN      
VARS  11
XLAT           1  0  LATITUDE, SOUTH IS NEGATIVE (degree_north)
XLONG          1  0  LONGITUDE, WEST IS NEGATIVE (degree_east)
HGT            1  0  Terrain Height (m)
T2             1  0  TEMP at 2 M (deg C)
U10            1  0  U at 10 M (m s-1)
V10            1  0  V at 10 M (m s-1)
SLP            1  0  SLP
PBLH           1  0  PBL height (m)
T              1  0  temperature (deg C)
U              1  0  x-wind component (m s-1)
V              1  0  y-wind component (m s-1)
ENDVARS
EOF

rm -f test

cat <<EOF>comm.scr
'open ${file_prefix}.ctl '
'exec plot_tempc_eta1.exec'
'draw title T (deg C) at lowest eta level: $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif tempc_eta1_${time_index}.gif
# =================
cat <<EOF>comm.scr
'open ${file_prefix}.ctl '
'exec plot_pblh.exec'
'draw title PBLH (m): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif pblh_${time_index}.gif

end

