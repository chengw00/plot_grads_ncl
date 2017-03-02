#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads
set EXP = 2009060912_init_SmallDomain

set DOMAIN = $1   # d03

if ( $1 == "" ) then
 echo "MISSING INPUT"
 echo "EXIT"
 exit
endif

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

cd $DOMAIN
#tempc_plev0850_2008-07-19_06.gdat
set filelist = ` ls -1 pcpac*.gdat `
cd ..

foreach ff ( $filelist )
 set time_index = ` echo $ff | cut -c16-28 `
 echo $ff $time_index
end

# ================================
set ncount = 0
STEP1:
@ ncount = $ncount + 1
if ( $ncount > $#filelist ) goto STEP1b
set time_index  = ` echo $filelist[$ncount] | cut -c16-28 `
set file_prefix = ` echo $filelist[$ncount] | cut -c1-28 `


echo $filelist[$ncount] $time_index $file_prefix


if ( $DOMAIN == d01 ) then

cat <<EOF> ${file_prefix}.ctl
dset ./${DOMAIN}/${file_prefix}.gdat
options little_endian
undef -999999
pdef  279  194 lcc 9.226906 65.54395  1.  1.  10.0  40.0 100.0 27000. 27000.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  56 linear 1  1
tdef 1 linear 12z12dec2008 1hr
vars 1
accum    1 99  precip accum (mm)
endvars
EOF

endif

if ( $DOMAIN == d02 ) then

cat <<EOF> ${file_prefix}.ctl
dset ./${DOMAIN}/${file_prefix}.gdat
options little_endian
undef -999999
pdef  651  483 lcc 15.52342 73.7645  1.  1.  10.0  40.0 100.0 9000. 9000.
* --------- original ----------------------
xdef  651 linear   68.01399  0.1048811
ydef  483 linear   15.16394  0.07852066
zdef  56 linear 1  1
tdef 1 linear 12z12dec2008 1hr
vars 1
accum    1 99  precip accum (mm)
endvars
EOF

endif

if ( $DOMAIN == d03 ) then

cat <<EOF> ${file_prefix}.ctl
dset ./${DOMAIN}/${file_prefix}.gdat
options little_endian
undef -999999
pdef  513 486 lcc  29.85936 117.885 1.  1.  30.00000  60.00000 127.0 2700. 2700.
xdef 513 linear   116.9581  0.0309388
ydef 486 linear    30.00962 0.02463833
zdef  56 linear 1  1
tdef 1 linear 12z12dec2008 1hr
vars 1
accum    1 99  precip accum (mm)
endvars
EOF

endif

goto STEP1

# ============
STEP1b:

set ncount = 1
STEP2:
@ ncount = $ncount + 1
@ ncount_m1 = $ncount - 1
if ( $ncount > $#filelist ) exit

set file_prefix1 = ` echo $filelist[$ncount_m1] | cut -c1-28 `
set time_index  = ` echo $filelist[$ncount] | cut -c16-28 `
set file_prefix2 = ` echo $filelist[$ncount] | cut -c1-28 `

echo $filelist[$ncount_m1] $filelist[$ncount] $time_index $file_prefix

cat <<EOF>comm.scr
'open ${file_prefix1}.ctl '
'open ${file_prefix2}.ctl '
'exec plot_3h_precip_accum.exec'
'draw title 3-h accum precip (mm): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF


$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif pcp3h_sigz0000_${time_index}_${DOMAIN}.gif


goto STEP2


