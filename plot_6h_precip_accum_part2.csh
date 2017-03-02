#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

#set filelist = ` ls -1 precip_6h_2009041600.gdat `
set filelist = ` ls -1 precip_6h_*.gdat `


set ncount = 0
STEP1:
@ ncount = $ncount + 1
if ( $ncount > $#filelist ) goto STEP1b

set time_index  = ` echo $filelist[$ncount] | cut -c11-20 `
set file_prefix = ` echo $filelist[$ncount] | cut -c1-20 `
 
echo $filelist[$ncount] $time_index $file_prefix

cat <<EOF> ${file_prefix}.ctl
dset $filelist[$ncount]
options  byteswapped
undef 1.e30
title  OUTPUT FROM WRF V3.0 MODEL
pdef  135 255 lcc  28.24144 32.33243 1.  1.  30.00000  60.00000 35.00000   3333.300   3333.300
xdef 135 linear   32.17249  0.03579462
ydef 255 linear   28.25152  0.03020343
zdef    1 levels 300
tdef    1 linear 00Z18APR2009      60MN      
VARS   27
LU_INDEX       1  0  LAND USE CATEGORY (-)
U              1  0  x-wind component (m s-1)
V              1  0  y-wind component (m s-1)
W              1  0  z-wind component (m s-1)
PH             1  0  perturbation geopotential (m2 s-2)
PHB            1  0  base-state geopotential (m2 s-2)
T              1  0  perturbation potential temperature (theta-t0) (K)
MU             1  0  perturbation dry air mass in column (Pa)
MUB            1  0  base state dry air mass in column (Pa)
P              1  0  perturbation pressure (Pa)
PB             1  0  BASE STATE PRESSURE (Pa)
T2             1  0  TEMP at 2 M (K)
TH2            1  0  POT TEMP at 2 M (K)
PSFC           1  0  SFC PRESSURE (Pa)
U10            1  0  U at 10 M (m s-1)
V10            1  0  V at 10 M (m s-1)
LANDMASK       1  0  LAND MASK (1 FOR LAND, 0 FOR WATER) (-)
VEGFRA         1  0  VEGETATION FRACTION (-)
SST            1  0  SEA SURFACE TEMPERATURE (K)
HGT            1  0  Terrain Height (m)
TSK            1  0  SURFACE SKIN TEMPERATURE (K)
RAINC          1  0  ACCUMULATED TOTAL CUMULUS PRECIPITATION (mm)
RAINNC         1  0  ACCUMULATED TOTAL GRID SCALE PRECIPITATION (mm)
XLAT           1  0  LATITUDE, SOUTH IS NEGATIVE (degree_north)
XLONG          1  0  LONGITUDE, WEST IS NEGATIVE (degree_east)
TMN            1  0  SOIL TEMPERATURE AT LOWER BOUNDARY (K)
PBLH           1  0  PBL HEIGHT (m)
ENDVARS
@ global String comment TITLE =  OUTPUT FROM WRF V3.0 MODEL
@ global String comment SIMULATION_START_DATE = 2009-04-16_00:00:00
@ global String comment GRIDTYPE = C
@ global String comment MMINLU = USGS
@ global String comment WEST-EAST_GRID_DIMENSION =   541
@ global String comment SOUTH-NORTH_GRID_DIMENSION =   571
@ global String comment BOTTOM-TOP_GRID_DIMENSION =    37
@ global String comment MAP_PROJ =     1
@ global String comment DX =      3333.30
@ global String comment DY =      3333.30
@ global String comment CEN_LAT =        40.60
@ global String comment CEN_LON =      -101.86
@ global String comment TRUELAT1 =        30.00
@ global String comment TRUELAT2 =        60.00
@ global String comment MOAD_CEN_LAT =        39.00
@ global String comment STAND_LON =      -104.00
@ global String comment DIFF_OPT =     1
@ global String comment KM_OPT =     4
@ global String comment DAMP_OPT =     3
@ global String comment KHDIF =         0.00
@ global String comment KVDIF =         0.00
@ global String comment MP_PHYSICS =     2
@ global String comment RA_LW_PHYSICS =     1
@ global String comment RA_SW_PHYSICS =     1
@ global String comment SF_SFCLAY_PHYSICS =     1
@ global String comment SF_SURFACE_PHYSICS =     2
@ global String comment BL_PBL_PHYSICS =     1
@ global String comment CU_PHYSICS =     0
@ global String comment SURFACE_INPUT_SOURCE =     1
@ global String comment SST_UPDATE =     0
@ global String comment GRID_FDDA =     0
@ global String comment UCMCALL =     0
@ global String comment FEEDBACK =     0
@ global String comment SMOOTH_OPTION =     0
@ global String comment W_DAMPING =     0
@ global String comment PD_MOIST =     1
@ global String comment PD_SCALAR =     1
@ global String comment PD_TKE =     0
@ global String comment OBS_NUDGE_OPT =     1
@ global String comment GRID_ID =     3
@ global String comment PARENT_ID =     2
@ global String comment I_PARENT_START =    38
@ global String comment J_PARENT_START =    21
@ global String comment PARENT_GRID_RATIO =     3
@ global String comment DT =        16.67
@ global String comment ISWATER =    16
@ global String comment ISICE =    24
@ global String comment ISURBAN =     1
@ global String comment ISOILWATER =    14
EOF

goto STEP1

# ============
STEP1b:

set ncount = 1
STEP2:
@ ncount = $ncount + 1
@ ncount_m1 = $ncount - 1
if ( $ncount > $#filelist ) exit

set file_prefix1 = ` echo $filelist[$ncount_m1] | cut -c1-20 `
set time_index  = ` echo $filelist[$ncount] | cut -c11-20 `
set file_prefix2 = ` echo $filelist[$ncount] | cut -c1-20 `
 
echo $filelist[$ncount_m1] $filelist[$ncount] $time_index $file_prefix

cat <<EOF>comm.scr
'open ${file_prefix1}.ctl '
'open ${file_prefix2}.ctl '
'exec plot_6h_precip_accum.exec'
'draw title 6-h accum precip (mm): $time_index '
'enable print test'
'print'
'disable print'
'quit'
EOF

$GRADS_PATH/bin/grads -lbc comm.scr

$GRADS_PATH/bin/gxgif -vr -x 1100 -y 950 test
rm -f test
mv test.gif 6h_precip_accum_${time_index}.gif


goto STEP2

