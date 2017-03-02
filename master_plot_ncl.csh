#!/bin/csh

set DOMAIN = $1

set TIME_INDEX = $2

set NCARG_ROOT = $3

# ====== plot total condensate path ======

if ( -e store_2d_${TIME_INDEX}_${DOMAIN}.gdat ) then
 nrec = 17
 plot_tcp_ncl.csh $DOMAIN $TIME_INDEX $nrec $NCARG_ROOT
endif

