#!/bin/csh

set GMID   = $1
set MEMBER = $2
set filelist = $3  # cycle list

#set DOMAIN = ( d01 d02 d03 )
set DOMAIN = ( d01 d02  )

set NCARG_ROOT = /opt/ncl
set NCARG_LIB =  /opt/ncl/lib

set NCARG_ROOT_precip = /opt/ncl
set NCARG_LIB_precip =  /opt/ncl/lib

#set MEMBER = GFS_WCTRL
set HOME = /data1/chengw/cepri/reanl/plot_v2_part1/scripts_ptop_10_ncl_v4

set PATH = /data1/chengw/cepri/reanl/plot_v2_part1/$MEMBER

cd $PATH
#if ( $filelist == "" ) then
# set filelist = ( 2009111618 2009111700 )
#endif
#set filelist = ` /bin/ls -1 $PATH `
#set filelist = ` cat /is-hydro/chengw/kma/floodcase/plot_v2/cycle.dat `
set filelist = ( 2012012700 )

foreach ff ( $filelist )

 set first_digit = ` echo $ff | cut -c1-1 `
 if ( $first_digit == 2 ) then
  cd $HOME
  echo $ff
  tar cf $PATH/$ff/abc.tar *
  cd $PATH/$ff
  tar xvf abc.tar
  rm -f abc.tar *dat *.gif test

# plot_6h_precip_accum_part1.csh
# plot_6h_precip_accum_part2.csh
# plot_eta_levels_part1.csh
# plot_eta_levels_part2.csh

# ==========================================
# UNCOMMENT THIS PART WHEN RUNNING

  foreach dd ( $DOMAIN )

   ln -sf diag_wrf_v2_${dd}.exe diag_wrf_v2.exe
   plot_pressure_levels_part1.csh $dd $MEMBER $NCARG_ROOT
   plot_pressure_levels_part2.csh $dd $MEMBER
# =============================================
#   extract_precip_accum_3h.csh $NCARG_ROOT_precip $NCARG_LIB_precip $dd
#   plot_precip_accum_3h.csh $dd

   extract_precip_accum_1h.csh $NCARG_ROOT_precip $NCARG_LIB_precip $dd
   plot_precip_accum_1h.csh $dd

  end # end dd loop
# plot_sounding.csh
 endif

end


#rm -f *.gdat
#rm -f */*.gdat
