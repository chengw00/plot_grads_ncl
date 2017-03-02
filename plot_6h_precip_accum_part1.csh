#!/bin/csh

set GRADS_PATH = /home/chengw/bin/grads

setenv GADDIR $GRADS_PATH/lib
setenv GASCRP $GRADS_PATH/lib/scripts

set filelist_wrf_d03 = ` ls -1 wrfout_d03* `

foreach ff3 ( $filelist_wrf_d03 )
 set year  = ` echo $ff3 | cut -c12-15 `
 set month = ` echo $ff3 | cut -c17-18 `
 set day   = ` echo $ff3 | cut -c20-21 `
 set hour  = ` echo $ff3 | cut -c23-24 `
 if ( $hour == 00 || $hour == 06 || $hour == 12 || $hour == 18 ) then
  set time_index = ` echo $ff3 | cut -c12-30 ` 
  echo $ff3 $year $month $day $hour $time_index
  set file_prefix = ` echo $ff3 | cut -c1-30 `

cat <<EOF>namelist.ARWpost
&datetime
 start_date = '${time_index}',
 end_date   = '${time_index}',
 interval_seconds = 3600,
 tacc = 0,
 debug_level = 0,
/

&io
 io_form_input  = 2,
 input_root_name = './${file_prefix}.nc'
 output_root_name = './test'
 plot = 'basic'
 fields = 'rainc'
 output_type = 'grads'
 mercator_defs = .true.
/
 split_output = .true.
 frames_per_outfile = 1

 output_type = 'grads'

 plot = 'list' 
! Below is a list of all available diagnostics
 fields = 'rainnc'
 

&interp
 interp_method = 1,
 interp_levels = 300.,
/
EOF

rm -f test.dat
ARWpost.exe
mv test.dat precip_6h_$year$month$day${hour}.gdat

 endif
end

cat <<EOF>comm.scr
'open 42h.ctl'
'open 48h.ctl'
'exec plot_6h_precip_accum.exec'
'enable print test'
'print'
'disable print'
'quit'
EOF

#$GRADS_PATH/bin/grads -lbc comm.scr

#$GRADS_PATH/bin/gxgif -vr test
