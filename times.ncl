;********************************************************
; WRF: latitude-z cross section.
;********************************************************
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
;load "/home/chengw/bin/ncl4.2.0/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "/home/chengw/bin/ncl4.2.0/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "/usr/local/ncarg-5.0.0/lib/ncarg/nclscripts/csm/contributed.ncl"   
;load "./roux/StationModel.ncl"
begin
  ; ==================
  setfileoption("bin","WriteByteOrder","LittleEndian")
  file_in = systemfunc ("/bin/ls -1 wrfout_d03_2009*.nc" )
   
  fall      = addfiles (file_in, "r")  
  ;========================

  Times = fall[:]->Times

  print(Times)
  end
