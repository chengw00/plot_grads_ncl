;*************************************************
; raster_3.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.

  ; ==== WC: define my own colors ========
  cmap = RGBtoCmap("color_yellow_green_rh.txt")
  gsn_define_colormap(wks,cmap)
  ; ============================


  gsn_merge_colormaps(wks,cmap,"BlAqGrYeOrRe")
  gsn_draw_colormap(wks)

end
