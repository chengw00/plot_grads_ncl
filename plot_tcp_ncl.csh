#!/bin/csh

set DOMAIN = $1

set TIME_INDEX = $2

set NREC = $3

set NCARG_ROOT = $4

cat <<EOF>test.ncl
;*************************************************
; raster_t.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  setfileoption("bin","ReadByteOrder","LittleEndian")
  ff = addfile("wrfout_${DOMAIN}_${TIME_INDEX}:00:00.nc","r")

  xlat      = ff->XLAT
  xlon      = ff->XLONG

  tlat1     = ff@TRUELAT1
  tlat2     = ff@TRUELAT2
  clat      = ff@MOAD_CEN_LAT
  clon      = ff@STAND_LON

  ; ===============
  dsizes_xlat = dimsizes(xlat)
  n_y = dsizes_xlat(1)
  n_x = dsizes_xlat(2)

  ; ===== read 2d file =======
  template = new((/n_y,n_x/),float,-999999)
  qh_total_path = template

  dims_2d = (/n_y,n_x/)
  nrec = $NREC - 1
  file_2d = "store_2d_${TIME_INDEX}_${DOMAIN}.gdat"
  qh_total_path = fbindirread(file_2d, nrec, dims_2d, "float")

  qh_total_path@lat2d = xlat(0,:,:)
  qh_total_path@lon2d = xlon(0,:,:)
  qh_total_path@units = " "

  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  ;gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.

  ; ==== WC: define my own colors ========
  cmap = RGBtoCmap("color_ybl_precip.txt")
  gsn_define_colormap(wks,cmap)
  ; ============================

  res                 = True                ; Plot mods desired.

  ; ====== WC: define my own contour levels =========
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels  = (/0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45/)

  ; ====== WC : change map thickness
  res@mpGeophysicalLineThicknessF = 3
  res@mpNationalLineThicknessF = 3
  ; ====== WC: define title ====
  res@tiMainFontHeightF = 0.02
  res@tiMainString = "${TIME_INDEX}"
  ; ============================================
  res@gsnMaximize     = True                ; Maximize plot in frame.

  res@gsnStringFontHeightF         = 0.013
  res@gsnRightStringOrthogonalPosF = 0.02
  res@gsnLeftStringOrthogonalPosF  = 0.02

  res@cnFillOn        = True               ; Color plot desired.
  res@cnLinesOn       = False              ; Turn off contour lines      .
  res@cnLineLabelsOn  = False              ; Turn off contour labels.

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = tlat1               ; two parallels
  res@mpLambertParallel2F    = tlat2
  res@mpLambertMeridianF     = clon             ; central meridian
  res@mpLimitMode            = "LatLon"

  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = xlat(0,0,0)
  res@mpLeftCornerLonF  = xlon(0,0,0)
  res@mpRightCornerLatF = xlat(0,n_y-1,n_x-1)
  res@mpRightCornerLonF = xlon(0,n_y-1,n_x-1)

  res@gsnSpreadColors     = False            ; Use full range of colormap.
  ;res@gsnSpreadColorStart = 16              ; Start color in colormap.
  ;res@gsnSpreadColorEnd   = 190             ; End color in colormap.
  res@cnInfoLabelOn       = False           ; Turn off info label.

  res@mpGeophysicalLineColor = "black"     ; color of continental outlines
  res@mpUSStateLineColor     = "black" 
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed
 
  ;res@mpOutlineBoundarySets  = "GeophysicalAndUSStates" 
 
  ; ====== modified by WC: 2010-01-27 =====
  res@mpOutlineBoundarySets  = "National"
  ; ======================================

  ; ====
  res@mpDataBaseVersion = "MediumRes"
  ;res@mpDataBaseVersion = "HighRes"
  res@mpDataSetName = "Earth..4" 
  ; ======
  res@pmTickMarkDisplayMode = "Always"       ; Turn on map tickmarks.

  res@pmLabelBarWidthF     = 0.6
  res@lbLabelFontHeightF   = 0.013
  res@lbLabelStride        = 1
  res@lbBoxLinesOn         = False

  map = gsn_csm_contour_map(wks,qh_total_path,res) 

end
EOF

ncl test.ncl

convert -geometry 1200x1200 -density 300 -trim raster.pdf raster.jpg
mv raster.jpg tot_cond_path_${TIME_INDEX}_${DOMAIN}.jpg
