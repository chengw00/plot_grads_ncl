#!/bin/csh

set DOMAIN = $1

set TIME_INDEX = $2

set NREC = $3

set NCARG_ROOT = $4

set LEV = $5

cat <<EOF>test.ncl
;*************************************************
; raster_3.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  ff = addfile("wrfout_${DOMAIN}_${TIME_INDEX}:00:00.nc","r")

  xlat      = ff->XLAT
  xlon      = ff->XLONG
  terr      = ff->HGT

  terr@lat2d = xlat(0,:,:)
  terr@lon2d = xlon(0,:,:)
  terr@units = "m"
  terr@description = " "

  printVarSummary(terr)

  tlat1     = ff@TRUELAT1
  tlat2     = ff@TRUELAT2
  clat      = ff@MOAD_CEN_LAT
  clon      = ff@STAND_LON

  ; ===============
  dsizes_xlat = dimsizes(xlat)
  n_y = dsizes_xlat(1)
  n_x = dsizes_xlat(2)

  ;do j=1,n_y
  ; do i=1,n_x
  ;  ii=i-1
  ;  jj=j-1
  ;  veg(jj,ii)=i
  ; end do
  ;end do

 ; =========
 rh_zeta = wrf_user_getvar(ff,"rh",0)    ; slp (Pa)

 printVarSummary(rh_zeta)

 ; ===== read 3d pressure file =======
  n_p = 1
 
  ; ==== fake level, just one, doesn't matter 
  template_1d = new(n_p,float,-999999)
  pres_1d = (/1000 /)

  print(pres_1d)

  kplot = 0

  template = new((/n_p,n_y,n_x/),float,-999999)
  a_3d = template

  dims_3d = (/n_p,n_y,n_x/)
  nrec = 13 - 1     ; RH
  file_3d = "store_2d_${TIME_INDEX}_d03.gdat"
  a_3d = fbindirread(file_3d, nrec, dims_3d, "float")

  a_3d@lat2d = xlat(0,:,:)
  a_3d@lon2d = xlon(0,:,:)
  a_3d@units = ""

  print(a_3d(:,50,50))

  u_3d = template
  nrec = 14 - 1
  u_3d = fbindirread(file_3d, nrec, dims_3d, "float")
  u_3d@lat2d = xlat(0,:,:)
  u_3d@lon2d = xlon(0,:,:)
  u_3d@units = "m/s"

  v_3d = template
  nrec = 15 - 1
  v_3d = fbindirread(file_3d, nrec, dims_3d, "float")
  v_3d@lat2d = xlat(0,:,:)
  v_3d@lon2d = xlon(0,:,:)
  v_3d@units = "m/s"

  ;kplot = 6 - 1
  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.

  ; ==== WC: define my own colors ========
  ;cmap = RGBtoCmap("color_yellow_green_rh.txt")
  ;gsn_define_colormap(wks,cmap)
  ; ============================

  ; ======= merge color map =======
  ;gsn_merge_colormaps(wks,cmap,"BlAqGrYeOrRe")
  ; ==============================

  res                 = True                ; Plot mods desired.
  res@cnFillOn         = True

  ;res@gsnMaximize     = True                ; Maximize plot in frame.

  res@gsnStringFontHeightF         = 0.013
  res@gsnRightStringOrthogonalPosF = 0.02
  res@gsnLeftStringOrthogonalPosF  = 0.02

  res@cnLinesOn       = False              ; Turn off contour lines      .
  res@cnLineLabelsOn  = False              ; Turn off contour labels.

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = tlat1               ; two parallels
  res@mpLambertParallel2F    = tlat2
  res@mpLambertMeridianF     = clon             ; central meridian

  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = xlat(0,0,0)
  res@mpLeftCornerLonF  = xlon(0,0,0)
  res@mpRightCornerLatF = xlat(0,n_y-1,n_x-1)
  res@mpRightCornerLonF = xlon(0,n_y-1,n_x-1)

  ; =======used for Native LC projection ========
  ; not needed since lat/lon defined with field
  ;res@tfDoNDCOverlay         = True

  ; ======= box 1 ==========
  ;res@mpLeftCornerLatF  = 34
  ;res@mpLeftCornerLonF  = 32
  ;res@mpRightCornerLatF = 36
  ;res@mpRightCornerLonF = 35

  ; ======= box 2 =========
  ;res@mpLeftCornerLatF  = 32
  ;res@mpLeftCornerLonF  = 34.5
  ;res@mpRightCornerLatF = 34.5
  ;res@mpRightCornerLonF = 36.5

  ; ====== box 3 ========
  ;res@mpLeftCornerLatF  = 30.2
  ;res@mpLeftCornerLonF  = 33.5
  ;res@mpRightCornerLatF = 32.2
  ;res@mpRightCornerLonF = 35.5

  ; ====== box 4 ========
  ;res@mpLeftCornerLatF  = 30.0
  ;res@mpLeftCornerLonF  = 31.8
  ;res@mpRightCornerLatF = 32.2
  ;res@mpRightCornerLonF = 33.5

  res@gsnSpreadColors     = True            ; Use full range of colormap.
  ;res@gsnSpreadColorStart = 16              ; Start color in colormap.
  ;res@gsnSpreadColorEnd   = 190             ; End color in colormap.
  res@cnInfoLabelOn       = False           ; Turn off info label.

  res@mpGeophysicalLineColor = "blue"     ; color of continental outlines
  res@mpUSStateLineColor     = "blue" 
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed
  res@mpOutlineBoundarySets  = "GeophysicalAndUSStates" 
 
  ; ====== modified by WC: 2010-01-27 =====
  ;res@mpOutlineBoundarySets  = "National"
  ; ======================================

  ; ====
  res@mpDataBaseVersion = "MediumRes"
  ;res@mpDataBaseVersion = "HighRes"
  res@mpDataSetName = "Earth..4" 
  ; ======
  res@pmTickMarkDisplayMode = "Always"       ; Turn on map tickmarks.

  ;res@cnFillMode           = "RasterFill"
  ;res@cnFillMode           = "CellFill"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -200              ; min contour level
  res@cnMaxLevelValF       = 1500              ; max contour level
  res@cnLevelSpacingF      =  100            ; contour spacing

  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels  = (/70, 90/)

  res@pmLabelBarWidthF     = 0.6
  res@lbLabelFontHeightF   = 0.013
  res@lbLabelStride        = 3
  res@lbBoxLinesOn         = False

  res@gsnFrame     =  False                   ; do not advance the frame
  ; ================
  ; wind vectors
  vecres                  =  True             ; vector only resources
  vecres@gsnDraw          = False           ; don't draw
  vecres@gsnFrame         = False           ; don't advance frame
  ;vecres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
  vecres@vcGlyphStyle     = "WindBarb"
  vecres@vcWindBarbLineThicknessF= 3.5              ; set the wind barb thickness
  vecres@vcRefMagnitudeF  = 20              ; define vector ref mag
  vecres@vcMinDistanceF   = 0.045
  vecres@vcRefLengthF     = 0.045           ; define length of vec ref
  vecres@gsnRightString   = " "             ; turn off right string
  vecres@gsnLeftString    = " "             ; turn off left string
  vecres@tiXAxisString    = " "             ; turn off axis label 


  vecres@vcMonoWindBarbColor   = True
  vecres@vcWindBarbColor = "purple"
  
  vecres@vcRefAnnoOn = False                ; reference vector box off

  ;vecres@vcRefAnnoOrthogonalPosF = -.535    ; move ref vector into plot

  ;vecres@gsnMaximize     = True

  vecres@mpProjection           = "LambertConformal"
  vecres@mpLambertParallel1F    = tlat1               ; two parallels
  vecres@mpLambertParallel2F    = tlat2
  vecres@mpLambertMeridianF     = clon             ; central meridian

  vecres@mpLimitMode       = "Corners"            ; choose range of map
  vecres@mpLeftCornerLatF  = xlat(0,0,0)
  vecres@mpLeftCornerLonF  = xlon(0,0,0)
  vecres@mpRightCornerLatF = xlat(0,n_y-1,n_x-1)
  vecres@mpRightCornerLonF = xlon(0,n_y-1,n_x-1)

  ; ====== WC: define title ====
  res@tiMainFontHeightF = 0.02
  res@tiMainString = "${TIME_INDEX}"
  ; ============================================

  vecres@mpGeophysicalLineColor = "white"     ; color of continental outli
  vecres@mpNationalLineColor    = "white"
  vecres@mpUSStateLineColor     = "white"
  vecres@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed

  ; ====== modified by WC: 2010-01-27 =====
  ;vecres@vcRefAnnoOrthogonalPosF = -0.06          ; move ref vector
  ;vecres@vcRefAnnoParallelPosF   = 1.            ; move ref vector

  vecres@mpFillOn               = False   ; no fill

  ; =======used for Native LC projection ========
  ; not needed since lat/lon defined with field
  ;vecres@tfDoNDCOverlay         = True
  ; =======================================
  vecres@mpOutlineBoundarySets  = "National"

  vecres@mpGeophysicalLineThicknessF = 3
  vecres@mpNationalLineThicknessF = 3
  ; ======================================

  ; ====
  vecres@mpDataBaseVersion = "MediumRes"
  ;vecres@mpDataBaseVersion = "HighRes"
  vecres@mpDataSetName = "Earth..4"

  ; =================

  conres                 = True                ; Plot mods desired.
  conres@cnFillOn         = False
  ;conres@gsnMaximize     = True                ; Maximize plot in frame.

  conres@cnLineThicknessF  = 3.0              ; contour line thickness
  conres@cnLineColor = "black"                ; contour line color
  conres@cnLineLabelFontColor = "black"       ; label font color
 
  conres@gsnStringFontHeightF         = 0.013
  conres@gsnRightStringOrthogonalPosF = 0.02
  conres@gsnLeftStringOrthogonalPosF  = 0.02

  conres@cnLinesOn       = True              ; Turn off contour lines      .
  conres@cnLineLabelsOn  = False              ; Turn off contour labels.

  conres@mpProjection           = "LambertConformal"
  conres@mpLambertParallel1F    = tlat1               ; two parallels
  conres@mpLambertParallel2F    = tlat2
  conres@mpLambertMeridianF     = clon             ; central meridian

  conres@mpLimitMode       = "Corners"            ; choose range of map
  conres@mpLeftCornerLatF  = xlat(0,0,0)
  conres@mpLeftCornerLonF  = xlon(0,0,0)
  conres@mpRightCornerLatF = xlat(0,n_y-1,n_x-1)
  conres@mpRightCornerLonF = xlon(0,n_y-1,n_x-1)

  ; ====== WC: define title ====
  ;conres@tiMainFontHeightF = 0.02
  ;conres@tiMainString = "2009-05-31_17"
  ; ============================================

  ; =======used for Native LC projection ========
  ; not needed since lat/lon defined with field
  ;conres@tfDoNDCOverlay         = True

  ; ======= box 1 ==========
  ;conres@mpLeftCornerLatF  = 34
  ;conres@mpLeftCornerLonF  = 32
  ;conres@mpRightCornerLatF = 36
  ;conres@mpRightCornerLonF = 35

  ; ======= box 2 =========
  ;conres@mpLeftCornerLatF  = 32
  ;conres@mpLeftCornerLonF  = 34.5
  ;conres@mpRightCornerLatF = 34.5
  ;conres@mpRightCornerLonF = 36.5

  ; ====== box 3 ========
  ;conres@mpLeftCornerLatF  = 30.2
  ;conres@mpLeftCornerLonF  = 33.5
  ;conres@mpRightCornerLatF = 32.2
  ;conres@mpRightCornerLonF = 35.5

  ; ====== box 4 ========
  ;conres@mpLeftCornerLatF  = 30.0
  ;conres@mpLeftCornerLonF  = 31.8
  ;conres@mpRightCornerLatF = 32.2
  ;conres@mpRightCornerLonF = 33.5

  ;conres@gsnSpreadColors     = False            ; Use full range of colormap.
  ;conres@gsnSpreadColorStart = 16              ; Start color in colormap.
  ;conres@gsnSpreadColorEnd   = 190             ; End color in colormap.
  conres@cnInfoLabelOn       = False           ; Turn off info label.

  conres@mpGeophysicalLineColor = "green"     ; color of continental outlines
  conres@mpNationalLineColor    = "green"
  conres@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed

  conres@mpGeophysicalLineThicknessF = 3.
  conres@mpNationalLineThicknessF = 3.

  ;conres@mpOutlineBoundarySets  = "GeophysicalAndUSStates" 
 
  ; ====== modified by WC: 2010-01-27 =====
  conres@mpOutlineBoundarySets  = "National"
  ; ======================================

  ; ====
  conres@mpDataBaseVersion = "MediumRes"
  ;conres@mpDataBaseVersion = "HighRes"
  conres@mpDataSetName = "Earth..4" 
  ; ======
  conres@pmTickMarkDisplayMode = "Always"       ; Turn on map tickmarks.

  ;conres@cnFillMode           = "RasterFill"
  ;conres@cnFillMode           = "CellFill"
  ;conres@cnLevelSelectionMode = "ManualLevels"
  ;conres@cnMinLevelValF       = a_min              ; min contour level
  ;conres@cnMaxLevelValF       = a_max              ; max contour level
  ;conres@cnLevelSpacingF      = a_int             ; contour spacing

  conres@cnLevelSelectionMode = "ExplicitLevels"
  conres@cnLevels  = (/90 /)

  conres@pmLabelBarWidthF     = 0.6
  conres@lbLabelFontHeightF   = 0.013
  conres@lbLabelStride        = 2
  conres@lbBoxLinesOn         = False

  conres@gsnFrame     =  False                   ; do not advance the frame
; =====================
  ; ===================
  ;map  = gsn_csm_contour_map(wks,terr(0,:,:),res)
  ;map2 = gsn_csm_contour_map(wks,a_3d(kplot,:,:),conres)
  plot = gsn_csm_contour_map_overlay(wks,terr(0,:,:),a_3d(kplot,:,:), \
                                     res,conres)

  plot_vec =  gsn_csm_vector_map(wks,u_3d(kplot,:,:),v_3d(kplot,:,:),vecres)
  ;plot_con =  gsn_csm_contour_map(wks,terr(0,:,:),conres)

  draw(plot_vec)

  frame(wks) 
end
EOF

ncl test.ncl

convert -geometry 1200x1200 -density 300 -trim raster.pdf raster.jpg
mv raster.jpg terrain_rh_${LEV}_${TIME_INDEX}_${DOMAIN}.jpg
