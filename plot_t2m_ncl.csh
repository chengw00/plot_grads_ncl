#!/bin/csh

set DOMAIN = $1

set TIME_INDEX = $2

set NREC = $3

set NCARG_ROOT = $4

set LEV = $5

cat <<EOF>test.ncl
;*************************************************
; raster_t.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  ff = addfile("wrfout_${DOMAIN}_${TIME_INDEX}:00:00.nc","r")
  ;veg       = ff->LANDUSEF
  ;veg@lat2d = ff->XLAT_M
  ;veg@lon2d = ff->XLONG_M

  landuse   = ff->T2
  ;landuse   = ff->LANDUSEF
  xlat      = ff->XLAT
  xlon      = ff->XLONG

  veg       = landuse(0,:,:) - 273.15 ; convert to deg C

  u10       = ff->U10
  v10       = ff->V10

  ;veg       = veg + 0.001
  ;veg       = landuse(0,0,:,:)
  veg@lat2d = xlat(0,:,:)
  veg@lon2d = xlon(0,:,:)
 
  u10@lat2d = xlat(0,:,:)
  u10@lon2d = xlon(0,:,:)

  v10@lat2d = xlat(0,:,:)
  v10@lon2d = xlon(0,:,:)
 
  tlat1     = ff@TRUELAT1
  tlat2     = ff@TRUELAT2
  clat      = ff@MOAD_CEN_LAT
  clon      = ff@STAND_LON

  u10 = 2*u10  ; scale
  v10 = 2*v10
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
  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  ;gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.
  gsn_define_colormap(wks,"wh-bl-gr-ye-re")
 
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

  res@mpGeophysicalLineColor = "black"     ; color of continental outlines
  res@mpUSStateLineColor     = "black" 
  res@mpNationalLineColor    = "black"
  res@mpGeophysicalLineThicknessF = 1.5
  res@mpNationalLineThicknessF = 1.5

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
  res@cnMinLevelValF       = 0.0              ; min contour level
  res@cnMaxLevelValF       = 34.0              ; max contour level
  res@cnLevelSpacingF      = 2.             ; contour spacing

  res@pmLabelBarWidthF     = 0.6
  res@lbLabelFontHeightF   = 0.013
  res@lbLabelStride        = 2
  res@lbBoxLinesOn         = False

  res@gsnFrame     =  False                   ; do not advance the frame

  ; ====== WC: define title ====
  res@tiMainFontHeightF = 0.02
  res@tiMainString = "${TIME_INDEX}"
  ; ============================================

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

  ; ======= box 1 ==========
  ;vecres@mpLeftCornerLatF  = 34
  ;vecres@mpLeftCornerLonF  = 32
  ;vecres@mpRightCornerLatF = 36
  ;vecres@mpRightCornerLonF = 35

  ; ======= box 2 =========
  ;vecres@mpLeftCornerLatF  = 32
  ;vecres@mpLeftCornerLonF  = 34.5
  ;vecres@mpRightCornerLatF = 34.5
  ;vecres@mpRightCornerLonF = 36.5

  ; ====== box 3 ========
  ;vecres@mpLeftCornerLatF  = 30.2
  ;vecres@mpLeftCornerLonF  = 33.5
  ;vecres@mpRightCornerLatF = 32.2
  ;vecres@mpRightCornerLonF = 35.5

  ; ====== box 4 ========
  ;vecres@mpLeftCornerLatF  = 30.0
  ;vecres@mpLeftCornerLonF  = 31.8
  ;vecres@mpRightCornerLatF = 32.2
  ;vecres@mpRightCornerLonF = 33.5

  vecres@mpGeophysicalLineColor = "black"     ; color of continental outli
  vecres@mpUSStateLineColor     = "black"
  vecres@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed

  ; ====== modified by WC: 2010-01-27 =====
  ;vecres@vcRefAnnoOrthogonalPosF = -0.06          ; move ref vector
  ;vecres@vcRefAnnoParallelPosF   = 1.            ; move ref vector

  ;vecres@vcMonoWindBarbColor   = True
  ;vecres@vcWindBarbColor = "white"

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

  ; ===================
  map = gsn_csm_contour_map(wks,veg,res)
  plot_vec =  gsn_csm_vector_map(wks,u10(0,:,:),v10(0,:,:),vecres)

  draw(plot_vec)

  ;overlay(map,plot_vec)
  ;draw(map)
  ;draw(plot_vec)

  ;map = gsn_vector_scalar_map(wks,u10(0,:,:),v10(0,:,:),veg,res)

  frame(wks) 
end
EOF

ncl test.ncl

convert -geometry 1200x1200 -density 300 -trim raster.pdf raster.jpg
mv raster.jpg tmp_${LEV}_${TIME_INDEX}_${DOMAIN}.jpg
