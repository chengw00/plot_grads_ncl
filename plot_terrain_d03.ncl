;*************************************************
; raster_3.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  ff = addfile("wrfout_d03_2009-06-01_02:00:00.nc","r")
  xlat      = ff->XLAT
  xlon      = ff->XLONG
  hgt       = ff->HGT

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

  hgt@lat2d = xlat(0,:,:)
  hgt@lon2d = xlon(0,:,:)

  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.

  ; ==== WC: define my own colors ========
  ;cmap = RGBtoCmap("color_yellow_green_rh.txt")
  ;gsn_define_colormap(wks,cmap)
  ; ============================

  res                 = True                ; Plot mods desired.
  res@cnFillOn         = True
  ;res@gsnMaximize     = True                ; Maximize plot in frame.

  res@gsnSpreadColors     = True          

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

  ;res@gsnSpreadColors     = True            ; Use full range of colormap.
  ;res@gsnSpreadColorStart = 16              ; Start color in colormap.
  ;res@gsnSpreadColorEnd   = 190             ; End color in colormap.
  res@cnInfoLabelOn       = False           ; Turn off info label.

  res@mpGeophysicalLineColor = "black"     ; color of continental outlines
  res@mpNationalLineColor     = "black"
  res@mpUSStateLineColor     = "black" 
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
  res@cnMinLevelValF       = -500              ; min contour level
  res@cnMaxLevelValF       = 2500              ; max contour level
  res@cnLevelSpacingF      = 250.             ; contour spacing

  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels  = (/15, 30, 60, 70, 80, 90, 100/)


  res@pmLabelBarWidthF     = 0.6
  res@lbLabelFontHeightF   = 0.013
  res@lbLabelStride        = 3
  res@lbBoxLinesOn         = False

  ;res@gsnFrame     =  False                   ; do not advance the frame
  ; ===================
  map = gsn_csm_contour_map(wks,hgt(0,:,:),res)

  frame(wks) 
end
