;*************************************************
; raster_t.ncl
;************************************************
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "/usr/local/ncarg-5.1.1/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  setfileoption("bin","ReadByteOrder","LittleEndian")
  ff = addfile("wrfout_d03_2009-06-01_05:00:00.nc","r")

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

  ; =========================
  qv = ff->QVAPOR
  qcloud = ff->QCLOUD
  qrain  = ff->QRAIN
  qice   = ff->QICE
  qsnow  = ff->QSNOW
  qgraup = ff->QGRAUP

  qhydromet = qcloud + qrain + qice + qsnow + qgraup

  ; ===== read 2d file =======
  kplot = 11 - 1

  template = new((/n_y,n_x/),float,-999999)
  qh_total_path = template

  qh_total_path = 1000*qhydromet(0,kplot,:,:)

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
  res@cnLevels  = (/0.01, 0.1, 0.2, 0.4, 0.8, 1.0, 1.2, 2.0/)

  ; ====== WC : change map thickness
  res@mpGeophysicalLineThicknessF = 3
  res@mpNationalLineThicknessF = 3
  ; ====== WC: define title ====
  res@tiMainFontHeightF = 0.02
  res@tiMainString = "2009-06-01_12"
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
