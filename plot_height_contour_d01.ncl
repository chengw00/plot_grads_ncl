;*************************************************
; raster_3.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
  ff = addfile("wrfout_d01_2009-05-31_13:00:00.nc","r")
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

  ;do j=1,n_y
  ; do i=1,n_x
  ;  ii=i-1
  ;  jj=j-1
  ;  veg(jj,ii)=i
  ; end do
  ;end do

 ; ===== read 3d pressure file =======
  n_p = 26

  template_1d = new(n_p,float,-999999)
  pres_1d = (/1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,    \
                  600, 550, 500, 450, 400, 350, 300, 250, 200, 150, \
                  100,  70,  50,  30,  20,  10/)

  print(pres_1d)

  ;plev = 1000
  ;a_min = -6
  ;a_max = 20
  ;a_int = 2

  ;plev = 850
  ;a_min = 120
  ;a_max = 168
  ;a_int = 3

  plev = 500
  a_min = 522
  a_max = 592
  a_int = 3

  do k=0,n_p-1
   if (plev.eq.pres_1d(k)) then
    kplot = k
   end if
  end do

  print(kplot)
  ; ==========================
  template = new((/n_p,n_y,n_x/),float,-999999)
  a_3d = template

  dims_3d = (/n_p,n_y,n_x/)
  nrec = 4 - 1
  file_3d = "store_3d_2009-05-31_13_d01.gdat"
  a_3d = fbindirread(file_3d, nrec, dims_3d, "float")

  a_3d@lat2d = xlat(0,:,:)
  a_3d@lon2d = xlon(0,:,:)
  a_3d@units = ""

  a_3d = a_3d/10 ; convert to dm
  print(a_3d(:,50,50))

  u_3d = template
  nrec = 1 - 1
  u_3d = fbindirread(file_3d, nrec, dims_3d, "float")
  u_3d@lat2d = xlat(0,:,:)
  u_3d@lon2d = xlon(0,:,:)
  u_3d@units = "m/s"

  v_3d = template
  nrec = 2 - 1
  v_3d = fbindirread(file_3d, nrec, dims_3d, "float")
  v_3d@lat2d = xlat(0,:,:)
  v_3d@lon2d = xlon(0,:,:)
  v_3d@units = "m/s"

  ; =================
  wks = gsn_open_wks("pdf","raster")         ; Open a workstation.
  gsn_define_colormap(wks,"BlAqGrYeOrRe")   ; Choose colormap.

  res                 = True                ; Plot mods desired.
  res@cnFillOn         = False
  ;res@gsnMaximize     = True                ; Maximize plot in frame.

  res@cnLineThicknessF  = 3.0              ; contour line thickness
  res@cnLineColor = "red"                ; contour line color
  res@cnLineLabelFontColor = "red"       ; label font color
 
  res@gsnStringFontHeightF         = 0.013
  res@gsnRightStringOrthogonalPosF = 0.02
  res@gsnLeftStringOrthogonalPosF  = 0.02

  res@cnLinesOn       = True              ; Turn off contour lines      .
  res@cnLineLabelsOn  = True              ; Turn off contour labels.

  res@mpProjection           = "LambertConformal"
  res@mpLambertParallel1F    = tlat1               ; two parallels
  res@mpLambertParallel2F    = tlat2
  res@mpLambertMeridianF     = clon             ; central meridian

  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = xlat(0,0,0)
  res@mpLeftCornerLonF  = xlon(0,0,0)
  res@mpRightCornerLatF = xlat(0,n_y-1,n_x-1)
  res@mpRightCornerLonF = xlon(0,n_y-1,n_x-1)

  ; ====== WC: define title ====
  res@tiMainFontHeightF = 0.02
  res@tiMainString = "2009-05-31_17"
  ; ============================================

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
  res@gsnSpreadColorStart = 16              ; Start color in colormap.
  res@gsnSpreadColorEnd   = 190             ; End color in colormap.
  res@cnInfoLabelOn       = False           ; Turn off info label.

  res@mpGeophysicalLineColor = "green"     ; color of continental outlines
  res@mpUSStateLineColor     = "green" 
  res@mpNationalLineColor    = "green"
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed

  res@mpGeophysicalLineThicknessF = 3.
  res@mpNationalLineThicknessF = 3.

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

  ;res@cnFillMode           = "RasterFill"
  ;res@cnFillMode           = "CellFill"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = a_min              ; min contour level
  res@cnMaxLevelValF       = a_max              ; max contour level
  res@cnLevelSpacingF      = a_int             ; contour spacing

  res@pmLabelBarWidthF     = 0.6
  res@lbLabelFontHeightF   = 0.013
  res@lbLabelStride        = 2
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

  vecres@mpGeophysicalLineColor = "green"     ; color of continental outli
  vecres@mpUSStateLineColor     = "green"
  vecres@mpNationalLineColor       = "green"
  vecres@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed

  vecres@mpOutlineBoundarySets  = "National"
  ; ====== modified by WC: 2010-01-27 =====
  ;vecres@vcRefAnnoOrthogonalPosF = -0.06          ; move ref vector
  ;vecres@vcRefAnnoParallelPosF   = 1.            ; move ref vector

  vecres@mpFillOn               = False   ; no fill

  ; =======used for Native LC projection ========
  ; not needed since lat/lon defined with field
  ;vecres@tfDoNDCOverlay         = True
  ; =======================================
  vecres@mpOutlineBoundarySets  = "National"

  vecres@mpGeophysicalLineThicknessF = 1.3
  vecres@mpNationalLineThicknessF = 1.3
  ; ======================================

  ; ====
  vecres@mpDataBaseVersion = "MediumRes"
  ;vecres@mpDataBaseVersion = "HighRes"
  vecres@mpDataSetName = "Earth..4"

  ; ===================
  map = gsn_csm_contour_map(wks,a_3d(kplot,:,:),res)
  plot_vec =  gsn_csm_vector_map(wks,u_3d(kplot,:,:),v_3d(kplot,:,:),vecres)
  draw(plot_vec)

  ;overlay(map,plot_vec)
  ;draw(map)
  ;draw(plot_vec)

  ;map = gsn_vector_scalar_map(wks,u10(0,:,:),v10(0,:,:),veg,res)

  frame(wks) 
end
