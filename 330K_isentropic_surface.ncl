; ===========================================
; Plots 330 K PV, wind, pressure from NCEP CFSR
; Last updated: 26 December 2017
; Created by: Alicia M. Bentley
; Email: abentley1110@gmail.com
; ===========================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "~/ncl/lib/tTick.ncl"
load "~/ncl/lib/ut_string.ncl"
load "~/ncl/lib/latlon.ncl"
load "~/ncl/lib/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "330K_isentropic_PV"   ; image name

foldername = "superstorm93"  ; event name

folder = (/"northamer","conus"/)	; domain(s) to plot 
enddomain = dimsizes(folder)

do x = 0,enddomain-1       ; loop through all the domains you selected (only 1 domain if 1 domain selected)

domain = folder(x)

img_create = "/lb13/abentley/scripts/github"   ; Define directory where images will be initially made (and create it)
system("mkdir "+img_create+"/"+foldername)
system("mkdir "+img_create+"/"+foldername+"/"+name)
img_folder = img_create+"/"+foldername+"/"+name+"/"

img_loc = "/lb11/abentley/www/research_images"   ; Define directory where images will ultimately be located (and create it)
system("mkdir "+img_loc+"/"+foldername)
system("mkdir "+img_loc+"/"+foldername+"/"+domain)
system("mkdir "+img_loc+"/"+foldername+"/"+domain+"/"+name)
dir_out = img_loc+"/"+foldername+"/"+domain+"/"+name+"/"

print(" ")
print("Plotting '"+name+"' over '"+domain+"' for '"+foldername+"'")  ; Print what you are plotting
print("Images initially made: "+img_folder)
print("Images finally located: "+dir_out)

;*******************************************************
; Define what dates to plot
;*******************************************************

byyyy = 1993		; Define first date to plot (yyyy mm dd hh)
bmm = 3
bdd = 13
bhh = 6

eyyyy = 1993		; Define last date to plot (yyyy mm dd hh)
emm = 3
edd = 14
ehh = 0

timeUnits = "hours since 1800-01-01 00:00:00"
bdate = cd_inv_calendar(byyyy,bmm,bdd,bhh,00,00,timeUnits,0)
edate = cd_inv_calendar(eyyyy,emm,edd,ehh,00,00,timeUnits,0)

print(" ")
print("bdate: "+ut_string(bdate,""))
print("edate: "+ut_string(edate,""))
print(" ")

number = 1   		; Number of the first image to be created

do date = bdate,edate,6       ; Create images from start date to end date (with 6 h increments)

; ******************************************************
; Load files from NCEP CFSR
; ******************************************************
    print("Date:  "+ut_string(date,"")+" | Creating image #"+number) 	; Print date currently being plotted
	yyyy = ut_string(date,"%Y")    			; Determine year of current date in yyyy form
	
    pv_file = addfile("/cfsr/data/"+yyyy+"/pv_isen."+yyyy+".0p5.anl.nc","r")
;    p_file = addfile("/cfsr/data/"+yyyy+"/pres_pv."+yyyy+".0p5.anl.nc","r")
    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    t_file = addfile("/cfsr/data/"+yyyy+"/t."+yyyy+".0p5.anl.nc","r")
    uisen_file = addfile("/cfsr/data/"+yyyy+"/u_isen."+yyyy+".0p5.anl.nc","r")
    visen_file = addfile("/cfsr/data/"+yyyy+"/v_isen."+yyyy+".0p5.anl.nc","r")

    isen_lev = (/330/)

    pv = pv_file->pv_isen({date},{isen_lev},:,:)
    pv = pv * (10^6)

    uisen = uisen_file->u_isen({date},{isen_lev},:,:)
    visen = visen_file->v_isen({date},{isen_lev},:,:)
    uisen = uisen * 1.94384449 
    visen = visen * 1.94384449 
    
;*************************************    
  
    t = t_file->t({date},{1000:100},:,:)
    u = u_file->u({date},{1000:100},:,:)
    v = v_file->v({date},{1000:100},:,:)
  
  do il = 0,dimsizes(isen_lev)-1
  
  lats = t_file->lat    
  lons = t_file->lon
  lev = t_file->lev({1000:100})

  nlat = dimsizes(lats)
  nlon = dimsizes(lons)
  nlev = dimsizes(lev)
  
  p = t
  do l=0,nlev-1
    p(l,:,:) = lev(l)
  end do
  p@units = "hPa"

  ;Make potential temp. (all levs)
  theta = t
  theta = t*((1000./p)^(287./1004.))


;------------- Create Isentropic Pressure --------------

  th_lev = (/isen_lev(il)-1,isen_lev(il),isen_lev(il)+1/) ; Only calculate desired levels!
  nth_lev = dimsizes(th_lev)
  
  levs_todo = (/isen_lev(il)/)
  nlevs_todo = dimsizes(levs_todo)

  p_th = new((/nth_lev,nlat,nlon/),"float")

  do j=0,nth_lev-1,1
    p_th(j,:,:)  = (/wrf_user_intrp3d(p(:,:,:),theta(:,:,:),"h",th_lev(j),0.,False)/)
  end do

  p_th!0 = "lev"
  p_th!1 = "lat"
  p_th!2 = "lon"
    
  p_th = p_th

  p_th&lev = th_lev
  p_th&lat = lats
  p_th&lon = lons
  p_th@_FillValue = -99999.

	do d=0,2
		p_th=smth9(p_th,.5,.25,True)  	; Smoothes 330 K pressure field
	end do

;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)
gsn_define_colormap(wks,"lance_isenpv")
;**********************************************

res = True

res@gsnFrame = False                            ;DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"             ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"
res@mpFillOn = False
res@mpOutlineOn = True

res@cnFillOn = True                                              
res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn = False                          	; Shaded 330K PV information
res@cnLevelSpacingF = 1                            
res@cnLevels := (/0,.5,1,1.5,2,4,6,8,10/)
res@cnFillColors := (/6,-1,14,13,12,16,19,22,33,26/)

;**********************************************

res@lbAutoManage = False                   ; Label Bar Information
res@pmLabelBarOrthogonalPosF = .02  
res@pmLabelBarHeightF = 0.08
res@pmLabelBarWidthF = 0.73
res@lbLabelAutoStride = True
res@lbLabelFontHeightF = 0.0124
res@lbBoxLinesOn  = True                    
res@lbTopMarginF = .3

;**********************************************

res@tiMainString = " "						; Image labels
res@tiMainJust = "CenterCenter"
res@tiMainPosition = "Center" 
res@tiMainOffsetYF = .02
res@tiXAxisFontHeightF = .015

res@gsnLeftString = "PV, pressure, and wind (330 K isen. sfc.)"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 4.0		; Map information
res@mpNationalLineThicknessF = 4.0
res@mpUSStateLineThicknessF = 3.5
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@mpGeophysicalLineColor = "navyblue"
res@mpNationalLineColor = "navyblue"
res@mpUSStateLineColor = "navyblue"

; *********************************************
res@mpGridAndLimbOn = True					; Grid lines and tick marks
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 20.
res@mpGridLonSpacingF = 20.
res@mpGridLineThicknessF = 3
res@mpGridLineColor = "grey29"
res@pmTickMarkDisplayMode = "Always"

tmres = True
tmres@tmXBLabelFontHeightF = .013
tmres@tmYLLabelFontHeightF = .013
tmres@tmYLLabelFontHeightF = .013

tmres@tmXBMajorLengthF = .01
tmres@tmXBMinorLengthF = .006

tmres@tmYRMinorValues = (/5,10,15,25,30,35,45,50,55,65/)
tmres@tmYLMinorValues = (/5,10,15,25,30,35,45,50,55,65/)
tmres@tmXBMinorValues = (/-155,-150,-145,-135,-130,-125,-120,-115,-110,-105,-95,-90,-85,-80,-75,-70,-65,-55,-50,-45,-35,-30,-25,-15,-10,-5,5,10,15,25,30,35/)


; *************************************************************

if domain .eq. "northamer"					; Define domains	
res@mpProjection = "LambertConformal"
res@mpLambertParallel2F = 31
res@mpLambertParallel1F = 11
res@mpLambertMeridianF = -100
res@mpLimitMode = "LatLon"
res@mpMinLatF = 8
res@mpMaxLatF = 69
res@mpMinLonF = -151
res@mpMaxLonF = -49
tmres@tmXBValues := ispan(-180,179,40)
tmres@tmYLValues := ispan(0,80,20)
tmres@tmYRValues := ispan(0,80,20)
end if


if domain .eq. "conus" then
res@mpProjection = "LambertConformal"
res@mpLambertParallel2F = 11
res@mpLambertParallel1F = 5
res@mpLambertMeridianF = -100
res@mpLimitMode = "LatLon"
res@mpMinLatF = 20
res@mpMaxLatF = 52
res@mpMinLonF = -125
res@mpMaxLonF = -67
tmres@tmXBValues := ispan(-180,179,20)
tmres@tmYLValues := ispan(0,80,20)
tmres@tmYRValues := ispan(0,80,20)
end if

; *************************************************************

cnres = True						; 330 K pressure contour information

cnres@gsnFrame = False                                        
cnres@gsnDraw = False

cnres@gsnLeftString = " "
cnres@gsnRightString = " "

cnres@cnFillMode = "AreaFill"
cnres@cnInfoLabelOn = False               
cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLabelMasking = True
cnres@cnLineThicknessF = 4.0                   
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 30
cnres@cnMinLevelValF = 90
cnres@cnMaxLevelValF = 600    
cnres@cnLineLabelAngleF = 0.0
cnres@cnLineLabelInterval = 2
cnres@cnLineLabelDensityF = .70
cnres@cnLineLabelBackgroundColor = -1
cnres@cnLineLabelFontHeightF = 0.008

cnres@gsnAddCyclic =  True

; *************************************************************


    bres = True						; Wind barb information

    bres@gsnFrame = False                                        
    bres@gsnDraw = False

    bres@gsnLeftString = " "
    bres@gsnRightString = " "

    bres@vcGlyphStyle = "WindBarb"
    bres@vcWindBarbColor = "grey14"
    bres@vcWindBarbCalmCircleSizeF = 0
    bres@vcWindBarbLineThicknessF = 3.0
    bres@vcRefAnnoOn = False
    bres@vcRefLengthF = 0.0242
    bres@vcRefMagnitudeF = 10.0
    bres@vcRefAnnoString2On = False
    bres@vcRefAnnoOrthogonalPosF = -1.08
    bres@vcMinDistanceF = 0.040
    bres@vcWindBarbTickLengthF = 0.37
    bres@vcWindBarbTickSpacingF = 0.130
    bres@vcWindBarbTickAngleF = 55.0
    bres@gsnFrame = False
    
    bres@gsnAddCyclic =  True


; *************************************************************
res@mpGridAndLimbDrawOrder = "PostDraw"		; Define variable draw order
res@cnFillDrawOrder = "PreDraw"
bres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"

plot = gsn_csm_contour_map(wks, pv, res)          			; Plot filled 330K PV map
contours = gsn_csm_contour(wks, p_th({330},:,:), cnres)		; Plot 330 K pressure contours
barbs = gsn_csm_vector(wks, uisen, visen, bres)				; Plot wind barbs
plot = add_map_tickmarks(wks,plot,tmres)					; Plot tick marks

;********************************************

overlay(plot,contours)
overlay(plot,barbs)

;********************************************

draw(plot)					; Draw final image
frame(wks)

;********************************************

system("convert "+img_folder+name+"_"+number+".png -trim +repage "+dir_out+name+"_"+number+".gif")
system("'rm' "+img_folder+name+"_"+number+".png")

number = number + 1

end do   	; end isen. pressure loop  
end do      ; end date loop
print(" ")
end do		; end domain loop

print(" ")
print("**************")
print(" END SCRIPT ")
print("**************")
print(" ")

end