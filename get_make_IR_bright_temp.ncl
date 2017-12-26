; ===========================================
; Plots IR brightness temp, etc. from NCEP CFSR
; IR brightness temps only available through 2016
; Last updated: 26 December 2017
; Created by: Alicia M. Bentley
; Email: abentley1110@gmail.com
; ===========================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "~/ncl/lib/tTick.ncl"
load "~/ncl/lib/ut_string.ncl"
load "~/ncl/lib/latlon.ncl"
load "~/ncl/lib/shea_util.ncl"
; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "IR_brightness_temps"

foldername = "superstorm93"  ; event name

folder = (/"northamer","conus"/)	; domain(s) to plot 
enddomain = dimsizes(folder)

do x = 0,enddomain-1       ; loop through all the domains you selected (only 1 domain if 1 domain selected)

domain = folder(x)

img_create = "/lb13/abentley/scripts/github"   ; Define directory where images will be initially made (and create it)
system("mkdir "+img_create+"/"+foldername)
system("mkdir "+img_create+"/"+foldername+"/"+name)
img_folder = img_create+"/"+foldername+"/"+name+"/"

img_loc = "/lb11/abentley/www/research_images"   ; Define directory where images will uldateately be located (and create it)
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

dateeUnits = "hours since 1800-01-01 00:00:00"
bdate = cd_inv_calendar(byyyy,bmm,bdd,bhh,00,00,dateeUnits,0)
edate = cd_inv_calendar(eyyyy,emm,edd,ehh,00,00,dateeUnits,0)

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
	    
	year = cd_string(date,"%Y")
   	mon = cd_string(date,"%N")
   	day = cd_string(date,"%D")
   	hr = cd_string(date,"%H")
   	dayn = cd_string(date,"%J")

  print("https://www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/"+year+"/GRIDSAT-B1."+year+"."+mon+"."+day+"."+hr+".v02r01.nc")
  
if domain .ne. "conus" then  
  system("'rm' GRIDSAT-B1."+year+"."+mon+"."+day+"."+hr+".v02r01.nc")
  system("wget https://www.ncei.noaa.gov/data/geostationary-ir-channel-brightness-temperature-gridsat-b1/access/"+year+"/GRIDSAT-B1."+year+"."+mon+"."+day+"."+hr+".v02r01.nc")
	print(" ")
	print("Successfully downloaded GRIDSAT-B1."+year+"."+mon+"."+day+"."+hr+".v02r01.nc!")
end if  

; ******************************************************
; Load files
; ******************************************************

    ir_file = addfile("GRIDSAT-B1."+year+"."+mon+"."+day+"."+hr+".v02r01.nc","r")
    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    t_file = addfile("/cfsr/data/"+yyyy+"/t."+yyyy+".0p5.anl.nc","r")
    w_file = addfile("/cfsr/data/"+yyyy+"/w."+yyyy+".0p5.anl.nc","r")
    
    ir = short2flt(ir_file->irwin_cdr(0,:,:))
    ir = ir - 273.15

    lev_pv  = (/200,250,300/)
    lev_pv@units = "hPa"
    upv = u_file->u({date},{lev_pv},:,:)
    vpv = v_file->v({date},{lev_pv},:,:)
    tpv = t_file->t({date},{lev_pv},:,:)

    pv_all = PotVortIsobaric((/20000,25000,30000/),upv,vpv,tpv,tpv&lat,1,0)
    pv = dim_avg_n_Wrap(pv_all,0)
    pv = pv*(10^6)
    copy_VarCoords(upv(0,:,:),pv)
    
    w_orig := w_file->w({date},{600:400},:,:)
    w_temp := w_orig/100          ; convert from Pa/s to hPa/s
    w := dim_avg_n_Wrap(w_temp,0) ; Averages over the specified layer
    copy_VarCoords(upv(0,:,:),w)

    uirro = u_file->u({date},{300:200},:,:)
    virro = v_file->v({date},{300:200},:,:)
    div = uv2dvF_Wrap(uirro,virro)
    irr = dv2uvF_Wrap(div)
    irro = dim_avg_n_Wrap(irr,1) ; Averages over the specified layer

	do d=0,3
		pv=smth9(pv,.5,.25,True)		; Smooth PV field
	end do

;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)
gsn_define_colormap(wks,"nhc_tc")
colors = gsn_retrieve_colormap(wks)  ; retrieve color map for editing. dimensioned (102,3)
colors(32,:) = (/ 255,255,255 /)/255 ; replace the 32nd color
gsn_reverse_colormap(wks) ; reverse colormap
;**********************************************

res = True

res@gsnFrame = False                                    ;DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"                     ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"

res@cnFillOn = True                                     ; IR brightness temp plotting information
;res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
res@cnMinLevelValF       = -80                ; set min contour level
res@cnMaxLevelValF       = 50                ; set max contour level
res@cnLevelSpacingF      = 2.0              ; set contour spacing
res@lbLabelStride        = 5
          
res@gsnSpreadColors = True
res@gsnSpreadColorStart = 2
res@gsnSpreadColorEnd = 60

res@cnLineLabelInterval = 10                 ; set label intervals to every other contour
res@cnLineLabelPlacementMode = 2

res@cnFillMode      = "RasterFill"            ; Use raster to draw map faster
res@cnRasterSmoothingOn = True
res@cnRasterCellSizeF = .0002
res@cnLinesOn       = False                 ; draw contour lines
res@cnLineThicknessF  = 1.0                 ; line thickness of the contour lines
res@cnNoDataLabelOn   = False
res@cnConstFLabelPerimOn = False
res@cnLowLabelsOn        = False
res@cnLineLabelsOn       = False

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

res@gsnLeftString = "IR brightness temperature; 300~F15~-~F21~200-hPa PV/irr. wind"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 3.0		; Map information
res@mpNationalLineThicknessF = 3.0
res@mpUSStateLineThicknessF = 3.0
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@pmTickMarkDisplayMode = "Always"

res@mpGeophysicalLineColor = "black"
res@mpNationalLineColor = "black"
res@mpUSStateLineColor = "black"
;green4"

res@mpDefaultFillColor = "transparent"

; *********************************************

res@mpGridAndLimbOn = True					; Grid line and tick mark information
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 20.
res@mpGridLonSpacingF = 20.
res@mpGridLineThicknessF = 3
res@mpGridLineColor = "grey12"
res@pmTickMarkDisplayMode = "Always"

tmres = True
tmres@tmYLValues := ispan(0,80,20)
tmres@tmYRValues := ispan(0,80,20)
tmres@tmXBValues := ispan(-180,40,40)
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

pvres = True									; PV contour information

pvres@gsnFrame = False                                        
pvres@gsnDraw = False

pvres@gsnLeftString = " "
pvres@gsnRightString = " "

pvres@cnInfoLabelOn = False                    
pvres@cnFillOn = False
pvres@cnLineLabelsOn = True
pvres@cnLinesOn = True
pvres@cnLabelMasking = True
pvres@cnLineThicknessF = 4.5                    
pvres@cnLevelSelectionMode = "ManualLevels"
pvres@cnLevelSpacingF = 1
pvres@cnMinLevelValF = 1
pvres@cnMaxLevelValF = 30  
pvres@cnLineLabelAngleF = 0.0
pvres@cnLineLabelInterval = 2
pvres@cnLineLabelDensityF = .77
pvres@cnLineLabelBackgroundColor = -1
pvres@cnLineLabelFontHeightF = 0.008
pvres@cnLineLabelFontColor = "red"
pvres@cnLineColor = "red"

pvres@gsnAddCyclic =  True

; *************************************************************

omres = True									; Vertical motion information

omres@gsnFrame = False                                        
omres@gsnDraw = False

omres@gsnLeftString = " "
omres@gsnRightString = " "
omres@cnInfoLabelOn = False                    
omres@cnFillOn = False
omres@cnLinesOn = True
omres@cnLabelMasking = True
omres@cnLineLabelsOn = False
omres@cnLineThicknessF = 4
omres@cnLineColor = "red"
omres@cnLevelSelectionMode = "ManualLevels"
omres@cnLevelSpacingF = 0.005
omres@cnMaxLevelValF = -0.005 
omres@cnLineLabelAngleF = 0.0
omres@cnLineLabelInterval = 1
omres@cnLineLabelDensityF = .77
omres@cnLineLabelBackgroundColor = -1
omres@cnLineLabelFontHeightF = 0.018

omres@gsnAddCyclic =  True


; *************************************************************

    vres = True									; Irrotational wind vector information

    vres@gsnFrame = False                                        
    vres@gsnDraw = False

    vres@gsnLeftString = " "
    vres@gsnRightString = " "

    vres@vcGlyphStyle = "FillArrow"
    vres@vcWindBarbColor = "black"
    vres@vcFillArrowEdgeColor = "black"
    vres@vcFillArrowFillColor = "black"
    vres@vcWindBarbCalmCircleSizeF = 0
    vres@vcWindBarbLineThicknessF = 7.5
    vres@vcMinDistanceF = 0.026
    vres@vcMinMagnitudeF = 3.0
    vres@gsnAddCyclic =  True
    
    if domain .eq. "northamer" then
    vres@vcRefAnnoOn = True 
    vres@vcRefAnnoString1 = "$VMG$ m s~S~~F15~-~F21~1"
    vres@vcRefLengthF = 0.035
    vres@vcRefMagnitudeF = 10.0
    vres@vcRefAnnoSide = "Top"
    vres@vcRefAnnoOrthogonalPosF = -1.046
    vres@vcRefAnnoParallelPosF = .999
    vres@vcRefAnnoString2On = False
    end if
    
    if domain .eq. "conus" then
    vres@vcRefAnnoOn = True 
    vres@vcRefAnnoString1 = "$VMG$ m s~S~~F15~-~F21~1"
    vres@vcRefLengthF = 0.035
    vres@vcRefMagnitudeF = 10.0
    vres@vcRefAnnoSide = "Top"
    vres@vcRefAnnoOrthogonalPosF = -1.047
    vres@vcRefAnnoParallelPosF = .999
    vres@vcRefAnnoString2On = False
    vres@vcMinDistanceF = 0.024
    end if
    
    vres@gsnAddCyclic =  True
; *************************************************************
res@mpFillOn = False
res@mpOutlineOn = True
;res@mpFillDrawOrder = "PreDraw"
res@cnFillDrawOrder = "PreDraw"
pvres@cnLineDrawOrder ="Draw"
omres@cnLineDrawOrder ="Draw" 
vres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpGridAndLimbDrawOrder = "PostDraw"

plot = gsn_csm_contour_map(wks, ir, res)          ; Plot filled contour map
contours = gsn_csm_contour(wks, pv, pvres)
other = gsn_csm_contour(wks, w, omres)
vect = gsn_csm_vector(wks, irro(0,:,:), irro(1,:,:), vres)
plot = add_map_tickmarks(wks,plot,tmres)

;********************************************

overlay(plot,contours)
;overlay(plot,other)
overlay(plot,vect)

;********************************************

draw(plot)					; Draw final image
frame(wks)

;********************************************

system("convert "+img_folder+name+"_"+number+".png -trim +repage "+dir_out+name+"_"+number+".gif")
system("'rm' "+img_folder+name+"_"+number+".png")

number = number + 1

end do     ; end date loop
print(" ")
end do     ; end domain loop

;  system("'rm' GRIDSAT-B1.*.v02r01.nc")

print(" ")
print("**************")
print(" END SCRIPT ")
print("**************")
print(" ")

end