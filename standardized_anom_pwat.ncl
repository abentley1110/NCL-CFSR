; ===========================================
; Plots standardized PWAT anom. from NCEP CFSR
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

name = "PW_standard_anom"   ; image name

foldername = "superstorm93"  ; event name

folder = (/"northamer","conus"/)	; domain(s) to plot 
enddomain = dimsizes(folder)

do x = 0,enddomain-1       ; loop through all the domains you selected (only 1 domain if 1 domain selected)

domain = folder(x)

img_create = "/path/to/image/creation/location"   ; Define directory where images will be initially made (and create it)
system("mkdir "+img_create+"/"+foldername)
system("mkdir "+img_create+"/"+foldername+"/"+name)
img_folder = img_create+"/"+foldername+"/"+name+"/"

img_loc = "/path/to/image/final/location"   ; Define directory where images will ultimately be located (and create it)
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
	
    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    g_file = addfile("/cfsr/data/"+yyyy+"/g."+yyyy+".0p5.anl.nc","r")
    pw_file = addfile("/cfsr/data/"+yyyy+"/pwat."+yyyy+".0p5.anl.nc","r")

    u = u_file->u({date},{700},:,:)
    v = v_file->v({date},{700},:,:)
    
    g = g_file->g({date},{700},:,:)
    g = g/10 ; dam

  fmean = addfile("/cfsr/data/climo/mean/pwat.0p5.mean.climo.nc","r")
  fstdev = addfile("/cfsr/data/climo/stdev/pwat.0p5.stdev.climo.nc","r")

	dne_yyyy = stringtoint(ut_string(date,"%Y"))
	dne_mm = stringtoint(ut_string(date,"%N"))
	dne_dd = stringtoint(ut_string(date,"%D"))
	dne_hh = stringtoint(ut_string(date,"%H"))

	print("Anom. Info.: "+dne_yyyy+" "+dne_mm+" "+dne_dd+" "+dne_hh)
    sdoy = day_of_year(dne_yyyy,dne_mm,dne_dd)

	g_all = pw_file->pwat({date},:,:)
;    g_all = lonFlip(g_all)
    g_mean := fmean->pwat({sdoy},{dne_hh},:,:)
    g_stdev := fstdev->pwat({sdoy},{dne_hh},:,:)
    
  	std_anom2 := ((g_all - g_mean)/g_stdev)
  	copy_VarMeta( g_all, std_anom2 )
  	
  	u = u * 1.94384449 
  	v = v * 1.94384449 
  	
  	do d=0,1
		g=smth9(g,.5,.25,True)		; Smoothes geopotential height
	end do 

;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)

cmap =  (/(/1.00,1.00,1.00/),\
          (/0.00,0.00,0.00/),\
          (/0.92,0.72,0.91/),\
          (/0.92,0.74,0.91/),\
          (/0.92,0.52,0.91/),\
          (/0.87,0.31,0.87/),\
          (/0.72,0.07,0.82/),\
          (/0.41,0.16,0.69/),\
          (/0.24,0.22,0.65/),\
          (/0.22,0.42,0.72/),\
          (/0.43,0.60,0.81/),\
          (/0.67,0.81,0.91/),\ first blue
          (/0.99,0.93,0.45/),\
          (/0.96,0.78,0.37/),\
          (/0.96,0.63,0.32/),\
          (/0.95,0.50,0.25/),\
          (/0.96,0.31,0.18/),\
          (/0.91,0.05,0.05/),\
          (/0.63,0.02,0.02/),\
          (/0.87,0.49,0.60/),\
          (/0.43,0.17,0.17/),\
          (/0.46,0.41,0.41/)/)
          
gsn_define_colormap(wks,cmap)  
;**********************************************

res = True

res@gsnFrame = False                                             ;DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"                              ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"
res@mpFillOn = True
res@mpOutlineOn = True

res@cnFillOn = True                                              
res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ExplicitLevels"             		; Shading standardized anom. information
res@cnLinesOn = False                               
res@cnLevelSpacingF = 1                            
res@cnLevels = (/-6.0,-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-.5,0,.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0/)
res@cnFillColors = (/3,4,5,6,7,8,9,10,11,-1,-1,12,13,14,15,16,17,18,20,21/)

;**********************************************

res@lbAutoManage = False                   	; Label Bar Information
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

res@gsnLeftString = "700-hPa height, wind; PW, std. PW anom. (1979-2009 CFSR)"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 2.0		; Map Information
res@mpNationalLineThicknessF = 2.0
res@mpUSStateLineThicknessF = 2.0
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@mpGeophysicalLineColor = "grey29"
res@mpNationalLineColor = "grey29"
res@mpUSStateLineColor = "grey29"

; *********************************************

res@mpGridAndLimbOn = True					; Grid lines and tick mark information
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 20.
res@mpGridLonSpacingF = 20.
res@mpGridLineThicknessF = 3
res@mpGridLineColor = "grey12"
res@pmTickMarkDisplayMode = "Always"

tmres = True
tmres@tmXBLabelFontHeightF = .013
tmres@tmYLLabelFontHeightF = .013
tmres@tmYLLabelFontHeightF = .013

tmres@tmXBMajorLengthF = .01
tmres@tmXBMinorLengthF = .006

tmres@tmYRMinorValues = (/5,10,15,25,30,35,45,50,55,65/)
tmres@tmYLMinorValues = (/5,10,15,25,30,35,45,50,55,65/)
tmres@tmXBMinorValues = (/-175,-170,-165,-160,-155,-150,-145,-135,-130,-125,-120,-115,-110,-105,-95,-90,-85,-80,-75,-70,-65,-55,-50,-45,-35,-30,-25,-15,-10,-5,5,10,15,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,145,150,155,160,165,170,175/)

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

cnres = True								; Geopotential height information

cnres@gsnFrame = False                                        
cnres@gsnDraw = False

cnres@gsnLeftString = " "
cnres@gsnRightString = " "

cnres@cnInfoLabelOn = False                    
cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLabelMasking = True
cnres@cnLineThicknessF = 4                     
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 3
cnres@cnMinLevelValF = 90
cnres@cnMaxLevelValF = 420   
cnres@cnLineLabelAngleF = 0.0
cnres@cnLineLabelInterval = 2
cnres@cnLineLabelDensityF = .90
;cnres@cnSmoothingOn = True
;cnres@cnSmoothingDistanceF = 0.009
cnres@cnLineLabelBackgroundColor = -1
cnres@cnLineLabelFontHeightF = 0.008
cnres@gsnAddCyclic =  True

; *************************************************************

    bres = True								; Wind barb information

    bres@gsnFrame = False                                        
    bres@gsnDraw = False

    bres@gsnLeftString = " "
    bres@gsnRightString = " "

    bres@vcGlyphStyle = "WindBarb"
    bres@vcWindBarbColor = "grey12"
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
res@mpFillDrawOrder = "PreDraw"
res@cnFillDrawOrder = "PreDraw"
cnres@cnLineDrawOrder ="Draw"
bres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpGridAndLimbDrawOrder = "PostDraw"

plot = gsn_csm_contour_map(wks, std_anom2, res)          	; Plot filled standardized anom. map
contours = gsn_csm_contour(wks, g, cnres)					; Plot geopotential height
barbs = gsn_csm_vector(wks, u, v, bres)						; Plot wind barbs
plot = add_map_tickmarks(wks,plot,tmres)					; Plot tick marks

;********************************************

overlay(plot,contours)		; Overlay variable order
overlay(plot,barbs)

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

print(" ")
print("**************")
print(" END SCRIPT ")
print("**************")
print(" ")

end
