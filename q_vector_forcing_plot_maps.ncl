; ===========================================
; Plots 850-hPa q-vectors and forcing for ascent from NCEP CFSR
; Last updated: 26 December 2017
; Created by: Alicia M. Bentley
; Email: abentley1110@gmail.com
; ===========================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "~/ncl/lib/ut_string.ncl"
load "~/ncl/lib/latlon.ncl"
load "~/ncl/lib/shea_util.ncl"
load "~/ncl/lib/tTick.ncl"
; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "850hPa_q_vector"   ; image name

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
	
    q1_file = addfile("q1_"+foldername+".nc","r")
    q2_file = addfile("q2_"+foldername+".nc","r")
    dv_file = addfile("dv_"+foldername+".nc","r")
    
    q1 = q1_file->q1({date},{850},:,:)     	; Selecting 850 hPa here, but could load all levels and calculate
    q2 = q2_file->q2({date},{850},:,:)		; a layer average as well (using dim_avg_n_Wrap)
    dv = dv_file->dv({date},{850},:,:)
    
    u_file = addfile("/cfsr/data/"+byyyy+"/u."+byyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+byyyy+"/v."+byyyy+".0p5.anl.nc","r")
    g_file = addfile("/cfsr/data/"+byyyy+"/g."+byyyy+".0p5.anl.nc","r")
    t_file = addfile("/cfsr/data/"+byyyy+"/t."+byyyy+".0p5.anl.nc","r")
    
    u = u_file->u({date},{850},:,:)
    v = v_file->v({date},{850},:,:)
    u = u * 1.94384449 
  	v = v * 1.94384449 
  	
  	g = g_file->g({date},{850},:,:)
	g = g/10
  	
	t = t_file->t({date},{850},:,:)
	t = t - 273.15
  	copy_VarCoords(g, t)

;-------------------------------------------------

 	do d=0,1
 		g=smth9(g,.5,.25,True)			; Smooth 850-hPa geo. height and temp fields
    	t=smth9(t,.5,.25,True)
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
          (/0.92,0.52,0.91/),\
          (/0.87,0.31,0.87/),\
          (/0.85,0.64,0.88/),\
          (/0.76,0.44,0.83/),\
          (/0.36,0.35,0.75/),\
          (/0.44,0.57,0.86/),\
          (/0.62,0.72,0.90/),\
          (/0.76,0.86,0.90/),\ first blue
          (/0.99,0.95,0.60/),\
          (/0.99,0.79,0.45/),\
          (/0.96,0.65,0.39/),\
          (/0.96,0.58,0.35/),\
          (/0.96,0.38,0.26/),\
          (/1.00,0.68,0.69/),\
          (/0.63,0.02,0.02/),\
          (/0.77,0.48,0.56/)/)      
          
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

res@cnFillOn             = True        			; Q-vector divergence shading information      
res@cnFillDrawOrder = 	"Draw"
res@cnLinesOn            = True           
res@cnLineThicknessF 	= 1
res@cnLineColor 			= "grey30"
res@cnLineLabelsOn       = False          
res@gsnSpreadColors      = False                
res@cnLevelSelectionMode = "ExplicitLevels"   
res@cnFillColors = (/4,5,7,8,9,-1,10,11,13,14,15/)
res@cnLevels = (/-4,-3,-2,-1,-0.5,0.5,1,2,3,4/)
res@cnSmoothingOn = True
res@cnSmoothingTensionF = -2.
res@cnSmoothingDistanceF = 0.001

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

res@tiMainString = " "							; Image labels
res@tiMainJust = "CenterCenter"
res@tiMainPosition = "Center" 
res@tiMainOffsetYF = .02
res@tiXAxisFontHeightF = .015

res@gsnLeftString = "850-hPa geo. height, temp., wind, Q-vectors, and QG forcing"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 2.0			; Map information
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

gres = True								; Geopotential height information

gres@gsnFrame = False                                        
gres@gsnDraw = False

gres@gsnLeftString = " "
gres@gsnRightString = " "

gres@cnInfoLabelOn = False                   
gres@cnFillOn = False
gres@cnLineLabelsOn = True
gres@cnLinesOn = True
gres@cnLabelMasking = True
gres@cnLineThicknessF = 3.5                   
gres@cnLevelSelectionMode = "ManualLevels"
gres@cnLevelSpacingF = 3
gres@cnMinLevelValF = 90
gres@cnMaxLevelValF = 420 
gres@cnLineLabelAngleF = 0.0
gres@cnLineLabelInterval = 2
gres@cnLineLabelDensityF = .6
gres@cnLineLabelBackgroundColor = -1
gres@cnLineLabelFontHeightF = 0.008
gres@cnLineLabelFontColor = "mediumblue"
gres@cnLineColor = "mediumblue"
gres@cnSmoothingOn = True
gres@cnSmoothingDistanceF = 0.001
gres@gsnAddCyclic =  True

; *************************************************************

tres = True								; Temperature information

tres@gsnFrame = False                                        
tres@gsnDraw = False

tres@gsnLeftString = " "
tres@gsnRightString = " "

tres@cnInfoLabelOn = False                    
tres@cnFillOn = False
tres@cnLineLabelsOn = True
tres@cnLinesOn = True
tres@cnLabelMasking = True
tres@cnLineThicknessF = 3.5                   
tres@cnLevelSelectionMode = "ManualLevels"
tres@cnLevelSpacingF = 4
tres@cnMinLevelValF = -40
tres@cnMaxLevelValF = 40  
tres@cnLineLabelAngleF = 0.0
tres@cnLineLabelInterval = 2
tres@cnLineLabelDensityF = .6
tres@cnLineLabelBackgroundColor = -1
tres@cnLineLabelFontHeightF = 0.008
tres@cnLineLabelFontColor = "red"
tres@cnLineColor = "red"
tres@cnSmoothingOn = True
tres@cnSmoothingDistanceF = 0.001
tres@cnLineDashPattern = 5
tres@gsnAddCyclic =  True
    
; *************************************************************

    vres = True					; Q-vector information

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
    vres@vcMinDistanceF = 0.030
    vres@vcMinMagnitudeF = 1
    vres@vcMaxMagnitudeF = 16
    vres@gsnAddCyclic =  True
    
    vres@vcRefAnnoOn = True 
    vres@vcRefMagnitudeF = 5
    vres@vcRefAnnoString1 = "$VMG$ ~F34~4~F21~ 10~S~~F19~@~F21~10~N~ K m~S~~F19~@~F21~1~N~ "+"s~S~~F19~@~F21~1~N~"
    vres@vcRefLengthF = 0.035
    vres@vcRefAnnoSide = "Top"
    vres@vcRefAnnoString2On = False
    
    if domain .eq. "northamer"
        vres@vcRefAnnoOrthogonalPosF = -1.0465
    	vres@vcRefAnnoParallelPosF = 0.1890    
	end if
	
	if domain .eq. "conus"
        vres@vcRefAnnoOrthogonalPosF = -1.0465
    	vres@vcRefAnnoParallelPosF = 0.1890  
	end if
	
; *************************************************************

    bres = True						; Wind barb information

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
    bres@vcMinDistanceF = 0.045
    bres@vcWindBarbTickLengthF = 0.37
    bres@vcWindBarbTickSpacingF = 0.130
    bres@vcWindBarbTickAngleF = 55.0
    bres@gsnFrame = False
    bres@gsnAddCyclic = False

; *************************************************************
res@mpFillDrawOrder = "PreDraw"
res@cnFillDrawOrder = "PreDraw" 
gres@cnLineDrawOrder ="Draw" 
tres@cnLineDrawOrder ="Draw" 
vres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpGridAndLimbDrawOrder = "PostDraw"

plot = gsn_csm_contour_map(wks, dv, res)          ; Plot filled contour map
other = gsn_csm_contour(wks, t, tres)
contours = gsn_csm_contour(wks, g, gres)
vect = gsn_csm_vector(wks, q1 ,q2, vres)
barbs = gsn_csm_vector(wks, u, v, bres)
plot = add_map_tickmarks(wks,plot,tmres)				; Plot tick marks
;********************************************


overlay(plot,other)
overlay(plot,contours)
overlay(plot,barbs)
overlay(plot,vect)

;********************************************

draw(plot)					; Draw final image
frame(wks)

;********************************************

system("convert "+img_folder+name+"_"+number+".png -trim +repage "+dir_out+name+"_"+number+".gif")
system("'rm' "+img_folder+name+"_"+number+".png")

number = number + 1

end do     ; end times loop
print(" ")
end do     ; end domain loop

print(" ")
print("**************")
print(" END SCRIPT ")
print("**************")
print(" ")

end