; ===========================================
; Plots DT (2-PVU) theta and winds
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
; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "dt_2pvu_relvort"   ; image name

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

    p_file = addfile("/cfsr/data/"+yyyy+"/pres_pv."+yyyy+".0p5.anl.nc","r")
    t_file = addfile("/cfsr/data/"+yyyy+"/t_pv."+yyyy+".0p5.anl.nc","r")
    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    upv_file = addfile("/cfsr/data/"+yyyy+"/u_pv."+yyyy+".0p5.anl.nc","r")
    vpv_file = addfile("/cfsr/data/"+yyyy+"/v_pv."+yyyy+".0p5.anl.nc","r")

    temp = t_file->t_pv({date},{2e-06},:,:)
    pres = p_file->pres_pv({date},{2e-06},:,:)

    theta = temp ; Keep the metadata
    theta = temp * (100000/pres)^(0.286)         ; Calculates DT theta from pressure and temp.

    u = u_file->u({date},{925:850},:,:)
    v = v_file->v({date},{925:850},:,:)

    upv = upv_file->u_pv({date},{2e-06},:,:)
    vpv = vpv_file->v_pv({date},{2e-06},:,:)
    upv = upv * 1.94384449 
    vpv = vpv * 1.94384449 

    vortdiv = uv2vrG_Wrap(u,v) ; Calculate relative vorticity using spherical harmonics
    relvort = dim_avg_n_Wrap(vortdiv,0) ; Makes layer averaged vorticity

;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)
gsn_define_colormap(wks,"alicia")
;**********************************************

res = True

res@gsnFrame = False                                             ;DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"                              ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"

res@cnFillOn = True                                              ; Turn on color fill
res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnLinesOn = False                               ; Shading DT theta information
res@cnLevelSpacingF = 1
res@cnLevelSpacingF = 6                             
res@cnMinLevelValF = 264
res@cnMaxLevelValF = 384                            
res@cnFillColors = (/24,25,26,27,7,9,10,11,12,13,14,16,17,18,19,20,21,22,23,35,37,38,34,31,0/)

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

res@tiMainString = " "						; Create image labels
res@tiMainJust = "CenterCenter"
res@tiMainPosition = "Center" 
res@tiMainOffsetYF = .02
res@tiXAxisFontHeightF = .015

res@gsnLeftString = "Theta and wind (DT = 2 PVU); 925~F15~-~F21~850-hPa rel. vort."
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 5.5		; Define map outlines
res@mpNationalLineThicknessF = 5.5
res@mpUSStateLineThicknessF = 5.5
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@mpGeophysicalLineColor = "grey29"
res@mpNationalLineColor = "grey29"
res@mpUSStateLineColor = "grey29"

; *********************************************

res@mpGridAndLimbOn = True					; Define lat lon labels and grid lines
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
relres = True							; Relative vorticity information

relres@gsnFrame = False                                        
relres@gsnDraw = False

relres@cnFillMode = "AreaFill"
relres@gsnLeftString = " "
relres@gsnRightString = " "
relres@cnInfoLabelOn = False                    
relres@cnFillOn = False
relres@cnLinesOn = True
relres@cnLabelMasking = True
relres@cnLineThicknessF = 4
relres@cnLineColor = "black"
;relres@cnLineColor = "navy"
relres@cnLevelSelectionMode = "ManualLevels"
relres@cnLevelSpacingF = .5*10^-4
relres@cnMinLevelValF = .5*10^-4
relres@cnMaxLevelValF = 10*10^-4
relres@cnLineLabelsOn = False
relres@cnLineLabelAngleF = 0.0
relres@cnLineLabelInterval = 1
relres@cnLineLabelDensityF = .77
relres@cnLineLabelBackgroundColor = -1
relres@cnLineLabelFontHeightF = 0.018

relres@gsnAddCyclic = True
; *************************************************************

    bres = True							; Wind barb information

    bres@gsnFrame = False                                        
    bres@gsnDraw = False

    bres@gsnLeftString = " "
    bres@gsnRightString = " "

    bres@vcGlyphStyle = "WindBarb"
    bres@vcWindBarbColor = "grey20"
    bres@vcWindBarbCalmCircleSizeF = 0
    bres@vcWindBarbLineThicknessF = 3.0
    bres@vcRefAnnoOn = False
    bres@vcRefLengthF = 0.0242
    bres@vcRefMagnitudeF = 10.0
    bres@vcRefAnnoString2On = False
    bres@vcRefAnnoOrthogonalPosF = -1.08
    bres@vcMinDistanceF = 0.035
    bres@vcWindBarbTickLengthF = 0.37
    bres@vcWindBarbTickSpacingF = 0.130
    bres@vcWindBarbTickAngleF = 55.0
    bres@gsnFrame = False
    
    bres@gsnAddCyclic = True

; *************************************************************
res@mpGridAndLimbDrawOrder = "PostDraw"					; Define draw order
res@cnFillDrawOrder = "PreDraw"
relres@cnLineDrawOrder ="Draw" 
bres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpFillOn = False
res@mpOutlineOn = True

plot = gsn_csm_contour_map(wks, theta, res)         	; Plot filled contour map
other = gsn_csm_contour(wks, relvort, relres)		 	; Plot relative vorticity
barbs = gsn_csm_vector(wks, upv, vpv, bres)				; Plot DT winds
plot = add_map_tickmarks(wks,plot,tmres)				; Plot tickmarks

;********************************************

overlay(plot,other)			; Overlay variables
overlay(plot,barbs)

;********************************************

draw(plot)	; Draw final image
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
