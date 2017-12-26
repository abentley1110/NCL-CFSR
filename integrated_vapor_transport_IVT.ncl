; ===========================================
; Plots integrated water vapor transport (IVT), etc. from NCEP CFSR
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
load "~/ncl/lib/IVT_util.ncl"

; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "IVT_700hPa"   ; image name

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

    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    q_file = addfile("/cfsr/data/"+yyyy+"/q."+yyyy+".0p5.anl.nc","r")
    g_file = addfile("/cfsr/data/"+yyyy+"/g."+yyyy+".0p5.anl.nc","r")

	g = g_file->g({date},{700},:,:)  
	g = g/10
	
;---------------------------------------------------------

	; Calculating IVT and moisture flux convergence

	q = q_file->q({date},{1000:300},:,:)  
	u = u_file->u({date},{1000:300},:,:)  
	v = v_file->v({date},{1000:300},:,:) 

	IVT = uvIVT(u,v,q)  ;calculate IVT
	uIVT = IVT(0,:,:)   ;u-comp of IVT
	vIVT = IVT(1,:,:)    ; v-comp of IVT

	mag_IVT = (uIVT^2.+vIVT^2.)^0.5   ;calculate magnitude of IVT
	copy_VarCoords(uIVT,mag_IVT)
	mag_IVT@units = "kg/m/s"
	mag_IVT@long_name = "magnitude of IVT"

	;IVT_conv = -1.0*uv2dv_cfd(uIVT,vIVT,uIVT&lat,vIVT&lon,0)
	;IVT_conv = IVT_conv*(10^3)
	;copy_VarCoords(uIVT,IVT_conv)
    ;printMinMax(IVT_conv,True)
	
;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)
gsn_define_colormap(wks,"ncl_default_with_gray")  
;**********************************************
res = True

res@gsnFrame = False                                             ; DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"                              ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"
res@mpFillOn = True
res@mpOutlineOn = True

res@cnFillOn = True                                              
res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLinesOn = True                               				; Shading IVT information
res@cnLineThicknessF 	= 1
res@cnLineColor 		= "grey30"
res@cnLevels = (/250,300,400,500,600,700,800,1000,1200,1400,1600/)
res@cnFillColors = (/-1,140,152,163,173,183,194,202,240,7,8,9,10,0/)

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

res@gsnLeftString = "Integrated Water Vapor Transport (IVT); 700-hPa geo. height"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 2.0		; Map information
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

gres = True									; Geopotential height information

gres@gsnFrame = False                                        
gres@gsnDraw = False

gres@gsnLeftString = " "
gres@gsnRightString = " "

gres@cnInfoLabelOn = False                    
gres@cnFillOn = False
gres@cnLineLabelsOn = True
gres@cnLinesOn = True
gres@cnLabelMasking = True
gres@cnLineThicknessF = 4.5                   
gres@cnLevelSelectionMode = "ManualLevels"
gres@cnLevelSpacingF = 3
gres@cnMinLevelValF = 90
gres@cnMaxLevelValF = 420  
gres@cnLineLabelAngleF = 0.0
gres@cnLineLabelInterval = 2
gres@cnLineLabelDensityF = .8
gres@cnLineLabelBackgroundColor = -1
gres@cnLineLabelFontHeightF = 0.008
gres@cnLineColor = "grey15"
gres@cnLineLabelFontColor = gres@cnLineColor
gres@cnSmoothingOn = True
gres@cnSmoothingDistanceF = 0.001

;gres@gsnAddCyclic =  True
    
; *************************************************************

    vres = True							; IVT vector information

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
    vres@vcMinMagnitudeF = 0.0
    vres@gsnAddCyclic =  True
    
    if domain .eq. "northamer" then
    vres@vcRefAnnoOn = True 
    vres@vcRefAnnoString1 = "1000 kg m~S~~F19~@~F21~1"+" ~N~"+"s~S~~F19~@~F21~1"
    vres@vcRefLengthF = 0.035
    vres@vcRefMagnitudeF = 1000.0 
    vres@vcRefAnnoSide = "Top"
    vres@vcRefAnnoOrthogonalPosF = -1.046
    vres@vcRefAnnoParallelPosF = .999
    vres@vcRefAnnoString2On = False
    end if
    
    if domain .eq. "conus" then
    vres@vcRefAnnoOn = True 
    vres@vcRefAnnoString1 = "1000 kg m~S~~F19~@~F21~1"+" ~N~"+"s~S~~F19~@~F21~1"
    vres@vcRefLengthF = 0.033
    vres@vcRefMagnitudeF = 1000.0 
    vres@vcRefAnnoSide = "Top"
    vres@vcRefAnnoOrthogonalPosF = -1.047
    vres@vcRefAnnoParallelPosF = .999
    vres@vcRefAnnoString2On = False
    vres@vcMinDistanceF = 0.035
    end if
    
    vres@gsnAddCyclic =  True


; *************************************************************
res@mpFillDrawOrder = "PreDraw"					; Variable draw order
res@cnFillDrawOrder = "PreDraw"
res@cnLineDrawOrder ="PreDraw" 
gres@cnLineDrawOrder ="Draw" 
vres@vcVectorDrawOrder = "PostDraw"
res@mpOutlineDrawOrder = "PostDraw"
res@mpGridAndLimbDrawOrder = "PostDraw"



plot = gsn_csm_contour_map(wks, mag_IVT, res)          ; Plot filled contour map
contours = gsn_csm_contour(wks, g, gres)
vect = gsn_csm_vector(wks, uIVT, vIVT, vres)
plot = add_map_tickmarks(wks,plot,tmres)

;********************************************

overlay(plot,contours)
overlay(plot,vect)

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
