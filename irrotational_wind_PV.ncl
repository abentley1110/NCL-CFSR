; ===========================================
; Plots PV, irrotational wind, etc. from NCEP CFSR
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

; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "irro_wind_PV"

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
	
    pwat_file = addfile("/cfsr/data/"+yyyy+"/pwat."+yyyy+".0p5.anl.nc","r")
    w_file = addfile("/cfsr/data/"+yyyy+"/w."+yyyy+".0p5.anl.nc","r")
    u_file = addfile("/cfsr/data/"+yyyy+"/u."+yyyy+".0p5.anl.nc","r")
    v_file = addfile("/cfsr/data/"+yyyy+"/v."+yyyy+".0p5.anl.nc","r")
    t_file = addfile("/cfsr/data/"+yyyy+"/t."+yyyy+".0p5.anl.nc","r")

    pwat = pwat_file->pwat({date},:,:)

    lev_pv  = (/200,250,300/)
    lev_pv@units = "hPa"
    upv = u_file->u({date},{lev_pv},:,:)
    vpv = v_file->v({date},{lev_pv},:,:)
    tpv = t_file->t({date},{lev_pv},:,:)

    pv_all = PotVortIsobaric((/20000,25000,30000/),upv,vpv,tpv,tpv&lat,1,0)
    pv = dim_avg_n_Wrap(pv_all,0)
    pv = pv*(10^6)
    copy_VarCoords(upv(0,:,:),pv)

    u = u_file->u({date},{250},:,:)
    v = v_file->v({date},{250},:,:)
    wind = ((u^2)+(v^2))^(0.5)
    copy_VarCoords(u,wind)

    w_orig := w_file->w({date},{600:400},:,:)
    w_temp := w_orig/100          ; convert from Pa/s to hPa/s
    w := dim_avg_n_Wrap(w_temp,0) ; Averages over the specified layer
    copy_VarCoords(u,w)

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
gsn_define_colormap(wks,"alicia_pwat")
;**********************************************

res = True

res@gsnFrame = False                                    ;DO NOT PLOT FIRST MAP
res@gsnDraw = False

res@mpDataBaseVersion = "MediumRes"                     ; Country Outlines
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "National"
res@mpFillOn = False
res@mpOutlineOn = True

res@cnFillOn = True                                     ; Preciptable water information 			
res@cnFillMode = "AreaFill"
res@cnRasterSmoothingOn = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnLinesOn = False                               
res@cnLevelSpacingF = 5                             
res@cnMinLevelValF = 20
res@cnMaxLevelValF = 70                          
res@cnFillColors = (/-1,116,103,84,75,129,139,144,149,157,185,192/)

;**********************************************

res@lbAutoManage = False                  	; Label Bar Information
res@pmLabelBarOrthogonalPosF = .02
res@pmLabelBarParallelPosF = .19  
res@pmLabelBarHeightF = 0.08
res@pmLabelBarWidthF = 0.36
res@lbLabelAutoStride = True
res@lbLabelFontHeightF = 0.0124
res@lbBoxLinesOn  = True                   
res@lbTopMarginF = .3

;**********************************************

res@tiMainString = " "						; Map label information
res@tiMainJust = "CenterCenter"
res@tiMainPosition = "Center" 
res@tiMainOffsetYF = .02
res@tiXAxisFontHeightF = .015

res@gsnLeftString = "250-hPa jet; 600~F15~-~F21~400-hPa ascent; 300~F15~-~F21~200-hPa PV/irr. wind; PW"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 4.5		; Map information
res@mpNationalLineThicknessF = 4.5
res@mpUSStateLineThicknessF = 4.5
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@pmTickMarkDisplayMode = "Always"

res@mpGeophysicalLineColor = "grey31"
res@mpNationalLineColor = "grey31"
res@mpUSStateLineColor = "grey31"

res@mpDefaultFillColor = "green2"

; *********************************************

res@mpGridAndLimbOn = True					; Grid point and tickmark information
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

jres = True							; Shaded wind speed information

jres@gsnFrame = False                       
jres@gsnDraw = False

jres@cnFillOn = True                            
jres@cnFillMode = "AreaFill"
jres@cnRasterSmoothingOn = True
jres@cnLevelSelectionMode = "ExplicitLevels"
jres@cnLinesOn = False                               
jres@cnLevelSpacingF = 1                           
jres@cnLevels = (/30,40,50,60,70,80/)
jres@cnFillColors = (/-1,213,212,211,231,232,233/)
jres@lbLabelBarOn = True
jres@lbAutoManage = False
jres@cnLineLabelsOn = False
jres@cnInfoLabelOn = False 


if domain .eq. "northamer" then					; Define wind speed color bar (northamer)
jres@pmLabelBarOrthogonalPosF = .0213
jres@pmLabelBarParallelPosF = .81    
jres@pmLabelBarHeightF = 0.11
jres@pmLabelBarWidthF = 0.36
jres@lbLabelAutoStride = True
jres@lbLabelFontHeightF = 0.0124
jres@lbBoxLinesOn  = True                            
jres@lbTopMarginF = .3
jres@gsnLeftString = " "
jres@gsnRightString = " "
end if

if domain .eq. "conus"							; Define wind speed color bar (CONUS)
jres@pmLabelBarOrthogonalPosF = .0180
jres@pmLabelBarParallelPosF = .81    
jres@pmLabelBarHeightF = 0.118
jres@pmLabelBarWidthF = 0.36
jres@lbLabelAutoStride = True
jres@lbLabelFontHeightF = 0.0124
jres@lbBoxLinesOn  = True                            	 
jres@lbTopMarginF = .3
jres@gsnLeftString = " "
jres@gsnRightString = " "
end if

jres@gsnAddCyclic =  True

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
pvres@cnLineLabelFontColor = "grey10"
pvres@cnLineColor = "grey10"

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
res@cnFillDrawOrder = "PreDraw"						; Draw order information
jres@cnFillDrawOrder = "PreDraw"
pvres@cnLineDrawOrder ="Draw"
omres@cnLineDrawOrder ="Draw" 
vres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpGridAndLimbDrawOrder = "PostDraw"

plot = gsn_csm_contour_map(wks, pwat, res)          ; Plot filled contour map
jet = gsn_csm_contour(wks, wind, jres)
contours = gsn_csm_contour(wks, pv, pvres)
other = gsn_csm_contour(wks, w, omres)
vect = gsn_csm_vector(wks, irro(0,:,:), irro(1,:,:), vres)
plot = add_map_tickmarks(wks,plot,tmres)

;********************************************

overlay(plot,jet)
overlay(plot,contours)
overlay(plot,other)
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
