; ===========================================
; Plots 850-hPa theta-e, wind, etc. from NCEP CFSR
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
load "~/ncl/lib/fgen_globe.ncl"

; ===========================================

begin

print(" ")
print("**************")
print(" START SCRIPT ")
print("**************")

;*******************************************************
; Define what variables, domains, and directories
;*******************************************************

name = "thetae_850hPa"   ; image name

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
    
	lev = 850
    
    g = g_file->g({date},{lev},:,:)
    g = g/10

    u = u_file->u({date},{lev},:,:)
    v = v_file->v({date},{lev},:,:)
    ubarb = u * 1.94384449
    copy_VarCoords(u,ubarb) 
    vbarb = v * 1.94384449
    copy_VarCoords(v,vbarb)  
    
; ******************************************************    
    
; Calculating potential temperature
fil1 = addfile("/cfsr/data/"+yyyy+"/t."+yyyy+".0p5.anl.nc","r")
fil2 = addfile("/cfsr/data/"+yyyy+"/q."+yyyy+".0p5.anl.nc","r")
t0 = fil1->t({date},{1000:200},:,:)
q0 = fil2->q({date},{1000:200},:,:)

t = t0
t@_FillValue = -999.

q = q0
q@_FillValue = -999.

w = q   
w = q/(1.-q)  ; convert specific humidity to mixing ratio
w@_FillValue = -999.

c_pd = 1004. ; specific heat at constant pressure for air [approx 1004 J/(kg-K)]
R = 287. ; specific gas constant for air [J/(kg-K)]
Rcpd = R/c_pd
p0 = 100000.  ; Pa
L_v = 2400000.  ; latent heat of vaporization at 25C
	
pc = conform(t,t&lev*100.,0) ; array containing pressure
copy_VarCoords(t,pc)
pc@_FillValue = -999.

;print("computing RH")
rh = t	
rh = relhum(t,w,pc)  ; compute relative humidity
rh@_FillValue = -999.

;print("computing dewpoint")
td = t
td = dewtemp_trh(t,rh) ; compute dewpoint temperature
td@_FillValue = -999.
td = where(td.eq.0,.01,td)  ; make sure no dewpoint values = 0; for some reason the dewpoint calculation in ncl yields some 0 K dewpoints 


;Calculate theta-e from Bolton (1980) that yields very accurate theta-e values

tL = (1./((1./(td-56.))+(log(t/td)/800.)))+56.
term2 = (3.376/tL-.00254)*w*1000.*(1.+.81*w)
et = (.2854*(1.-(0.28*w)))
term1 =  t*(p0/pc)^et

theta_e = term1*exp(term2)
copy_VarCoords(t,theta_e)

;==============================================
; Create Shaded Map
;==============================================

wks_type = "png"

wks_type@wkWidth = 1600
wks_type@wkHeight = 1600

wks = gsn_open_wks(wks_type,img_folder+name+"_"+number)
gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")
;gsn_reverse_colormap(wks)
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
res@cnLevelSelectionMode = "ManualLevels"
res@cnLinesOn = False                               
res@cnLevelSpacingF = 6											; Shading theta-e information
res@cnMinLevelValF = 250
res@cnMaxLevelValF = 370
res@cnFillColors = (/0,21,31,35,42,47,55,62,71,80,90,110,120,130,138,154,162,170,186,192,196,199/)

;**********************************************

res@lbAutoManage = False					; Label bar information
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

res@gsnLeftString = lev+"-hPa theta-e, heights, and wind"
res@gsnLeftStringOrthogonalPosF = .012
res@gsnLeftStringFontHeightF = .008
res@gsnRightString = "(CFSR) "+ut_string(date,"")
res@gsnRightStringOrthogonalPosF = .0075
res@gsnRightStringFontHeightF = .013

; *********************************************

res@mpGeophysicalLineThicknessF = 3.5		; Map information
res@mpNationalLineThicknessF = 3.5
res@mpUSStateLineThicknessF = 3.5
res@mpOutlineSpecifiers = (/"United States"/)+" : States"

res@mpGeophysicalLineColor = "navyblue"
res@mpNationalLineColor = "navyblue"
res@mpUSStateLineColor = "navyblue"

; *********************************************

res@mpGridAndLimbOn = True					; Grid lines and tick mark information
res@mpGridLineDashPattern = 5
res@mpGridLatSpacingF = 20.
res@mpGridLonSpacingF = 20.
res@mpGridLineThicknessF = 3
res@mpGridLineColor = "grey12"

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

cnres = True								; Geopotential height information

cnres@gsnFrame = False                                        
cnres@gsnDraw = False

cnres@gsnLeftString = " "
cnres@gsnRightString = " "

cnres@cnFillMode = "AreaFill"
cnres@cnInfoLabelOn = False                   
cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLabelMasking = True
cnres@cnLineThicknessF = 4.5                    
cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 3
cnres@cnMinLevelValF = 90
cnres@cnMaxLevelValF = 420  
cnres@cnLineLabelAngleF = 0.0
cnres@cnLineLabelInterval = 2
cnres@cnLineLabelDensityF = 1.1

cnres@cnLineLabelBackgroundColor = -1
cnres@cnLineLabelFontHeightF = 0.008

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
    bres@gsnAddCyclic = True

; *************************************************************
res@mpFillDrawOrder = "PreDraw"
res@cnFillDrawOrder = "PreDraw"
cnres@cnLineDrawOrder ="Draw"
bres@vcVectorDrawOrder = "Draw"
res@mpOutlineDrawOrder = "Draw"
res@mpGridAndLimbDrawOrder = "PostDraw"

plot = gsn_csm_contour_map(wks, theta_e({850},:,:), res)          	; Plot filled thetae map
contours = gsn_csm_contour(wks, g, cnres)							; Plot geopotential height
barbs = gsn_csm_vector(wks, ubarb, vbarb, bres)						; Plot wind barbs
plot = add_map_tickmarks(wks,plot,tmres)							; Plot tickmarks

;********************************************

overlay(plot,contours)		; Define overlay order
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
