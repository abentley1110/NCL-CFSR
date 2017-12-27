; ===========================================
; Creates files of q-vectors and forcing for ascent from NCEP CFSR
; Last updated: 26 December 2017
; Created by: Alicia M. Bentley
; Email: abentley1110@gmail.com
; ===========================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "~/ncl/lib/ut_string.ncl"
load "~/ncl/lib/latlon.ncl"
load "~/ncl/lib/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "~/ncl/lib/q_vector.ncl"
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

folder = (/"northamer"/)	; domain(s) to plot 
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

; ******************************************************
; Load files
; ******************************************************

ntimes = toint(((edate - bdate)/6) + 1)   	; Find the number of times from the start time to the end time.
                                        	; E.g., if startTime = 700024 and endTime = 700048, and the data is 
                                            ; available every 6 hours, we'd expect ntimes = 5. 
                                            ; We see that (700048 - 700024)/6 + 1 = 4+1 = 5

timestepfactor = 4 					;in 6 hr intervals

times = new( (/ntimes/), "float")  	; Create a new 1-dimensional float array of size nTimes
times@units = timeUnits            	; Give the new times array the same units as timeUnits (in this case: 
                                   	; "hours since 1800-01-01 00:00:00")

do n=0,ntimes-1                           	; This do loop fills the times array with the times from the startTime to the endTime
   times(n) = tofloat(bdate + (6*n))  		; in 6 hour increments, in units of "hours since 1800-01-01 00:00:00"
end do

;-------------------------------------------------

if domain .eq. "northamer" then      ; define domain (must not be global)
;domain of interest
    lW = -179.
    lE = -10.
    lN = 80.
    lS = -10.
else if domain .eq. "conus" then
    lW = -150.
    lE = -35.
    lN = 65.
    lS = 5.
end if
end if

;----------------------------------------------------

u_file = addfile("/cfsr/data/"+byyyy+"/u."+byyyy+".0p5.anl.nc","r")		; download data
v_file = addfile("/cfsr/data/"+byyyy+"/v."+byyyy+".0p5.anl.nc","r")
g_file = addfile("/cfsr/data/"+byyyy+"/g."+byyyy+".0p5.anl.nc","r")
t_file = addfile("/cfsr/data/"+byyyy+"/t."+byyyy+".0p5.anl.nc","r")

lev = (/900,850,800/) 		; define layer for q-vector calculation

fileName = "/cfsr/data/"+byyyy+"/t."+byyyy+".0p5.anl.nc"
inFile_t = addfile( fileName, "r" )             
t =  inFile_t->t({times},{lev},{lS-5:lN+5},{lW-5:lE+5})      	; Get data a little bigger than what you want
																; so that function can calculate q-vector at edges
fileName = "/cfsr/data/"+byyyy+"/g."+byyyy+".0p5.anl.nc"
inFile_g = addfile( fileName, "r" )             
z =  inFile_g->g({times},{lev},{lS-5:lN+5},{lW-5:lE+5})   

q1 := new(dimsizes(z),float) ;predefine arrays to hold Q-vector information
q2 := new(dimsizes(z),float)
qn1 := new(dimsizes(z),float)
qn2 := new(dimsizes(z),float)
qs1 := new(dimsizes(z),float)
qs2 := new(dimsizes(z),float)

dv := new(dimsizes(z),float)
dvn := new(dimsizes(z),float)
dvs := new(dimsizes(z),float)

np = 50       ;number of smoothing passes

Q_vector_natural(t,z,np,q1,q2,qn1,qn2,qs1,qs2,dv,dvn,dvs) ;call Q-vector in natural coordinates procedure

;q1_avg = dim_avg_n_Wrap(q1,1)
q1 = q1*(10.^10.)					; give q1, q2, and dv correct units
;printMinMax(q1, True)

;q2_avg = dim_avg_n_Wrap(q2,1)
q2 = q2*(10.^10.)
;printMinMax(q2, True)

;dv_avg = dim_avg_n_Wrap(dv,1)
dv = dv*(10.^17.)
;printMinMax(dv, True)

;-----------------------------------------------
; Create dummy global file and insert q1, q2, dv
;-----------------------------------------------

fileName = "/cfsr/data/"+byyyy+"/g."+byyyy+".0p5.anl.nc"
inFile_g = addfile( fileName, "r" )             
g = g_file->g({times},{lev},:,:)

q1_final = g
q1_final(:,:,:,:) = 0.0
q1_final(:,:,{lS:lN},{lW:lE}) = q1(:,:,{lS:lN},{lW:lE})
system ("'rm' q1_"+foldername+".nc")
ncdf_q1 = addfile("q1_"+foldername+".nc" , "c")  
ncdf_q1->q1 = q1_final

q2_final = g
q2_final(:,:,:,:) = 0.0
q2_final(:,:,{lS:lN},{lW:lE}) = q2(:,:,{lS:lN},{lW:lE})
system ("'rm' q2_"+foldername+".nc")
ncdf_q2 = addfile("q2_"+foldername+".nc" , "c")  
ncdf_q2->q2 = q2_final

dv_final = g
dv_final(:,:,:,:) = 0.0
dv_final(:,:,{lS:lN},{lW:lE}) = dv(:,:,{lS:lN},{lW:lE})
system ("'rm' dv_"+foldername+".nc")
ncdf_dv = addfile("dv_"+foldername+".nc" , "c")  
ncdf_dv->dv = dv_final

print("Q-vector data saved to files!")

;-------------------------------------------------

end do 		; end domain

print(" ")
print("**************")
print(" END SCRIPT ")
print("**************")
print(" ")

end