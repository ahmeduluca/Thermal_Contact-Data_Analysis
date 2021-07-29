# This code calculates load-displacement-area from
# analytical results of indentation for various geometries.

## After above calculation, determines A/F functions to optimize the most
## area for in a range of load.

### At last, for thermal side of the problem, by borrowing fin problem solution
### for 3 types of pinned fins- cylindrical, conical and parabolic cases.
### fin effectiveness eta*area functions wrt max heat transfer possible
### through fin. optimizing that to %99.9 of qmax for extracting length or
### diameter/length. moreover, other variables' effect will be shown via graphs
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi,sqrt,tan,tanh
import scipy
from scipy.special import iv
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

## Indenter Area Functions:

## ideal theta values:
vick_th=68
berk_th=65.27
knoop_th1=86.25
knoop_th2=65
cube_th=35.26

def conic_area(h,alpha):
    alpha=math.radians(alpha)
    return pi*(h**2)*tan(alpha)**2

def spherical_area(h,r):
    return pi*(2*r*h-h**2)

def vickers_area(h, theta): ## 4 sided pyramids
    theta=math.radians(theta)
    return 4*(h**2)*tan(theta)**2

def berko_area(h, theta): ## 3 sided pyramids
    theta=math.radians(theta)
    return 3*sqrt(3)*(h**2)*tan(theta)**2

def knoop_area(h, theta1, theta2):
    theta1=math.radians(theta1)
    theta2=math.radians(theta2)
    return 2*(h**2)*tan(theta1)*tan(theta2)

def cube_area(h): ## a special case of 3 sided pyramids
    return 3*sqrt(3)*(h**2)*tan(35.26)**2

def flat_area(a):
    return pi*a**2

## Indenter Load Functions:

## Equivalent cone angles for pyramids for using in conic_load function:
vickers=70.3
berkovich=70.296
knoop=77.64
cube=42.778

def conic_load(E, h, alpha):
    alpha=math.radians(alpha)
    return (2*E*tan(alpha)*h**2)/pi

def spherical_load(E, h, r):
    return (4/3)*E*sqrt(r)*h**1.5

def flat_load(E, h, a):
    return 2*a*E*h

## 1 important thing is 'h' of load function is indentation depth h_max
## but, h of area functions is contact depth h_c

def h_c(h,eps,gam):
    return h-eps*gam
## gam is a gamma function that Fmax/S ; S is contact stiffness obtained through
## retraction from Fmax slope of F-h at Fmax

## E is indentation modulus Young's(Y) modulus div by poisson's ratio:

def Ei(Y, poi):
    return Y/(1-poi**2)

## if indenter is not very high stiffness then a reduced modulus:

def Er(Y1, Y2, poi1, poi2):
    return 1/(((1-poi1**2)/Y1)+((1-poi2**2)/Y2))

## Inverse Hardness optimization: A/F function extraction;

## conic angle effect: 
depth=np.linspace(0.000000001,0.000010,num=10000)
area=np.zeros((89,10000))
force=np.zeros((89,10000))
a_f=np.zeros(89)
for j in range(0,89):
    for i in range(len(depth)):
        area[j][i]=(conic_area(depth[i],j+1))
        force[j][i]=(conic_load(Er(130e9,1220e9,0.34,0.2),depth[i],j+1))
    if(force[j][5000]==0):
        a_f[j]=0
    else:
        a_f[j]=(area[j][5000]/force[j][5000])
##        if(a_f[j]>=1.01e-10):
##           print(j)
##           break
figure, axis_1=plt.subplots()
axis_1.plot(depth, np.transpose(area),'.')
axis_1.yaxis.set_label_text('Area')
axis_1.xaxis.set_label_text('Depth')
axis_2=axis_1.twinx()
axis_2.plot(depth, np.transpose(force),':')
axis_2.yaxis.set_label_text('Force')
##plt.plot(force,'o-', area, label="f vs A")
lines_1, labels_1 = axis_1.get_legend_handles_labels()
lines_2, labels_2 = axis_2.get_legend_handles_labels()
lines = lines_1 + lines_2
labels = labels_1 + labels_2
axis_1.legend()
plt.show()
plt.axes(yscale='linear')
plt.ylabel("Area/Force (m2/N)")
plt.xlabel("Conic Half Angle")
plt.title("Diamond-Copper Area/Force")
plt.plot(a_f,'o-',label="Conic Angles vs Area/Force")
plt.legend()
plt.xticks(np.arange(0,91,10))
plt.yticks(np.linspace(0,np.amax(a_f),10))
plt.minorticks_on()
plt.grid(True)
plt.show()

## Well-defined shapes comparison:
vick_ar=np.zeros(10000)
berk_ar=np.zeros(10000)
cube_ar=np.zeros(10000)
knoop_ar=np.zeros(10000)
spher_ar=np.zeros(10000)
flat_ar=np.zeros(10000)
vick_f=np.zeros(10000)
berk_f=np.zeros(10000)
cube_f=np.zeros(10000)
knoop_f=np.zeros(10000)
spher_f=np.zeros(10000)
flat_f=np.zeros(10000)

for i in range(len(depth)):
    vick_ar[i]=vickers_area(depth[i],vick_th)
    berk_ar[i]=berko_area(depth[i],berk_th)
    knoop_ar[i]=knoop_area(depth[i],knoop_th1,knoop_th2)
    cube_ar[i]=cube_area(depth[i])
    spher_ar[i]=spherical_area(depth[i],0.00001)
    flat_ar[i]=flat_area(0.000001)
    vick_f[i]=conic_load(Er(130e9,1220e9,0.34,0.2),depth[i],vickers)
    berk_f[i]=conic_load(Er(130e9,1220e9,0.34,0.2),depth[i],berkovich)
    knoop_f[i]=conic_load(Er(130e9,1220e9,0.34,0.2),depth[i],knoop)
    cube_f[i]=conic_load(Er(130e9,1220e9,0.34,0.2),depth[i],cube)
    spher_f[i]=spherical_load(Er(130e9,1220e9,0.34,0.2),depth[i],0.00001)
    flat_f[i]=flat_load(Er(130e9,1220e9,0.34,0.2),depth[i],0.000001)
vickers_af=vick_ar/vick_f
berko_af=berk_ar/berk_f
knoop_af=knoop_ar/knoop_f
cube_af=cube_ar/cube_f
spher_af=spher_ar/spher_f    
flat_af=flat_ar/flat_f
plt.plot(vick_f,vick_ar,'.-',label="Vickers")
plt.plot(berk_f,berk_ar,',--', label="Berkovich")
plt.plot(knoop_f,knoop_ar, '.-',label="Knoop")
plt.plot(cube_f,cube_ar, '.-',label="Cube Corner")
plt.plot(spher_f,spher_ar, '.-',label="Spherical R=10um")
plt.plot(flat_f,flat_ar,'.-', label="Flat Punch a=1um")
plt.xlabel("Force (N)")
plt.ylabel("Area (m2)")
plt.title("Force vs Area up to 10um Indentation of Copper")
plt.grid(True)
plt.legend()
plt.show()

plt.plot(depth,vickers_af, '.-',label="Vickers")
plt.plot(depth,berko_af, ',--',label="Berkovich")
plt.plot(depth,knoop_af,'.-',label="Knoop")
plt.plot(depth,cube_af,'.-',label="Cube Corner")
plt.plot(depth,spher_af,'.-',label="Spherical R=10um")
plt.plot(depth,flat_af,'.-',label="Flat Punch a=1um")
plt.legend()
plt.grid(True)
plt.xlabel("Depth")
plt.ylabel("Area/Force")
plt.title("Depth vs Area/Force Indentation of Copper")
plt.show()

## Pinned Fin effectiveness functions:

def cylinder_etha(h, k, d, l):
    m=2*sqrt(h/(k*d))
    return tanh((l+d/4)*m)/((l+d/4)*m)
def cylinder_surf(d,l):
    return pi*d*(l+d/4)
def cylinder_vol(d,l):
    return pi*d**2*l/4

def conic_etha(h, k, d, l):
    m=2*sqrt(h/(k*d))
    return (2/(m*l))*(iv(2,2*m*l)/iv(1,2*m*l))
def conic_surf(d, l):
    return pi*d*sqrt((l**2+(d/2)**2))/2
def conic_vol(d, l):
    return pi*d**2/(l*12)

def parabol_etha(h, k, d, l):
    m=2*sqrt(h/(k*d))
    return 2/(sqrt((4/9)*(m*l)**2+1)+1)
def parabol_area(d, l):
    c3=1+2*(d/l)**2
    c4=sqrt(1+(d/l)**2)
    return (pi*l**3/(8*d))*(c3*c4-(l/(2*d))*np.log((2*d*c4/l)+c3))
def parabol_vol(d, l):
    return pi*d**2*l/20

##efficiency functions wrt length of fin:

def Spreading(contact_dia,k1,k2):
    return contact_dia*(2*k1*k2/(k1+k2))
k=2200
k_sample=380
fin_len=np.linspace(0.000000001,0.00001,10000)
fin_dia=0.000002
#h_sp=Spreading(fin_dia,k,k_sample)
h_sp=1e6
eff_cylinder=np.zeros(len(fin_len))
eff_conic=np.zeros(len(fin_len))
eff_parabol=np.zeros(len(fin_len))
angle=np.zeros(len(fin_len))
tot_eff_cal=np.linspace(0,90,91)
tot_eff_ang=np.zeros(91)

for i in range(len(fin_len)):
    eff_cylinder[i]=cylinder_etha(h_sp,k,fin_dia,fin_len[i])
    eff_conic[i]=conic_etha(h_sp,k,fin_dia,fin_len[i])
    eff_parabol[i]=parabol_etha(h_sp,k,fin_dia,fin_len[i])
    angle[i]=np.degrees(np.arctan(fin_dia/(2*fin_len[i])))
    if(int(angle[i]) in tot_eff_cal):
        tot_eff_cal[int(angle[i])]=-1
        tot_eff_ang[int(angle[i])]=eff_conic[i]*a_f[int(angle[i])-1]
        
plt.plot(fin_len,eff_cylinder/max(eff_cylinder),'.',label="Cylinder Efficiency")
plt.plot(fin_len,eff_conic/max(eff_conic),'.',label="Conic Efficiency")
plt.plot(fin_len,eff_parabol/max(eff_parabol),'.',label="Parabola Efficiency")
plt.xticks(np.linspace(0,np.amax(fin_len),10))
plt.yticks(np.linspace(np.amin(eff_cylinder),np.amax(eff_cylinder),10))
plt.minorticks_on()
plt.grid(True)
plt.title("Surface Extension Efficiency wrt Length")
plt.annotate(" h=2*a*k_harmonic(diamond&copper) \n k=k_diamond \n diameter=2um",
         xy=(0,0.5),xycoords='axes fraction')
plt.xlabel("Length (m)")
plt.ylabel("Efficiency")
plt.legend()
plt.show()
plt.plot(angle,eff_cylinder,'.',label="Cylinder Efficiency")
plt.plot(angle,eff_conic,'.',label="Conic Efficiency")
plt.plot(angle,eff_parabol,'.',label="Parabola Efficiency")
plt.legend()
plt.xticks(np.linspace(0,np.amax(angle),10))
plt.yticks(np.linspace(np.amin(eff_cylinder),np.amax(eff_cylinder),10))
plt.minorticks_on()
plt.grid(True)
plt.title("Surface Extension Efficiency wrt Cone Angle")
plt.annotate(" h=2*a*k_harmonic(diamond&copper) \n k=k_diamond \n diameter=2um\n angle=arctan(diameter/(2*length))",
         xy=(0.3,0.5),xycoords='axes fraction')
plt.xlabel("Conic Half Angle")
plt.ylabel("Efficiency")
plt.show()
def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]
tot_eff_ang=tot_eff_ang/max(tot_eff_ang)
tot_eff_ang=zero_to_nan(tot_eff_ang)
plt.plot(tot_eff_ang, 'o-',label="Conic Total Efficiency")
plt.legend()
plt.title("Normalized Efficency for Conic Fin Indent")
plt.xlabel("Conic Half Angle")
plt.ylabel("Efficiency")
#plt.xticks(np.linspace(0,np.amax(angle),10))
#plt.yticks(np.linspace(np.amin(tot_eff_ang),np.amax(tot_eff_ang),10))
plt.minorticks_on()
plt.grid(True)
plt.annotate(" Conic Fin Efficiency * Area/Force = Q/F",
         xy=(0,0.5),xycoords='axes fraction')
plt.show()
## Qmax calculations for a certain diameter || diameter/length || length

def Q_unit(h,A,etha):
    return h*A*etha


##  calculation:


print(h_sp)
q_calc=np.zeros(10000)
tot_eff_cal=np.linspace(0,90,91)
tot_eff_ang=np.zeros(91)
tot_eff_force=np.zeros(91)
eff_for=np.zeros(91)
eff_angle=np.zeros(91)
print(tot_eff_cal[-1])
for i in range(0,10000):
    q_calc[i]=h_sp*cylinder_etha(h_sp,k,fin_dia,fin_len[i])*cylinder_surf(fin_dia,fin_len[i])
    angle[i]=np.degrees(np.arctan(fin_dia/(2*fin_len[i])))
#plt.axes(yscale="log")
q_max=max(q_calc)
#plt.plot(angle,q_calc/q_max, 'o-',label="Cylinder Pin Fin")
plt.plot(fin_len,q_calc/q_max, '.',label="Cylinder Pin Fin")

for i in range(0,10000):
    q_calc[i]=h_sp*parabol_etha(h_sp,k,fin_dia,fin_len[i])*parabol_area(fin_dia,fin_len[i])
q_max=max(q_calc)
#plt.plot(angle,q_calc/q_max, 'o-',label="Parabola Pin Fin")
plt.plot(fin_len,q_calc/q_max, '.',label="Parabola Pin Fin")

for i in range(0,10000):
    q_calc[i]=h_sp*conic_etha(h_sp,k,fin_dia,fin_len[i])*conic_surf(fin_dia,fin_len[i])
##    if(int(angle[i]) in tot_eff_cal):
##        tot_eff_cal[int(angle[i])]=-1
##        try:
##            tot_eff_force[int(angle[i])]=force[int(angle[i])][i]
##            eff_for[int(angle[i])]=q_calc[i]
##            eff_angle[int(angle[i])]=angle[i]
##        except:
##            print(angle[i])
##            continue
q_max=max(q_calc)
#plt.plot(angle,q_calc/q_max, 'o-',label="Conic Pin Fin")
plt.plot(fin_len,q_calc/q_max, '.',label="Conic Pin Fin")
#plt.plot(a_f,'o-',label="Conic Area/Force")
plt.xlabel("Length (m)")
#plt.xlabel("Conic Half Angle")
#plt.xticks(np.linspace(0,np.amax(fin_len),10))
#plt.yticks(np.linspace(0,1,10))
plt.minorticks_on()
plt.grid(True)
plt.ylabel("Q/Q_max")
plt.title("Heat Transfer Rate wrt Length")
#plt.annotate(" Q=h*A*etha*dT \n h=2*a*k_harmonic(diamond&copper) \n k=k_diamond \n constant diameter=2um\n angle=arctan(diameter/(2*length))",
        # xy=(0.3,0.5),xycoords='axes fraction')
plt.annotate(" Q=h*A*etha*dT \n h=2*a*k_harmonic(diamond&copper) \n k=k_diamond \n constant diameter=2um",
         xy=(0.3,0.6),xycoords='axes fraction')
##plt.plot(angle,fin_len)
plt.legend()
plt.show()

##plt.plot(tot_eff_ang, 'o-',label="Conic Total Efficiency")
##plt.legend()
##plt.show()
##plt.plot(tot_eff_force,eff_for/q_max,'o-',label="Conic Total Force Efficiency")
##plt.legend()
##plt.show()
##plt.plot(eff_angle,eff_for/q_max, 'o-',label="anglevQ")
##plt.legend()
##plt.show()
##plt.plot(eff_angle,tot_eff_force, 'o-',label="anglevForce")
##plt.legend()
##plt.show()
##plt.plot(fin_len,q_calc/q_max, 'o-')
##plt.show()

