
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt
from math import inf,pi,sqrt,tan
from scipy.signal import medfilt


readpath=r"D:\SEDA\09-02-2022\process1\09-02-2022_16-50-13\Step 7 Indentation\ForceDisplacement.txt"
with open (readpath,'r',encoding='utf-8') as f:
    lines=f.readlines()



def area_func(d,R_0=0.7,alpha=1.22):
    return pi*(R_0**2+2*R_0*np.tan(alpha)*d+tan(alpha)**2*d**2)

    
def fit_func(h,a,c,h_0,m):
    return  a*(h-h_0)**m+c
    

def der_func(h,a,c,h_0,m):
    return m*a*(h-h_0)**(m-1)





F=[]
disp=[]
for k in range(0,len(lines)):
    try:
        disp.append(float(lines[k].split('\t')[1]))
        F.append(float(lines[k].split('\t')[2]))

    except:
        print(lines[k])
F=medfilt(F,11)
F1=np.array(F)
disp1=np.array(disp)

f0=5.
i=np.where(F1<min(F1)+f0)[0]
h_0=disp1[i][0]

fmin=min(F)
for k in range(len(disp)):
    if (np.isfinite(disp[k])!=True):
        print(disp[k],k)
    F[k]=F[k]-fmin
##disp=np.flipud(disp)
##F=np.flipud(F)


p0=[1,1,h_0,1.5]
bounds=([0.,-200,h_0-5,1],[1000,200,h_0+5,2])
para,cov=curve_fit(fit_func,disp[0:i[0]],F[0:i[0]],check_finite=False,p0=p0,maxfev = 1000000)

hmax=max(disp)
S=der_func(hmax,*para)*10**3#contact stifness


E=sqrt(pi)*S/(2*sqrt(area_func(hmax)))#reduced young modulus
print('S=',S,'N/m','E=',E,'MPa')
plt.plot(disp[0:i[0]],F[0:i[0]])
plt.show()
#para,cov=curve_fit(fit_func,disp,F,p0=[1,h_0,1.5],bounds=([0.01,h_0-1,1],[1000,h_0+1,2]),check_finite=False)
## #
#print('a=',para[0],',h_0=',para[1],',m=',para[2])
labell='%s*(h-%s)**%s+%s'%(str(para[0]),str(para[2]),str(para[3]),str(para[1]))
plt.plot(disp,F,disp,fit_func(disp,*para),label=labell)
plt.legend()
plt.show()
