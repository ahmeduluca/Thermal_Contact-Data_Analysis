from io import open
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import scipy.interpolate as interp
from scipy.signal import savgol_filter
from scipy.fftpack import rfft, irfft, rfftfreq
from nptdms import TdmsFile
import numpy as np
import glob
import os
import importlib
from scipy.signal import medfilt
import math
from scipy.special import gamma
from natsort import natsorted

"""
This code is for reading experimental data of Indenter/RC Thermal setup obtained by
NI-DAQ written/converted txt data files are being searched as txt. 
"""
##Special Characters##
degsin=u"\N{DEGREE SIGN}"
degcel=degsin+"C"
deltasig=u"\N{GREEK CAPITAL LETTER DELTA}"
mikron=u"\N{GREEK SMALL LETTER MU}"+"m"
###

## Main file of Expriment Data of text files -that converted from NI TDMS-
#main=r"D:\SEDA\18-02-2022\process3" #WindowsFilePath
main="/Volumes/AhmedUluca/Indenter_Data/main/Berkovich/Au/Indent Sweep Different Locations"
 #MacFilePath
###

tip_chname="Temperature9"
sample_chname='Temperature12'
gage_chname="Voltage0"
pos_chname="Actuator_Voltage"

## Calibration constants of Experiment ##
lc_mass_cal= 4.75 #1000
grave= 9.81 #1
bit_to_load=lc_mass_cal*grave
if("Cu" in main):
    tip_chname="Temperature2"
    sample_chname='Temperature3'
loadcell_stiffness=536.5 ##mN/um
gage2um=10 ##um/V
load_threshold=2.5 ##mN loadcell sensitivity for touch
terr=0.08#error for temperature
disp_tresh=.5
sg_tol=0.2
gradSayı=1000
gradStep=0.007

###
exp_per=4000 ##millisec
daq_frq=0.2 ##kHz
osc_per=exp_per*daq_frq
load_files=[]
position_files=[]
temp_files=[]
sg_files=[]
lvsdis_file=[[]]

date=[[]]
n=0

## Fit functions for experiment data
def rise(t,N1,tau1,c):
    """The model function for the rise of the curve"""
    return N1*(1-np.exp(-t/tau1))+c


def decay(t,N2,tau2,c):
    """Model function for the decay part of the curve"""
    return N2*(np.exp(-t/tau2))+c

def linear(x,m,c):
    return m*x+c

def inverse(x,m,c):
    return m*x**-1+c

def quad(x,a,b,c):
    return a*x**2+b*x+c

def oliver(h,A,f,m):
    return A*(h-f)**m

def gradien(arr,inter,thresh):
    counter=0
    sayac=0
    for i in range(0,len(arr)):
        if(arr[i]+thresh<arr[i+1]):
            if(counter>=0):
                counter+=1
            else:
                sayac+=1
                if(sayac>3):
                    counter=0
                    sayac=0
                    
        elif(arr[i]-thresh>arr[i+1]):
            if(counter<=0):
                counter-=1
            else:
                sayac+=1
                if(sayac>3):
                    counter=0
                    sayac=0
        if(counter>inter and inter>0):
            return i-inter+1
            break
        elif(counter<inter and inter<0):
            return i+inter+1
            break
        elif(i==len(arr)-2):
            return int(-1)
            break

def indexFinder(arr,condition,value,start=0,stop=0):
    if(stop==0):
        stop=len(arr)
    for i in range(start,stop):
        if(condition==0):
            if(arr[i]==value):
                return i
        elif(condition==1):
            if(arr[i]>value):
                return i
        elif(condition==2):
            if(arr[i]<value):
                return i
        elif(condition==3):
            if(arr[i]>=value):
                return i
        elif(condition==4):
            if(arr[i]<=value):
                return i



def area_projected(d,R_0=0.7,alpha=1.16):
    return np.pi*(R_0**2+2*R_0*np.tan(alpha)*d+np.tan(alpha)**2*d**2)

    
def oliverPharr(h,a,c,h_0,m):
    return  a*(h-h_0)**m+c
    

def stiffness(h,a,c,h_0,m):
    return m*a*(h-h_0)**(m-1)

def martens(h,a):
    return a*h**2#a=1/(ms^2)


def area_surface(h,R0,h0=0.13*10**-6,alpha=1.16):
    return np.pi*((h+h0)**2*np.tan(alpha)-R0*h0)/np.cos(alpha)+np.pi*R0**2

def epsilon(m):
    return (1-(gamma(m/(2*m-2))/(gamma((2*m-1)/(2*m-2))))*(1/math.sqrt(np.pi)))*m
###



file_list=[f.path for f in os.scandir(main) if f.is_dir()]
file_list.sort()
expdates=[f.name for f in os.scandir(main) if f.is_dir()]
expdates.sort()
m=0
tcount=0
ch=0
starter=-1
lstart=-1
final=-1
lfinal=-1
tpt0=[]
tpt1=[]
transit=0
ltrans=0
t1=[]
t2=[]
tempo=[]
sg=[]
file_list.sort()
salla=0
load=[]
position=[]
volt=[]
voltaj=[]
temp=[[],[]]
time=[]
min_vol=0
loadcount=0
## Search for Experiment Files that include Step Files ##
for i in file_list:
    print(file_list[m])
    print(expdates[m])
    sub_dirs=[f.path for f in os.scandir(i) if f.is_dir()]
    sub_dirs=natsorted(sub_dirs)
    stepnames=[f.name for f in os.scandir(i) if f.is_dir()]
    stepnames=natsorted(stepnames)
    results_files=['']*len(stepnames)   
    n=0
    taular=[]
    taulad=[]
    lups=np.zeros((len(stepnames)))
    ldws=np.zeros((len(stepnames)))
    ## Search for txt files of Load - Strain Gage - Temperatures ##    
    for j in sub_dirs:
        load=[]
        volt=[]
        voltaj=[]
        tempTip=[]
        tempSample=[]
        time=[]
        tip=np.zeros(int(osc_per))
        sample=np.zeros(int(osc_per))
        voltMean=np.zeros(int(osc_per))
        txt_files=glob.glob(os.path.join(j+'/*.txt'))
        ch=0
        load_time=[]
        for k in txt_files:
            if("Load" in k):
                load=[]
                load_time=[]
                loadcount=0
                f=open(k,"r",encoding='ascii')
                spl=0
                if(len(f.readline().split())>1):
                   spl=1
                for x in f:
                    load.append(float(x.split()[spl])*bit_to_load/1000)
                    if(spl==1):
                        load_time.append(float(x.split()[0])/1000)
                    else:
                        load_time.append(loadcount*400/1000)
                        loadcount+=1
                f.close()
                if(len(load)>1):
                    print(len(load_time)/len(load))
                else:
                    print("NoLoad")
                    load.append(1.)
                    load_time.append(1.)
                    ##            elif(pos_chname in k):
##                f=open(k,"r")
##                for x in f:
##                    position[n].append(float(x))
            elif(tip_chname in k):
                with open(k,'r',encoding='ascii') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time.append(float(x.split()[0]))
                        tempTip.append(float(x.split()[1]))
                f.close()
                tcount+=1
            elif(sample_chname in k):
                with open(k,'r',encoding='ascii') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time.append(float(x.split()[0]))
                        tempSample.append(float(x.split()[1]))
                f.close()
                tcount+=1
            elif(gage_chname in k):
                with open(k,'r',encoding='ascii') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time.append(float(x.split()[0]))
                        voltaj.append(float(x.split()[1])*gage2um)##um conversion
                tcount+=1
                f.close()
        if(len(load)<1):
            load.append(0)
        if(len(time)<1):
            time.append(0)
        else:
            loadStPt=time[-1]-len(load)*0.05
            print(loadStPt)
#            if(spl!=1):
#                load_time=np.linspace(loadStPt,time[-1],len(load))
        if(len(tempTip)<1):
            tempTip=np.zeros(len(time))
        if(len(tempSample)<1):
            tempSample=np.zeros(len(time))
        if(len(voltaj)<1):
            voltaj=np.zeros(len(time))
            volt=voltaj
        else:
            volt=medfilt(voltaj,21)
##Indentation load displacement graphs & timeconstant extraction trial through indentation stage
        if('Indentation' in stepnames[n] or 'Calibration' in stepnames[n] or 'Approach' in stepnames[n] or 'Retract' in stepnames[n]):
            if(len(voltaj)>101):
                volt=medfilt(volt,13)
            zz=0
            #######   Interpolation ###################################################################

            if(len(load)>1):               
                kernel=int(len(load)/100)
                if(kernel%2==0):
                    kernel+=1
                load=medfilt(load,kernel)
                min_load_t=min(load_time)
                for i3 in range(len(load_time)):
                    load_time[i3]=load_time[i3]-min_load_t
                extdata=(time[-1]-load_time[-1])*daq_frq
                interpol=interp.interp1d(load_time,load,fill_value="extrapolate")
                load=interpol(time)
            #########################
                ### OFFSET ###
            try:
                #load=medfilt(load,31)
                a=np.where(load>min(load)+load_threshold)[0][0]
            except:
                if('Approach' in stepnames[n]):
                    a=len(volt)-1
                else:
                    a=np.where(volt==min(volt))[0][0]
            if((n==0) or (n>0 and ('Indent' in stepnames[n]) and ('Approach' in stepnames[n-1]))):
                min_vol=min(volt)
            for t in range(len(volt)):###Touch point offset substraction
                volt[t]=volt[t]-min_vol

        ##################
                
            tempTip=medfilt(tempTip,21)
            figs,axd=plt.subplots()
            axd.set_title(stepnames[n]+" "+expdates[m].replace('_',' '))
            if(len(load)>1):
                axl=axd.twinx()
                axl.plot(time,load,'.',color='cyan')
                axl.set_ylabel("Load (mN)",color="cyan")
                axl.tick_params(axis='y',labelcolor="cyan")
            axt=axd.twinx()#
            axt.spines['right'].set_position(("outward",75))#
            axd.set_ylabel("Displacement (%s)"%mikron,color='orange')
            axd.set_xlabel("Time (s)",color="black")
            axt.set_ylabel("Tip Temperature (%s)"%degcel,color="red")#
            axd.tick_params(axis='y',labelcolor="orange")
            axt.tick_params(axis='y',labelcolor="red")#
            axd.plot(time,volt,'.',color='orange')
            axt.plot(time,tempTip,'r.')#
            figs.tight_layout()
            plt.savefig(os.path.join(j,"TimevsLoad&Disp.png"),dpi=128)
            plt.close()                            
            figs,axt=plt.subplots()
            #axt.plot(position[zz:int(len(position))],reduceT[zz:int(len(position))],'r.')
            axt.get_yaxis().set_visible(False)
            if(len(load)>1):
                axd=axt.twinx()
                axd.plot(volt[zz:int(len(volt))],load[zz:int(len(volt))],'b.')
            axd.set_title(stepnames[n]+" "+expdates[m].replace('_',' '))
            plt.tight_layout()
            plt.savefig(os.path.join(j,"LoadvsDisp.png"),dpi=128)
            plt.close()
            maxT=max(tempTip)
            minT=min(tempTip)
            #steps of actuator are at each 100ms...for starting time/100 maybe tried
            dTemp=0.2
            number=(maxT-minT)/dTemp
##            lastpt=np.where(tempTip<=minT)
##            lastpt=int(lastpt[0][-1])
            lastpt=a
            tim=[]
            say=0
            tautr=[]
            tautd=[]
            starttime=time[a]
            for t in range(lastpt,len(tempTip)):
                if(tempTip[t]==maxT):
                    starttime=time[t]
                    say=0
                    break
                if(time[t]==starttime+2*(say+1)):
                    ared=tempTip[lastpt:t]
                    tims=time[lastpt:t]
                    if(len(tims)>1):
                        mt=min(tims)
                        tim=[(tims[a]-mt)*1000 for a in range(len(tims))]
                        try:
                            para, covs =curve_fit(rise,tim,ared,maxfev=100000,p0=[dTemp,500,minT])
                            tautr.append(para[1])
                        except:
                            say=say
                        Risefit=rise(np.array(tim), *para)
                        #plt.plot(tim,Risefit)
                        #plt.plot(tim,ared)
                        #plt.show()
                    lastpt=t+1
                    say+=1
            for t in range(lastpt,len(tempTip)):
                if(tempTip[t]==minT):
                            break
                if(time[t]==starttime+2*(say+1)):
                    ared=tempTip[lastpt:t]
                    tims=time[lastpt:t]
                    if(len(tims)>1):
                        mt=min(tims)
                        tim=[(tims[a]-mt)*1000 for a in range(len(tims))]
                        try:
                            para, covs =curve_fit(decay,tim,ared,maxfev=100000,p0=[dTemp,500,maxT])
                            tautd.append(para[1])
                        except:
                            say=say
                        DecayFit=decay(np.array(tim), *para)
                        #plt.plot(tim,DecayFit)
                        #plt.plot(tim,ared)
                        #plt.show()
                    lastpt=t+1
                    say+=1
 #           plt.show()
            #plt.plot(tautr,'.')
            #plt.plot(tautd,'.')
            #plt.show()
### FD Analysis Oliver-Pharr
            if(len(load)>1):
                try:
                    #plt.rcParams['text.usetex']=True
                    load=medfilt(load,11)
                    F1=np.array(load)
                    disp1=np.array(volt)
                    f0=min(load)+load_threshold
                    inop=np.where(F1<min(F1)+load_threshold)[0]#index oliver pharr
                    if(len(inop)>0):
                        h_0=disp1[inop][0]
                    else:
                        h_0=0
                        disp_tresh=(max(volt)-min(volt))/2
                        print('h_0 cannot be found h_0=0')
                    fmin=min(load)
                    for k in range(len(volt)):
                        if (np.isfinite(volt[k])!=True):
                            print(volt[k],k)
                        load[k]=load[k]-fmin
                    ##disp=np.flipud(disp)
                    ##F=np.flipud(F)
                    loaderr=2.5
                    p0=[1,1,round(h_0,3),1.5]
                    bounds=([0.,-100,h_0-disp_tresh,1],[1000,100,h_0+disp_tresh,2])
                    indexxx=gradien(load,gradSayı,gradStep)
                    indexx=len(load)-gradien(np.flipud(load),-gradSayı,gradStep)
                    if(indexxx==-1):#geri çekme
                        indexxx=gradien(load,-gradSayı,gradStep)
                        indexx=len(load)-gradien(np.flipud(load),gradSayı,gradStep)
                    para,cov=curve_fit(oliverPharr,volt[indexxx:indexx],load[indexxx:indexx],bounds=bounds,p0=p0,check_finite=False,maxfev = 10000,sigma=loaderr*np.ones(len(load[indexxx:indexx])),absolute_sigma=True)
                except:
                    print("error on FD")
                    indexxx=gradien(load,gradSayı,gradStep)
                    indexx=len(load)-gradien(np.flipud(load),-gradSayı,gradStep)
                    if(indexxx==-1 or  indexxx>=indexx):
                        indexxx=gradien(load,-gradSayı,gradStep)
 #                       indexx=len(load)-gradien(np.flipud(load),gradSayı,gradStep)
                        if(indexxx==-1 or indexxx==len(load)-1):
                            indexxx=0
                            indexx=len(load)-gradien(np.flipud(load),gradSayı,gradStep)
                    para=[1,1,1,1.1]
                    cov=[[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]
 #                   para,cov=curve_fit(oliverPharr,volt[indexxx:indexx],load[indexxx:indexx],check_finite=False,maxfev = 10000,sigma=loaderr*np.ones(len(load[indexxx:indexx])),absolute_sigma=True)
                #imax=np.where(load<max(load))[0][0]
                hmax=max(volt)
                Fmax=oliverPharr(hmax,*para)
                S=stiffness(hmax,*para)*10**3#contact volt
                c=Fmax-S*hmax/1000
                hr=-c*1000/S
                print('hr=',hr)
                E=math.sqrt(np.pi)*S/(2*math.sqrt(area_projected(hmax)))#reduced young modulus
                eps=epsilon(para[3])
                hc=hmax-eps*(hmax-hr)
                H_ıt=Fmax/area_projected(hc)#GPa
                print('S=',S,'N/m','E=',E,'MPa')
        ############ Uncertanity calculation for S 
                da_=cov[0][0]
                dh0_=cov[2][2]
                dm_=cov[3][3]
                dh_=0.05**2
                a_=para[0]
                h_0=para[2]
                m_=para[3]
                der1=(m_*(hmax-h_0)**(m_-1))**2
                der2=(-a_*m_*(m_-1)*(hmax-h_0)**(m_-2))**2
                der3=(a_*(hmax-h_0)**(m_-1)+(a_*m_*(hmax-h_0)**(m_-1)))**2
                der4=(a_*m_*(m_-1)*(hmax-h_0)**(m_-2))**2
    ##            uncert_S=math.sqrt(der1*da_+der2*dh0_+der3*dm_+der4*dh_)*10**3#in N/m
    ##            print('DeltaS=',uncert_S,'N/m')
        ####### End of uncertanity calculation
                plt.figure()
                labell='$%s(h-(%s))^{%s}+(%s)$'%(str(round(para[0],3)),str(round(para[2],3)),str(round(para[3],3)),str(round(para[1],3)))
                plt.plot(volt,load,label=r'$\varepsilon(m)$=%s'%eps)
                plt.plot(volt,oliverPharr(volt,*para),label=labell)
                #plt.plot(volt[indexxx:indexx],stiffness(volt[indexxx:indexx],*para),label='Stiffness=%s,$E_r$=%s,$h_r$=%s'%(str(np.round(S,3)),str(np.round(E,3)),str(np.round(hr,3))))
                plt.plot([hr,hmax],[0,Fmax],label='$h_r=$%.2f,S=%.2f,$h_c$=%.2f,$H_{IT}$=%.2f'%(hr,S,hc,H_ıt))
                #plt.ylim(0,max(load)*1.5)
                plt.legend()
                plt.xlabel('Displacement(%s)'%mikron)
                plt.ylabel('Load (mN)')
                plt.savefig(os.path.join(j,"Oliver-Pharr.png"),dpi=128)
                plt.close()
                print(j)
            
######################
############ Martens
##            if(load[0]<load[-1]):
##                i90=np.where(load>max(load)*0.9)[0][0]
##                i50=np.where(load>max(load)*0.5)[0][0]
##                F2=[]
####                for k5 in range(len(load)):
####                    up=max(load)*0.9
####                    down=max(load)*0.5
####                    if(load[k5]<up and load[k5]>down):
####                        F2.append(load[k5])
####                F2=np.array()
##                paraM,covM=curve_fit(martens,volt[i50:i90],load[i50:i90],check_finite=False,maxfev=100000)
##                plt.plot(volt[i50:i90],load[i50:i90],volt[i50:i90],martens(volt[i50:i90],*paraM))
##                print(j)
##                print(i50,i90,paraM[0])
##                plt.savefig(os.path.join(j,"Martens.png"),dpi=128)
##                plt.close()
##

            
                
#            plt.plot(position[zz:int(len(position))],load[zz:int(len(position))])
## Finish of FD analysis
            results_files[n]=j+"/ForceDisplacement.txt"
            with open (results_files[n],"w",encoding="utf-8") as dr:
                dr.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\tDisplacement ("+ mikron+ ")\tLoad (mN)\tTip Temperature ("+degcel+")\n" )
                for p in range(len(load)):
                    dr.write(str(time[p])+"\t"+str(volt[p])+"\t"+str(load[p])+"\t"+str(tempTip[p])+"\n")
            dr.close()
            with open(j+"/RisingTCs.txt","w",encoding='utf-8') as q:
                q.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\t\tDepth (um)\t\tLoad (mN)\t\tTime Constant (ms)\n")
                for t in range (len(tautr)):
                    q.write(str(tautr[t])+"\n")
            q.close()
            with open(j+"/DecayingTCs.txt","w",encoding='utf-8') as q:
                q.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\t\tDepth (um)\t\tLoad (mN)\t\tTime Constant (ms)\n")
                for t in range (len(tautd)):
                    q.write(str(tautd[t])+"\n")
            q.close()
            if('Retract' in stepnames[n]):
                if(len(load)>3):
                    load=medfilt(load,3)
                    load_threshold+=load[-1]

                    
## Data Process for oscillation cycles -> averaging & time constant extract
        elif('Oscillation' in stepnames[n]):
            ##Median filter of load to filter incident changes
            if(len(load)>300):
                load=medfilt(load,31)
            for t in range(len(volt)): ## Touch point offset subtraction
                volt[t]=volt[t]-min_vol
            ###
            ## Find tolerances and averages of position and load to extract each cyle through data
## strain anchor points:
            sg_av=np.mean(volt)
            sg_max=max(volt)
            sg_min=min(volt)
## load anchor points:
            loadmax=max(load)
            loadmin=min(load)
            loaddel=loadmax-loadmin
            loadav=np.mean(load)
            loadtol=15
###
            lup=[]
            ldw=[]
            ups=[]
            dwn=[]
            ltrans=0

            ### Search for oscillation on load itself find "each cycle" wrt change around average
            ## in a tolerance
            for t in range(len(load)):
                if(loadmax-loadtol<load[t]<loadmax+loadtol):
                    if(ltrans==2):
                        ups.append(t)
                    if(lstart==-1):
                        lstart=t
                    ltrans=1
                    final=t
                    lup.append(load[t])
                elif(loadmin-loadtol<load[t]<loadmin+loadtol):
                    if(lstart==-1):
                        lstart=t
                    final=t
                    if(ltrans==1):
                        dwn.append(t)
                    ltrans=2
                    ldw.append(load[t])
            try:
                if not lup:
                    lups[n]=np.mean(load)+loadtol
                else:
                    lups[n]=np.mean(lup)
            except:
                lups[n]=np.mean(load)+loadtol
            try:
                if not ldw:
                    ldws[n]=np.mean(load)-loadtol
                else:
                    ldws[n]=np.mean(ldw)
            except:
                ldws[n]=np.mean(load)-loadtol
            ##
            tpt0=np.floor(np.arange(int(osc_per/2),len(volt),int(osc_per)))
            tpt1=np.floor(np.arange(int(osc_per),len(volt),int(osc_per)))
            ### Look for gage voltage whether its behavior of oscillation resoluble within tolerance
            if(sg_max-sg_min-sg_tol<(sg_av-sg_min)*2<sg_max-sg_min+sg_tol):
                print("seems good-av-min")
                starter=-1
            else:
                #print("gage problem-av-min")
                final=len(volt)
                starter=-1
            if(sg_max-sg_min-sg_tol<(sg_max-sg_av)*2<sg_max-sg_min+sg_tol):
                print("seems good-av-max")
                starter=-1
            else:
                #print("gage problem-av-max")
                starter=-1
                final=len(volt)
            ###
            ## Walk through gage data to find change points of oscillation
            if(starter==-1):
                tpt0=[]
                tpt1=[]
                transit=0
                for t in range(len(volt)):
                    if(volt[t]<=sg_av-sg_tol):
                        if(starter==-1):
                            starter=t
                        final=t
                        if(transit==2):
                            tpt0.append(t)
                        transit=1
                    elif(sg_av+sg_tol<volt[t]):
                        if(starter==-1):
                            starter=t
                        final=t
                        if(transit==1):
                            tpt1.append(t)
                        transit=2
##                print(final)
##                print(starter)
##                for t in range(len(tpt1)-1):
##                    print(tpt1[t+1]-tpt1[t])
##                for t in range(len(tpt0)-1):
##                    print(tpt0[t+1]-tpt0[t])                

                    
## Find each cycle wrt before found tpt0/1 arrays of temperature
            ### take temperature vals in arrays to sum later
            ### Find each cycle's decay and rise exponential fits
            sg=[[]]
            ld=[[]]
            tempo=[]
            tfitTip=[[]]
            tfitSample=[[]]
            fitsay=0
            taur=[]
            taud=[]            
            for t in range(starter, final):
                tfitTip[fitsay].append(tempTip[t])
                tfitSample[fitsay].append(tempSample[t])
                sg[fitsay].append(volt[t])
                if t in tpt0:
                    ared=tfitTip[fitsay]
                    tim=np.linspace(0,len(ared)/daq_frq,len(ared))
                    amax=max(ared)
                    amin=min(ared)
                    delta=amax-amin
                    ##curve fit t1 & t2 exp rise
                    try:
                        popt1, pcov1 =curve_fit(decay,tim,ared,maxfev=100000,p0=[-delta,500,amin])
                        Decayfit=decay(np.array(tim), *popt1)
                        if(popt1[1]>2000 or popt1[1]<0):
                            taud.append(float('nan'))
                        else:
                            taud.append(popt1[1])
                    except:
                        taud.append(float('nan'))
                    tfitTip.append([])
                    tfitSample.append([])
                    sg.append([])
                    fitsay+=1
                elif t in tpt1:
                    ared=tfitTip[fitsay]
                    tim=np.linspace(0,len(ared)/daq_frq,len(ared))
                    amax=max(ared)
                    amin=min(ared)
                    delta=amax-amin
                    ##curve fit t1 & t2 exp rise
                    try:
                        popt1, pcov1 =curve_fit(rise,tim,ared,maxfev=100000,p0=[delta,500,amin])
                        Risefit=rise(np.array(tim), *popt1)
                        if(popt1[1]>2000 or popt1[1]<0):
                            taur.append(float('nan'))
                        else:
                            taur.append(popt1[1])
                            #taular.append(popt1[1])
                    except:
                        taur.append(float('nan'))
                    tfitTip.append([])
                    tfitSample.append([])
                    sg.append([])
                    fitsay+=1
            taular.append(taur)
            taulad.append(taud)
            
            ## writing each cycle's TC to txt file
            with open(j+"/CyclesTCs.txt","w") as q:
                q.write(expdates[m]+" "+stepnames[n]+"\n"+"Rise TC (ms)\t\t\t Decay TC (ms)\n")
                for t in range (len(taur)):
                    if(t<len(taud)):
                        q.write(str(taur[t])+"\t"+str(taud[t])+"\n")
                    else:
                        q.write(str(taur[t])+"\n")
        
            ###
            sayacs=0
            sayacs1=0
            nor1=0
            nor2=1
            ## Find Averages of Cycles
            try:
                if(min(tpt0)<min(tpt1)):
                    nor1=1
                    nor2=0
            except:
                nor1=0
            for p in range(0,int(osc_per/2)):
                sayacs=0
                for t in range(nor1,fitsay,2):
                    if(p<len(tfitTip[t])):
                       tip[p]=tip[p]+tfitTip[t][p]
                       sample[p]=sample[p]+tfitSample[t][p]
                       voltMean[p]=voltMean[p]+sg[t][p]
                       sayacs=sayacs+1
                if(sayacs>0):
                    sample[p]=sample[p]/sayacs
                    tip[p]=tip[p]/sayacs
                    voltMean[p]=voltMean[p]/sayacs
##                    if(np.mean(tempTip)+5<tip[p] or tip[p]<np.mean(tempTip)-5):
##                        print(tip[p])
##                        tip[p]=np.mean(tempTip)
##                    if(np.mean(tempSample)+5<sample[p] or sample[p]<np.mean(tempSample)-5):
##                        print(sample[p])
##                        sample[p]=np.mean(tempSample)
##                    if(np.mean(volt)+20<voltMean[p] or voltMean[p]<np.mean(volt)-20):
##                        print(voltMean[p])
##                        voltMean[p]=np.mean(volt)
            b=int(osc_per/2)
            for p in range(b,int(osc_per)):
                sayacs1=0
                for t in range(nor2,fitsay,2):
                    if((p-b)<len(tfitTip[t])):
                       tip[p]=tip[p]+tfitTip[t][p-b]
                       sample[p]=sample[p]+tfitSample[t][p-b]
                       voltMean[p]=voltMean[p]+sg[t][p-b]
                       sayacs1=sayacs1+1           
                if(sayacs1>0):
                    sample[p]=sample[p]/sayacs1
                    tip[p]=tip[p]/sayacs1
                    voltMean[p]=voltMean[p]/sayacs1
##                    if(np.mean(tempTip)+5<tip[p] or tip[p]<np.mean(tempTip)-5):
##                        print(tip[p])
##                        tip[p]=np.mean(tempTip)
##                    if(np.mean(tempSample)+5<sample[p] or sample[p]<np.mean(tempSample)-5):
##                        print(sample[p])
##                        sample[p]=np.mean(tempSample)
##                    if(np.mean(volt)+20<voltMean[p] or voltMean[p]<np.mean(volt)-20):
##                        print(voltMean[p])
##                        voltMean[p]=np.mean(volt)
            t=0
            try:
                upss=lups[n]
                downss=ldws[n]
            except:
                upss=0
                downss=0                
            lengvol=int(len(voltMean)/2)
            if (math.isnan(upss)):
                upss=0
            if (math.isnan(downss)):
                downss=0
##            for p in range(len(voltMean)):
##                if(p<lengvol):
##                    voltMean[p]=voltMean[p]-downss/loadcell_stiffness
##                else:
##                    voltMean[p]=voltMean[p]-upss/loadcell_stiffness
            ## Write Averages of Cycles of Temperatures and Strain results to txt file !
            results_files[n]=j+"/AverageResults.txt"
            with open (results_files[n],"w",encoding="utf-8") as dr:
                dr.write(expdates[m]+" "+stepnames[n]+"\n"+"Sample Temperature ("+ degcel+ ")\tTip Temperature ("+ degcel+ ")\tStrain ("+ mikron+")\n")
                for p in range(len(sample)):
                    dr.write(str(sample[p])+"\t"+str(tip[p])+"\t"+str(voltMean[p])+"\n")
            dr.close()
        tcount=0
        n+=1
        
### Average of cycles of Steps together in 1 Graph with expRise & Decay fits
        ##separate with labels indicate time constants; !! SHOW on inset change of TC
        ##through indent or etc..
        ## ++ Load average - Position average of steps
    loadup=[]
    loaddw=[]
    depthLoad=[]
    depthUnlod=[]
    fig,ax1=plt.subplots()
    ax1.set_title("Averages of Cycles: "+expdates[m].replace('_',' '))
    fig.set_size_inches(10.0,7.5)
    fig.subplots_adjust(right=0.6)
    ax2=ax1.twinx()
    ax3=ax1.twinx()
    ax3.spines['right'].set_position(("outward",75))
    ax4=ax1.twinx()
    ax4.spines['right'].set_position(("outward",150))
    fittime=np.linspace(0,int(exp_per/2),int((osc_per/2)))
    ax1.set_ylabel("Tip Temperature (%s)"%degcel,color='r')
    ax4.set_ylabel("Load (mN)",color="cyan")
    ax4.tick_params(axis='y',labelcolor="cyan")
    ax3.tick_params(axis='y',labelcolor="orange")
    ax1.tick_params(axis='y',labelcolor="red")
    ax2.set_ylabel("Sample Temperature (%s)"%degcel,color="olive")
    ax2.tick_params(axis='y',labelcolor="olive")
    ax3.set_ylabel("Position (%s)"%mikron,color="orange")
##Uncomment/Comment out below 9 lines for tau plot inset/separate state
#    left, bottom, width, height = [0.07, 0.8, 0.15, 0.15] 
#    inset = fig.add_axes([left, bottom, width, height])
#    inset.set_xlabel("Depth (%s)"%(mikron),color="red")
#    inset.tick_params(axis='x',labelcolor="red")
#    inset.set_ylabel("Time Constant (ms)")
#    inset.grid("True",color='red')
#    inset2=inset.twiny()
#    inset2.set_xlabel("Load (mN)",color="blue")
#    inset2.tick_params(axis='x',labelcolor="blue")
    n=0
    p=0
    riseTcs=[]
    riseTcsErr=[]
    decayTcs=[]
    decayTcsErr=[]
    tempMax=[]
    tempMin=[]
    yerrR=[]
    yerrD=[]
    with open(i+"/AverTCs.txt","w") as q:
        q.write(expdates[m]+"\tAverages of Cycles Time Constants\nStep\tRise TC (ms)\t\t\tDecay TC (ms)\n")
    q.close()
    for t in stepnames:
        if('Oscillation' in t):
            salla=1
            tTip=[]
            tSam=[]
            vGag=[]
            with open(results_files[n], 'r',encoding='utf-8') as res:
                rows=res.readlines()[2:]
                for x in rows:
                    tSam.append(float(x.split()[0]))
                    tTip.append(float(x.split()[1]))
                    vGag.append(float(x.split()[2]))
            res.close()
            t0=np.linspace(0+(p*int(exp_per+1/daq_frq)),int(exp_per/2+(p*(exp_per+1/daq_frq))),int(osc_per/2))
            t02=np.linspace(int((exp_per/2+1/daq_frq)+(p*(exp_per+1/daq_frq))),int((exp_per+1/daq_frq)+(p*(exp_per+1/daq_frq))),int(osc_per/2))
            t3=np.concatenate((t0,t02))
            areddecay=tTip[0:int(len(tTip)/2)]
            aredrise=tTip[int(len(tTip)/2):]
            amax=max(tTip)
            amin=min(tTip)
            delta=amax-amin
            ax1.plot(t3,tTip,'.',color='red',zorder=2)
        ##curve fit t1 & t2 exp rise
            try:
                popt1, pcov1 =curve_fit(rise,fittime,aredrise,maxfev=100000,p0=[delta,500,amin],sigma=terr*np.ones(len(aredrise)),absolute_sigma=True)
                Risefit=rise(fittime, *popt1)
                yr1 = rise(fittime, popt1[0], max(popt1[1] - pcov1[1,1]**0.5,10), popt1[2])
                yr2 = rise(fittime, popt1[0], max(popt1[1] + pcov1[1,1]**0.5,10), popt1[2])
                riseTcs.append(popt1[1])
                riseTcsErr.append(pcov1[1,1]**0.5)
                tempMax.append(max(aredrise))
                yerrR.append(np.floor(pcov1[1,1]**0.5))
                popt2, pcov2 =curve_fit(decay,fittime,areddecay,maxfev=100000,p0=[-delta,500,amin],sigma=terr*np.ones(len(areddecay)),absolute_sigma=True)
                yerrD.append(np.floor(pcov2[1,1]**0.5))
                decayTcs.append(popt2[1])
                decayTcsErr.append(pcov2[1,1]**0.5)
                tempMin.append(min(areddecay))
                Decayfit=decay(fittime, *popt2)
                yd1 = decay(fittime, popt2[0], max(popt2[1] - pcov2[1,1]**0.5,10), popt2[2])
                yd2 = decay(fittime, popt2[0], max(popt2[1] + pcov2[1,1]**0.5,10), popt2[2])
                with open(i+"/AverTCs.txt","a",encoding="utf-8") as q:
                    q.write(stepnames[n].replace('Oscillation',"")+"\t"+str(popt1[1])
                            +" +/- "+str(yerrR[p])+" ms \t"+str(popt2[1])+" +/- "+str(yerrD[p])+" ms\n")
                q.close()
                ax1.plot(t02,Risefit,'-',color="blue",label=stepnames[n].replace("Oscillation","")+"Rise TC ="+str(np.floor(riseTcs[p]))+" +/- "+str(yerrR[p]) +"ms",zorder=2)
                ax1.plot(t02,yr1,'--',color="blue",zorder=2)
                ax1.plot(t02,yr2,'--',color="blue",zorder=2)
                ax1.fill_between(t02, yr1, yr2, facecolor="gray", alpha=0.75)
                ax1.plot(t0,Decayfit,'-',color="black",label=stepnames[n].replace("Oscillation","")+"Decay TC ="+str(np.floor(decayTcs[p]))+" +/- "+str(yerrD[p]) +"ms",zorder=2)
                ax1.plot(t0,yd1,'--',color="black",zorder=2)
                ax1.plot(t0,yd2,'--',color="black",zorder=2)
                ax1.fill_between(t0, yd1, yd2, facecolor="gray", alpha=0.75)
            except:
                if(len(riseTcs)<p+1):
                    riseTcs.append(0.)
                if(len(decayTcs)<p+1):
                    decayTcs.append(0.)
                print("Exp Fit Problems Occurred!!")
            if(riseTcs[p]>10000 or riseTcs[p]<0):
                riseTcs[p]=0.
            elif(decayTcs[p]>10000 or decayTcs[p]<0):
                decayTcs[p]=0.
            ax2.plot(t3,tSam,'.',ms=1,color="olive",zorder=1)
            ax1.set_xlabel("Average Oscillation Time (ms)")
            ax4.plot(t0,np.ones(len(t0))*ldws[n],'.',color='cyan',zorder=4)
            ax4.plot(t02,np.ones(len(t02))*lups[n],'.',color='cyan',zorder=4)
            loadup.append(lups[n])
            loaddw.append(ldws[n])
            ax3.plot(t3,vGag,'.',color="orange",zorder=2)
            p+=1
            depthLoad.append(np.nanmean(vGag[int(len(vGag)/2):]))
            depthUnlod.append(np.nanmean(vGag[:int(len(vGag)/2)]))
        n+=1
    if(salla==1):
        deltadep=[]
        deltaload=[]
        ax4.grid('True',axis='y',color='green')
        ax1.grid('True',color='red')
        ax1.legend(loc=2,fontsize='x-small')##Uncomment/Comment out for tau plot inset/separate state
        try:
            fig.tight_layout()##Uncomment/Comment out for tau plot inset/separate state
        except:
            continue
        plt.savefig(os.path.join(i,"AverageGraph.png"),dpi=1024)##Uncomment/Comment out for tau plot inset/separate state
        plt.close()    ##Uncomment/Comment out for tau plot inset/separate state
        fig=plt.figure()##Uncomment/Comment out for tau plot inset/separate state
        inset=fig.add_subplot(111,label="1")
        inDec=fig.add_subplot(111,label="2",frame_on=False)
        inDec2=fig.add_subplot(111,label="3",frame_on=False)
        inset2=fig.add_subplot(111,label="4",frame_on=False)
        inset.set_xlabel("Depth (%s)"%(mikron),color="cyan")
        inset.tick_params(axis='x',labelcolor="cyan")
        inset.tick_params(axis='y',labelcolor="red")
        inset.set_ylabel("Rise Time Constant (ms)",color="red")
        inset.grid("True",axis='y',color='red')
        inset2.set_xlabel("Load (mN)",color="blue")
        inset2.tick_params(axis='x',labelcolor="blue")
        inDec2.xaxis.tick_top()
        inDec2.yaxis.tick_right()
        inDec.yaxis.tick_right()
        inset2.xaxis.tick_top()
        inDec2.set_xticks([])
        inset2.xaxis.set_label_position('top')
        inDec.yaxis.set_label_position('right')
        for t in range(len(depthLoad)):
            deltadep.append(depthLoad[t]-depthUnlod[t])
        ytik1=riseTcs
        ytik1=[np.floor(ytik1[a]) for a in range (len(ytik1))]
        ytik1.sort()
        ytik1=np.linspace(min(ytik1)-5,max(ytik1)+5,10)
        inset.set_yticks(ytik1)
        inDec.set_ylabel("Decay Time Constant (ms)",color='orange')
        inDec.tick_params(axis='y',labelcolor="orange")
        xtik=np.concatenate((depthLoad,depthUnlod))
        xtik.sort()
        xtik=[np.floor(xtik[a]) for a in range (xtik.size)]
        xtik=np.around(np.linspace(min(xtik)-1,max(xtik)+1,10),decimals=1)
        inset.set_xticks(xtik)
        inDec.set_xticks(xtik)
        inset.set_ylim(ymin=min(ytik1)-10, ymax=max(ytik1)+10)
        inset.set_xlim(xmin=min(xtik)-sg_tol, xmax=max(xtik)+sg_tol)
        inDec.set_xlim(xmin=min(xtik)-sg_tol, xmax=max(xtik)+sg_tol)
        ytik2=decayTcs
        ytik2=[np.floor(ytik2[a]) for a in range (len(ytik2))]
        ytik2.sort()
        ytik2=np.linspace(min(ytik2)-5,max(ytik2)+5,10)
        inDec.set_yticks(ytik2)
        inDec2.set_yticks(ytik2)
        inDec2.set_ylim(ymin=min(ytik2)-10, ymax=max(ytik2)+10)
        inDec.set_ylim(ymin=min(ytik2)-10, ymax=max(ytik2)+10)
        xtik2=np.concatenate((loadup,loaddw))
        xtik2.sort()
        xtik2=[np.floor(xtik2[a]) for a in range (xtik2.size)]
        xtik2=np.around(np.linspace(min(xtik2)-5,max(xtik2)+5,10),decimals=1)
        inDec2.set_xticks(xtik2)
        inset2.set_xticks(xtik2)
        inset2.set_yticks(ytik1)
        inset2.set_ylim(ymin=min(ytik1)-10, ymax=max(ytik1)+10)
        inset2.set_xlim(right=max(xtik2)+loadtol,left=min(xtik2)-loadtol)
        inDec2.set_xlim(right=max(xtik2)+loadtol,left=min(xtik2)-loadtol)
        inDec2.axes.xaxis.set_visible(False)
        inDec2.axes.yaxis.set_visible(False)
        inDec.axes.xaxis.set_visible(False)
        inset2.axes.yaxis.set_visible(False)
        deptErr=np.ones(len(riseTcs))*sg_tol
        loadErr=np.ones(len(riseTcs))*load_threshold
        l1=inset.errorbar(depthLoad,riseTcs, yerr=yerrR, xerr=deptErr, label="Rise TC vs Depth",color='red',ls=':',marker='o',capsize=5, capthick=2)
        l2=inDec.errorbar(depthUnlod,decayTcs, yerr=yerrD, xerr=deptErr, label="Decay TC vs Depth",color='orange',ls=' ',marker='o',capsize=5, capthick=2)
        l3=inset2.errorbar(loadup,riseTcs, yerr=yerrR, xerr=loadErr, label="Rise TC vs Load",color='blue',ls=':',marker='o',capsize=5, capthick=2)
        l4=inDec2.errorbar(loaddw,decayTcs, yerr=yerrD, xerr=loadErr, label="Decay TC vs Load",color='green',ls=' ',marker='o',capsize=5, capthick=2)
        [caps.set_alpha(0.3) for caps in l1[2]]
        [caps.set_alpha(0.3) for caps in l2[2]]
        [caps.set_alpha(0.3) for caps in l3[2]]
        [caps.set_alpha(0.3) for caps in l4[2]]
        [bars.set_alpha(0.3) for bars in l1[1]]
        [bars.set_alpha(0.3) for bars in l2[1]]
        [bars.set_alpha(0.3) for bars in l3[1]]
        [bars.set_alpha(0.3) for bars in l4[1]]
        fig.legend(handles=[l1,l2,l3,l4],fontsize='x-small',ncol=2,mode='expand',loc=2)
        fig.suptitle(" Time Constants of Averages of the Cycles: "+expdates[m].replace('_',' '),fontsize='small')
        try:
            fig.tight_layout()
        except:
            continue
        plt.savefig(os.path.join(i,"avTCGraph.png"),dpi=512)##Uncomment/Comment out for tau plot inset/separate state
        plt.close()
    ###
        ## Writing rise tau vs depth & load to file
        with open(i+"/riseTC_disp_load.txt","w",encoding="utf-8") as q:
            q.write("Time Constant (ms)\tSigma (ms)\t Temperature ("+degcel+")\tDepth ("+mikron+")\tLoad (mN)\n")
            for t in range(len(riseTcs)):
                q.write(str(riseTcs[t])+"\t"+str(riseTcsErr[t])+"\t"+str(tempMax[t])+"\t"+str(depthLoad[t])+"\t"+str(loadup[t])+"\n")
        q.close()
        ## Writing decay tau vs depth & load to file
        with open(i+"/decayTC_disp_load.txt","w",encoding="utf-8") as q:
            q.write("Time Constant (ms)\tSigma (ms)\t Temperature ("+degcel+")\tDepth ("+mikron+")\tLoad (mN)\n")
            for t in range(len(decayTcs)):
                q.write(str(decayTcs[t])+"\t"+str(decayTcsErr[t])+"\t"+str(tempMin[t])+"\t"+str(depthUnlod[t])+"\t"+str(loaddw[t])+"\n")
        q.close()
    ## Plot of each cycle's tau of whole experiment's oscillations
        ax=plt.gca()
        fig,plto=plt.subplots()
        for t in range(len(taular)):
            print(t)
            color=next(ax._get_lines.prop_cycler)['color']
            risedots=np.arange(1+t*(len(taular[t])),1+len(taular[t])+(len(taular[t])*t))
            plto.plot(risedots,taular[t],'.',color=color)
            riseav=np.nanmean(taular[t])*np.ones(len(taular[t]))
            plto.plot(risedots,riseav,'-',label="Average Rise Step "+str(t+1)+" ="+str(np.floor(np.nanmean(taular[t])))+" ms",color=color)
            color=next(ax._get_lines.prop_cycler)['color']
            decaydots=np.arange(1+t*(len(taulad[t])),1+len(taulad[t])+(len(taulad[t])*t))
            plto.plot(decaydots,taulad[t],'.',color=color)
            decayav=np.nanmean(taulad[t])*np.ones(len(taulad[t]))
            plto.plot(decaydots,decayav,'-',label="Average Decay Step "+str(t+1)+" ="+str(np.floor(np.nanmean(taulad[t])))+" ms",color=color)
        plto.legend()
        plto.set_xlabel("Number of Cycle")
        plto.set_ylabel("Time Constant (ms)")
        plto.set_title(expdates[m]+" Time Constants of Each Cycle")
        plt.savefig(os.path.join(i,"allTCGraph.png"),dpi=512)
        plt.close()
    else:
        plt.close()
    salla=0
#all-in-one displacement vs load&temperature graphs double yaxis multi x axis-subplots
#all-in-one time vs disp&load&temperature graph triple yaxis sole x axis
    n=0
    p=0
    inden=0
    tDis=[]
    tLoad=[]
    tTip=[]
    axt=[]
    for t in stepnames:
        if('Indentation' in t or 'Retract' in t):
            p+=1
    figs =plt.figure()
    figs.suptitle("Temperature and Load vs Displacement: "+expdates[m].replace('_',' '),fontsize=16,y=1.)
    common=figs.add_subplot(111,label='1')
    com1=common.twinx()
    common.set_xlabel('Displacement (%s)'%mikron)
    common.tick_params(color='white',top=False, bottom=False, right=False, left=False)
    com1.spines['bottom'].set_color('white')
    com1.tick_params(color='white',top=False, bottom=False, left=False,right=False)
    for t in range(p):
        if(9>p>4):
            axt.append(figs.add_subplot(2,4,t+1))
##            if(4>t>0):
##                axt[t].sharey(axt[0])
##            elif(t>4):
##                axt[t].sharey(axt[4])
        elif(13>p>8):
            axt.append(figs.add_subplot(3,4,t+1))
##            if(4>t>0):
##                axt[t].sharey(axt[0])
##            elif(9>t>4):
##                axt[t].sharey(axt[4])
##            elif(t>8):
##                axt[t].sharey(axt[8])
        elif(17>p>12):
            axt.append(figs.add_subplot(4,4,t+1))
##            if(4>t>0):
##                axt[t].sharey(axt[0])
##            elif(9>t>4):
##                axt[t].sharey(axt[4])
##            elif(13>t>8):
##                axt[t].sharey(axt[8])
##            elif(17>t>12):
##                axt[t].sharey(axt[12])
        else:            
            axt.append(figs.add_subplot(1,p,t+1))
##            if(t>0):
##                axt[t].sharey(axt[0])
    figs.set_size_inches(p*2.0,6*int(np.ceil(p/4)))
    axd=[]
    for t in range(p):
        axd.append(axt[t].twinx())
##        if(t>0):
##            axd[t].sharey(axd[0])
###        axt[t].get_yaxis().set_visible(False)
        axt[t].spines['left'].set_color('red')
        axd[t].spines['right'].set_color('blue')
        axd[t].tick_params(color='blue', left=False, right=True)
        axt[t].tick_params(color='red', right=False, left=True)
        axt[t].grid('True',axis='y',color='blue')
        axd[t].grid('True',axis='y',color='red')
    common.set_xticks([])
    com1.set_xticks([])
    common.set_yticks([])
    com1.set_yticks([])
    common.set_ylabel("Tip Temperature (%s)"%degcel,color='r')
    com1.set_ylabel("Load (mN)",color="blue")
    common.spines['left'].set_position(("outward",50))
    common.spines['bottom'].set_position(("outward",50))
    com1.spines['right'].set_position(("outward",50))
#    common.sharey(axt[0])
 #   com1.sharey(axd[0])
 #   com1.grid('True',axis='y',color='blue')
 #   common.grid('True',axis='y',color='red')
    ##21.01.2022-single graph addition:
    singleTemp=[]
    singleDis=[]
    singleTim=[]
    singleLoad=[]
    multiTim=[]
    multiTemp=[]
    multiDis=[]
    multiLoad=[]
    tim0=0
    p=0
    for t in stepnames:
        if('Indentation' in t or 'Approach' in t or 'Retract' in t):
            inden=1
            if('Indentation' in t and n>0 and 'Approach' in stepnames[n-1]):
                zz=0
            else:
                tDis=[]
                tLoad=[]
                tTip=[]
            zz=0
            multiTim.append([])
            multiDis.append([])
            multiTemp.append([])
            multiLoad.append([])
            with open(results_files[n], 'r',encoding='utf-8') as res:
                rows=res.readlines()[2:]
                for x in rows:
                    singleTim.append(float(x.split()[0])+tim0)
                    multiTim[p].append(float(x.split()[0]))
                    tDis.append(float(x.split()[1]))
                    singleDis.append(float(x.split()[1]))
                    multiDis[p].append(float(x.split()[1]))
                    tLoad.append(float(x.split()[2]))
                    singleLoad.append(float(x.split()[2]))
                    multiLoad[p].append(float(x.split()[2]))
                    tTip.append(float(x.split()[3]))
                    singleTemp.append(float(x.split()[3]))
                    multiTemp[p].append(float(x.split()[3]))
            res.close()
            try:
                tim0=singleTim[-1]
            except:
                print("No TIME")
            if('Approach' in t):
                zz=0
            else:
                try:                       
                    axt[p].plot(tDis[zz:int(len(tDis))],tTip[zz:int(len(tDis))],'r.')
                    xtik=np.around(np.linspace(np.floor(min(tDis)-0.5),np.floor(max(tDis)+0.5),4),decimals=1)
                    axt[p].set_xticks(xtik)
                    axt[p].set_title(t)
                    if(len(tLoad)>0):
                        axd[p].plot(tDis[zz:int(len(tDis))],tLoad[zz:int(len(tDis))],'b.')
                except:
                    print("load-disp-problem on"+t)
                p+=1
        n+=1
    if(inden==1):
        if(p>1):
            singleLoad=medfilt(singleLoad,101)
            figs.tight_layout()
            figs.subplots_adjust(top=(0.88-0.1*(p/9)))
            plt.savefig(os.path.join(i,"LoadvsDisp.png"),dpi=256)
            plt.figure()
            singleDis=np.array(singleDis)
            singleTemp=np.array(singleTemp)
            plt.scatter(singleDis,singleTemp)
            plt.xlabel('Displacement (%s)'%mikron)
            plt.ylabel('Temperature (%s)'%degcel)
            plt.savefig(os.path.join(i,"TemperatureAll_Depth.png"),dpi=256)
            plt.figure()
            singleTim=np.array(singleTim)
            plt.scatter(singleTim,singleTemp)
            plt.xlabel('Time (s)')
            plt.ylabel('Temperature (%s)'%degcel)
            plt.savefig(os.path.join(i,"TemperatureAll_Time.png"),dpi=256)
            fig,ax=plt.subplots()
            ax1=ax.twinx()
            ax.scatter(singleTim,singleDis,label="Displacement",color="orange")
            ax1.scatter(singleTim,singleLoad,label="Load",color="cyan")
            ax1.yaxis.set_label("Load (mN)")
            ax.set_ylabel("Displacement (%s)"%mikron)
            ax1.set_ylabel("Load (mN)")
            ax.set_xlabel('Time (s)')
            fig.legend()
            plt.savefig(os.path.join(i,"DispAll.png"),dpi=256)
            plt.figure()
            plt.scatter(singleDis,singleLoad)
            plt.xlabel('Displacement (%s)'%mikron)
            plt.ylabel('Load (mN)')
            plt.savefig(os.path.join(i,"LoadDisp-Single.png"),dpi=256)
            plt.figure()
            for inc in range(p):
                plt.scatter(multiTim[inc],medfilt(multiTemp[inc],101))
            plt.xlabel('Time (s)')
            plt.ylabel('Temperature (%s)'%degcel)
            plt.savefig(os.path.join(i,"TemperatureAll_Mtime.png"),dpi=256)
            plt.figure()
            for inc in range(p):
                plt.scatter(multiTim[inc],medfilt(multiLoad[inc],101))
            plt.xlabel('Time (s)')
            plt.ylabel('Load (mN)')
            plt.savefig(os.path.join(i,"TemperatureAll_Mload.png"),dpi=256)
            plt.figure()
            for inc in range(p):
                plt.scatter(multiTim[inc],medfilt(multiDis[inc],101))
            plt.xlabel('Time (s)')
            plt.ylabel('Displacement (%s)'%mikron)
            plt.savefig(os.path.join(i,"TemperatureAll_Mdis.png"),dpi=256)                                        
        plt.close()
    m+=1
