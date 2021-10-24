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

"""
This code is for reading experimental data of Indenter/RC Thermal setup obtained by
NI-DAQ written/converted txt data files are being searched as txt. 
"""
## Calibration constants of Experiment ##
lc_mass_cal=4.75
grave=9.81
bit_to_load=lc_mass_cal*grave
loadcell_stiffness=536.5 ##mN/um
gage2um=10 ##um/V
load_threshold=15 ##mN loadcell sensitivity for touch
terr=1.1#error for temperature

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
###

##Special Characters##
degsin=u"\N{DEGREE SIGN}"
degcel=degsin+"C"
deltasig=u"\N{GREEK CAPITAL LETTER DELTA}"
mikron=u"\N{GREEK SMALL LETTER MU}"+"m"
###

## Main file of Expriment Data of text files -that converted from NI TDMS-
#main="D:/ahmed/RC experiments/Al/Al-Process/September/25-09-2021" #WindowsFilePath
main="/Users/ahmeduluca/Desktop/download/process"
 #MacFilePath
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
tip_chname="Temperature2"
sample_chname='Temperature3'
gage_chname="Voltage0"
pos_chname="Actuator_Voltage"
salla=0
load=[]
position=[]
volt=[]
voltaj=[]
temp=[[],[]]
time=[]
min_vol=0
sg_tol=0.25
## Search for Experiment Files that include Step Files ##
for i in file_list:
    print(file_list[m])
    print(expdates[m])
    sub_dirs=[f.path for f in os.scandir(i) if f.is_dir()]
    sub_dirs.sort()
    stepnames=[f.name for f in os.scandir(i) if f.is_dir()]
    stepnames.sort()
    results_files=['']*len(stepnames)   
    n=0
    taular=[]
    taulad=[]
    lups=np.zeros((len(stepnames)))
    ldws=np.zeros((len(stepnames)))
    ## Search for txt files of Load - Strain Gage - Temperatures ##    
    for j in sub_dirs:
        load=[]
        position=[]
        reduceT=[]
        redtim=[]
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
        for k in txt_files:
            if("Load" in k):
                load=[]
                f=open(k,"r",encoding='ascii')
                for x in f:
                    load.append(float(x)*bit_to_load/1000)
                f.close()
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
        volt=medfilt(voltaj,21)
        
##Indentation load displacement graphs & timeconstant extraction trial through indentation stage
        if('Indentation' in stepnames[n] or 'Calibration' in stepnames[n] or 'Approach' in stepnames[n] or 'Retract' in stepnames[n]):
            volt=medfilt(volt,501)
            if(len(load)>0):
                load=medfilt(load,31)
                a=0
                for kuvvet in load:
                    a+=1
                    if(kuvvet>load_threshold):
                        break
                oran=1-(len(load)-a)/len(load)
                b=int(oran*len(volt))-1
                orta=int(len(volt)/len(load))
            else:
                b=0
                orta=1
            tempTip=medfilt(tempTip,21)
            if(n==0 or (n>0 and 'Retract' in stepnames[n-1]) or 'Approach' in stepnames[n]):
                min_vol=volt[b]
            zz=0
            if("Approach" in stepnames[n]):
                mxpos=max(volt)
                mnpos=min(volt)
                dpoz=mxpos-mnpos
                if(volt[b]>mnpos):
                    for t in range(len(volt)-1,-1,-1):
                        if(volt[t]-sg_tol<mnpos):
                            zz=t
                            break
            zz=0
            for t in range(len(volt)):
                volt[t]=volt[t]-min_vol
            figs,axd=plt.subplots()
            axd.set_title(stepnames[n]+" "+expdates[m].replace('_',' '))
            if(len(load)>0):
                axl=axd.twinx()
                axl.plot(np.linspace(0,time[-1],len(load)),load,'.',color='cyan')
                axl.set_ylabel("Load (mN)",color="cyan")
                axl.tick_params(axis='y',labelcolor="cyan")
            axt=axd.twinx()
            axt.spines['right'].set_position(("outward",75))
            axd.set_ylabel("Displacement (%s)"%mikron,color='orange')
            axd.set_xlabel("Time (s)",color="black")
            axt.set_ylabel("Tip Temperature (%s)"%degcel,color="red")
            axd.tick_params(axis='y',labelcolor="orange")
            axt.tick_params(axis='y',labelcolor="red")
            axd.plot(time,volt,'.',color='orange')
            axt.plot(time,tempTip,'r.')
            figs.tight_layout()
            plt.savefig(os.path.join(j,"TimevsLoad&Disp.png"),dpi=128)
            plt.close()
            
            ## Reduced data for load vs graphs:
            p2=[]
            p3=[]
            p4=[]
            for t in range(0,int(len(volt)-orta/2),orta):
                p2.append(volt[int(t+(orta/2))])
                p3.append(tempTip[int(t+(orta/2))])
                p4.append(time[int(t+(orta/2))])
            if(len(load)>0):               
                leng=len(p2)-1
                artik=leng%len(load)
                if(artik==0):
                    artik=1
                out=int(leng/artik)
                if(out==0):
                    out=1
                bas=leng%artik
            else:
                bas=0
                leng=len(p2)-1
                out=1
                artik=1
            for t in range(bas,leng,out):
                position=np.concatenate((np.array(position),p2[t:t+out-1]))
                reduceT=np.concatenate((np.array(reduceT),p3[t:t+out-1]))
                redtim=np.concatenate((np.array(redtim),p4[t:t+out-1]))
            inc=0
            gg=0
            for t in range(len(position)):
                position[t]=position[t]-load[t]/loadcell_stiffness
            position=medfilt(position,101)
            figs,axt=plt.subplots()
            axt.set_title(stepnames[n]+" "+expdates[m].replace('_',' '))
            axt.plot(position[zz:int(len(position))],reduceT[zz:int(len(position))],'r.')
            if(len(load)>0):
                axd=axt.twinx()
                axd.plot(position[zz:int(len(position))],load[zz:int(len(position))],'b.')
            plt.savefig(os.path.join(j,"LoadvsDisp.png"),dpi=128)
            plt.close()
            maxT=max(tempTip)
            minT=min(tempTip)
            #steps of actuator are at each 100ms...for starting time/100 maybe tried
            dTemp=0.2
            number=(maxT-minT)/dTemp
##            lastpt=np.where(tempTip<=minT)
##            lastpt=int(lastpt[0][-1])
            lastpt=b
            tim=[]
            say=0
            tautr=[]
            tautd=[]
            starttime=time[b]
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
            results_files[n]=j+"/ForceDisplacement.txt"
            with open (results_files[n],"w+",encoding="utf-8") as dr:
                dr.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\tDisplacement ("+ mikron+ ")\tLoad (mN)\tTip Temperature ("+degcel+")\n" )
                for p in range(len(position)):
                    dr.write(str(redtim[p])+"\t"+str(position[p])+"\t"+str(load[p])+"\t"+str(reduceT[p])+"\n")
            dr.close()
            with open(j+"/RisingTCs.txt","w+",encoding='utf-8') as q:
                q.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\t\tDepth (um)\t\tLoad (mN)\t\tTime Constant (ms)\n")
                for t in range (len(tautr)):
                    q.write(str(tautr[t])+"\n")
            q.close()
            with open(j+"/DecayingTCs.txt","w+",encoding='utf-8') as q:
                q.write(expdates[m]+" "+stepnames[n]+"\n"+"Time (s)\t\tDepth (um)\t\tLoad (mN)\t\tTime Constant (ms)\n")
                for t in range (len(tautd)):
                    q.write(str(tautd[t])+"\n")
            q.close()
            if('Retract' in stepnames[n]):
                if(len(load)>11):
                    load=medfilt(load,21)
                    load_threshold+=load[-1]

                    
## Data Process for oscillation cycles -> averaging & time constant extract
        elif('Oscillation' in stepnames[n]):
            ##Median filter of load to filter incident changes
            if(len(load)>300):
                load=medfilt(load,31)
            for t in range(len(volt)):
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
            with open(j+"/CyclesTCs.txt","w+") as q:
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
            for p in range(len(voltMean)):
                if(p<lengvol):
                    voltMean[p]=voltMean[p]-downss/loadcell_stiffness
                else:
                    voltMean[p]=voltMean[p]-upss/loadcell_stiffness
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
    decayTcs=[]
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
                popt1, pcov1 =curve_fit(rise,fittime,aredrise,maxfev=100000,p0=[delta,500,amin],sigma=terr*np.ones(len(aredrise)),absolute_sigma=False)
                Risefit=rise(fittime, *popt1)
                yr1 = rise(fittime, popt1[0], max(popt1[1] - pcov1[1,1]**0.5,10), popt1[2])
                yr2 = rise(fittime, popt1[0], max(popt1[1] + pcov1[1,1]**0.5,10), popt1[2])
                riseTcs.append(popt1[1])
                yerrR.append(np.floor(pcov1[1,1]**0.5))
                popt2, pcov2 =curve_fit(decay,fittime,areddecay,maxfev=100000,p0=[-delta,500,amin],sigma=terr*np.ones(len(areddecay)),absolute_sigma=False)
                yerrD.append(np.floor(pcov2[1,1]**0.5))
                decayTcs.append(popt2[1])
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
        loadErr=np.ones(len(riseTcs))*loadtol
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
            q.write("Time Constant (ms)\tDepth ("+mikron+")\tLoad (mN)\n")
            for t in range(len(riseTcs)):
                q.write(str(riseTcs[t])+"\t"+str(depthLoad[t])+"\t"+str(loadup[t])+"\n")
        q.close()
        ## Writing decay tau vs depth & load to file
        with open(i+"/decayTC_disp_load.txt","w",encoding="utf-8") as q:
            q.write("Time Constant (ms)\tDepth (%s)\tLoad (mN)\n"%mikron)
            for t in range(len(decayTcs)):
                q.write(str(decayTcs[t])+"\t"+str(depthUnlod[t])+"\t"+str(loaddw[t])+"\n")
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
    common=figs.add_subplot(111,label='1')
    com1=figs.add_subplot(111,label='2')
    com1.sharex(common)
    common.set_xlabel('Displacement (%s)'%mikron)
    common.spines['top'].set_color('none')
    common.spines['bottom'].set_color('none')
    common.tick_params(top=False, bottom=False, right=False)
    com1.spines['top'].set_color('none')
    com1.spines['bottom'].set_color('none')
    com1.tick_params(top=False, bottom=False, left=False)
    for t in range(p):
        axt.append(figs.add_subplot(1,p,t+1))
        if(t>0):
            axt[t].sharey(axt[0])
    figs.set_size_inches(p*2.0,6)
    axd=[]
    for t in range(p):
        axd.append(axt[t].twinx())
        if(t>0):
            axd[t].sharey(axd[0])
            axt[t].spines['left'].set_color('none')
            axt[t].spines['right'].set_color('none')
        axt[t].spines['top'].set_color('none')
        if(t<p-1):
            axd[t].spines['left'].set_color('none')
            axd[t].spines['right'].set_color('none')
        axd[t].spines['top'].set_color('none')
        axd[t].tick_params(left=False, right=False)
        axt[t].tick_params(right=False, left=False)
    figs.suptitle("Temperature and Load vs Displacement: "+expdates[m].replace('_',' '))
    common.set_ylabel("Tip Temperature (%s)"%degcel,color='r')
    com1.set_ylabel("Load (mN)",color="blue")
    common.sharey(axt[0])
    com1.sharey(axd[0])
    com1.grid('True',axis='y',color='blue')
    common.grid('True',axis='y',color='red')
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
            with open(results_files[n], 'r',encoding='utf-8') as res:
                rows=res.readlines()[2:]
                for x in rows:
                    tDis.append(float(x.split()[1]))
                    tLoad.append(float(x.split()[2]))
                    tTip.append(float(x.split()[3]))
            res.close()
            if('Approach' in t):
                zz=0
            else:
                axt[p].plot(tDis[zz:int(len(tDis))],tTip[zz:int(len(tDis))],'r.')
                if(len(tLoad)>0):
                    axd[p].plot(tDis[zz:int(len(tDis))],tLoad[zz:int(len(tDis))],'b.')
                p+=1
        n+=1
    if(inden==1):
        figs.tight_layout()
        plt.savefig(os.path.join(i,"LoadvsDisp.png"),dpi=256)
        plt.close()
    m+=1
