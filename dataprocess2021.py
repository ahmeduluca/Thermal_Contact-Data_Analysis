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
from scipy.signal import medfilt

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
###

load_files=[]
position_files=[]
temp_files=[]
sg_files=[]
lvsdis_file=[[]]

date=[[]]
n=0
a=[]
b=[]
c=[]
#t1=np.linspace(1,1000,1000)
#t2=np.linspace(999,1999,1000)


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
main="/Users/ahmeduluca/Desktop/Cu/Cu-Process/22-06-2021"
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
oscsay=0
## Search for Experiment Files that include Step Files ##
for i in file_list:
    load=[[]]
    position=[[]]
    voltage=[[]]
    temp=[[[],[]]]
    time=[[]]
    print(file_list[m])
    print(expdates[m])
    sub_dirs=[f.path for f in os.scandir(i) if f.is_dir()]
    sub_dirs.sort()
    stepnames=[f.name for f in os.scandir(i) if f.is_dir()]
    stepnames.sort()
    lups=np.zeros((len(stepnames)))
    ldws=np.zeros((len(stepnames)))
    tip=np.empty((len(stepnames),842))
    sample=np.empty((len(stepnames),842))
    voltMean=np.empty((len(stepnames),842))    
    n=0
    taular=[]
    taulad=[]
    ## Search for txt files of Load - Strain Gage - Temperatures ##    
    for j in sub_dirs:

        txt_files=glob.glob(os.path.join(j+'/*.txt'))
        ch=0
        for k in txt_files:
            if("Load" in k):
                load.append([])
                f=open(k,"r")
                for x in f:
                    load[n].append(float(x)*bit_to_load/1000)
            elif("Position" in k):
                position.append([])
                f=open(k,"r")
                for x in f:
                    position[n].append(float(x))
            elif("Temperature" in k):
                if(ch==0 and n>0):
                    temp.append([[],[]])
                if(tcount==0):
                    time.append([])
                with open(k,'r') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time[n].append(float(x.split()[0]))
                        temp[n][ch].append(float(x.split()[1]))
                tcount+=1
                ch+=1
            elif("Voltage" in k):
                voltage.append([])
                if(tcount==0):
                    time.append([])
                with open(k,'r') as f:
                    rows=f.readlines()
                    for x in rows:
                        if(tcount==0):
                            time[n].append(float(x.split()[0]))
                        voltage[n].append(float(x.split()[1])*gage2um)##um conversion
                tcount+=1
        voltage[n]=medfilt(voltage[n],21)
        if(n==0):
            min_vol=min(voltage[0])
        for t in range(len(voltage[n])):
            voltage[n][t]=voltage[n][t]-min_vol
            
        ## Data Process for oscillation cycles -> averaging & time constant extract
        if('Oscillation' in stepnames[n]):
            ##Median filter of load to filter incident changes
            load[n]=medfilt(load[n],101)
            ###
            ## Find tolerances and averages of position and load to extract each cyle through data
## strain anchor points:
            sg_tol=0
            #plt.plot(voltage[n])
            #plt.show()
#            for t in range(1000):
#                sg_tol+=abs(voltage[n][t]-voltage[n][t+1])
 #           sg_tol=sg_tol/1000
            #print(sg_tol)
            sg_tol=0.4
            sg_av=np.mean(voltage[n])
            sg_max=max(voltage[n])
            sg_min=min(voltage[n])
## load anchor points:
            loadmax=max(load[n])
            loadmin=min(load[n])
            loaddel=loadmax-loadmin
            loadav=np.mean(load[n])
            loadtol=15
###
            lup=[]
            ldw=[]
            ups=[]
            dwn=[]
            ltrans=0

            ### Search for oscillation on load itself find "each cycle" wrt change around average
            ## in a tolerance
            for t in range(len(load[n])):
                if(loadmax-loadtol<load[n][t]<loadmax+loadtol):
                    if(ltrans==2):
                        ups.append(t)
                    if(lstart==-1):
                        lstart=t
##                        else:
##                            plt.plot(load[n][:lup[-1]])
##                            plt.show()
                    ltrans=1
                    final=t
                    lup.append(load[n][t])
                elif(loadmin-loadtol<load[n][t]<loadmin+loadtol):
                    if(lstart==-1):
                        lstart=t
                    final=t
                    if(ltrans==1):
                        dwn.append(t)
                    ltrans=2
                    ldw.append(load[n][t])
            lups[n]=np.mean(lup)
            ldws[n]=np.mean(ldw)
            ##
            
            ### Look for gage voltage whether its behavior of oscillation resoluble within tolerance
            if(sg_max-sg_min-sg_tol<(sg_av-sg_min)*2<sg_max-sg_min+sg_tol):
                #print("seems good-av-min")
                starter=-1
            else:
                #print("gage problem-av-min")
                final=len(voltage[n])
                starter=-1
            if(sg_max-sg_min-sg_tol<(sg_max-sg_av)*2<sg_max-sg_min+sg_tol):
                #print("seems good-av-max")
                starter=-1
            else:
                #print("gage problem-av-max")
                starter=-1
                final=len(voltage[n])
            ###
            ## Walk through gage data to find change points of oscillation
            if(starter==-1):
                tpt0=[]
                tpt1=[]
                transit=0
                t1=[]
                t2=[]
                tempo=[]
                for t in range(len(voltage[n])):
                    if(sg_av<voltage[n][t]):
                        if(starter==-1):
                            starter=t
                        final=t
                        if(transit==2):
                            tpt0.append(t)
                        transit=1
                    elif(sg_av>voltage[n][t]):
                        if(starter==-1):
                            starter=t
                        final=t
                        if(transit==1):
                            tpt1.append(t)
                        transit=2
                if(min(tpt1)-starter>842):
                    starter=min(tpt1)-842
#            else: ## If it is not well enough use a generic data pts-its for 250mHz osc & 200Hz daq
#                starter=0
#                final=len(voltage[n])
#                tpt0=[835, 1676, 2518, 3360, 4202, 5044, 5885, 6727, 7569, 8410, 9253, 10095, 10937, 11779, 12621, 13464, 14306, 15148, 15990, 16832, 17674, 18516, 19358, 20200, 21041, 21884, 22726, 23567, 24409, 25251, 26093, 26935, 27776, 28618, 29460, 30301, 31143, 31985, 32826, 33668, 34510, 35352, 36193, 37035, 37877, 38719, 39561, 40403, 41245, 42087, 42929, 43771, 44613, 45455, 46297, 47139, 47981, 48823, 49665, 50506, 51349, 52191, 53033, 53875, 54717, 55559, 56401, 57243, 58085, 58927, 59769, 60611, 61453, 62295, 63137, 63979, 64822, 65664, 66506, 67348, 68190, 69032, 69874, 70716, 71558, 72400, 73242, 74084, 74926, 75768]
#                tpt1=[1255, 2097, 2939, 3781, 4623, 5464, 6306, 7148, 7989, 8831, 9674, 10516, 11358, 12200, 13043, 13885, 14727, 15569, 16411, 17253, 18095, 18937, 19779, 20621, 21463, 22305, 23147, 23988, 24830, 25672, 26514, 27356, 28197, 29039, 29880, 30722, 31564, 32405, 33247, 34089, 34931, 35773, 36615, 37456, 38298, 39140, 39982, 40824, 41666, 42508, 43350, 44192, 45034, 45876, 46718, 47560, 48402, 49244, 50086, 50927, 51770, 52612, 53454, 54296, 55138, 55980, 56822, 57664, 58506, 59348, 60190, 61033, 61874, 62716, 63558, 64401, 65243, 66085, 66927, 67769, 68611, 69453, 70295, 71137, 71979, 72821, 73663, 74505, 75347]
##                for t in range(len(tpt0)):
##                    tpt0[t]=2*tpt0[t]
##                for t in range(len(tpt1)):
##                    tpt1[t]=2*tpt1[t]
            #print(len(tpt0))
            #print(len(tpt1))
            #print(starter)
            #print(final)
            
## Find each cycle wrt before found tpt0/1 arrays of temperature
            ### take temperature vals in arrays to sum later
            ### Find each cycle's decay and rise exponential fits
            sg=[[]]
            ld=[[]]
            tempo=[]
            tfit1=[[]]
            tfit2=[[]]
            fitsay=0
            taur=[]
            taud=[]            
            for t in range(starter, final):
                t1.append(temp[n][0][t])
                t2.append(temp[n][1][t])
                tempo.append(time[n][t])
                tfit1[fitsay].append(temp[n][0][t])
                tfit2[fitsay].append(temp[n][1][t])
                sg[fitsay].append(voltage[n][t])
                if t in tpt0:
                    tim=np.linspace(0,len(t2)*5,len(t2))
                    ared=t2
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
                    t1=[]
                    t2=[]
                    tempo=[]
                elif t in tpt1:
                    tim=np.linspace(0,len(t2)*5,len(t2))
                    ared=t2
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
                    tfit1.append([])
                    tfit2.append([])
                    sg.append([])
                    fitsay+=1
                    t1=[]
                    t2=[]
                    tempo=[]
            taular.append(taur)
            taulad.append(taud)
##            linpara, lincov=curve_fit(linear,np.arange(0,len(taur)),taur,check_finite=False,maxfev=100000,p0=[1,np.nanmean(taur)])
##            Linearfit=linear(np.arange(0,len(taur),0.1),*linpara)
##            plt.plot(np.arange(0,len(taur),0.1),Linearfit,label="Rise Linear Fit m="+str(round(linpara[0])))
##            plt.plot(np.arange(0,len(taur)),np.nanmean(taur)*np.ones(len(taur)),'r-',label="Average Rise ="+str(round(np.nanmean(taur)))+" ms")
##            plt.plot(np.arange(0,len(taud)),np.nanmean(taud)*np.ones(len(taud)),'b-',label="Average Decay ="+str(round(np.nanmean(taud)))+" ms")
##            plt.plot(taur,'r.',label="Rise TC")
##            plt.plot(taud,'b.',label="Decay TC")
##            plt.legend()
##            plt.show()
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

            ## Find Averages of Cycles of
            #for t in range(fitsay):3
#                plt.plot(tfit2[fitsay])
                #print(len(tfit2[fitsay]))
#            plt.show()
            for p in range(0,842):
                sayacs=0
                for t in range(fitsay):
                    #print("Length of average:"+str(len(tfit2[t])))
                    if(len(tfit1[t])>p):
                        sample[n][p]+=tfit1[t][p]
                        tip[n][p]+=tfit2[t][p]
                        voltMean[n][p]+=sg[t][p]
                        sayacs+=1                         
                sample[n][p]=sample[n][p]/sayacs
                tip[n][p]=tip[n][p]/sayacs
                #print(tip[n][p])
                voltMean[n][p]=voltMean[n][p]/sayacs
            t=0
#            print(lups[n])
            upss=lups[n]
            downss=ldws[n]
            lengvol=int(len(voltMean[n])/2)
#            print("Length:"+str(lengvol))
            for p in range(len(voltMean[n])):
                if(t<lengvol):
                    voltMean[n][p]-=downss/536.5
                else:
                    voltMean[n][p]-=upss/536.5
            oscsay+=1
            ## Write Averages of Cycles of Temperatures and Strain results to txt file !
            results=j+"/AverageResults.txt"
            with open (results,"w") as dr:
                dr.write(expdates[m]+" "+stepnames[n]+"\n"+"Sample Temperature (C)\tTip Temperature (C)\tStrain (um)\n")
                for p in range(len(sample[n])):
                    dr.write(str(sample[n][p])+"\t"+str(tip[n][p])+"\t"+str(voltMean[n][p])+"\n")
            dr.close()
        #print(stepnames[n])
        #print(j)
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
    fig.subplots_adjust(right=0.6)
    ax2=ax1.twinx()
    ax3=ax1.twinx()
    ax3.spines['right'].set_position(("axes",1.25))
    ax4=ax1.twinx()
    ax4.spines['right'].set_position(("axes",1.5))
    fittime=np.linspace(0,2105,421)
    ax1.set_ylabel("Tip Temperature ("+degcel+")",color='r')
    ax1.tick_params(axis='y',labelcolor="red")
    ax2.set_ylabel("Sample Temperature ("+ degcel+")",color="olive")
    ax2.tick_params(axis='y',labelcolor="olive")
    ax3.set_ylabel("Position ("+mikron+")",color="orange")
    ax3.tick_params(axis='y',labelcolor="orange")
    ax4.set_ylabel("Load (mN)",color="green")
    ax4.tick_params(axis='y',labelcolor="green")
    left, bottom, width, height = [0.07, 0.8, 0.15, 0.15]
    inset = fig.add_axes([left, bottom, width, height])
    inset.set_xlabel("Amplitude_pp (%s)"%(mikron),color="red")
    inset.tick_params(axis='x',labelcolor="red")
    inset.set_ylabel("Time Constant (ms)")
    inset.grid("True",color='red')
    inset2=inset.twiny()
    inset2.set_xlabel("Load (mN)",color="blue")
    inset2.tick_params(axis='x',labelcolor="blue")
    ax1.set_title(expdates[m]+" Averages of Cycles")
    n=0
    p=0
    riseTcs=[]
    decayTcs=[]
    with open(i+"/AverTCs.txt","w") as q:
        q.write(expdates[m]+"\tAverages of Cycles Time Constants\n")
    q.close()
    for t in stepnames:
        if('Oscillation' in t):
            ared=tip[n]
            t0=np.linspace(0+(p*4215),2105+(p*4215),421)
            t02=np.linspace(2110+(p*4215),4215+(p*4215),421)
            t03=np.linspace(0+(p*4215),4215+(p*4215),42)
            t3=np.concatenate((t0,t02))
            areddecay=ared[0:int(len(ared)/2)]
            aredrise=ared[int(len(ared)/2):]
            amax=max(ared)
            amin=min(ared)
            delta=amax-amin
            ax1.plot(t3,tip[n],'r.')
            ##curve fit t1 & t2 exp rise
            try:
                popt1, pcov1 =curve_fit(rise,fittime,aredrise,maxfev=100000,p0=[delta,500,amin])
                Risefit=rise(fittime, *popt1)
                riseTcs.append(popt1[1])
                popt2, pcov2 =curve_fit(decay,fittime,areddecay,maxfev=100000,p0=[-delta,500,amin])
                decayTcs.append(popt2[1])
                Decayfit=decay(fittime, *popt2)
                with open(i+"/AverTCs.txt","a") as q:
                    q.write("Step\tRise TC\t\t\tDecay TC\t\t\tCov Rise\t\t\tCov Decay\n"
                            +stepnames[n]+"\t"+str(popt1[1])
                            +" ms\t"+str(popt2[1])+" ms"+str(pcov1[1])+"\t"+str(pcov2[1])+"\n")
                q.close()
                ax1.plot(t02,Risefit,'-',color="blue",label=stepnames[n].replace("Oscillation","")+"Rise TC ="+str(round(riseTcs[p]))+"ms")
                ax1.plot(t0,Decayfit,'-',color="black",label=stepnames[n].replace("Oscillation","")+"Decay TC ="+str(round(decayTcs[p]))+"ms")
#          plt.plot(t1,Decayfit,'yellow',label="T_decay ="+str(round(popt2[1]))+"ms")
            except:
                if(len(riseTcs)<p+1):
                    riseTcs.append(0.)
                if(len(decayTcs)<p+1):
                    decayTcs.append(0.)
                print("Exp Fit Problems Occured!!")
            ax2.plot(t3,sample[n],'.',ms=1,color="olive")
            ax1.set_xlabel("Time (ms)")
            ax4.plot(t0,np.ones(len(t0))*ldws[n],'g.')
            ax4.plot(t02,np.ones(len(t02))*lups[n],'g.')
            loadup.append(lups[n])
            loaddw.append(ldws[n])
            ax3.plot(t3,voltMean[n],".",color="orange")
            p+=1
            depthLoad.append(np.mean(voltMean[n][int(len(voltMean[n])/2):]))
            depthUnlod.append(np.mean(voltMean[n][:int(len(voltMean[n])/2)]))
        n+=1
    deltadep=[]
    deltaload=[]
    for t in range(len(depthLoad)):
        deltadep.append(depthLoad[t]-depthUnlod[t])
    ax4.grid('True',axis='y',color='green')
    ax1.grid('True',color='red')
    inset.set_yticks(np.linspace(round(min(riseTcs)),round(max(riseTcs)),5))
    l1,=inset.plot(deltadep,riseTcs,'r.-',label="Rise TC vs Depth")
    inset.set_xticks(deltadep)
    inset2.set_xticks(np.arange(0,450,50))
    #l2,=inset.plot(depthUnlod,decayTcs,'b.-',label="Decay TC vs Depth")
    l3,=inset2.plot(loadup,riseTcs,'b.-',label="Rise TC vs Load")
    #l4,=inset2.plot(loaddw,decayTcs,'g.-',label="Decay TC vs Load")
    inset.legend(handles=[l1,l3],fontsize='x-small',frameon='False')
#    ax1.legend()
    fig.tight_layout()
    plt.show()
    ###
##    for t in range(len(taular)):
##        print(t)
##        plt.plot(np.arange(1+t*(len(taular[t])),1+len(taular[t])+(len(taular[t])*t)),taular[t],'.')
##        plt.plot(np.arange(1+t*(len(taulad[t])),1+len(taulad[t])+(len(taulad[t])*t)),taulad[t],'.')
##        plt.plot(np.arange(1+t*(len(taular[t])),1+len(taular[t])+(len(taular[t])*t)),np.nanmean(taular[t])*np.ones(len(taular[t])),'-',label="Average Rise Step "+str(t+1)+" ="+str(round(np.nanmean(taular[t])))+" ms")
##        plt.plot(np.arange(1+t*(len(taulad[t])),1+len(taulad[t])+(len(taulad[t])*t)),np.nanmean(taulad[t])*np.ones(len(taulad[t])),'-',label="Average Decay Step "+str(t+1)+" ="+str(round(np.nanmean(taulad[t])))+" ms")
##    plt.legend()
##    plt.xlabel("Number of Cycle")
##    plt.ylabel("Time Constant (ms)")
##    plt.title(expdates[m]+" Time Constants of Each Cycle")
##    plt.show()
    m+=1
